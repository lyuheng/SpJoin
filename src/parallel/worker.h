#pragma once

#include <iostream>
#include "comper.h"
#include <unistd.h>
#include <omp.h>

#include <src/mpi/mpi_global.h>
#include <src/mpi/communication.h>

struct steal_plan
{
	int id; // steal from/to which server. negative number means receive
	int num; // steal batch size

    steal_plan(int id_, int num_): id(id_), num(num_) {}
    steal_plan() {}

};

obinstream & operator>>(obinstream & m, steal_plan & s)
{
	m >> s.id;
	m >> s.num;
	return m;
}

ibinstream & operator<<(ibinstream & m, const steal_plan & s)
{
	m << s.id;
	m << s.num;
	return m;
}

ofbinstream & operator>>(ofbinstream & m, steal_plan & s)
{
	m >> s.id;
	m >> s.num;
	return m;
}

ifbinstream & operator<<(ifbinstream & m, const steal_plan & s)
{
	m << s.id;
	m << s.num;
	return m;
}

template <class ComperT>
class Worker
{
public:
    typedef typename ComperT::TaskType TaskT;
    typedef typename ComperT::DataType DataT;
    typedef std::deque<TaskT *> TaskQ;
    typedef std::stack<DataT *> DataStack;

    // Dynamic array of compers
    ComperT *compers = nullptr;
    // Contains all data loaded from file
    std::vector<DataT *> data_array;
    // Contains pointers to data-array without initialized tasks pointers
    DataStack *data_stack;
    // Regular tasks queue
    TaskQ *Qreg;
    // Big tasks queue
    TaskQ *Qbig;

    typedef std::chrono::_V2::steady_clock::time_point timepoint;
    // Worker's start time
    timepoint start_time = steady_clock::now();
    // To save latest elapsed time in seconds, in order to find exact duration for each method
    float latest_elapsed_time = 0;

    float init_time;
    float load_data_time;

    // Save files sequence number when spill big tasks to disk
    std::vector<size_t> big_files_seq;
    // Save files sequence number when spill reg tasks to disk
    std::vector<size_t> reg_files_seq;

    std::vector<std::string> big_files_names;
    std::vector<std::string> reg_files_names;

    // Init files seq with 1 for each thread.
    Worker(int comper_num) : 
        big_files_seq(comper_num, 1), 
        reg_files_seq(comper_num, 1), 
        big_files_names(comper_num), 
        reg_files_names(comper_num)
    {
        // Create disk buffer dir
        TASK_DISK_BUFFER_DIR = "buffered_tasks" + std::to_string(_my_rank);
        recursive_mkdir(TASK_DISK_BUFFER_DIR.c_str());

        global_data_stack = data_stack = new std::stack<DataT *>;

        global_Qreg = Qreg = new TaskQ;
        global_Qbig = Qbig = new TaskQ;
        num_compers = comper_num;
        global_end_label = false;
    }

    virtual ~Worker()
    {
        for (int i = 0; i < data_array.size(); i++)
        {
            delete data_array[i];
        }

        if (compers)
            delete[] compers;
        delete Qreg;
        delete Qbig;
    }

    // UDF1: read data from file_path, and insert into data_array
    // virtual void load_data(const std::string &file_path) {}

    // UDF2
    virtual bool task_spawn(DataT &data) = 0;

    // UDF3
    virtual bool is_bigTask(TaskT *task)
    {
        return false;
    }

    // Insert some tasks into Qreg before spawn compers, so compers have some tasks to work on
    void initialize_tasks()
    {

        // 1- Spawn some tasks from data_array
        size_t _size = std::min(Qreg_capacity, data_array.size());

#pragma omp parallel for schedule(dynamic, 1) num_threads(num_compers)
        for (int i = 0; i < _size; i++)
        {
            task_spawn(*(data_array[i]));
        }
        // 2- Add the >>rest<< of data_array to a data stack, to be used by comper when spawn
        for (int i = data_array.size() - 1; i >= _size; i--)
        {
            data_stack->push(data_array[i]);
        }
    }

    bool add_task(TaskT *task)
    {
        if (is_bigTask(task))
        {
            add_bigTask(task);
            return true;
        }

        add_regTask(task);
        return false;
    }

    void add_bigTask(TaskT *task)
    {
        Qbig_mtx.lock();
        // Check if spill is needed
        if (Qbig->size() == Qbig_capacity)
        {
            spill_Qbig();
            Qbig_mtx.lock();
        }
        Qbig->push_back(task);
        Qbig_mtx.unlock();
    }

    void add_regTask(TaskT *task)
    {
        Qreg_mtx.lock();
        // Check if spill is needed.
        if (Qreg->size() == Qreg_capacity)
        {
            spill_Qreg();
            Qreg_mtx.lock();
        }
        Qreg->push_back(task);
        Qreg_mtx.unlock();
    }

    void spill_Qbig()
    {
        int i = 0;
        std::queue<TaskT *> collector;
        while (i < BT_TASKS_PER_FILE && !Qbig->empty())
        {
            // Get task at the tail
            TaskT *t = Qbig->back();
            Qbig->pop_back();
            collector.push(t);
            i++;
        }
        Qbig_mtx.unlock();

        if (!collector.empty())
        {
            int thread_id = omp_get_thread_num();
            set_bigTask_fname(thread_id);
            ifbinstream bigTask_out(big_files_names[thread_id].c_str());

            while (!collector.empty())
            {
                TaskT *t = collector.front();
                collector.pop();
                // Stream to file
                bigTask_out << t;
                // Release from memory
                delete t;
            }

            bigTask_out.close();
            global_Lbig.enqueue(big_files_names[thread_id]);
            global_Lbig_num++;
        }
    }

    void spill_Qreg()
    {
        int i = 0;
        std::queue<TaskT *> collector;
        while (i < RT_TASKS_PER_FILE && !Qreg->empty())
        {
            // Get task at the tail
            TaskT *t = Qreg->back();
            Qreg->pop_back();
            collector.push(t);
            i++;
        }
        Qreg_mtx.unlock();

        if (!collector.empty())
        {
            int thread_id = omp_get_thread_num();
            set_regTask_fname(thread_id);
            ifbinstream regTask_out(reg_files_names[thread_id].c_str());

            while (!collector.empty())
            {
                TaskT *t = collector.front();
                collector.pop();
                // Stream to file
                regTask_out << t;
                // Release from memory
                delete t;
            }
            regTask_out.close();
            global_Lreg.enqueue(reg_files_names[thread_id]);
            global_Lreg_num++;
        }
    }

    void set_bigTask_fname(const int thread_id)
    {
        // Reset filename
        big_files_names[thread_id] = "";
        big_files_names[thread_id] += TASK_DISK_BUFFER_DIR + "/w_" + std::to_string(thread_id) + 
                                        "_" + std::to_string(big_files_seq[thread_id]) + "_bt";
        big_files_seq[thread_id]++;
    }

    void set_regTask_fname(const int thread_id)
    {
        // Reset filename
        reg_files_names[thread_id] = "";
        reg_files_names[thread_id] += TASK_DISK_BUFFER_DIR + "/w_" + std::to_string(thread_id) + 
                                        "_" + std::to_string(reg_files_seq[thread_id]) + "_rt";
        reg_files_seq[thread_id]++;
    }

    size_t get_remaining_task_num()
	//not counting number of active tasks in memory (for simplicity)
	{
		Qreg_mtx.lock();
		int q_regtask_size = Qreg->size();
		Qreg_mtx.unlock();
		return q_regtask_size + global_Lreg_num.load() * RT_TASKS_PER_FILE;
	}

    bool regTask_file2vec(std::vector<TaskT *> & tvec)
	{
		std::string file;
		bool succ = global_Lreg.dequeue(file);
		//@@@@@@@@@@@@@@@@@
		//if(succ) cout<<"!!!!!!!!!!! there is big file !!!!!!!!"<<endl;
		if (!succ) return false; //"global_file_list" is empty
		else
		{
			global_Lreg_num --;
			ofbinstream in(file.c_str());
			while (!in.eof())
			{
				TaskT * task = new TaskT;
				in >> *task;
				tvec.push_back(task);
			}
			in.close();
			//------
			if (remove(file.c_str()) != 0) {
				std::cout << "Error removing file: " << file << std::endl;
				perror("Error printed by perror");
			}
			return true;
		}
	}

    bool file2regTask_queue()
	{
		std::string file;
		bool succ = global_Lreg.dequeue(file);
		if (!succ) return false; //"global_bigTask_fileList" is empty
		else
		{
			global_Lreg_num --;
			ofbinstream in(file.c_str());
			while(!in.eof())
			{
				TaskT* task;
				in >> task;
				Qreg->push_back(task);
			}
			in.close();

			if (remove(file.c_str()) != 0) {
				std::cout << "Error removing file: " << file << std::endl;
				perror("Error printed by perror");
			}
			return true;
		}
	}

    struct max_heap_entry
	{
		size_t num_remain;
		int rank;

		bool operator<(const max_heap_entry& o) const
		{
			return num_remain < o.num_remain;
		}
	};

	struct min_heap_entry
	{
		size_t num_remain;
		int rank;

		bool operator<(const min_heap_entry& o) const
		{
			return num_remain > o.num_remain;
		}
	};

    bool steal_planning()
    {
        std::vector<steal_plan> my_single_steal_list; // if my_single_steal_list[i] = 3, I need to steal a tasks from Worker 3's global queue
		std::vector<int> my_batch_steal_list; // if my_batch_steal_list[i] = 3, I need to steal a task-file from Worker 3's global file list
        int avg_num;
		//====== set my_steal_list
		if (_my_rank != MASTER_RANK)
		{
			send_data(get_remaining_task_num(), MASTER_RANK, STATUS_CHANNEL);
			recv_data<std::vector<steal_plan>>(MASTER_RANK, STATUS_CHANNEL, my_single_steal_list);
			recv_data<std::vector<int>>(MASTER_RANK, STATUS_CHANNEL, my_batch_steal_list);
            recv_data<int>(MASTER_RANK, STATUS_CHANNEL, avg_num);
		}
		else
        {
            std::vector<size_t> remain_vec(_num_workers);
			for (int i = 0; i < _num_workers; i++)
			{
				if(i != MASTER_RANK)
					remain_vec[i] = recv_data<size_t>(i, STATUS_CHANNEL);
				else
					remain_vec[i] = get_remaining_task_num();
			}
			//------
			std::priority_queue<max_heap_entry> max_heap;
			std::priority_queue<min_heap_entry> min_heap;
			avg_num = std::ceil((float)std::accumulate(remain_vec.begin(), remain_vec.end(), 0)/remain_vec.size());

            for (int i = 0; i < _num_workers; i++)
			{
				if (remain_vec[i] > avg_num)
				{
					max_heap_entry en;
					en.num_remain = remain_vec[i];
					en.rank = i;
					max_heap.push(en);
				}
				else if (remain_vec[i] < avg_num)
				{
					min_heap_entry en;
					en.num_remain = remain_vec[i];
					en.rank = i;
					min_heap.push(en);
				}
			}

            std::vector<int> steal_num(_num_workers, 0); // steal_num[i] = Worker i currently will steal how many tasks
			std::vector<std::vector<steal_plan>> single_steal_lists(_num_workers); //steal_list[i] = {j{id,num}, ...} means Worker_i will steal a j.num tasks.. from Worker j.id, ...
			std::vector<std::vector<int>> batch_steal_lists(_num_workers); //steal_list[i] = {j, k, ...} means Worker_i will steal a task-file from each Worker j, k, ...
			int total_plans_num = 0;

            while (!max_heap.empty() && !min_heap.empty())
			{
				max_heap_entry max = max_heap.top();
				max_heap.pop();
				min_heap_entry min = min_heap.top();
				min_heap.pop();
				if (avg_num - min.num_remain > RT_TASKS_PER_FILE
					&& max.num_remain - avg_num > RT_TASKS_PER_FILE)
				{
                    // both has a gap >= task-batchsize, steal file
					max.num_remain -= RT_TASKS_PER_FILE;
					min.num_remain += RT_TASKS_PER_FILE;
					steal_num[min.rank] += RT_TASKS_PER_FILE;
					total_plans_num += RT_TASKS_PER_FILE; // only for printing

					//a negative tag (-x-1) means receiving
					batch_steal_lists[min.rank].push_back(-max.rank - 1);
					batch_steal_lists[max.rank].push_back(min.rank);
					//---
					if (max.num_remain > avg_num) max_heap.push(max);
					if (steal_num[min.rank] < MAX_STEAL_TASK_NUM && min.num_remain < avg_num)
						min_heap.push(min);
				}
				else
				{
                    //steal n task; n < BIG_TASK_FLUSH_BATCH
					int steal_batch = std::min((max.num_remain - avg_num), (avg_num - min.num_remain));
					max.num_remain -= steal_batch;
					min.num_remain += steal_batch;
					steal_num[min.rank] += steal_batch;
					total_plans_num += steal_batch;

					single_steal_lists[min.rank].push_back(steal_plan(-max.rank - 1, steal_batch));
					single_steal_lists[max.rank].push_back(steal_plan(min.rank, steal_batch));
					//---
					if(max.num_remain > avg_num) max_heap.push(max);
					if(steal_num[min.rank] < MAX_STEAL_TASK_NUM &&
							min.num_remain < avg_num)
						min_heap.push(min);
				}
			}
            // if (total_plans_num > 0)
            //     std::cout << total_plans_num << " stealing plans generated at the master\n";
			
            for (int i = 0; i < _num_workers; i++)
			{
				if (i == _my_rank)
				{
					single_steal_lists[i].swap(my_single_steal_list);
					batch_steal_lists[i].swap(my_batch_steal_list);
				}
				else
				{
					send_data(single_steal_lists[i], i, STATUS_CHANNEL);
					send_data(batch_steal_lists[i], i, STATUS_CHANNEL);
                    send_data(avg_num, i, STATUS_CHANNEL);
				}
			}
        }
        if (my_single_steal_list.empty() && my_batch_steal_list.empty()) 
            return false;
        
        if (!my_batch_steal_list.empty())
		{
			for (int i = 0; i < my_batch_steal_list.size(); i++)
			{
				int other = my_batch_steal_list[i];
				if(other < 0)
				{
					std::vector<TaskT> tvec;
					recv_data<std::vector<TaskT>>(-other-1, STATUS_CHANNEL, tvec);
					if(!tvec.empty())
					{
                        int thread_id = omp_get_thread_num();
						set_regTask_fname(thread_id);
						ifbinstream out(reg_files_names[thread_id].c_str());
						//------
						for(int i = 0; i < tvec.size(); i++)
						{
							out << tvec[i];
						}
						out.close();
						//------
						//register with "global_file_list"
						global_Lreg.enqueue(reg_files_names[thread_id].c_str());
						global_Lreg_num ++;
					}
				}
				else
				{
					std::vector<TaskT *> tvec;
					if(get_remaining_task_num() > avg_num)
					//check this since time has passed, and more tasks may have been processed
					//send empty tvec if no longer a task heavy-hitter
						regTask_file2vec(tvec);
					send_data(tvec, other, STATUS_CHANNEL); //send even if it's empty
					for(int i = 0; i < tvec.size(); i++)
						delete tvec[i];
				}
			}
		}


        if (!my_single_steal_list.empty())
		{
			for (int i = 0; i < my_single_steal_list.size(); i++)
			{
				int other = my_single_steal_list[i].id;
				if (other < 0)
				{
					std::vector<TaskT *> tvec;
					recv_data<std::vector<TaskT* >>(-other-1, STATUS_CHANNEL, tvec);
					for(int i = 0; i < tvec.size(); i++)
						add_regTask(tvec[i]);
				}
				else
				{
					int steal_num = my_single_steal_list[i].num;
					std::vector<TaskT *> tvec;
					for (int i = 0; i < steal_num; i++)
					{
						if (get_remaining_task_num() <= avg_num) break;
						//check this since time has passed, and more tasks may have been processed
						//send empty task-vec if no longer a task heavy-hitter

						TaskT * task = NULL;
						
						Qreg_mtx.lock();
						if(Qreg->size() <= MAX_STEAL_TASK_NUM)
							file2regTask_queue();

						if(!Qreg->empty())
						{
							task = Qreg->front();
							Qreg->pop_front();
						}
						Qreg_mtx.unlock();
						if (task != NULL) tvec.push_back(task);
					}
					send_data(tvec, other, STATUS_CHANNEL); //send even if it's empty
					for(int i = 0; i < tvec.size(); i++)
						delete tvec[i];
				}
			}
		}
        return true;
    }

    void create_compers()
    {
        compers = new ComperT[num_compers];
        for (int i = 0; i < num_compers; i++)
        {
            compers[i].start(i);
        }
    }

    void status_sync(bool sth2steal)
    {
        bool worker_idle = (sth2steal == false) && 
            (global_num_idle.load(std::memory_order_relaxed) == num_compers);
        
        if(_my_rank != MASTER_RANK)
        {
            send_data(worker_idle, MASTER_RANK, STATUS_CHANNEL);
            bool all_idle = recv_data<bool>(MASTER_RANK, STATUS_CHANNEL);
            if (all_idle)
                global_end_label = true;
        }
        else
        {
            bool all_idle = worker_idle;
            for (int i = 0; i < _num_workers; i++)
            {
                if(i != MASTER_RANK) 
                    all_idle = (recv_data<bool>(i, STATUS_CHANNEL) && all_idle);
            }
            if (all_idle) 
                global_end_label = true;
            for (int i = 0; i < _num_workers; i++)
            {
                if(i != MASTER_RANK) 
                    send_data(all_idle, i, STATUS_CHANNEL);
            }
        }
    }

    // Program entry point
    void run()
    {
        assert(RT_TASKS_PER_FILE <= Qreg_capacity);
        assert(BT_TASKS_PER_FILE <= Qbig_capacity);
        // Kernel_app's load_data method add tasks directly to both queues not to data_array.
        if (data_array.size() > 0)
        {
            // Initialize some tasks
            auto start = steady_clock::now();
            initialize_tasks();
            auto end = steady_clock::now();
            init_time = (float)duration_cast<milliseconds>(end - start).count() / 1000;
            std::cout << "initialize_tasks() execution time:" << init_time << std::endl;
        }

        // Setup computing threads
        create_compers();

        // Call status_sync() periodically
        while (global_end_label == false)
        {
            bool sth2steal = steal_planning();
            status_sync(sth2steal);
            //------
            //reset idle status of Worker, compers will add back if idle
            mtx_go.lock();
            ready_go = true;
            cv_go.notify_all(); //release threads to compute tasks
            mtx_go.unlock();
            // Avoid busy-checking
            usleep(WAIT_TIME_WHEN_IDLE);
        }


        double avg = 0;
        for (int i = 0; i < 32; i++) {
            avg += th_wait_time_g[i];
        }
        std::cout << "thread avg wait time = " << avg / 32 << std::endl; 
    }
};