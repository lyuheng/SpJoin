#include "src/parallel/spatialjoin.h"

int main(int argc, char *argv[])
{
    if (argc != 4) {
        exit(-1);
    }
    init_worker(&argc, &argv);
    std::cout << "Rank: " << _my_rank << std::endl; 

    SJWorker worker(atoi(argv[1]));

    std::string rtree_name_1 = std::string(argv[2]);
    std::string rtree_name_2 = std::string(argv[3]);

    worker.load_data(rtree_name_1, rtree_name_2);

    auto start_time = std::chrono::steady_clock::now();

    worker.run();

    auto end_time = std::chrono::steady_clock::now();

    std::cout << "Elasped Time: " << (float)std::chrono::duration_cast<ms>(end_time - start_time).count() / 1000 << "s\n";

    if (_my_rank == MASTER_RANK)
    {
        size_t ttl_result = std::accumulate(thread_counter.begin(), thread_counter.end(), (size_t)0);
        std::cout << "Results by each worker: {" << ttl_result; 
        for (int i = 0; i < _num_workers; ++i)
        {
            if (i != MASTER_RANK)
            {
                size_t found = recv_data<size_t>(i, RESULT_CHANNEL);
                std::cout << ", " << found;
                ttl_result += found;
            }
        }
        std::cout << "}\n";
        std::cout << "Join Result: " << ttl_result << std::endl;
    }
    else
    {
        for (int i = 0; i < _num_workers; ++i)
        {
            if (i != MASTER_RANK)
            {
                send_data(std::accumulate(thread_counter.begin(), thread_counter.end(), (size_t)0), MASTER_RANK, RESULT_CHANNEL);
            }
        }
    }

    worker_finalize();
}