#pragma once

#include <src/rtree/Node.h>
#include <src/rtree/Leaf.h>
#include <src/rtree/Index.h>
#include <src/rtree/BulkLoader.h>
#include <src/rtree/RTree.h>

#include <spatialindex/SpatialIndex.h>

#include "task.h"
#include "comper.h"
#include "worker.h"

#define TIMEOUT 1 // 0.1s

using namespace SpatialIndex;
using namespace std::literals;

enum OP
{
    DELETE = 0,
    INSERT,
    QUERY
};


// SpatialIndex::RTree::RTree *Srtree;
// SpatialIndex::RTree::RTree *Rrtree;


// IStorageManager *diskfile1;
// IStorageManager *diskfile2;

std::vector<size_t> thread_counter(64, 0);

std::string baseName1;
std::string baseName2;


class MyDataStream : public IDataStream
{
public:
    MyDataStream(std::string inputFile) : m_pNext(nullptr)
    {
        m_fin.open(inputFile.c_str());

        if (!m_fin)
            throw Tools::IllegalArgumentException("Input file not found.");

        readNextEntry();
    }

    ~MyDataStream() override
    {
        if (m_pNext != nullptr)
            delete m_pNext;
    }

    IData *getNext() override
    {
        if (m_pNext == nullptr)
            return nullptr;

        RTree::Data *ret = m_pNext;
        m_pNext = nullptr;
        readNextEntry();
        return ret;
    }

    bool hasNext() override
    {
        return (m_pNext != nullptr);
    }

    uint32_t size() override
    {
        throw Tools::NotSupportedException("Operation not supported.");
    }

    void rewind() override
    {
        if (m_pNext != nullptr)
        {
            delete m_pNext;
            m_pNext = nullptr;
        }

        m_fin.seekg(0, std::ios::beg);
        readNextEntry();
    }

    void readNextEntry()
    {
        id_type id;
        uint32_t op;
        double low[2], high[2];

        m_fin >> op >> id >> low[0] >> low[1] >> high[0] >> high[1];

        if (m_fin.good())
        {
            if (op != INSERT)
                throw Tools::IllegalArgumentException(
                    "The data input should contain insertions only.");

            Region r(low, high, 2);
            // m_pNext = new RTree::Data(sizeof(double), reinterpret_cast<uint8_t *>(low), r, id);
            m_pNext = new RTree::Data(0, nullptr, r, id);
            // Associate a bogus data array with every entry for testing purposes.
            // Once the data array is given to RTRee:Data a local copy will be created.
            // Hence, the input data array can be deleted a	fter this operation if not
            // needed anymore.
        }
    }

    std::ifstream m_fin;
    RTree::Data *m_pNext;
};


struct ContextValue
{
    id_type id1;
    id_type id2; 
    // Region r;
};

ofbinstream & operator>>(ofbinstream & m, Region & r)
{
    r.m_dimension = 2;
    r.m_pLow = new double[2];
    r.m_pHigh = new double[2];
    m >> r.m_pLow[0];
    m >> r.m_pLow[1];
    m >> r.m_pHigh[0];
    m >> r.m_pHigh[1];
    return m;
}
ifbinstream & operator<<(ifbinstream & m, const Region & r)
{
	m << r.m_pLow[0];
    m << r.m_pLow[1];
    m << r.m_pHigh[0];
    m << r.m_pHigh[1];
    return m;
}

ofbinstream & operator>>(ofbinstream & m, ContextValue & c)
{
	m >> c.id1;
    m >> c.id2;
    // m >> c.r;
    return m;
}

ifbinstream & operator<<(ifbinstream & m, const ContextValue & c)
{
	m << c.id1;
    m << c.id2;
    // m << c.r;
    return m;
}

obinstream & operator>>(obinstream & m, ContextValue & c)
{
	m >> c.id1;
    m >> c.id2;
    // m >> c.r;
    return m;
}

ibinstream & operator<<(ibinstream & m, const ContextValue & c)
{
	m << c.id1;
    m << c.id2;
    // m << c.r;
    return m;
}

typedef Task<ContextValue> SJTask; // spatial join task
typedef std::pair<id_type, id_type> SJData;

class SJComper: public Comper<SJTask, SJData>
{
public:

    SpatialIndex::RTree::RTree *Rrtree;
    SpatialIndex::RTree::RTree *Srtree;


    SJComper()
    {
        // auto start_time = std::chrono::steady_clock::now();

        // std::string baseName1 = "tree1"s;
        // std::string baseName2 = "tree2"s;

        IStorageManager *diskfile1 = StorageManager::createNewDiskStorageManager(baseName1, 4096);
        StorageManager::IBuffer *file1 = StorageManager::createNewRandomEvictionsBuffer(*diskfile1, 100, false);
        // applies a main memory random buffer on top of the persistent storage manager
        // (LRU buffer, etc can be created the same way).
        // MyDataStream stream1(str1);
        // Create and bulk load a new RTree with dimensionality 2, using "file" as
        // the StorageManager and the RSTAR splitting policy.
        // id_type indexIdentifier1;
        // ISpatialIndex *tree1 = RTree::createAndBulkLoadNewRTree(
        //     RTree::BLM_STR, stream1, *file1, utilization, 50, 50, 2, SpatialIndex::RTree::RV_RSTAR, indexIdentifier1);
        ISpatialIndex *tree1 = RTree::loadRTree(*file1, 1);
        Rrtree = dynamic_cast<RTree::RTree *>(tree1);

        // Create a new storage manager with the provided base name and a 4K page size.
        IStorageManager *diskfile2 = StorageManager::createNewDiskStorageManager(baseName2, 4096);
        StorageManager::IBuffer *file2 = StorageManager::createNewRandomEvictionsBuffer(*diskfile2, 100, false);
        // applies a main memory random buffer on top of the persistent storage manager
        // (LRU buffer, etc can be created the same way).

        // MyDataStream stream2(str2);
        // Create and bulk load a new RTree with dimensionality 2, using "file" as
        // the StorageManager and the RSTAR splitting policy.
        // id_type indexIdentifier2;
        // ISpatialIndex *tree2 = RTree::createAndBulkLoadNewRTree(
        //     RTree::BLM_STR, stream2, *file2, utilization, 20, 20, 2, SpatialIndex::RTree::RV_RSTAR, indexIdentifier2);

        ISpatialIndex *tree2 = RTree::loadRTree(*file2, 1);
        Srtree = dynamic_cast<RTree::RTree *>(tree2);

        // auto end_time = std::chrono::steady_clock::now();
        // std::cout << "Loading Rtree Time: " << (float)std::chrono::duration_cast<ms>(end_time - start_time).count() / 1000 << "s\n";


        // IStorageManager *diskfile1 = StorageManager::createNewMemoryStorageManager();
        // StorageManager::IBuffer *file1 = StorageManager::createNewRandomEvictionsBuffer(*diskfile1, 10, false);
        // MyDataStream stream1("./data3.txt");
        // id_type indexIdentifier1;
        // ISpatialIndex *tree1 = RTree::createAndBulkLoadNewRTree(
		// 	RTree::BLM_STR, stream1, *file1, 0.7, 100, 100, 2, SpatialIndex::RTree::RV_RSTAR, indexIdentifier1);
        // Rrtree = dynamic_cast<RTree::RTree *>(tree1);

        // IStorageManager *diskfile2 = StorageManager::createNewMemoryStorageManager();
        // StorageManager::IBuffer *file2 = StorageManager::createNewRandomEvictionsBuffer(*diskfile2, 10, false);
        // MyDataStream stream2("./data3.txt");
        // id_type indexIdentifier2;
        // ISpatialIndex *tree2 = RTree::createAndBulkLoadNewRTree(
		// 	RTree::BLM_STR, stream2, *file2, 0.7, 100, 100, 2, SpatialIndex::RTree::RV_RSTAR, indexIdentifier2);
        // Srtree = dynamic_cast<RTree::RTree *>(tree2);

        // auto end_time = std::chrono::steady_clock::now();
        // std::cout << "Loading Rtree Time: " << (float)std::chrono::duration_cast<ms>(end_time - start_time).count() / 1000 << "s\n";
    }

    bool isTimeout(auto start_time)
    {
        auto cur_time = std::chrono::steady_clock::now();
        float elapsed_time = (float)std::chrono::duration_cast<ms>(cur_time - start_time).count() / 1000;
        if (elapsed_time > TIMEOUT)
            return true;
        else 
            return false;
    }

    /*
    void SpatialJoinQuery(id_type id1, id_type id2, const Region& r, size_t &joinRes, auto start_time)
    {
        if (id1 < 0)
        {
            RTree::NodePtr n2 = Srtree->readNode(id2);
            for (uint32_t cChild2 = 0; cChild2 < n2->m_children; ++cChild2)
            {
                if (
                    !r.intersectsRegion(*(n2->m_ptrMBR[cChild2]))
                    )
                    continue;
                if (n2->m_level == 0)
                {
                    joinRes++;
                }
                else
                {
                    // Region rr = r.getIntersectingRegion(*(n2->m_ptrMBR[cChild2]));
                    // if (rr.m_pLow[0] != std::numeric_limits<double>::max())
                    if (!isTimeout(start_time))
                        SpatialJoinQuery(-1, n2->m_pIdentifier[cChild2], r, joinRes, start_time);
                    else
                    {
                        SJTask * t = new SJTask {-1, n2->m_pIdentifier[cChild2]};
                        add_task(t);
                    }
                }
            }
        }
        else if (id2 < 0)
        {
            RTree::NodePtr n1 = Rrtree->readNode(id1);
            for (uint32_t cChild1 = 0; cChild1 < n1->m_children; ++cChild1)
            {
                if (
                    !r.intersectsRegion(*(n1->m_ptrMBR[cChild1]))
                    )
                    continue;
                if (n1->m_level == 0)
                {
                    joinRes++;
                }
                else
                {
                    // Region rr = r.getIntersectingRegion(*(n1->m_ptrMBR[cChild1]));
                    // if (rr.m_pLow[0] != std::numeric_limits<double>::max())

                    if (!isTimeout(start_time))
                        SpatialJoinQuery(n1->m_pIdentifier[cChild1], -1, r, joinRes, start_time);
                    else
                    {
                        SJTask * t = new SJTask {n1->m_pIdentifier[cChild1], -1};
                        add_task(t);
                    }
                }
            }
        }
        else
        {
            RTree::NodePtr n1 = Rrtree->readNode(id1);
            RTree::NodePtr n2 = Srtree->readNode(id2);

            for (uint32_t cChild1 = 0; cChild1 < n1->m_children; ++cChild1)
            {
                // if (!r.intersectsRegion(*(n1->m_ptrMBR[cChild1])))
                // 	continue;

                for (uint32_t cChild2 = 0; cChild2 < n2->m_children; ++cChild2)
                {
                    if (
                        !n1->m_ptrMBR[cChild1]->intersectsRegion(*(n2->m_ptrMBR[cChild2]))
                        )
                        continue;
                    
                    if (n1->m_level == 0 && n2->m_level == 0)
                    {
                        joinRes++;
                    }
                    else if (n1->m_level == 0) // n2->m_level != 0
                    {
                        // Region rr = n1->m_ptrMBR[cChild1]->getIntersectingRegion(*(n2->m_ptrMBR[cChild2]));
                        // if (rr.m_pLow[0] != std::numeric_limits<double>::max())
                        if (!isTimeout(start_time))
                            SpatialJoinQuery(-1, n2->m_pIdentifier[cChild2], *n1->m_ptrMBR[cChild1], joinRes, start_time);
                        else
                        {
                            SJTask * t = new SJTask {-1, n2->m_pIdentifier[cChild2]};
                            add_task(t);
                        }
                    }
                    else if (n2->m_level == 0) // n1->m_level != 0
                    {
                        // Region rr = n1->m_ptrMBR[cChild1]->getIntersectingRegion(*(n2->m_ptrMBR[cChild2]));
                        // if (rr.m_pLow[0] != std::numeric_limits<double>::max())
                        if (!isTimeout(start_time))
                            SpatialJoinQuery(n1->m_pIdentifier[cChild1], -1, *n2->m_ptrMBR[cChild2], joinRes, start_time);
                        else
                        {
                            SJTask * t = new SJTask {n1->m_pIdentifier[cChild1], -1};
                            add_task(t);
                        }
                    }
                    else
                    {
                        // Region rr = n1->m_ptrMBR[cChild1]->getIntersectingRegion(*(n2->m_ptrMBR[cChild2]));
                        // if (rr.m_pLow[0] != std::numeric_limits<double>::max())

                        if (!isTimeout(start_time))
                            SpatialJoinQuery(n1->m_pIdentifier[cChild1], n2->m_pIdentifier[cChild2], r, joinRes, start_time);
                        else
                        {
                            SJTask * t = new SJTask {n1->m_pIdentifier[cChild1], 
                                                    n2->m_pIdentifier[cChild2]};
                            add_task(t);
                        }
                    }
                }
            }
        }
    }
    */


    void SpatialJoinQuery(id_type id1, id_type id2, const Region& r, size_t &joinRes, auto start_time)
    {
        if (id1 < 0)
        {
            th_wait_time -= get_time();
            RTree::NodePtr n2 = Srtree->readNode(id2);
            th_wait_time += get_time();

            uint32_t k = 0;

            while (k < n2->m_children && n2->m_ptrMBR[k]->m_pLow[0] <= r.m_pHigh[0])
            {
                if (r.m_pLow[0] <= n2->m_ptrMBR[k]->m_pHigh[0] &&
                    r.m_pLow[1] <= n2->m_ptrMBR[k]->m_pHigh[1] && 
                    r.m_pHigh[1] >= n2->m_ptrMBR[k]->m_pLow[1])
                {
                    if (n2->m_level == 0)
                    {
                        joinRes ++;
                    }
                    else
                    {
                        Region rr = r.getIntersectingRegion(*(n2->m_ptrMBR[k]));
                        SpatialJoinQuery(-1, n2->m_pIdentifier[k], rr, joinRes, start_time);
                    }
                }
                k++;
            }
        }
        else if (id2 < 0)
        {
            th_wait_time -= get_time();
            RTree::NodePtr n1 = Rrtree->readNode(id1);
            th_wait_time += get_time();

            uint32_t k = 0;

            while (k < n1->m_children && n1->m_ptrMBR[k]->m_pLow[0] <= r.m_pHigh[0])
            {
                if (r.m_pLow[0] <= n1->m_ptrMBR[k]->m_pHigh[0] &&
                    r.m_pLow[1] <= n1->m_ptrMBR[k]->m_pHigh[1] && 
                    r.m_pHigh[1] >= n1->m_ptrMBR[k]->m_pLow[1])
                {
                    if (n1->m_level == 0)
                    {
                        joinRes ++;
                    }
                    else
                    {
                        Region rr = r.getIntersectingRegion(*(n1->m_ptrMBR[k]));
                        SpatialJoinQuery(n1->m_pIdentifier[k], -1, rr, joinRes, start_time);
                    }
                }
                k++;
            }
        }
        else 
        {
            th_wait_time -= get_time();
            RTree::NodePtr n1 = Rrtree->readNode(id1);
            RTree::NodePtr n2 = Srtree->readNode(id2);
            th_wait_time += get_time();

            uint32_t cur_level = std::max(n1->m_level, n2->m_level);

            uint32_t i = 0, j = 0, k;
            while (i < n1->m_children && j < n2->m_children)
            {
                if (n1->m_ptrMBR[i]->m_pLow[0] < n2->m_ptrMBR[j]->m_pLow[0])
                {
                    RegionPtr region = n1->m_ptrMBR[i];
                    k = j;
                    while (k < n2->m_children && n2->m_ptrMBR[k]->m_pLow[0] <= region->m_pHigh[0])
                    {
                        if (region->m_pLow[1] <= n2->m_ptrMBR[k]->m_pHigh[1] && 
                            region->m_pHigh[1] >= n2->m_ptrMBR[k]->m_pLow[1])
                        {
                            if (n1->m_level == 0 && n2->m_level == 0)
                            {
                                joinRes ++;
                            }
                            else if (n1->m_level == 0)
                            {
                                Region rr = region->getIntersectingRegion(*(n2->m_ptrMBR[k]));
                                SpatialJoinQuery(-1, n2->m_pIdentifier[k], rr, joinRes, start_time);
                            }
                            else if (n2->m_level == 0)
                            {
                                Region rr = region->getIntersectingRegion(*(n2->m_ptrMBR[k]));
                                SpatialJoinQuery(n1->m_pIdentifier[i], -1, rr, joinRes, start_time);
                            }
                            else 
                            {
                                if (!isTimeout(start_time))
                                    SpatialJoinQuery(n1->m_pIdentifier[i], n2->m_pIdentifier[k], r, joinRes, start_time);
                                else
                                {
                                    SJTask * t = new SJTask {n1->m_pIdentifier[i], n2->m_pIdentifier[k]};
                                    add_task(t);
                                }
                            }
                        }
                        k++;
                    }
                    i++;
                }
                else
                {
                    RegionPtr region = n2->m_ptrMBR[j];
                    k = i;
                    while (k < n1->m_children && n1->m_ptrMBR[k]->m_pLow[0] <= region->m_pHigh[0])
                    {
                        if (region->m_pLow[1] <= n1->m_ptrMBR[k]->m_pHigh[1] && 
                            region->m_pHigh[1] >= n1->m_ptrMBR[k]->m_pLow[1])
                        {
                            if (n1->m_level == 0 && n2->m_level == 0)
                            {
                                joinRes ++;
                            }
                            else if (n1->m_level == 0)
                            {
                                Region rr = region->getIntersectingRegion(*(n1->m_ptrMBR[k]));
                                SpatialJoinQuery(-1, n2->m_pIdentifier[j], rr, joinRes, start_time);
                            }
                            else if (n2->m_level == 0)
                            {
                                Region rr = region->getIntersectingRegion(*(n1->m_ptrMBR[k]));
                                SpatialJoinQuery(n1->m_pIdentifier[k], -1, rr, joinRes, start_time);
                            }
                            else 
                            {
                                // SpatialJoinQueryWithSorting2(Stree, n1->m_pIdentifier[k], n2->m_pIdentifier[j], r, joinRes);
                                if (!isTimeout(start_time))
                                    SpatialJoinQuery(n1->m_pIdentifier[k], n2->m_pIdentifier[j], r, joinRes, start_time);
                                else
                                {
                                    SJTask * t = new SJTask {n1->m_pIdentifier[k], n2->m_pIdentifier[j]};
                                    add_task(t);
                                }
                            }
                        }
                        k++;
                    }
                    j++;
                }
            }
        }
    }

    virtual void compute(ContextT &context) override
    {
        auto start_time =  std::chrono::steady_clock::now();
        SpatialJoinQuery(context.id1, context.id2, Region{}, thread_counter[thread_id], start_time);
    }

    virtual bool task_spawn(SJData &data) override
	{
		SJTask *task = new SJTask;
        task->context.id1 = data.first;
        task->context.id2 = data.second;

        // double low[2] = {0, 0};
        // double high[2] = {100, 100};
        // Region r(low, high, 2);

        // task->context.r = r;
		return add_task(task);
	}

    virtual bool is_bigTask(SJTask *task) override
	{
		return false;
	}
};


class SJWorker: public Worker<SJComper>
{
public:
    SJWorker(int num_compers) : Worker(num_compers)
    {
    }

    virtual bool task_spawn(SJData &data) override
	{
		SJTask *task = new SJTask;
        task->context.id1 = data.first;
        task->context.id2 = data.second;

        // double low[2] = {0, 0};
        // double high[2] = {100, 100};
        // Region r(low, high, 2);

        // task->context.r = r;

		return add_task(task);
	}    

    void load_data(std::string const& rtree_name_1, std::string const& rtree_name_2)
    {
        baseName1 = rtree_name_1;
        baseName2 = rtree_name_2;
        double utilization = 0.7;
        uint32_t capacity1 = 100;
        uint32_t capacity2 = 100;

        // =========  Build First R-tree ========
        IStorageManager *diskfile1 = StorageManager::createNewDiskStorageManager(baseName1, 4096);
        // Create a new storage manager with the provided base name and a 4K page size.

        // StorageManager::IBuffer *file1 = StorageManager::createNewRandomEvictionsBuffer(*diskfile1, 10, false);
        ISpatialIndex *tree1 = RTree::loadRTree(*diskfile1, 1);
        SpatialIndex::RTree::RTree * Rrtree = dynamic_cast<RTree::RTree *>(tree1);

        // =========  Build Second R-tree ========
        IStorageManager *diskfile2 = StorageManager::createNewDiskStorageManager(baseName2, 4096);
        // StorageManager::IBuffer *file2 = StorageManager::createNewRandomEvictionsBuffer(*diskfile2, 10, false);
        ISpatialIndex *tree2 = RTree::loadRTree(*diskfile2, 1);
        SpatialIndex::RTree::RTree * Srtree = dynamic_cast<RTree::RTree *>(tree2);
        
        // do a few levels (=1) of spatial join
        // to generate some preliminary tasks
        // generate some subtasks for processing

        auto init_hash = [] (const id_type a, const id_type b) {
            return std::hash<id_type>{}(a) + std::hash<id_type>{}(b);
        };

        SpatialIndex::RTree::NodePtr n1 = Rrtree->readNode(Rrtree->m_rootID);
		SpatialIndex::RTree::NodePtr n2 = Srtree->readNode(Srtree->m_rootID);

        uint32_t i = 0, j = 0, k;
		while (i < n1->m_children && j < n2->m_children)
		{
			if (n1->m_ptrMBR[i]->m_pLow[0] < n2->m_ptrMBR[j]->m_pLow[0])
			{
				RegionPtr region = n1->m_ptrMBR[i];
				k = j;
				while (k < n2->m_children && n2->m_ptrMBR[k]->m_pLow[0] <= region->m_pHigh[0])
				{
					if (region->m_pLow[1] <= n2->m_ptrMBR[k]->m_pHigh[1] && 
						region->m_pHigh[1] >= n2->m_ptrMBR[k]->m_pLow[1])
					{
						if (init_hash(n1->m_pIdentifier[i], n2->m_pIdentifier[k]) % _num_workers == _my_rank)
                        {
                            data_array.push_back(new SJData{n1->m_pIdentifier[i], 
                                                    n2->m_pIdentifier[k]});
                        }
					}
					k++;
				}
				i++;
			}
			else
			{
				RegionPtr region = n2->m_ptrMBR[j];
				k = i;
				while (k < n1->m_children && n1->m_ptrMBR[k]->m_pLow[0] <= region->m_pHigh[0])
				{
					if (region->m_pLow[1] <= n1->m_ptrMBR[k]->m_pHigh[1] && 
						region->m_pHigh[1] >= n1->m_ptrMBR[k]->m_pLow[1])
					{
                        if (init_hash(n1->m_pIdentifier[k], n2->m_pIdentifier[j]) % _num_workers == _my_rank)
                        {
                            data_array.push_back(new SJData{n1->m_pIdentifier[k], 
                                                    n2->m_pIdentifier[j]});
                        }
					}
					k++;
				}
				j++;
			}
		}
    }

    virtual bool is_bigTask(SJTask *task) override
	{
		return false;
	}
};