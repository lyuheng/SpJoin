#pragma once

#include "meta.h"
#include "distsys.h"
#include "engine.h"
#include "aggregator.h"
#include "workload.h"

#include <numa.h>

using namespace SpatialIndex;
typedef RTree::RTree RTreeType;

class DistributedApplication {
public:
    constexpr static size_t NODE_PAIR_SIZE  = sizeof(id_type) * 2;
    int num_threads_;
    int num_sockets_;

    int max_depth_;
    int max_depth_R_;
    int max_depth_S_; 
    std::string rtree_name_1_;
    std::string rtree_name_2_;

    int node_id_;
    int num_nodes_;

    int num_partitions_;

    EmbeddingExplorationEngine * engines_[MAX_NUM_SOCKETS];
    WorkloadDistributer * workload_distributer_;


    size_t memory_size_per_level_; // the amount of memory allocated for the chunk of a level, default: 2GB

    void init_execution_engines();
    void finalize_execution_engines();

    void partition_workload();
    void partition_workload_sequential();

    DistributedApplication( std::string const& rtree_name_1,
                            std::string const& rtree_name_2,
                            int num_threads,
                            size_t memory_size_per_level = (size_t) 8 * 1024 * 1024 * 1024
                        );
    virtual ~DistributedApplication();

    void run();
    void load_Rtrees();

    // helper functions
    int get_num_threads() {
        return num_threads_;
    }
    int get_num_sockets() {
        return num_sockets_;
    }

    Aggregator<size_t> * count;
    RTreeType * Rrtree[MAX_NUM_SOCKETS + 1];
    RTreeType * Srtree[MAX_NUM_SOCKETS + 1];
    uint32_t * num_local_nodepair_per_partition_[MAX_NUM_SOCKETS];
    id_type * local_nodepair_[MAX_NUM_SOCKETS];
};

DistributedApplication::DistributedApplication( std::string const& rtree_name_1,
                                                std::string const& rtree_name_2,
                                                int num_threads, 
                                                size_t memory_size_per_level
                                            )
{
    node_id_ = DistributedSys::get_instance()->get_node_id();
    num_nodes_ = DistributedSys::get_instance()->get_num_nodes();

    rtree_name_1_ = rtree_name_1;
    rtree_name_2_ = rtree_name_2;

    int num_threads_limit = numa_num_configured_cpus() - 2 * numa_num_configured_nodes();
    if (num_threads == -1) {
        num_threads_ = num_threads_limit;
    } else {
        if (num_threads <= 0 || num_threads > num_threads_limit) {
            fprintf(stderr, "The number of threads should be between 1 and %d\n",
                    num_threads_limit);
            exit(-1);
        }
        num_threads_ = num_threads;
    }
    // int num_threads_limit_per_socket = numa_num_configured_cpus() / numa_num_configured_nodes() - 2;
    // if (num_threads_ <= num_threads_limit_per_socket) {
    //     num_sockets_ = 1;
    // } else {
        num_sockets_ = 8; // 8; //numa_num_configured_nodes();
        num_threads_ = 32 - 2 * num_sockets_;
        // assert(num_threads_ % num_sockets_ == 0);
    // }

    memory_size_per_level_ = memory_size_per_level / num_sockets_;
    printf("Chunk size per level: %.3f MB\n", memory_size_per_level / 1024. / 1024.);
    printf("Number of used sockets: %d\n", num_sockets_);


    num_partitions_ = num_nodes_ * num_sockets_;

    std::cout << "num_partitions_ = " << num_partitions_ << std::endl;


    // * load Rtree and Stree
    load_Rtrees();
    max_depth_ = std::max(Rrtree[0]->m_stats.getTreeHeight(), 
                            Srtree[0]->m_stats.getTreeHeight()) - 2;
    max_depth_R_ = Rrtree[0]->m_stats.getTreeHeight() - 2;
    max_depth_S_ = Srtree[0]->m_stats.getTreeHeight() - 2;
    
    std::cout << "max depth = " << max_depth_ << std::endl;

    partition_workload();
    // partition_workload_sequential();

    count = new Aggregator<size_t>(num_sockets_, num_threads_);

    init_execution_engines();

}

DistributedApplication::~DistributedApplication() {

    finalize_execution_engines();
    delete count;

    for (int s_i = 0; s_i < num_sockets_; ++s_i) {
        int partition_id = node_id_ * num_sockets_ + s_i;
        numa_free(local_nodepair_[s_i], NODE_PAIR_SIZE * num_local_nodepair_per_partition_[s_i][partition_id]);
        numa_free(num_local_nodepair_per_partition_[s_i], sizeof(uint32_t) * num_partitions_);
    }

    // TODO: destory Rrtree and Srtree
    // * maybe later
}


void DistributedApplication::init_execution_engines()
{
    lli max_num_embeddings_global_queue = (memory_size_per_level_/ 16) / sizeof(CompactExtendableEmbedding);
    std::cout << "          max_num_embeddings_global_queue = " << max_num_embeddings_global_queue << std::endl;
    lli global_disk_data_buffer_size = memory_size_per_level_/ 2 / 16 * 15; // we need 2 individual buffers for S-Rtree and R-Rtree

    for (int s_i = 0; s_i < num_sockets_; ++ s_i) {
        printf("Allocating the execution engine (size: %u) on node %d...\n", 
                sizeof(EmbeddingExplorationEngine), s_i);
        int partition_id = node_id_ * num_sockets_ + s_i;
        void * ptx = numa_alloc_onnode(sizeof(EmbeddingExplorationEngine), s_i/ENGINE_PER_SOCKET);

        printf("Check point A.\n");
        engines_[s_i] = new(ptx) EmbeddingExplorationEngine(
                max_depth_,
                max_depth_R_,
                max_depth_S_,
                max_num_embeddings_global_queue,
                global_disk_data_buffer_size,
                num_partitions_,
                s_i,
                num_sockets_,
                num_threads_ / num_sockets_,
                CHUNK_SIZE,
                Rrtree[s_i + 1],
                Srtree[s_i + 1],
                count
            );
        printf("Check point B.\n");
    }

    for (int s_i = 0; s_i < num_sockets_; ++ s_i) {
        engines_[s_i]->set_engines_on_the_same_node(engines_);
    }
    printf("Finish initializing engines...\n");

    workload_distributer_ = new WorkloadDistributer(node_id_, 
                                                    num_nodes_, 
                                                    num_sockets_,
                                                    num_local_nodepair_per_partition_,
                                                    local_nodepair_
                                                );
    printf("Finish initializing Workload Distributer ...\n");
}

void DistributedApplication::finalize_execution_engines()
{
    for (int s_i = 0; s_i < num_sockets_; ++ s_i) {
        void * ptx = engines_[s_i];
        engines_[s_i]->~EmbeddingExplorationEngine();
        numa_free(ptx, sizeof(EmbeddingExplorationEngine));
    }
    delete workload_distributer_;
}

void DistributedApplication::partition_workload_sequential()
{
    SpatialIndex::RTree::NodePtr n1 = Rrtree[0]->readNode(Rrtree[0]->m_rootID);
    SpatialIndex::RTree::NodePtr n2 = Srtree[0]->readNode(Srtree[0]->m_rootID);

    uint32_t i = 0, j = 0, k;
    uint32_t num_ttl_nodepair = 0;

    while (i < n1->m_children && j < n2->m_children) {
        if (n1->m_ptrMBR[i]->m_pLow[0] < n2->m_ptrMBR[j]->m_pLow[0]) {
            RegionPtr region = n1->m_ptrMBR[i];
            k = j;
            while (k < n2->m_children && n2->m_ptrMBR[k]->m_pLow[0] <= region->m_pHigh[0]) {
                if (region->m_pLow[1] <= n2->m_ptrMBR[k]->m_pHigh[1] && 
                    region->m_pHigh[1] >= n2->m_ptrMBR[k]->m_pLow[1]) {

                    num_ttl_nodepair ++;
                    
                }
                k++;
            }
            i++;
        }
        else {
            RegionPtr region = n2->m_ptrMBR[j];
            k = i;
            while (k < n1->m_children && n1->m_ptrMBR[k]->m_pLow[0] <= region->m_pHigh[0]) {
                if (region->m_pLow[1] <= n1->m_ptrMBR[k]->m_pHigh[1] &&
                    region->m_pHigh[1] >= n1->m_ptrMBR[k]->m_pLow[1]) {

                    num_ttl_nodepair++;
                }
                k++;
            }
            j++;
        }
    }

    uint32_t num_avg_nodepair = num_ttl_nodepair % num_partitions_ == 0 ? 
                                num_ttl_nodepair / num_partitions_ :
                                num_ttl_nodepair / num_partitions_ + 1;
    uint32_t num_nodepair_last_partition = num_ttl_nodepair - num_avg_nodepair * (num_partitions_ - 1);
    // uint32_t prefix_sum[MAX_NUM_PARTITIONS];
    // for (int i = 0; i < num_partitions_; ++i) {
    //     prefix_sum[i] = i * num_avg_nodepair;
    // }
    // prefix_sum[num_partitions_] = num_ttl_nodepair;

    // std::cout << num_avg_nodepair << " " << num_nodepair_last_partition << std::endl;

    uint32_t counter[MAX_NUM_SOCKETS];
    memset(counter, 0, sizeof(uint32_t) * num_sockets_);

    for (int s_i = 0; s_i < num_sockets_; ++s_i) {

        num_local_nodepair_per_partition_[s_i] = (uint32_t *) numa_alloc_onnode(
                        sizeof(uint32_t) * num_partitions_, s_i/ENGINE_PER_SOCKET);

        for (int idx = 0; idx < num_partitions_; ++idx) {
            if (idx != num_partitions_ - 1)
                num_local_nodepair_per_partition_[s_i][idx] = num_avg_nodepair;
            else
                num_local_nodepair_per_partition_[s_i][idx] = num_nodepair_last_partition;
        }

        int partition_id = node_id_ * num_sockets_ + s_i;
        if (partition_id != num_partitions_ - 1) {
            local_nodepair_[s_i] = (id_type *) numa_alloc_onnode(
                                NODE_PAIR_SIZE * num_avg_nodepair, s_i/ENGINE_PER_SOCKET);

        } else {
            local_nodepair_[s_i] = (id_type *) numa_alloc_onnode(
                                NODE_PAIR_SIZE * num_nodepair_last_partition, s_i/ENGINE_PER_SOCKET);

        }
        
        i = 0;
        j = 0;
        num_ttl_nodepair = 0;
        while (i < n1->m_children && j < n2->m_children) {
            if (n1->m_ptrMBR[i]->m_pLow[0] < n2->m_ptrMBR[j]->m_pLow[0]) {
                RegionPtr region = n1->m_ptrMBR[i];
                k = j;
                while (k < n2->m_children && n2->m_ptrMBR[k]->m_pLow[0] <= region->m_pHigh[0]) {
                    if (region->m_pLow[1] <= n2->m_ptrMBR[k]->m_pHigh[1] && 
                        region->m_pHigh[1] >= n2->m_ptrMBR[k]->m_pLow[1]) {

                        int p_i = num_ttl_nodepair / num_avg_nodepair;
                        if (p_i % num_sockets_ == s_i && p_i / num_sockets_ == node_id_) {
                            local_nodepair_[s_i][counter[s_i]++] = n1->m_pIdentifier[i];
                            local_nodepair_[s_i][counter[s_i]++] = n2->m_pIdentifier[k];
                        }
                        num_ttl_nodepair ++;
                    }
                    k++;
                }
                i++;
            }
            else {
                RegionPtr region = n2->m_ptrMBR[j];
                k = i;
                while (k < n1->m_children && n1->m_ptrMBR[k]->m_pLow[0] <= region->m_pHigh[0]) {
                    if (region->m_pLow[1] <= n1->m_ptrMBR[k]->m_pHigh[1] &&
                        region->m_pHigh[1] >= n1->m_ptrMBR[k]->m_pLow[1]) {

                        int p_i = num_ttl_nodepair / num_avg_nodepair;
                        if (p_i % num_sockets_ == s_i && p_i / num_sockets_ == node_id_) {
                            local_nodepair_[s_i][counter[s_i]++] = n1->m_pIdentifier[k];
                            local_nodepair_[s_i][counter[s_i]++] = n2->m_pIdentifier[j];
                        }
                        num_ttl_nodepair++;
                    }
                    k++;
                }
                j++;
            }
        }
    }
}

void DistributedApplication::partition_workload()
{
    auto INIT_HASH = [] (const id_type a, const id_type b) {
        return std::hash<id_type>{}(a) + std::hash<id_type>{}(b);
    };

    SpatialIndex::RTree::NodePtr n1 = Rrtree[0]->readNode(Rrtree[0]->m_rootID);
    SpatialIndex::RTree::NodePtr n2 = Srtree[0]->readNode(Srtree[0]->m_rootID);

    for (int s_i = 0; s_i < num_sockets_; ++s_i) {

        int partition_id = node_id_ * num_sockets_ + s_i;

        num_local_nodepair_per_partition_[s_i] = (uint32_t *) numa_alloc_onnode(
                        sizeof(uint32_t) * num_partitions_, s_i/ENGINE_PER_SOCKET);
        
        memset(num_local_nodepair_per_partition_[s_i], 0, sizeof(uint32_t) * num_partitions_);

        uint32_t i = 0, j = 0, k;
        while (i < n1->m_children && j < n2->m_children) {
            if (n1->m_ptrMBR[i]->m_pLow[0] < n2->m_ptrMBR[j]->m_pLow[0]) {
                RegionPtr region = n1->m_ptrMBR[i];
                k = j;
                while (k < n2->m_children && n2->m_ptrMBR[k]->m_pLow[0] <= region->m_pHigh[0]) {
                    if (region->m_pLow[1] <= n2->m_ptrMBR[k]->m_pHigh[1] && 
                        region->m_pHigh[1] >= n2->m_ptrMBR[k]->m_pLow[1]) {

                        int p_i = INIT_HASH(n1->m_pIdentifier[i], n2->m_pIdentifier[k]) % num_partitions_;
                        num_local_nodepair_per_partition_[s_i][p_i] ++;
                        
                    }
                    k++;
                }
                i++;
            }
            else {
                RegionPtr region = n2->m_ptrMBR[j];
                k = i;
                while (k < n1->m_children && n1->m_ptrMBR[k]->m_pLow[0] <= region->m_pHigh[0]) {
                    if (region->m_pLow[1] <= n1->m_ptrMBR[k]->m_pHigh[1] &&
                        region->m_pHigh[1] >= n1->m_ptrMBR[k]->m_pLow[1]) {

                        int p_i = INIT_HASH(n1->m_pIdentifier[k], n2->m_pIdentifier[j]) % num_partitions_;
                        num_local_nodepair_per_partition_[s_i][p_i] ++;
                    }
                    k++;
                }
                j++;
            }
        }

        local_nodepair_[s_i] = (id_type *) numa_alloc_onnode(
                                NODE_PAIR_SIZE * num_local_nodepair_per_partition_[s_i][partition_id], s_i/ENGINE_PER_SOCKET);

        uint32_t num_local_cnt = 0;
        i = 0;
        j = 0;
        while (i < n1->m_children && j < n2->m_children) {
            if (n1->m_ptrMBR[i]->m_pLow[0] < n2->m_ptrMBR[j]->m_pLow[0]) {
                RegionPtr region = n1->m_ptrMBR[i];
                k = j;
                while (k < n2->m_children && n2->m_ptrMBR[k]->m_pLow[0] <= region->m_pHigh[0]) {
                    if (region->m_pLow[1] <= n2->m_ptrMBR[k]->m_pHigh[1] && 
                        region->m_pHigh[1] >= n2->m_ptrMBR[k]->m_pLow[1]) {

                        int p_i = INIT_HASH(n1->m_pIdentifier[i], n2->m_pIdentifier[k]) % num_partitions_;
                        if (p_i == partition_id) {
                            local_nodepair_[s_i][num_local_cnt ++] = n1->m_pIdentifier[i];
                            local_nodepair_[s_i][num_local_cnt ++] = n2->m_pIdentifier[k];   
                        }
                    }
                    k++;
                }
                i++;
            }
            else {
                RegionPtr region = n2->m_ptrMBR[j];
                k = i;
                while (k < n1->m_children && n1->m_ptrMBR[k]->m_pLow[0] <= region->m_pHigh[0]) {
                    if (region->m_pLow[1] <= n1->m_ptrMBR[k]->m_pHigh[1] &&
                        region->m_pHigh[1] >= n1->m_ptrMBR[k]->m_pLow[1]) {
                        
                        int p_i = INIT_HASH(n1->m_pIdentifier[k], n2->m_pIdentifier[j]) % num_partitions_;
                        if (p_i == partition_id) {
                            local_nodepair_[s_i][num_local_cnt ++] = n1->m_pIdentifier[k];
                            local_nodepair_[s_i][num_local_cnt ++] = n2->m_pIdentifier[j];
                        }
                    }
                    k++;
                }
                j++;
            }
        }
        assert(num_local_cnt == num_local_nodepair_per_partition_[s_i][partition_id] * 2);
    }
}

void DistributedApplication::load_Rtrees()
{
    std::string baseName1 = rtree_name_1_;
    std::string baseName2 = rtree_name_2_;
    
    double utilization = 0.7;
    uint32_t capacity1 = 100;
    uint32_t capacity2 = 100;


    // =========  Build  Rrtree ========
    for (int s_i = 0; s_i <= num_sockets_; ++ s_i) {
        IStorageManager *diskfile1 = StorageManager::createNewDiskStorageManager(baseName1, 4096);
        // Create a new storage manager with the provided base name and a 4K page size.
        // StorageManager::IBuffer *file1 = StorageManager::createNewRandomEvictionsBuffer(*diskfile1, 10, false);
        ISpatialIndex *tree1 = RTree::loadRTree(*diskfile1, 1);
        Rrtree[s_i] = dynamic_cast<RTreeType *>(tree1);
    }

    // =========  Build Second Srtree ========
    for (int s_i = 0; s_i <= num_sockets_; ++ s_i) {
        IStorageManager *diskfile2 = StorageManager::createNewDiskStorageManager(baseName2, 4096);
        // StorageManager::IBuffer *file2 = StorageManager::createNewRandomEvictionsBuffer(*diskfile2, 10, false);
        ISpatialIndex *tree2 = RTree::loadRTree(*diskfile2, 1);
        Srtree[s_i] = dynamic_cast<RTreeType *>(tree2); 
    }
}


void DistributedApplication::run()
{
    std::thread * main_threads[MAX_NUM_SOCKETS];

    auto init_hash = [] (const id_type a, const id_type b) {
        return std::hash<id_type>{}(a) + std::hash<id_type>{}(b);
    };

    SpatialIndex::RTree::NodePtr n1 = Rrtree[0]->readNode(Rrtree[0]->m_rootID);
    SpatialIndex::RTree::NodePtr n2 = Srtree[0]->readNode(Srtree[0]->m_rootID);

    double ttl_wait_time[8];
    size_t total_data_read_r[8];
    size_t total_data_read_s[8];

    for (int s_i = 0; s_i < num_sockets_; ++ s_i) {
        engines_[s_i]->clear_all_buffers();
        
        main_threads[s_i] = new std::thread([&](int socket_id) {

            auto th_start_t = std::chrono::steady_clock::now();
            
            assert(numa_run_on_node(socket_id/ENGINE_PER_SOCKET) == 0);
            
            // cpu_set_t cpu_mask;
            // CPU_ZERO(&cpu_mask);
            // for (int i = 4; i < cpu_topo[socket_id/2].size(); ++i) {
            //     CPU_SET(i, &cpu_mask);
            // }
            // bind_to_core(cpu_mask);

            Workload workload;
            while (true) {
                workload_distributer_->fetch_next_batch(socket_id, workload);
                // workload_distributer_->fetch_next_batch_no_contention(socket_id, workload);

                if (workload.num == 0) { // no more workload
                    break;
                }
                for (uint32_t idx = 0; idx < workload.num; ++ idx) {
                    id_type id1 = workload.list[2 * idx];
                    id_type id2 = workload.list[2 * idx + 1];
                    engines_[socket_id]->scatter_vertex_extendable_embedding(id1, id2); 
                }
                engines_[socket_id]->flush_all_extendable_embeddings(); 
            }

            auto th_end_t = std::chrono::steady_clock::now();

            std::cout << socket_id << " takes " << 
                (float)std::chrono::duration_cast<std::chrono::milliseconds>(th_end_t - th_start_t).count() / 1000 
                << "s, wait time = " << (engines_[socket_id]->th_wait_time[0] + engines_[socket_id]->th_wait_time[1])/2 << std::endl;
            
            ttl_wait_time[socket_id] = (engines_[socket_id]->th_wait_time[0] + engines_[socket_id]->th_wait_time[1])/2;

            total_data_read_r[socket_id] = engines_[socket_id]->disk_reader1_->total_data_read;
            total_data_read_s[socket_id] = engines_[socket_id]->disk_reader2_->total_data_read;

        }, s_i);
    }



    for (int s_i = 0; s_i < num_sockets_; ++ s_i) {
        main_threads[s_i]->join();
        delete main_threads[s_i];
    }

    double avg = 0;
    for (int i=0; i<num_sockets_; ++i) {
        avg += ttl_wait_time[i];
    }

    double ttl_read_r = 0.0, ttl_read_s = 0.0;
    for (int i=0; i<num_sockets_; ++i) {
        ttl_read_r += (double) total_data_read_r[i] / 1024.0 / 1024.0 / 1024.0;
        ttl_read_s += (double) total_data_read_s[i] / 1024.0 / 1024.0 / 1024.0;
    }

    std::cout << "Average wait time = " << avg/num_sockets_ << std::endl;

    std::cout << "R read = " << ttl_read_r << std::endl;
    std::cout << "S read = " << ttl_read_s << std::endl;
}
