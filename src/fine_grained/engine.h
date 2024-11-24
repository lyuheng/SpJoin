#pragma once

#include <spatialindex/SpatialIndex.h>
#include <src/rtree/Node.h>
#include <src/rtree/RTree.h>

#include <src/rtree/Leaf.h>
#include <src/rtree/Index.h>
#include <src/rtree/BulkLoader.h>


#include <numa.h>

#include <thread>
#include <mutex>
#include <condition_variable>

#include "aggregator.h"
#include "nodeptr.h"


using namespace SpatialIndex;
typedef RTree::RTree RTreeType;

typedef SpatialIndex::id_type id_type;
typedef SpatialIndex::RTree::NodePtr NodePtr;


class EmbeddingExplorationEngine;

struct CompactExtendableEmbedding {
    id_type id1;                        // 8 byte
    id_type id2;                        // 8 byte
    
    uint8_t * new_disk_data_1;          // 8 byte
    uint8_t * new_disk_data_2;          // 8 byte
    lli disk_data_length_1;             // 8 byte
    lli disk_data_length_2;             // 8 byte

    // reserved space for unequal tree heights scenario
    double low0;
    double low1;
    double high0;
    double high1;
} __attribute__((packed));

enum GlobalEmbeddingQueueState {
    InValid = 1,
    PartialReady = 2
};

struct NumaAwareBuffer {
    lli buffer_size;
    uint8_t ** buffer_pointer;
};

class NumaAwareBufferManager {
private:
    std::vector<NumaAwareBuffer> ** registered_buffers; // [num_sockets_]
    lli * total_registered_buffer_size; // [num_sockets_]
    int num_sockets_;
    void ** allocated_ptx_;

public:
    NumaAwareBufferManager(int num_sockets)
        :num_sockets_(num_sockets)
    {
        total_registered_buffer_size = new lli [num_sockets];
        registered_buffers = new std::vector<NumaAwareBuffer>* [num_sockets];
        allocated_ptx_ = new void* [num_sockets];
        for (int s_i = 0; s_i < num_sockets; ++ s_i) {
            total_registered_buffer_size[s_i] = 0;
            registered_buffers[s_i] = new std::vector<NumaAwareBuffer>();
            registered_buffers[s_i]->clear();
        }
    }
    
    ~NumaAwareBufferManager()
    {
        delete [] total_registered_buffer_size;
        for (int s_i = 0; s_i < num_sockets_; ++ s_i) {
            delete registered_buffers[s_i];
        }
        delete [] registered_buffers;
        delete [] allocated_ptx_;
    }

    void register_numa_aware_buffer(lli buffer_size, uint8_t ** buffer_pointer, int socket_id)
    {
        assert(socket_id >= 0 && socket_id < num_sockets_);
        total_registered_buffer_size[socket_id] += buffer_size;
        NumaAwareBuffer buffer;
        buffer.buffer_size = buffer_size;
        buffer.buffer_pointer = buffer_pointer;
        registered_buffers[socket_id]->push_back(buffer);
    }
    void allocate_numa_aware_buffers()
    {
        for (int s_i = 0; s_i < num_sockets_; ++ s_i) {
            if (total_registered_buffer_size[s_i] == 0) {
                continue;
            }

            uint8_t * ptx = (uint8_t*) numa_alloc_onnode(total_registered_buffer_size[s_i], s_i/ENGINE_PER_SOCKET);
            allocated_ptx_[s_i] = ptx;

            if (ptx == NULL) {
                fprintf(stderr, "Failed to allocate %lld bytes memory on socket %d.\n", 
                        total_registered_buffer_size[s_i], s_i);
                exit(-1);
            }

            // *@@@@@@@@@@@@@@@@@@@@@@@
            // memset(ptx, 0, total_registered_buffer_size[s_i]);
            int num_buffers = registered_buffers[s_i]->size();
            lli size_sum = 0;
            for (int i = 0; i < num_buffers; ++ i) {
                NumaAwareBuffer buffer = registered_buffers[s_i]->at(i);
                *buffer.buffer_pointer = ptx;
                ptx += buffer.buffer_size;
                size_sum += buffer.buffer_size;
            }
            assert(size_sum == total_registered_buffer_size[s_i]);
        }
    }
    void deallocate_numa_aware_buffers() // must be called after allocate_numa_aware_buffers()
    {
        for (int s_i = 0; s_i < num_sockets_; ++ s_i) {
            if (total_registered_buffer_size[s_i] > 0) {
                numa_free(allocated_ptx_[s_i], total_registered_buffer_size[s_i]);
            }
        }
    }
    lli get_total_allocated_size()
    {
        lli total_buffer_size = 0;
        for (int s_i = 0; s_i < num_sockets_; ++ s_i) {
            total_buffer_size += total_registered_buffer_size[s_i];
        }
        return total_buffer_size;
    }

    template<typename T>
        T* alloc_numa_aware_array(int array_size, int socket_id) {
            void * ptx = numa_alloc_onnode(sizeof(T) * array_size, socket_id/ENGINE_PER_SOCKET);
            assert(ptx != NULL);
            return (T*) ptx;
        }
};

struct DiskRef {
    uint8_t * addr;
    id_type node_id;
}__attribute__((packed));

class DiskDataReader_R
{
public:
    const lli disk_data_size_per_request = DATA_SIZE_PER_REQUEST;
    const lli num_nodes_per_request = MAX_NUM_NODES_PER_REQUEST;

    EmbeddingExplorationEngine * engine_;

    int socket_id_;
    int num_sockets_;
    std::thread * thread_;
    volatile bool is_terminated_;

    const int num_confused_bits = 0;
    uint64_t hash_set_size_;
    lli node_id_mask_;

    NumaAwareBufferManager * buff_manager_;
    uint32_t curr_generation_[MAX_DEPTH]; // using the generation method to avoid bitmap resetting cost
    uint32_t * previous_request_generation_[MAX_DEPTH]; 
    id_type * previous_request_node_id_[MAX_DEPTH]; 
    uint8_t ** previous_request_disk_buffer_[MAX_DEPTH];

    // * for alternative queues
    uint32_t * previous_request_generation_alternative_[MAX_DEPTH]; 
    id_type * previous_request_node_id_alternative_[MAX_DEPTH]; 
    uint8_t ** previous_request_disk_buffer_alternative_[MAX_DEPTH];

    DiskRef * requested_node_data_address[MAX_DEPTH];


    size_t total_data_read = 0;
    
    DiskDataReader_R(EmbeddingExplorationEngine * engine);
    ~DiskDataReader_R();

    int get_next_request_batch( lli & num_requested_nodes,
                                lli & requested_disk_data_size,
                                lli & requested_memcpy_size,
                                lli & curr_embedding_idx
                            );
    void thread_main();
};


//=================================================================================
//=================================================================================
//=================================================================================
//=================================================================================
//=================================================================================

class DiskDataReader_S
{
public:
    const lli disk_data_size_per_request = DATA_SIZE_PER_REQUEST;
    const lli num_nodes_per_request = MAX_NUM_NODES_PER_REQUEST;

    EmbeddingExplorationEngine * engine_;

    int socket_id_;
    int num_sockets_;
    std::thread * thread_;
    volatile bool is_terminated_;

    const int num_confused_bits = 2;
    uint64_t hash_set_size_;
    lli node_id_mask_;

    NumaAwareBufferManager * buff_manager_;
    uint32_t curr_generation_[MAX_DEPTH]; // using the generation method to avoid bitmap resetting cost
    uint32_t * previous_request_generation_[MAX_DEPTH]; 
    id_type * previous_request_node_id_[MAX_DEPTH]; 
    uint8_t ** previous_request_disk_buffer_[MAX_DEPTH];

    // * for alternative queues
    uint32_t * previous_request_generation_alternative_[MAX_DEPTH]; 
    id_type * previous_request_node_id_alternative_[MAX_DEPTH]; 
    uint8_t ** previous_request_disk_buffer_alternative_[MAX_DEPTH];

    DiskRef * requested_node_data_address[MAX_DEPTH];


    size_t total_data_read = 0;
    
    DiskDataReader_S(EmbeddingExplorationEngine * engine);
    ~DiskDataReader_S();

    int get_next_request_batch( lli & num_requested_nodes,
                                lli & requested_disk_data_size,
                                lli & requested_memcpy_size,
                                lli & curr_embedding_idx
                            );
    void thread_main();
};


//=================================================================================
//=================================================================================
//=================================================================================
//=================================================================================
//=================================================================================

class EmbeddingExplorationEngine {

public:
    EmbeddingExplorationEngine ** engines_same_node_ = nullptr;

    RTreeType * Rrtree_;
    RTreeType * Srtree_;

    Aggregator<size_t> * count_ = nullptr;

    int num_sockets_;
    int socket_id_;
    volatile int next_depth_;

    int num_partitions_;
    int partition_id_;

    int chunk_size_;

    DiskDataReader_R * disk_reader1_;
    DiskDataReader_S * disk_reader2_;

    int num_computation_threads_;
    int max_depth_; // depth = 0...max_depth
    int max_depth_R_; 
    int max_depth_S_;

    lli max_num_embeddings_global_queue_;
    lli global_disk_data_buffer_size_;

    GlobalEmbeddingQueueState global_queue_states_[MAX_DEPTH];

    NumaAwareBufferManager * buff_manager_;

    CompactExtendableEmbedding * global_embedding_queues_[MAX_DEPTH];
    // * main queues
    uint8_t * global_disk_data_buffer1_[MAX_DEPTH];
    uint8_t * global_disk_data_buffer2_[MAX_DEPTH];
    // * alternative queues
    uint8_t * global_disk_data_buffer1_alternative_[MAX_DEPTH];
    uint8_t * global_disk_data_buffer2_alternative_[MAX_DEPTH];

    lli global_embedding_queue_size_[MAX_DEPTH];

    lli global_num_ready_embeddings1_[MAX_DEPTH];
    lli global_num_ready_embeddings2_[MAX_DEPTH];

    lli global_used_disk_data_buffer_size1_[MAX_DEPTH];
    lli global_used_disk_data_buffer_size2_[MAX_DEPTH];


    std::mutex global_embedding_queue_mutex_;

    // thread-local data structures
    lli max_num_embeddings_local_queue_; // this paramter should be carefully choosen
    lli local_disk_data_buffer_size_;

    // padding to avoid L1/L2-cache-line false sharing
    CompactExtendableEmbedding * local_embedding_queues_[MAX_NUM_THREADS][MAX_DEPTH + 128];
    lli local_embedding_queue_size_[MAX_NUM_THREADS][MAX_DEPTH + 128];
    lli local_needed_disk_data_buffer_size1_[MAX_NUM_THREADS][MAX_DEPTH + 128];
    lli local_needed_disk_data_buffer_size2_[MAX_NUM_THREADS][MAX_DEPTH + 128];

    // CV and mutex used to managing computation threads
    std::mutex num_suspended_comp_threads_mutex_[MAX_DEPTH];
    std::condition_variable num_suspended_comp_threads_cv_[MAX_DEPTH];
    volatile int num_suspended_threads_[MAX_DEPTH];
    volatile int global_phases_[MAX_DEPTH];
    // there could be some false sharing concerning this array => however it is fine since the operations updating this array are rare
    int local_phases_[MAX_DEPTH][MAX_NUM_THREADS][128];
    // we do not use std::barrier since it is only supported after C++20
    pthread_barrier_t comp_thread_barrier_[MAX_DEPTH];

    EmbeddingExplorationEngine(
        lli max_depth,
        lli max_depth_R,
        lli max_depth_S,
        lli max_num_embeddings_global_queue,
        lli global_disk_data_buffer_size,
        int num_partitions,
        int socket_id,
        int num_sockets,
        int num_computation_threads,           /* computation threads per socket */
        int chunk_size,
        RTreeType * Rrtree,
        RTreeType * Srtree,
        Aggregator<size_t> * count
    );
    ~EmbeddingExplorationEngine();

    void flush_all_extendable_embeddings();
    void change_global_queue_state_to_partial_ready(int depth);
    void extend_embeddings(int depth);
    void shuffle_global_embedding_queue(int depth);
    void shuffle_global_embedding_alternative_queue(int depth);
    int flush_local_embeddings(int thread_id, int depth);
    void clear_global_queue(int depth);
    void scatter(id_type id1, id_type id2, int thread_id_, int depth);
    void scatter_imbalance(id_type id1, id_type id2, int thread_id_, int depth, SimpleRegion const& region);
    void suspend_thread(const int thread_id, int depth);

    void init_buffers();
    void release_buffer();
    void clear_all_buffers();
    void scatter_vertex_extendable_embedding(id_type id1, id_type id2);

    void set_engines_on_the_same_node(EmbeddingExplorationEngine ** engines) {
        engines_same_node_ = engines;
    }

    // spatial-join functions
    void spatial_join_func( id_type id1, 
                            id_type id2, 
                            NodePtr & n1, 
                            NodePtr & n2,
                            int thread_id,
                            int depth
                        );
    void spatial_join_func_simplenode(  id_type id1, 
                                        id_type id2, 
                                        SimpleNodePtr * n1, 
                                        SimpleNodePtr * n2,
                                        int thread_id,
                                        int depth);

    void spatial_join_func_simplenode_oneside( CompactExtendableEmbedding * e,
                                                SimpleNodePtr * n,
                                                int thread_id,
                                                int depth);


    double th_wait_time[2] {0.0, 0.0}; 

};