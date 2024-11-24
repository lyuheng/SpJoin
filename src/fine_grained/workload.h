#pragma once

#include "meta.h"
#include <mpi.h>

#include <numa.h>

struct Workload {
    id_type * list;
    uint32_t num;
} __attribute__((packed));

class WorkloadDistributer {
public:

    constexpr static size_t NODE_PAIR_SIZE  = sizeof(id_type) * 2;

    uint32_t batch_size = MAX_WORKLOAD_BATCH_SIZE;
    id_type * stolen_nodepair_buff_[MAX_NUM_SOCKETS]; // one per socket

    int node_id_;
    int num_nodes_;
    int num_sockets_;
    int num_partitions_;

    // the progress of all execution engines on this node
    uint32_t curr_progress_[MAX_NUM_SOCKETS]; 
    MPI_Win progress_win_;
    std::mutex progress_win_mutex_;

    // tracking the node && socket being stolen
    int remote_partition_id_[MAX_NUM_SOCKETS];
    uint32_t num_local_nodepair_[MAX_NUM_SOCKETS];
    uint32_t num_local_nodepair_per_partition_[MAX_NUM_PARTITIONS];

    id_type * local_nodepair_[MAX_NUM_SOCKETS];
    MPI_Win local_nodepair_win_[MAX_NUM_SOCKETS];

    void fetch_next_batch_local(int s_i, Workload &workload);
    void fetch_next_batch_remote(int s_i, Workload &workload);
    void next_remote_partition(int s_i);
    void fetch_next_batch(int s_i, Workload &workload);
    void fetch_next_batch_no_contention(int s_i, Workload &workload); // * used when no work stealing is applied

    WorkloadDistributer(int node_id, 
                        int num_nodes, 
                        int num_sockets,
                        uint32_t ** num_local_nodepair_per_partition,
                        id_type ** local_nodepair
                    );
    ~WorkloadDistributer();
};

void WorkloadDistributer::fetch_next_batch_local(int s_i, Workload &workload)
{
    uint32_t boundary = num_local_nodepair_[s_i];
    uint32_t pos;
    uint32_t delta = batch_size;

    progress_win_mutex_.lock();
    int ret = MPI_Fetch_and_op(&delta, &pos, MPI_UNSIGNED,
                                node_id_, s_i, MPI_SUM, progress_win_);
    assert(ret == MPI_SUCCESS);
    MPI_Win_flush(node_id_, progress_win_);
    progress_win_mutex_.unlock();
    
    if (pos >= boundary) {
        workload.num = 0;
        workload.list = nullptr;
    }
    else {
        uint32_t begin = pos;
        uint32_t end = std::min(pos + batch_size, boundary);
        workload.num = end - begin;
        workload.list = local_nodepair_[s_i] + 2 * begin;
    }
}

void WorkloadDistributer::fetch_next_batch_remote(int s_i, Workload &workload)
{
    workload.num = 0;
    workload.list = nullptr;

    int partition_id = node_id_ * num_sockets_ + s_i;

    while (remote_partition_id_[s_i] != partition_id) {
        int remote_p_i = remote_partition_id_[s_i];
        int remote_node_id = remote_p_i / num_sockets_;
        int remote_socket_id = remote_p_i % num_sockets_;

        uint32_t boundary = num_local_nodepair_per_partition_[remote_p_i]; 
        uint32_t pos;
        uint32_t delta = batch_size;

        progress_win_mutex_.lock();
        int ret = MPI_Fetch_and_op(
            &delta, &pos, MPI_UNSIGNED,
            remote_node_id, remote_socket_id, MPI_SUM, progress_win_);
        assert(ret == MPI_SUCCESS);
        MPI_Win_flush(remote_node_id, progress_win_);
        progress_win_mutex_.unlock();

        if (pos < boundary) {
            uint32_t begin = pos;
            uint32_t end = std::min(pos + batch_size, boundary);
            workload.num = end - begin;
            workload.list = stolen_nodepair_buff_[s_i];

            MPI_Get(
                stolen_nodepair_buff_[s_i], 2 * (end - begin),
                MPI_LONG,
                remote_node_id,
                2 * begin, 2 * (end - begin),
                MPI_LONG,
                local_nodepair_win_[remote_socket_id]);
            MPI_Win_flush(remote_node_id, local_nodepair_win_[remote_socket_id]);
            break;
        }
        else {
            next_remote_partition(s_i);
        }
    }
}

void WorkloadDistributer::next_remote_partition(int s_i)
{
    remote_partition_id_[s_i] += 1;
    remote_partition_id_[s_i] %= num_partitions_;
}

void WorkloadDistributer::fetch_next_batch(int s_i, Workload &workload)
{
    fetch_next_batch_local(s_i, workload);
    if (workload.num != 0) return; 
    fetch_next_batch_remote(s_i, workload);
}

void WorkloadDistributer::fetch_next_batch_no_contention(int s_i, Workload &workload)
{
    uint32_t boundary = num_local_nodepair_[s_i];
    uint32_t pos;
    uint32_t delta = batch_size;

    pos = curr_progress_[s_i];
    curr_progress_[s_i] += delta;
    
    if (pos >= boundary) {
        workload.num = 0;
        workload.list = nullptr;
    }
    else {
        uint32_t begin = pos;
        uint32_t end = std::min(pos + batch_size, boundary);
        workload.num = end - begin;
        workload.list = local_nodepair_[s_i] + 2 * begin;
    }
}


WorkloadDistributer::WorkloadDistributer(int node_id, 
                                        int num_nodes, 
                                        int num_sockets,
                                        uint32_t ** num_local_nodepair_per_partition,
                                        id_type ** local_nodepair
                                    )
{
    node_id_ = node_id;
    num_nodes_ = num_nodes;
    num_sockets_ = num_sockets;
    num_partitions_ = num_nodes_ * num_sockets_;

    for (int s_i = 0; s_i < num_sockets_; ++s_i) {
        stolen_nodepair_buff_[s_i] = (id_type *)numa_alloc_onnode(
            NODE_PAIR_SIZE * MAX_WORKLOAD_BATCH_SIZE, s_i/ENGINE_PER_SOCKET);
        memset(stolen_nodepair_buff_[s_i], 0, NODE_PAIR_SIZE * MAX_WORKLOAD_BATCH_SIZE);
    }

    for (int s_i = 0; s_i < num_sockets_; ++s_i) {
        curr_progress_[s_i] = 0;
        remote_partition_id_[s_i] = node_id_ * num_sockets_ + s_i;
        next_remote_partition(s_i);
    }
    MPI_Win_create(curr_progress_, sizeof(uint32_t) * num_sockets_,
                    sizeof(uint32_t), MPI_INFO_NULL, MPI_COMM_WORLD, &progress_win_);
    // passive synchronization
    for (int n_i = 0; n_i < num_nodes_; ++n_i) {
        MPI_Win_lock(MPI_LOCK_SHARED, n_i, 0, progress_win_);
    }
    for (int s_i = 0; s_i < num_sockets_; ++s_i) {
        int p_i = node_id_ * num_sockets_ + s_i;
        num_local_nodepair_[s_i] = num_local_nodepair_per_partition[s_i][p_i];
    }

    MPI_Allgather(
        num_local_nodepair_, num_sockets_,
        MPI_UNSIGNED,
        num_local_nodepair_per_partition_, num_sockets_,
        MPI_UNSIGNED, MPI_COMM_WORLD
    );
    for (int s_i = 0; s_i < num_sockets_; ++s_i) {
        assert(num_local_nodepair_[s_i] == num_local_nodepair_per_partition_[node_id_ * num_sockets_ + s_i]);
    }
    for (int s_i = 0; s_i < num_sockets_; ++s_i) {
        local_nodepair_[s_i] = local_nodepair[s_i];
        MPI_Win_create(
            local_nodepair_[s_i], NODE_PAIR_SIZE * num_local_nodepair_[s_i],
            sizeof(id_type), MPI_INFO_NULL, MPI_COMM_WORLD, &local_nodepair_win_[s_i]);

        for (int n_i = 0; n_i < num_nodes_; ++n_i) {
            MPI_Win_lock(MPI_LOCK_SHARED, n_i, 0, local_nodepair_win_[s_i]);
        }
    }
}


WorkloadDistributer::~WorkloadDistributer()
{
    for (int n_i = 0; n_i < num_nodes_; ++n_i) {
        MPI_Win_unlock(n_i, progress_win_);
    }
    MPI_Win_free(&progress_win_);

    for (int s_i = 0; s_i < num_sockets_; ++s_i) {
        for (int n_i = 0; n_i < num_nodes_; ++n_i) {
            MPI_Win_unlock(n_i, local_nodepair_win_[s_i]);
        }
        MPI_Win_free(&local_nodepair_win_[s_i]);
    }

    for (int s_i = 0; s_i < num_sockets_; ++s_i) {
        numa_free(stolen_nodepair_buff_[s_i], NODE_PAIR_SIZE * MAX_WORKLOAD_BATCH_SIZE);
    }
}