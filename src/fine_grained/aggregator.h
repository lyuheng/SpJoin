#pragma once

#include "meta.h"
#include "mpi.h"
#include "distsys.h"

template<typename T>
class Aggregator {
private:
    const int padding = 128;
    int num_threads_;
    int num_sockets_;
    T * local_reducers_[MAX_NUM_SOCKETS][MAX_NUM_THREADS];
public:
    inline void clear() {
        for (int s_i = 0; s_i < num_sockets_; ++ s_i) {
            for (int t_i = 0; t_i < num_threads_; ++ t_i) {
                memset(local_reducers_[s_i][t_i], 0, sizeof(T));
            }
        }
    }
    inline T evaluate() {
        T reducer = 0;
        for (int s_i = 0; s_i < num_sockets_; ++ s_i) {
            for (int t_i = 0; t_i < num_threads_; ++ t_i) {
                reducer += *local_reducers_[s_i][t_i];
            }
        }

        std::cout << "############# reducer = " << reducer << std::endl;

        T global_reducer = 0;
        MPI_Allreduce(
                &reducer, &global_reducer, 1,
                DistributedSys::get_mpi_data_type<T>(),
                MPI_SUM, MPI_COMM_WORLD
                );
        return global_reducer;
    }
    inline void add(int s_i, int t_i, T delta) {
        *(local_reducers_[s_i][t_i]) += delta;
    }
    inline void list_all_result() {
        for (int i = 0; i < num_sockets_; ++i) {
            for (int j = 0; j < num_threads_; ++j) {
                std::cout << *local_reducers_[i][j] << " ";
            }
        }
        std::cout << std::endl;
    }

    Aggregator(int num_sockets, int num_threads) {
        num_sockets_ = num_sockets;
        num_threads_ = num_threads;
        std::cout << "              num_sockets = " << num_sockets << ", num_threads_ = " << num_threads_ << std::endl; 
        for (int s_i = 0; s_i < num_sockets_; ++ s_i) {
            for (int t_i = 0; t_i < num_threads_; ++ t_i) {
                local_reducers_[s_i][t_i] = (T*) numa_alloc_onnode(
                        sizeof(T) * (1 + padding / sizeof(T)), s_i/ENGINE_PER_SOCKET
                        );
            }
        }
        clear();
    }
    ~Aggregator() {
        for (int s_i = 0; s_i < num_sockets_; ++ s_i) {
            for (int t_i = 0; t_i < num_threads_; ++ t_i) {
                numa_free(local_reducers_[s_i][t_i], 
                        sizeof(T) * (1 + padding / sizeof(T)));
            }
        }
    }
};