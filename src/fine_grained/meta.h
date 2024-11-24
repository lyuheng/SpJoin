#pragma once

#include <spatialindex/SpatialIndex.h>
#include <vector>
#include <sched.h>
#include <iostream>

#include <sys/time.h>

// * Hardward related
#define MAX_NUM_SOCKETS 8
#define MAX_NUM_THREADS 64 
#define MAX_NUM_NODES 16 // at most 16 machines/nodes used in experiments

#define MAX_NUM_PARTITIONS (MAX_NUM_SOCKETS * MAX_NUM_NODES)

inline const std::vector<std::vector<int>> cpu_topo = { /*socket 0*/ {0,1,2,3,4,5,6,7,           /*hyper-thread*/ 32,33,34,35,36,37,38,39},
                                                        /*socket 1*/ {8,9,10,11,12,13,14,15,     /*hyper-thread*/ 40,41,42,43,44,45,46,47},
                                                        /*socket 2*/ {16,17,18,19,20,21,22,23,   /*hyper-thread*/ 48,49,50,51,52,53,54,55},
                                                        /*socket 3*/ {24,25,26,27,28,29,30,31,   /*hyper-thread*/ 56,57,58,59,60,61,62,63} };

// * Cache size related (on Polaris, ALCF for example)
#define REAL_SOCKET_NODE 4
#define L1_DCACHE_SIZE (1 * 1024 * 1024)    // 1MB
// #define L2_CACHE_SIZE (16 * 1024 * 1024)    // 16MB 
// #define LLC_CACHE_SIZE (256 * 1024 * 1024)  // 256MB

#define ENGINE_PER_SOCKET 2 // specify how many engines executed on a socket node


// * typedef related
typedef long long int lli;


// * R-tree related
#define MAX_DEPTH 6     // maximum height for R-tree


// * Workload control related
#define CHUNK_SIZE 4
#define DATA_SIZE_PER_REQUEST (4 * 1024 * 1024)  // 4 MB
#define MAX_NUM_NODES_PER_REQUEST 100 // read 100 R-tree nodes once at most 

// * Workload stealing related
#define MAX_WORKLOAD_BATCH_SIZE 1


inline void bind_to_core(size_t core)
{
    cpu_set_t mask;
    CPU_ZERO(&mask);
    CPU_SET(core, &mask);
    if (sched_setaffinity(0, sizeof(mask), &mask) != 0)
        std::cout << "Failed to set affinity (core: " << core << ")" << std::endl;
}

inline void bind_to_core(cpu_set_t mask)
{
    if (sched_setaffinity(0, sizeof(mask), &mask) != 0)
        std::cout << "Fail to set affinity!" << std::endl;
}

// socket id can be greater than real # socket CPU node 
inline int which_core_to_use(int socket_id, int thread_id)
{
    if (socket_id >= REAL_SOCKET_NODE * 2) {
        throw std::runtime_error("Not implemented!");
    }

    if (ENGINE_PER_SOCKET == 2) {
        if (socket_id < REAL_SOCKET_NODE) {
            switch (thread_id)
            {
            case 0:
                return 4;
                break;
            case 1: 
                return 5;
                break;
            default:
                throw std::runtime_error("Not implemented!");
            }
        }
        else {
            switch (thread_id)
            {
            case 0:
                return 6;
                break;
            case 1: 
                return 7;
                break;
            default:
                throw std::runtime_error("Not implemented!");
            }
        }
    }
    else if (ENGINE_PER_SOCKET == 1) {
        if (thread_id < 6) return 2 + thread_id;
        else throw std::runtime_error("Not implemented!");
    }

    throw std::runtime_error("Not implemented!");

    return -1;
}


inline double get_time() {
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return tv.tv_sec + (tv.tv_usec / 1e6);
}