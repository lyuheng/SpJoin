#pragma once

#include <mpi.h>
#include <iostream>
#include <type_traits>

class DistributedSys {
private:
	static DistributedSys * instance_;

	int node_id_;
	int num_nodes_;

	DistributedSys();

public:
	static void init_distributed_sys();
	static void finalize_distributed_sys();
	static DistributedSys * get_instance();

	inline int get_node_id() {
	    return node_id_;
	}
	inline int get_num_nodes() {
	    return num_nodes_;
	}
    inline bool is_master_node() {
        return node_id_ == 0;
    }
	template<typename T>
    static MPI_Datatype get_mpi_data_type() {
        if (std::is_same<T, char>::value) {
            return MPI_CHAR;
        } else if (std::is_same<T, unsigned char>::value) {
            return MPI_UNSIGNED_CHAR;
        } else if (std::is_same<T, int>::value) {
            return MPI_INT;
        } else if (std::is_same<T, unsigned>::value) {
            return MPI_UNSIGNED;
        } else if (std::is_same<T, long>::value) {
            return MPI_LONG;
        } else if (std::is_same<T, unsigned long>::value) {
            return MPI_UNSIGNED_LONG;
        } else if (std::is_same<T, float>::value) {
            return MPI_FLOAT;
        } else if (std::is_same<T, double>::value) {
            return MPI_DOUBLE;
        } else {
            printf("type not supported\n");
            exit(-1);
        }
    }
};

inline DistributedSys * DistributedSys::instance_ = nullptr;

inline DistributedSys::DistributedSys() {
    int provided;
    MPI_Init_thread(NULL, NULL, MPI_THREAD_MULTIPLE, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &node_id_);
    MPI_Comm_size(MPI_COMM_WORLD, &num_nodes_);
    char host_name[128];
    int host_name_len;
    MPI_Get_processor_name(host_name, &host_name_len);
    host_name[host_name_len] = 0;
    printf("Nodename: %s\n", host_name);
}

inline void DistributedSys::init_distributed_sys() {
    assert(instance_ == nullptr);
    instance_ = new DistributedSys();
}

inline void DistributedSys::finalize_distributed_sys() {
    if (instance_ != nullptr) {
        MPI_Barrier(MPI_COMM_WORLD);
	    MPI_Finalize();
    }
}

inline DistributedSys * DistributedSys::get_instance() {
    if (instance_ == nullptr) {
	    init_distributed_sys();
    }
    return instance_;
}
