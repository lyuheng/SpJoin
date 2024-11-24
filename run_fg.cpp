#include "src/fine_grained/distributed.h"
#include "src/fine_grained/distsys.h"

#include <iostream>
#include <chrono>


int main(int argc, char *argv[])
{
    DistributedSys::init_distributed_sys();

    int node_id = DistributedSys::get_instance()->get_node_id();
    std::cout << "My node_id = " << node_id << std::endl;

    std::string rtree_name_1;
    std::string rtree_name_2;
    if (node_id == 0) {
        rtree_name_1 = std::string(argv[1]);
        rtree_name_2 = std::string(argv[2]);
    } else {
        rtree_name_1 = std::string(argv[3]);
        rtree_name_2 = std::string(argv[4]);
    }

    {
        DistributedApplication app(rtree_name_1, rtree_name_2, -1); // thread_num = -1
        

        auto start_time = std::chrono::steady_clock::now();

        app.count->clear();
        app.run();
        auto end_time = std::chrono::steady_clock::now();
        std::cout << "Elasped Time: " << (float)std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count() / 1000 << "s\n";
        
        std::cout << "Result = " << app.count->evaluate() << std::endl;
        // app.count->list_all_result();
    }

    DistributedSys::finalize_distributed_sys();
    
    return 0;
}