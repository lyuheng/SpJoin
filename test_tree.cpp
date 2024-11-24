#include <spatialindex/SpatialIndex.h>

#include "src/rtree/RTree.h"
#include "src/rtree/Node.h"
#include "src/fine_grained/nodeptr.h"

#include <chrono>
#include <cstring>
#include <iostream>

using namespace SpatialIndex;
using namespace std::literals;

#define DELETE 0
#define INSERT 1
#define QUERY 2

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

int main(int argc, char **argv)
{
    try
    {
        std::string baseName1 = "tree1"s;
        std::string baseName2 = "tree2"s;
        double utilization = 0.7;
        uint32_t capacity = 100;

        // =========  Build First R-tree ========
        IStorageManager *diskfile1 = StorageManager::createNewDiskStorageManager(baseName1, 4096);
        // Create a new storage manager with the provided base name and a 4K page size.

        // StorageManager::IBuffer *file1 = StorageManager::createNewRandomEvictionsBuffer(*diskfile1, 10, false);
        // applies a main memory random buffer on top of the persistent storage manager
        // (LRU buffer, etc can be created the same way).

        // MyDataStream stream1(argv[1]);

        // Create and bulk load a new RTree with dimensionality 2, using "file" as
        // the StorageManager and the RSTAR splitting policy.
        // id_type indexIdentifier1;
        ISpatialIndex *tree1 = RTree::loadRTree(*diskfile1, 1);

        RTree::RTree *Rrtree = dynamic_cast<RTree::RTree *>(tree1);

        // =========  Build Second R-tree ========
        IStorageManager *diskfile2 = StorageManager::createNewDiskStorageManager(baseName2, 4096);
        
        // Create a new storage manager with the provided base name and a 4K page size.

        // StorageManager::IBuffer *file2 = StorageManager::createNewRandomEvictionsBuffer(*diskfile2, 10, false);
        // applies a main memory random buffer on top of the persistent storage manager
        // (LRU buffer, etc can be created the same way).

        ISpatialIndex *tree2 = RTree::loadRTree(*diskfile2, 1);

        RTree::RTree *Srtree = dynamic_cast<RTree::RTree *>(tree2);

        // ============= Initialize explict stack, etc ================
        uint32_t Rtree_height = Rrtree->m_stats.getTreeHeight();
        std::cout << "R tree_height = " << Rtree_height << " " << Rrtree->m_headerID << '\n';
        std::cout << "R tree nodes = " << Rrtree->m_stats.getNumberOfNodes() << '\n';

        Rrtree->readNode(Rrtree->m_rootID);
        std::cout << Rrtree->readNodeLength(Rrtree->m_rootID) << std::endl;

        uint32_t Stree_height = Srtree->m_stats.getTreeHeight();
        std::cout << "S tree_height = " << Stree_height << " " << Srtree->m_headerID << '\n'; 
        std::cout << "S tree nodes = " << Srtree->m_stats.getNumberOfNodes() << '\n';

        Srtree->readNode(Srtree->m_rootID);
        std::cout << Srtree->readNodeLength(Srtree->m_rootID) << std::endl;


        uint32_t joinRes = 0;
        double low[2] = {0, 0};
        double high[2] = {100, 100};
        Region r(low, high, 2);
        Rrtree->SpatialJoinQuery1(Srtree, Rrtree->m_rootID, Srtree->m_rootID, r, joinRes);
        std::cout << joinRes << "\n";

        RTree::NodePtr n1 = Rrtree->readNode(20);
        // std::cout << n1->m_children << std::endl;

        auto len = Rrtree->readNodeLength(20);
        uint8_t *arr = new uint8_t[len];
        Rrtree->readRawData(20, arr);
        SimpleNodePtr sn;
        sn.fill_data(20, arr);

        // std::cout << sn.returnChildRegion(5).low0 << " " << n1->m_ptrMBR[5]->m_pLow[0] << std::endl;
        // std::cout << sn.returnChildRegion(5).low[1] << " " << n1->m_ptrMBR[5]->m_pLow[1] << std::endl;

        // std::cout << sn.returnChildRegion(5).high[0] << " " << n1->m_ptrMBR[5]->m_pHigh[0] << std::endl;
        // std::cout << sn.returnChildRegion(5).high[1] << " " << n1->m_ptrMBR[5]->m_pHigh[1] << std::endl;

        // std::cout << sn.returnChildIdentifier(5) << " " << n1->m_pIdentifier[5] << std::endl;


        

        // RTree::NodePtr n2 = Srtree->readNode(25);
        // std::cout << n2->m_children << std::endl;


        delete tree1;
		// delete file1;
		delete diskfile1;


        delete tree2;
		// delete file2;
		delete diskfile2;

        // delete everything     
    }
    catch (Tools::Exception &e)
    {
        std::cerr << "******ERROR******" << std::endl;
        std::string s = e.what();
        std::cerr << s << std::endl;
        return -1;
    }

    return 0;
}
