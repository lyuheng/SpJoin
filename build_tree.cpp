#include <spatialindex/SpatialIndex.h>

#include "src/rtree/RTree.h"

#include <chrono>
#include <cstring>

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
    // try
    {
        if (argc != 3)
        {
            std::cerr << "Usage: " << argv[0] << " input_file output_file" << std::endl;
            return -1;
        }

        std::string baseName1 = std::string(argv[2]);
        double utilization = 0.7;
        uint32_t capacity = 100;


        // =========  Build First R-tree ========
        IStorageManager *diskfile1 = StorageManager::createNewDiskStorageManager(baseName1, 4096);
        // Create a new storage manager with the provided base name and a 4K page size.

        StorageManager::IBuffer *file1 = StorageManager::createNewRandomEvictionsBuffer(*diskfile1, 10, false);
        // applies a main memory random buffer on top of the persistent storage manager
        // (LRU buffer, etc can be created the same way).

        MyDataStream stream1(argv[1]); 

        // Create and bulk load a new RTree with dimensionality 2, using "file" as
        // the StorageManager and the RSTAR splitting policy.
        id_type indexIdentifier1;
        ISpatialIndex *tree1 = RTree::createAndBulkLoadNewRTree(
            RTree::BLM_STR, stream1, *file1, utilization, capacity, capacity, 2, SpatialIndex::RTree::RV_RSTAR, indexIdentifier1);

        RTree::RTree *Rrtree = dynamic_cast<RTree::RTree *>(tree1);

        // =========  Build Second R-tree ========
        // IStorageManager *diskfile2 = StorageManager::createNewDiskStorageManager(baseName2, 4096);
        
        // // Create a new storage manager with the provided base name and a 4K page size.

        // StorageManager::IBuffer *file2 = StorageManager::createNewRandomEvictionsBuffer(*diskfile2, 10, false);
        // // applies a main memory random buffer on top of the persistent storage manager
        // // (LRU buffer, etc can be created the same way).

        // MyDataStream stream2(argv[2]);

        // // Create and bulk load a new RTree with dimensionality 2, using "file" as
        // // the StorageManager and the RSTAR splitting policy.
        // id_type indexIdentifier2;
        // ISpatialIndex *tree2 = RTree::createAndBulkLoadNewRTree(
        //     RTree::BLM_STR, stream2, *file2, utilization, capacity, capacity, 2, SpatialIndex::RTree::RV_RSTAR, indexIdentifier2);

        // RTree::RTree *Srtree = dynamic_cast<RTree::RTree *>(tree2);

        // ============= Initialize explict stack, etc ================
        // uint32_t Rtree_height = Rrtree->m_stats.getTreeHeight();
        // std::cout << "R tree_height = " << Rtree_height << " " << Rrtree->m_headerID << " " << indexIdentifier1 << '\n';

        // uint32_t Stree_height = Srtree->m_stats.getTreeHeight();
        // std::cout << "S tree_height = " << Stree_height << " " << Srtree->m_headerID << " " << indexIdentifier2 << '\n';   

        // * need reorder Tree's ID


        // uint32_t joinRes = 0;
        // double low[2] = {0, 0};
        // double high[2] = {100, 100};
        // Region r(low, high, 2);
        // Rrtree->SpatialJoinQuery1(Srtree, Rrtree->m_rootID, Srtree->m_rootID, r, joinRes);
        // std::cout << joinRes << "\n";


        delete tree1;
		delete file1;
		delete diskfile1;


        // delete tree2;
		// delete file2;
		// delete diskfile2;

        // delete everything     
    }
    // catch (Tools::Exception &e)
    // {
    //     std::cerr << "******ERROR******" << std::endl;
    //     std::string s = e.what();
    //     std::cerr << s << std::endl;
    //     return -1;
    // }

    return 0;
}
