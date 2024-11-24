/******************************************************************************
 * Project:  libspatialindex - A C++ library for spatial indexing
 * Author:   Marios Hadjieleftheriou, mhadji@gmail.com
 ******************************************************************************
 * Copyright (c) 2002, Marios Hadjieleftheriou
 *
 * All rights reserved.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
******************************************************************************/

#include <cstring>
#include <cmath>
#include <limits>
#include <algorithm> // std::sort

#include <spatialindex/SpatialIndex.h>
#include <spatialindex/capi/IdVisitor.h>
#include "Node.h"
#include "Leaf.h"
#include "Index.h"
#include "BulkLoader.h"
#include "RTree.h"

using namespace SpatialIndex::RTree;
using namespace SpatialIndex;

SpatialIndex::RTree::Data::Data(uint32_t len, uint8_t* pData, Region& r, id_type id)
	: m_id(id), m_region(r), m_pData(nullptr), m_dataLength(len)
{
	if (m_dataLength > 0)
	{
		m_pData = new uint8_t[m_dataLength];
		memcpy(m_pData, pData, m_dataLength);
	}
}

SpatialIndex::RTree::Data::~Data()
{
	delete[] m_pData;
}

SpatialIndex::RTree::Data* SpatialIndex::RTree::Data::clone()
{
	return new Data(m_dataLength, m_pData, m_region, m_id);
}

id_type SpatialIndex::RTree::Data::getIdentifier() const
{
	return m_id;
}

void SpatialIndex::RTree::Data::getShape(IShape** out) const
{
	*out = new Region(m_region);
}

void SpatialIndex::RTree::Data::getData(uint32_t& len, uint8_t** data) const
{
	len = m_dataLength;
	*data = nullptr;

	if (m_dataLength > 0)
	{
		*data = new uint8_t[m_dataLength];
		memcpy(*data, m_pData, m_dataLength);
	}
}

uint32_t SpatialIndex::RTree::Data::getByteArraySize()
{
	return
		sizeof(id_type) +
		sizeof(uint32_t) +
		m_dataLength +
		m_region.getByteArraySize();
}

void SpatialIndex::RTree::Data::loadFromByteArray(const uint8_t* ptr)
{
	memcpy(&m_id, ptr, sizeof(id_type));
	ptr += sizeof(id_type);

	delete[] m_pData;
	m_pData = nullptr;

	memcpy(&m_dataLength, ptr, sizeof(uint32_t));
	ptr += sizeof(uint32_t);

	if (m_dataLength > 0)
	{
		m_pData = new uint8_t[m_dataLength];
		memcpy(m_pData, ptr, m_dataLength);
		ptr += m_dataLength;
	}

	m_region.loadFromByteArray(ptr);
}

void SpatialIndex::RTree::Data::storeToByteArray(uint8_t** data, uint32_t& len)
{
	// it is thread safe this way.
	uint32_t regionsize;
	uint8_t* regiondata = nullptr;
	m_region.storeToByteArray(&regiondata, regionsize);

	len = sizeof(id_type) + sizeof(uint32_t) + m_dataLength + regionsize;

	*data = new uint8_t[len];
	uint8_t* ptr = *data;

	memcpy(ptr, &m_id, sizeof(id_type));
	ptr += sizeof(id_type);
	memcpy(ptr, &m_dataLength, sizeof(uint32_t));
	ptr += sizeof(uint32_t);

	if (m_dataLength > 0)
	{
		memcpy(ptr, m_pData, m_dataLength);
		ptr += m_dataLength;
	}

	memcpy(ptr, regiondata, regionsize);
	delete[] regiondata;
	// ptr += regionsize;
}

SpatialIndex::ISpatialIndex* SpatialIndex::RTree::returnRTree(SpatialIndex::IStorageManager& sm, Tools::PropertySet& ps)
{
	SpatialIndex::ISpatialIndex* si = new SpatialIndex::RTree::RTree(sm, ps);
	return si;
}

SpatialIndex::ISpatialIndex* SpatialIndex::RTree::createNewRTree(
	SpatialIndex::IStorageManager& sm,
	double fillFactor,
	uint32_t indexCapacity,
	uint32_t leafCapacity,
	uint32_t dimension,
	RTreeVariant rv,
	id_type& indexIdentifier)
{
	Tools::Variant var;
	Tools::PropertySet ps;

	var.m_varType = Tools::VT_DOUBLE;
	var.m_val.dblVal = fillFactor;
	ps.setProperty("FillFactor", var);

	var.m_varType = Tools::VT_ULONG;
	var.m_val.ulVal = indexCapacity;
	ps.setProperty("IndexCapacity", var);

	var.m_varType = Tools::VT_ULONG;
	var.m_val.ulVal = leafCapacity;
	ps.setProperty("LeafCapacity", var);

	var.m_varType = Tools::VT_ULONG;
	var.m_val.ulVal = dimension;
	ps.setProperty("Dimension", var);

	var.m_varType = Tools::VT_LONG;
	var.m_val.lVal = rv;
	ps.setProperty("TreeVariant", var);

	ISpatialIndex* ret = returnRTree(sm, ps);

	var.m_varType = Tools::VT_LONGLONG;
	var = ps.getProperty("IndexIdentifier");
	indexIdentifier = var.m_val.llVal;

	return ret;
}

SpatialIndex::ISpatialIndex* SpatialIndex::RTree::createAndBulkLoadNewRTree(
	BulkLoadMethod m,  		// BLM_STR
	IDataStream& stream,
	SpatialIndex::IStorageManager& sm,  // cache
	double fillFactor, 		// 0.7
	uint32_t indexCapacity, // 100
	uint32_t leafCapacity,  // 100
	uint32_t dimension,     // 2
	SpatialIndex::RTree::RTreeVariant rv, // RV_RSTAR
	id_type& indexIdentifier)
{
	SpatialIndex::ISpatialIndex* tree = createNewRTree(sm, fillFactor, indexCapacity, leafCapacity, dimension, rv, indexIdentifier);

	uint32_t bindex = static_cast<uint32_t>(std::floor(static_cast<double>(indexCapacity * fillFactor))); // 70
	uint32_t bleaf = static_cast<uint32_t>(std::floor(static_cast<double>(leafCapacity * fillFactor)));   // 70

	SpatialIndex::RTree::BulkLoader bl;

	switch (m)
	{
	case BLM_STR:
		bl.bulkLoadUsingSTR(static_cast<RTree*>(tree), stream, bindex, bleaf, 10000, 100);
		break;
	default:
		throw Tools::IllegalArgumentException("createAndBulkLoadNewRTree: Unknown bulk load method.");
		break;
	}

	return tree;
}

SpatialIndex::ISpatialIndex* SpatialIndex::RTree::createAndBulkLoadNewRTree(
	BulkLoadMethod m,
	IDataStream& stream,
	SpatialIndex::IStorageManager& sm,
	Tools::PropertySet& ps,
	id_type& indexIdentifier)
{
	Tools::Variant var;
	RTreeVariant rv(RV_LINEAR);
	double fillFactor(0.0);
	uint32_t indexCapacity(0);
	uint32_t leafCapacity(0);
	uint32_t dimension(0);
	uint32_t pageSize(0);
	uint32_t numberOfPages(0);

	// tree variant
	var = ps.getProperty("TreeVariant");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (
			var.m_varType != Tools::VT_LONG ||
			(var.m_val.lVal != RV_LINEAR &&
			var.m_val.lVal != RV_QUADRATIC &&
			var.m_val.lVal != RV_RSTAR))
			throw Tools::IllegalArgumentException("createAndBulkLoadNewRTree: Property TreeVariant must be Tools::VT_LONG and of RTreeVariant type");

		rv = static_cast<RTreeVariant>(var.m_val.lVal);
	}

	// fill factor
	// it cannot be larger than 50%, since linear and quadratic split algorithms
	// require assigning to both nodes the same number of entries.
	var = ps.getProperty("FillFactor");
	if (var.m_varType != Tools::VT_EMPTY)
	{
	    if (var.m_varType != Tools::VT_DOUBLE)
            throw Tools::IllegalArgumentException("createAndBulkLoadNewRTree: Property FillFactor was not of type Tools::VT_DOUBLE");

        if (var.m_val.dblVal <= 0.0)
            throw Tools::IllegalArgumentException("createAndBulkLoadNewRTree: Property FillFactor was less than 0.0");

        if (((rv == RV_LINEAR || rv == RV_QUADRATIC) && var.m_val.dblVal > 0.5))
            throw Tools::IllegalArgumentException( "createAndBulkLoadNewRTree: Property FillFactor must be in range (0.0, 0.5) for LINEAR or QUADRATIC index types");
        if ( var.m_val.dblVal >= 1.0)
            throw Tools::IllegalArgumentException("createAndBulkLoadNewRTree: Property FillFactor must be in range (0.0, 1.0) for RSTAR index type");
		fillFactor = var.m_val.dblVal;
	}

	// index capacity
	var = ps.getProperty("IndexCapacity");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (var.m_varType != Tools::VT_ULONG || var.m_val.ulVal < 4)
			throw Tools::IllegalArgumentException("createAndBulkLoadNewRTree: Property IndexCapacity must be Tools::VT_ULONG and >= 4");

		indexCapacity = var.m_val.ulVal;
	}

	// leaf capacity
	var = ps.getProperty("LeafCapacity");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (var.m_varType != Tools::VT_ULONG || var.m_val.ulVal < 4)
			throw Tools::IllegalArgumentException("createAndBulkLoadNewRTree: Property LeafCapacity must be Tools::VT_ULONG and >= 4");

		leafCapacity = var.m_val.ulVal;
	}

	// dimension
	var = ps.getProperty("Dimension");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (var.m_varType != Tools::VT_ULONG)
			throw Tools::IllegalArgumentException("createAndBulkLoadNewRTree: Property Dimension must be Tools::VT_ULONG");
		if (var.m_val.ulVal <= 1)
			throw Tools::IllegalArgumentException("createAndBulkLoadNewRTree: Property Dimension must be greater than 1");

		dimension = var.m_val.ulVal;
	}

	// page size
	var = ps.getProperty("ExternalSortBufferPageSize");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (var.m_varType != Tools::VT_ULONG)
			throw Tools::IllegalArgumentException("createAndBulkLoadNewRTree: Property ExternalSortBufferPageSize must be Tools::VT_ULONG");
		if (var.m_val.ulVal <= 1)
			throw Tools::IllegalArgumentException("createAndBulkLoadNewRTree: Property ExternalSortBufferPageSize must be greater than 1");

		pageSize = var.m_val.ulVal;
	}

	// number of pages
	var = ps.getProperty("ExternalSortBufferTotalPages");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (var.m_varType != Tools::VT_ULONG)
			throw Tools::IllegalArgumentException("createAndBulkLoadNewRTree: Property ExternalSortBufferTotalPages must be Tools::VT_ULONG");
		if (var.m_val.ulVal <= 1)
			throw Tools::IllegalArgumentException("createAndBulkLoadNewRTree: Property ExternalSortBufferTotalPages must be greater than 1");

		numberOfPages = var.m_val.ulVal;
	}

	SpatialIndex::ISpatialIndex* tree = createNewRTree(sm, fillFactor, indexCapacity, leafCapacity, dimension, rv, indexIdentifier);

	uint32_t bindex = static_cast<uint32_t>(std::floor(static_cast<double>(indexCapacity * fillFactor)));
	uint32_t bleaf = static_cast<uint32_t>(std::floor(static_cast<double>(leafCapacity * fillFactor)));

	SpatialIndex::RTree::BulkLoader bl;

	switch (m)
	{
	case BLM_STR:
		bl.bulkLoadUsingSTR(static_cast<RTree*>(tree), stream, bindex, bleaf, pageSize, numberOfPages);
		break;
	default:
		throw Tools::IllegalArgumentException("createAndBulkLoadNewRTree: Unknown bulk load method.");
		break;
	}

	return tree;
}

SpatialIndex::ISpatialIndex* SpatialIndex::RTree::loadRTree(IStorageManager& sm, id_type indexIdentifier)
{
	Tools::Variant var;
	Tools::PropertySet ps;

	var.m_varType = Tools::VT_LONGLONG;
	var.m_val.llVal = indexIdentifier;
	ps.setProperty("IndexIdentifier", var);

	return returnRTree(sm, ps);
}

SpatialIndex::RTree::RTree::RTree(IStorageManager& sm, Tools::PropertySet& ps) :
	m_pStorageManager(&sm),
	m_rootID(StorageManager::NewPage),
	m_headerID(StorageManager::NewPage),
	m_treeVariant(RV_RSTAR),
	m_fillFactor(0.7),
	m_indexCapacity(100),
	m_leafCapacity(100),
	m_nearMinimumOverlapFactor(32),
	m_splitDistributionFactor(0.4),
	m_reinsertFactor(0.3),
	m_dimension(2),
	m_bTightMBRs(true),
	m_pointPool(500),
	m_regionPool(1000),
	m_indexPool(100),
	m_leafPool(100)
{
	Tools::Variant var = ps.getProperty("IndexIdentifier");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (var.m_varType == Tools::VT_LONGLONG) m_headerID = var.m_val.llVal;
		else if (var.m_varType == Tools::VT_LONG) m_headerID = var.m_val.lVal;
			// for backward compatibility only.
		else throw Tools::IllegalArgumentException("RTree: Property IndexIdentifier must be Tools::VT_LONGLONG");

		initOld(ps);
	}
	else
	{
		initNew(ps);
		var.m_varType = Tools::VT_LONGLONG;
		var.m_val.llVal = m_headerID;
		ps.setProperty("IndexIdentifier", var);
	}

	// // ============= Generate Z-order ================
	// uint32_t n = 16;
	// std::vector<uint32_t> array1(n * n);
	// float zorder_margin = (float)100/n;
	// for(uint32_t i = 0; i < n; i++)
	// 	for(uint32_t j = 0; j < n; j++)
	// 		array1[i * n + j] = i * n + j;
	// zorder.generate_z_order_curve(array1.data(), pow(n, 2));
	// zorder.print_t_z_order();
	// //============= Generate Z-order Done ================
}

SpatialIndex::RTree::RTree::~RTree()
{
	storeHeader();
}

//
// ISpatialIndex interface
//

void SpatialIndex::RTree::RTree::insertData(uint32_t len, const uint8_t* pData, const IShape& shape, id_type id)
{
	if (shape.getDimension() != m_dimension) 
		throw Tools::IllegalArgumentException("insertData: Shape has the wrong number of dimensions.");

	// convert the shape into a Region (R-Trees index regions only; i.e., approximations of the shapes).
	RegionPtr mbr = m_regionPool.acquire();
	shape.getMBR(*mbr);

	uint8_t* buffer = nullptr;

	if (len > 0)
	{
		buffer = new uint8_t[len];
		memcpy(buffer, pData, len);
	}

	insertData_impl(len, buffer, *mbr, id);
		// the buffer is stored in the tree. Do not delete here.
}

bool SpatialIndex::RTree::RTree::deleteData(const IShape& shape, id_type id)
{
	if (shape.getDimension() != m_dimension) throw Tools::IllegalArgumentException("deleteData: Shape has the wrong number of dimensions.");

	RegionPtr mbr = m_regionPool.acquire();
	shape.getMBR(*mbr); // shape copy to mbr
	bool ret = deleteData_impl(*mbr, id);

	return ret;
}

/*
void SpatialIndex::RTree::RTree::internalNodesQuery(const IShape& query, IVisitor& v)
{
	if (query.getDimension() != m_dimension) throw Tools::IllegalArgumentException("containsWhatQuery: Shape has the wrong number of dimensions.");

#ifdef HAVE_PTHREAD_H
	Tools::LockGuard lock(&m_lock);
#endif

	try
	{
		std::stack<NodePtr> st;
		NodePtr root = readNode(m_rootID);
		st.push(root);

		while (! st.empty())
		{
			NodePtr n = st.top(); st.pop();

			if(query.containsShape(n->m_nodeMBR))
			{
				IdVisitor vId = IdVisitor();
				visitSubTree(n, vId);
				const uint64_t nObj = vId.GetResultCount();
				uint64_t *obj = new uint64_t[nObj];
				std::copy(vId.GetResults().begin(), vId.GetResults().end(), obj);

				Data data = Data((uint32_t)(sizeof(uint64_t) * nObj), (uint8_t *) obj, n->m_nodeMBR, n->getIdentifier());
				v.visitData(data);
				++(m_stats.m_u64QueryResults);
			}
			else
			{
				if(n->m_level == 0)
				{
					for (uint32_t cChild = 0; cChild < n->m_children; ++cChild)
					{
						if(query.containsShape(*(n->m_ptrMBR[cChild])))
						{
							Data data = Data(sizeof(id_type), (uint8_t *) &n->m_pIdentifier[cChild], *(n->m_ptrMBR[cChild]), n->getIdentifier());
							v.visitData(data);
							++(m_stats.m_u64QueryResults);
						}
					}
				}
				else //not a leaf
				{
					if(query.intersectsShape(n->m_nodeMBR))
					{
						for (uint32_t cChild = 0; cChild < n->m_children; ++cChild)
						{
							st.push(readNode(n->m_pIdentifier[cChild]));
						}
					}
				}
			}
		}

	}
	catch (...)
	{
		throw;
	}
}
*/

void SpatialIndex::RTree::RTree::containsWhatQuery(const IShape& query, IVisitor& v)
{
	if (query.getDimension() != m_dimension) throw Tools::IllegalArgumentException("containsWhatQuery: Shape has the wrong number of dimensions.");

	try
	{
		std::stack<NodePtr> st;
		NodePtr root = readNode(m_rootID);
		st.push(root);

		while (! st.empty())
		{
			NodePtr n = st.top(); st.pop();

			if(n->m_level == 0)
			{
				v.visitNode(*n);

				for (uint32_t cChild = 0; cChild < n->m_children; ++cChild)
				{
					if(query.containsShape(*(n->m_ptrMBR[cChild])))
					{
						Data data = Data(n->m_pDataLength[cChild], n->m_pData[cChild], *(n->m_ptrMBR[cChild]), n->m_pIdentifier[cChild]);
						v.visitData(data);
						++(m_stats.m_u64QueryResults);
					}
				}
			}
			else //not a leaf
			{
				if(query.containsShape(n->m_nodeMBR))
				{
					visitSubTree(n, v);
				}
				else if(query.intersectsShape(n->m_nodeMBR))
				{
					v.visitNode(*n);

					for (uint32_t cChild = 0; cChild < n->m_children; ++cChild)
					{
						st.push(readNode(n->m_pIdentifier[cChild]));
					}
				}
			}
		}

	}
	catch (...)
	{
		throw;
	}
}
void SpatialIndex::RTree::RTree::intersectsWithQuery(const IShape& query, IVisitor& v)
{
	if (query.getDimension() != m_dimension) throw Tools::IllegalArgumentException("intersectsWithQuery: Shape has the wrong number of dimensions.");
	rangeQuery(IntersectionQuery, query, v);
}

void SpatialIndex::RTree::RTree::pointLocationQuery(const Point& query, IVisitor& v)
{
	if (query.m_dimension != m_dimension) throw Tools::IllegalArgumentException("pointLocationQuery: Shape has the wrong number of dimensions.");
	Region r(query, query);
	rangeQuery(IntersectionQuery, r, v);
}

void SpatialIndex::RTree::RTree::nearestNeighborQuery(uint32_t k, const IShape& query, IVisitor& v, INearestNeighborComparator& nnc)
{
	if (query.getDimension() != m_dimension) throw Tools::IllegalArgumentException("nearestNeighborQuery: Shape has the wrong number of dimensions.");

	auto ascending = [](const NNEntry* lhs, const NNEntry* rhs) { return lhs->m_minDist > rhs->m_minDist;  };
	std::priority_queue<NNEntry*, std::vector<NNEntry*>, decltype(ascending)> queue(ascending);

	queue.push(new NNEntry(m_rootID, nullptr, 0.0));

	uint32_t count = 0;
	double knearest = 0.0;

	while (! queue.empty())
	{
		NNEntry* pFirst = queue.top();

		// report all nearest neighbors with equal greatest distances.
		// (neighbors can be more than k, if many happen to have the same greatest distance).
		if (count >= k && pFirst->m_minDist > knearest)	break;

		queue.pop();

		if (pFirst->m_pEntry == nullptr)
		{
			// n is a leaf or an index.
			NodePtr n = readNode(pFirst->m_id);
			v.visitNode(*n);

			for (uint32_t cChild = 0; cChild < n->m_children; ++cChild)
			{
				if (n->m_level == 0)
				{
					Data* e = new Data(n->m_pDataLength[cChild], n->m_pData[cChild], *(n->m_ptrMBR[cChild]), n->m_pIdentifier[cChild]);
					// we need to compare the query with the actual data entry here, so we call the
					// appropriate getMinimumDistance method of NearestNeighborComparator.
					queue.push(new NNEntry(n->m_pIdentifier[cChild], e, nnc.getMinimumDistance(query, *e)));
				}
				else
				{
					queue.push(new NNEntry(n->m_pIdentifier[cChild], nullptr, nnc.getMinimumDistance(query, *(n->m_ptrMBR[cChild]))));
				}
			}
		}
		else
		{
			v.visitData(*(static_cast<IData*>(pFirst->m_pEntry)));
			++(m_stats.m_u64QueryResults);
			++count;
			knearest = pFirst->m_minDist;
			delete pFirst->m_pEntry;
		}

		delete pFirst;
	}

	while (! queue.empty())
	{
		NNEntry* e = queue.top(); queue.pop();
		if (e->m_pEntry != nullptr) delete e->m_pEntry;
		delete e;
	}
}

void SpatialIndex::RTree::RTree::nearestNeighborQuery(uint32_t k, const IShape& query, IVisitor& v)
{
	if (query.getDimension() != m_dimension) throw Tools::IllegalArgumentException("nearestNeighborQuery: Shape has the wrong number of dimensions.");
	NNComparator nnc;
	nearestNeighborQuery(k, query, v, nnc);
}


void SpatialIndex::RTree::RTree::selfJoinQuery(const IShape& query, IVisitor& v)
{
	if (query.getDimension() != m_dimension)
		throw Tools::IllegalArgumentException("selfJoinQuery: Shape has the wrong number of dimensions.");

	RegionPtr mbr = m_regionPool.acquire();
	query.getMBR(*mbr);
	selfJoinQuery(m_rootID, m_rootID, *mbr, v);
}

void SpatialIndex::RTree::RTree::queryStrategy(IQueryStrategy& qs)
{
	id_type next = m_rootID;
	bool hasNext = true;

	while (hasNext)
	{
		NodePtr n = readNode(next);
		qs.getNextEntry(*n, next, hasNext);
	}
}

void SpatialIndex::RTree::RTree::getIndexProperties(Tools::PropertySet& out) const
{
	Tools::Variant var;

	// dimension
	var.m_varType = Tools::VT_ULONG;
	var.m_val.ulVal = m_dimension;
	out.setProperty("Dimension", var);

	// index capacity
	var.m_varType = Tools::VT_ULONG;
	var.m_val.ulVal = m_indexCapacity;
	out.setProperty("IndexCapacity", var);

	// leaf capacity
	var.m_varType = Tools::VT_ULONG;
	var.m_val.ulVal = m_leafCapacity;
	out.setProperty("LeafCapacity", var);

	// R-tree variant
	var.m_varType = Tools::VT_LONG;
	var.m_val.lVal = m_treeVariant;
	out.setProperty("TreeVariant", var);

	// fill factor
	var.m_varType = Tools::VT_DOUBLE;
	var.m_val.dblVal = m_fillFactor;
	out.setProperty("FillFactor", var);

	// near minimum overlap factor
	var.m_varType = Tools::VT_ULONG;
	var.m_val.ulVal = m_nearMinimumOverlapFactor;
	out.setProperty("NearMinimumOverlapFactor", var);

	// split distribution factor
	var.m_varType = Tools::VT_DOUBLE;
	var.m_val.dblVal = m_splitDistributionFactor;
	out.setProperty("SplitDistributionFactor", var);

	// reinsert factor
	var.m_varType = Tools::VT_DOUBLE;
	var.m_val.dblVal = m_reinsertFactor;
	out.setProperty("ReinsertFactor", var);

	// tight MBRs
	var.m_varType = Tools::VT_BOOL;
	var.m_val.blVal = m_bTightMBRs;
	out.setProperty("EnsureTightMBRs", var);

	// index pool capacity
	var.m_varType = Tools::VT_ULONG;
	var.m_val.ulVal = m_indexPool.getCapacity();
	out.setProperty("IndexPoolCapacity", var);

	// leaf pool capacity
	var.m_varType = Tools::VT_ULONG;
	var.m_val.ulVal = m_leafPool.getCapacity();
	out.setProperty("LeafPoolCapacity", var);

	// region pool capacity
	var.m_varType = Tools::VT_ULONG;
	var.m_val.ulVal = m_regionPool.getCapacity();
	out.setProperty("RegionPoolCapacity", var);

	// point pool capacity
	var.m_varType = Tools::VT_ULONG;
	var.m_val.ulVal = m_pointPool.getCapacity();
	out.setProperty("PointPoolCapacity", var);

	var.m_varType = Tools::VT_LONGLONG;
	var.m_val.llVal = m_headerID;
	out.setProperty("IndexIdentifier", var);

}

void SpatialIndex::RTree::RTree::addCommand(ICommand* pCommand, CommandType ct)
{
	switch (ct)
	{
		case CT_NODEREAD:
			m_readNodeCommands.push_back(std::shared_ptr<ICommand>(pCommand));
			break;
		case CT_NODEWRITE:
			m_writeNodeCommands.push_back(std::shared_ptr<ICommand>(pCommand));
			break;
		case CT_NODEDELETE:
			m_deleteNodeCommands.push_back(std::shared_ptr<ICommand>(pCommand));
			break;
	}
}

bool SpatialIndex::RTree::RTree::isIndexValid()
{
	bool ret = true;
	std::stack<ValidateEntry> st;
	NodePtr root = readNode(m_rootID);

	if (root->m_level != m_stats.m_u32TreeHeight - 1)
	{
		std::cerr << "Invalid tree height." << std::endl;
		return false;
	}

	std::map<uint32_t, uint32_t> nodesInLevel;
	nodesInLevel.insert(std::pair<uint32_t, uint32_t>(root->m_level, 1));

	ValidateEntry e(root->m_nodeMBR, root);
	st.push(e);

	while (! st.empty())
	{
		e = st.top(); st.pop();

		Region tmpRegion;
		tmpRegion = m_infiniteRegion;

		for (uint32_t cDim = 0; cDim < tmpRegion.m_dimension; ++cDim)
		{
			tmpRegion.m_pLow[cDim] = std::numeric_limits<double>::max();
			tmpRegion.m_pHigh[cDim] = -std::numeric_limits<double>::max();

			for (uint32_t cChild = 0; cChild < e.m_pNode->m_children; ++cChild)
			{
				tmpRegion.m_pLow[cDim] = std::min(tmpRegion.m_pLow[cDim], e.m_pNode->m_ptrMBR[cChild]->m_pLow[cDim]);
				tmpRegion.m_pHigh[cDim] = std::max(tmpRegion.m_pHigh[cDim], e.m_pNode->m_ptrMBR[cChild]->m_pHigh[cDim]);
			}
		}

		if (! (tmpRegion == e.m_pNode->m_nodeMBR))
		{
			std::cerr << "Invalid parent information." << std::endl;
			ret = false;
		}
		else if (! (tmpRegion == e.m_parentMBR))
		{
			std::cerr << "Error in parent." << std::endl;
			ret = false;
		}

		if (e.m_pNode->m_level != 0)
		{
			for (uint32_t cChild = 0; cChild < e.m_pNode->m_children; ++cChild)
			{
				NodePtr ptrN = readNode(e.m_pNode->m_pIdentifier[cChild]);
				ValidateEntry tmpEntry(*(e.m_pNode->m_ptrMBR[cChild]), ptrN);

				std::map<uint32_t, uint32_t>::iterator itNodes = nodesInLevel.find(tmpEntry.m_pNode->m_level);

				if (itNodes == nodesInLevel.end())
				{
					nodesInLevel.insert(std::pair<uint32_t, uint32_t>(tmpEntry.m_pNode->m_level, 1l));
				}
				else
				{
					// equals nodesInLevel[tmpEntry.m_pNode->m_level]++;
					nodesInLevel[tmpEntry.m_pNode->m_level] = nodesInLevel[tmpEntry.m_pNode->m_level] + 1; 
				}

				st.push(tmpEntry);
			}
		}
	}

	uint32_t nodes = 0;
	for (uint32_t cLevel = 0; cLevel < m_stats.m_u32TreeHeight; ++cLevel)
	{
		if (nodesInLevel[cLevel] != m_stats.m_nodesInLevel[cLevel])
		{
			std::cerr << "Invalid nodesInLevel information." << std::endl;
			ret = false;
		}

		nodes += m_stats.m_nodesInLevel[cLevel];
	}

	if (nodes != m_stats.m_u32Nodes)
	{
		std::cerr << "Invalid number of nodes information." << std::endl;
		ret = false;
	}

	return ret;
}

void SpatialIndex::RTree::RTree::getStatistics(IStatistics** out) const
{
	*out = new Statistics(m_stats);
}

void SpatialIndex::RTree::RTree::flush()
{
	storeHeader();
}

void SpatialIndex::RTree::RTree::initNew(Tools::PropertySet& ps)
{
	Tools::Variant var;

	// tree variant
	var = ps.getProperty("TreeVariant");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (
			var.m_varType != Tools::VT_LONG ||
			(var.m_val.lVal != RV_LINEAR &&
			var.m_val.lVal != RV_QUADRATIC &&
			var.m_val.lVal != RV_RSTAR))
			throw Tools::IllegalArgumentException("initNew: Property TreeVariant must be Tools::VT_LONG and of RTreeVariant type");

		m_treeVariant = static_cast<RTreeVariant>(var.m_val.lVal);
	}

	// fill factor
	// it cannot be larger than 50%, since linear and quadratic split algorithms
	// require assigning to both nodes the same number of entries.
	var = ps.getProperty("FillFactor");
	if (var.m_varType != Tools::VT_EMPTY)
	{
	    if (var.m_varType != Tools::VT_DOUBLE)
            throw Tools::IllegalArgumentException("initNew: Property FillFactor was not of type Tools::VT_DOUBLE");

        if (var.m_val.dblVal <= 0.0)
            throw Tools::IllegalArgumentException("initNew: Property FillFactor was less than 0.0");

        if (((m_treeVariant == RV_LINEAR || m_treeVariant == RV_QUADRATIC) && var.m_val.dblVal > 0.5))
            throw Tools::IllegalArgumentException(  "initNew: Property FillFactor must be in range "
                                                    "(0.0, 0.5) for LINEAR or QUADRATIC index types");
        if ( var.m_val.dblVal >= 1.0)
            throw Tools::IllegalArgumentException(  "initNew: Property FillFactor must be in range "
                                                    "(0.0, 1.0) for RSTAR index type");
		m_fillFactor = var.m_val.dblVal;
	}

	// index capacity
	var = ps.getProperty("IndexCapacity");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (var.m_varType != Tools::VT_ULONG || var.m_val.ulVal < 4)
			throw Tools::IllegalArgumentException("initNew: Property IndexCapacity must be Tools::VT_ULONG and >= 4");

		m_indexCapacity = var.m_val.ulVal;
	}

	// leaf capacity
	var = ps.getProperty("LeafCapacity");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (var.m_varType != Tools::VT_ULONG || var.m_val.ulVal < 4)
			throw Tools::IllegalArgumentException("initNew: Property LeafCapacity must be Tools::VT_ULONG and >= 4");

		m_leafCapacity = var.m_val.ulVal;
	}

	// near minimum overlap factor
	var = ps.getProperty("NearMinimumOverlapFactor");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (
			var.m_varType != Tools::VT_ULONG ||
			var.m_val.ulVal < 1 ||
			var.m_val.ulVal > m_indexCapacity ||
			var.m_val.ulVal > m_leafCapacity)
			throw Tools::IllegalArgumentException("initNew: Property NearMinimumOverlapFactor must be Tools::VT_ULONG and less than both index and leaf capacities");

		m_nearMinimumOverlapFactor = var.m_val.ulVal;
	}

	// split distribution factor
	var = ps.getProperty("SplitDistributionFactor");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (
			var.m_varType != Tools::VT_DOUBLE ||
			var.m_val.dblVal <= 0.0 ||
			var.m_val.dblVal >= 1.0)
			throw Tools::IllegalArgumentException("initNew: Property SplitDistributionFactor must be Tools::VT_DOUBLE and in (0.0, 1.0)");

		m_splitDistributionFactor = var.m_val.dblVal;
	}

	// reinsert factor
	var = ps.getProperty("ReinsertFactor");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (
			var.m_varType != Tools::VT_DOUBLE ||
			var.m_val.dblVal <= 0.0 ||
			var.m_val.dblVal >= 1.0)
			throw Tools::IllegalArgumentException("initNew: Property ReinsertFactor must be Tools::VT_DOUBLE and in (0.0, 1.0)");

		m_reinsertFactor = var.m_val.dblVal;
	}

	// dimension
	var = ps.getProperty("Dimension");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (var.m_varType != Tools::VT_ULONG)
			throw Tools::IllegalArgumentException("initNew: Property Dimension must be Tools::VT_ULONG");
		if (var.m_val.ulVal <= 1)
			throw Tools::IllegalArgumentException("initNew: Property Dimension must be greater than 1");

		m_dimension = var.m_val.ulVal;
	}

	// tight MBRs
	var = ps.getProperty("EnsureTightMBRs");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (var.m_varType != Tools::VT_BOOL)
			throw Tools::IllegalArgumentException("initNew: Property EnsureTightMBRs must be Tools::VT_BOOL");

		m_bTightMBRs = var.m_val.blVal;
	}

	// index pool capacity
	var = ps.getProperty("IndexPoolCapacity");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (var.m_varType != Tools::VT_ULONG)
			throw Tools::IllegalArgumentException("initNew: Property IndexPoolCapacity must be Tools::VT_ULONG");

		m_indexPool.setCapacity(var.m_val.ulVal);
	}

	// leaf pool capacity
	var = ps.getProperty("LeafPoolCapacity");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (var.m_varType != Tools::VT_ULONG)
			throw Tools::IllegalArgumentException("initNew: Property LeafPoolCapacity must be Tools::VT_ULONG");

		m_leafPool.setCapacity(var.m_val.ulVal);
	}

	// region pool capacity
	var = ps.getProperty("RegionPoolCapacity");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (var.m_varType != Tools::VT_ULONG)
			throw Tools::IllegalArgumentException("initNew: Property RegionPoolCapacity must be Tools::VT_ULONG");

		m_regionPool.setCapacity(var.m_val.ulVal);
	}

	// point pool capacity
	var = ps.getProperty("PointPoolCapacity");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (var.m_varType != Tools::VT_ULONG)
			throw Tools::IllegalArgumentException("initNew: Property PointPoolCapacity must be Tools::VT_ULONG");

		m_pointPool.setCapacity(var.m_val.ulVal);
	}

	m_infiniteRegion.makeInfinite(m_dimension);

	m_stats.m_u32TreeHeight = 1;
	m_stats.m_nodesInLevel.push_back(0);

	Leaf root(this, -1); //@@: this -1 in writeNode(.) will replaced with a new valid ID
	m_rootID = writeNode_complete(&root);

	storeHeader();
}

void SpatialIndex::RTree::RTree::initOld(Tools::PropertySet& ps)
{
	loadHeader();

	// only some of the properties may be changed.
	// the rest are just ignored.

	Tools::Variant var;

	// tree variant
	var = ps.getProperty("TreeVariant");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (
			var.m_varType != Tools::VT_LONG ||
			(var.m_val.lVal != RV_LINEAR &&
			 var.m_val.lVal != RV_QUADRATIC &&
			 var.m_val.lVal != RV_RSTAR))
			throw Tools::IllegalArgumentException("initOld: Property TreeVariant must be Tools::VT_LONG and of RTreeVariant type");

		m_treeVariant = static_cast<RTreeVariant>(var.m_val.lVal);
	}

	// near minimum overlap factor
	var = ps.getProperty("NearMinimumOverlapFactor");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if ( 
			var.m_varType != Tools::VT_ULONG ||
			var.m_val.ulVal < 1 ||
			var.m_val.ulVal > m_indexCapacity ||
			var.m_val.ulVal > m_leafCapacity)
			throw Tools::IllegalArgumentException("initOld: Property NearMinimumOverlapFactor must be Tools::VT_ULONG and less than both index and leaf capacities");

		m_nearMinimumOverlapFactor = var.m_val.ulVal;
	}

	// split distribution factor
	var = ps.getProperty("SplitDistributionFactor");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (var.m_varType != Tools::VT_DOUBLE || var.m_val.dblVal <= 0.0 || var.m_val.dblVal >= 1.0)
			throw Tools::IllegalArgumentException("initOld: Property SplitDistributionFactor must be Tools::VT_DOUBLE and in (0.0, 1.0)");

		m_splitDistributionFactor = var.m_val.dblVal;
	}

	// reinsert factor
	var = ps.getProperty("ReinsertFactor");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (var.m_varType != Tools::VT_DOUBLE || var.m_val.dblVal <= 0.0 || var.m_val.dblVal >= 1.0)
			throw Tools::IllegalArgumentException("initOld: Property ReinsertFactor must be Tools::VT_DOUBLE and in (0.0, 1.0)");

		m_reinsertFactor = var.m_val.dblVal;
	}

	// tight MBRs
	var = ps.getProperty("EnsureTightMBRs");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (var.m_varType != Tools::VT_BOOL) throw Tools::IllegalArgumentException("initOld: Property EnsureTightMBRs must be Tools::VT_BOOL");

		m_bTightMBRs = var.m_val.blVal;
	}

	// index pool capacity
	var = ps.getProperty("IndexPoolCapacity");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (var.m_varType != Tools::VT_ULONG) throw Tools::IllegalArgumentException("initOld: Property IndexPoolCapacity must be Tools::VT_ULONG");

		m_indexPool.setCapacity(var.m_val.ulVal);
	}

	// leaf pool capacity
	var = ps.getProperty("LeafPoolCapacity");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (var.m_varType != Tools::VT_ULONG) throw Tools::IllegalArgumentException("initOld: Property LeafPoolCapacity must be Tools::VT_ULONG");

		m_leafPool.setCapacity(var.m_val.ulVal);
	}

	// region pool capacity
	var = ps.getProperty("RegionPoolCapacity");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (var.m_varType != Tools::VT_ULONG) throw Tools::IllegalArgumentException("initOld: Property RegionPoolCapacity must be Tools::VT_ULONG");

		m_regionPool.setCapacity(var.m_val.ulVal);
	}

	// point pool capacity
	var = ps.getProperty("PointPoolCapacity");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (var.m_varType != Tools::VT_ULONG) throw Tools::IllegalArgumentException("initOld: Property PointPoolCapacity must be Tools::VT_ULONG");

		m_pointPool.setCapacity(var.m_val.ulVal);
	}

	m_infiniteRegion.makeInfinite(m_dimension);
}

void SpatialIndex::RTree::RTree::storeHeader()
{
	const uint32_t headerSize =
		sizeof(id_type) +						// m_rootID
		sizeof(RTreeVariant) +					// m_treeVariant
		sizeof(double) +						// m_fillFactor
		sizeof(uint32_t) +						// m_indexCapacity
		sizeof(uint32_t) +						// m_leafCapacity
		sizeof(uint32_t) +						// m_nearMinimumOverlapFactor
		sizeof(double) +						// m_splitDistributionFactor
		sizeof(double) +						// m_reinsertFactor
		sizeof(uint32_t) +						// m_dimension
		sizeof(char) +							// m_bTightMBRs
		sizeof(uint32_t) +						// m_stats.m_nodes
		sizeof(uint64_t) +						// m_stats.m_data
		sizeof(uint32_t) +						// m_stats.m_treeHeight
		m_stats.m_u32TreeHeight * sizeof(uint32_t);	// m_stats.m_nodesInLevel

	uint8_t* header = new uint8_t[headerSize];
	uint8_t* ptr = header;

	memcpy(ptr, &m_rootID, sizeof(id_type));
	ptr += sizeof(id_type);
	memcpy(ptr, &m_treeVariant, sizeof(RTreeVariant));
	ptr += sizeof(RTreeVariant);
	memcpy(ptr, &m_fillFactor, sizeof(double));
	ptr += sizeof(double);
	memcpy(ptr, &m_indexCapacity, sizeof(uint32_t));
	ptr += sizeof(uint32_t);
	memcpy(ptr, &m_leafCapacity, sizeof(uint32_t));
	ptr += sizeof(uint32_t);
	memcpy(ptr, &m_nearMinimumOverlapFactor, sizeof(uint32_t));
	ptr += sizeof(uint32_t);
	memcpy(ptr, &m_splitDistributionFactor, sizeof(double));
	ptr += sizeof(double);
	memcpy(ptr, &m_reinsertFactor, sizeof(double));
	ptr += sizeof(double);
	memcpy(ptr, &m_dimension, sizeof(uint32_t));
	ptr += sizeof(uint32_t);
	char c = (char) m_bTightMBRs;
	memcpy(ptr, &c, sizeof(char));
	ptr += sizeof(char);
	memcpy(ptr, &(m_stats.m_u32Nodes), sizeof(uint32_t));
	ptr += sizeof(uint32_t);
	memcpy(ptr, &(m_stats.m_u64Data), sizeof(uint64_t));
	ptr += sizeof(uint64_t);
	memcpy(ptr, &(m_stats.m_u32TreeHeight), sizeof(uint32_t));
	ptr += sizeof(uint32_t);

	for (uint32_t cLevel = 0; cLevel < m_stats.m_u32TreeHeight; ++cLevel)
	{
		memcpy(ptr, &(m_stats.m_nodesInLevel[cLevel]), sizeof(uint32_t));
		ptr += sizeof(uint32_t);
	}

	m_pStorageManager->storeByteArray(m_headerID, headerSize, header);

	delete[] header;
}

void SpatialIndex::RTree::RTree::loadHeader()
{
	uint32_t headerSize;
	uint8_t* header = nullptr;
	m_pStorageManager->loadByteArray(m_headerID, headerSize, &header);

	uint8_t* ptr = header;

	memcpy(&m_rootID, ptr, sizeof(id_type));
	ptr += sizeof(id_type);
	memcpy(&m_treeVariant, ptr, sizeof(RTreeVariant));
	ptr += sizeof(RTreeVariant);
	memcpy(&m_fillFactor, ptr, sizeof(double));
	ptr += sizeof(double);
	memcpy(&m_indexCapacity, ptr, sizeof(uint32_t));
	ptr += sizeof(uint32_t);
	memcpy(&m_leafCapacity, ptr, sizeof(uint32_t));
	ptr += sizeof(uint32_t);
	memcpy(&m_nearMinimumOverlapFactor, ptr, sizeof(uint32_t));
	ptr += sizeof(uint32_t);
	memcpy(&m_splitDistributionFactor, ptr, sizeof(double));
	ptr += sizeof(double);
	memcpy(&m_reinsertFactor, ptr, sizeof(double));
	ptr += sizeof(double);
	memcpy(&m_dimension, ptr, sizeof(uint32_t));
	ptr += sizeof(uint32_t);
	char c;
	memcpy(&c, ptr, sizeof(char));
	m_bTightMBRs = (c != 0);
	ptr += sizeof(char);
	memcpy(&(m_stats.m_u32Nodes), ptr, sizeof(uint32_t));
	ptr += sizeof(uint32_t);
	memcpy(&(m_stats.m_u64Data), ptr, sizeof(uint64_t));
	ptr += sizeof(uint64_t);
	memcpy(&(m_stats.m_u32TreeHeight), ptr, sizeof(uint32_t));
	ptr += sizeof(uint32_t);

	for (uint32_t cLevel = 0; cLevel < m_stats.m_u32TreeHeight; ++cLevel)
	{
		uint32_t cNodes;
		memcpy(&cNodes, ptr, sizeof(uint32_t));
		ptr += sizeof(uint32_t);
		m_stats.m_nodesInLevel.push_back(cNodes);
	}

	delete[] header;
}

void SpatialIndex::RTree::RTree::insertData_impl(uint32_t dataLength, uint8_t* pData, Region& mbr, id_type id)
{
	assert(mbr.getDimension() == m_dimension);

	std::stack<id_type> pathBuffer;
	uint8_t* overflowTable = nullptr;

	try
	{
		NodePtr root = readNode(m_rootID);

		overflowTable = new uint8_t[root->m_level];
		memset(overflowTable, 0, root->m_level);

		NodePtr l = root->chooseSubtree(mbr, 0, pathBuffer);
		if (l.get() == root.get())
		{
			assert(root.unique());
			root.relinquish();
		}
		l->insertData(dataLength, pData, mbr, id, pathBuffer, overflowTable);

		delete[] overflowTable;
		++(m_stats.m_u64Data);
	}
	catch (...)
	{
		delete[] overflowTable;
		throw;
	}
}

void SpatialIndex::RTree::RTree::insertData_impl(uint32_t dataLength, uint8_t* pData, Region& mbr, id_type id, uint32_t level, uint8_t* overflowTable)
{
	assert(mbr.getDimension() == m_dimension);

	std::stack<id_type> pathBuffer;
	NodePtr root = readNode(m_rootID);
	NodePtr n = root->chooseSubtree(mbr, level, pathBuffer);

	assert(n->m_level == level);

	if (n.get() == root.get())
	{
		assert(root.unique());
		root.relinquish();
	}
	n->insertData(dataLength, pData, mbr, id, pathBuffer, overflowTable);
}

bool SpatialIndex::RTree::RTree::deleteData_impl(const Region& mbr, id_type id)
{
	assert(mbr.m_dimension == m_dimension);

	std::stack<id_type> pathBuffer;
	NodePtr root = readNode(m_rootID);
	NodePtr l = root->findLeaf(mbr, id, pathBuffer);
	if (l.get() == root.get())
	{
		assert(root.unique());
		root.relinquish();
	}

	if (l.get() != nullptr)
	{
		Leaf* pL = static_cast<Leaf*>(l.get());
		pL->deleteData(mbr, id, pathBuffer);
		--(m_stats.m_u64Data);
		return true;
	}

	return false;
}

SpatialIndex::id_type SpatialIndex::RTree::RTree::writeNode_complete(Node* n)
{
	uint8_t* buffer;
	uint32_t dataLength;
	n->storeToByteArray(&buffer, dataLength);

	id_type page;
	// if (n->m_identifier < 0) page = StorageManager::NewPage;
	// else page = n->m_identifier;

	page = StorageManager::NewPage;

	try
	{
		m_pStorageManager->storeByteArray(page, dataLength, buffer);
		delete[] buffer;
	}
	catch (InvalidPageException& e)
	{
		delete[] buffer;
		std::cerr << e.what() << std::endl;
		throw;
	}

	if (n->m_identifier < 0)
	{
		n->m_identifier = page;

		++(m_stats.m_u32Nodes);

#ifndef NDEBUG
		try
		{
			m_stats.m_nodesInLevel[n->m_level] = m_stats.m_nodesInLevel.at(n->m_level) + 1;
		}
		catch(...)
		{
			throw Tools::IllegalStateException("writeNode: writing past the end of m_nodesInLevel.");
		}
#else
		m_stats.m_nodesInLevel[n->m_level] = m_stats.m_nodesInLevel[n->m_level] + 1;
#endif
	}

	++(m_stats.m_u64Writes);

	for (size_t cIndex = 0; cIndex < m_writeNodeCommands.size(); ++cIndex)
	{
		m_writeNodeCommands[cIndex]->execute(*n);
	}

	return page;
}

SpatialIndex::id_type SpatialIndex::RTree::RTree::writeNode_real(Node* n)
{
	uint8_t* buffer;
	uint32_t dataLength;
	n->storeToByteArray(&buffer, dataLength);


	id_type page;
	// if (n->m_identifier < 0) page = StorageManager::NewPage;
	// else page = n->m_identifier;

	page = StorageManager::NewPage;

	try
	{

		m_pStorageManager->storeByteArray(page, dataLength, buffer);

		delete[] buffer;
	}
	catch (InvalidPageException& e)
	{
		delete[] buffer;
		std::cerr << e.what() << std::endl;
		throw;
	}

	n->m_identifier = page;
}

SpatialIndex::id_type SpatialIndex::RTree::RTree::writeNode(Node* n)
{
	// id_type page = -1;

	uint8_t* buffer;
	uint32_t dataLength;
	n->storeToByteArray(&buffer, dataLength);

	id_type page;
	// if (n->m_identifier < 0) page = StorageManager::NewPage;
	// else page = n->m_identifier;

	page = StorageManager::NewPage;

	try
	{
		m_pStorageManager->storeByteArray(page, dataLength, buffer);
		delete[] buffer;
	}
	catch (InvalidPageException& e)
	{
		delete[] buffer;
		std::cerr << e.what() << std::endl;
		throw;
	}

	if (n->m_identifier < 0)
	{
		n->m_identifier = page;

		++(m_stats.m_u32Nodes);

#ifndef NDEBUG
		try
		{
			m_stats.m_nodesInLevel[n->m_level] = m_stats.m_nodesInLevel.at(n->m_level) + 1;
		}
		catch(...)
		{
			throw Tools::IllegalStateException("writeNode: writing past the end of m_nodesInLevel.");
		}
#else
		m_stats.m_nodesInLevel[n->m_level] = m_stats.m_nodesInLevel[n->m_level] + 1;
#endif
	}

	++(m_stats.m_u64Writes);

	for (size_t cIndex = 0; cIndex < m_writeNodeCommands.size(); ++cIndex)
	{
		m_writeNodeCommands[cIndex]->execute(*n);
	}

	return page;
}



SpatialIndex::RTree::NodePtr SpatialIndex::RTree::RTree::readNode(id_type page)
{
	uint32_t dataLength;
	uint8_t* buffer;

	try
	{
		m_pStorageManager->loadByteArray(page, dataLength, &buffer);
	}
	catch (InvalidPageException& e)
	{
		std::cerr << e.what() << std::endl;
		throw;
	}

	try
	{
		uint32_t nodeType;
		memcpy(&nodeType, buffer, sizeof(uint32_t));

		NodePtr n;

		if (nodeType == PersistentIndex) n = m_indexPool.acquire();
		else if (nodeType == PersistentLeaf) n = m_leafPool.acquire();
		else throw Tools::IllegalStateException("readNode: failed reading the correct node type information");

		if (n.get() == nullptr)
		{
			if (nodeType == PersistentIndex) n = NodePtr(new Index(this, -1, 0), &m_indexPool);
			else if (nodeType == PersistentLeaf) n = NodePtr(new Leaf(this, -1), &m_leafPool);
		}

		//n->m_pTree = this;
		n->m_identifier = page;
		n->loadFromByteArray(buffer);

		++(m_stats.m_u64Reads);

		for (size_t cIndex = 0; cIndex < m_readNodeCommands.size(); ++cIndex)
		{
			m_readNodeCommands[cIndex]->execute(*n);
		}

		delete[] buffer;
		return n;
	}
	catch (...)
	{
		delete[] buffer;
		throw;
	}
}

// *=========================================================
// * In memory version
// *=========================================================
void SpatialIndex::RTree::RTree::readRawInMemData(id_type page, uint8_t ** addr) {

	uint32_t dataLength;
	m_pStorageManager->loadByteArray(page, dataLength, addr);
}


// *=========================================================
// * self-implemented readNode()
// * read node data and write in a given address
// *=========================================================
void SpatialIndex::RTree::RTree::readRawData(id_type page, uint8_t * addr)
{
	try
	{
		// m_pStorageManager->loadByteArray_nomalloc_cstyle(page, addr);
		m_pStorageManager->loadByteArray_nomalloc(page, addr);

		// uint32_t dataLength;
		// m_pStorageManager->loadByteArray(page, dataLength, &addr);
	}
	catch (InvalidPageException& e)
	{
		std::cerr << e.what() << std::endl;
		throw;
	}
}

void SpatialIndex::RTree::RTree::readNode(id_type page, uint8_t * addr, NodePtr & n)
{
	// uint32_t dataLength;
	// uint8_t* buffer;
	// try
	// {
	// 	m_pStorageManager->loadByteArray_nomalloc(page, addr);
	// 	// m_pStorageManager->loadByteArray(page, dataLength, &buffer);
	// }
	// catch (InvalidPageException& e)
	// {
	// 	std::cerr << e.what() << std::endl;
	// 	throw;
	// }

	try
	{
		uint32_t nodeType;
		memcpy(&nodeType, addr, sizeof(uint32_t));
		// memcpy(&nodeType, buffer, sizeof(uint32_t));

		if (nodeType == PersistentIndex) n = m_indexPool.acquire();
		else if (nodeType == PersistentLeaf) n = m_leafPool.acquire();
		else throw Tools::IllegalStateException("readNode: failed reading the correct node type information");

		if (n.get() == nullptr)
		{
			if (nodeType == PersistentIndex) n = NodePtr(new Index(this, -1, 0), &m_indexPool);
			else if (nodeType == PersistentLeaf) n = NodePtr(new Leaf(this, -1), &m_leafPool);
		}

		//n->m_pTree = this;
		n->m_identifier = page;
		n->loadFromByteArray(addr);
		// n->loadFromByteArray(buffer);

		++(m_stats.m_u64Reads);

		for (size_t cIndex = 0; cIndex < m_readNodeCommands.size(); ++cIndex)
		{
			m_readNodeCommands[cIndex]->execute(*n);
		}
	}
	catch (...)
	{
		throw;
	}
}

uint32_t SpatialIndex::RTree::RTree::readNodeLength(id_type page)
{
	return m_pStorageManager->loadByteArrayLength(page);
}


void SpatialIndex::RTree::RTree::deleteNode(Node* n)
{
	try
	{
		m_pStorageManager->deleteByteArray(n->m_identifier);
	}
	catch (InvalidPageException& e)
	{
		std::cerr << e.what() << std::endl;
		throw;
	}

	--(m_stats.m_u32Nodes);
	m_stats.m_nodesInLevel[n->m_level] = m_stats.m_nodesInLevel[n->m_level] - 1;

	for (size_t cIndex = 0; cIndex < m_deleteNodeCommands.size(); ++cIndex)
	{
		m_deleteNodeCommands[cIndex]->execute(*n);
	}
}

void SpatialIndex::RTree::RTree::rangeQuery(RangeQueryType type, const IShape& query, IVisitor& v)
{
	std::stack<NodePtr> st;
	NodePtr root = readNode(m_rootID);

	if (root->m_children > 0 && query.intersectsShape(root->m_nodeMBR)) st.push(root);

	while (! st.empty())
	{
		NodePtr n = st.top(); st.pop();

		if (n->m_level == 0)
		{
			v.visitNode(*n);

			for (uint32_t cChild = 0; cChild < n->m_children; ++cChild)
			{
				bool b;
				if (type == ContainmentQuery) b = query.containsShape(*(n->m_ptrMBR[cChild]));
				else b = query.intersectsShape(*(n->m_ptrMBR[cChild]));

				if (b)
				{
					Data data = Data(n->m_pDataLength[cChild], n->m_pData[cChild], *(n->m_ptrMBR[cChild]), n->m_pIdentifier[cChild]);
					v.visitData(data);
					++(m_stats.m_u64QueryResults);
				}
			}
		}
		else
		{
			v.visitNode(*n);

			for (uint32_t cChild = 0; cChild < n->m_children; ++cChild)
			{
				if (query.intersectsShape(*(n->m_ptrMBR[cChild]))) 
					st.push(readNode(n->m_pIdentifier[cChild]));
			}
		}
	}
}

void SpatialIndex::RTree::RTree::selfJoinQuery(id_type id1, id_type id2, const Region& r, IVisitor& vis)
{
	NodePtr n1 = readNode(id1);
	NodePtr n2 = readNode(id2);
	vis.visitNode(*n1);
	vis.visitNode(*n2);

	for (uint32_t cChild1 = 0; cChild1 < n1->m_children; ++cChild1)
	{
		if (r.intersectsRegion(*(n1->m_ptrMBR[cChild1])))
		{
			for (uint32_t cChild2 = 0; cChild2 < n2->m_children; ++cChild2)
			{
				if (r.intersectsRegion(*(n2->m_ptrMBR[cChild2])) &&
					n1->m_ptrMBR[cChild1]->intersectsRegion(*(n2->m_ptrMBR[cChild2])))
				{
					if (n1->m_level == 0)
					{
						if (n1->m_pIdentifier[cChild1] != n2->m_pIdentifier[cChild2])
						{
							assert(n2->m_level == 0);

							std::vector<const IData*> v;
							Data e1(n1->m_pDataLength[cChild1], n1->m_pData[cChild1], *(n1->m_ptrMBR[cChild1]), n1->m_pIdentifier[cChild1]);
							Data e2(n2->m_pDataLength[cChild2], n2->m_pData[cChild2], *(n2->m_ptrMBR[cChild2]), n2->m_pIdentifier[cChild2]);
							v.push_back(&e1);
							v.push_back(&e2);
							vis.visitData(v);
						}
					}
					else
					{
						Region rr = r.getIntersectingRegion(n1->m_ptrMBR[cChild1]->getIntersectingRegion(*(n2->m_ptrMBR[cChild2])));
						selfJoinQuery(n1->m_pIdentifier[cChild1], n2->m_pIdentifier[cChild2], rr, vis);
					}
				}
			}
		}
	}
}

/**
 * @brief 			spatial join 2 separate r-trees
 * @param Stree 	another r-tree 
 * @param id1 		id in this r-tree
 * @param id2 		id in Stree
 * @param r  		not used in this version
 * @param joinRes 	number of results
 */
void SpatialIndex::RTree::RTree::SpatialJoinQuery1(SpatialIndex::RTree::RTree *Stree, id_type id1, id_type id2, const Region& r,
													uint32_t &joinRes)
{
	NodePtr n1 = readNode(id1);
	NodePtr n2 = Stree->readNode(id2);

	// std::cout << n1->m_children << ", " << n2->m_children << std::endl;

	for (uint32_t cChild1 = 0; cChild1 < n1->m_children; ++cChild1)
	{
		for (uint32_t cChild2 = 0; cChild2 < n2->m_children; ++cChild2)
		{
			if (!n1->m_ptrMBR[cChild1]->intersectsRegion(*(n2->m_ptrMBR[cChild2])))
				continue;
			
			if (n1->m_level == 0 && n2->m_level == 0)
			{
				joinRes++;
			}
			else if (n1->m_level == 0 && n2->m_level > 0) // n2->m_level != 0
			{
				// SpatialJoinQuery1(Stree, n1->m_identifier, n2->m_pIdentifier[cChild2], r, joinRes);
			}
			else if (n1->m_level > 0 && n2->m_level == 0) // n1->m_level != 0
			{
				// SpatialJoinQuery1(Stree, n1->m_pIdentifier[cChild1], n2->m_identifier, r, joinRes);
			}
			else
			{
				SpatialJoinQuery1(Stree, n1->m_pIdentifier[cChild1], n2->m_pIdentifier[cChild2], r, joinRes);
			}
		}

		if (n1->m_level == 0 && n2->m_level > 0) // n2->m_level != 0
		{
			// SpatialJoinQuery1(Stree, n1->m_identifier, n2->m_pIdentifier[cChild2], r, joinRes);
		}
		else if (n1->m_level > 0 && n2->m_level == 0) // n1->m_level != 0
		{
			SpatialJoinQuery1(Stree, n1->m_pIdentifier[cChild1], n2->m_identifier, r, joinRes);
		}
	}
}

void SpatialIndex::RTree::RTree::SpatialJoinQuery1Iterative(SpatialIndex::RTree::RTree *Stree, id_type id1, id_type id2, const Region& r,
													uint32_t &joinRes)
{
	cur_level = 0;
	id1_lvl[cur_level] = id1;
	id2_lvl[cur_level] = id2;
	child1[cur_level] = 0;
	child2[cur_level] = 0;
	n1_lvl[cur_level] = readNode(id1_lvl[cur_level]);
	n2_lvl[cur_level] = Stree->readNode(id2_lvl[cur_level]);

	while (true)
	{
		while (child1[cur_level] < n1_lvl[cur_level]->m_children)
		{
			if (!n1_lvl[cur_level]->m_ptrMBR[child1[cur_level]]->intersectsRegion(*(n2_lvl[cur_level]->m_ptrMBR[child2[cur_level]])))
			{
				child2[cur_level]++;
				if (child2[cur_level] == n2_lvl[cur_level]->m_children)
				{
					child2[cur_level] = 0;
					child1[cur_level] ++;
				}
				continue;
			}
		
			if (n1_lvl[cur_level]->m_level == 0 && n2_lvl[cur_level]->m_level == 0)
			{
				joinRes++;
				child2[cur_level]++;
				if (child2[cur_level] == n2_lvl[cur_level]->m_children)
				{
					child2[cur_level] = 0;
					child1[cur_level] ++;
				}
			}
			else
			{
				id1_lvl[cur_level + 1] = n1_lvl[cur_level]->m_pIdentifier[child1[cur_level]];
				id2_lvl[cur_level + 1] = n2_lvl[cur_level]->m_pIdentifier[child2[cur_level]];
				child2[cur_level]++;
				if (child2[cur_level] == n2_lvl[cur_level]->m_children)
				{
					child2[cur_level] = 0;
					child1[cur_level] ++;
				}
				cur_level ++;
				child1[cur_level] = 0;
				child2[cur_level] = 0;
				n1_lvl[cur_level] = readNode(id1_lvl[cur_level]);
				n2_lvl[cur_level] = Stree->readNode(id2_lvl[cur_level]);
			}
		}
		cur_level--;
		if (cur_level < 0)
			break;
	}
}

/**
 * @brief 		spatial join 2 separate r-trees
 * @param Stree another r-tree 
 * @param id1 	id in this r-tree
 * @param r  	used
 * @param id2 	id in Stree
 */
// get rid of recursions.
void SpatialIndex::RTree::RTree::SpatialJoinQuery2(SpatialIndex::RTree::RTree *Stree, id_type id1, id_type id2, const Region& r,
													uint32_t &joinRes)
{
	if (id1 < 0) // slow
	{
		NodePtr n2 = Stree->readNode(id2);
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
					SpatialJoinQuery2(Stree, -1, n2->m_pIdentifier[cChild2], r, joinRes);
			}
		}
	}
	else if (id2 < 0) // fast
	{
		NodePtr n1 = this->readNode(id1);
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
					SpatialJoinQuery2(Stree, n1->m_pIdentifier[cChild1], -1, r, joinRes);
			}
		}
	}
	else
	{
		NodePtr n1 = readNode(id1);
		NodePtr n2 = Stree->readNode(id2);

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
						SpatialJoinQuery2(Stree, -1, n2->m_pIdentifier[cChild2], *n1->m_ptrMBR[cChild1], joinRes);
				}
				else if (n2->m_level == 0) // n1->m_level != 0
				{
					// Region rr = n1->m_ptrMBR[cChild1]->getIntersectingRegion(*(n2->m_ptrMBR[cChild2]));
					// if (rr.m_pLow[0] != std::numeric_limits<double>::max())
						SpatialJoinQuery2(Stree, n1->m_pIdentifier[cChild1], -1, *n2->m_ptrMBR[cChild2], joinRes);
				}
				else
				{
					// Region rr = n1->m_ptrMBR[cChild1]->getIntersectingRegion(*(n2->m_ptrMBR[cChild2]));
					// if (rr.m_pLow[0] != std::numeric_limits<double>::max())
						SpatialJoinQuery2(Stree, n1->m_pIdentifier[cChild1], n2->m_pIdentifier[cChild2], r, joinRes);
				}
			}
		}
	}
}

void SpatialIndex::RTree::RTree::SpatialJoinQueryWithSorting1(SpatialIndex::RTree::RTree *Stree, id_type id1, id_type id2, const Region& r,
													uint32_t &joinRes)
{
	if (id1 < 0)
	{
		NodePtr n2 = Stree->readNode(id2);
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
					outputs[n2->m_level].push_back({{RegionPtr{}, -1, 0},
													{n2->m_ptrMBR[k], n2->m_pIdentifier[k], k}});
			}
			k++;
		}

		if (n2->m_level != 0)
		{
			for (uint32_t k = 0; k < outputs[n2->m_level].size(); ++k)
			{
				Region rr = r.getIntersectingRegion(*(outputs[n2->m_level][k].second.regionPtr));
				// std::cout << outputs[n1->m_level][k].first.id << " " << outputs[n1->m_level][k].second.id << std::endl;
				SpatialJoinQueryWithSorting1(Stree, -1, outputs[n2->m_level][k].second.id, rr, joinRes);
			}
		}
		outputs[n2->m_level].resize(0);
	}
	else if (id2 < 0)
	{
		NodePtr n1 = readNode(id1);
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
					outputs[n1->m_level].push_back({{n1->m_ptrMBR[k], n1->m_pIdentifier[k], k},
													{RegionPtr{}, -1, 0}});
			}
			k++;
		}

		if (n1->m_level != 0)
		{
			for (uint32_t k = 0; k < outputs[n1->m_level].size(); ++k)
			{
				Region rr = r.getIntersectingRegion(*(outputs[n1->m_level][k].first.regionPtr));
				// std::cout << outputs[n1->m_level][k].first.id << " " << outputs[n1->m_level][k].second.id << std::endl;
				SpatialJoinQueryWithSorting1(Stree, outputs[n1->m_level][k].first.id, -1, rr, joinRes);
			}
		}
		outputs[n1->m_level].resize(0);
	}
	else 
	{
		NodePtr n1 = readNode(id1);
		NodePtr n2 = Stree->readNode(id2);

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
						else
							outputs[cur_level].push_back({{region, n1->m_pIdentifier[i], i},
															{n2->m_ptrMBR[k], n2->m_pIdentifier[k], k}});
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
						else
							outputs[cur_level].push_back({{n1->m_ptrMBR[k], n1->m_pIdentifier[k], k},
															{region, n2->m_pIdentifier[j], j}});
					}
					k++;
				}
				j++;
			}
		}
		if (n1->m_level != 0 && n2->m_level != 0)
		{
			for (uint32_t k = 0; k < outputs[cur_level].size(); ++k)
			{
				Region rr = outputs[cur_level][k].first.regionPtr->getIntersectingRegion(*(outputs[cur_level][k].second.regionPtr));
				// std::cout << outputs[n1->m_level][k].first.id << " " << outputs[n1->m_level][k].second.id << std::endl;
				// if (rr.m_pLow[0] != std::numeric_limits<double>::max())
				SpatialJoinQueryWithSorting1(Stree, outputs[cur_level][k].first.id, outputs[cur_level][k].second.id, rr, joinRes);
			}
		}
		else if (n1->m_level != 0)
		{
			for (uint32_t k = 0; k < outputs[cur_level].size(); ++k)
			{
				Region rr = outputs[cur_level][k].first.regionPtr->getIntersectingRegion(*(outputs[cur_level][k].second.regionPtr));
				// std::cout << outputs[n1->m_level][k].first.id << " " << outputs[n1->m_level][k].second.id << std::endl;
				// if (rr.m_pLow[0] != std::numeric_limits<double>::max())
				SpatialJoinQueryWithSorting1(Stree, outputs[cur_level][k].first.id, -1, rr, joinRes);
			}
		}
		else if (n2->m_level != 0)
		{
			for (uint32_t k = 0; k < outputs[cur_level].size(); ++k)
			{
				Region rr = outputs[cur_level][k].first.regionPtr->getIntersectingRegion(*(outputs[cur_level][k].second.regionPtr));
				// std::cout << outputs[n1->m_level][k].first.id << " " << outputs[n1->m_level][k].second.id << std::endl;
				// if (rr.m_pLow[0] != std::numeric_limits<double>::max())
				SpatialJoinQueryWithSorting1(Stree, -1, outputs[cur_level][k].second.id, rr, joinRes);
			}
		}
		outputs[cur_level].resize(0);
	}
}


/************************** not used anymore ***************************/
void SpatialIndex::RTree::RTree::InternalLoop(const SortingIndex &region, uint32_t k, 
											const std::vector<SortingIndex> &Seq, 
											std::vector<std::pair<SortingIndex, SortingIndex>> &outputs,
											bool Rfist)
{
	while (k < Seq.size() && Seq[k].regionPtr->m_pLow[0] <= region.regionPtr->m_pHigh[0])
	{
		if (region.regionPtr->m_pLow[1] <= Seq[k].regionPtr->m_pHigh[1] && 
			region.regionPtr->m_pHigh[1] >= Seq[k].regionPtr->m_pLow[1])
		{
			if (Rfist)
			{
				outputs.push_back({region, Seq[k]});
			}
			else 
			{
				outputs.push_back({Seq[k], region});
			}
		}
		k++;
	}
}

void SpatialIndex::RTree::RTree::SpatialJoinQueryWithSorting2(SpatialIndex::RTree::RTree *Stree, id_type id1, id_type id2, const Region& r,
													uint32_t &joinRes)
{
	if (id1 < 0)
	{
		NodePtr n2 = Stree->readNode(id2);
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
					SpatialJoinQueryWithSorting2(Stree, -1, n2->m_pIdentifier[k], rr, joinRes);
				}
			}
			k++;
		}
	}
	else if (id2 < 0)
	{
		NodePtr n1 = readNode(id1);
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
					SpatialJoinQueryWithSorting2(Stree, n1->m_pIdentifier[k], -1, rr, joinRes);
				}
			}
			k++;
		}
	}
	else 
	{
		NodePtr n1 = readNode(id1);
		NodePtr n2 = Stree->readNode(id2);

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
							SpatialJoinQueryWithSorting2(Stree, -1, n2->m_pIdentifier[k], rr, joinRes);
						}
						else if (n2->m_level == 0)
						{
							Region rr = region->getIntersectingRegion(*(n2->m_ptrMBR[k]));
							SpatialJoinQueryWithSorting2(Stree, n1->m_pIdentifier[i], -1, rr, joinRes);
						}
						else 
							SpatialJoinQueryWithSorting2(Stree, n1->m_pIdentifier[i], n2->m_pIdentifier[k], r, joinRes);
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
							SpatialJoinQueryWithSorting2(Stree, -1, n2->m_pIdentifier[j], rr, joinRes);
						}
						else if (n2->m_level == 0)
						{
							Region rr = region->getIntersectingRegion(*(n1->m_ptrMBR[k]));
							SpatialJoinQueryWithSorting2(Stree, n1->m_pIdentifier[k], -1, rr, joinRes);
						}
						else 
							SpatialJoinQueryWithSorting2(Stree, n1->m_pIdentifier[k], n2->m_pIdentifier[j], r, joinRes);
					}
					k++;
				}
				j++;
			}
		}
	}
}

void SpatialIndex::RTree::RTree::SpatialJoinQueryWithSortingPinning(SpatialIndex::RTree::RTree *Stree, id_type id1, id_type id2, const Region& r,
													uint32_t &joinRes)
{
	NodePtr n1 = readNode(id1);
	NodePtr n2 = Stree->readNode(id2);

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
					if (n1->m_level == 0)
					{
						assert(n2->m_level == 0);
						joinRes ++;
					}
					else
					{
						degree[n1->m_level][i][k] = true;
						outputs[n1->m_level].push_back({{region, n1->m_pIdentifier[i], i},
														{n2->m_ptrMBR[k], n2->m_pIdentifier[k], k}});
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
					if (n1->m_level == 0)
					{
						assert(n2->m_level == 0);
						joinRes ++;
					}
					else
					{
						degree[n1->m_level][k][j] = true;
						outputs[n1->m_level].push_back({{n1->m_ptrMBR[k], n1->m_pIdentifier[k], k},
														{region, n2->m_pIdentifier[j], j}});
					}	
				}
				k++;
			}
			j++;
		}
	}

	if (n1->m_level != 0)
	{
		// perform pinning
		uint32_t cursor = 0;
		while (cursor < outputs[n1->m_level].size())
		{
			const auto &inner = outputs[n1->m_level][cursor];
			if (!degree[n1->m_level][inner.first.child_id][inner.second.child_id]) 
			{
				cursor++;
				continue;
			}

			SpatialJoinQueryWithSortingPinning(Stree, inner.first.id, inner.second.id, r, joinRes);

			degree[n1->m_level][inner.first.child_id][inner.second.child_id] = false;

			// compute degree
			uint32_t r_cur_deg = 0;
			for (uint32_t i = 0; i < n2->m_children; ++i)
			{
				if (degree[n1->m_level][inner.first.child_id][i]) // TODO: SIMD???
					r_cur_deg++;
			}
			uint32_t s_cur_deg = 0;
			for (uint32_t i = 0; i < n1->m_children; ++i)
			{
				if (degree[n1->m_level][i][inner.second.child_id]) 
					s_cur_deg++;
			}
			if (r_cur_deg < s_cur_deg)
			{
				// pin inner.second
				for (uint32_t i = 0; i < n1->m_children; ++i)
				{
					if (!degree[n1->m_level][i][inner.second.child_id])
						continue;
					SpatialJoinQueryWithSortingPinning(Stree, n1->m_pIdentifier[i], 
														n2->m_pIdentifier[inner.second.child_id], r, joinRes);
					degree[n1->m_level][i][inner.second.child_id] = false;
				}
			}
			else 
			{
				// pin inner.first
				for (uint32_t i = 0; i < n2->m_children; ++i)
				{
					if (!degree[n1->m_level][inner.first.child_id][i])
						continue;
					SpatialJoinQueryWithSortingPinning(Stree, n1->m_pIdentifier[inner.first.child_id],
														n2->m_pIdentifier[i], r, joinRes);
					degree[n1->m_level][inner.first.child_id][i] = false;
				}
			}
			cursor++;
		}
	}
	outputs[n1->m_level].resize(0);
	// for (uint32_t i = 0; i < n1->m_children; ++i)
	// 	memset(degree[n1->m_level][i].get(), 0, m_indexCapacity);
}


void SpatialIndex::RTree::RTree::visitSubTree(NodePtr subTree, IVisitor& v)
{
	std::stack<NodePtr> st;
	st.push(subTree);

	while (! st.empty())
	{
		NodePtr n = st.top(); st.pop();
		v.visitNode(*n);

		if(n->m_level == 0)
		{
			for (uint32_t cChild = 0; cChild < n->m_children; ++cChild)
			{
				Data data = Data(n->m_pDataLength[cChild], n->m_pData[cChild], *(n->m_ptrMBR[cChild]), n->m_pIdentifier[cChild]);
				v.visitData(data);
				++(m_stats.m_u64QueryResults);
			}
		}
		else
		{
			for (uint32_t cChild = 0; cChild < n->m_children; ++cChild)
			{
				st.push(readNode(n->m_pIdentifier[cChild]));
			}
		}
	}
}

std::ostream& SpatialIndex::RTree::operator<<(std::ostream& os, const RTree& t)
{
	os	<< "Dimension: " << t.m_dimension << std::endl
		<< "Fill factor: " << t.m_fillFactor << std::endl
		<< "Index capacity: " << t.m_indexCapacity << std::endl
		<< "Leaf capacity: " << t.m_leafCapacity << std::endl
		<< "Tight MBRs: " << ((t.m_bTightMBRs) ? "enabled" : "disabled") << std::endl;

	if (t.m_treeVariant == RV_RSTAR)
	{
		os	<< "Near minimum overlap factor: " << t.m_nearMinimumOverlapFactor << std::endl
			<< "Reinsert factor: " << t.m_reinsertFactor << std::endl
			<< "Split distribution factor: " << t.m_splitDistributionFactor << std::endl;
	}

	if (t.m_stats.getNumberOfNodesInLevel(0) > 0)
		os	<< "Utilization: " << 100 * t.m_stats.getNumberOfData() / (t.m_stats.getNumberOfNodesInLevel(0) * t.m_leafCapacity) << "%" << std::endl
			<< t.m_stats;

	#ifndef NDEBUG
	os	<< "Leaf pool hits: " << t.m_leafPool.m_hits << std::endl
		<< "Leaf pool misses: " << t.m_leafPool.m_misses << std::endl
		<< "Index pool hits: " << t.m_indexPool.m_hits << std::endl
		<< "Index pool misses: " << t.m_indexPool.m_misses << std::endl
		<< "Region pool hits: " << t.m_regionPool.m_hits << std::endl
		<< "Region pool misses: " << t.m_regionPool.m_misses << std::endl
		<< "Point pool hits: " << t.m_pointPool.m_hits << std::endl
		<< "Point pool misses: " << t.m_pointPool.m_misses << std::endl;
	#endif

	return os;
}
