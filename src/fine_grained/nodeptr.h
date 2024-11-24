#pragma once

#include <cstdint>
#include <cstring>
#include <vector>
#include <climits>

#include "meta.h"

typedef SpatialIndex::id_type id_type;

#define DIMENSION 2

struct SimpleRegion
{
    double low0, low1, high0, high1;

    SimpleRegion getIntersectingArea(SimpleRegion const& r) {
        assert (!(low0 > r.high0 || high0 < r.low0));
        SimpleRegion ret;
        ret.low0 = std::max(low0, r.low0);
		ret.high0 = std::min(high0, r.high0);
        ret.low1 = std::max(low1, r.low1);
		ret.high1 = std::min(high1, r.high1);
        return ret;
    }

    SimpleRegion getIntersectingArea(double rlow0,
                                    double rlow1,
                                    double rhigh0,
                                    double rhigh1)
    {
        assert (!(low0 > rhigh0 || high0 < rlow0));
        SimpleRegion ret;
        ret.low0 = std::max(low0, rlow0);
		ret.high0 = std::min(high0, rhigh0);
        ret.low1 = std::max(low1, rlow1);
		ret.high1 = std::min(high1, rhigh1);
        return ret;
    }

    void initialize() {
        low0 = std::numeric_limits<double>::max();
        low1 = std::numeric_limits<double>::max();
        high0 = std::numeric_limits<double>::min();
        high1 = std::numeric_limits<double>::min();
    }
};

struct SimpleNodePtr
{
    constexpr static uint32_t PER_CHILD_DATA_LEN = sizeof(SimpleRegion) + sizeof(id_type) + sizeof(uint32_t);               

    uint32_t nodeType;
    id_type m_identifier{-1};
    uint32_t m_level{0};
    uint32_t m_children{0};

    // id_type* m_pIdentifier {nullptr};
    // SimpleRegion* m_ptrMBR {nullptr};
    // SimpleRegion m_nodeMBR;

    uint8_t *m_rawData {nullptr};
    uint8_t *m_childData {nullptr};

    SimpleNodePtr() {}

    void fill_data(id_type page, uint8_t * rawData)
    {    
        m_rawData = rawData;
        uint8_t * ptr = rawData;
        m_identifier = page;

        memcpy(&nodeType, ptr, sizeof(uint32_t));
        ptr += sizeof(uint32_t);

        memcpy(&m_level, ptr, sizeof(uint32_t));
        ptr += sizeof(uint32_t);

        memcpy(&m_children, ptr, sizeof(uint32_t));
        ptr += sizeof(uint32_t);

        m_childData = ptr;
    }
    SimpleRegion returnChildRegion(uint32_t childIndex)
    {
        assert(childIndex < m_children);
        uint8_t * ptr = m_childData + childIndex * PER_CHILD_DATA_LEN;

        SimpleRegion region;
        memcpy(&region.low0, ptr, sizeof(double));
        ptr += sizeof(double);
        memcpy(&region.low1, ptr, sizeof(double));
        ptr += sizeof(double);
        memcpy(&region.high0, ptr, sizeof(double));
        ptr += sizeof(double);
        memcpy(&region.high1, ptr, sizeof(double));
        return region;
    }
    id_type returnChildIdentifier(uint32_t childIndex)
    {
        assert(childIndex < m_children);
        uint8_t * ptr = m_childData + childIndex * PER_CHILD_DATA_LEN + sizeof(SimpleRegion);
        id_type identifier;
        memcpy(&identifier, ptr, sizeof(id_type));
        return identifier;
    }
};