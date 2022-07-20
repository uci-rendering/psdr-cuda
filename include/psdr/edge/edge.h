#pragma once

#include <psdr/psdr.h>
#include <psdr/core/ray.h>
#include <psdr/core/pmf.h>

namespace psdr
{

/********************************************
* Primary edge info types
********************************************/

struct PrimaryEdgeSample {
    FloatD  x_dot_n;
    IntC    idx;
    RayC    ray_n, ray_p;

#ifdef PSDR_PRIMARY_EDGE_VIS_CHECK
    RayC    ray_c;
#endif

    FloatC  pdf;
};


template <typename Float_>
struct PrimaryEdgeInfo_ {
    static constexpr bool ad = std::is_same_v<Float_, FloatD>;

#ifdef PSDR_PRIMARY_EDGE_VIS_CHECK
    Vector3f<ad>            p0, p1;
#else
    Vector2f<ad>            p0, p1;
#endif
    Vector2f<ad>            edge_normal;
    Float<ad>               edge_length;

    ENOKI_STRUCT(PrimaryEdgeInfo_, p0, p1, edge_normal, edge_length)
};

using PrimaryEdgeInfo = PrimaryEdgeInfo_<FloatD>;


/********************************************
 * Secondary edge info
 ********************************************/

template <typename Float_>
struct SecondaryEdgeInfo_ {
    static constexpr bool ad = std::is_same_v<Float_, FloatD>;

    // p0 and (p0 + e1) are the two endpoints of the edge
    Vector3f<ad>            p0, e1;

    // n0 and n1 are the normals of the two faces sharing the edge
    Vector3f<ad>            n0, n1;

    // p2 is the third vertex of the face with normal n0
    Vector3f<ad>            p2;

    Mask<ad>                is_boundary;

    ENOKI_STRUCT(SecondaryEdgeInfo_, p0, e1, n0, n1, p2, is_boundary)
};

using SecondaryEdgeInfo = SecondaryEdgeInfo_<FloatD>;

} // namespace psdr

ENOKI_STRUCT_SUPPORT(psdr::PrimaryEdgeInfo_, p0, p1, edge_normal, edge_length)
ENOKI_STRUCT_SUPPORT(psdr::SecondaryEdgeInfo_, p0, e1, n0, n1, p2, is_boundary)

