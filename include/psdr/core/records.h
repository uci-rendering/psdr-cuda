#pragma once

#include <enoki/array_macro.h>
#include "ray.h"
#include "intersection.h"

namespace psdr
{

template <typename Float_>
struct SampleRecord_ {
    static constexpr bool ad = std::is_same_v<Float_, FloatD>;

    Float<ad>       pdf;
    Mask<ad>        is_valid;

    ENOKI_STRUCT(SampleRecord_, pdf, is_valid)
};


template <typename Float_>
struct PositionSample_ : public SampleRecord_<Float_> {
    PSDR_IMPORT_BASE(SampleRecord_<Float_>, ad, pdf, is_valid)

    Vector3f<ad>    p, n;
    Float<ad>       J;

    ENOKI_DERIVED_STRUCT(PositionSample_, Base,
        ENOKI_BASE_FIELDS(pdf, is_valid),
        ENOKI_DERIVED_FIELDS(p, n, J)
    )
};


struct BoundarySegSampleDirect : public SampleRecord_<FloatC> {
    PSDR_IMPORT_BASE(SampleRecord_<FloatC>, pdf, is_valid)
    
    // Sample point on a face edge
    Vector3fD           p0;
    Vector3fC           edge, edge2;

    // Sample point on an emitter
    Vector3fC           p2, n;
};

} // namespace psdr

ENOKI_STRUCT_SUPPORT(psdr::SampleRecord_, pdf, is_valid)
ENOKI_STRUCT_SUPPORT(psdr::PositionSample_, pdf, is_valid, p, n, J)
