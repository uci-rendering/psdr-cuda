#pragma once

#include "frame.h"

namespace psdr
{

template <typename Float_>
struct Interaction_ {
    static constexpr bool ad = std::is_same_v<Float_, FloatD>;

    template <typename Float1_>
    inline Interaction_(const Interaction_<Float1_> &in) : wi(in.wi), p(in.p), t(in.t) {}

    virtual Mask<ad> is_valid() const = 0;

    Vector3f<ad> wi, p;
    Float<ad> t = Infinity;

    ENOKI_STRUCT(Interaction_, wi, p, t)
};


template <typename Float_>
struct Intersection_ : public Interaction_<Float_> {
    PSDR_IMPORT_BASE(Interaction_<Float_>, ad, wi, p, t)

    Mask<ad> is_valid() const override {
        return neq(shape, nullptr);
    }

    inline Mask<ad> is_emitter(Mask<ad> active) const {
        return neq(shape->emitter(active), nullptr);
    }

    inline Spectrum<ad> Le(Mask<ad> active) const {
        return shape->emitter(active)->eval(*this, active);
    }

    MeshArray<ad>       shape;

    Vector3f<ad>        n;                  // geometric normal
    Frame<ad>           sh_frame;           // shading frame

    Vector2f<ad>        uv;
    Float<ad>           J;                  // Jacobian determinant for material-form reparam

    Int<ad>             v0_idx,
                        v1_idx,
                        v2_idx;             // indices of the three vertices into the mesh's vertex array
    Vector2f<ad>        barycentric_uv;     // barycentric coordinates of the intersection point

    ENOKI_DERIVED_STRUCT(Intersection_, Interaction_<Float_>,
        ENOKI_BASE_FIELDS(wi, p, t),
        ENOKI_DERIVED_FIELDS(shape, n, sh_frame, uv, J, v0_idx, v1_idx, v2_idx, barycentric_uv)
    )
};


//inline IntersectionC detach(const IntersectionD &its) {
//    IntersectionC result;
//    result.wi       = detach(its.wi);
//    result.p        = detach(its.p);
//    result.t        = detach(its.t);
//    result.shape    = detach(its.shape);
//    result.n        = detach(its.n);
//    result.sh_frame = detach(its.sh_frame);
//    result.uv       = detach(its.uv);
//    result.J        = detach(its.J);
//    return result;
//}

} // namespace psdr

ENOKI_STRUCT_SUPPORT(psdr::Interaction_, wi, p, t)
ENOKI_STRUCT_SUPPORT(psdr::Intersection_, wi, p, t, shape, n, sh_frame, uv, J, v0_idx, v1_idx, v2_idx, barycentric_uv)
