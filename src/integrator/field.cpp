#include <misc/Exception.h>
#include <psdr/core/ray.h>
#include <psdr/core/intersection.h>
#include <psdr/scene/scene.h>
#include <psdr/integrator/field.h>

namespace psdr
{

FieldExtractionIntegrator::FieldExtractionIntegrator(const char *field) : m_field(field) {
    PSDR_ASSERT_MSG(
        m_field == "silhouette" ||
        m_field == "position"   ||
        m_field == "depth"      ||
        m_field == "geoNormal"  ||
        m_field == "shNormal"   ||
        m_field == "uv",
        "Unsupported field: " + m_field
    );
}


SpectrumC FieldExtractionIntegrator::Li(const Scene &scene, Sampler &sampler, const RayC &ray, MaskC active) const {
    return __Li<false>(scene, ray, active);
}


SpectrumD FieldExtractionIntegrator::Li(const Scene &scene, Sampler &sampler, const RayD &ray, MaskD active) const {
    return __Li<true>(scene, ray, active);
}


template <bool ad>
Spectrum<ad> FieldExtractionIntegrator::__Li(const Scene &scene, const Ray<ad> &ray, Mask<ad> active) const {
    Intersection<ad> its = scene.ray_intersect<ad>(ray);
    Vector3f<ad> result;

    if ( m_field == "silhouette" ) {
        result = full<Spectrum<ad>>(1.f);
    } else if ( m_field == "position" ) {
        result = its.p;
    } else if ( m_field == "depth" ) {
        result = its.t;
    } else if ( m_field == "geoNormal" ) {
        result = its.n;
    } else if ( m_field == "shNormal" ) {
        result = its.sh_frame.n;
    } else if ( m_field == "uv" ) {
        result = concat(its.uv, 0.f);
    } else {
        PSDR_ASSERT(false);
    }
    return result & (active && its.is_valid());
}

} // namespace psdr
