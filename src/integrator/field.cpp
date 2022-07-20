#include <misc/Exception.h>
#include <psdr/core/ray.h>
#include <psdr/core/intersection.h>
#include <psdr/scene/scene.h>
#include <psdr/integrator/field.h>
#include <psdr/bsdf/bsdf.h>
#include <psdr/shape/mesh.h>

namespace psdr
{

FieldExtractionIntegrator::FieldExtractionIntegrator(char *field) {
    std::string field_name = strtok(field, " ");
    PSDR_ASSERT_MSG(
        field_name == "segmentation" ||
        field_name == "silhouette" ||
        field_name == "position"   ||
        field_name == "depth"      ||
        field_name == "geoNormal"  ||
        field_name == "shNormal"   ||
        field_name == "uv",
        "Unsupported field: " + field_name
    );
    m_field = field_name;

   char *obj_name;
   obj_name = strtok(NULL, " ");
   if( obj_name != NULL ) {
      m_object = obj_name;
   } else {
      m_object = "";
   }

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

    BSDFArray<ad> bsdf_array = its.shape->bsdf(active);
    if ( scene.m_emitter_env != nullptr ) {
        // Skip reflectance computations for intersections on the bounding mesh
        active &= neq(bsdf_array, nullptr);
    }
    Mask<ad> valid_obj(1);
    if (m_object != "") {
        if constexpr ( !ad ) { 
            valid_obj = its.shape->get_obj_mask(m_object);
        } else {
            valid_obj = detach(its.shape)->get_obj_mask(m_object);
        }
    }

    if ( m_field == "segmentation" ) {
        if constexpr ( !ad ) { 
            result = its.shape->get_obj_id();
        } else {
            result = detach(its.shape)->get_obj_id();
        }
    } else if ( m_field == "silhouette" ) {
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
    return result & (active && (its.is_valid()) && valid_obj);
}

} // namespace psdr
