#include <psdr/core/bitmap.h>
#include <psdr/core/warp.h>
#include <psdr/bsdf/diffuse.h>

namespace psdr
{

Diffuse::Diffuse(const char *refl_file) : m_reflectance(refl_file) {}


Diffuse::Diffuse(const Bitmap3fD &reflectance) : m_reflectance(reflectance) {}


SpectrumC Diffuse::eval(const IntersectionC &its, const Vector3fC &wo, MaskC active) const {
    return __eval<false>(its, wo, active);
}


SpectrumD Diffuse::eval(const IntersectionD &its, const Vector3fD &wo, MaskD active) const {
    return __eval<true>(its, wo, active);
}


template <bool ad>
Spectrum<ad> Diffuse::__eval(const Intersection<ad> &its, const Vector3f<ad> &wo, Mask<ad> active) const {
    Float<ad> cos_theta_i = Frame<ad>::cos_theta(its.wi),
              cos_theta_o = Frame<ad>::cos_theta(wo);

    active &= (cos_theta_i > 0.f && cos_theta_o > 0.f);

    Spectrum<ad> value = m_reflectance.eval<ad>(its.uv) * InvPi * cos_theta_o;
    return value & active;
}


BSDFSampleC Diffuse::sample(const IntersectionC &its, const Vector3fC &sample, MaskC active) const {
    return __sample<false>(its, sample, active);
}


BSDFSampleD Diffuse::sample(const IntersectionD &its, const Vector3fD &sample, MaskD active) const {
    return __sample<true>(its, sample, active);
}

BSDFSampleDualC Diffuse::sampleDual(const IntersectionC &its, const Vector3fC &sample, MaskC active) const {
    PSDR_ASSERT(0);
    BSDFSampleDualC bs;
    return bs;
}


BSDFSampleDualD Diffuse::sampleDual(const IntersectionD &its, const Vector3fD &sample, MaskD active) const {
    PSDR_ASSERT(0);
    BSDFSampleDualD bs;
    return bs;
}

template <bool ad>
BSDFSample<ad> Diffuse::__sample(const Intersection<ad> &its, const Vector3f<ad> &sample, Mask<ad> active) const {
    Float<ad> cos_theta_i = Frame<ad>::cos_theta(its.wi);
    BSDFSample<ad> bs;
    bs.wo = warp::square_to_cosine_hemisphere<ad>(tail<2>(sample));
    bs.eta = 1.0f;
    bs.pdf = warp::square_to_cosine_hemisphere_pdf<ad>(bs.wo);
    bs.is_valid = active && (cos_theta_i > 0.f);
    return detach(bs);
}


FloatC Diffuse::pdf(const IntersectionC &its, const Vector3fC &wo, MaskC active) const {
    return __pdf<false>(its, wo, active);
}


FloatD Diffuse::pdf(const IntersectionD &its, const Vector3fD &wo, MaskD active) const {
    return __pdf<true>(its, wo, active);
}


template <bool ad>
Float<ad> Diffuse::__pdf(const Intersection<ad> &its, const Vector3f<ad> &wo, Mask<ad> active) const {
    FloatC cos_theta_i, cos_theta_o;
    if constexpr ( ad ) {
        cos_theta_i = FrameC::cos_theta(detach(its.wi));
        cos_theta_o = FrameC::cos_theta(detach(wo));
    } else {
        cos_theta_i = FrameC::cos_theta(its.wi);
        cos_theta_o = FrameC::cos_theta(wo);
    }
    active &= (cos_theta_i > 0.f && cos_theta_o > 0.f);

    Float<ad> value = InvPi * cos_theta_o;
    return value & active;
}

} // namespace psdr
