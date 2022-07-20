#include <psdr/core/bitmap.h>
#include <psdr/core/warp.h>
#include <psdr/bsdf/ggx.h>
#include <psdr/bsdf/roughconductor.h>

namespace psdr
{

SpectrumC RoughConductor::eval(const IntersectionC& its, const Vector3fC& wo, MaskC active) const {
    return __eval<false>(its, wo, active);
}


SpectrumD RoughConductor::eval(const IntersectionD& its, const Vector3fD& wo, MaskD active) const {
    return __eval<true>(its, wo, active);
}


BSDFSampleC RoughConductor::sample(const IntersectionC& its, const Vector3fC& sample, MaskC active) const {
    return __sample<false>(its, sample, active);
}


BSDFSampleD RoughConductor::sample(const IntersectionD& its, const Vector3fD& sample, MaskD active) const {
    return __sample<true>(its, sample, active);
}

BSDFSampleDualC RoughConductor::sampleDual(const IntersectionC &its, const Vector3fC &sample, MaskC active) const {
    return __sampleDual<false>(its, sample, active);
}


BSDFSampleDualD RoughConductor::sampleDual(const IntersectionD &its, const Vector3fD &sample, MaskD active) const {
    return __sampleDual<true>(its, sample, active);
}

template <bool ad>
BSDFSampleDual<ad> RoughConductor::__sampleDual(const Intersection<ad>& its, const Vector3f<ad>& sample, Mask<ad> active) const {
    BSDFSampleDual<ad> bs;
    Float<ad> cos_theta_i = Frame<ad>::cos_theta(its.wi);
    Float<ad> alpha_u = m_alpha_u.eval<ad>(its.uv);
    Float<ad> alpha_v = m_alpha_v.eval<ad>(its.uv);
    GGXDistribution distr(alpha_u, alpha_v);
    Mask<ad> isDual = alpha_u < 0.05;

    auto m_pair = distr.sample<ad>(its.wi, sample);
    Vector3f<ad> m = std::get<0>(m_pair);
    Vector3f<ad> m_flip(-m.x(), -m.y(), m.z());

    bs.pdf1 = std::get<1>(m_pair);
    bs.pdf2 = distr.smith_g1<ad>(its.wi, m_flip) * abs(dot(its.wi, m_flip)) * distr.eval<ad>(m_flip) / abs(Frame<ad>::cos_theta(its.wi));

    bs.pdf1 = 0.5f*(bs.pdf1+bs.pdf2);
    bs.pdf2 = bs.pdf1;

    bs.wo1 = fmsub(Vector3f<ad>(m), 2.f * dot(its.wi, m), its.wi);
    bs.pdf1 /= (4.f * dot(bs.wo1, m));
    bs.is_valid1 = (cos_theta_i > 0.f && neq(bs.pdf1, 0.f) && Frame<ad>::cos_theta(bs.wo1) > 0.f) & active;

    bs.wo2 = fmsub(Vector3f<ad>(m_flip), 2.f * dot(its.wi, m_flip), its.wi);
    bs.pdf2 /= (4.f * dot(bs.wo2, m_flip));
    bs.is_valid2 = (cos_theta_i > 0.f && neq(bs.pdf2, 0.f) && Frame<ad>::cos_theta(bs.wo2) > 0.f) & active;
    bs.isDual = isDual;
    return detach(bs);

}

FloatC RoughConductor::pdf(const IntersectionC& its, const Vector3fC& wo, MaskC active) const {
    return __pdf<false>(its, wo, active);
}


FloatD RoughConductor::pdf(const IntersectionD& its, const Vector3fD& wo, MaskD active) const {
    return __pdf<true>(its, wo, active);
}


template <bool ad>
Spectrum<ad> RoughConductor::__eval(const Intersection<ad>& its, const Vector3f<ad>& wo, Mask<ad> active) const {
    Float<ad> cos_theta_i = Frame<ad>::cos_theta(its.wi),
              cos_theta_o = Frame<ad>::cos_theta(wo);
    active &= cos_theta_i > 0.f && cos_theta_o > 0.f;
    Float<ad> alpha_u = m_alpha_u.eval<ad>(its.uv);
    Float<ad> alpha_v = m_alpha_v.eval<ad>(its.uv);
    GGXDistribution m_distr(alpha_u, alpha_v);
    Vector3f<ad> H = normalize(wo + its.wi);
    Float<ad> D = m_distr.eval<ad>(H);
    active &= neq(D, 0.f);
    Float<ad> G = m_distr.G<ad>(its.wi, wo, H);
    Spectrum<ad> result = D * G / (4.f * Frame<ad>::cos_theta(its.wi));
    Spectrum<ad> F = fresnel<ad>(m_eta.eval<ad>(its.uv), m_k.eval<ad>(its.uv), dot(its.wi, H));
    Spectrum<ad> specular_reflectance = m_specular_reflectance.eval<ad>(its.uv);

    return (F * result * specular_reflectance) & active;
}


template <bool ad>
Float<ad> RoughConductor::__pdf(const Intersection<ad>& its, const Vector3f<ad>& wo, Mask<ad> active) const {
    Float<ad> cos_theta_i = Frame<ad>::cos_theta(its.wi),
              cos_theta_o = Frame<ad>::cos_theta(wo);

    Spectrum<ad> m = normalize(wo + its.wi);
    
    active &= cos_theta_i > 0.f && cos_theta_o > 0.f &&
              dot(its.wi, m) > 0.f && dot(wo, m) > 0.f;

    Float<ad> alpha_u = m_alpha_u.eval<ad>(its.uv);
    Float<ad> alpha_v = m_alpha_v.eval<ad>(its.uv);
    GGXDistribution distr(alpha_u, alpha_v);
    Float<ad> result = distr.eval<ad>(m) * distr.smith_g1<ad>(its.wi, m) /
                       (4.f * cos_theta_i);
    result = detach(result);   
    return result;
}


template <bool ad>
BSDFSample<ad> RoughConductor::__sample(const Intersection<ad>& its, const Vector3f<ad>& sample, Mask<ad> active) const {
    BSDFSample<ad> bs;
    Float<ad> cos_theta_i = Frame<ad>::cos_theta(its.wi);
    Float<ad> alpha_u = m_alpha_u.eval<ad>(its.uv);
    Float<ad> alpha_v = m_alpha_v.eval<ad>(its.uv);
    GGXDistribution distr(alpha_u, alpha_v);

    auto m_pair = distr.sample<ad>(its.wi, sample);
    Vector3f<ad> m = std::get<0>(m_pair);
    Float<ad> m_pdf = std::get<1>(m_pair);
    bs.wo = fmsub(Vector3f<ad>(m), 2.f * dot(its.wi, m), its.wi);
    bs.eta = 1.0f; // TODO: Fix later
    bs.pdf = m_pdf / (4.f * dot(bs.wo, m));
    bs.is_valid = (cos_theta_i > 0.f && neq(bs.pdf, 0.f) && Frame<ad>::cos_theta(bs.wo) > 0.f) & active;
    return detach(bs);
}

} // namespace psdr
