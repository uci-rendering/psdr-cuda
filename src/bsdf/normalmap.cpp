#include <psdr/core/bitmap.h>
#include <psdr/core/warp.h>
#include <psdr/core/sampler.h>

#include <psdr/bsdf/normalmap.h>

namespace psdr
{

NormalMap::NormalMap(const Bitmap3fD &n_map) : m_nmap(n_map) {}


SpectrumC NormalMap::eval(const IntersectionC &its, const Vector3fC &wo, MaskC active) const {
    return __eval<false>(its, wo, active);
}


SpectrumD NormalMap::eval(const IntersectionD &its, const Vector3fD &wo, MaskD active) const {
    return __eval<true>(its, wo, active);
}
template <bool ad>
Vector3f<ad> wt(Vector3f<ad> wp) {
    return normalize(Vector3f<ad>(-wp.x(), -wp.y(), 0.f));
}

template <bool ad>
Float<ad> pdot(Vector3f<ad> a, Vector3f<ad> b) {
    return max(Float<ad>(0.f), dot(a, b));
}

template <bool ad>
Float<ad> G1(Vector3f<ad> wp, Vector3f<ad> w) {
    return min(1.f,
        max(0.f, Frame<ad>::cos_theta(w)) * max(0.f, Frame<ad>::cos_theta(wp))
        / (pdot<ad>(w, wp) + pdot<ad>(w, wt<ad>(wp)) * Frame<ad>::sin_theta(wp))
    );
}

template <bool ad>
Float<ad> lambda_p(Vector3f<ad> wp, Vector3f<ad> wi) {
    Float<ad> i_dot_p = pdot<ad>(wp, wi);
    return i_dot_p / (i_dot_p + pdot<ad>(wt<ad>(wp), wi) * Frame<ad>::sin_theta(wp));
}


template <bool ad>
Spectrum<ad> NormalMap::__eval(const Intersection<ad> &its, const Vector3f<ad> &wo, Mask<ad> active) const {
    Float<ad> cos_theta_i = Frame<ad>::cos_theta(its.wi),
              cos_theta_o = Frame<ad>::cos_theta(wo);
    active &= cos_theta_i > 0.f && cos_theta_o > 0.f;

    Vector3f<ad> wp = normalize(fmadd(m_nmap.eval<ad>(its.uv), 2, -1.f));    
    Frame<ad> p_frame(wp, normalize(fnmadd(wp, dot(wp, its.dp_du), its.dp_du)));
    Intersection<ad> perturbed_its = its;
    perturbed_its.wi = p_frame.to_local(its.wi);

    Vector3f<ad> perturbed_wo = p_frame.to_local(wo);
 
    Float<ad> shadowing = G1<ad>(wp, wo);
    Float<ad> lambda_p_ = lambda_p<ad>(wp, its.wi);
    Vector3f<ad> wt_ = wt<ad>(wp);

    // i -> p -> o
    Spectrum<ad> value = m_bsdf->eval(perturbed_its, perturbed_wo, active) * lambda_p_ * shadowing;

    // i -> p -> t -> o
    // Vector3f<ad> wo_reflected = normalize(wo - 2.0f * dot(wo, wt_) * wt_);
    // Vector3f<ad> reflected_perturbed_wo = p_frame.to_local(wo_reflected);
    // Float<ad> notShadowedWpMirror = 1.f - G1<ad>(wp, wo_reflected);
    // value[dot(wo, wt_) > 0] += m_bsdf->eval(perturbed_its, reflected_perturbed_wo)  * (lambda_p_ * notShadowedWpMirror * shadowing);

    // i -> t -> p -> o
    Vector3f<ad> wi_reflected = normalize(its.wi - 2.0f * dot(its.wi, wt_) * wt_);
    Intersection<ad> reflected_perturbed_its = perturbed_its;
    reflected_perturbed_its.wi = p_frame.to_local(wi_reflected);
    value[dot(its.wi, wt_) > 0] += m_bsdf->eval(reflected_perturbed_its, perturbed_wo, active) * (1.f - lambda_p_) * shadowing;

    return value & active;
}


BSDFSampleC NormalMap::sample(const IntersectionC &its, const Vector3fC &sample, MaskC active) const {
    return __sample<false>(its, sample, active);
}


BSDFSampleD NormalMap::sample(const IntersectionD &its, const Vector3fD &sample, MaskD active) const {
    return __sample<true>(its, sample, active);
}

BSDFSampleDualC NormalMap::sampleDual(const IntersectionC &its, const Vector3fC &sample, MaskC active) const {
    PSDR_ASSERT(0);
    BSDFSampleDualC bs;
    return bs;
}


BSDFSampleDualD NormalMap::sampleDual(const IntersectionD &its, const Vector3fD &sample, MaskD active) const {
    PSDR_ASSERT(0);
    BSDFSampleDualD bs;
    return bs;
}

FloatC NormalMap::pdf(const IntersectionC &its, const Vector3fC &wo, MaskC active) const {
    return __pdf<false>(its, wo, active);
}


FloatD NormalMap::pdf(const IntersectionD &its, const Vector3fD &wo, MaskD active) const {
    return __pdf<true>(its, wo, active);
}


template <bool ad>
Float<ad> NormalMap::__pdf(const Intersection<ad> &its, const Vector3f<ad> &wo, Mask<ad> active) const {    
    Float<ad> cos_theta_i = Frame<ad>::cos_theta(its.wi),
              cos_theta_o = Frame<ad>::cos_theta(wo);
    active &= cos_theta_i > 0.f && cos_theta_o > 0.f;

    Vector3f<ad> wp = normalize(fmadd(m_nmap.eval<ad>(its.uv), 2, -1.f));    
    Frame<ad> p_frame(wp, normalize(fnmadd(wp, dot(wp, its.dp_du), its.dp_du)));

    Vector3f<ad> perturbed_wo = p_frame.to_local(wo);
    Float<ad> probability_wp = lambda_p<ad>(wp, its.wi);
    Vector3f<ad> wt_ = wt<ad>(wp);

    // i -> p -> o
    Intersection<ad> perturbed_its = its;
    perturbed_its.wi = p_frame.to_local(its.wi);
    // Float<ad> value = probability_wp * m_bsdf->pdf(perturbed_its, perturbed_wo, active);

    // i -> t -> p -> o
    Vector3f<ad> wi_reflected = normalize(its.wi - 2.0f * dot(its.wi, wt_) * wt_);
    Intersection<ad> reflected_perturbed_its = perturbed_its;
    reflected_perturbed_its.wi = p_frame.to_local(wi_reflected);
    Float<ad> value = probability_wp * m_bsdf->pdf(perturbed_its, perturbed_wo, active) + (1.f-probability_wp) * m_bsdf->pdf(reflected_perturbed_its, perturbed_wo, active);

    return detach(value) & active;
}

template <bool ad>
BSDFSample<ad> NormalMap::__sample(const Intersection<ad> &its, const Vector3f<ad> &sample, Mask<ad> active) const {
    Vector3f<ad> wp = normalize(fmadd(m_nmap.eval<ad>(its.uv), 2, -1.f));    
    Frame<ad> p_frame(wp, normalize(fnmadd(wp, dot(wp, its.dp_du), its.dp_du)));

    Intersection<ad> perturbed_its = its;
    perturbed_its.wi = p_frame.to_local(its.wi);
    
    Float<ad> probability_wp = lambda_p<ad>(wp, its.wi);
    Vector3f<ad> wt_ = wt<ad>(wp);
    Mask<ad> itpo_mask = sample.z() >= probability_wp;

    // i -> p -> o
    BSDFSample<ad> bs = m_bsdf->sample(perturbed_its, sample, active && !itpo_mask);
    Vector3f<ad> perturbed_wo = p_frame.to_world(bs.wo);

    // i -> t -> p -> o
    Vector3f<ad> wi_reflected = normalize(its.wi - 2.0f * dot(its.wi, wt_) * wt_);
    Intersection<ad> reflected_perturbed_its = its;
    reflected_perturbed_its.wi = p_frame.to_local(wi_reflected);
    BSDFSample<ad> bs_itpo = m_bsdf->sample(reflected_perturbed_its, sample, active && itpo_mask);
    
    bs.wo[itpo_mask] = bs_itpo.wo;

    Float<ad> pdf1 = m_bsdf->pdf(perturbed_its, bs.wo, active);
    Float<ad> pdf2 = m_bsdf->pdf(reflected_perturbed_its, bs.wo, active);

    bs.pdf = probability_wp * pdf1 + (1.f-probability_wp)*pdf2;
    bs.wo = p_frame.to_world(bs.wo);
    bs.is_valid = active && (bs.is_valid || bs_itpo.is_valid);

    return detach(bs);
}

} // namespace psdr
