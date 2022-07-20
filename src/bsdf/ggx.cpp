#include <misc/Exception.h>
#include <psdr/core/warp.h>
#include <psdr/bsdf/ggx.h>

namespace psdr
{

template <bool ad>
Float<ad> GGXDistribution::G(const Vector3f<ad>& wi, const Vector3f<ad>& wo, const Vector3f<ad>& m) const {
    return smith_g1<ad>(wi, m) * smith_g1<ad>(wo, m);
}


template <bool ad>
Float<ad> GGXDistribution::eval(const Vector3f<ad>& m) const {
    Float<ad> alpha_uv;
    if constexpr (ad) {
        alpha_uv = m_alpha_u * m_alpha_v;
    }
    else {
        alpha_uv = detach(m_alpha_u) * detach(m_alpha_v);
    }
    Float<ad> cos_theta = Frame<ad>::cos_theta(m),
        cos_theta_2 = sqr(cos_theta),
        result;
    if constexpr (ad) {
        result = rcp(Pi * alpha_uv * sqr(sqr(m.x() / m_alpha_u) + sqr(m.y() / m_alpha_v) + sqr(m.z())));
    }
    else {
        result = rcp(Pi * alpha_uv * sqr(sqr(m.x() / detach(m_alpha_u)) + sqr(m.y() / detach(m_alpha_v)) + sqr(m.z())));
    }
    return select(result * cos_theta > 1e-5f, result, 0.f);
}


template <bool ad>
std::pair<Vector3f<ad>, Float<ad>> GGXDistribution::sample(const Vector3f<ad>& wi, const Vector3f<ad>& sample) const {
    Float<ad> sin_phi, cos_phi, cos_theta;
    Vector3f<ad> wi_p;

    if constexpr (ad) {
        wi_p = normalize(Vector3f<ad>(
            m_alpha_u * wi.x(),
            m_alpha_v * wi.y(),
            wi.z()
        ));
    }
    else {
        wi_p = normalize(Vector3f<ad>(
            detach(m_alpha_u) * wi.x(),
            detach(m_alpha_v) * wi.y(),
            wi.z()
        ));
    }

    sin_phi = Frame<ad>::sin_phi(wi_p);
    cos_phi = Frame<ad>::cos_phi(wi_p);
    cos_theta = Frame<ad>::cos_theta(wi_p);
    Vector2f<ad> sample2(sample.x(), sample.y());
    Vector2f<ad> slope = sample_visible_11<ad>(cos_theta, sample2);

    if constexpr (ad) {
        slope = Vector2f<ad>(
            fmsub(cos_phi, slope.x(), sin_phi * slope.y()) * m_alpha_u,
            fmadd(sin_phi, slope.x(), cos_phi * slope.y()) * m_alpha_v
        );
    }
    else {
        slope = Vector2f<ad>(
            fmsub(cos_phi, slope.x(), sin_phi * slope.y()) * detach(m_alpha_u),
            fmadd(sin_phi, slope.x(), cos_phi * slope.y()) * detach(m_alpha_v)
        );
    }
    Vector3f<ad> m = normalize(Vector3f<ad>(-slope.x(), -slope.y(), 1));

    // Compute probability density of the sampled position
    Float<ad> pdf = smith_g1<ad>(wi, m) * abs(dot(wi, m)) * eval<ad>(m) / abs(Frame<ad>::cos_theta(wi));
    pdf = detach(pdf);
    return {m, pdf};
}


template <bool ad>
Float<ad> GGXDistribution::smith_g1(const Vector3f<ad>& v, const Vector3f<ad>& m) const {
    Float<ad> xy_alpha_2;
    if constexpr (ad) {
        xy_alpha_2 = sqr(m_alpha_u * v.x()) + sqr(m_alpha_v * v.y());
    }
    else {
        xy_alpha_2 = sqr(detach(m_alpha_u) * v.x()) + sqr(detach(m_alpha_v) * v.y());
    }
    Float<ad> tan_theta_alpha_2 = xy_alpha_2 / sqr(v.z()), result;
    result = 2.f / (1.f + sqrt(1.f + tan_theta_alpha_2));
    masked(result, eq(xy_alpha_2, 0.f)) = 1.f;
    masked(result, dot(v, m) * Frame<ad>::cos_theta(v) <= 0.f) = 0.f;
    return result;
}


template <bool ad>
Vector2f<ad> GGXDistribution::sample_visible_11(const Float<ad>& cos_theta_i, const Vector2f<ad>& sample) const {
    Vector2f<ad> p = warp::square_to_uniform_disk_concentric<ad>(sample);
    Float<ad> s = .5f * (1.f + cos_theta_i);
    p.y() = lerp(safe_sqrt(1.f - sqr(p.x())), p.y(), s);
    Float<ad> x = p.x(), y = p.y(),
              z = safe_sqrt(1.f - squared_norm(p));
    Float<ad> sin_theta_i = safe_sqrt(1.f - sqr(cos_theta_i));
    Float<ad> norm = rcp(fmadd(sin_theta_i, y, cos_theta_i * z));
    return Vector2f<ad>(fmsub(cos_theta_i, y, sin_theta_i * z), x) * norm;
}

// Explicit instancitations

template FloatC GGXDistribution::eval<false>(const Vector3fC& m) const;
template FloatD GGXDistribution::eval<true >(const Vector3fD& m) const;
template FloatC GGXDistribution::G<false>(const Vector3fC& wi, const Vector3fC& wo, const Vector3fC& m) const;
template FloatD GGXDistribution::G<true >(const Vector3fD& wi, const Vector3fD& wo, const Vector3fD& m) const;
template std::pair<Vector3fC, FloatC> GGXDistribution::sample<false>(const Vector3fC& wi, const Vector3fC& sample) const;
template std::pair<Vector3fD, FloatD> GGXDistribution::sample<true >(const Vector3fD& wi, const Vector3fD& sample) const;
template FloatC GGXDistribution::smith_g1<false>(const Vector3fC& v, const Vector3fC& m) const;
template FloatD GGXDistribution::smith_g1<true >(const Vector3fD& v, const Vector3fD& m) const;

} // namespace psdr
