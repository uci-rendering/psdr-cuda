#include <misc/Exception.h>
#include <psdr/core/frame.h>
#include <psdr/core/intersection.h>
#include <psdr/core/transform.h>
#include <psdr/emitter/envmap.h>

namespace psdr
{

void EnvironmentMap::configure() {
    int width = m_radiance.m_resolution.x(), height = m_radiance.m_resolution.y();
    PSDR_ASSERT(width > 1 && height > 1);

    width = (width - 1) << 1; height = (height - 1) << 1;
    m_cell_distrb.set_resolution(ScalarVector2i(width, height));

    Vector2fC uv = (m_cell_distrb.m_cells + Vector2fC(.5f, .5f))*m_cell_distrb.m_unit;
    SpectrumC val = m_radiance.eval<false>(uv, false);

    FloatC theta = ((arange<IntC>(width*height) % height) + .5f)*(Pi/static_cast<float>(height));  
    m_cell_distrb.set_mass(rgb2luminance<false>(val)*sin(theta));

    m_to_world = m_to_world_left*m_to_world_raw;
    m_from_world = inverse(m_to_world);
    m_ready = true;
}


SpectrumC EnvironmentMap::eval(const IntersectionC &its, MaskC active) const {
    Vector3fC wi_world = its.sh_frame.to_world(its.wi);
    return eval_direction<false>(-wi_world, active);
}


SpectrumD EnvironmentMap::eval(const IntersectionD &its, MaskD active) const {
    Vector3fD wi_world = its.sh_frame.to_world(its.wi);
    return eval_direction<true>(-wi_world, active);
}


template <bool ad>
Spectrum<ad> EnvironmentMap::eval_direction(const Vector3f<ad> &wi, Mask<ad> active) const {
    PSDR_ASSERT(m_ready);
    Vector3f<ad> v;
    if constexpr ( ad ) {
        v = transform_dir<FloatD>(m_from_world, wi);
    } else {
        v = transform_dir<FloatC>(detach(m_from_world), wi);
    }

    Vector2f<ad> uv(atan2(v.x(), -v.z())*InvTwoPi, safe_acos(v.y())*InvPi);
    uv -= floor(uv);
    if constexpr ( ad ) {
        return m_radiance.eval<true>(uv, false)*m_scale;
    } else {
        return m_radiance.eval<false>(uv, false)*detach(m_scale);
    }
}


PositionSampleC EnvironmentMap::sample_position(const Vector3fC &ref_p, const Vector2fC &sample2, MaskC active) const {
    return __sample_position<false>(ref_p, sample2, active);
}


PositionSampleD EnvironmentMap::sample_position(const Vector3fD &ref_p, const Vector2fD &sample2, MaskD active) const {
    return __sample_position<true>(ref_p, sample2, active);
}


template <bool ad>
PositionSample<ad> EnvironmentMap::__sample_position(const Vector3f<ad> &ref_p, const Vector2f<ad> &_sample2, Mask<ad> active) const {
    PSDR_ASSERT(m_ready);
    PositionSample<ad> result;

    RayC ray;
    Vector2fC sample2;
    FloatC pdf;
    if constexpr ( ad ) {
        ray.o = detach(ref_p);
        sample2 = detach(_sample2);
    } else {
        ray.o = ref_p;
        sample2 = _sample2;
    }
    std::tie(ray.d, pdf) = sample_direction(sample2);
    auto [t, n, G] = ray_intersect_scene_aabb<false>(ray, m_lower, m_upper);

    result.is_valid = active;
    result.p = ray(t);
    result.n = n;
    result.pdf = pdf*G;
    result.J = 1.f;
    return result;
}


std::pair<Vector3fC, FloatC> EnvironmentMap::sample_direction(Vector2fC &uv) const {
    PSDR_ASSERT(m_ready);
    FloatC pdf = m_cell_distrb.sample_reuse(uv);

    FloatC theta = uv.y()*Pi, phi = uv.x()*TwoPi;
    Vector3fC d = sphdir<false>(theta, phi);
    d = Vector3fC(d.y(), d.z(), -d.x());

    FloatC inv_sin_theta = safe_rsqrt(max(sqr(d.x()) + sqr(d.z()), sqr(Epsilon)));
    masked(pdf, pdf > Epsilon) *= inv_sin_theta*(.5f/sqr(Pi));

    d = transform_dir<FloatC>(detach(m_to_world), d);
    return { d, pdf };
}


FloatC EnvironmentMap::sample_position_pdf(const Vector3fC &ref_p, const IntersectionC &its, MaskC active) const {
    return __sample_position_pdf<false>(ref_p, its, active);
}


FloatD EnvironmentMap::sample_position_pdf(const Vector3fD &ref_p, const IntersectionD &its, MaskD active) const {
    return __sample_position_pdf<true>(ref_p, its, active);
}


template <bool ad>
Float<ad> EnvironmentMap::__sample_position_pdf(const Vector3f<ad> &ref_p, const Intersection<ad> &its, Mask<ad> active) const {
    Vector3fC d;
    FloatC G, dist2;
    if constexpr ( ad ) {
        d = detach(its.p) - detach(ref_p);
        dist2 = squared_norm(d); d /= safe_sqrt(dist2);
        G = abs(dot(d, detach(its.n)))/dist2;
    } else {
        d = its.p - ref_p;
        dist2 = squared_norm(d); d /= safe_sqrt(dist2);
        G = abs(dot(d, its.n))/dist2;
    }

    d = transform_dir<FloatC>(detach(m_from_world), d);
    FloatC factor = G*safe_rsqrt(max(sqr(d.x()) + sqr(d.z()), sqr(Epsilon)))*(.5f/sqr(Pi));
    Vector2fC uv(atan2(d.x(), -d.z())*InvTwoPi, safe_acos(d.y())*InvPi);
    uv -= floor(uv);
    return m_cell_distrb.pdf(uv)*factor;
}


std::string EnvironmentMap::to_string() const {
    std::ostringstream oss;
    oss << "EnvironmentMap[sampling_weight = " << m_sampling_weight << "]";
    return oss.str();
}

// Explicit instantiations
template SpectrumC EnvironmentMap::eval_direction<false>(const Vector3fC&, MaskC) const;
template SpectrumD EnvironmentMap::eval_direction<true >(const Vector3fD&, MaskD) const;

} // namespace psdr
