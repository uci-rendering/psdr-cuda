#include <misc/Exception.h>
#include <psdr/core/cube_distrb.h>
#include <psdr/core/ray.h>
#include <psdr/core/intersection.h>
#include <psdr/core/sampler.h>
#include <psdr/core/transform.h>
#include <psdr/bsdf/bsdf.h>
#include <psdr/emitter/emitter.h>
#include <psdr/shape/mesh.h>
#include <psdr/scene/scene.h>
#include <psdr/sensor/perspective.h>
#include <psdr/integrator/direct.h>

namespace psdr
{

template <bool ad>
static inline Float<ad> mis_weight(const Float<ad> &pdf1, const Float<ad> &pdf2) {
    Float<ad> w1 = sqr(pdf1), w2 = sqr(pdf2);
    return w1/(w1 + w2);
}


DirectIntegrator::~DirectIntegrator() {
    for ( auto *item : m_warpper ) {
        if ( item != nullptr ) delete item;
    }
}


DirectIntegrator::DirectIntegrator(int bsdf_samples, int light_samples) : m_bsdf_samples(bsdf_samples), m_light_samples(light_samples) {
    PSDR_ASSERT((bsdf_samples >= 0) && (light_samples >= 0) && (bsdf_samples + light_samples > 0));
}


SpectrumC DirectIntegrator::Li(const Scene &scene, Sampler &sampler, const RayC &ray, MaskC active) const {
    return __Li<false>(scene, sampler, ray, active);
}


SpectrumD DirectIntegrator::Li(const Scene &scene, Sampler &sampler, const RayD &ray, MaskD active) const {
    return __Li<true>(scene, sampler, ray, active);
}


template <bool ad>
Spectrum<ad> DirectIntegrator::__Li(const Scene &scene, Sampler &sampler, const Ray<ad> &ray, Mask<ad> active) const {
    Intersection<ad> its = scene.ray_intersect<ad>(ray, active);
    active &= its.is_valid();

    Spectrum<ad> result = m_hide_emitters ? zero<Spectrum<ad>>() : its.Le(active);

    BSDFArray<ad> bsdf_array = its.shape->bsdf(active);
    if ( scene.m_emitter_env != nullptr ) {
        // Skip reflectance computations for intersections on the bounding mesh
        active &= neq(bsdf_array, nullptr);
    }

    const BSDF *bsdf = nullptr;
    if ( scene.m_bsdfs.size() == 1U || scene.m_meshes.size() == 1U ) {
        bsdf = scene.m_meshes[0]->m_bsdf;
    }

    if ( m_bsdf_samples > 0 ) {
        // BSDF sampling

        for ( int i = 0; i < m_bsdf_samples; ++i ) {
            BSDFSample<ad> bs;
            if ( bsdf != nullptr ) {
                bs = bsdf->sample(its, sampler.next_nd<3, ad>(), active);
            } else {
                bs = bsdf_array->sample(its, sampler.next_nd<3, ad>(), active);
            }
            Mask<ad> active1 = active && bs.is_valid;

            Ray<ad> ray1(its.p, its.sh_frame.to_world(bs.wo));
            Intersection<ad> its1 = scene.ray_intersect<ad, ad>(ray1, active1);
            active1 &= its1.is_valid();
            active1 &= neq(its1.shape->emitter(active1), nullptr);

            Spectrum<ad> bsdf_val;
            Float<ad> pdf0;
            if constexpr ( ad ) {
                Vector3fD wo = its1.p - its.p;
                wo /= its1.t;

                if ( bsdf != nullptr ) {
                    bsdf_val = bsdf->eval(its, its.sh_frame.to_local(wo), active1);
                } else {
                    bsdf_val = bsdf_array->eval(its, its.sh_frame.to_local(wo), active1);
                }
                FloatD cos_val = dot(its1.n, -wo);
                FloatD G_val = abs(cos_val)/sqr(its1.t);
                pdf0 = bs.pdf*detach(G_val);
                bsdf_val *= G_val*its1.J/pdf0;
            } else {
                if ( bsdf != nullptr ) {
                    bsdf_val = bsdf->eval(its, bs.wo, active1);
                } else {
                    bsdf_val = bsdf_array->eval(its, bs.wo, active1);
                }
                FloatC cos_val = dot(its1.n, -ray1.d);
                FloatC G_val = abs(cos_val)/sqr(its1.t);
                pdf0 = bs.pdf*G_val;
                bsdf_val /= bs.pdf;
            }

            Float<ad> weight = 1.f/static_cast<float>(m_bsdf_samples);
            if ( m_light_samples > 0 ) {
                weight *= mis_weight<ad>(pdf0, scene.emitter_position_pdf<ad>(its.p, its1, active1));
            }
            masked(result, active1) += its1.Le(active1)*bsdf_val*weight;
        }
    }

    if ( m_light_samples > 0 ) {
        // Light sampling

        for ( int i = 0; i < m_light_samples; ++i ) {
            PositionSample<ad> ps = scene.sample_emitter_position<ad>(its.p, sampler.next_2d<ad>(), active);
            Mask<ad> active1 = active && ps.is_valid;

            Vector3f<ad> wo = ps.p - its.p;
            Float<ad> dist_sqr = squared_norm(wo);
            Float<ad> dist = safe_sqrt(dist_sqr);
            wo /= dist;

            Ray<ad> ray1(its.p, wo);
            Intersection<ad> its1 = scene.ray_intersect<ad, ad>(ray1, active1);
            active1 &= its1.is_valid();
            active1 &= (its1.t > dist - ShadowEpsilon) && its1.is_emitter(active1);
            //ps.pdf = scene.emitter_position_pdf<ad>(its.p, its1, active1);

            Float<ad> cos_val = dot(its1.n, -wo);
            Float<ad> G_val = abs(cos_val) / dist_sqr;
            Spectrum<ad> bsdf_val;
            Float<ad> pdf1;
            Vector3f<ad> wo_local = its.sh_frame.to_local(wo);
            if ( bsdf != nullptr ) {
                bsdf_val = bsdf->eval(its, wo_local, active1);
                pdf1 = bsdf->pdf(its, wo_local, active1);
            } else {
                bsdf_val = bsdf_array->eval(its, wo_local, active1);
                pdf1 = bsdf_array->pdf(its, wo_local, active1);
            }
            bsdf_val *= G_val*ps.J/ps.pdf;

            if constexpr ( ad ) {
                pdf1 *= detach(G_val);
            } else {
                pdf1 *= G_val;
            }

            Float<ad> weight = 1.f/static_cast<float>(m_light_samples);
            if ( m_bsdf_samples > 0 ) {
                weight *= mis_weight<ad>(ps.pdf, pdf1);
            }
            masked(result, active1) += its1.Le(active1)*bsdf_val*weight;
        }
    }

    return result;
}


void DirectIntegrator::preprocess_secondary_edges(const Scene &scene, int sensor_id, const ScalarVector4i &reso, int nrounds) {
    PSDR_ASSERT(nrounds > 0);
    PSDR_ASSERT_MSG(scene.is_ready(), "Scene needs to be configured!");

    if ( static_cast<int>(m_warpper.size()) != scene.m_num_sensors )
        m_warpper.resize(scene.m_num_sensors, nullptr);

    if ( m_warpper[sensor_id] == nullptr )
        m_warpper[sensor_id] = new HyperCubeDistribution3f();
    auto warpper = m_warpper[sensor_id];

    warpper->set_resolution(head<3>(reso));
    int num_cells = warpper->m_num_cells;
    const int64_t num_samples = static_cast<int64_t>(num_cells)*reso[3];
    PSDR_ASSERT(num_samples <= std::numeric_limits<int>::max());

    IntC idx = divisor<int>(reso[3])(arange<IntC>(num_samples));
    Vector3iC sample_base = gather<Vector3iC>(warpper->m_cells, idx);

    Sampler sampler;
    sampler.seed(arange<UInt64C>(num_samples));

    FloatC result = zero<FloatC>(num_cells);
    for ( int j = 0; j < nrounds; ++j ) {
        SpectrumC value0;
        std::tie(std::ignore, value0) = eval_secondary_edge<false>(scene, *scene.m_sensors[sensor_id],
                                                                   (sample_base + sampler.next_nd<3, false>())*warpper->m_unit);
        masked(value0, ~enoki::isfinite<SpectrumC>(value0)) = 0.f;
        if ( likely(reso[3] > 1) ) {
            value0 /= static_cast<float>(reso[3]);
        }
        //PSDR_ASSERT(all(hmin(value0) > -Epsilon));
        scatter_add(result, hmax(value0), idx);
    }
    if ( nrounds > 1 ) result /= static_cast<float>(nrounds);
    warpper->set_mass(result);

    cuda_eval(); cuda_sync();
}


void DirectIntegrator::render_secondary_edges(const Scene &scene, int sensor_id, SpectrumD &result) const {
    const RenderOption &opts = scene.m_opts;

    Vector3fC sample3 = scene.m_samplers[2].next_nd<3, false>();
    FloatC pdf0 = (m_warpper.empty() || m_warpper[sensor_id] == nullptr) ?
                  1.f : m_warpper[sensor_id]->sample_reuse(sample3);

    auto [idx, value] = eval_secondary_edge<true>(scene, *scene.m_sensors[sensor_id], sample3);
    masked(value, ~enoki::isfinite<SpectrumD>(value)) = 0.f;
    masked(value, pdf0 > Epsilon) /= pdf0;
    if ( likely(opts.sppse > 1) ) {
        value /= static_cast<float>(opts.sppse);
    }
    scatter_add(result, value, IntD(idx), idx >= 0);
}


template <bool ad>
std::pair<IntC, Spectrum<ad>> DirectIntegrator::eval_secondary_edge(const Scene &scene, const Sensor &sensor, const Vector3fC &sample3) const {
    BoundarySegSampleDirect bss = scene.sample_boundary_segment_direct(sample3);
    MaskC valid = bss.is_valid;

    // _p0 on a face edge, _p2 on an emitter
    const Vector3fC &_p0    = detach(bss.p0);
    Vector3fC       &_p2    = bss.p2,
                    _dir    = normalize(_p2 - _p0);

    // check visibility between _p0 and _p2
    IntersectionC _its2;
    TriangleInfoD tri_info;
    if constexpr ( ad ) {
        _its2 = scene.ray_intersect<false>(RayC(_p0, _dir), valid, &tri_info);
    } else {
        _its2 = scene.ray_intersect<false>(RayC(_p0, _dir), valid);
    }
    valid &= _its2.is_valid() && norm(_its2.p - _p2) < ShadowEpsilon;

    // trace another ray in the opposite direction to complete the boundary segment (_p1, _p2)
    IntersectionC _its1 = scene.ray_intersect<false>(RayC(_p0, -_dir), valid);
    valid &= _its1.is_valid();
    Vector3fC &_p1 = _its1.p;

    // project _p1 onto the image plane and compute the corresponding pixel id
    SensorDirectSampleC sds = sensor.sample_direct(_p1);
    valid &= sds.is_valid;

    // trace a camera ray toward _p1 in a differentiable fashion
    Ray<ad> camera_ray;
    Intersection<ad> its1;
    if constexpr ( ad ) {
        camera_ray = sensor.sample_primary_ray(Vector2fD(sds.q));
        its1 = scene.ray_intersect<true, false>(camera_ray, valid);
        valid &= its1.is_valid() && norm(detach(its1.p) - _p1) < ShadowEpsilon;
    } else {
        camera_ray = sensor.sample_primary_ray(sds.q);
        its1 = scene.ray_intersect<false>(camera_ray, valid);
        valid &= its1.is_valid() && norm(its1.p - _p1) < ShadowEpsilon;
    }

    // calculate base_value
    FloatC      dist    = norm(_p2 - _p1),
                cos2    = abs(dot(bss.n, -_dir));
    Vector3fC   e       = cross(bss.edge, _dir);
    FloatC      sinphi  = norm(e);
    Vector3fC   proj    = normalize(cross(e, bss.n));
    FloatC      sinphi2 = norm(cross(_dir, proj));
    FloatC      base_v  = (_its1.t/dist)*(sinphi/sinphi2)*cos2;
    valid &= (sinphi > Epsilon) && (sinphi2 > Epsilon);

    // evaluate BSDF at _p1
    SpectrumC bsdf_val;
    Vector3fC d0;
    if constexpr ( ad ) {
        d0 = -detach(camera_ray.d);
    } else {
        d0 = -camera_ray.d;
    }
    Vector3fC d0_local = _its1.sh_frame.to_local(d0);
    if ( scene.m_bsdfs.size() == 1U || scene.m_meshes.size() == 1U ) {
        const BSDF *bsdf = scene.m_meshes[0]->m_bsdf;
        bsdf_val = bsdf->eval(_its1, d0_local, valid);
    } else {
        BSDFArrayC bsdf_array = _its1.shape->bsdf(valid);
        bsdf_val = bsdf_array->eval(_its1, d0_local, valid);
    }
    // accounting for BSDF's asymmetry caused by shading normals
    FloatC correction = abs((_its1.wi.z()*dot(d0, _its1.n))/(d0_local.z()*dot(_dir, _its1.n)));
    masked(bsdf_val, valid) *= correction;

    SpectrumC value0 = (bsdf_val*_its2.Le(valid)*(base_v*sds.sensor_val/bss.pdf)) & valid;
    if constexpr ( ad ) {
        Vector3fC n = normalize(cross(bss.n, proj));
        value0 *= sign(dot(e, bss.edge2))*sign(dot(e, n));

        const Vector3fD &v0 = tri_info.p0,
                        &e1 = tri_info.e1,
                        &e2 = tri_info.e2;

        RayD shadow_ray(its1.p, normalize(bss.p0 - its1.p));
        Vector2fD uv;
        std::tie(uv, std::ignore) = ray_intersect_triangle<true>(v0, e1, e2, shadow_ray);
        Vector3fD u2 = bilinear<true>(detach(v0), detach(e1), detach(e2), uv);

        SpectrumD result = (SpectrumD(value0)*dot(Vector3fD(n), u2)) & valid;
        return { select(valid, sds.pixel_idx, -1), result - detach(result) };
    } else {
        // returning the value without multiplying normal velocity for guiding
        return { -1, value0 };
    }
}

} // namespace psdr
