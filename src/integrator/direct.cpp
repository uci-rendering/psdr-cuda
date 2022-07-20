#include <misc/Exception.h>
#include <psdr/core/cube_distrb.h>
#include <psdr/core/ray.h>
#include <psdr/core/intersection.h>
#include <psdr/core/sampler.h>
#include <psdr/core/transform.h>
#include <psdr/bsdf/bsdf.h>
#include <psdr/emitter/emitter.h>
#include <psdr/emitter/envmap.h>

#include <psdr/shape/mesh.h>
#include <psdr/scene/scene.h>
#include <psdr/sensor/perspective.h>
#include <psdr/integrator/direct.h>

#include <psdr/core/AQ_distrb.h>

#include <fstream>

namespace psdr
{

DirectIntegrator::~DirectIntegrator() {
    for ( auto *item : m_warpper ) {
        if ( item != nullptr ) delete item;
    }

    for ( auto *item : m_aq ) {
        if ( item != nullptr ) delete item;
    }
}

FloatC DirectIntegrator::get_grid() {
#if 1 // For debug
    // m_aq m_aq_emitter
    if (!(m_aq.empty() || m_aq[0] == nullptr)) {
        Sampler sampler;
        sampler.seed(arange<UInt64C>(100*100*100));
        Histogram3D<100, 100, 100> histo;
        for (int i=0; i<64; ++i) {
            Vector3fC sample3 = sampler.next_3d<false>();
            FloatC pdf0 = 0.0f;
            Vector3fC aq_sample = m_aq[0]->sample(sample3, pdf0);;
            histo.update(aq_sample);
        }
        histo.normalize();
        return histo.get_data();

    } else if (m_warpper[0] != nullptr) {
        return m_warpper[0]->m_distrb.m_pmf;
    }
#endif
    return FloatC(0.0f);
}

DirectIntegrator::DirectIntegrator(int bsdf_samples, int light_samples, int edge_direct) : m_bsdf_samples(bsdf_samples), m_light_samples(light_samples), m_edge_direct(edge_direct) {
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

void DirectIntegrator::preprocess_secondary_edges(const Scene &scene, const std::vector<int> &sensor_id, const std::vector<float> &config, int option) {
    PSDR_ASSERT_MSG(scene.is_ready(), "Scene needs to be configured!");
    preprocess_option = option;
    if (option == 0 || option == 2) {

        std::cout << "Guiding using AQ Distribution with option: " << option << std::endl;

        if (option == 2 || !m_edge_direct) {
            if ( static_cast<int>(m_aq.size()) != scene.m_num_sensors )
                m_aq.resize(scene.m_num_sensors, nullptr);
            if ( m_aq[sensor_id[0]] == nullptr ) {
                m_aq[sensor_id[0]] = new AdaptiveQuadratureDistribution3f();
                m_aq[sensor_id[0]]->aq_edge_direct         = false;
            }
        }

        if (option == 2 || m_edge_direct) {
            if ( static_cast<int>(m_aq_emitter.size()) != scene.m_num_sensors )
                m_aq_emitter.resize(scene.m_num_sensors, nullptr);
            if ( m_aq_emitter[sensor_id[0]] == nullptr ) {
                m_aq_emitter[sensor_id[0]] = new AdaptiveQuadratureDistribution3f();
                m_aq_emitter[sensor_id[0]]->aq_edge_direct = true;
            }
        }
        std::cout << "AQ Distribution init" << std::endl;

        if (option == 0) {
            PSDR_ASSERT_MSG(config.size()==10, "Incorrect size of config");
            FloatC cdfx;
            bool enable_edge_sort = scene.m_edge_sort.enable_sort;
            if (enable_edge_sort == false) {
                cdfx = scene.m_sec_edge_distrb->cmf();
            } else {
                FloatC raw_cdf = scene.m_sec_edge_distrb->cmf();
                cdfx = gather<FloatC>(raw_cdf, scene.m_edge_cut);
            }
            // FloatC cdfx = (arange<FloatC>(int(config[0]))+FloatC(1.0)) / FloatC(config[0]);
            FloatC cdfy = (arange<FloatC>(int(config[1]))+FloatC(1.0)) / FloatC(config[1]);
            FloatC cdfz = (arange<FloatC>(int(config[2]))+FloatC(1.0)) / FloatC(config[2]);
            cuda_eval(); cuda_sync();

            AQ_Option aqconfig = AQ_Option(config, option);

            if (m_edge_direct) {
                m_aq_emitter[sensor_id[0]]->setup(scene, sensor_id, cdfx, cdfy, cdfz, aqconfig);
            } else {
                m_aq[sensor_id[0]]->setup(scene, sensor_id, cdfx, cdfy, cdfz, aqconfig);
            }
        } else if (option == 2) {
            std::cout << "option 2 is for MIS" << std::endl;
            PSDR_ASSERT_MSG(config.size()==20, "Incorrect size of config");
            std::vector<float> config1(&config[0],&config[10]); 
            std::vector<float> config2(&config[10],&config[20]); 


            FloatC cdfx1;
            bool enable_edge_sort = scene.m_edge_sort.enable_sort;
            if (enable_edge_sort == false) {
                cdfx1 = scene.m_sec_edge_distrb->cmf();
            } else {
                FloatC raw_cdf = scene.m_sec_edge_distrb->cmf();
                cdfx1 = gather<FloatC>(raw_cdf, scene.m_edge_cut);
            }
            // FloatC cdfx1 = (arange<FloatC>(int(config1[0]))+FloatC(1.0)) / FloatC(config1[0]);
            FloatC cdfy1 = (arange<FloatC>(int(config1[1]))+FloatC(1.0)) / FloatC(config1[1]);
            FloatC cdfz1 = (arange<FloatC>(int(config1[2]))+FloatC(1.0)) / FloatC(config1[2]);
            cuda_eval(); cuda_sync();
            AQ_Option aqconfig1 = AQ_Option(config1, 0);

            FloatC cdfx2;
            if (enable_edge_sort == false) {
                cdfx2 = scene.m_sec_edge_distrb->cmf();
            } else {
                FloatC raw_cdf = scene.m_sec_edge_distrb->cmf();
                cdfx2 = gather<FloatC>(raw_cdf, scene.m_edge_cut);
            }
            // FloatC cdfx2 = (arange<FloatC>(int(config2[0]))+FloatC(1.0)) / FloatC(config2[0]);
            FloatC cdfy2 = (arange<FloatC>(int(config2[1]))+FloatC(1.0)) / FloatC(config2[1]);
            FloatC cdfz2 = (arange<FloatC>(int(config2[2]))+FloatC(1.0)) / FloatC(config2[2]);
            cuda_eval(); cuda_sync();
            AQ_Option aqconfig2 = AQ_Option(config2, 0);
            m_aq_emitter[sensor_id[0]]->setup(scene, sensor_id, cdfx1, cdfy1, cdfz1, aqconfig1);
            m_aq[sensor_id[0]]->setup(scene, sensor_id, cdfx2, cdfy2, cdfz2, aqconfig2);

        } else {
            PSDR_ASSERT(0);
        }
        // PSDR_ASSERT(0);
        cuda_eval(); cuda_sync();
    } else if (option == 1) {
        std::cout << "Guiding using MC regular Cube Distribution: " << config[0] << " " << config[1] << " " << config[2] << std::endl;
        int nrounds = int(config[4]);
        if ( static_cast<int>(m_warpper.size()) != scene.m_num_sensors )
            m_warpper.resize(scene.m_num_sensors, nullptr);

        if ( m_warpper[sensor_id[0]] == nullptr )
            m_warpper[sensor_id[0]] = new HyperCubeDistribution3f();
        auto warpper = m_warpper[sensor_id[0]];

        Array<int, 3> wconfig{int(config[0]), int(config[1]), int(config[2])};

        warpper->set_resolution(wconfig);
        int num_cells = warpper->m_num_cells;
        const int64_t num_samples = static_cast<int64_t>(num_cells)*int(config[3]);
        PSDR_ASSERT(num_samples <= std::numeric_limits<int>::max());

        IntC idx = divisor<int>(int(config[3]))(arange<IntC>(num_samples));
        Vector3iC sample_base = gather<Vector3iC>(warpper->m_cells, idx);

        Sampler sampler;
        sampler.seed(arange<UInt64C>(num_samples));

        FloatC result = zero<FloatC>(num_cells);
        for ( int j = 0; j < nrounds; ++j ) {
            // SpectrumC value0;

            BoundaryMISRecordD BMR1 = __eval_secondary_edgeMIS<false>(scene, *scene.m_sensors[sensor_id[0]],
                                                                        (sample_base + sampler.next_nd<3, false>())*warpper->m_unit, m_edge_direct == 1);
            SpectrumC value0 = detach(BMR1.value);
            // std::tie(std::ignore, value0) = eval_secondary_edge<false>(scene, *scene.m_sensors[sensor_id[0]],
            //                                                            (sample_base + sampler.next_nd<3, false>())*warpper->m_unit);
            // std::cout << value0 << std::endl;
            cuda_eval(); cuda_sync();

            masked(value0, ~enoki::isfinite<SpectrumC>(value0)) = 0.000f;
            if ( likely(config[3] > 1) ) {
                value0 /= static_cast<float>(config[3]);
            }
            PSDR_ASSERT(all(hmin(value0) > -Epsilon));
            scatter_add(result, hmax(value0), idx);
        }
        if ( nrounds > 1 ) result /= static_cast<float>(nrounds);
        result /= hsum(result); // new code
        float ueps = 0.0f;
        result = ueps + (1.0f - ueps)*result;
        warpper->set_mass(result);
        cuda_eval(); cuda_sync();
    } else {
        PSDR_ASSERT_MSG(0, "ERROR no such config for guiding SE");
    }
}


void DirectIntegrator::render_secondary_edges(const Scene &scene, int sensor_id, SpectrumD &result) const {
    const RenderOption &opts = scene.m_opts;
    Vector3fC sample3 = scene.m_samplers[2].next_nd<3, false>();

    if (preprocess_option == 0) {
        if ((preprocess_option == 0 && m_edge_direct == 0) || preprocess_option == 2) {
            FloatC pdf1 = 1.f;
            Vector3fC aq_sample1;

            aq_sample1 = m_aq[sensor_id]->sample(sample3, pdf1);

            auto [idx1, value1] = __eval_secondary_edge<true>(scene, *scene.m_sensors[sensor_id], aq_sample1, false);
            masked(value1, ~enoki::isfinite<SpectrumD>(value1)) = 0.f;
            if ( likely(opts.sppse > 1) ) {
                value1 /= static_cast<float>(opts.sppse);
            }
            masked(value1, pdf1 > Epsilon) /= pdf1;
            cuda_eval(); cuda_sync();
            if (preprocess_option == 2) {
                FloatC pdf1_1 = m_aq[sensor_id]->pdf_mis(scene, sensor_id, aq_sample1);
                FloatC pdf1_2 = m_aq_emitter[sensor_id]->pdf_mis(scene, sensor_id, aq_sample1);
                masked(value1, pdf1 > Epsilon) *= mis_weight<false>(pdf1_1, pdf1_2);
            }
            scatter_add(result, value1, IntD(idx1), idx1 >= 0);

        } else {
            if (preprocess_option == 2 || !m_edge_direct) {
                PSDR_ASSERT_MSG(0, "m_aq is empty!");
            }
        }

        sample3 = scene.m_samplers[2].next_nd<3, false>();
        if ((preprocess_option == 0 && m_edge_direct == 1) || preprocess_option == 2) {
            FloatC pdf2 = 1.f;
            Vector3fC aq_sample2;

            aq_sample2 = m_aq_emitter[sensor_id]->sample(sample3, pdf2);
            auto [idx2, value2] = __eval_secondary_edge<true>(scene, *scene.m_sensors[sensor_id], aq_sample2, true);
            masked(value2, ~enoki::isfinite<SpectrumD>(value2)) = 0.f;
            if ( likely(opts.sppse > 1) ) {
                value2 /= static_cast<float>(opts.sppse);
            }
            masked(value2, pdf2 > Epsilon) /= pdf2;
            cuda_eval(); cuda_sync();
            if (preprocess_option == 2) {
                FloatC pdf2_1 = m_aq_emitter[sensor_id]->pdf_mis(scene, sensor_id, aq_sample2);
                FloatC pdf2_2 = m_aq[sensor_id]->pdf_mis(scene, sensor_id, aq_sample2);
                masked(value2, pdf2 > Epsilon) *= mis_weight<false>(pdf2_1, pdf2_2);
            }
            scatter_add(result, value2, IntD(idx2), idx2 >= 0);
        } else {
            if (preprocess_option == 2 || m_edge_direct) {
                PSDR_ASSERT_MSG(0, "m_aq_emitter is empty!");
            }
        }

    } else {
        if (preprocess_option == 2) {
            PSDR_ASSERT_MSG(m_edge_direct == 2, "Not using MIS but guide MIS!");

            FloatC pdf0 = 1.f;
            Vector3fC aq_sample0;
            aq_sample0 = m_aq[sensor_id]->sample(sample3, pdf0);

            BoundaryMISRecordD BMR0 = __eval_secondary_edgeMIS<true>(scene, *scene.m_sensors[sensor_id], aq_sample0, false);
            masked(BMR0.value, ~enoki::isfinite<SpectrumD>(BMR0.value)) = 0.f;
            if ( likely(opts.sppse > 1) ) {
                BMR0.value /= static_cast<float>(opts.sppse);
            }
            masked(BMR0.value, pdf0 > Epsilon) /= pdf0;
            Vector3fC radiance1 = scene.m_emitter_env->eval_direction<false>(BMR0.dir);
            FloatC pdf0_1 = hmax(radiance1) / scene.m_emitter_env->m_cube_distrb.m_distrb.m_sum;
            pdf0_1 /= 3.0f;
            masked(BMR0.value, BMR0.pdf > Epsilon) *= mis_weight<false>(BMR0.pdf, pdf0_1);
            scatter_add(result, BMR0.value, IntD(BMR0.idx), BMR0.idx >= 0);


            // emitter sampling
            sample3 = scene.m_samplers[2].next_nd<3, false>();
            FloatC pdf1 = 1.f;
            Vector3fC aq_sample1;
            aq_sample1 = m_aq_emitter[sensor_id]->sample(sample3, pdf1);
            BoundaryMISRecordD BMR1 = __eval_secondary_edgeMIS<true>(scene, *scene.m_sensors[sensor_id], aq_sample1, true);

            masked(BMR1.value, ~enoki::isfinite<SpectrumD>(BMR1.value)) = 0.f;
            if ( likely(opts.sppse > 1) ) {
                BMR1.value /= static_cast<float>(opts.sppse);
            }
            masked(BMR1.value, pdf1 > Epsilon) /= pdf1;
            Vector3fC radiance2 = scene.m_emitter_env->eval_direction<false>(BMR1.dir);
            FloatC pdf1_0 = hmax(radiance2) / scene.m_emitter_env->m_cube_distrb.m_distrb.m_sum;
            pdf1_0 /= 3.0f;
            BoundarySegSampleDirect BMR1_1 = scene.sample_edge_ray(aq_sample1);
            BMR1.value *= mis_weight<false>(pdf1_0, BMR1_1.pdf);
            scatter_add(result, BMR1.value, IntD(BMR1.idx), BMR1.idx >= 0);

            // PSDR_ASSERT_MSG(0, "Guiding MIS not ready");

        } else if (m_edge_direct == 2 && preprocess_option == -1) {
            // direction sampling
            FloatC pdf0 = 1.f;
            BoundaryMISRecordD BMR0 = __eval_secondary_edgeMIS<true>(scene, *scene.m_sensors[sensor_id], sample3, false);
            masked(BMR0.value, ~enoki::isfinite<SpectrumD>(BMR0.value)) = 0.f;
            masked(BMR0.value, BMR0.pdf > Epsilon) /= pdf0;
            if ( likely(opts.sppse > 1) ) {
                BMR0.value /= static_cast<float>(opts.sppse);
            }
            Vector3fC radiance1 = scene.m_emitter_env->eval_direction<false>(BMR0.dir);
            FloatC pdf0_1 = hmax(radiance1) / scene.m_emitter_env->m_cube_distrb.m_distrb.m_sum;
            pdf0_1 /= 3.0f;
            masked(BMR0.value, BMR0.pdf > Epsilon) *= mis_weight<false>(BMR0.pdf, pdf0_1);
            scatter_add(result, BMR0.value, IntD(BMR0.idx), BMR0.idx >= 0);

            // emitter sampling
            sample3 = scene.m_samplers[2].next_nd<3, false>();
            FloatC pdf1 = 1.f;
            BoundaryMISRecordD BMR1 = __eval_secondary_edgeMIS<true>(scene, *scene.m_sensors[sensor_id], sample3, true);
            masked(BMR1.value, ~enoki::isfinite<SpectrumD>(BMR1.value)) = 0.f;
            masked(BMR1.value, pdf1 > Epsilon) /= pdf1;
            if ( likely(opts.sppse > 1) ) {
                BMR1.value /= static_cast<float>(opts.sppse);
            }
            Vector3fC radiance2 = scene.m_emitter_env->eval_direction<false>(BMR1.dir);
            FloatC pdf1_0 = hmax(radiance2) / scene.m_emitter_env->m_cube_distrb.m_distrb.m_sum;
            pdf1_0 /= 3.0f;

            BoundarySegSampleDirect BMR1_1 = scene.sample_edge_ray(sample3);

            BMR1.value *= mis_weight<false>(pdf1_0, BMR1_1.pdf);
            scatter_add(result, BMR1.value, IntD(BMR1.idx), BMR1.idx >= 0);

        } else {

            FloatC pdf0 = (m_warpper.empty() || m_warpper[sensor_id] == nullptr) ?
                           1.f : m_warpper[sensor_id]->sample_reuse(sample3);
            auto [idx, value] = eval_secondary_edge<true>(scene, *scene.m_sensors[sensor_id], sample3);
            masked(value, ~enoki::isfinite<SpectrumD>(value)) = 0.f;
            masked(value, pdf0 > Epsilon) /= pdf0;
            if ( likely(opts.sppse > 1) ) {
                value /= static_cast<float>(opts.sppse);
            }
            scatter_add(result, value, IntD(idx), idx >= 0);

            // FloatC pdf0 = (m_warpper.empty() || m_warpper[sensor_id] == nullptr) ?
            //                1.f : m_warpper[sensor_id]->sample_reuse(sample3);
            // // std::cout << "debug MIS on: " << m_edge_direct << std::endl;
            // BoundaryMISRecordD BMR1 = __eval_secondary_edgeMIS<true>(scene, *scene.m_sensors[sensor_id], sample3, m_edge_direct == 1);
            // masked(BMR1.value, ~enoki::isfinite<SpectrumD>(BMR1.value)) = 0.f;
            // masked(BMR1.value, pdf0 > Epsilon) /= pdf0;
            // if ( likely(opts.sppse > 1) ) {
            //     BMR1.value /= static_cast<float>(opts.sppse);
            // }
            // scatter_add(result, BMR1.value, IntD(BMR1.idx), BMR1.idx >= 0);
        }
    }
}

template <bool ad>
BoundaryMISRecord<ad> DirectIntegrator::__eval_secondary_edgeMIS(const Scene &scene, const Sensor &sensor, const Vector3fC &sample3, bool emitter_sampling) const {
    BoundarySegSampleDirect bss;
    if (emitter_sampling) {
        bss = scene.sample_emitter_ray(sample3);

    } else {
        bss = scene.sample_edge_ray(sample3);
    }

    MaskC valid = bss.is_valid;

    // _p0 on a face edge, _p2 on an emitter
    const Vector3fC &_p0    = detach(bss.p0);
    Vector3fC       _p2, _dir;

    _dir = bss.p2;

    // check visibility between _p0 and _p2
    IntersectionC _its2;
    TriangleInfoD tri_info;
    if constexpr ( ad ) {
        _its2 = scene.ray_intersect<false>(RayC(_p0, _dir), valid, &tri_info);
    } else {
        _its2 = scene.ray_intersect<false>(RayC(_p0, _dir), valid); 
    }

    _p2 = _its2.p;

    valid &= _its2.is_emitter(valid) && _its2.is_valid() && norm(_its2.p - _p2) < ShadowEpsilon;

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

    BSDFArray<ad> bsdf_array = its1.shape->bsdf(valid);
    if ( scene.m_emitter_env != nullptr ) {
        valid &= neq(bsdf_array, nullptr);
    }

    // calculate base_value
    FloatC      dist    = norm(_p2 - _p1),
                cos2    = abs(dot(_its2.n, -_dir));
    Vector3fC   e       = cross(bss.edge, _dir);
    FloatC      sinphi  = norm(e);
    Vector3fC   proj    = normalize(cross(e, _its2.n));
    FloatC      sinphi2 = norm(cross(_dir, proj));
    FloatC      base_v  = (_its1.t/dist)*(sinphi/sinphi2)*cos2;
    // valid &= (sinphi > Epsilon) && (sinphi2 > Epsilon);

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
        Vector3fC n = normalize(cross(_its2.n, proj));
        value0 *= sign(dot(e, bss.edge2))*sign(dot(e, n));

        const Vector3fD &v0 = tri_info.p0,
                        &e1 = tri_info.e1,
                        &e2 = tri_info.e2;

        RayD shadow_ray(its1.p, normalize(bss.p0 - its1.p));
        Vector2fD uv;
        std::tie(uv, std::ignore) = ray_intersect_triangle<true>(v0, e1, e2, shadow_ray);
        Vector3fD u2 = bilinear<true>(detach(v0), detach(e1), detach(e2), uv);

        SpectrumD result = (SpectrumD(value0)*dot(Vector3fD(n), u2)) & valid;

        BoundaryMISRecord<ad> resultMIS;
        resultMIS.p0 = _p0;
        resultMIS.dir = _dir;
        resultMIS.idx = select(valid, sds.pixel_idx, -1);
        resultMIS.value = result - detach(result);
        resultMIS.pdf = bss.pdf;
        resultMIS.is_valid = valid;
        resultMIS.bss = bss;
        return resultMIS;
    } else {
        BoundaryMISRecord<ad> resultMIS;
        resultMIS.p0 = _p0;
        resultMIS.dir = _dir;
        resultMIS.idx = -1;
        resultMIS.value = value0;
        resultMIS.pdf = bss.pdf;
        resultMIS.is_valid = valid;
        resultMIS.bss = bss;
        // returning the value without multiplying normal velocity for guiding
        return resultMIS;
    }
}


template <bool ad>
std::pair<IntC, Spectrum<ad>> DirectIntegrator::__eval_secondary_edge(const Scene &scene, const Sensor &sensor, const Vector3fC &sample3, bool emitter_sampling) const {
    BoundarySegSampleDirect bss;
    if (emitter_sampling) {
        bss = scene.sample_boundary_segment_direct(sample3);
        // std::cout << "emitter_sampling pdf: " << bss.pdf << std::endl;
    } else {
        bss = scene.sample_edge_ray(sample3);
        // std::cout << "sample_edge_ray pdf: " << bss.pdf << std::endl;
    }

    MaskC valid = bss.is_valid;

    // _p0 on a face edge, _p2 on an emitter
    const Vector3fC &_p0    = detach(bss.p0);
    Vector3fC       _p2, _dir;

    if (emitter_sampling) {
        _p2  = bss.p2;
        _dir = normalize(_p2 - _p0);
    } else {
        _dir = bss.p2;
    }

    // check visibility between _p0 and _p2
    IntersectionC _its2;
    TriangleInfoD tri_info;
    if constexpr ( ad ) {
        _its2 = scene.ray_intersect<false>(RayC(_p0, _dir), valid, &tri_info);
    } else {
        _its2 = scene.ray_intersect<false>(RayC(_p0, _dir), valid); 
    }

    if (!emitter_sampling) {
        _p2 = _its2.p;
    }

    valid &= _its2.is_emitter(valid) && _its2.is_valid() && norm(_its2.p - _p2) < ShadowEpsilon;

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

    BSDFArray<ad> bsdf_array = its1.shape->bsdf(valid);
    if ( scene.m_emitter_env != nullptr ) {
        valid &= neq(bsdf_array, nullptr);
    }

    // calculate base_value
    FloatC      dist    = norm(_p2 - _p1),
                cos2    = abs(dot(_its2.n, -_dir));
    Vector3fC   e       = cross(bss.edge, _dir);
    FloatC      sinphi  = norm(e);
    Vector3fC   proj    = normalize(cross(e, _its2.n));
    FloatC      sinphi2 = norm(cross(_dir, proj));
    FloatC      base_v  = (_its1.t/dist)*(sinphi/sinphi2)*cos2;
    // valid &= (sinphi > Epsilon) && (sinphi2 > Epsilon);

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
        Vector3fC n = normalize(cross(_its2.n, proj));
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


template <bool ad>
std::pair<IntC, Spectrum<ad>> DirectIntegrator::eval_secondary_edge(const Scene &scene, const Sensor &sensor, const Vector3fC &sample3) const {
    BoundarySegSampleDirect bss;
    if (m_edge_direct) {
        bss = scene.sample_boundary_segment_direct(sample3);
    } else {
        bss = scene.sample_edge_ray(sample3);
    }

    MaskC valid = bss.is_valid;

    // _p0 on a face edge, _p2 on an emitter
    const Vector3fC &_p0    = detach(bss.p0);
    Vector3fC       _p2, _dir;

    if (m_edge_direct) {
        _p2  = bss.p2;
        _dir = normalize(_p2 - _p0);
    } else {
        _dir = bss.p2;
    }

    // check visibility between _p0 and _p2
    IntersectionC _its2;
    TriangleInfoD tri_info;
    if constexpr ( ad ) {
        _its2 = scene.ray_intersect<false>(RayC(_p0, _dir), valid, &tri_info);
    } else {
        _its2 = scene.ray_intersect<false>(RayC(_p0, _dir), valid); 
    }

    if (!m_edge_direct) {
        _p2 = _its2.p;
    }

    valid &= _its2.is_emitter(valid) && _its2.is_valid() && norm(_its2.p - _p2) < ShadowEpsilon;

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

    BSDFArray<ad> bsdf_array = its1.shape->bsdf(valid);
    if ( scene.m_emitter_env != nullptr ) {
        valid &= neq(bsdf_array, nullptr);
    }

    // calculate base_value
    FloatC      dist    = norm(_p2 - _p1),
                cos2    = abs(dot(_its2.n, -_dir));
    Vector3fC   e       = cross(bss.edge, _dir);
    FloatC      sinphi  = norm(e);
    Vector3fC   proj    = normalize(cross(e, _its2.n));
    FloatC      sinphi2 = norm(cross(_dir, proj));
    FloatC      base_v  = (_its1.t/dist)*(sinphi/sinphi2)*cos2;
    // valid &= (sinphi > Epsilon) && (sinphi2 > Epsilon);

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

    // masked(value0, !valid) = 1.f;

    if constexpr ( ad ) {
        Vector3fC n = normalize(cross(_its2.n, proj));
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
