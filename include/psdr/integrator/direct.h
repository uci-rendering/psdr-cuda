#pragma once

#include "integrator.h"

namespace psdr
{

PSDR_CLASS_DECL_BEGIN(DirectIntegrator, final, Integrator)
public:
    DirectIntegrator(int bsdf_samples = 1, int light_samples = 1, int edge_direct = 1);
    virtual ~DirectIntegrator();

    FloatC get_grid();
    
    void preprocess_secondary_edges(const Scene &scene, const std::vector<int> &sensor_id, const std::vector<float> &config, int option = 0) override;
    
    bool m_hide_emitters = false;

protected:
    SpectrumC Li(const Scene &scene, Sampler &sampler, const RayC &ray, MaskC active = true) const override;
    SpectrumD Li(const Scene &scene, Sampler &sampler, const RayD &ray, MaskD active = true) const override;

    template <bool ad>
    Spectrum<ad> __Li(const Scene &scene, Sampler &sampler, const Ray<ad> &ray, Mask<ad> active) const;

    void render_secondary_edges(const Scene &scene, int sensor_id, SpectrumD &result) const override;

    template <bool ad>
    std::pair<IntC, Spectrum<ad>> eval_secondary_edge(const Scene &scene, const Sensor &sensor, const Vector3fC &sample3) const;
    
    template <bool ad>
    std::pair<IntC, Spectrum<ad>> __eval_secondary_edge(const Scene &scene, const Sensor &sensor, const Vector3fC &sample3, bool emitter_sampling) const;

    template <bool ad>
    BoundaryMISRecord<ad> __eval_secondary_edgeMIS(const Scene &scene, const Sensor &sensor, const Vector3fC &sample3, bool emitter_sampling) const;

    int m_bsdf_samples, m_light_samples, m_edge_direct;
    std::vector<HyperCubeDistribution3f*> m_warpper;
    std::vector<AdaptiveQuadratureDistribution3f*>      m_aq;
    std::vector<AdaptiveQuadratureDistribution3f*>      m_aq_emitter;

    AdaptiveQuadratureDistribution3f*                   m_global_aq = nullptr;
    int preprocess_option = -1;

PSDR_CLASS_DECL_END(DirectIntegrator)

} // namespace psdr
