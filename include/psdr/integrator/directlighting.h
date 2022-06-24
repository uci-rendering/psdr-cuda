#pragma once

#include "direct.h"

namespace psdr
{

PSDR_CLASS_DECL_BEGIN(DirectLightingIntegrator, final, DirectIntegrator)
public:
    DirectLightingIntegrator(int bsdf_samples = 1, int light_samples = 1);

protected:
    SpectrumC Li(const Scene &scene, Sampler &sampler, const RayC &ray, MaskC active = true) const override;
    SpectrumD Li(const Scene &scene, Sampler &sampler, const RayD &ray, MaskD active = true) const override;

    template <bool ad>
    Spectrum<ad> __Li(const Scene &scene, Sampler &sampler, const Ray<ad> &ray, Mask<ad> active) const;

    std::pair<IntC, SpectrumC> eval_secondary_edgeC(const Scene &scene, const Sensor &sensor, const Vector3fC &sample3) const override;
    std::pair<IntC, SpectrumD> eval_secondary_edgeD(const Scene &scene, const Sensor &sensor, const Vector3fC &sample3) const override;

    template <bool ad>
    std::pair<IntC, Spectrum<ad>> __eval_secondary_edge(const Scene &scene, const Sensor &sensor, const Vector3fC &sample3) const;
PSDR_CLASS_DECL_END(DirectLightingIntegrator)

} // namespace psdr