#pragma once

#include <string>
#include <psdr/psdr.h>

namespace psdr
{

PSDR_CLASS_DECL_BEGIN(Integrator,, Object)
public:
    virtual ~Integrator() {}

    SpectrumC renderC(const Scene &scene, int sensor_id = 0, int npass=1) const;
    SpectrumD renderD(const Scene &scene, int sensor_id = 0) const;

    virtual void preprocess_secondary_edges(const Scene &scene, const std::vector<int> &sensor_id, const std::vector<float> &config, int option = 0) {}

protected:
    virtual SpectrumC Li(const Scene &scene, Sampler &sampler, const RayC &ray, MaskC active = true) const = 0;
    virtual SpectrumD Li(const Scene &scene, Sampler &sampler, const RayD &ray, MaskD active = true) const = 0;

    virtual void render_primary_edges(const Scene &scene, int sensor_id, SpectrumD &result) const;

    virtual void render_secondary_edges(const Scene &scene, int sensor_id, SpectrumD &result) const {}

    template <bool ad>
    Spectrum<ad> __render(const Scene &scene, int sensor_id) const;
PSDR_CLASS_DECL_END(SamplingIntegrator)

} // namespace psdr
