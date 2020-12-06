#pragma once

#include <string>
#include <psdr/psdr.h>

namespace psdr
{

PSDR_CLASS_DECL_BEGIN(Integrator,, Object)
public:
    virtual ~Integrator() {}

    SpectrumC renderC(const Scene &scene, int sensor_id = 0) const;
    SpectrumD renderD(const Scene &scene, int sensor_id = 0) const;

    virtual void preprocess_secondary_edges(const Scene &scene, int sensor_id, const ScalarVector4i &reso, int nrounds = 1) {}

protected:
    virtual SpectrumC Li(const Scene &scene, Sampler &sampler, const RayC &ray, MaskC active = true) const = 0;
    virtual SpectrumD Li(const Scene &scene, Sampler &sampler, const RayD &ray, MaskD active = true) const = 0;

    virtual void render_primary_edges(const Scene &scene, int sensor_id, SpectrumD &result) const;

    virtual void render_secondary_edges(const Scene &scene, int sensor_id, SpectrumD &result) const {}

    template <bool ad>
    Spectrum<ad> __render(const Scene &scene, int sensor_id) const;
PSDR_CLASS_DECL_END(SamplingIntegrator)

} // namespace psdr
