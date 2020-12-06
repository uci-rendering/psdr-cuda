#pragma once

#include <psdr/psdr.h>
#include "emitter.h"

namespace psdr
{

PSDR_CLASS_DECL_BEGIN(AreaLight, final, Emitter)
public:
    AreaLight(const ScalarVector3f &radiance, const Mesh *mesh) : m_radiance(radiance), m_mesh(mesh) {}

    void configure() override;

    SpectrumC eval(const IntersectionC &its, MaskC active = true) const override;
    SpectrumD eval(const IntersectionD &its, MaskD active = true) const override;

    PositionSampleC sample_position(const Vector3fC &ref_p, const Vector2fC &sample2, MaskC active = true) const override;
    PositionSampleD sample_position(const Vector3fD &ref_p, const Vector2fD &sample2, MaskD active = true) const override;

    FloatC sample_position_pdf(const Vector3fC &ref_p, const IntersectionC &its, MaskC active = true) const override;
    FloatD sample_position_pdf(const Vector3fD &ref_p, const IntersectionD &its, MaskD active = true) const override;

    std::string to_string() const override;

    SpectrumD m_radiance;
    const Mesh *m_mesh;

    ENOKI_PINNED_OPERATOR_NEW(FloatD)

protected:
    template <bool ad>
    PositionSample<ad> __sample_position(const Vector2f<ad>&, Mask<ad>) const;

    template <bool ad>
    Float<ad> __sample_position_pdf(const Vector3f<ad>&, const Intersection<ad>&, Mask<ad>) const;
PSDR_CLASS_DECL_END(AreaLight)

} // namespace psdr
