#pragma once

#include <psdr/psdr.h>
#include <psdr/core/records.h>

namespace psdr
{

PSDR_CLASS_DECL_BEGIN(Emitter,, Object)
public:
    virtual ~Emitter() override {}

    virtual void configure() = 0;

    // Returns the emitted radiance at its.p in direction its.wi

    virtual SpectrumC eval(const IntersectionC &its, MaskC active = true) const = 0;
    virtual SpectrumD eval(const IntersectionD &its, MaskD active = true) const = 0;

    virtual PositionSampleC sample_position(const Vector3fC &ref_p, const Vector2fC &sample2, MaskC active = true) const = 0;
    virtual PositionSampleD sample_position(const Vector3fD &ref_p, const Vector2fD &sample2, MaskD active = true) const = 0;

    virtual FloatC sample_position_pdf(const Vector3fC &ref_p, const IntersectionC &its, MaskC active = true) const = 0;
    virtual FloatD sample_position_pdf(const Vector3fD &ref_p, const IntersectionD &its, MaskD active = true) const = 0;

    bool m_ready = false;
    float m_sampling_weight = 1.f;

    ENOKI_PINNED_OPERATOR_NEW(FloatD)

PSDR_CLASS_DECL_END(Emitter)

} // namespace psdr

ENOKI_CALL_SUPPORT_BEGIN(psdr::Emitter)
    ENOKI_CALL_SUPPORT_METHOD(eval)
    ENOKI_CALL_SUPPORT_METHOD(sample_position)
    ENOKI_CALL_SUPPORT_METHOD(sample_position_pdf)
ENOKI_CALL_SUPPORT_END(psdr::Emitter)
