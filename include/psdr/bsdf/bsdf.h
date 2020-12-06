#pragma once

#include <psdr/psdr.h>
#include <psdr/core/intersection.h>
#include <psdr/core/records.h>

namespace psdr
{

template <typename Float_>
struct BSDFSample_ : public SampleRecord_<Float_> {
    PSDR_IMPORT_BASE(SampleRecord_<Float_>, ad, pdf, is_valid)

    Vector3f<ad> wo;

    ENOKI_DERIVED_STRUCT(BSDFSample_, Base,
        ENOKI_BASE_FIELDS(pdf, is_valid),
        ENOKI_DERIVED_FIELDS(wo)
    )
};


PSDR_CLASS_DECL_BEGIN(BSDF,, Object)
public:
    virtual ~BSDF() override {}

    virtual SpectrumC eval(const IntersectionC &its, const Vector3fC &wo, MaskC active = true) const = 0;
    virtual SpectrumD eval(const IntersectionD &its, const Vector3fD &wo, MaskD active = true) const = 0;

    virtual BSDFSampleC sample(const IntersectionC &its, const Vector3fC &sample, MaskC active = true) const = 0;
    virtual BSDFSampleD sample(const IntersectionD &its, const Vector3fD &sample, MaskD active = true) const = 0;

    virtual FloatC pdf(const IntersectionC &its, const Vector3fC &wo, MaskC active) const = 0;
    virtual FloatD pdf(const IntersectionD &its, const Vector3fD &wo, MaskD active) const = 0;

    virtual bool anisotropic() const = 0;
PSDR_CLASS_DECL_END(BSDF)

} // namespace psdr

ENOKI_STRUCT_SUPPORT(psdr::BSDFSample_, pdf, is_valid, wo)

ENOKI_CALL_SUPPORT_BEGIN(psdr::BSDF)
    ENOKI_CALL_SUPPORT_METHOD(eval)
    ENOKI_CALL_SUPPORT_METHOD(sample)
    ENOKI_CALL_SUPPORT_METHOD(pdf)
    ENOKI_CALL_SUPPORT_METHOD(anisotropic)
ENOKI_CALL_SUPPORT_END(psdr::BSDF)
