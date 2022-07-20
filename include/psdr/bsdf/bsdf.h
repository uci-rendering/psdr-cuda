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
    Float<ad> eta;

    ENOKI_DERIVED_STRUCT(BSDFSample_, Base,
        ENOKI_BASE_FIELDS(pdf, is_valid),
        ENOKI_DERIVED_FIELDS(wo, eta)
    )
};

template <typename Float_>
struct BSDFSampleDual_ : public SampleRecordDual_<Float_> {
    PSDR_IMPORT_BASE(SampleRecordDual_<Float_>, ad, pdf1, pdf2, is_valid1, is_valid2)
    Vector3f<ad> wo1, wo2;
    Mask<ad>     isDual;
    ENOKI_DERIVED_STRUCT(BSDFSampleDual_, Base,
        ENOKI_BASE_FIELDS(pdf1, pdf2, is_valid1, is_valid2),
        ENOKI_DERIVED_FIELDS(wo1, wo2, isDual)
    )
};

PSDR_CLASS_DECL_BEGIN(BSDF,, Object)
public:
    virtual ~BSDF() override {}

    virtual SpectrumC eval(const IntersectionC &its, const Vector3fC &wo, MaskC active = true) const = 0;
    virtual SpectrumD eval(const IntersectionD &its, const Vector3fD &wo, MaskD active = true) const = 0;

    virtual BSDFSampleC sample(const IntersectionC &its, const Vector3fC &sample, MaskC active = true) const = 0;
    virtual BSDFSampleD sample(const IntersectionD &its, const Vector3fD &sample, MaskD active = true) const = 0;

    virtual BSDFSampleDualC sampleDual(const IntersectionC &its, const Vector3fC &sample, MaskC active = true) const = 0;
    virtual BSDFSampleDualD sampleDual(const IntersectionD &its, const Vector3fD &sample, MaskD active = true) const = 0;

    virtual FloatC pdf(const IntersectionC &its, const Vector3fC &wo, MaskC active) const = 0;
    virtual FloatD pdf(const IntersectionD &its, const Vector3fD &wo, MaskD active) const = 0;

    virtual bool anisotropic() const = 0;
PSDR_CLASS_DECL_END(BSDF)

} // namespace psdr

ENOKI_STRUCT_SUPPORT(psdr::BSDFSample_, pdf, is_valid, wo, eta)
ENOKI_STRUCT_SUPPORT(psdr::BSDFSampleDual_, pdf1, pdf2, is_valid1, is_valid2, wo1, wo2, isDual)

ENOKI_CALL_SUPPORT_BEGIN(psdr::BSDF)
    ENOKI_CALL_SUPPORT_METHOD(eval)
    ENOKI_CALL_SUPPORT_METHOD(sample)
    ENOKI_CALL_SUPPORT_METHOD(sampleDual)
    ENOKI_CALL_SUPPORT_METHOD(pdf)
    ENOKI_CALL_SUPPORT_METHOD(anisotropic)
ENOKI_CALL_SUPPORT_END(psdr::BSDF)
