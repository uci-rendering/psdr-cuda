#pragma once

#include <psdr/core/bitmap.h>
#include "bsdf.h"

namespace psdr
{

PSDR_CLASS_DECL_BEGIN(Diffuse, final, BSDF)
public:
    Diffuse() : m_reflectance(0.5f) {}
    Diffuse(const ScalarVector3f &ref) : m_reflectance(ref) {}
    Diffuse(const char *refl_file);
    Diffuse(const Bitmap3fD &reflectance);

    SpectrumC eval(const IntersectionC &its, const Vector3fC &wo, MaskC active = true) const override;
    SpectrumD eval(const IntersectionD &its, const Vector3fD &wo, MaskD active = true) const override;

    BSDFSampleC sample(const IntersectionC &its, const Vector3fC &sample, MaskC active = true) const override;
    BSDFSampleD sample(const IntersectionD &its, const Vector3fD &sample, MaskD active = true) const override;

    FloatC pdf(const IntersectionC &its, const Vector3fC &wo, MaskC active) const override;
    FloatD pdf(const IntersectionD &its, const Vector3fD &wo, MaskD active) const override;

    bool anisotropic() const override { return false; }

    SpectrumC albedo(const IntersectionC &its, MaskC active) const override;
    SpectrumD albedo(const IntersectionD &its, MaskD active) const override;

    std::string to_string() const override { return std::string("Diffuse[id=") + m_id + "]"; }

    Bitmap3fD m_reflectance;

protected:
    template <bool ad>
    Spectrum<ad> __eval(const Intersection<ad>&, const Vector3f<ad>&, Mask<ad>) const;

    template <bool ad>
    BSDFSample<ad> __sample(const Intersection<ad>&, const Vector3f<ad>&, Mask<ad>) const;

    template <bool ad>
    Float<ad> __pdf(const Intersection<ad> &, const Vector3f<ad> &, Mask<ad>) const;

    template <bool ad>
    Spectrum<ad> __albedo(const Intersection<ad>&, Mask<ad>) const;
PSDR_CLASS_DECL_END(Diffuse)

} // namespace psdr
