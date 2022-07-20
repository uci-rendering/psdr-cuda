#pragma once

#include <psdr/psdr.h>
#include <psdr/core/bitmap.h>
#include <psdr/core/cube_distrb.h>
#include "emitter.h"

namespace psdr
{

PSDR_CLASS_DECL_BEGIN(EnvironmentMap, final, Emitter)
public:
    inline EnvironmentMap(const char *file_name) {
        m_radiance.load_openexr(file_name);
    }

    void configure() override;

    inline void set_transform(const Matrix4fD &mat) {
        m_to_world_left = mat;
        m_ready = false;
    }

    SpectrumC eval(const IntersectionC &its, MaskC active = true) const override;
    SpectrumD eval(const IntersectionD &its, MaskD active = true) const override;

    template <bool ad>
    Spectrum<ad> eval_direction(const Vector3f<ad> &wi, Mask<ad> active = true) const;

    PositionSampleC sample_position(const Vector3fC &ref_p, const Vector2fC &sample2, MaskC active = true) const override;
    PositionSampleD sample_position(const Vector3fD &ref_p, const Vector2fD &sample2, MaskD active = true) const override;

    FloatC sample_position_pdf(const Vector3fC &ref_p, const IntersectionC &its, MaskC active = true) const override;
    FloatD sample_position_pdf(const Vector3fD &ref_p, const IntersectionD &its, MaskD active = true) const override;

    std::string to_string() const override;

    Bitmap3fD               m_radiance;     // stored in latitude-longitude format
    FloatD                  m_scale = 1.f;

    Matrix4fD               m_to_world_raw = identity<Matrix4fD>(),
                            m_to_world_left = identity<Matrix4fD>(),
                            m_to_world, m_from_world;

    Vector3fC               m_lower, m_upper;
    // HyperCubeDistribution2f m_cell_distrb;

    CubeDistribution        m_cube_distrb;

    ENOKI_PINNED_OPERATOR_NEW(FloatD)

    template <bool ad>
    PositionSample<ad> __sample_position(const Vector3f<ad> &ref_p, const Vector2f<ad>&, Mask<ad>) const;

    std::pair<Vector3fC, FloatC> sample_direction(Vector2fC&) const;

    template <bool ad>
    Float<ad> __sample_position_pdf(const Vector3f<ad> &, const Intersection<ad> &, Mask<ad>) const;
PSDR_CLASS_DECL_END(EnvironmentMap)

} // namespace psdr
