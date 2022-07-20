#pragma once

#include <psdr/core/bitmap.h>
#include "bsdf.h"

namespace psdr
{

struct GGXDistribution {
public:
    GGXDistribution() : m_alpha_u(0.1f), m_alpha_v(0.1f) {}
    GGXDistribution(const FloatD &alpha) : m_alpha_u(alpha), m_alpha_v(alpha) {}
    GGXDistribution(const FloatD &alpha_u, const FloatD &alpha_v) : m_alpha_u(alpha_u), m_alpha_v(alpha_v) {}

    template <bool ad>
    Float<ad> G(const Vector3f<ad>& wi, const Vector3f<ad>& wo, const Vector3f<ad>& m) const;

    template <bool ad>
    Float<ad> eval(const Vector3f<ad>& m) const;

    template <bool ad>
    std::pair<Vector3f<ad>, Float<ad>> sample(const Vector3f<ad>& wi, const Vector3f<ad>& sample) const;

    template <bool ad>
    Float<ad> smith_g1(const Vector3f<ad>& v, const Vector3f<ad>& m) const;

    template <bool ad>
    Vector2f<ad> sample_visible_11(const Float<ad>& cos_theta_i, const Vector2f<ad>& sample) const;

    FloatD m_alpha_u, m_alpha_v;
};

} // namespace psdr
