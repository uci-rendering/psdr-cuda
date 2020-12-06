#pragma once

#include <psdr/psdr.h>
#include "frame.h"

namespace psdr
{

namespace warp
{

/// Low-distortion concentric square to disk mapping by Peter Shirley
template <bool ad>
inline Vector2f<ad> square_to_uniform_disk_concentric(const Vector2f<ad> &sample) {
    Float<ad> x = fmsub(2.f, sample.x(), 1.f),
              y = fmsub(2.f, sample.y(), 1.f);

    /* Modified concentric map code with less branching (by Dave Cline), see
       http://psgraphics.blogspot.ch/2011/01/improved-code-for-concentric-map.html

      Original non-vectorized version:

        Value phi, r;
        if (x == 0 && y == 0) {
            r = phi = 0;
        } else if (x * x > y * y) {
            r = x;
            phi = (math::Pi / 4.f) * (y / x);
        } else {
            r = y;
            phi = (math::Pi / 2.f) - (x / y) * (math::Pi / 4.f);
        }
    */

    Mask<ad> is_zero         = eq(x, zero<Float<ad>>()) &&
                               eq(y, zero<Float<ad>>()),
             quadrant_1_or_3 = abs(x) < abs(y);

    Float<ad> r  = select(quadrant_1_or_3, y, x),
              rp = select(quadrant_1_or_3, x, y);

    Float<ad> phi = .25f * Pi * rp / r;
    masked(phi, quadrant_1_or_3) = .5f * Pi - phi;
    masked(phi, is_zero) = zero<Float<ad>>();

    auto [s, c] = sincos(phi);
    return { r * c, r * s };
}


/// Sample a cosine-weighted vector on the unit hemisphere with respect to solid angles
template <bool ad>
inline Vector3f<ad> square_to_cosine_hemisphere(const Vector2f<ad> &sample) {
    // Low-distortion warping technique based on concentric disk mapping
    Vector2f<ad> p = square_to_uniform_disk_concentric<ad>(sample);

    // Guard against numerical imprecisions
    Float<ad> z = safe_sqrt(1.f - squared_norm(p));

    return { p.x(), p.y(), z };
}


/// Density of \ref square_to_cosine_hemisphere() with respect to solid angles
template <bool ad, bool TestDomain = false>
inline Float<ad> square_to_cosine_hemisphere_pdf(const Vector3f<ad> &v) {
    if constexpr (TestDomain)
        return select(abs(squared_norm(v) - 1.f) > RayEpsilon ||
                      v.z() < 0.f, zero<Float<ad>>(), InvPi * v.z());
    else
        return InvPi * v.z();
}


/// Convert an uniformly distributed square sample into barycentric coordinates
template <bool ad>
inline Vector2f<ad> square_to_uniform_triangle(const Vector2f<ad> &sample) {
    Float<ad> t = safe_sqrt(1.f - sample.x());
    return { 1.f - t, t * sample.y() };
}

/// Density of \ref square_to_uniform_triangle per unit area.
template <bool ad, bool TestDomain = false>
inline Float<ad> square_to_uniform_triangle_pdf(const Vector2f<ad> &p) {
    if constexpr (TestDomain) {
        return select(
            p.x() < zero<Float<ad>>() || p.y() < zero<Float<ad>>()
                                  || (p.x() + p.y() > 1.f),
            zero<Float<ad>>(),
            2.f
        );
    } else {
        return 2.f;
    }
}


} // namespace warp

} // namespace psdr
