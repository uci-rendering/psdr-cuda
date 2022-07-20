#pragma once

#include <limits>

namespace psdr
{

constexpr float Epsilon         = 1e-5f;
constexpr float RayEpsilon      = 1e-3f;
constexpr float ShadowEpsilon   = 1e-3f;
constexpr float EdgeEpsilon     = 1e-5f;
constexpr float AQEpsilon       = 1e-5f;
constexpr float DeepEpsilon     = 1e-8f;

constexpr float E               = 2.71828182845904523536f;
constexpr float Pi              = 3.14159265358979323846f;
constexpr float TwoPi           = 6.28318530717958647692f;
constexpr float InvPi           = 0.31830988618379067154f;
constexpr float InvTwoPi        = 0.15915494309189533577f;
constexpr float InvFourPi       = 0.07957747154594766788f;
constexpr float SqrtPi          = 1.77245385090551602793f;
constexpr float InvSqrtPi       = 0.56418958354775628695f;
constexpr float SqrtTwo         = 1.41421356237309504880f;
constexpr float InvSqrtTwo      = 0.70710678118654752440f;
constexpr float SqrtTwoPi       = 2.50662827463100050242f;
constexpr float InvSqrtTwoPi    = 0.39894228040143267794f;

constexpr float Infinity = std::numeric_limits<float>::infinity();

} // namespace psdr
