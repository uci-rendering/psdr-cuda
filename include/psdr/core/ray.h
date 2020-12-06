#pragma once

#include <psdr/psdr.h>

namespace psdr
{

template <bool ad>
struct Ray {
    inline Ray() = default;

    inline Ray(const Vector3f<ad> &o, const Vector3f<ad> &d, const Float<ad> &tmax) : o(o), d(d), tmax(tmax) {
        assert(slices(o) == slices(d) && slices(o) == slices(tmax));
    }

    inline Ray(const Vector3f<ad> &o, const Vector3f<ad> &d) : o(o), d(d) {
        assert(slices(o) == slices(d));
        tmax = full<Float<ad>>(Infinity, slices(o));
    }

    template <bool ad1> inline Ray(const Ray<ad1> &ray) : o(ray.o), d(ray.d), tmax(ray.tmax) {}

    Vector3f<ad> o, d;
    Float<ad> tmax;

    Vector3f<ad> operator() (const Float<ad> &t) const { return fmadd(d, t, o); }

    inline Ray reversed() const { return Ray(o, -d, tmax); }
};


inline RayC detach(const RayD &ray) {
    return RayC(detach(ray.o), detach(ray.d), detach(ray.tmax));
}

} // namespace psdr
