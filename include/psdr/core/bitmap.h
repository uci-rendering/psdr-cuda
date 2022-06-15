#pragma once

#include <psdr/psdr.h>
#include <misc/Exception.h>
#include <psdr/core/intersection.h>

namespace psdr
{

template <int channels>
struct Bitmap {
    static_assert(channels == 1 || channels == 3);
    using ScalarValue = typename std::conditional<channels == 1, float, const ScalarVector3f&>::type;

    template <bool ad>
    using Value = typename std::conditional<channels == 1, Float<ad>, Vector3f<ad>>::type;

    using ValueC = Value<false>;
    using ValueD = Value<true>;

    Bitmap();
    Bitmap(const char *file_name);
    Bitmap(ScalarValue value);
    Bitmap(int width, int height, const ValueD &data);

    void load_openexr(const char *file_name);

    inline void fill(ScalarValue value) {
        m_resolution = ScalarVector2i(1, 1);
        m_data = ValueD(value);
    }

    template <bool ad> Value<ad> eval(Vector2f<ad> uv, bool flip_v = true) const;

    template <bool ad> Value<ad> at(Int<ad> idx) const;

    template <bool ad> Value<ad> sample(const Intersection<ad>& its) const;

    ScalarVector2i m_resolution;
    ValueD m_data;
};

using Bitmap1fD = Bitmap<1>;
using Bitmap3fD = Bitmap<3>;

} // namespace psdr
