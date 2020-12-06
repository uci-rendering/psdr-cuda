#pragma once

namespace psdr
{
	
template <typename ArrayD>
ENOKI_INLINE ArrayD compressD(const ArrayD &array, const MaskD &active) {
    PSDR_ASSERT(slices(array) == slices(active));
    IntD idx = IntD(compress(arange<IntC>(slices(array)), detach(active)));
    return gather<ArrayD>(array, idx);
}


template <typename T, int n, bool async = false>
ENOKI_INLINE void copy_cuda_array(const Array<CUDAArray<T>, n> &cuda_array, std::array<std::vector<T>, n> &cpu_array) {
    if constexpr ( !async ) {
        cuda_eval();
    }
    size_t m = slices(cuda_array);
    for ( int i = 0; i < n; ++i ) {
        cpu_array[i].resize(m);
        if constexpr ( async ) {
            cuda_memcpy_from_device_async(cpu_array[i].data(), cuda_array[i].data(), m*sizeof(T));
        } else {
            cuda_memcpy_from_device(cpu_array[i].data(), cuda_array[i].data(), m*sizeof(T));
        }
    }
}


template <bool ad>
ENOKI_INLINE Int<ad> sign(const Float<ad> &x, float eps) {
    Int<ad> result = zero<Int<ad>>(slices(x));
    masked(result, x > eps) = 1;
    masked(result, x < -eps) = -1;
    return result;
}


template <bool ad>
ENOKI_INLINE Vector3f<ad> sphdir(const Float<ad> &theta, const Float<ad> &phi) {
    auto [sin_theta, cos_theta] = sincos(theta);
    auto [sin_phi,   cos_phi]   = sincos(phi);
    return Vector3f<ad>(cos_phi*sin_theta, sin_phi*sin_theta, cos_theta);
}


template <bool ad>
ENOKI_INLINE Vector3f<ad> bilinear(const Vector3f<ad> &p0, const Vector3f<ad> &e1, const Vector3f<ad> &e2, const Vector2f<ad> &st) {
    return fmadd(e1, st.x(), fmadd(e2, st.y(), p0));
}


template <bool ad>
ENOKI_INLINE Vector2f<ad> bilinear2(const Vector2f<ad> &p0, const Vector2f<ad> &e1, const Vector2f<ad> &e2, const Vector2f<ad> &st) {
    return fmadd(e1, st.x(), fmadd(e2, st.y(), p0));
}


template <bool ad>
ENOKI_INLINE Float<ad> rgb2luminance(const Vector3f<ad> &rgb) {
    return rgb.x()*.2126f + rgb.y()*.7152f + rgb.z()*.0722f;
}


template <bool ad>
ENOKI_INLINE auto ray_intersect_triangle(const Vector3f<ad> &p0, const Vector3f<ad> &e1, const Vector3f<ad> &e2, const Ray<ad> &ray) {
    Vector3f<ad>    h = cross(ray.d, e2);
    Float<ad>       a = dot(e1, h);
    Float<ad>       f = rcp(a); // 1.0f/a_d;
    Vector3f<ad>    s = ray.o - p0;
    Float<ad>       u = f*dot(s, h);
    Vector3f<ad>    q = cross(s, e1);
    Float<ad>       v = f*dot(ray.d, q);
    Float<ad>       t = f*dot(e2, q);
    return std::make_pair(Vector2f<ad>(u, v), t);
}


template <int ndim, bool ad>
ENOKI_INLINE auto argmin(const Vectorf<ndim, ad> &vector) {
    Float<ad> result = vector[0];
    Int<ad> idx = zero<Int<ad>>(slices(vector));
    for ( int i = 1; i < ndim; ++i ) {
        Mask<ad> less = vector[i] < result;
        masked(result, less) = vector[i];
        masked(idx, less) = i;
    }
    return std::make_pair(result, idx);
}


template <int ndim, bool ad>
ENOKI_INLINE std::pair<Float<ad>, Int<ad>> argmax(const Vectorf<ndim, ad> &vector) {
    Float<ad> result = vector[0];
    Int<ad> idx = zero<Int<ad>>(slices(vector));
    for ( int i = 1; i < ndim; ++i ) {
        Mask<ad> greater = vector[i] > result;
        masked(result, greater) = vector[i];
        masked(idx, greater) = i;
    }
    return { result, idx };
}


template <bool ad>
ENOKI_INLINE auto ray_intersect_box(const Ray<ad> &ray, const Vector3f<ad> &lower, const Vector3f<ad> &upper) {
    /* First, ensure that the ray either has a nonzero slope on each axis,
       or that its origin on a zero-valued axis is within the box bounds */
    Mask<ad> active = all(neq(ray.d, zero<Vector3f<ad>>()) || ((ray.o > lower) || (ray.o < upper)));

    // Compute intersection intervals for each axis
    Vector3f<ad> t1 = (lower - ray.o)/ray.d,
                 t2 = (upper - ray.o)/ray.d;

    // Ensure proper ordering
    Vector3f<ad> t1p = enoki::min(t1, t2),
                 t2p = enoki::max(t1, t2);

    // Intersect intervals
    Float<ad> mint = hmax(t1p), maxt = hmin(t2p);
    active &= (maxt >= mint);

    return std::make_tuple(active, mint, maxt);
}


template <bool ad>
ENOKI_INLINE auto ray_intersect_scene_aabb(const Ray<ad> &ray, const Vector3f<ad> &lower, const Vector3f<ad> &upper) {
    // Compute intersection intervals for each axis
    Vector3f<ad> t1 = (lower - ray.o)/ray.d,
                 t2 = (upper - ray.o)/ray.d;

    // Ensure proper ordering
    Vector3f<ad> t2p = enoki::max(t1, t2);

    // Intersect intervals
    auto [t, idx] = argmin<3, ad>(t2p);
    Vector3f<ad> n = zero<Vector3f<ad>>(slices(ray));
    for ( int i = 0; i < 3; ++i ) {
        masked(n[i], eq(idx, i)) = -enoki::sign(ray.d[i]);
    }
    Float<ad> G = dot(n, -ray.d)*rcp(sqr(t));
    return std::make_tuple(t, n, G);
}


template <bool ad>
Spectrum<ad> fresnel(const Spectrum<ad> &eta_r, const Spectrum<ad> &eta_i, const Float<ad> &cos_theta_i) {
    // Modified from "Optics" by K.D. Moeller, University Science Books, 1988
    Float<ad> cos_theta_i_2 = sqr(cos_theta_i),
              sin_theta_i_2 = 1.f - cos_theta_i_2,
              sin_theta_i_4 = sqr(sin_theta_i_2);
    Spectrum<ad>     temp_1 = sqr(eta_r) - sqr(eta_i) - sin_theta_i_2,
                   a_2_pb_2 = safe_sqrt(sqr(temp_1) + 4.f * sqr(eta_i * eta_r)),
                          a = safe_sqrt(.5f * (a_2_pb_2 + temp_1));
    Spectrum<ad>     term_1 = a_2_pb_2 + cos_theta_i_2,
                     term_2 = 2.f * cos_theta_i * a;
    Spectrum<ad>        r_s = (term_1 - term_2) / (term_1 + term_2);
    Spectrum<ad>     term_3 = a_2_pb_2 * cos_theta_i_2 + sin_theta_i_4,
                     term_4 = term_2 * sin_theta_i_2;
    Spectrum<ad>        r_p = r_s * (term_3 - term_4) / (term_3 + term_4);
    return .5f * (r_s + r_p);
}

} // namespace psdr
