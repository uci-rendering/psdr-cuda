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

template <bool ad>
std::tuple<Float<ad>, Float<ad>, Float<ad>, Float<ad>> fresnel_dielectric(const Float<ad> &eta, const Float<ad> &cos_theta_i) {
    auto outside_mask = cos_theta_i >= 0.f;
    Float<ad> rcp_eta = rcp(eta),
                 eta_it = select(outside_mask, eta, rcp_eta),
                 eta_ti = select(outside_mask, rcp_eta, eta);

    Float<ad> cos_theta_t_sqr =
        fnmadd(fnmadd(cos_theta_i, cos_theta_i, 1.f), eta_ti * eta_ti, 1.f);

    Float<ad> cos_theta_i_abs = abs(cos_theta_i);
    Float<ad> cos_theta_t_abs = safe_sqrt(cos_theta_t_sqr);

    auto index_matched = eq(eta, 1.f),
         special_case  = index_matched || eq(cos_theta_i_abs, 0.f);

    Float<ad> r_sc = select(index_matched, Float<ad>(0.f), Float<ad>(1.f));

    Float<ad> a_s = fnmadd(eta_it, cos_theta_t_abs, cos_theta_i_abs) /
                       fmadd(eta_it, cos_theta_t_abs, cos_theta_i_abs);

    Float<ad> a_p = fnmadd(eta_it, cos_theta_i_abs, cos_theta_t_abs) /
                       fmadd(eta_it, cos_theta_i_abs, cos_theta_t_abs);

    Float<ad> r = .5f * (sqr(a_s) + sqr(a_p));

    masked(r, special_case) = r_sc;

    Float<ad> cos_theta_t = mulsign_neg(cos_theta_t_abs, cos_theta_i);
    return { r, cos_theta_t, eta_it, eta_ti };
}

ENOKI_INLINE Vector3fC squareToEdgeRayDirection_NB(const Vector2fC &sample, const Vector3fC &n0, const Vector3fC &n1, FloatC &pdf) {
    FloatC tmp = dot(n0, n1);
    FloatC phi0 = acos(tmp);
    pdf = 1.0f/(4.0f*phi0);

    Vector3fC Z = normalize(n0 + n1);
    Vector3fC Y = normalize(cross(n0, Z));
    Vector3fC X = cross(Y, Z);

    FloatC phi = (sample[0] - 0.5f)*phi0;

    phi = clamp(phi, -0.5f*phi0 + Epsilon, 0.5f*phi0 - Epsilon);

    Vector3fC X1 = X*cos(phi) + Z*sin(phi);

    FloatC b = select(sample[1] > 0.5f, 4.0f*sample[1] - 3.0f, 4.0f*sample[1] - 1.0f);
    FloatC a = select(sample[1] > 0.5f, -sqrt(1.0f - b*b), sqrt(1.0f - b*b));

    Vector3fC ret = normalize(X1*a + Y*b);

    return ret;
}

ENOKI_INLINE std::pair<Vector3fC, Vector3fC> coordinate_system2(const Vector3fC &n) {
    FloatC sign = enoki::sign(n.z()),
          a    = -rcp(sign + n.z()),
          b    = n.x() * n.y() * a;

    return {
        Vector3fC(mulsign(sqr(n.x()) * a, n.z()) + 1.f,
                 mulsign(b, n.z()),
                 mulsign_neg(n.x(), n.z())),
        Vector3fC(b, sign + sqr(n.y()) * a, -n.y())
    };
}

ENOKI_INLINE Vector3fC squareToEdgeRayDirection_B(const Vector2fC &sample, const Vector3fC &n0, const Vector3fC &e0, const Vector3fC &e1) {
    Vector3fC tmp = cross(n0, e0 - e1);
    FloatC tmpNorm = norm(tmp);

    Vector3fC s, t, n;

    n = tmp/tmpNorm;
    std::tie(s, t) = coordinate_system2(n);

    FloatC z = 1.0f - 2.0f*sample.y();

    // Avoid sampling directions that are almost within the plane
    z = min(z, FloatC(1.0f - EdgeEpsilon));

    FloatC r = sqrt(max(FloatC(0.0f), FloatC(1.0f - z*z)));

    FloatC Phi = FloatC(2.0f*M_PI)*sample[0];
    FloatC sinPhi = sin(Phi);
    FloatC cosPhi = cos(Phi);
    Vector3fC ret = s*r*cosPhi + t*r*sinPhi + n*z;

    return ret;
}

template <bool ad>
ENOKI_INLINE Float<ad> mis_weight(const Float<ad> &pdf1, const Float<ad> &pdf2) {
    Float<ad> w1 = sqr(pdf1), w2 = sqr(pdf2);
    return w1/(w1 + w2);
}

template <bool ad>
ENOKI_INLINE Float<ad> mis_weight2(const Float<ad> &pdf1, const Float<ad> &pdf2) {
    Float<ad> w1 = pdf1, w2 = pdf2;
    return w1/(w1 + w2);
}


} // namespace psdr
