#pragma once

#include <psdr/psdr.h>

namespace psdr
{


template <typename Vector3f> std::pair<Vector3f, Vector3f> coordinate_system(const Vector3f &n) {
    static_assert(Vector3f::Size == 3, "coordinate_system() expects a 3D vector as input!");
    using Float = value_t<Vector3f>;

    /* Based on "Building an Orthonormal Basis, Revisited" by
       Tom Duff, James Burgess, Per Christensen,
       Christophe Hery, Andrew Kensler, Max Liani,
       and Ryusuke Villemin (JCGT Vol 6, No 1, 2017) */

    Float sign = enoki::sign(n.z()),
          a    = -rcp(sign + n.z()),
          b    = n.x() * n.y() * a;

    return {
        Vector3f(mulsign(sqr(n.x()) * a, n.z()) + 1.f,
                 mulsign(b, n.z()),
                 mulsign_neg(n.x(), n.z())),
        Vector3f(b, sign + sqr(n.y()) * a, -n.y())
    };
}


template <typename Float_>
struct Frame_
{
    static constexpr bool ad = std::is_same_v<Float_, FloatD>;

    inline Frame_(const Vector3f<ad> &v) : n(v) {
        std::tie(s, t) = coordinate_system(v);
    }

    inline Frame_(const Vector3f<ad> &s, const Vector3f<ad> &t, const Vector3f<ad> &n) : s(s), t(t), n(n) {}

    inline Frame_(const Frame_ &frame) : s(frame.s), t(frame.t), n(frame.n) {}

    /// Convert from world coordinates to local coordinates
    Vector3f<ad> to_local(const Vector3f<ad> &v) const {
        return Vector3f<ad>(dot(v, s), dot(v, t), dot(v, n));
    }

    /// Convert from local coordinates to world coordinates
    Vector3f<ad> to_world(const Vector3f<ad> &v) const {
        return s*v.x() + t*v.y() + n*v.z();
    }

    /** \brief Give a unit direction, this function returns the cosine of the
     * elevation angle in a reference spherical coordinate system (see the \ref
     * Frame description)
     */
    static Float<ad> cos_theta(const Vector3f<ad> &v) { return v.z(); }

    /** \brief Give a unit direction, this function returns the square cosine
     * of the elevation angle in a reference spherical coordinate system (see
     * the \ref Frame description)
     */
    static Float<ad> cos_theta_2(const Vector3f<ad> &v) { return sqr(v.z()); }

    /** \brief Give a unit direction, this function returns the sine
     * of the elevation angle in a reference spherical coordinate system (see
     * the \ref Frame description)
     */
    static Float<ad> sin_theta(const Vector3f<ad> &v) { return safe_sqrt(sin_theta_2(v)); }

    /** \brief Give a unit direction, this function returns the square sine
     * of the elevation angle in a reference spherical coordinate system (see
     * the \ref Frame description)
     */
    static Float<ad> sin_theta_2(const Vector3f<ad> &v) { return fmadd(v.x(), v.x(), sqr(v.y())); }

    /** \brief Give a unit direction, this function returns the tangent
     * of the elevation angle in a reference spherical coordinate system (see
     * the \ref Frame description)
     */
    static Float<ad> tan_theta(const Vector3f<ad> &v) {
        Float<ad> temp = fnmadd(v.z(), v.z(), 1.f);
        return safe_sqrt(temp) / v.z();
    }

    /** \brief Give a unit direction, this function returns the square tangent
     * of the elevation angle in a reference spherical coordinate system (see
     * the \ref Frame description)
     */
    static Float<ad> tan_theta_2(const Vector3f<ad> &v) {
        Float<ad> temp = fnmadd(v.z(), v.z(), 1.f);
        return max(temp, 0.f) / sqr(v.z());
    }

    /** \brief Give a unit direction, this function returns the sine of the
     * azimuth in a reference spherical coordinate system (see the \ref Frame
     * description)
     */
    static Float<ad> sin_phi(const Vector3f<ad> &v) {
        Float<ad> sin_theta_2 = Frame_::sin_theta_2(v),
               inv_sin_theta = rsqrt(Frame_::sin_theta_2(v));
        return select(abs(sin_theta_2) <= 4.f * Epsilon, 0.f,
                      clamp(v.y() * inv_sin_theta, -1.f, 1.f));
    }

    /** \brief Give a unit direction, this function returns the cosine of the
     * azimuth in a reference spherical coordinate system (see the \ref Frame
     * description)
     */
    static Float<ad> cos_phi(const Vector3f<ad> &v) {
        Float<ad> sin_theta_2 = Frame_::sin_theta_2(v),
               inv_sin_theta = rsqrt(Frame_::sin_theta_2(v));
        return select(abs(sin_theta_2) <= 4.f * Epsilon, 1.f,
                      clamp(v.x() * inv_sin_theta, -1.f, 1.f));
    }

    /** \brief Give a unit direction, this function returns the sine and cosine
     * of the azimuth in a reference spherical coordinate system (see the \ref
     * Frame description)
     */
    static std::pair<Float<ad>, Float<ad>> sincos_phi(const Vector3f<ad> &v) {
        Float<ad> sin_theta_2 = Frame_::sin_theta_2(v),
               inv_sin_theta = rsqrt(Frame_::sin_theta_2(v));

        Vector2fD result = head<2>(v) * inv_sin_theta;

        result = select(abs(sin_theta_2) <= 4.f * Epsilon,
                        Vector2fD(1.f, 0.f),
                        clamp(result, -1.f, 1.f));

        return { result.y(), result.x() };
    }

    /** \brief Give a unit direction, this function returns the squared sine of
     * the azimuth in a reference spherical coordinate system (see the \ref
     * Frame description)
     */
    static Float<ad> sin_phi_2(const Vector3f<ad> &v) {
        Float<ad> sin_theta_2 = Frame_::sin_theta_2(v);
        return select(abs(sin_theta_2) <= 4.f * Epsilon, 0.f,
                      clamp(sqr(v.y()) / sin_theta_2, -1.f, 1.f));
    }

    /** \brief Give a unit direction, this function returns the squared cosine of
     * the azimuth in a reference spherical coordinate system (see the \ref
     * Frame description)
     */
    static Float<ad> cos_phi_2(const Vector3f<ad> &v) {
        Float<ad> sin_theta_2 = Frame_::sin_theta_2(v);
        return select(abs(sin_theta_2) <= 4.f * Epsilon, 1.f,
                      clamp(sqr(v.x()) / sin_theta_2, -1.f, 1.f));
    }

    /** \brief Give a unit direction, this function returns the squared sine
     * and cosine of the azimuth in a reference spherical coordinate system
     * (see the \ref Frame description)
     */
    static std::pair<Float<ad>, Float<ad>> sincos_phi_2(const Vector3f<ad> &v) {
        Float<ad> sin_theta_2 = Frame_::sin_theta_2(v),
               inv_sin_theta_2 = rcp(sin_theta_2);

        Vector2fD result = sqr(head<2>(v)) * inv_sin_theta_2;

        result = select(abs(sin_theta_2) <= 4.f * Epsilon,
                        Vector2fD(1.f, 0.f), clamp(result, -1.f, 1.f));

        return { result.y(), result.x() };
    }

    /// Equality test
    Mask<ad> operator==(const Frame_ &frame) const {
        return all(eq(frame.s, s) && eq(frame.t, t) && eq(frame.n, n));
    }

    /// Inequality test
    Mask<ad> operator!=(const Frame_ &frame) const {
        return any(neq(frame.s, s) || neq(frame.t, t) || neq(frame.n, n));
    }

    Vector3f<ad> s, t, n;

    ENOKI_STRUCT(Frame_, s, t, n)
};

//inline FrameC detach(const FrameD &frame) {
//    return FrameC(detach(frame.s), detach(frame.t), detach(frame.n));
//}

} // namespace psdr

ENOKI_STRUCT_SUPPORT(psdr::Frame_, s, t, n)
