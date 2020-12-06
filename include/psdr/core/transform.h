#pragma once

#include <psdr/psdr.h>
#include <enoki/transform.h>

namespace psdr
{

namespace transform
{

/// Create a translation transformation
template <typename Float>
static Matrix<Float, 4> translate(const Array<Float, 3> &v) {
    return enoki::translate<Matrix<Float, 4>>(v);
}

/// Create a scale transformation
template <typename Float>
static Matrix<Float, 4> scale(const Array<Float, 3> &v) {
    return enoki::scale<Matrix<Float, 4>>(v);
}

/// Create a rotation transformation around an arbitrary axis in 3D. The angle is specified in degrees
template <typename Float>
static Matrix<Float, 4> rotate(const Array<Float, 3> &axis, float angle) {
    return enoki::rotate<Matrix<Float, 4>>(axis, deg_to_rad(angle));
}

/** \brief Create a perspective transformation.
 *   (Maps [near, far] to [0, 1])
 *
 *  Projects vectors in camera space onto a plane at z=1:
 *
 *  x_proj = x / z
 *  y_proj = y / z
 *  z_proj = (far * (z - near)) / (z * (far-near))
 *
 *  Camera-space depths are not mapped linearly!
 *
 * \param fov Field of view in degrees
 * \param near Near clipping plane
 * \param far  Far clipping plane
 */
static inline ScalarMatrix4f perspective(float fov, float near_, float far_) {
    float recip = 1.f / (far_ - near_);

    /* Perform a scale so that the field of view is mapped
       to the interval [-1, 1] */
    float tan = enoki::tan(deg_to_rad(fov * .5f)),
          cot = 1.f / tan;

    ScalarMatrix4f trafo = diag<ScalarMatrix4f>(ScalarVector4f(cot, cot, far_ * recip, 0.f));
    trafo(2, 3) = -near_ * far_ * recip;
    trafo(3, 2) = 1.f;

    return trafo;
}

/** \brief Create a look-at camera transformation
 *
 * \param origin Camera position
 * \param target Target vector
 * \param up     Up vector
 */
template <typename Float>
static Matrix<Float, 4> look_at(const Array<Float, 3> &origin, const Array<Float, 3> &target, const Array<Float, 3> &up) {
    Array<Float, 3> dir = normalize(target - origin);
    Array<Float, 3> left = normalize(cross(up, dir));
    Array<Float, 3> new_up = cross(dir, left);

    return Matrix<Float, 4>::from_cols(
        concat(left, 0.0f),
        concat(new_up, 0.0f),
        concat(dir, 0.0f),
        concat(origin, 1.0f)
    );
}


} // namespace transform


template <typename Float>
static Array<Float, 3> transform_pos(const Matrix<Float, 4> &mat, const Array<Float, 3> &vec) {
    Array<Float, 4> tmp = mat*concat(vec, 1.f);
    return head<3>(tmp)/tmp.w();
}


template <typename Float>
static Array<Float, 3> transform_dir(const Matrix<Float, 4> &mat, const Array<Float, 3> &vec) {
    return head<3>(mat*concat(vec, 0.f));
}

} // namespace psdr
