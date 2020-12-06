#pragma once

#include <psdr/core/frame.h>
#include "sensor.h"

namespace psdr
{

PSDR_CLASS_DECL_BEGIN(PerspectiveCamera, final, Sensor)
public:
    PerspectiveCamera(float fov_x, float near, float far) : m_fov_x(fov_x), m_near_clip(near), m_far_clip(far) {}

    void configure() override;

    RayC sample_primary_ray(const Vector2fC &samples) const override;
    RayD sample_primary_ray(const Vector2fD &samples) const override;

    SensorDirectSampleC sample_direct(const Vector3fC &p) const override;

    PrimaryEdgeSample sample_primary_edge(const FloatC &sample1) const override;

    std::string to_string() const override;

    float       m_fov_x,
                m_near_clip,
                m_far_clip;

    Matrix4fD   m_sample_to_camera, m_camera_to_sample,
                m_world_to_sample, m_sample_to_world;

    Vector3fD   m_camera_pos, m_camera_dir;
    FloatD      m_inv_area;
PSDR_CLASS_DECL_END(PerspectiveCamera)

} // namespace psdr
