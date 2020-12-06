#pragma once

#include <psdr/psdr.h>
#include <psdr/edge/edge.h>
#include <psdr/core/pmf.h>
#include <psdr/core/records.h>

namespace psdr
{

template <typename Float_>
struct SensorDirectSample_ : public SampleRecord_<Float_> {
    PSDR_IMPORT_BASE(SampleRecord_<Float_>, ad, pdf, is_valid)

    Vector2f<ad>    q;
    Int<ad>         pixel_idx;
    Float<ad>       sensor_val;

    ENOKI_DERIVED_STRUCT(SensorDirectSample_, Base,
        ENOKI_BASE_FIELDS(pdf, is_valid),
        ENOKI_DERIVED_FIELDS(q, pixel_idx, sensor_val)
    )
};

PSDR_CLASS_DECL_BEGIN(Sensor,, Object)
public:
    virtual ~Sensor() override {}

    virtual void configure();

    virtual RayC sample_primary_ray(const Vector2fC &samples) const = 0;
    virtual RayD sample_primary_ray(const Vector2fD &samples) const = 0;

    virtual SensorDirectSampleC sample_direct(const Vector3fC &p) const = 0;

    virtual PrimaryEdgeSample sample_primary_edge(const FloatC &sample1) const = 0;

    ScalarVector2i          m_resolution;
    float                   m_aspect;

    Matrix4fD               m_to_world = identity<Matrix4fD>();

    const Scene             *m_scene = nullptr;

    // Properties for primary edge sampling

    bool                    m_enable_edges = false;
    PrimaryEdgeInfo         m_edge_info;
    DiscreteDistribution    m_edge_distrb;
PSDR_CLASS_DECL_END(Sensor)

} // namespace psdr

ENOKI_STRUCT_SUPPORT(psdr::SensorDirectSample_, pdf, is_valid, q, pixel_idx, sensor_val)

ENOKI_CALL_SUPPORT_BEGIN(psdr::Sensor)
    ENOKI_CALL_SUPPORT_METHOD(sample_primary_ray)
    ENOKI_CALL_SUPPORT_METHOD(sample_primary_edge)
ENOKI_CALL_SUPPORT_END(psdr::Sensor)
