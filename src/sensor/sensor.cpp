#include <misc/Exception.h>
#include <psdr/core/ray.h>
#include <psdr/sensor/sensor.h>

namespace psdr
{

void Sensor::configure() {
    m_aspect = static_cast<float>(m_resolution.x())/m_resolution.y();
    PSDR_ASSERT_MSG(std::abs(det(Matrix3fD(m_to_world))[0] - 1.f) < Epsilon,
                    "Sensor transformation should not involve scaling!");
}

} // namespace psdr
