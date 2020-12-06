#pragma once

#include <psdr/psdr.h>

struct PathTracerState;

namespace psdr
{

struct Intersection_OptiX {
    void reserve(int64_t size);

    int64_t m_size = 0;
    IntC triangle_id;
    IntC shape_id;
    Vector2fC uv;
};


class Scene_OptiX {
    friend class Scene;
    //friend std::unique_ptr<Scene_OptiX> std::make_unique<Scene_OptiX>();

public:
    ~Scene_OptiX();

protected:
    Scene_OptiX();
    void configure(const std::vector<Mesh*> &meshes);
    bool is_ready() const;

    template <bool ad>
    Vector2i<ad> ray_intersect(const Ray<ad> &ray, Mask<ad> &active) const;

    PathTracerState                 *m_accel;
    mutable Intersection_OptiX      m_its;
};

} // namespace psdr
