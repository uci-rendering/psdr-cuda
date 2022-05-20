#pragma once

#include <unordered_map>
#include <psdr/psdr.h>
#include <psdr/edge/edge.h>

namespace psdr
{

PSDR_CLASS_DECL_BEGIN(Scene, final, Object)
    friend class SceneLoader;

public:
    using ParamMap = std::unordered_map<std::string, const Object&>;

    Scene();
    ~Scene() override;

    void load_file(const char *file_name, bool auto_configure = true);
    void load_string(const char *scene_xml, bool auto_configure = true);

    void reload_mesh(Mesh& mesh, const char *file_name, bool verbose=false);
    void reload_mesh_mem(Mesh& mesh, const Vector3fD &vertex_positions, const Vector3iD &face_indices,
                         const Vector2fD &vertex_uv = Vector2fD(), const Vector3iD &face_uv_indices = Vector3iD(),
                         bool verbose=false);

    void configure();
    bool is_ready() const;

    template <bool ad, bool path_space = false>
    Intersection<ad> ray_intersect(const Ray<ad> &ray, Mask<ad> active = true, TriangleInfoD *out_info = nullptr) const;

    template <bool ad>
    Spectrum<ad> Lenv(const Vector3f<ad> &wi, Mask<ad> active = true) const;

    template <bool ad>
    PositionSample<ad> sample_emitter_position(const Vector3f<ad> &ref_p, const Vector2f<ad> &sample, Mask<ad> active = true) const;

    template <bool ad>
    Float<ad> emitter_position_pdf(const Vector3f<ad> &ref_p, const Intersection<ad> &its, Mask<ad> active = true) const;

    BoundarySegSampleDirect sample_boundary_segment_direct(const Vector3fC &sample3, MaskC active = true) const;

    std::string to_string() const override;

    int                     m_num_sensors;
    std::vector<Sensor*>    m_sensors;

    std::vector<Emitter*>   m_emitters;
    EnvironmentMap          *m_emitter_env;
    EmitterArrayD           m_emitters_cuda;
    DiscreteDistribution    *m_emitters_distrb;

    std::vector<BSDF*>      m_bsdfs;

    int                     m_num_meshes;
    std::vector<Mesh*>      m_meshes;
    MeshArrayD              m_meshes_cuda;

    // Scene bounding box
    Vector3fC               m_lower, m_upper;

    ParamMap                m_param_map;

    RenderOption            m_opts;
    mutable Sampler         *m_samplers;

protected:
    TriangleInfoD           m_triangle_info;
    TriangleUVD             m_triangle_uv;
    MaskD                   m_triangle_face_normals;
    bool                    m_has_bound_mesh;

    SecondaryEdgeInfo       m_sec_edge_info;
    DiscreteDistribution    *m_sec_edge_distrb;

    bool                    m_loaded;
    Scene_OptiX             *m_optix;

PSDR_CLASS_DECL_END(Scene)

} // namespace psdr
