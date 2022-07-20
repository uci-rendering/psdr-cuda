#pragma once


#include <psdr/psdr.h>
#include <psdr/core/records.h>
#include <psdr/edge/edge.h>
#include <map>

namespace psdr
{

struct Edge_Graph {
    std::map<int, std::vector<int>>           vertices;
    std::map<std::pair<int,int>, std::pair<bool, int>>        edge_map;

    std::vector<int> edge_id;

    Vector2iD                                 sorted_indices;
    IntD                                      sorted_edge_id;

    std::array<std::vector<float>, 3> m_vertex_positions;
    std::array<std::vector<int  >, 2> m_edge_indices;

    std::map<std::pair<int,int>, int>         m_edge_indices_map;

    std::vector<std::pair<int,int>>           draw;
    std::vector<int>                          jump;

    bool                                      m_ready = false;

    std::vector<double>                        cos_data;
    std::vector<double>                        length_data;

    float                                     cut_thold = 1.0f;
    float                                     global_cut_thold = 1.0f;

    EdgeSortOption      m_edge_sort;


    Edge_Graph() {};
    void config(const Vector3fD& vertex_positions, const Vector2iD& valid_edge_indices, int vertex_size, const EdgeSortOption &sort_option);
    void Greedy(int edge_key, int prev);
    void print();

};


PSDR_CLASS_DECL_BEGIN(Mesh, final, Object)
public:
    Mesh() = default;
    ~Mesh() override;

    void load(const char *fname, bool verbose = false);
    void configure();
    void prepare_optix_buffers();

    inline void set_transform(const Matrix4fD &mat, bool set_left = true) {
        if ( set_left ) {
            m_to_world_left = mat;
        } else {
            m_to_world_right = mat;
        }
        m_ready = false;
    }

    inline void append_transform(const Matrix4fD &mat, bool append_left = true) {
        if ( append_left ) {
            m_to_world_left = mat*m_to_world_left;
        } else {
            m_to_world_right *= mat;
        }
        m_ready = false;
    }

#ifdef PSDR_MESH_ENABLE_1D_VERTEX_OFFSET
    void shift_vertices();
#endif

    PositionSampleC sample_position(const Vector2fC &sample2, MaskC active = true) const;
    PositionSampleD sample_position(const Vector2fD &sample2, MaskD active = true) const;

    FloatC sample_position_pdf(const IntersectionC &its, MaskC active = true) const;
    FloatD sample_position_pdf(const IntersectionD &its, MaskD active = true) const;

    MaskC get_obj_mask(std::string obj_name) const{
        if (m_id == obj_name) {
            return MaskC(1);
        } else {
            return MaskC(0);
        }
    };

    IntC get_obj_id() const{
        return IntC(m_mesh_id);
    };

    void dump(const char *fname) const;

    std::string to_string() const override;

    int                 m_mesh_id = -1;

    bool                m_ready = false;

    bool                m_use_face_normals = false,
                        m_has_uv = false;

    bool                m_enable_edges = true; 

    EdgeSortOption      m_edge_sort;

    // Indicates if the mesh creates primiary and secondary edges

    Matrix4fD           m_to_world_raw   = identity<Matrix4fD>(),
                        m_to_world_left  = identity<Matrix4fD>(),
                        m_to_world_right = identity<Matrix4fD>();

    const BSDF*         m_bsdf = nullptr;
    const Emitter*      m_emitter = nullptr;

    int                 m_num_vertices = 0,
                        m_num_faces = 0;

    IntC                m_cut_position;

    Vector2iD           m_valid_edge_indices;

    Vector3fD           m_vertex_positions_raw,
                        m_vertex_normals_raw;

#ifdef PSDR_MESH_ENABLE_1D_VERTEX_OFFSET
    FloatD              m_vertex_offset;
#endif
    Vector2fD           m_vertex_uv;

    // When PSDR_MESH_ENABLE_1D_VERTEX_OFFSET is defined:
    //   m_vertex_positions = m_to_world * (m_vertex_positions_raw +
    //                                      m_vertex_offset * m_vertex_normals_raw)
    // Otherwise:
    //   m_vertex_positions = m_to_world * m_vertex_positions_raw
    //
    Vector3fD           m_vertex_positions;

    Vector3iD           m_face_indices,
                        m_face_uv_indices;

    // 0, 1: storing vertex indices of the end points
    // 2, 3: storing face indices sharing each edge
    // 4   : storing the third vertex of face with index stored in [2]
    Vectori<5, true>    m_edge_indices;

    float               m_total_area, m_inv_total_area;

    // For position sampling
    DiscreteDistribution *m_face_distrb = nullptr;

    // Temporary triangle info for Scene::configure()
    TriangleInfoD       *m_triangle_info = nullptr;
    TriangleUVD         *m_triangle_uv = nullptr;
    SecondaryEdgeInfo   *m_sec_edge_info = nullptr;

    // For OptiX ray tracing
    FloatC              m_vertex_buffer;
    IntC                m_face_buffer;

    Edge_Graph          edge_graph;

    ENOKI_PINNED_OPERATOR_NEW(FloatD)

protected:
    template <bool ad>
    PositionSample<ad> __sample_position(const Vector2f<ad>&, Mask<ad>) const;
PSDR_CLASS_DECL_END(Mesh)

} // namespace psdr

ENOKI_CALL_SUPPORT_BEGIN(psdr::Mesh)
    ENOKI_CALL_SUPPORT_GETTER(bsdf, m_bsdf)
    ENOKI_CALL_SUPPORT_GETTER(emitter, m_emitter)
    ENOKI_CALL_SUPPORT_METHOD(get_obj_mask)
    ENOKI_CALL_SUPPORT_METHOD(get_obj_id)
    ENOKI_CALL_SUPPORT_METHOD(sample_position)
    ENOKI_CALL_SUPPORT_METHOD(sample_position_pdf)
ENOKI_CALL_SUPPORT_END(psdr::Mesh)
