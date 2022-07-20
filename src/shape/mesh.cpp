#include <vector>
#include <map>

#define TINYOBJLOADER_IMPLEMENTATION
#include <tiny_obj_loader/tiny_obj_loader.h>

#include <misc/Exception.h>
#include <psdr/core/transform.h>
#include <psdr/core/pmf.h>
#include <psdr/core/warp.h>
#include <psdr/bsdf/bsdf.h>
#include <psdr/emitter/emitter.h>
#include <psdr/shape/mesh.h>

#include <chrono>
#include <numeric>

namespace psdr
{

template <bool ad>
static std::pair<TriangleInfo<ad>, Vector3f<ad>> process_mesh(const Vector3f<ad> &vertex_positions, const Vector3i<ad> &face_indices) {
    const int num_vertices = static_cast<int>(slices(vertex_positions));

    TriangleInfo<ad> triangles;
    triangles.p0 = gather<Vector3f<ad>>(vertex_positions, face_indices[0]);
    triangles.e1 = gather<Vector3f<ad>>(vertex_positions, face_indices[1]) - triangles.p0;
    triangles.e2 = gather<Vector3f<ad>>(vertex_positions, face_indices[2]) - triangles.p0;

    Vector3f<ad> &face_normals = triangles.face_normal;
    Float<ad>    &face_areas = triangles.face_area;

    face_normals = cross(triangles.e1, triangles.e2);
    face_areas   = norm(face_normals);

    Vector3f<ad> vertex_normals = zero<Vector3f<ad>>(num_vertices);
    Float<ad>    vertex_weights = zero<Float<ad>>(num_vertices);
    for ( int i = 0; i < 3; ++i ) {
        scatter_add(vertex_normals, face_normals, face_indices[i]);
        scatter_add(vertex_weights, face_areas  , face_indices[i]);
    }
    vertex_normals = normalize(vertex_normals/vertex_weights);

    triangles.n0 = gather<Vector3f<ad>>(vertex_normals, face_indices[0]);
    triangles.n1 = gather<Vector3f<ad>>(vertex_normals, face_indices[1]);
    triangles.n2 = gather<Vector3f<ad>>(vertex_normals, face_indices[2]);

    // Normalize the face normals
    face_normals /= face_areas;
    face_areas *= 0.5f;
    
    return { triangles, vertex_normals };
}


Mesh::~Mesh() {
    if ( m_face_distrb ) delete m_face_distrb;
    if ( m_triangle_info ) delete m_triangle_info;
    if ( m_triangle_uv ) delete m_triangle_uv;
    if ( m_sec_edge_info ) delete m_sec_edge_info;
}


void Mesh::load(const char *fname, bool verbose) {
    tinyobj::attrib_t attrib;
    std::vector<tinyobj::shape_t> shapes;
    std::vector<tinyobj::material_t> materials;
    std::string warn, err;

    PSDR_ASSERT_MSG(
        tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, fname),
        std::string("Failed to load OBJ from: ") + fname
    );

    assert(attrib.vertices.size() % 3 == 0);
    m_num_vertices = static_cast<int>(attrib.vertices.size())/3;

    // Loading vertex positions
    {
        std::vector<float> buffers[3];
        buffers[0].resize(m_num_vertices);
        buffers[1].resize(m_num_vertices);
        buffers[2].resize(m_num_vertices);

        for ( int i = 0; i < m_num_vertices; ++i )
            for ( int j = 0; j < 3; ++j )
                buffers[j][i] = attrib.vertices[3*i + j];

        m_vertex_positions_raw = Vector3fD(FloatD::copy(buffers[0].data(), m_num_vertices),
                                           FloatD::copy(buffers[1].data(), m_num_vertices),
                                           FloatD::copy(buffers[2].data(), m_num_vertices));
    }

#ifdef PSDR_MESH_ENABLE_1D_VERTEX_OFFSET
    m_vertex_offset = zero<FloatD>(m_num_vertices);
#endif

    // Loading vertex uv coordinates
    if ( (m_has_uv = !attrib.texcoords.empty()) ) {
        int n = static_cast<int>(attrib.texcoords.size())/2;

        std::vector<float> buffers[2];
        buffers[0].resize(n);
        buffers[1].resize(n);

        for ( int i = 0; i < n; ++i )
            for ( int j = 0; j < 2; ++j )
                buffers[j][i] = attrib.texcoords[2*i + j];

        m_vertex_uv = Vector2fD(FloatD::copy(buffers[0].data(), n),
                                FloatD::copy(buffers[1].data(), n));
    }

    m_num_faces = 0;
    for ( size_t s = 0; s < shapes.size(); s++ )
        m_num_faces += static_cast<int>(shapes[s].mesh.num_face_vertices.size());

    std::vector<int> v[6];
    for ( int i = 0; i < 6; ++i ) v[i].reserve(m_num_faces);

    // Loading face indices
    for ( size_t s = 0; s < shapes.size(); ++s ) {
        for ( size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); ++f ) {
            int fv = shapes[s].mesh.num_face_vertices[f];
            assert(fv == 3);
            for ( int i = 0; i < fv; ++i ) {
                auto idx = shapes[s].mesh.indices[3*f + i];
                v[i].push_back(idx.vertex_index);
                if ( m_has_uv )
                    v[3 + i].push_back(idx.texcoord_index);
            }
        }
    }
    assert(static_cast<int>(v[0].size()) == m_num_faces);

    m_face_indices = Vector3iD(IntD::copy(v[0].data(), m_num_faces),
                               IntD::copy(v[1].data(), m_num_faces),
                               IntD::copy(v[2].data(), m_num_faces));
    if ( m_has_uv ) {
        m_face_uv_indices = Vector3iD(IntD::copy(v[3].data(), m_num_faces),
                                      IntD::copy(v[4].data(), m_num_faces),
                                      IntD::copy(v[5].data(), m_num_faces));
    }

    // Constructing edge list

    int m_num_edges = 0;
    if ( m_enable_edges ) {
        std::vector<int> buffers[5];
        buffers[0].reserve(3*m_num_faces);
        buffers[1].reserve(3*m_num_faces);
        buffers[2].reserve(3*m_num_faces);
        buffers[3].reserve(3*m_num_faces);
        buffers[4].reserve(3*m_num_faces);

        std::map<std::pair<int, int>, std::vector<int>> edge_map;
        for ( size_t s = 0; s < shapes.size(); ++s ) {
            for ( size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); ++f ) {
                int fv = shapes[s].mesh.num_face_vertices[f];
                assert(fv == 3);
                for ( int i = 0; i < fv; ++i ) {
                    int k1 = i, k2 = (i + 1) % 3, k3 = (i + 2) % 3;
                    int idx1 = shapes[s].mesh.indices[3*f + k1].vertex_index;
                    int idx2 = shapes[s].mesh.indices[3*f + k2].vertex_index;
                    int idx3 = shapes[s].mesh.indices[3*f + k3].vertex_index;
                    auto key = idx1 < idx2 ?
                        std::make_pair(idx1, idx2) : std::make_pair(idx2, idx1);
                    if (edge_map.find(key) == edge_map.end()) {
                        auto it = edge_map.insert(edge_map.end(), { key, std::vector<int>() });
                        it->second.push_back(idx3);
                    }
                    edge_map[key].push_back(static_cast<int>(f));
                }
            }
        }

        for ( auto it: edge_map ) {
            if ( it.second.size() > 3 ) {
                // Non-manifold mesh is not allowed
                PSDR_ASSERT_MSG(false, std::string("Edge shared by more than 2 faces: ") + fname);
            } else {
                buffers[0].push_back(it.first.first);
                buffers[1].push_back(it.first.second);
                if ( it.second.size() == 3 ) {
                    PSDR_ASSERT_MSG(it.second[1] != it.second[2], std::string("Duplicated faces: ") + fname);
                    buffers[2].push_back(it.second[1]);
                    buffers[3].push_back(it.second[2]);
                    buffers[4].push_back(it.second[0]);
                    ++m_num_edges;
                } else {
                    PSDR_ASSERT_MSG(it.second.size() == 2, std::string("Edge should be boundary: ") + fname);
                    buffers[2].push_back(it.second[1]);
                    buffers[3].push_back(-1);
                    buffers[4].push_back(it.second[0]);
                    ++m_num_edges;
                }
            }
        }

        m_edge_indices = Vectori<5, true>(IntD::copy(buffers[0].data(), m_num_edges),
                                          IntD::copy(buffers[1].data(), m_num_edges),
                                          IntD::copy(buffers[2].data(), m_num_edges),
                                          IntD::copy(buffers[3].data(), m_num_edges),
                                          IntD::copy(buffers[4].data(), m_num_edges));
    }

    if ( verbose ) {
        std::cout << "Loaded " << m_num_vertices << " vertices, "
                               << m_num_faces    << " faces, "
                               << m_num_edges    << " edges. " << std::endl;
    }

    m_ready = false;
}


void Mesh::configure() {
    if ( m_bsdf != nullptr ) {
        PSDR_ASSERT(!m_bsdf->anisotropic() || !m_use_face_normals);
    }

    // Calculating the "raw" (i.e., object-space) vertex normals
    std::tie(std::ignore, m_vertex_normals_raw) = process_mesh<true>(m_vertex_positions_raw, m_face_indices);

    Matrix4fD to_world = m_to_world_left * m_to_world_raw * m_to_world_right;

    // Calculating the world-space vertex positions
#ifdef PSDR_MESH_ENABLE_1D_VERTEX_OFFSET
    m_vertex_positions = transform_pos(
        to_world, fmadd(m_vertex_normals_raw, m_vertex_offset, m_vertex_positions_raw)
    );
#else
    m_vertex_positions = transform_pos(to_world, m_vertex_positions_raw);
#endif

    m_triangle_info = new TriangleInfoD();
    TriangleInfoD &triangle_info = *m_triangle_info;
    std::tie(triangle_info, std::ignore) = process_mesh<true>(m_vertex_positions, m_face_indices);

    const FloatD &face_areas = triangle_info.face_area;
    m_total_area = hsum(face_areas)[0];
    m_inv_total_area = 1.f/m_total_area;

    if ( m_has_uv ) {
        m_triangle_uv = new TriangleUVD();
        for ( int i = 0; i < 3; ++i )
            (*m_triangle_uv)[i] = gather<Vector2fD>(m_vertex_uv, m_face_uv_indices[i]);
    }

    if ( m_face_distrb == nullptr ) m_face_distrb = new DiscreteDistribution();
    m_face_distrb->init(detach(face_areas));

    if ( m_enable_edges ) {
        if (m_edge_sort.enable_sort == false) {
            m_cut_position = 0;
            if ( m_sec_edge_info == nullptr ) m_sec_edge_info = new SecondaryEdgeInfo();
            SecondaryEdgeInfo secEdgeInfo;
            secEdgeInfo.is_boundary = (m_edge_indices[3] < 0);
            secEdgeInfo.p0 = gather<Vector3fD>(m_vertex_positions, m_edge_indices[0]);
            secEdgeInfo.e1 = gather<Vector3fD>(m_vertex_positions, m_edge_indices[1]) - secEdgeInfo.p0;
            secEdgeInfo.n0 = gather<Vector3fD>(m_triangle_info->face_normal, m_edge_indices[2]);
            secEdgeInfo.n1 = gather<Vector3fD>(m_triangle_info->face_normal, m_edge_indices[3], ~secEdgeInfo.is_boundary);
            secEdgeInfo.p2 = gather<Vector3fD>(m_vertex_positions, m_edge_indices[4]);

            MaskD keep = (dot(secEdgeInfo.n0, secEdgeInfo.n1) < 1.f - EdgeEpsilon);
            *m_sec_edge_info = compressD<SecondaryEdgeInfo>(secEdgeInfo, keep);
        } else {
            if ( m_sec_edge_info == nullptr ) m_sec_edge_info = new SecondaryEdgeInfo();

            Vector3fD n0 = gather<Vector3fD>(m_triangle_info->face_normal, m_edge_indices[2]);
            Vector3fD n1 = gather<Vector3fD>(m_triangle_info->face_normal, m_edge_indices[3], ~(m_edge_indices[3] < 0));

            MaskD keep = (dot(n0, n1) < 1.f - EdgeEpsilon);

            IntD v0 = compressD<IntD>((m_edge_indices[0]), keep);
            IntD v1 = compressD<IntD>((m_edge_indices[1]), keep);
            IntD v2 = compressD<IntD>((m_edge_indices[2]), keep);
            IntD v3 = compressD<IntD>((m_edge_indices[3]), keep);
            IntD v4 = compressD<IntD>((m_edge_indices[4]), keep);

            m_valid_edge_indices[0] = v0;
            m_valid_edge_indices[1] = v1;

            // TODO: should not have this
            // if (edge_graph.m_ready == false) {
                edge_graph.config(m_vertex_positions, m_valid_edge_indices, m_num_vertices, m_edge_sort);
                // edge_graph.m_ready = true;
            // }
            m_cut_position = IntC::copy(edge_graph.jump.data(), edge_graph.jump.size());

            IntD nv0 = edge_graph.sorted_indices[0];
            IntD nv1 = edge_graph.sorted_indices[1];
            IntD nv2 = zero<IntC>(slices(v2));
            IntD nv3 = zero<IntC>(slices(v3));
            IntD nv4 = zero<IntC>(slices(v4));

            scatter(nv2, v2, edge_graph.sorted_edge_id);
            scatter(nv3, v3, edge_graph.sorted_edge_id);
            scatter(nv4, v4, edge_graph.sorted_edge_id);

            m_sec_edge_info->is_boundary = (nv3 < 0);
            m_sec_edge_info->p0 = gather<Vector3fD>(m_vertex_positions, nv0);
            m_sec_edge_info->e1 = gather<Vector3fD>(m_vertex_positions, nv1) - m_sec_edge_info->p0;
            m_sec_edge_info->n1 = gather<Vector3fD>(m_triangle_info->face_normal, nv2);
            m_sec_edge_info->n0 = gather<Vector3fD>(m_triangle_info->face_normal, nv3, ~m_sec_edge_info->is_boundary);
            m_sec_edge_info->p2 = gather<Vector3fD>(m_vertex_positions, nv4);
        }
    } else {
        if ( m_sec_edge_info != nullptr ) {
            delete m_sec_edge_info;
            m_sec_edge_info = nullptr;
        }
    }

    m_ready = true;

    prepare_optix_buffers();
}


void Mesh::prepare_optix_buffers() {
    PSDR_ASSERT(m_ready);
    IntC idx;

    m_vertex_buffer = empty<FloatC>(m_num_vertices*3);
    idx = arange<IntC>(m_num_vertices)*3;
    for ( int i = 0; i < 3; ++i ) {
        scatter(m_vertex_buffer, detach(m_vertex_positions[i]), idx + i);
    }

    m_face_buffer = empty<IntC>(m_num_faces*3);
    idx = arange<IntC>(m_num_faces)*3;
    for ( int i = 0; i < 3; ++i ) {
        scatter(m_face_buffer, detach(m_face_indices[i]), idx + i);
    }
}


PositionSampleC Mesh::sample_position(const Vector2fC &sample2, MaskC active) const {
    return __sample_position<false>(sample2, active);
}


PositionSampleD Mesh::sample_position(const Vector2fD &sample2, MaskD active) const {
    return __sample_position<true>(sample2, active);
}


template <bool ad>
PositionSample<ad> Mesh::__sample_position(const Vector2f<ad> &_sample2, Mask<ad> active) const {
    PSDR_ASSERT(m_ready && m_emitter != nullptr);
    PSDR_ASSERT(m_triangle_info != nullptr);

    PositionSample<ad> result;
    Vector2f<ad> sample2 = _sample2;

    IntC idx;
    std::tie(idx, std::ignore) = m_face_distrb->sample_reuse<ad>(sample2.x());
    sample2 = warp::square_to_uniform_triangle<ad>(sample2);

    TriangleInfo<ad> tri_info;
    if constexpr ( ad ) {
        tri_info = gather<TriangleInfoD>(*m_triangle_info, IntD(idx), active);
        result.J = tri_info.face_area/detach(tri_info.face_area);
    } else {
        tri_info = gather<TriangleInfoC>(detach(*m_triangle_info), idx, active);
        result.J = 1.f;
    }
    result.p = bilinear<ad>(tri_info.p0, tri_info.e1, tri_info.e2, sample2);
    result.n = tri_info.face_normal;
    result.pdf = m_inv_total_area;
    result.is_valid = true;
    return result;
}


FloatC Mesh::sample_position_pdf(const IntersectionC &its, MaskC active) const {
    active &= eq(its.shape, this);
    return FloatC(m_inv_total_area) & active;
}


FloatD Mesh::sample_position_pdf(const IntersectionD &its, MaskD active) const {
    active &= eq(its.shape, this);
    return FloatD(m_inv_total_area) & active;
}


#ifdef PSDR_MESH_ENABLE_1D_VERTEX_OFFSET
void Mesh::shift_vertices() {
    m_vertex_positions_raw += detach(m_vertex_normals_raw)*detach(m_vertex_offset);
    m_vertex_offset = zero<FloatD>(m_num_vertices);
    m_ready = false;
}
#endif


void Mesh::dump(const char *fname) const {
    std::array<std::vector<float>, 3> vertex_positions, vertex_normals;
    {
#ifdef PSDR_MESH_ENABLE_1D_VERTEX_OFFSET
        Vector3fC vertex_positions_ = fmadd(detach(m_vertex_normals_raw), detach(m_vertex_offset), detach(m_vertex_positions_raw));
#else
        const Vector3fC &vertex_positions_ = detach(m_vertex_positions_raw);
#endif
        copy_cuda_array<float, 3>(vertex_positions_, vertex_positions);

        if ( !m_use_face_normals ) {
            Vector3fC vertex_normals_;
            std::tie(std::ignore, vertex_normals_) = process_mesh<false>(vertex_positions_, detach(m_face_indices));
            copy_cuda_array<float, 3>(vertex_normals_, vertex_normals);
        }
    }

    FILE *fout = fopen(fname, "wt");
    for ( int i = 0; i < m_num_vertices; ++i ) {
        fprintf(fout, "v %.6e %.6e %.6e\n", vertex_positions[0][i], vertex_positions[1][i], vertex_positions[2][i]);
        if ( !m_use_face_normals ) {
            fprintf(fout, "vn %.6e %.6e %.6e\n", vertex_normals[0][i], vertex_normals[1][i], vertex_normals[2][i]);
        }
    }

    std::array<std::vector<int32_t>, 3> face_indices;
    copy_cuda_array<int32_t, 3>(detach(m_face_indices), face_indices);
    if ( m_has_uv ) {
        std::array<std::vector<float>, 2> vertex_uv;
        std::array<std::vector<int32_t>, 3> face_uv_indices;
        copy_cuda_array<float, 2>(detach(m_vertex_uv), vertex_uv);
        copy_cuda_array<int32_t, 3>(detach(m_face_uv_indices), face_uv_indices);

        int m = static_cast<int>(slices(m_vertex_uv));
        for ( int i = 0; i < m; ++i )
            fprintf(fout, "vt %.6le %.6le\n", vertex_uv[0][i], vertex_uv[1][i]);

        for ( int i = 0; i < m_num_faces; ++i ) {
            int v0 = face_indices[0][i] + 1, v1 = face_indices[1][i] + 1, v2 = face_indices[2][i] + 1;
            if ( m_use_face_normals ) {
                fprintf(fout, "f %d/%d %d/%d %d/%d\n",
                        v0, face_uv_indices[0][i] + 1,
                        v1, face_uv_indices[1][i] + 1,
                        v2, face_uv_indices[2][i] + 1);

            } else {
                fprintf(fout, "f %d/%d/%d %d/%d/%d %d/%d/%d\n",
                        v0, face_uv_indices[0][i] + 1, v0,
                        v1, face_uv_indices[1][i] + 1, v1,
                        v2, face_uv_indices[2][i] + 1, v2);
            }
        }
    } else {
        for ( int i = 0; i < m_num_faces; ++i ) {
            int v0 = face_indices[0][i] + 1, v1 = face_indices[1][i] + 1, v2 = face_indices[2][i] + 1;
            if ( m_use_face_normals ) {
                fprintf(fout, "f %d %d %d\n", v0, v1, v2);
            } else {
                fprintf(fout, "f %d//%d %d//%d %d//%d\n", v0, v0, v1, v1, v2, v2);
            }
        }
    }

    fclose(fout);
}


std::string Mesh::to_string() const {
    std::stringstream oss;
    oss << "Mesh[nv=" << m_num_vertices << ", nf=" << m_num_faces;
    if ( m_id != "" ) oss << ", id=" << m_id;
    oss << ", bsdf=" << m_bsdf->to_string() << "]";
    return oss.str();
}


void Edge_Graph::print() {
    // std::cout << "total Vertex: " << vertices.size() << std::endl;
    // std::cout << "total edge: " << edge_map.size() << std::endl;
    std::cout << "total Segment: " << jump.size() << std::endl;
}

void Edge_Graph::config(const Vector3fD& vertex_positions, const Vector2iD& valid_edge_indices, int vertex_size, const EdgeSortOption& sort_option) {
    PSDR_ASSERT(sort_option.enable_sort == true);
    edge_id.clear();
    draw.clear();
    jump.clear();
    vertices.clear();
    edge_map.clear();
    m_edge_indices_map.clear();
    cos_data.clear();
    length_data.clear();


    m_edge_sort = sort_option;
    // auto start = std::chrono::high_resolution_clock::now();

    const Vector3fC &m_vertex_positions_ = detach(vertex_positions);
    const Vector2iC &m_valid_edge_indices_ = detach(valid_edge_indices);
    copy_cuda_array<float, 3>(m_vertex_positions_, m_vertex_positions);
    copy_cuda_array<int  , 2>(m_valid_edge_indices_, m_edge_indices);

    int index_size = m_edge_indices[0].size();

    // auto stop = std::chrono::high_resolution_clock::now();
    // std::cout << "Copy from GPU to CPU: " << std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count() / 1000000.0f << " second" << std::endl;
    // start = stop;


    for (int i=0; i<vertex_size; ++i) {
        std::vector<int> temp;
        vertices[i] = temp;
    }

    for (int i=0; i<index_size; ++i) {
        int e1 = m_edge_indices[0][i];
        int e2 = m_edge_indices[1][i];
        vertices[e1].push_back(e2);
        vertices[e2].push_back(e1);
    }

    for (int i=0; i<index_size; ++i) {
        std::pair<int,int> ind_key(m_edge_indices[0][i], m_edge_indices[1][i]);
        m_edge_indices_map[ind_key] = i;
    }

    for (int i=0; i<vertex_size; ++i) {
        for (int j=0; j<vertices[i].size(); ++j) {
            std::pair<int,int> edge_key(i, vertices[i][j]);
            if (edge_key.first > edge_key.second) {
                std::swap(edge_key.first, edge_key.second);
            }
            edge_map[edge_key].first = false;
            std::pair<int,int> sec_key(edge_key.first, edge_key.second);
            edge_map[edge_key].second = m_edge_indices_map[sec_key];
        }
    }

    // stop = std::chrono::high_resolution_clock::now();
    // std::cout << "Building data structure: " << std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count() / 1000000.0f << " second" << std::endl;
    // start = stop;
    cut_thold = cos(sort_option.local_angle * M_PI / 180.f);
    global_cut_thold = cos(sort_option.global_angle * M_PI / 180.f);

    // TODO: start at any node and give up odd vertices?
    for(auto iter = vertices.begin(); iter != vertices.end(); ++iter)
    {
        for (auto &v : iter->second) {
            std::pair<int,int> temp(v, iter->first);
            if (temp.first > temp.second) {
                std::swap(temp.first, temp.second);
            }
            if (edge_map[temp].first == false) {
                Greedy(v, iter->first);
                jump.push_back(draw.size()-1);
            }
        }
    }

    // stop = std::chrono::high_resolution_clock::now();
    // std::cout << "Greedy search edges: " << std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count() / 1000000.0f << " second" << std::endl;
    // start = stop;

    double cos_sum = std::accumulate(cos_data.begin(), cos_data.end(), 0.0);
    double cos_mean = cos_sum / cos_data.size();

    double cos_sq_sum = std::inner_product(cos_data.begin(), cos_data.end(), cos_data.begin(), 0.0);
    double cos_stdev = std::sqrt(cos_sq_sum / cos_data.size() - cos_mean * cos_mean);

    // std::cout << cos_data.size() << " Average edge angle: mean = " << cos_mean << " degree +- " << cos_stdev << std::endl;

    double length_sum = std::accumulate(length_data.begin(), length_data.end(), 0.0);
    double length_mean = length_sum / length_data.size();

    double length_sq_sum = std::inner_product(length_data.begin(), length_data.end(), length_data.begin(), 0.0);
    double length_stdev = std::sqrt(length_sq_sum / length_data.size() - length_mean * length_mean);

    // std::cout << length_data.size() << " Average edge length: mean = " << length_mean << " +- " << length_stdev << std::endl;

    // stop = std::chrono::high_resolution_clock::now();
    // start = stop;

    std::vector<int> ind1;
    std::vector<int> ind2;

    ind1.resize(draw.size());
    ind2.resize(draw.size());

    for (int i=0; i<draw.size(); ++i) {
            ind1[i] = draw[i].first;
            ind2[i] = draw[i].second;
    }

    std::map<std::pair<int,int>, int>        draw_map;
    for (int i=0; i<draw.size(); ++i) {
        std::pair<int, int> draw_temp(ind1[i], ind2[i]);
        if (draw_temp.first > draw_temp.second) {
            std::swap(draw_temp.first, draw_temp.second);
        }
        draw_map[draw_temp] = i;
    }

    sorted_indices = Vector2iD(IntD::copy(ind1.data(), draw.size()),
                               IntD::copy(ind2.data(), draw.size()));
    std::vector<int> sorted_edge;
    for (int i=0; i<edge_id.size(); ++i) {
        std::pair<int, int> edge_temp_map(m_edge_indices[0][i], m_edge_indices[1][i]);
        sorted_edge.push_back(draw_map[edge_temp_map]);
    }
    sorted_edge_id = IntD::copy(sorted_edge.data(), sorted_edge.size());

    // stop = std::chrono::high_resolution_clock::now();
    // std::cout << "Copy to GPU: " << std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count() / 1000000.0f << " second" << std::endl;
    // start = stop;

    print();
}

void Edge_Graph::Greedy(int init_curr, int init_prev) {
    int curr = init_curr;
    int prev = init_prev;
    int steps = 0;
    float best_cos = -1.0f;
    float best_global_cos = -1.0f;

    float p0_x_init = m_vertex_positions[0][init_prev];
    float p0_y_init = m_vertex_positions[1][init_prev];
    float p0_z_init = m_vertex_positions[2][init_prev];

    while(true) {
        draw.push_back(std::pair<int,int>(prev, curr));
        std::pair<int,int> temp(prev, curr);
        if (temp.first > temp.second) {
            std::swap(temp.first, temp.second);
        }
        // PSDR_ASSERT(edge_map[temp].first == false);
        edge_map[temp].first = true;
        edge_id.push_back(edge_map[temp].second);
        steps++;
        if (steps >= m_edge_sort.max_depth) {
            break;
        }

        // use a loop to find the best cos:
        int best_id = -1;
        int best_global_id = -1;

        best_cos = 1.1f;
        best_global_cos = 1.1f;


        for (int i=0; i<vertices[curr].size(); ++i) {
            std::pair<int,int> temp2(curr, vertices[curr][i]);
            if (temp2.first > temp2.second) {
                std::swap(temp2.first, temp2.second);
            }
            if (edge_map[temp2].first == false) {


                float p0_x = m_vertex_positions[0][prev];
                float p0_y = m_vertex_positions[1][prev];
                float p0_z = m_vertex_positions[2][prev];

                float p1_x = m_vertex_positions[0][curr];
                float p1_y = m_vertex_positions[1][curr];
                float p1_z = m_vertex_positions[2][curr];
                float p2_x = m_vertex_positions[0][vertices[curr][i]];
                float p2_y = m_vertex_positions[1][vertices[curr][i]];
                float p2_z = m_vertex_positions[2][vertices[curr][i]];

                float e1_x = p0_x - p1_x;
                float e1_y = p0_y - p1_y;
                float e1_z = p0_z - p1_z;
                float e1_norm = sqrt(e1_x*e1_x + e1_y*e1_y + e1_z*e1_z);
                float e2_x = p2_x - p1_x;
                float e2_y = p2_y - p1_y;
                float e2_z = p2_z - p1_z;
                float e2_norm = sqrt(e2_x*e2_x + e2_y*e2_y + e2_z*e2_z);

                float cos_val = (e1_x*e2_x + e1_y*e2_y + e1_z*e2_z) / (e1_norm * e2_norm);

                e1_x = p0_x_init - p1_x;
                e1_y = p0_y_init - p1_y;
                e1_z = p0_z_init - p1_z;
                e1_norm = sqrt(e1_x*e1_x + e1_y*e1_y + e1_z*e1_z);
                e2_x = p2_x - p1_x;
                e2_y = p2_y - p1_y;
                e2_z = p2_z - p1_z;
                e2_norm = sqrt(e2_x*e2_x + e2_y*e2_y + e2_z*e2_z);

                float cos_global_val = (e1_x*e2_x + e1_y*e2_y + e1_z*e2_z) / (e1_norm * e2_norm);

                if (cos_val < best_cos) {
                    best_cos = cos_val;
                    best_id = i;
                }

                if (cos_global_val < best_global_cos) {
                    best_global_cos = cos_global_val;
                    best_global_id = i;
                }

            }
        }

        if  (best_global_id != -1) {
            std::pair<int,int> temp2(curr, vertices[curr][best_global_id]);
            if (temp2.first > temp2.second) {
                std::swap(temp2.first, temp2.second);
            }
            if (edge_map[temp2].first == false) {
                prev = curr;
                curr = vertices[prev][best_global_id];
            } else {
                break;
            }
        } else {
            break;
        }

        if (best_cos > cut_thold || (best_global_cos > global_cut_thold && steps > m_edge_sort.min_global_step)) {
            break;
        } else {
            cos_data.push_back(best_global_cos);
        }
        
    }

    length_data.push_back(static_cast<double>(steps));
}

} // namespace psdr
