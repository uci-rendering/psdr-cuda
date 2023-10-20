#include <chrono>
#include <misc/Exception.h>

#include <psdr/core/ray.h>
#include <psdr/core/intersection.h>
#include <psdr/core/sampler.h>
#include <psdr/core/pmf.h>
#include <psdr/bsdf/bsdf.h>
#include <psdr/emitter/envmap.h>
#include <psdr/sensor/perspective.h>
#include <psdr/shape/mesh.h>

#include <psdr/scene/scene_optix.h>
#include <psdr/scene/scene_loader.h>
#include <psdr/scene/scene.h>

namespace psdr
{

Scene::Scene() {
    m_has_bound_mesh = false;
    m_loaded = false;
    m_samplers = new Sampler[3];
    m_optix = new Scene_OptiX();
    m_emitter_env = nullptr;
    m_emitters_distrb = new DiscreteDistribution();
    m_sec_edge_distrb = new DiscreteDistribution();
}


Scene::~Scene() {
    for ( Sensor*  s : m_sensors  ) delete s;
    for ( Emitter* e : m_emitters ) delete e;
    for ( BSDF*    b : m_bsdfs    ) delete b;
    for ( Mesh*    m : m_meshes   ) delete m;

    delete[]    m_samplers;
    delete      m_optix;
    delete      m_emitters_distrb;
    delete      m_sec_edge_distrb;
}


void Scene::load_file(const char *file_name, bool auto_configure) {
    SceneLoader::load_from_file(file_name, *this);
    if ( auto_configure ) configure();
}


void Scene::load_string(const char *scene_xml, bool auto_configure) {
    SceneLoader::load_from_string(scene_xml, *this);
    if ( auto_configure ) configure();
}


void Scene::reload_mesh(Mesh& mesh, const char *file_name, bool verbose) {
    mesh.load(file_name, verbose);
    if ( m_has_bound_mesh ) {
        Mesh* bound_mesh = m_meshes.back();
        delete bound_mesh;
        m_meshes.pop_back();
        m_num_meshes--;
        m_has_bound_mesh = false;
    }
    delete m_optix;
    m_optix = new Scene_OptiX();
}


void Scene::reload_mesh_mem(Mesh &mesh, const Vector3fD &vertex_positions, const Vector3iD &face_indices,
                            const Vector2fD &vertex_uv, const Vector3iD &face_uv_indices, bool verbose) {
    mesh.load_mem(vertex_positions, face_indices, vertex_uv, face_uv_indices, verbose);
    if ( m_has_bound_mesh ) {
        Mesh* bound_mesh = m_meshes.back();
        delete bound_mesh;
        m_meshes.pop_back();
        m_num_meshes--;
        m_has_bound_mesh = false;
    }
    delete m_optix;
    m_optix = new Scene_OptiX();
}


void Scene::configure() {
    PSDR_ASSERT_MSG(m_loaded, "Scene not loaded yet!");
    PSDR_ASSERT(m_num_sensors == static_cast<int>(m_sensors.size()));
    PSDR_ASSERT(m_num_meshes == static_cast<int>(m_meshes.size()));

    using namespace std::chrono;
    auto start_time = high_resolution_clock::now();

    // Seed samplers
    if ( m_opts.spp  > 0 ) {
        int64_t sample_count = static_cast<int64_t>(m_opts.height)*m_opts.width*m_opts.spp;
        if ( m_samplers[0].m_sample_count != sample_count )
            m_samplers[0].seed(arange<UInt64C>(sample_count));
    }
    if ( m_opts.sppe > 0 ) {
        int64_t sample_count = static_cast<int64_t>(m_opts.height)*m_opts.width*m_opts.sppe;
        if ( m_samplers[1].m_sample_count != sample_count )
            m_samplers[1].seed(arange<UInt64C>(sample_count));
    }
    if ( m_opts.sppse > 0 ) {
        int64_t sample_count = static_cast<int64_t>(m_opts.height)*m_opts.width*m_opts.sppse;
        if ( m_samplers[2].m_sample_count != sample_count )
            m_samplers[2].seed(arange<UInt64C>(sample_count));
    }

    // Preprocess meshes
    PSDR_ASSERT_MSG(!m_meshes.empty(), "Missing meshes!");
    std::vector<int> face_offset, edge_offset;
    face_offset.reserve(m_num_meshes + 1);
    face_offset.push_back(0);
    edge_offset.reserve(m_num_meshes + 1);
    edge_offset.push_back(0);
    m_lower = full<Vector3fC>(std::numeric_limits<float>::max());
    m_upper = full<Vector3fC>(std::numeric_limits<float>::min());
    for ( Mesh *mesh : m_meshes ) {
        mesh->configure();
        face_offset.push_back(face_offset.back() + mesh->m_num_faces);
        if ( m_opts.sppse > 0 && mesh->m_enable_edges )
            edge_offset.push_back(edge_offset.back() + static_cast<int>(slices(*(mesh->m_sec_edge_info))));
        else
            edge_offset.push_back(edge_offset.back());
        for ( int i = 0; i < 3; ++i ) {
            m_lower[i] = enoki::min(m_lower[i], hmin(detach(mesh->m_vertex_positions[i])));
            m_upper[i] = enoki::max(m_upper[i], hmax(detach(mesh->m_vertex_positions[i])));
        }
    }

    // Preprocess sensors
    PSDR_ASSERT_MSG(!m_sensors.empty(), "Missing sensor!");
    std::vector<size_t> num_edges;
    if ( m_opts.sppe > 0 ) num_edges.reserve(m_sensors.size());
    for ( Sensor *sensor : m_sensors ) {
        sensor->m_resolution = ScalarVector2i(m_opts.width, m_opts.height);
        sensor->m_scene = this;
        sensor->configure();
        if ( m_opts.sppe > 0 ) num_edges.push_back(sensor->m_edge_distrb.m_size);

        for ( int i = 0; i < 3; ++i ) {
            if ( PerspectiveCamera *camera = dynamic_cast<PerspectiveCamera *>(sensor) ) {
                m_lower[i] = enoki::min(m_lower[i], detach(camera->m_camera_pos[i]));
                m_upper[i] = enoki::max(m_upper[i], detach(camera->m_camera_pos[i]));
            }
        }
    }
    if ( m_opts.log_level > 0 ) {
        std::stringstream oss;
        oss << "AABB: [lower = " << m_lower << ", upper = " << m_upper << "]";
        log(oss.str().c_str());

        if ( m_opts.sppe > 0 ) {
            std::stringstream oss;
            oss << "(" << num_edges[0];
            for ( size_t i = 1; i < num_edges.size(); ++i ) oss << ", " << num_edges[i];
            oss << ") primary edges initialized.";
            log(oss.str().c_str());
        }
    }

    // Handling env. lighting
    if ( m_emitter_env != nullptr && !m_has_bound_mesh ) {
        FloatC margin = hmin((m_upper - m_lower)*0.05f);
        m_lower -= margin; m_upper += margin;

        m_emitter_env->m_lower = m_lower;
        m_emitter_env->m_upper = m_upper;

        // Adding a bounding mesh
        std::array<std::vector<float>, 3> lower_, upper_;
        copy_cuda_array<float, 3>(m_lower, lower_);
        copy_cuda_array<float, 3>(m_upper, upper_);

        float vtx_data[3][8];
        for ( int i = 0; i < 8; ++i )
            for ( int j = 0; j < 3; ++j )
                vtx_data[j][i] = (i & (1 << j)) ? upper_[j][0] : lower_[j][0];

        const int face_data[3][12] = {
            0, 0, 1, 1, 2, 2, 0, 0, 0, 0, 4, 4,
            1, 3, 5, 7, 3, 7, 5, 4, 2, 6, 7, 6,
            3, 2, 7, 3, 7, 6, 1, 5, 6, 4, 5, 7
        };

        Mesh *bound_mesh = new Mesh();
        bound_mesh->m_num_vertices = 8;
        bound_mesh->m_num_faces = 12;
        bound_mesh->m_use_face_normals = true;
        bound_mesh->m_enable_edges = false;
        bound_mesh->m_bsdf = nullptr;
        bound_mesh->m_emitter = m_emitter_env;
        for ( int i = 0; i < 3; ++i ) {
            bound_mesh->m_vertex_positions_raw[i] = FloatD::copy(vtx_data[i], 8);
            bound_mesh->m_face_indices[i] = IntD::copy(face_data[i], 12);
        }
        bound_mesh->configure();
        m_meshes.push_back(bound_mesh);

        face_offset.push_back(face_offset.back() + 12);
        edge_offset.push_back(edge_offset.back());

        ++m_num_meshes;
        m_has_bound_mesh = true;
        if ( m_opts.log_level > 0 ) {
            log("Bounding mesh added for environmental lighting.");
        }
    }

    // Preprocess emitters
    if ( !m_emitters.empty() ) {
        std::vector<float> weights;
        weights.reserve(m_emitters.size());
        for ( Emitter *e : m_emitters ) {
            e->configure();
            weights.push_back(e->m_sampling_weight);
        }
        m_emitters_distrb->init(FloatC::copy(weights.data(), weights.size()));

        float inv_total_weight = rcp(m_emitters_distrb->m_sum)[0];
        for ( Emitter *e : m_emitters ) {
            e->m_sampling_weight *= inv_total_weight;
        }
    }

    // Initialize CUDA arrays
    if ( !m_emitters.empty() ) {
        m_emitters_cuda = EmitterArrayD::copy(m_emitters.data(), m_emitters.size());
    }
    m_meshes_cuda = MeshArrayD::copy(m_meshes.data(), m_num_meshes);

    // Generate global triangle arrays
    m_triangle_info = empty<TriangleInfoD>(face_offset.back());
    m_triangle_uv = zero<TriangleUVD>(face_offset.back());
    m_triangle_face_normals = empty<MaskD>(face_offset.back());
    for ( int i = 0; i < m_num_meshes; ++i ) {
        const Mesh &mesh = *m_meshes[i];
        const IntD idx = arange<IntD>(mesh.m_num_faces) + face_offset[i];
        scatter(m_triangle_info, *mesh.m_triangle_info, idx);
        scatter(m_triangle_face_normals, MaskD(mesh.m_use_face_normals), idx);
        if ( mesh.m_has_uv ) {
            scatter(m_triangle_uv, *mesh.m_triangle_uv, idx);
        }
    }

    // Generate global sec. edge arrays
    if ( m_opts.sppse > 0 ) {
        m_sec_edge_info = empty<SecondaryEdgeInfo>(edge_offset.back());
        for ( int i = 0; i < m_num_meshes; ++i ) {
            const Mesh &mesh = *m_meshes[i];
            if ( mesh.m_enable_edges ) {
                PSDR_ASSERT(mesh.m_sec_edge_info != nullptr);
                const int m = static_cast<int>(slices(*mesh.m_sec_edge_info));
                const IntD idx = arange<IntD>(m) + edge_offset[i];
                scatter(m_sec_edge_info, *mesh.m_sec_edge_info, idx);
            }
        }
#if 0
        FloatC angle = select(detach(m_sec_edge_info.is_boundary), Pi,
                              safe_acos(dot(detach(m_sec_edge_info.n0), detach(m_sec_edge_info.n1))));
        m_sec_edge_distrb->init(norm(detach(m_sec_edge_info.e1))*angle);
#else
        m_sec_edge_distrb->init(norm(detach(m_sec_edge_info.e1)));
#endif
        if ( m_opts.log_level > 0 ) {
            std::stringstream oss;
            oss << edge_offset.back() << " secondary edges initialized.";
            log(oss.str().c_str());
        }
    } else {
        m_sec_edge_info = empty<SecondaryEdgeInfo>();
    }

    // Initialize OptiX
    cuda_eval();
    m_optix->configure(m_meshes);

    // Cleanup
    for ( int i = 0; i < m_num_meshes; ++i ) {
        Mesh *mesh = m_meshes[i];

        if ( mesh->m_emitter == nullptr ) {
            PSDR_ASSERT(mesh->m_triangle_info != nullptr);
            delete mesh->m_triangle_info;
            mesh->m_triangle_info = nullptr;

            if ( mesh->m_triangle_uv != nullptr ) {
                delete mesh->m_triangle_uv;
                mesh->m_triangle_uv = nullptr;
            }
        }

        if ( mesh->m_enable_edges ) {
            PSDR_ASSERT(mesh->m_sec_edge_info != nullptr);
            delete mesh->m_sec_edge_info;
            mesh->m_sec_edge_info = nullptr;
        }
    }

    auto end_time = high_resolution_clock::now();
    if ( m_opts.log_level > 0 ) {
        std::stringstream oss;
        oss << "Configured in " << duration_cast<duration<double>>(end_time - start_time).count() << " seconds.";
        log(oss.str().c_str());
    }
}


bool Scene::is_ready() const {
    return (m_opts.spp   == 0 || m_samplers[0].is_ready()) &&
           (m_opts.sppe  == 0 || m_samplers[1].is_ready()) &&
           (m_opts.sppse == 0 || m_samplers[2].is_ready()) &&
           m_optix->is_ready();
}


template <bool ad, bool path_space>
Intersection<ad> Scene::ray_intersect(const Ray<ad> &ray, Mask<ad> active, TriangleInfoD *out_info) const {
    static_assert(ad || !path_space);

    Intersection<ad> its;
    Vector2i<ad> idx = m_optix->ray_intersect<ad>(ray, active);

    TriangleInfo<ad>    tri_info;
    TriangleUV<ad>      tri_uv_info;
    Mask<ad>            face_normal_mask;
    if constexpr ( ad ) {
        tri_info = gather<TriangleInfoD>(m_triangle_info, idx[1], active);
        if ( out_info != nullptr ) *out_info = tri_info;

        its.v0_idx = tri_info.v0_idx;
        its.v1_idx = tri_info.v1_idx;
        its.v2_idx = tri_info.v2_idx;

        tri_uv_info = gather<TriangleUVD>(m_triangle_uv, idx[1], active);
        face_normal_mask = gather<MaskD>(m_triangle_face_normals, idx[1], active);
        if constexpr ( path_space ) {
            its.J = tri_info.face_area/detach(tri_info.face_area);
        } else {
            its.J = 1.f;
        }
    } else {
        if ( out_info != nullptr ) {
            *out_info = gather<TriangleInfoD>(m_triangle_info, IntD(idx[1]), active);
            tri_info = detach(*out_info);
        } else {
            tri_info = gather<TriangleInfoC>(detach(m_triangle_info), idx[1], active);
        }

        its.v0_idx = tri_info.v0_idx;
        its.v1_idx = tri_info.v1_idx;
        its.v2_idx = tri_info.v2_idx;

        tri_uv_info = gather<TriangleUVC>(detach(m_triangle_uv), idx[1], active);
        face_normal_mask = gather<MaskC>(detach(m_triangle_face_normals), idx[1], active);
        its.J = 1.f;
    }
    its.n = tri_info.face_normal;

    const Vector3f<ad> &vertex0 = tri_info.p0, &edge1 = tri_info.e1, &edge2 = tri_info.e2;

    if constexpr ( !ad || path_space ) {
        // Path-space formulation
        const Vector2fC &uv = m_optix->m_its.uv;

        its.barycentric_uv = uv;

        //Vector3f<ad> sh_n = normalize(fmadd(tri_info.n0, 1.f - u - v, fmadd(tri_info.n1, u, tri_info.n2*v)));
        Vector3f<ad> sh_n = normalize(bilinear<ad>(tri_info.n0,
                                                   tri_info.n1 - tri_info.n0,
                                                   tri_info.n2 - tri_info.n0,
                                                   uv));
        masked(sh_n, face_normal_mask) = its.n;

        if constexpr ( ad ) {
            its.shape = gather<MeshArrayD>(m_meshes_cuda, idx[0], active);
        } else {
            its.shape = gather<MeshArrayC>(detach(m_meshes_cuda), idx[0], active);
        }
        its.p = bilinear<ad>(vertex0, edge1, edge2, uv);

        Vector3f<ad> dir = its.p - ray.o;
        its.t = norm(dir);
        dir /= its.t;

        its.sh_frame = Frame<ad>(sh_n);
        its.wi = its.sh_frame.to_local(-dir);
        //its.uv = fmadd(tri_uv_info[0], 1.f - u - v, fmadd(tri_uv_info[1], u, tri_uv_info[2]*v));
        its.uv = bilinear2<ad>(tri_uv_info[0],
                               tri_uv_info[1] - tri_uv_info[0],
                               tri_uv_info[2] - tri_uv_info[0],
                               uv);
    } else {
        // Standard (solid-angle) formulation
        auto [uv, t] = ray_intersect_triangle<true>(vertex0, edge1, edge2, ray);

        its.barycentric_uv = uv;

        //Vector3f<ad> sh_n = normalize(fmadd(tri_info.n0, 1.f - u_d - v_d, fmadd(tri_info.n1, u_d, tri_info.n2*v_d)));
        Vector3fD sh_n = normalize(bilinear<true>(tri_info.n0,
                                                  tri_info.n1 - tri_info.n0,
                                                  tri_info.n2 - tri_info.n0,
                                                  uv));
        masked(sh_n, face_normal_mask) = its.n;

        its.shape = gather<MeshArrayD>(m_meshes_cuda, idx[0], active);
        its.p = ray(t);
        its.t = t;

        its.sh_frame = Frame<ad>(sh_n);
        its.wi = its.sh_frame.to_local(-ray.d);
        //its.uv = fmadd(tri_uv_info[0], 1.f - u_d - v_d, fmadd(tri_uv_info[1], u_d, tri_uv_info[2]*v_d));
        its.uv = bilinear2<true>(tri_uv_info[0],
                                 tri_uv_info[1] - tri_uv_info[0],
                                 tri_uv_info[2] - tri_uv_info[0],
                                 uv);

        // int num_samples = static_cast<int>(slices(ray.o));
        // FloatC tmp = zero<FloatC>(num_samples);
        // masked(tmp, detach(active)) = norm(Vector2fC(detach(u_d), detach(v_d)) - m_optix->m_its.uv);
        // std::cout << hmax(tmp) << std::endl;
    }
    return its;
}


template <bool ad>
Spectrum<ad> Scene::Lenv(const Vector3f<ad> &wi, Mask<ad> active) const {
    return m_emitter_env == nullptr ? 0.f : m_emitter_env->eval_direction<ad>(wi, active);
}


template <typename T>
static void print_to_string(const std::vector<T*>& arr, const char* name, std::stringstream& oss) {
    if ( !arr.empty() ) {
        oss << "  # " << name << "\n";
        for ( size_t i = 0; i < arr.size(); ++i )
            oss << "  " << arr[i]->to_string() << "\n";
    }
}


template <typename T>
static void print_to_string(const Type<T*, true> &arr, const char *name, std::stringstream &oss) {
    if ( !arr.empty() ) {
        oss << "  # " << name << "\n";
        for ( size_t i = 0; i < slices(arr); ++i )
            oss << "  " << arr[i]->to_string() << "\n";
    }
}


std::string Scene::to_string() const {
    std::stringstream oss;
    oss << "Scene[\n";
    print_to_string<Sensor>(m_sensors, "Sensors", oss);
    oss << "\n";
    print_to_string<BSDF  >(m_bsdfs  , "BSDFs"  , oss);
    oss << "\n";
    print_to_string<Mesh  >(m_meshes , "Meshes" , oss);
    oss << "]";
    return oss.str();
}


template <bool ad>
PositionSample<ad> Scene::sample_emitter_position(const Vector3f<ad> &ref_p, const Vector2f<ad> &_sample2, Mask<ad> active) const {
    PSDR_ASSERT_MSG(!m_emitters.empty(), "No Emitter!");

    PositionSample<ad> result;
    if ( m_emitters.size() == 1U ) {
        result = m_emitters[0]->sample_position(ref_p, _sample2, active);
    } else {
        Vector2f<ad> sample2 = _sample2;
        auto [emitter_index, emitter_pdf] = m_emitters_distrb->sample_reuse<ad>(sample2.y());

        EmitterArray<ad> emitter_arr;
        if constexpr ( ad ) {
            emitter_arr = gather<EmitterArrayD>(m_emitters_cuda, IntD(emitter_index), active);
        } else {
            emitter_arr = gather<EmitterArrayC>(detach(m_emitters_cuda), emitter_index, active);
        }
        result = emitter_arr->sample_position(ref_p, sample2, active);
        result.pdf *= emitter_pdf;
    }
    return result;
}


template <bool ad>
Float<ad> Scene::emitter_position_pdf(const Vector3f<ad> &ref_p, const Intersection<ad> &its, Mask<ad> active) const {
    return its.shape->emitter(active)->sample_position_pdf(ref_p, its, active);
}


BoundarySegSampleDirect Scene::sample_boundary_segment_direct(const Vector3fC &sample3, MaskC active) const {
    BoundarySegSampleDirect result;

    // Sample a point p0 on a face edge

    FloatC sample1 = sample3.x();
    auto [edge_idx, pdf0] = m_sec_edge_distrb->sample_reuse<false>(sample1);

    SecondaryEdgeInfo info = gather<SecondaryEdgeInfo>(m_sec_edge_info, IntD(edge_idx), active);
    result.p0 = fmadd(info.e1, sample1, info.p0); // p0 = info.p0 + info.e1*sample1;
    result.edge = normalize(detach(info.e1));
    result.edge2 = detach(info.p2) - detach(info.p0);
    const Vector3fC &p0 = detach(result.p0);
    pdf0 /= norm(detach(info.e1));

    // Sample a point ps2 on a emitter

    PositionSampleC ps2 = sample_emitter_position<false>(p0, tail<2>(sample3), active);
    result.p2 = ps2.p;
    result.n = ps2.n;

    // Construct the edge "ray" and check if it is valid

    Vector3fC e = result.p2 - p0;
    const FloatC distSqr = squared_norm(e);
    e /= safe_sqrt(distSqr);
    const FloatC cosTheta = dot(result.n, -e);

    IntC sgn0 = sign<false>(dot(detach(info.n0), e), EdgeEpsilon),
         sgn1 = sign<false>(dot(detach(info.n1), e), EdgeEpsilon);
    result.is_valid = active && (cosTheta > Epsilon) && (
        (detach(info.is_boundary) && neq(sgn0, 0)) || (~detach(info.is_boundary) && (sgn0*sgn1 < 0))
    );

    result.pdf = (pdf0*ps2.pdf*(distSqr/cosTheta)) & result.is_valid;
    return result;
}

// Explicit instantiations
template IntersectionC Scene::ray_intersect<false, false>(const RayC&, MaskC, TriangleInfoD*) const;
template IntersectionD Scene::ray_intersect<true , false>(const RayD&, MaskD, TriangleInfoD*) const;
template IntersectionD Scene::ray_intersect<true , true >(const RayD&, MaskD, TriangleInfoD*) const;

template SpectrumC Scene::Lenv<false>(const Vector3fC&, MaskC) const;
template SpectrumD Scene::Lenv<true >(const Vector3fD&, MaskD) const;

template PositionSampleC Scene::sample_emitter_position<false>(const Vector3fC&, const Vector2fC&, MaskC) const;
template PositionSampleD Scene::sample_emitter_position<true >(const Vector3fD&, const Vector2fD&, MaskD) const;

template FloatC Scene::emitter_position_pdf<false>(const Vector3fC&, const IntersectionC&, MaskC) const;
template FloatD Scene::emitter_position_pdf<true >(const Vector3fD&, const IntersectionD&, MaskD) const;

} // namespace psdr
