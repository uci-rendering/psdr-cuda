#include <misc/Exception.h>
#include <psdr/core/ray.h>
#include <psdr/scene/optix.h>
#include <psdr/scene/scene_optix.h>
#include <psdr/shape/mesh.h>

namespace psdr
{

void Intersection_OptiX::reserve(int64_t size) {
    PSDR_ASSERT(size > 0);
    if ( size != m_size ) {
        m_size = size;
        triangle_id = empty<IntC>(size);
        shape_id = empty<IntC>(size);
        uv = empty<Vector2fC>(size);
    }
}


Scene_OptiX::Scene_OptiX() {
    m_accel = nullptr;
}


Scene_OptiX::~Scene_OptiX() {
    if ( m_accel != nullptr ) {
        optix_release(*m_accel);
        delete m_accel;
    }
}


void Scene_OptiX::configure(const std::vector<Mesh *> &meshes) {
    PSDR_ASSERT(!meshes.empty());
    size_t num_meshes = meshes.size();

    if ( m_accel == nullptr ) {
        std::vector<int> face_offset(num_meshes + 1);
        face_offset[0] = 0;
        for ( size_t i = 0; i < num_meshes; ++i )
            face_offset[i + 1] = face_offset[i] + meshes[i]->m_num_faces;
        m_accel = new PathTracerState();
        optix_config(*m_accel, face_offset);
    }

    uint32_t triangle_input_flags[] = { OPTIX_GEOMETRY_FLAG_DISABLE_ANYHIT };
    std::vector<CUdeviceptr> vertex_buffer_ptrs(num_meshes);
    std::vector<OptixBuildInput> build_inputs(num_meshes);
    for ( size_t i = 0; i < num_meshes; ++i ) {
        const Mesh *mesh = meshes[i];

        PSDR_ASSERT(static_cast<int>(slices(mesh->m_vertex_buffer)) == mesh->m_num_vertices*3);
        PSDR_ASSERT(static_cast<int>(slices(mesh->m_face_buffer)) == mesh->m_num_faces*3);

        build_inputs[i].type                                = OPTIX_BUILD_INPUT_TYPE_TRIANGLES;
        build_inputs[i].triangleArray.vertexFormat          = OPTIX_VERTEX_FORMAT_FLOAT3;
        build_inputs[i].triangleArray.vertexStrideInBytes   = 3*sizeof(float);
        build_inputs[i].triangleArray.numVertices           = static_cast<uint32_t>(slices(mesh->m_vertex_buffer)/3);
        vertex_buffer_ptrs[i]                               = reinterpret_cast<CUdeviceptr>(mesh->m_vertex_buffer.data());
        build_inputs[i].triangleArray.vertexBuffers         = &vertex_buffer_ptrs[i];

        build_inputs[i].triangleArray.indexFormat           = OPTIX_INDICES_FORMAT_UNSIGNED_INT3;
        build_inputs[i].triangleArray.numIndexTriplets      = static_cast<uint32_t>(slices(mesh->m_face_buffer)/3);
        build_inputs[i].triangleArray.indexBuffer           = reinterpret_cast<CUdeviceptr>(mesh->m_face_buffer.data());
        build_inputs[i].triangleArray.indexStrideInBytes    = 3*sizeof(int);

        build_inputs[i].triangleArray.flags                 = triangle_input_flags;
        build_inputs[i].triangleArray.numSbtRecords         = 1;
    }
    build_accel(*m_accel, build_inputs);
}


bool Scene_OptiX::is_ready() const {
    return m_accel != nullptr;
}


template <bool ad>
Vector2i<ad> Scene_OptiX::ray_intersect(const Ray<ad> &ray, Mask<ad> &active) const {
    const int m = static_cast<int>(slices(ray.o));
    m_its.reserve(m);

    cuda_eval();

    m_accel->params.ray_o_x         = ray.o.x().data();
    m_accel->params.ray_o_y         = ray.o.y().data();
    m_accel->params.ray_o_z         = ray.o.z().data();

    m_accel->params.ray_d_x         = ray.d.x().data();
    m_accel->params.ray_d_y         = ray.d.y().data();
    m_accel->params.ray_d_z         = ray.d.z().data();

    m_accel->params.ray_tmax        = ray.tmax.data();
    m_accel->params.tri_index       = m_its.triangle_id.data();
    m_accel->params.shape_index     = m_its.shape_id.data();
    m_accel->params.barycentric_u   = m_its.uv.x().data();
    m_accel->params.barycentric_v   = m_its.uv.y().data();

    CUDA_CHECK(
        cudaMemcpyAsync(
            reinterpret_cast<void*>( m_accel->d_params ),
            &m_accel->params, sizeof( Params ),
            cudaMemcpyHostToDevice, m_accel->stream
        )
    );

    OPTIX_CHECK(
        optixLaunch(
            m_accel->pipeline,
            m_accel->stream,
            reinterpret_cast<CUdeviceptr>( m_accel->d_params ),
            sizeof( Params ),
            &m_accel->sbt,
            m,              // launch size
            1,              // launch height
            1               // launch depth
        )
    );

    CUDA_SYNC_CHECK();

    active &= (m_its.shape_id >= 0) && (m_its.triangle_id >= 0);
    return Vector2i<ad>(m_its.shape_id, m_its.triangle_id);
}


// Explicit instantiations
template Vector2iC Scene_OptiX::ray_intersect(const RayC &ray, MaskC &active) const;
template Vector2iD Scene_OptiX::ray_intersect(const RayD &ray, MaskD &active) const;

} // namespace psdr
