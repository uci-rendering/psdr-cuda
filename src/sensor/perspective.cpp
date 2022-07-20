#include <misc/Exception.h>
#include <psdr/core/ray.h>
#include <psdr/core/transform.h>
#include <psdr/shape/mesh.h>
#include <psdr/scene/scene.h>
#include <psdr/sensor/perspective.h>

namespace psdr
{

void PerspectiveCamera::configure() {
    Sensor::configure();

    ScalarMatrix4f camera_to_sample =
        transform::scale(ScalarVector3f(-0.5f, -0.5f * m_aspect, 1.f)) *
        transform::translate(ScalarVector3f(-1.f, -1.f / m_aspect, 0.f)) *
        transform::perspective(m_fov_x, m_near_clip, m_far_clip);

    m_camera_to_sample = Matrix4fD(camera_to_sample);
    m_sample_to_camera = Matrix4fD(inverse(camera_to_sample));

    Matrix4fD m_to_world = m_to_world_left * m_to_world_raw * m_to_world_right;

    m_world_to_sample = m_camera_to_sample*inverse(m_to_world);
    m_sample_to_world = m_to_world*m_sample_to_camera;

    m_camera_pos = transform_pos(m_to_world, zero<Vector3fD>());
    m_camera_dir = transform_dir(m_to_world, Vector3fD(0.f, 0.f, 1.f));

    Vector3fD v00 = transform_pos(m_sample_to_camera, Vector3fD(0.f, 0.f, 0.f)), //  bottom-left corner of the image at the near plane
              v10 = transform_pos(m_sample_to_camera, Vector3fD(1.f, 0.f, 0.f)), // bottom-right corner of the image at the near plane
              v11 = transform_pos(m_sample_to_camera, Vector3fD(1.f, 1.f, 0.f)), //    top-right corner of the image at the near plane
              vc  = transform_pos(m_sample_to_camera, Vector3fD(.5f, .5f, 0.f)); //              center of the image at the near plane
    m_inv_area = rcp(norm(v00 - v10)*norm(v11 - v10))*squared_norm(vc);
    //std::cout << m_inv_pixel_area << std::endl;

    // Construct list of primary edges

    PSDR_ASSERT_MSG(m_scene, "Missing scene data!");

    if ( m_scene->m_opts.sppe > 0 ) {
        const std::vector<Mesh*> &meshes = m_scene->m_meshes;
        std::vector<Vectori<5, true>> per_mesh_info(meshes.size());
        std::vector<int> offset(1, 0);

        for ( size_t i = 0; i < meshes.size(); ++i ) {
            const Mesh *mesh = meshes[i];
            if ( mesh->m_enable_edges ) {
                auto &info = per_mesh_info[i];

                PSDR_ASSERT(all(mesh->m_edge_indices[2] >= 0));
                MaskD valid = (mesh->m_edge_indices[3] >= 0);

                Vector3fD e0 = normalize(m_camera_pos - gather<Vector3fD>(mesh->m_triangle_info->p0, mesh->m_edge_indices[2])),
                          e1 = normalize(m_camera_pos - gather<Vector3fD>(mesh->m_triangle_info->p0, mesh->m_edge_indices[3], valid)),
                          n0 = gather<Vector3fD>(mesh->m_triangle_info->face_normal, mesh->m_edge_indices[2]),
                          n1 = gather<Vector3fD>(mesh->m_triangle_info->face_normal, mesh->m_edge_indices[3], valid);

                // add UV edge
                MaskD uv_mask;

                if (mesh->m_has_uv) {
                    Vector3iD fuv1 = gather<Vector3iD>(mesh->m_face_uv_indices, mesh->m_edge_indices[2]);
                    Vector3iD fuv2 = gather<Vector3iD>(mesh->m_face_uv_indices, mesh->m_edge_indices[3]);
                    IntC uv_cut = zero<IntC>(slices(valid));
                    uv_cut = select( (eq(fuv1[0], fuv2[0]) || eq(fuv1[0], fuv2[1]) || eq(fuv1[0], fuv2[2])), IntC(1), 0);
                    uv_cut = select( (eq(fuv1[1], fuv2[0]) || eq(fuv1[1], fuv2[1]) || eq(fuv1[1], fuv2[2])), uv_cut+1, uv_cut);
                    uv_cut = select( (eq(fuv1[2], fuv2[0]) || eq(fuv1[2], fuv2[1]) || eq(fuv1[2], fuv2[2])), uv_cut+1, uv_cut);

                    uv_mask = neq(uv_cut, 2);
                    PSDR_ASSERT(any(neq(uv_cut, 3)));
                }

                if ( mesh->m_use_face_normals ) {
                    MaskD skip = valid;
                    skip &= (dot(e0, n0) < Epsilon && dot(e1, n1) < Epsilon) ||
                            (dot(n0, n1) > 1.f - Epsilon);
                    if (mesh->m_has_uv) {
                        info = compressD<Vectori<5, true>>(mesh->m_edge_indices, ~skip || uv_mask);
                    } else {
                        info = compressD<Vectori<5, true>>(mesh->m_edge_indices, ~skip);
                    }
                } else {
                    MaskD active = ~valid;
                    active |= (dot(e0, n0) > Epsilon)^(dot(e1, n1) > Epsilon);
                    if (mesh->m_has_uv) {
                        info = compressD<Vectori<5, true>>(mesh->m_edge_indices, active || uv_mask);
                    } else {
                        info = compressD<Vectori<5, true>>(mesh->m_edge_indices, active);
                    }
                }
                PSDR_ASSERT(slices(info) > 0);

                // std::cout << info << std::endl;
                offset.push_back(offset.back() + static_cast<int>(slices(info)));
            } else {
                offset.push_back(offset.back());
            }
        }

        if ( offset.back() > 0 ) {
            m_edge_info = empty<PrimaryEdgeInfo>(offset.back());
            for ( size_t i = 0; i < meshes.size(); ++i ) {
                const Mesh *mesh = meshes[i];
                if ( mesh->m_enable_edges ) {
                    const auto &info = per_mesh_info[i];
                    const int m = static_cast<int>(slices(info));
                    const IntD idx = arange<IntD>(m) + offset[i];

                    Vector3fD p0 = gather<Vector3fD>(mesh->m_vertex_positions, info[0]),
                              p1 = gather<Vector3fD>(mesh->m_vertex_positions, info[1]);

                    Vector3fD q0 = transform_pos(m_world_to_sample, p0),
                              q1 = transform_pos(m_world_to_sample, p1);

#ifdef PSDR_PRIMARY_EDGE_VIS_CHECK
                    scatter(m_edge_info.p0, q0, idx);
                    scatter(m_edge_info.p1, q1, idx);
#else
                    scatter(m_edge_info.p0, head<2>(q0), idx);
                    scatter(m_edge_info.p1, head<2>(q1), idx);
#endif

                    Vector2fD e = head<2>(detach(q1) - detach(q0));
                    FloatD len = norm(e);
                    e /= len;
                    scatter(m_edge_info.edge_normal, Vector2fD(-e.y(), e.x()), idx);
                    scatter(m_edge_info.edge_length, len, idx);
                }
            }
            m_edge_distrb.init(detach(m_edge_info.edge_length));
            m_enable_edges = true;
        } else {
            m_enable_edges = false;
        }
    }
}


std::string PerspectiveCamera::to_string() const {
    return "PerspectiveCamera";
}


RayC PerspectiveCamera::sample_primary_ray(const Vector2fC &samples) const {
    Vector3fC d = normalize(transform_pos<FloatC>(detach(m_sample_to_camera), concat(samples, 0.f)));
    Matrix4fD m_to_world = m_to_world_left * m_to_world_raw * m_to_world_right;
    Matrix4fC to_world = detach(m_to_world);
    return RayC(
        transform_pos<FloatC>(to_world, zero<Vector3fC>(slices(samples))),
        transform_dir<FloatC>(to_world, d)
    );
}


RayD PerspectiveCamera::sample_primary_ray(const Vector2fD &samples) const {
    Vector3fD d = normalize(transform_pos<FloatD>(m_sample_to_camera, concat(samples, 0.f)));
    Matrix4fD m_to_world = m_to_world_left * m_to_world_raw * m_to_world_right;
    return RayD(
        transform_pos<FloatD>(m_to_world, zero<Vector3fD>(slices(samples))),
        transform_dir<FloatD>(m_to_world, d)
    );
}


SensorDirectSampleC PerspectiveCamera::sample_direct(const Vector3fC &p) const {
    SensorDirectSampleC result;
    result.q = head<2>(transform_pos<FloatC>(detach(m_world_to_sample), p));

    Vector2iC iq = floor2int<Vector2iC, Vector2fC>(result.q*m_resolution);
    result.is_valid = iq.x() >= 0 && iq.x() < m_resolution.x() &&
                      iq.y() >= 0 && iq.y() < m_resolution.y();
    result.pixel_idx = select(result.is_valid, iq.y()*m_resolution.x() + iq.x(), -1);

    Vector3fC dir = p - detach(m_camera_pos);
    FloatC dist2 = squared_norm(dir);
    dir /= safe_sqrt(dist2);

    FloatC cosTheta = dot(detach(m_camera_dir), dir);
    result.sensor_val = rcp(dist2)*pow(rcp(cosTheta), 3.f)*detach(m_inv_area);
    return result;
}


PrimaryEdgeSample PerspectiveCamera::sample_primary_edge(const FloatC &_sample1) const {
    FloatC sample1 = _sample1;

    PrimaryEdgeSample result;
    const int m = static_cast<int>(slices(sample1));

    IntC edge_idx;
    std::tie(edge_idx, result.pdf) = m_edge_distrb.sample_reuse<false>(sample1);

    PrimaryEdgeInfo edge_info = gather<PrimaryEdgeInfo>(m_edge_info, IntD(edge_idx));
    result.pdf /= detach(edge_info.edge_length);

    Vector2fC edge_normal = detach(edge_info.edge_normal);
#ifdef PSDR_PRIMARY_EDGE_VIS_CHECK
    const Vector3fD &p0 = edge_info.p0, &p1 = edge_info.p1;
    Vector3fD p_3   = fmadd(p0, 1.0f - sample1, p1*sample1);
    Vector2fD p_    = Vector2fD(p_3.x(), p_3.y());
#else
    const Vector2fD &p0 = edge_info.p0, &p1 = edge_info.p1;
    Vector2fD p_    = fmadd(p0, 1.0f - sample1, p1*sample1);
#endif
    Vector2fC p     = detach(p_);
    result.x_dot_n  = dot(p_, edge_normal);

    Vector2iC ip    = floor2int<Vector2iC, Vector2fC>(p*m_resolution);
    MaskC valid     = ip.x() >= 0 && ip.x() < m_resolution.x() &&
                      ip.y() >= 0 && ip.y() < m_resolution.y();

    result.idx      = full<IntC>(-1, m);
    masked(result.idx, valid) = ip.y()*m_resolution.x() + ip.x();

    result.ray_p    = sample_primary_ray(p + EdgeEpsilon*edge_normal);
    result.ray_n    = sample_primary_ray(p - EdgeEpsilon*edge_normal);

#ifdef PSDR_PRIMARY_EDGE_VIS_CHECK
    RayC &ray_c     = result.ray_c;
    ray_c           = sample_primary_ray(p);
    Vector3fD q     = transform_pos(m_sample_to_world, p_3);
    ray_c.tmax      = norm(detach(q) - detach(m_camera_pos)) - 100.f*ShadowEpsilon;    // Being conservative here to avoid numerical issues
#endif

    return result;
}

} // namespace psdr
