import psdr_cuda
import enoki as ek
from enoki.cuda_autodiff import Float32 as FloatD, Vector3f as Vector3fD, Matrix4f as Matrix4fD

def mesh_transform(sc, mesh_ID, dir_vector):
    para = "Mesh[" + str(mesh_ID) + "]"
    sc.param_map[para].set_transform(Matrix4fD.translate(dir_vector))

def mesh_rotate(sc, mesh_ID, axis, angle):
    para = "Mesh[" + str(mesh_ID) + "]"
    sc.param_map[para].set_transform(Matrix4fD.rotate(axis, angle))

def vertex_transform(sc, mesh_ID, vertex_ID, dir_vector, orig_vtx_pos, P):
    assert isinstance(P, FloatD)

    para = "Mesh[" + str(mesh_ID) + "]"
    n = sc.param_map[para].num_vertices
    assert vertex_ID >= 0 and vertex_ID < n

    x_vals = [0.]*n
    y_vals = [0.]*n
    z_vals = [0.]*n
    x_vals[vertex_ID] = dir_vector[0]
    y_vals[vertex_ID] = dir_vector[1]
    z_vals[vertex_ID] = dir_vector[2]
    sc.param_map[para].vertex_positions = Vector3fD(orig_vtx_pos) + Vector3fD(x_vals, y_vals, z_vals)*P

def material_roughness(sc, bsdf_ID, alpha, P):
    para = "BSDF[" + str(bsdf_ID) + "]"
    sc.param_map[para].alpha_u.data = FloatD(alpha[0]) + P
    sc.param_map[para].alpha_v.data = FloatD(alpha[1]) + P

def envmap_rotate(sc, emitter_ID, axis, angle):
    para = "Emitter[" + str(emitter_ID) + "]"
    sc.param_map[para].set_transform(Matrix4fD.rotate(axis, angle))
