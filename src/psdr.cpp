#include <psdr/psdr.h>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
// #include <pybind11/stl_bind.h>
#include <enoki/python.h>

#include <psdr/core/ray.h>
#include <psdr/core/frame.h>
#include <psdr/core/bitmap.h>
#include <psdr/core/intersection.h>
#include <psdr/core/sampler.h>
#include <psdr/core/pmf.h>
#include <psdr/core/cube_distrb.h>
#include <psdr/core/records.h>

#include <psdr/bsdf/diffuse.h>
#include <psdr/bsdf/ggx.h>
#include <psdr/bsdf/roughconductor.h>

#include <psdr/emitter/area.h>
#include <psdr/emitter/envmap.h>

#include <psdr/sensor/perspective.h>

#include <psdr/shape/mesh.h>

#include <psdr/scene/scene.h>
#include <psdr/scene/scene_loader.h>

#include <psdr/integrator/integrator.h>
#include <psdr/integrator/field.h>
#include <psdr/integrator/direct.h>

namespace py = pybind11;
using namespace py::literals;
using namespace psdr;


PYBIND11_MODULE(psdr_cuda, m) {
    py::module::import("enoki");
    py::module::import("enoki.cuda");
    py::module::import("enoki.cuda_autodiff");

    m.doc() = "Path-space differentiable renderer";

    py::class_<Object>(m, "Object")
        .def("type_name", &Object::type_name)
        .def_readonly("id", &Object::m_id)
        .def("__repr__", &Object::to_string);

    py::class_<RenderOption>(m, "RenderOption")
        .def(py::init<>())
        .def(py::init<int, int, int>(), "width"_a, "height"_a, "spp/sppe"_a)
        .def(py::init<int, int, int, int>(), "width"_a, "height"_a, "spp"_a, "sppe"_a)
        .def(py::init<int, int, int, int, int>(), "width"_a, "height"_a, "spp"_a, "sppe"_a, "sppse"_a)
        .def_readwrite("width", &RenderOption::width)
        .def_readwrite("height", &RenderOption::height)
        .def_readwrite("spp", &RenderOption::spp)
        .def_readwrite("sppe", &RenderOption::sppe)
        .def_readwrite("sppse", &RenderOption::sppse)
        .def_readwrite("log_level", &RenderOption::log_level)
        .def("__repr__",
            [](const RenderOption &ro) {
                std::stringstream oss;
                oss << "[width: " << ro.width << ", height: " << ro.height
                    << ", spp: " << ro.spp << ", sppe: " << ro.sppe << ", sppse: " << ro.sppse
                    << ", log_level: " << ro.log_level << "]";
                return oss.str();
            }
        );

    // Core classes

    py::class_<RayC>(m, "RayC")
        .def(py::init<>())
        .def("reversed", &RayC::reversed)
        .def_readwrite("o", &RayC::o)
        .def_readwrite("d", &RayC::d);

    py::class_<RayD>(m, "RayD")
        .def(py::init<>())
        .def("reversed", &RayD::reversed)
        .def_readwrite("o", &RayD::o)
        .def_readwrite("d", &RayD::d);

    py::class_<FrameC>(m, "FrameC")
        .def(py::init<>())
        .def(py::init<const Vector3fC &>())
        .def_readwrite("s", &FrameC::s)
        .def_readwrite("t", &FrameC::t)
        .def_readwrite("n", &FrameC::n);

    py::class_<FrameD>(m, "FrameD")
        .def(py::init<>())
        .def(py::init<const Vector3fD &>())
        .def_readwrite("s", &FrameD::s)
        .def_readwrite("t", &FrameD::t)
        .def_readwrite("n", &FrameD::n);

    py::class_<Bitmap1fD>(m, "Bitmap1fD")
        .def(py::init<>())
        .def(py::init<float>())
        .def(py::init<int, int, const FloatD&>())
        .def(py::init<const char*>())
        .def("load_openexr", &Bitmap1fD::load_openexr)
        .def("eval", &Bitmap1fD::eval<true>, "uv"_a, "flip_v"_a = true)
        .def_readwrite("resolution", &Bitmap1fD::m_resolution)
        .def_readwrite("data", &Bitmap1fD::m_data);

    py::class_<Bitmap3fD>(m, "Bitmap3fD")
        .def(py::init<>())
        .def(py::init<const ScalarVector3f&>())
        .def(py::init<int, int, const SpectrumD&>())
        .def(py::init<const char*>())
        .def("load_openexr", &Bitmap3fD::load_openexr)
        .def("eval", &Bitmap3fD::eval<true>, "uv"_a, "flip_v"_a = true)
        .def_readwrite("resolution", &Bitmap3fD::m_resolution)
        .def_readwrite("data", &Bitmap3fD::m_data);

    /*
    py::class_<InteractionC>(m, "InteractionC")
        .def("is_valid", &InteractionC::is_valid)
        .def_readonly("wi", &InteractionC::wi)
        .def_readonly("p", &InteractionC::p)
        .def_readonly("t", &InteractionC::t);

    py::class_<InteractionD>(m, "InteractionD")
        .def("is_valid", &InteractionD::is_valid)
        .def_readonly("wi", &InteractionD::wi)
        .def_readonly("p", &InteractionD::p)
        .def_readonly("t", &InteractionD::t);

    py::class_<IntersectionC, InteractionC>(m, "IntersectionC")
        .def_readonly("shape", &IntersectionC::shape)
        .def_readonly("n", &IntersectionC::n)
        .def_readonly("sh_frame", &IntersectionC::sh_frame)
        .def_readonly("uv", &IntersectionC::uv)
        .def_readonly("J", &IntersectionC::J);

    py::class_<IntersectionD, InteractionD>(m, "IntersectionD")
        .def_readonly("shape", &IntersectionD::shape)
        .def_readonly("n", &IntersectionD::n)
        .def_readonly("sh_frame", &IntersectionD::sh_frame)
        .def_readonly("uv", &IntersectionD::uv)
        .def_readonly("J", &IntersectionD::J);

    py::class_<Sampler>(m, "Sampler")
        .def(py::init<>())
        .def("clone", &Sampler::clone)
        .def("seed", &Sampler::seed)
        .def("next_1d", &Sampler::next_1d<false>)
        .def("next_2d", &Sampler::next_2d<false>);
    */

    py::class_<DiscreteDistribution>(m, "DiscreteDistribution")
        .def(py::init<>())
        .def("init", &DiscreteDistribution::init)
        .def("sample", &DiscreteDistribution::sample)
        .def_readonly("sum", &DiscreteDistribution::m_sum)
        .def("pmf", &DiscreteDistribution::pmf);

    py::class_<HyperCubeDistribution2f>(m, "HyperCubeDistribution2f")
        .def(py::init<>())
        .def("set_resolution", &HyperCubeDistribution2f::set_resolution)
        .def("set_mass", &HyperCubeDistribution2f::set_mass)
        .def("sample_reuse", &HyperCubeDistribution2f::sample_reuse)
        .def("pdf", &HyperCubeDistribution2f::pdf)
        .def_readonly("cells", &HyperCubeDistribution2f::m_cells);

    py::class_<HyperCubeDistribution3f>(m, "HyperCubeDistribution3f")
        .def(py::init<>())
        .def("set_resolution", &HyperCubeDistribution3f::set_resolution)
        .def("set_mass", &HyperCubeDistribution3f::set_mass)
        .def("sample_reuse", &HyperCubeDistribution3f::sample_reuse)
        .def("pdf", &HyperCubeDistribution3f::pdf)
        .def_readonly("cells", &HyperCubeDistribution3f::m_cells);

    // Records

    py::class_<SampleRecordC>(m, "SampleRecordC")
        .def_readonly("pdf", &SampleRecordC::pdf)
        .def_readonly("is_valid", &SampleRecordC::is_valid);

    py::class_<SampleRecordD>(m, "SampleRecordD")
        .def_readonly("pdf", &SampleRecordD::pdf)
        .def_readonly("is_valid", &SampleRecordD::is_valid);

    py::class_<PositionSampleC, SampleRecordC>(m, "PositionSampleC")
        .def_readonly("p", &PositionSampleC::p)
        .def_readonly("J", &PositionSampleC::J);

    py::class_<PositionSampleD, SampleRecordD>(m, "PositionSampleD")
        .def_readonly("p", &PositionSampleD::p)
        .def_readonly("J", &PositionSampleD::J);

    // BSDFs

    py::class_<BSDF, Object>(m, "BSDF")
        .def("anisotropic", &BSDF::anisotropic);

    py::class_<Diffuse, BSDF>(m, "DiffuseBSDF")
        //.def(py::init<>())
        //.def(py::init<const ScalarVector3f&>())
        //.def(py::init<const char*>())
        //.def(py::init<const Bitmap3fD&>())
        .def_readwrite("reflectance", &Diffuse::m_reflectance);

    py::class_<RoughConductor, BSDF>(m, "RoughConductorBSDF")
        .def_readwrite("alpha_u", &RoughConductor::m_alpha_u)
        .def_readwrite("alpha_v", &RoughConductor::m_alpha_v)
        .def_readwrite("eta", &RoughConductor::m_eta)
        .def_readwrite("k", &RoughConductor::m_k)
        .def_readwrite("specular_reflectance", &RoughConductor::m_specular_reflectance);

    // Sensors

    py::class_<Sensor, Object>(m, "Sensor")
        .def_readwrite("to_world", &Sensor::m_to_world);

    py::class_<PerspectiveCamera, Sensor>(m, "PerspectiveCamera")
        .def(py::init<float, float, float>())
        .def_readwrite("to_world", &PerspectiveCamera::m_to_world);

    // Emitters

    py::class_<Emitter, Object>(m, "Emitter");

    py::class_<AreaLight, Emitter>(m, "AreaLight")
        .def(py::init<const ScalarVector3f&, const Mesh*>());

    py::class_<EnvironmentMap, Emitter>(m, "EnvironmentMap")
        .def(py::init<const char *>())
        .def_readonly("to_world", &EnvironmentMap::m_to_world_raw)
        .def("set_transform", &EnvironmentMap::set_transform)
        .def_readwrite("radiance", &EnvironmentMap::m_radiance)
        .def_readwrite("scale", &EnvironmentMap::m_scale);

    // Shapes

    py::class_<Mesh, Object>(m, "Mesh")
        .def(py::init<>())
        .def("load", &Mesh::load, "filename"_a, "verbose"_a = false)
        .def("configure", &Mesh::configure)
        .def("set_transform", &Mesh::set_transform, "mat"_a, "set_left"_a = true)
        .def("append_transform", &Mesh::append_transform, "mat"_a, "append_left"_a = true)
        .def("sample_position", py::overload_cast<const Vector2fC&, MaskC>(&Mesh::sample_position, py::const_), "sample2"_a, "active"_a = true)
        .def("sample_position", py::overload_cast<const Vector2fD&, MaskD>(&Mesh::sample_position, py::const_), "sample2"_a, "active"_a = true)
        .def_readonly("num_vertices", &Mesh::m_num_vertices)
        .def_readonly("num_faces", &Mesh::m_num_faces)
        .def_readonly("bsdf", &Mesh::m_bsdf)
        .def_readonly("to_world", &Mesh::m_to_world_raw)
        .def_readwrite("vertex_positions", &Mesh::m_vertex_positions_raw, "Object-space vertex positions")
        .def_readonly("vertex_normals", &Mesh::m_vertex_normals_raw, "Object-space vertex normals")
        .def_readwrite("vertex_uv", &Mesh::m_vertex_uv)
        .def_readwrite("face_indices", &Mesh::m_face_indices)
        .def_readwrite("face_uv_indices", &Mesh::m_face_uv_indices)
#ifdef PSDR_MESH_ENABLE_1D_VERTEX_OFFSET
        .def_readwrite("vertex_offset", &Mesh::m_vertex_offset)
        .def("shift_vertices", &Mesh::shift_vertices)
#endif
        .def_readwrite("enable_edges", &Mesh::m_enable_edges)
        .def("edge_indices", [](const Mesh &mesh) { return head<4>(mesh.m_edge_indices); })
        .def("dump", &Mesh::dump);

    // Scene

    py::class_<Scene, Object>(m, "Scene")
        .def(py::init<>())
        .def("load_file", &Scene::load_file, "file_name"_a, "auto_configure"_a = true)
        .def("load_string", &Scene::load_string, "scene_xml"_a, "auto_configure"_a = true)
        .def("reload_mesh", &Scene::reload_mesh, "mesh"_a, "file_name"_a, "verbose"_a = false)
        .def("configure", &Scene::configure)
        .def("sample_boundary_segment_direct", &Scene::sample_boundary_segment_direct, "sample3"_a, "active"_a = true)
        .def_readwrite("opts", &Scene::m_opts, "Render options")
        .def_readonly("num_sensors", &Scene::m_num_sensors)
        .def_readonly("num_meshes", &Scene::m_num_meshes)
        .def_readonly("param_map", &Scene::m_param_map, "Parameter map");

    // Integrator base and basic integrators

    py::class_<Integrator, Object>(m, "Integrator")
        .def("renderC", &Integrator::renderC, "scene"_a, "sensor_id"_a = 0)
        .def("renderD", &Integrator::renderD, "scene"_a, "sensor_id"_a = 0)
        .def("preprocess_secondary_edges", &Integrator::preprocess_secondary_edges, "scene"_a, "sensor_id"_a, "resolution"_a, "nrounds"_a = 1);

    py::class_<FieldExtractionIntegrator, Integrator>(m, "FieldExtractionIntegrator")
        .def(py::init<const char*>());

    // Direct-illumination Integrators

    py::class_<DirectIntegrator, Integrator>(m, "DirectIntegrator")
        .def(py::init<int, int>(), "bsdf_samples"_a = 1, "light_samples"_a = 1)
        .def_readwrite("hide_emitters", &DirectIntegrator::m_hide_emitters);
}
