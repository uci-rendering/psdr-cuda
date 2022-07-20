import psdr_cuda
import enoki as ek
from enoki.cuda_autodiff import Float32 as FloatD, Matrix4f as Matrix4fD, Vector3f as Vector3fD
import numpy as np
from matplotlib import pyplot as plt
import cv2
import time, os
from utils.differential import *
from config import scene_path, output_path

time_threshold = 0.2

def run_orig(integrator, sc, fname, args):
    global time_threshold

    ro = sc.opts

    npass = args["npass"] if "npass" in args else 1
    num_sensors = sc.num_sensors
    img_org = [None]*num_sensors

    t0 = time.process_time()
    t1 = t0
    for i in range(npass):
        for sensor_id in range(num_sensors):
            img = integrator.renderC(sc, sensor_id)
            if i == 0:
                img_org[sensor_id] = img.numpy()
            else:
                img_org[sensor_id] += img.numpy()
            del img

        t2 = time.process_time()
        if t2 - t1 > time_threshold:
            print("(%d/%d) done in %.2f seconds." % (i + 1, npass, t2 - t0), end="\r")
            t1 = t2
    print("(%d/%d) Total orig. rendering time: %.2f seconds." % (npass, npass, t2 - t0))

    for sensor_id in range(num_sensors):
        img = (img_org[sensor_id]/float(npass)).reshape((ro.height, ro.width, 3))
        output = cv2.cvtColor(img, cv2.COLOR_RGB2BGR)
        cv2.imwrite(fname[:-4] + "_" + str(sensor_id) + fname[-4:], output)


def run_ad(integrator, sc, fname, args):
    global time_threshold
    ad_config = args["AD"]

    if "spp" in ad_config:
        sc.opts.spp = ad_config["spp"]
    if "sppe" in ad_config:
        sc.opts.sppe = ad_config["sppe"]
    if "sppse" in ad_config:
        sc.opts.sppse = ad_config["sppse"]

    if "no_edge" in ad_config:
        for i in ad_config["no_edge"]:
            sc.param_map["Mesh[" + str(i) + "]"].enable_edges = False

    ro = sc.opts
    if ad_config["type"] == "mesh_transform":
        if len(ad_config["Mesh_ID"]) != len(ad_config["Mesh_dir"]):
            raise Exception("Mesh_ID and Mesh_dir have different sizes")
    elif ad_config["type"] == "mesh_rotate":
        if len(ad_config["Mesh_ID"]) != len(ad_config["axis"]):
            raise Exception("Mesh_ID and axis have different sizes")
    elif ad_config["type"] == "vertex_transform":
        if len(ad_config["Mesh_ID"]) != len(ad_config["Vertex_ID"]):
            raise Exception("Mesh_ID and Vertex_ID have different sizes")
        orig_vtx_pos = {}
        for j in ad_config["Mesh_ID"]:
            mesh_obj = sc.param_map["Mesh[" + str(j) + "]"]
            orig_vtx_pos[j] = ek.detach(mesh_obj.vertex_positions)
    elif ad_config["type"] == "rc_roughness":
        base_roughness = {}
        for j in ad_config["BSDF_ID"]:
            bsdf_obj = sc.param_map["BSDF[" + str(j) + "]"]
            base_roughness[j] = (ek.detach(bsdf_obj.alpha_u.data), ek.detach(bsdf_obj.alpha_v.data))
    elif ad_config["type"] == "material_roughness":
        base_roughness = {}
        for j in ad_config["BSDF_ID"]:
            bsdf_obj = sc.param_map["BSDF[" + str(j) + "]"]
            base_roughness[j] = ek.detach(bsdf_obj.roughness.data)
    elif ad_config["type"] == "material_diffuse":
        base_roughness = {}
        for j in ad_config["BSDF_ID"]:
            bsdf_obj = sc.param_map["BSDF[" + str(j) + "]"]
            base_roughness[j] = ek.detach(bsdf_obj.diffuseReflectance.data)

    elif ad_config["type"] == "envmap_rotate":
        if "Emitter_ID" not in ad_config:
            raise Exception("Missing Emitter_ID")
    else:
        raise Exception("Unknown transform")

    if "npass" in ad_config:
        npass = ad_config["npass"]
    elif "npass" in args:
        npass = args["npass"]
    else:
        npass = 1

    num_sensors = sc.num_sensors
    print("total sensor:", num_sensors)
    img_ad = [None]*num_sensors

    total_time = 0.0
    t0 = time.process_time()
    t1 = t0
    for i in range(npass):
        for sensor_id in range(num_sensors):
            # AD config
            P = FloatD(0.)
            ek.set_requires_gradient(P)

            if ad_config["type"] == "mesh_transform":
                for j in range(len(ad_config["Mesh_ID"])):
                    mesh_transform(sc, ad_config["Mesh_ID"][j], Vector3fD(ad_config["Mesh_dir"][j]) * P)
            elif ad_config["type"] == "mesh_rotate":
                for j in range(len(ad_config["Mesh_ID"])):
                    mesh_rotate(sc, ad_config["Mesh_ID"][j], Vector3fD(ad_config["axis"][j]), P)
            elif ad_config["type"] == "vertex_transform":
                for j in range(len(ad_config["Mesh_ID"])):
                    vertex_transform(sc, ad_config["Mesh_ID"][j], ad_config["Vertex_ID"][j], ad_config["Vertex_dir"][j], orig_vtx_pos[j], P)
            elif ad_config["type"] == "rc_roughness":
                for j in ad_config["BSDF_ID"]:
                    rc_roughness(sc, j, base_roughness[j], P)
            elif ad_config["type"] == "material_roughness":
                for j in ad_config["BSDF_ID"]:
                    material_roughness(sc, j, base_roughness[j], P)
            elif ad_config["type"] == "material_diffuse":
                for j in ad_config["BSDF_ID"]:
                    material_diffuse(sc, j, base_roughness[j], P)
            elif ad_config["type"] == "envmap_rotate":
                envmap_rotate(sc, ad_config["Emitter_ID"], ad_config["axis"], P)
            # End AD config
            # if i == 0:
            sc.configure()

            before = time.process_time()
            if i == 0 and "guide" in ad_config:
                t2 = time.process_time()
                guide_info = ad_config["guide"]

                if guide_info["option"] == 2:
                    if sensor_id == 0:
                        active_sensor = np.arange(0,num_sensors)
                        integrator.preprocess_secondary_edges(sc, active_sensor, np.array(guide_info["config"]), guide_info["option"])
                        print("Global guiding done in %.2f seconds." % (time.process_time() - t2))
                else:
                    integrator.preprocess_secondary_edges(sc, [sensor_id], np.array(guide_info["config"]), guide_info["option"])
                    print("Single Camera guiding done in %.2f seconds." % (time.process_time() - t2))

            img = integrator.renderD(sc, sensor_id)
            after = time.process_time()
            total_time += after-before

            ek.forward(P, free_graph=True)
            grad_img = ek.gradient(img).numpy()
            grad_img[np.logical_not(np.isfinite(grad_img))] = 0.
            if i == 0:
                img_ad[sensor_id] = grad_img
            else:
                img_ad[sensor_id] += grad_img
            del img
            del P

        t2 = time.process_time()
        if t2 - t1 > time_threshold:
            # print("(%d/%d) done in %.2f seconds." % (i + 1, npass, t2 - t0), end="\r")
            print("(%d/%d) done in %.2f seconds." % (i + 1, npass, total_time), end="\r")
            t1 = t2
    # print("(%d/%d) Total AD rendering time: %.2f seconds." % (npass, npass, t2 - t0))

    print("(%d/%d) Total AD rendering time: %.2f seconds." % (npass, npass, total_time))

    for sensor_id in range(num_sensors):
        img = (img_ad[sensor_id]/float(npass)).reshape((ro.height, ro.width, 3))
        output = cv2.cvtColor(img, cv2.COLOR_RGB2BGR)
        cv2.imwrite(fname[:-4] + "_" + str(sensor_id) + fname[-4:], output)


def run_fd(integrator, fname, args):
    global time_threshold

    sc1, sc2 = psdr_cuda.Scene(), psdr_cuda.Scene()
    sc1.load_file(scene_path + args["scene_file"], False)
    sc2.load_file(scene_path + args["scene_file"], False)
    sc1.opts.sppe, sc1.opts.sppse = 0, 0
    sc2.opts.sppe, sc2.opts.sppse = 0, 0
    sc1.opts.log_level = 0
    sc2.opts.log_level = 0

    ad_config = args["AD"]
    eps = args["FD"]["eps"]
    if ad_config["type"] == "mesh_transform":
        for tid in range(len(ad_config["Mesh_ID"])):
            mesh_transform(sc1, ad_config["Mesh_ID"][tid], Vector3fD(ad_config["Mesh_dir"][tid]) * FloatD(-eps))
            mesh_transform(sc2, ad_config["Mesh_ID"][tid], Vector3fD(ad_config["Mesh_dir"][tid]) * FloatD( eps))
    elif ad_config["type"] == "mesh_rotate":
        for tid in range(len(ad_config["Mesh_ID"])):
            mesh_rotate(sc1, ad_config["Mesh_ID"][tid], Vector3fD(ad_config["axis"][tid]), FloatD(-eps))
            mesh_rotate(sc2, ad_config["Mesh_ID"][tid], Vector3fD(ad_config["axis"][tid]), FloatD( eps))
    elif ad_config["type"] == "vertex_transform":
        assert len(ad_config["Mesh_ID"]) == len(ad_config["Vertex_ID"])
        for j in range(len(ad_config["Mesh_ID"])):
            mesh_obj = sc1.param_map["Mesh[" + str(j) + "]"]
            orig_vtx_pos = ek.detach(mesh_obj.vertex_positions)
            vertex_transform(sc1, ad_config["Mesh_ID"][j], ad_config["Vertex_ID"][j], ad_config["Vertex_dir"][j], orig_vtx_pos, FloatD(-eps))
            vertex_transform(sc2, ad_config["Mesh_ID"][j], ad_config["Vertex_ID"][j], ad_config["Vertex_dir"][j], orig_vtx_pos, FloatD( eps))
    
    elif ad_config["type"] == "rc_roughness":
        for j in ad_config["BSDF_ID"]:
            bsdf_obj = sc1.param_map["BSDF[" + str(j) + "]"]
            base_roughness = (ek.detach(bsdf_obj.alpha_u.data), ek.detach(bsdf_obj.alpha_v.data))
            rc_roughness(sc1, j, base_roughness, -eps)
            rc_roughness(sc2, j, base_roughness,  eps)


    elif ad_config["type"] == "material_roughness":
        for j in ad_config["BSDF_ID"]:
            bsdf_obj = sc1.param_map["BSDF[" + str(j) + "]"]
            base_roughness = ek.detach(bsdf_obj.roughness.data)
            material_roughness(sc1, j, base_roughness, -eps)
            material_roughness(sc2, j, base_roughness,  eps)
    elif ad_config["type"] == "material_diffuse":
        for j in ad_config["BSDF_ID"]:
            bsdf_obj = sc1.param_map["BSDF[" + str(j) + "]"]
            base_roughness = ek.detach(bsdf_obj.diffuseReflectance.data)
            material_diffuse(sc1, j, base_roughness, -eps)
            material_diffuse(sc2, j, base_roughness,  eps)
    elif ad_config["type"] == "envmap_rotate":
        envmap_rotate(sc1, ad_config["Emitter_ID"], ad_config["axis"], -eps)
        envmap_rotate(sc2, ad_config["Emitter_ID"], ad_config["axis"],  eps)
    else:
        raise Exception("Unknown transform")

    sc1.configure()
    sc2.configure()

    ro = sc1.opts

    if "npass" in args["FD"]:
        npass = args["FD"]["npass"]
    elif "npass" in args:
        npass = args["npass"]
    else:
        npass = 1

    num_sensors = sc1.num_sensors
    img1_org = [None]*num_sensors
    img2_org = [None]*num_sensors

    t0 = time.process_time()
    t1 = t0
    for i in range(npass):
        for sensor_id in range(num_sensors):
            img1 = integrator.renderC(sc1)
            img2 = integrator.renderC(sc2)
            if i == 0:
                img1_org[sensor_id] = img1.numpy()
                img2_org[sensor_id] = img2.numpy()
            else:
                img1_org[sensor_id] += img1.numpy()
                img2_org[sensor_id] += img2.numpy()
            del img1, img2

        t2 = time.process_time()
        if t2 - t1 > time_threshold:
            print("(%d/%d) done in %.2f seconds." % (i + 1, npass, t2 - t0), end="\r")
            t1 = t2
    print("(%d/%d) Total FD rendering time: %.2f seconds." % (npass, npass, t2 - t0))
    del sc1, sc2

    for sensor_id in range(num_sensors):
        img_FD = (img2_org[sensor_id] - img1_org[sensor_id])/(2.0*eps*float(npass))
        img = img_FD.reshape((ro.height, ro.width, 3))
        output = cv2.cvtColor(img, cv2.COLOR_RGB2BGR)
        cv2.imwrite(fname[:-4] + "_" + str(sensor_id) + fname[-4:], output)
