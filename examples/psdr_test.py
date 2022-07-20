import argparse
from config import scene_path, output_path, psdr_tests
from run_test import *

def create_output_dir(dirname):
    if not os.path.isdir(dirname):
        try:
            os.mkdir(dirname)
        except OSError:
            print ("Creation of the directory '%s' failed" % dirname)
            exit()

def field_test(name, args):
    dirname = output_path + args["fname"]
    create_output_dir(dirname)

    field_name = args["field_name"]
    print("# %s (field = %s):" % (name, field_name))
    if "orig" in args and args["orig"]:
        sc = psdr_cuda.Scene()
        sc.load_file(scene_path + args["scene_file"], False)
        integrator = psdr_cuda.FieldExtractionIntegrator(field_name)
        sc.opts.sppe, sc.opts.sppse = 0, 0
        sc.opts.log_level = 0
        sc.configure()
        run_orig(integrator, sc, dirname + "/field_orig.exr", args)
        del sc, integrator

    if "AD" in args:
        sc = psdr_cuda.Scene()
        sc.load_file(scene_path + args["scene_file"], False)
        integrator = psdr_cuda.FieldExtractionIntegrator(field_name)
        sc.opts.log_level = 0
        run_ad(integrator, sc, dirname + "/field_AD.exr", args)
        del sc, integrator

    if "FD" in args:
        integrator = psdr_cuda.FieldExtractionIntegrator(field_name)
        run_fd(integrator, dirname + "/field_FD.exr", args)
        del integrator

    print()

def direct_test(name, args):
    dirname = output_path + args["fname"]
    create_output_dir(dirname)

    print("# %s (direct):" % name)
    if "orig" in args and args["orig"]:
        sc = psdr_cuda.Scene()
        sc.load_file(scene_path + args["scene_file"], False)
        sc.opts.log_level = 0
        sc.configure()
        integrator = psdr_cuda.DirectIntegrator(bsdf_samples=args["bsdf_samples"], light_samples=args["light_samples"])
        # integrator.hide_emitters = True
        run_orig(integrator, sc, dirname + "/direct_orig.exr", args)
        del sc, integrator

    if "AD" in args:
        sc = psdr_cuda.Scene()
        sc.load_file(scene_path + args["scene_file"], False)
        sc.opts.log_level = 0
        if "edge_sort" in args["AD"]:
            sc.edge_sort_opt.enable_sort = True
            sc.edge_sort_opt.local_angle = args["AD"]["edge_sort"][0]
            sc.edge_sort_opt.global_angle = args["AD"]["edge_sort"][1]
            sc.edge_sort_opt.min_global_step = args["AD"]["edge_sort"][2]
            sc.edge_sort_opt.max_depth = args["AD"]["edge_sort"][3]

        else:
            sc.edge_sort_opt.enable_sort = False

        integrator = psdr_cuda.DirectIntegrator(bsdf_samples=args["bsdf_samples"], light_samples=args["light_samples"], edge_direct=args["AD"]["edge_direct"])
        run_ad(integrator, sc, dirname + "/direct_AD.exr", args)
        del sc, integrator

    if "FD" in args:
        print("FD")
        integrator = psdr_cuda.DirectIntegrator(bsdf_samples=args["bsdf_samples"], light_samples=args["light_samples"])
        run_fd(integrator, dirname + "/direct_FD.exr", args)
        del integrator
    print()

def directDual_test(name, args):
    dirname = output_path + args["fname"]
    create_output_dir(dirname)

    print("# %s (directDual):" % name)
    if "orig" in args and args["orig"]:
        sc = psdr_cuda.Scene()
        sc.load_file(scene_path + args["scene_file"], False)
        sc.opts.log_level = 0

        sc.configure()
        integrator = psdr_cuda.DirectIntegratorDual(bsdf_samples=args["bsdf_samples"], light_samples=args["light_samples"])
        run_orig(integrator, sc, dirname + "/direct_orig.exr", args)
        del sc, integrator

    if "AD" in args:
        sc = psdr_cuda.Scene()
        sc.load_file(scene_path + args["scene_file"], False)
        sc.opts.log_level = 0
        integrator = psdr_cuda.DirectIntegratorDual(bsdf_samples=args["bsdf_samples"], light_samples=args["light_samples"], edge_direct=args["AD"]["edge_direct"])
        run_ad(integrator, sc, dirname + "/direct_AD.exr", args)
        del sc, integrator

    if "FD" in args:
        print("FD")
        integrator = psdr_cuda.DirectIntegratorDual(bsdf_samples=args["bsdf_samples"], light_samples=args["light_samples"])
        run_fd(integrator, dirname + "/direct_FD.exr", args)
        del integrator
    print()

def path_test(name, args):
    dirname = output_path + args["fname"]
    create_output_dir(dirname)

    print("# %s (PathTracer):" % name)
    if "orig" in args and args["orig"]:
        sc = psdr_cuda.Scene()
        sc.load_file(scene_path + args["scene_file"], False)
        sc.opts.log_level = 0
        sc.configure()
        integrator = psdr_cuda.PathTracer(max_depth=args["max_depth"])
        run_orig(integrator, sc, dirname + "/path_orig" + str(args["max_depth"]) + ".exr", args)
        del sc, integrator
    if "AD" in args:
        sc = psdr_cuda.Scene()
        sc.load_file(scene_path + args["scene_file"], False)
        sc.opts.log_level = 0
        integrator = psdr_cuda.PathTracer(max_depth=args["max_depth"])
        run_ad(integrator, sc, dirname + "/path_AD" + str(args["max_depth"]) + ".exr", args)
        del sc, integrator

    if "FD" in args:
        integrator = psdr_cuda.PathTracer(max_depth=args["max_depth"])
        run_fd(integrator, dirname + "/path_FD" + str(args["max_depth"]) + ".exr", args)
        del integrator
    print()

def process(name, args):
    dirname = output_path + args["fname"]
    
    if args["test_type"] == "field":
        field_test(name, args)
    elif args["test_type"] == "directDual":
        directDual_test(name, args)
    elif args["test_type"] == "direct":
        direct_test(name, args)

    elif args["test_type"] == "path":
        path_test(name, args)
    else:
        raise Exception("Incorrect test type")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='psdr_test',
        description='psdr_cuda example for SIG22 PSDR Guiding Paper',
        epilog='Kai Yan (kyan8@uci.edu)'
    )
    parser.add_argument('--test', required=False, nargs=1)
    args = parser.parse_args()

    if not os.path.isdir(output_path):
        try:
            os.mkdir(output_path)
        except OSError:
            print ("Creation of the directory %s failed" % output_path)
            exit()

    if args.test is None:
        for name, value in psdr_tests.items():
            process(name, value)
    else:
        name = args.test[0]
        process(name, psdr_tests[name])
