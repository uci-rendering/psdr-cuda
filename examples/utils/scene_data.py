import psdr_cuda
from utils.system import *
from utils.direct import run_direct, run_direct_ad, run_direct_FD
from utils.collocated import run_collocated

def create_output_dir(dirname):
    if not os.path.isdir(dirname):
        try:
            os.mkdir(dirname)
        except OSError:
            print ("Creation of the directory '%s' failed" % dirname)
            exit()

def collocated_test(name, args):
    dirname = output_path + args["fname"]
    create_output_dir(dirname)

    print("# %s (collocated):" % name)
    sc = psdr_cuda.Scene()
    sc.load_file(scene_path + args["scene_file"])
    run_collocated(sc, dirname + "/collocated.exr", args)
    del sc
    print()

def direct_test(name, args):
    dirname = output_path + args["fname"]
    create_output_dir(dirname)

    sc = psdr_cuda.Scene()
    sc.load_file(scene_path+args["scene_file"], False)
    sc.opts.log_level = 0
    print("# %s (direct):" % name)

    if "orig" in args and args["orig"]:
        run_direct(sc, dirname + "/direct_org.exr", args)

    if "AD" in args:
        run_direct_ad(sc, dirname + "/direct.exr", args)

    if "FD" in args:
        run_direct_FD(dirname + "/direct_FD.exr", args)
    print()
    del sc
