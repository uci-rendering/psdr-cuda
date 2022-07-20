scene_path = "./"
output_path = "./outputs/"


# edge sort is enable by default
# "edge_sort": [170, 175, 5, 20] controls the greedy search for draw edge segments
# "guide" : { "config" : [1, 1, 1, 0.003, 100000, 100000, 0, 16, 0.0, 0.0], "option" : 0 } is adaptive guiding without MIS
# "guide"         : { "config" : [1, 1, 1, 0.01, 100000, 100000, 0, 48, 0.0, 0.0,
#                                 1, 1, 1, 0.01, 100000, 100000, 0, 48, 0.0, 0.0], "option" : 2 } is adaptive guiding with MIS
# "guide"         : { "config" : [100, 100, 100, 2, 8], "option" : 1 } is the [Zhang 2020] regular guiding
# "edge_direct"   : 0 is direcition sampling for boundary term
# "edge_direct"   : 1 is emitter sampling for boundary term
# "edge_direct"   : 2 is MIS sampling for boundary term


psdr_tests = {
    "example" : {
        "test_type"     : "direct",
        "scene_file"    : "example.xml",
        "bsdf_samples"  : 1,
        "light_samples" : 1,
        "orig"          : False,
        "npass"         : 8,
        "AD"            : {
            "edge_direct"   : 2, # 0 or 1 or 2
            "type"          : "mesh_rotate",
            "Mesh_ID"       : [0],
            "axis"          : [[0., 0., 1.]],
            "spp"           : 0,
            "sppe"          : 0,
            "sppse"         : 8,
            "npass"         : 8,
            # "edge_sort"     : [170, 175, 5, 20], # enable or not

            # "guide"         : { "config" : [1, 1, 1, 0.003, 100000, 100000, 0, 16, 0.0, 0.0], "option" : 0 }, # emi or dir guiding
            # "guide"         : { "config" : [1, 1, 1, 0.01, 100000, 100000, 0, 16, 0.0, 0.0], "option" : 0 }, # emi or dir guiding

             # "guide"         : { "config" : [1, 1, 1, 0.01, 100000, 100000, 0, 48, 0.0, 0.0,
             #                                 1, 1, 1, 0.01, 100000, 100000, 0, 48, 0.0, 0.0], "option" : 2 }, # MIS guiding

            # "guide"         : { "config" : [100, 100, 100, 2, 8], "option" : 1 }, # [Zhang 2020] baseline
            "no_edge"       : [1]
        },
        # "FD"            : { "npass": 10240, "eps": 0.01 },
        "fname"         : "example"
    }
}
