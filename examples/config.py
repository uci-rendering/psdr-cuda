scene_path = "./data/scenes/"
output_path = "./outputs/"

AD_config1 = {
    "type"          : "mesh_transform",
    "Mesh_ID"       : [1],
    "Mesh_dir"      : [[1.0, 0.0, 0.0]],
    "spp"           : 8,
    "sppe"          : 8,
    "sppse"         : 8,
    "guide"         : { "reso" : [40000, 5, 5, 2], "nround" : 16 }

}

AD_config2 = {
    "type"          : "mesh_transform",
    "Mesh_ID"       : [1, 2],
    "Mesh_dir"      : [[100.0, 0.0, 0.0], [0.0, 100.0, 0.0]]
}

AD_config3 = {
    "type"          : "vertex_transform",
    "Mesh_ID"       : [0],
    "Vertex_ID"     : [0],
    "Vertex_dir"    : [[-50.0, 0.0, 0.0]],
    "spp"           : 16,
    "sppe"          : 8,
    "sppse"         : 64,
    "guide"         : { "reso" : [40000, 5, 5, 2], "nround" : 32 }
}

AD_config4 = {
    "type"          : "mesh_rotate",
    "Mesh_ID"       : [1],
    "axis"          : [[0., 0., 1.]],
    "spp"           : 0,
    "sppe"          : 0,
    "sppse"         : 64,
    "guide"         : { "reso" : [40000, 5, 5, 2], "nround" : 16 }
}

FD_config1 = { "npass": 512, "eps": 0.1 }
FD_config2 = { "npass": 64, "eps": 0.01 }

psdr_tests = {
    "cbox_MIS" : {
        "test_type"     : "direct",
        "scene_file"    : "cbox_bunny.xml",
        "npass"         : 20,
        "bsdf_samples"  : 2,
        "light_samples" : 2,
        "orig"          : True,
        "AD"            : AD_config3,
        "FD"            : FD_config1,
        "fname"         : "cbox_MIS_sampling"
    },

    "cbox_bs" : {
        "test_type"     : "direct",
        "scene_file"    : "cbox_bunny.xml",
        "npass"         : 100,
        "bsdf_samples"  : 5,
        "light_samples" : 0,
        "orig"          : True,
        "AD"            : AD_config3,
        "fname"         : "cbox_bsdf_sampling"
    },

    "cbox_es" : {
        "test_type"     : "direct",
        "scene_file"    : "cbox_bunny.xml",
        "npass"         : 20,
        "bsdf_samples"  : 0,
        "light_samples" : 2,
        "orig"          : True,
        "AD"            : AD_config3,
        "fname"         : "cbox_emitter_sampling"
    },

    "cbox_mutie" : {
        "test_type"     : "direct",
        "scene_file"    : "cbox_bunny_mutiemitter.xml",
        "npass"         : 2,
        "bsdf_samples"  : 2,
        "light_samples" : 2,
        "orig"          : True,
        "fname"         : "cbox_muti_emitter"
    },

    "tree" : {
        "test_type"     : "direct",
        "scene_file"    : "tree.xml",
        "bsdf_samples"  : 0,
        "light_samples" : 2,
        "orig"          : True,
        "AD"            : {
            "type"          : "mesh_rotate",
            "Mesh_ID"       : [1],
            "axis"          : [[0., 0., 1.]],
            "spp"           : 0,
            "sppe"          : 0,
            "sppse"         : 64,
            "guide"         : { "reso" : [40000, 5, 5, 2], "nround" : 16 },
            "npass"         : 32,
            "no_edge"       : [0, 2]
        },
        "FD"            : { "npass": 64, "eps": 0.01 },
        "fname"         : "tree"
    },

    "bunny_silhouette" : {
        "test_type"     : "field",
        "field_name"    : "silhouette",
        "scene_file"    : "bunny.xml",
        "orig"          : False,
        "AD"            :  {
            "type"          : "mesh_rotate",
            "Mesh_ID"       : [0, 1],
            "axis"          : [[0., 0.1, 0.], [0., -0.1, 0.]],
            "spp"           : 64,
            "sppe"          : 64,
            "sppse"         : 0
        },
        "FD"            : { "npass": 20, "eps": 0.01 },
        "fname"         : "bunny_silhouette"
    },

    "bunny_env_1" : {
        "test_type"     : "direct",
        "scene_file"    : "bunny_env.xml",
        "bsdf_samples"  : 4,
        "light_samples" : 4,
        "orig"          : True,
        "AD"            :  {
            "type"          : "envmap_rotate",
            "Emitter_ID"    : 0,
            "axis"          : [0., 0.1, 0.],
            "spp"           : 64,
            "sppe"          : 0,
            "sppse"         : 0,
            "npass"         : 25
        },
        "FD"            : { "npass": 25, "eps": 0.01 },
        "fname"         : "bunny_env_1"
    },

    "bunny_env_2" : {
        "test_type"     : "direct",
        "scene_file"    : "bunny_env_2.xml",
        "bsdf_samples"  : 2,
        "light_samples" : 2,
        "orig"          : True,
        "npass"         : 8,
        "AD"            : {
            "type"          : "mesh_rotate",
            "Mesh_ID"       : [0],
            "axis"          : [[0., 0., 1.]],
            "spp"           : 4,
            "sppe"          : 4,
            "sppse"         : 64,
            "npass"         : 40,
            "no_edge"       : [1]
        },
        "FD"            : { "npass": 64, "eps": 0.01 },
        "fname"         : "bunny_env_2"
    }
}
