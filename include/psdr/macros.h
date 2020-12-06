#pragma once

#ifdef _WIN32
#	define likely(x)       (x)
#	define unlikely(x)     (x)
#else
#	define likely(x)       __builtin_expect((x),1)
#	define unlikely(x)     __builtin_expect((x),0)
#endif

// #define PSDR_OPTIX_DEBUG
// #define PSDR_MESH_ENABLE_1D_VERTEX_OFFSET
// #define PSDR_PRIMARY_EDGE_VIS_CHECK


#define PSDR_CLASS_DECL_BEGIN(_class_, _mode_, _parent_)    \
    class _class_ _mode_ : public _parent_ {


#define PSDR_CLASS_DECL_END(_class_)                        \
    public:                                                 \
        virtual std::string type_name() const override {    \
            return #_class_;                                \
        }                                                   \
    };


#define PSDR_IMPORT_BASE_HELPER(...) Base, ##__VA_ARGS__

#define PSDR_IMPORT_BASE(Name, ...)                         \
    using Base = Name;                                      \
    ENOKI_USING_MEMBERS(PSDR_IMPORT_BASE_HELPER(__VA_ARGS__))
