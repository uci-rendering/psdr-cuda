#pragma once

#include <psdr/psdr.h>

namespace pugi
{
    class xml_document;
    class xml_node;
} // namespace pugi


namespace psdr
{

class SceneLoader {
public:
    static void load_from_file(const char *file_name, Scene &scene);
    static void load_from_string(const char *scene_xml, Scene &scene);

protected:
    static void load_scene(const pugi::xml_document &doc, Scene &scene);
    static void load_sensor(const pugi::xml_node &node, Scene &scene);
    static void load_emitter(const pugi::xml_node &node, Scene &scene);
    static void load_bsdf(const pugi::xml_node &node, Scene &scene);
    static void load_shape(const pugi::xml_node &node, Scene &scene, int shape_id = -1);
};

} // namespace psdr
