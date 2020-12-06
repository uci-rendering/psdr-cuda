#pragma once

#include <psdr/psdr.h>

namespace psdr
{


struct BitmapLoader {
    static std::pair<Vector4fC, ScalarVector2i> load_openexr_rgba(const char *file_name);
};


} // namespace psdr
