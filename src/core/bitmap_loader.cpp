#include <psdr/psdr.h>
#include <psdr/core/bitmap_loader.h>
#include <misc/Exception.h>

#define TINYEXR_USE_MINIZ 0
#include <psdr/core/miniz.h>
#define TINYEXR_IMPLEMENTATION
#include <psdr/core/tinyexr.h>

namespace psdr
{

std::pair<Vector4fC, ScalarVector2i> BitmapLoader::load_openexr_rgba(const char *file_name) {
    {
        float* out; // width * height * RGBA
        int width;
        int height;
        const char* err = NULL; // or nullptr in C++11

        int ret = LoadEXR(&out, &width, &height, file_name, &err);
        int size = width*height;

        if (ret != TINYEXR_SUCCESS) {
            if (err) {
            fprintf(stderr, "ERR : %s\n", err);
               FreeEXRErrorMessage(err); // release memory of error message.
               PSDR_ASSERT(0);
            }
          } 
        std::vector<float> buf[4];
        for ( int i = 0; i < 4; ++i ) buf[i].resize(size);

        int offset = 0;
        for ( int i = 0; i < height; ++i )
            for ( int j = 0; j < width; ++j ) {
                buf[0][offset] = out[offset*4];
                buf[1][offset] = out[offset*4+1];
                buf[2][offset] = out[offset*4+2];
                buf[3][offset] = out[offset*4+3];
                ++offset;
            }
        free(out); // release memory of image data
        return {
            Vector4fC(
                FloatC::copy(buf[0].data(), size),
                FloatC::copy(buf[1].data(), size),
                FloatC::copy(buf[2].data(), size),
                FloatC::copy(buf[3].data(), size)
            ),
            ScalarVector2i(width, height)
        }; 
    }
}

} // namespace psdr
