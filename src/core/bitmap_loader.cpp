#ifdef _WIN32
#   define OPENEXR_DLL
#endif
#include <ImathBox.h>
#include <ImfArray.h>
#include <ImfRgba.h>
#include <ImfRgbaFile.h>
#include <psdr/psdr.h>
#include <psdr/core/bitmap_loader.h>

namespace psdr
{

std::pair<Vector4fC, ScalarVector2i> BitmapLoader::load_openexr_rgba(const char *file_name) {
    Imf::RgbaInputFile file(file_name);
    Imath::Box2i dw = file.dataWindow();

    Imf::Array2D<Imf::Rgba> pixels;

    int width = dw.max.x - dw.min.x + 1;
    int height = dw.max.y - dw.min.y + 1;
    pixels.resizeErase(height, width);
    file.setFrameBuffer(&pixels[0][0] - dw.min.x - dw.min.y * width, 1, width);
    file.readPixels(dw.min.y, dw.max.y);

    // std::cout << width << ' ' << height << std::endl;
    int size = width*height;

    std::vector<float> buf[4];
    for ( int i = 0; i < 4; ++i ) buf[i].resize(size);


    int offset = 0;
    for ( int i = 0; i < height; ++i )
        for ( int j = 0; j < width; ++j ) {
            buf[0][offset] = pixels[i][j].r;
            buf[1][offset] = pixels[i][j].g;
            buf[2][offset] = pixels[i][j].b;
            buf[3][offset] = pixels[i][j].a;
            ++offset;
        }

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

} // namespace psdr
