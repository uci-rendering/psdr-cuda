#include <vector_functions.h>
#include <vector_types.h>

namespace psdr_cuda {

void init_tree(const float* cdf_x, const float* cdf_y, const float* cdf_z, int dimx, int dimy, int dimz, float* p0_x, float* p0_y, float* p0_z, float* p1_x, float* p1_y, float* p1_z);

void generate_eval_point(int leaf_size, const float* p0_x, const float* p0_y, const float* p0_z, const float* p1_x, const float* p1_y, const float* p1_z, float* outx, float* outy, float* outz, int npass);

int cut_grid(const float* eval_value, float* p0_x, float* p0_y, float* p0_z, 
                                      float* p1_x, float* p1_y, float* p1_z, 
                                      float* eval_x, float* eval_y, float* eval_z, 
                                      int fix_size, float thold, float wt1);

}
