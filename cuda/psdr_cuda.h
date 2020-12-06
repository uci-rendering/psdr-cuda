#include <stdint.h>

struct Params
{
    const float  *ray_o_x, *ray_o_y, *ray_o_z;
    const float  *ray_d_x, *ray_d_y, *ray_d_z;
    const float  *ray_tmax;

    int          *tri_index;
    int          *shape_index;
    float        *barycentric_u, *barycentric_v;

    OptixTraversableHandle handle;
};

struct RayGenData
{

};

struct MissData
{

};

struct HitGroupData
{
	int shape_offset;
	int shape_id;
};
