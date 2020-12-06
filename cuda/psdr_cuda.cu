#include <optix.h>
#include <psdr/constants.h>
#include "psdr_cuda.h"

extern "C" {
__constant__ Params params;
}

extern "C" __global__ void __raygen__psdr_rg()
{
    const int tid = optixGetLaunchIndex().x;
    optixTrace(
        params.handle,
        make_float3(params.ray_o_x[tid], params.ray_o_y[tid], params.ray_o_z[tid]),
        make_float3(params.ray_d_x[tid], params.ray_d_y[tid], params.ray_d_z[tid]),
        psdr::RayEpsilon,
        params.ray_tmax[tid],
        0.0f,
        OptixVisibilityMask( 1 ),
        OPTIX_RAY_FLAG_NONE,
        0,
        1,
        0
    );
}

extern "C" __global__ void __miss__psdr_ms()
{
    const uint32_t image_index  = optixGetLaunchIndex().x;
    params.tri_index[ image_index ] = -1;
    params.shape_index[ image_index ] = -1;
    params.barycentric_u[ image_index ] = -1.0f;
    params.barycentric_v[ image_index ] = -1.0f;
}
 
extern "C" __global__ void __closesthit__psdr_ch()
{
    HitGroupData* rt_data = (HitGroupData*)optixGetSbtDataPointer();
    const uint32_t image_index  = optixGetLaunchIndex().x;
    params.tri_index[ image_index ] = optixGetPrimitiveIndex() + rt_data->shape_offset;
    params.shape_index[ image_index ] = rt_data->shape_id;
    float2 uv = optixGetTriangleBarycentrics();
    params.barycentric_u[ image_index ] = uv.x;
    params.barycentric_v[ image_index ] = uv.y;
}
