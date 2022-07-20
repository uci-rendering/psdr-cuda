#include <stdio.h>
#include <cuda.h>
#include <cublas.h>
#include <cub/cub.cuh>
#include <cuda/host/NEE.cuh>


template<unsigned int N>
static __host__ __device__ __inline__ unsigned int tea( unsigned int val0, unsigned int val1 )
{
  unsigned int v0 = val0;
  unsigned int v1 = val1;
  unsigned int s0 = 0;

  for( unsigned int n = 0; n < N; n++ )
  {
    s0 += 0x9e3779b9;
    v0 += ((v1<<4)+0xa341316c)^(v1+s0)^((v1>>5)+0xc8013ea4);
    v1 += ((v0<<4)+0xad90777d)^(v0+s0)^((v0>>5)+0x7e95761e);
  }

  return v0;
}

// Generate random unsigned int in [0, 2^24)
static __host__ __device__ __inline__ unsigned int lcg(unsigned int &prev)
{
  const unsigned int LCG_A = 1664525u;
  const unsigned int LCG_C = 1013904223u;
  prev = (LCG_A * prev + LCG_C);
  return prev & 0x00FFFFFF;
}

// Generate random float in [0, 1)
static __host__ __device__ __inline__ float rnd(unsigned int &prev)
{
  return ((float) lcg(prev) / (float) 0x01000000);
}

__global__ void build_tree(float* p0_x, float* p0_y, float* p0_z, float* p1_x, float* p1_y, float* p1_z, const float* cdf_x, int x_num, int y_num, int z_num) {
    int tid=blockIdx.x * blockDim.x + threadIdx.x;
    if (tid < x_num*y_num*z_num) {
        float dimy = 1.0f/y_num;
        float dimz = 1.0f/z_num;

        int i = tid % x_num;
        int j = tid / x_num % y_num;
        int k = tid / (x_num*y_num) % z_num;

        p0_x[tid] = (i!=0) ? cdf_x[i-1] : 0.0f;
        p0_y[tid] = j*dimy;
        p0_z[tid] = k*dimz;

        p1_x[tid] = cdf_x[i];
        p1_y[tid] = (j+1.f)*dimy;
        p1_z[tid] = (k+1.f)*dimz;
    }
}

__global__ void gen_eval_point( const float* p0_x, const float* p0_y, const float* p0_z, const float* p1_x, const float* p1_y, const float* p1_z, float* outx, float* outy, float* outz, int npass) {
    int bid=blockIdx.x;
    int tid=threadIdx.x;
    int i = bid * blockDim.x + tid;

    unsigned int seed = tea<4>( i, npass );

    float midx = (p1_x[bid]+p0_x[bid])/2.0f;
    float midy = (p1_y[bid]+p0_y[bid])/2.0f;
    float midz = (p1_z[bid]+p0_z[bid])/2.0f;

    float rndx = rnd(seed);
    float rndy = rnd(seed);
    float rndz = rnd(seed);

    switch(tid) {
        case 0:
            outx[i] = p0_x[bid];
            outy[i] = p0_y[bid];
            outz[i] = p0_z[bid];
            break;
        case 1:
            outx[i] = p0_x[bid];
            outy[i] = p0_y[bid];
            outz[i] = p1_z[bid];
            break;
        case 2:
            outx[i] = p0_x[bid];
            outy[i] = p0_y[bid];
            outz[i] = p0_z[bid]+midz;
            break;

        case 3:
            outx[i] = p0_x[bid];
            outy[i] = p1_y[bid];
            outz[i] = p0_z[bid];
            break;
        case 4:
            outx[i] = p0_x[bid];
            outy[i] = p1_y[bid];
            outz[i] = p1_z[bid];
            break;
        case 5:
            outx[i] = p0_x[bid];
            outy[i] = p1_y[bid];
            outz[i] = p0_z[bid]+midz;
            break;

        case 6:
            outx[i] = p0_x[bid];
            outy[i] = p0_y[bid]+midy;
            outz[i] = p0_z[bid];
            break;
        case 7:
            outx[i] = p0_x[bid];
            outy[i] = p0_y[bid]+midy;
            outz[i] = p1_z[bid];
            break;

        case 8:
            outx[i] = p1_x[bid];
            outy[i] = p0_y[bid];
            outz[i] = p0_z[bid];
            break;
        case 9:
            outx[i] = p1_x[bid];
            outy[i] = p0_y[bid];
            outz[i] = p1_z[bid];
            break;
        case 10:
            outx[i] = p1_x[bid];
            outy[i] = p0_y[bid];
            outz[i] = p0_z[bid]+midz;
            break;

        case 11:
            outx[i] = p1_x[bid];
            outy[i] = p1_y[bid];
            outz[i] = p0_z[bid];
            break;
        case 12:
            outx[i] = p1_x[bid];
            outy[i] = p1_y[bid];
            outz[i] = p1_z[bid];
            break;
        case 13:
            outx[i] = p1_x[bid];
            outy[i] = p1_y[bid];
            outz[i] = p0_z[bid]+midz;
            break;

        case 14:
            outx[i] = p1_x[bid];
            outy[i] = p0_y[bid]+midy;
            outz[i] = p0_z[bid];
            break;
        case 15:
            outx[i] = p1_x[bid];
            outy[i] = p0_y[bid]+midy;
            outz[i] = p1_z[bid];
            break;

        case 16:
            outx[i] = p0_x[bid]+midx;
            outy[i] = p0_y[bid];
            outz[i] = p0_z[bid];
            break;
        case 17:
            outx[i] = p0_x[bid]+midx;
            outy[i] = p0_y[bid];
            outz[i] = p1_z[bid];
            break;
        case 18:
            outx[i] = p0_x[bid]+midx;
            outy[i] = p1_y[bid];
            outz[i] = p0_z[bid];
            break;
        case 19:
            outx[i] = p0_x[bid]+midx;
            outy[i] = p1_y[bid];
            outz[i] = p1_z[bid];
            break;
        
        // random sample for grid
        case 20:
            outx[i] = p0_x[bid] + (midx-p0_x[bid])*rndx;
            outy[i] = p0_y[bid] + (midy-p0_y[bid])*rndy;
            outz[i] = p0_z[bid] + (midz-p0_z[bid])*rndz;
            break;
        case 21:
            outx[i] = p0_x[bid] + (midx-p0_x[bid])*rndx;
            outy[i] = p0_y[bid] + (midy-p0_y[bid])*rndy;
            outz[i] = midz      + (midz-p0_z[bid])*rndz;
            break;
        case 22:
            outx[i] = p0_x[bid] + (midx-p0_x[bid])*rndx;
            outy[i] = midy      + (midy-p0_y[bid])*rndy;
            outz[i] = p0_z[bid] + (midz-p0_z[bid])*rndz;
            break;
        case 23:
            outx[i] = p0_x[bid] + (midx-p0_x[bid])*rndx;
            outy[i] = midy      + (midy-p0_y[bid])*rndy;
            outz[i] = midz      + (midz-p0_z[bid])*rndz;
            break;

        case 24:
            outx[i] = midx      + (midx-p0_x[bid])*rndx;
            outy[i] = p0_y[bid] + (midy-p0_y[bid])*rndy;
            outz[i] = p0_z[bid] + (midz-p0_z[bid])*rndz;
            break;
        case 25:
            outx[i] = midx      + (midx-p0_x[bid])*rndx;
            outy[i] = p0_y[bid] + (midy-p0_y[bid])*rndy;
            outz[i] = midz      + (midz-p0_z[bid])*rndz;
            break;
        case 26:
            outx[i] = midx      + (midx-p0_x[bid])*rndx;
            outy[i] = midy      + (midy-p0_y[bid])*rndy;
            outz[i] = p0_z[bid] + (midz-p0_z[bid])*rndz;
            break;
        case 27:
            outx[i] = midx      + (midx-p0_x[bid])*rndx;
            outy[i] = midy      + (midy-p0_y[bid])*rndy;
            outz[i] = midz      + (midz-p0_z[bid])*rndz;
            break;

        case 28:
            outx[i] = p0_x[bid] + (midx-p0_x[bid])*(1.0f-rndx);
            outy[i] = p0_y[bid] + (midy-p0_y[bid])*(1.0f-rndy);
            outz[i] = p0_z[bid] + (midz-p0_z[bid])*(1.0f-rndz);
            break;
        case 29:
            outx[i] = p0_x[bid] + (midx-p0_x[bid])*(1.0f-rndx);
            outy[i] = p0_y[bid] + (midy-p0_y[bid])*(1.0f-rndy);
            outz[i] = midz      + (midz-p0_z[bid])*(1.0f-rndz);
            break;
        case 30:
            outx[i] = p0_x[bid] + (midx-p0_x[bid])*(1.0f-rndx);
            outy[i] = midy      + (midy-p0_y[bid])*(1.0f-rndy);
            outz[i] = p0_z[bid] + (midz-p0_z[bid])*(1.0f-rndz);
            break;
        case 31:
            outx[i] = p0_x[bid] + (midx-p0_x[bid])*(1.0f-rndx);
            outy[i] = midy      + (midy-p0_y[bid])*(1.0f-rndy);
            outz[i] = midz      + (midz-p0_z[bid])*(1.0f-rndz);
            break;

        case 32:
            outx[i] = midx      + (midx-p0_x[bid])*(1.0f-rndx);
            outy[i] = p0_y[bid] + (midy-p0_y[bid])*(1.0f-rndy);
            outz[i] = p0_z[bid] + (midz-p0_z[bid])*(1.0f-rndz);
            break;
        case 33:
            outx[i] = midx      + (midx-p0_x[bid])*(1.0f-rndx);
            outy[i] = p0_y[bid] + (midy-p0_y[bid])*(1.0f-rndy);
            outz[i] = midz      + (midz-p0_z[bid])*(1.0f-rndz);
            break;
        case 34:
            outx[i] = midx      + (midx-p0_x[bid])*(1.0f-rndx);
            outy[i] = midy      + (midy-p0_y[bid])*(1.0f-rndy);
            outz[i] = p0_z[bid] + (midz-p0_z[bid])*(1.0f-rndz);
            break;
        case 35:
            outx[i] = midx      + (midx-p0_x[bid])*(1.0f-rndx);
            outy[i] = midy      + (midy-p0_y[bid])*(1.0f-rndy);
            outz[i] = midz      + (midz-p0_z[bid])*(1.0f-rndz);
            break;
    }
}

__global__ void cut_tree(const float* eval_value, float* p0_x, float* p0_y, float* p0_z, float* p1_x, float* p1_y, float* p1_z,
                                                  float* eval_x, float* eval_y, float* eval_z, 
                                                  float thold, float wt1, int bound_size, int* global_index) {
    unsigned int seed = tea<4>( blockIdx.x, threadIdx.x);
    int bid = blockIdx.x * blockDim.x + threadIdx.x;
    if (bid < bound_size) {
        int eid = bid*36;
        float leaf_size =   (p1_x[bid]-p0_x[bid])*
                            (p1_y[bid]-p0_y[bid])*
                            (p1_z[bid]-p0_z[bid]);

        float lag_val = (eval_value[eid] + eval_value[eid+1] + eval_value[eid+3] + eval_value[eid+4] + eval_value[eid+8] + eval_value[eid+9] + eval_value[eid+11] + eval_value[eid+12]) / 8.0f;
        float MC_val = (eval_value[eid+20] + eval_value[eid+21] + eval_value[eid+22] + eval_value[eid+23] + 
                        eval_value[eid+24] + eval_value[eid+25] + eval_value[eid+26] + eval_value[eid+27] +
                        eval_value[eid+28] + eval_value[eid+29] + eval_value[eid+30] + eval_value[eid+31] +
                        eval_value[eid+32] + eval_value[eid+33] + eval_value[eid+34] + eval_value[eid+35]) / 16.0f;
        
        float error = abs(lag_val-MC_val) * leaf_size * wt1;
        if (error> thold) {
            int current_index = atomicAdd(&global_index[0], 1) + bound_size;
            float dim1a  = (eval_value[eid]+eval_value[eid+1]+eval_value[eid+3]+eval_value[eid+4]+eval_value[eid+16]+eval_value[eid+17]+eval_value[eid+18]+eval_value[eid+19]) / 8.0f;
            float dim1b  = (eval_value[eid+8]+eval_value[eid+9]+eval_value[eid+11]+eval_value[eid+12]+eval_value[eid+16]+eval_value[eid+17]+eval_value[eid+18]+eval_value[eid+19]) / 8.0f;

            float dim2a  = (eval_value[eid]+eval_value[eid+1]+eval_value[eid+6]+eval_value[eid+7]+eval_value[eid+8]+eval_value[eid+9]+eval_value[eid+14]+eval_value[eid+15]) / 8.0f;
            float dim2b  = (eval_value[eid+3]+eval_value[eid+4]+eval_value[eid+6]+eval_value[eid+7]+eval_value[eid+11]+eval_value[eid+12]+eval_value[eid+14]+eval_value[eid+15]) / 8.0f;

            float dim3a  = (eval_value[eid]+eval_value[eid+2]+eval_value[eid+3]+eval_value[eid+5]+eval_value[eid+8]+eval_value[eid+10]+eval_value[eid+11]+eval_value[eid+13]) / 8.0f;
            float dim3b  = (eval_value[eid+1]+eval_value[eid+2]+eval_value[eid+4]+eval_value[eid+5]+eval_value[eid+9]+eval_value[eid+10]+eval_value[eid+12]+eval_value[eid+13]) / 8.0f;

            float MC1a  = (eval_value[eid+20]+eval_value[eid+21]+eval_value[eid+22]+eval_value[eid+23]+eval_value[eid+20+8]+eval_value[eid+21+8]+eval_value[eid+22+8]+eval_value[eid+23+8]) / 8.0f;
            float MC1b  = (eval_value[eid+24]+eval_value[eid+25]+eval_value[eid+26]+eval_value[eid+27]+eval_value[eid+24+8]+eval_value[eid+25+8]+eval_value[eid+26+8]+eval_value[eid+27+8]) / 8.0f;

            float MC2a  = (eval_value[eid+20]+eval_value[eid+21]+eval_value[eid+24]+eval_value[eid+25]+eval_value[eid+20+8]+eval_value[eid+21+8]+eval_value[eid+24+8]+eval_value[eid+25+8]) / 8.0f;
            float MC2b  = (eval_value[eid+22]+eval_value[eid+23]+eval_value[eid+26]+eval_value[eid+27]+eval_value[eid+22+8]+eval_value[eid+23+8]+eval_value[eid+26+8]+eval_value[eid+27+8]) / 8.0f;

            float MC3a  = (eval_value[eid+20]+eval_value[eid+22]+eval_value[eid+24]+eval_value[eid+26]+eval_value[eid+20+8]+eval_value[eid+22+8]+eval_value[eid+24+8]+eval_value[eid+26+8]) / 8.0f;
            float MC3b  = (eval_value[eid+21]+eval_value[eid+23]+eval_value[eid+25]+eval_value[eid+27]+eval_value[eid+21+8]+eval_value[eid+23+8]+eval_value[eid+25+8]+eval_value[eid+27+8]) / 8.0f;;

            float error0 = abs(abs(dim1a-MC1a) - abs(dim1b-MC1b));
            float error1 = abs(abs(dim2a-MC2a) - abs(dim2b-MC2b));
            float error2 = abs(abs(dim3a-MC3a) - abs(dim3b-MC3b));
            if (error0 > error1 && error0 > error2) {
                p1_x[current_index] = p1_x[bid];
                p1_y[current_index] = p1_y[bid];
                p1_z[current_index] = p1_z[bid];
                p0_x[current_index] = (p1_x[bid] + p0_x[bid]) / 2.0f;
                p0_y[current_index] = p0_y[bid];
                p0_z[current_index] = p0_z[bid];
                p1_x[bid] = (p1_x[bid] + p0_x[bid]) / 2.0f;
            } else if (error1 > error0 && error1 > error2) {
                p1_x[current_index] = p1_x[bid];
                p1_y[current_index] = p1_y[bid];
                p1_z[current_index] = p1_z[bid];
                p0_y[current_index] = (p1_y[bid] + p0_y[bid]) / 2.0f;
                p0_x[current_index] = p0_x[bid];
                p0_z[current_index] = p0_z[bid];
                p1_y[bid] = (p1_y[bid] + p0_y[bid]) / 2.0f;
            } else if (error2 > error0 && error2 > error1) {
                p1_x[current_index] = p1_x[bid];
                p1_y[current_index] = p1_y[bid];
                p1_z[current_index] = p1_z[bid];
                p0_z[current_index] = (p1_z[bid] + p0_z[bid]) / 2.0f;
                p0_x[current_index] = p0_x[bid];
                p0_y[current_index] = p0_y[bid];
                p1_z[bid] = (p1_z[bid] + p0_z[bid]) / 2.0f;
            } else {
                float rnd_num = rnd(seed) * 3.0f;
                if (rnd_num < 1.0f) {
                    p1_x[current_index] = p1_x[bid];
                    p1_y[current_index] = p1_y[bid];
                    p1_z[current_index] = p1_z[bid];
                    p0_x[current_index] = (p1_x[bid] + p0_x[bid]) / 2.0f;
                    p0_y[current_index] = p0_y[bid];
                    p0_z[current_index] = p0_z[bid];
                    p1_x[bid] = (p1_x[bid] + p0_x[bid]) / 2.0f;
                } else if (rnd_num < 2.0f) {
                    p1_x[current_index] = p1_x[bid];
                    p1_y[current_index] = p1_y[bid];
                    p1_z[current_index] = p1_z[bid];
                    p0_y[current_index] = (p1_y[bid] + p0_y[bid]) / 2.0f;
                    p0_x[current_index] = p0_x[bid];
                    p0_z[current_index] = p0_z[bid];
                    p1_y[bid] = (p1_y[bid] + p0_y[bid]) / 2.0f;
                } else {
                    p1_x[current_index] = p1_x[bid];
                    p1_y[current_index] = p1_y[bid];
                    p1_z[current_index] = p1_z[bid];
                    p0_z[current_index] = (p1_z[bid] + p0_z[bid]) / 2.0f;
                    p0_x[current_index] = p0_x[bid];
                    p0_y[current_index] = p0_y[bid];
                    p1_z[bid] = (p1_z[bid] + p0_z[bid]) / 2.0f;
                }
            }

            // float error0 = abs(abs(dim1a-MC1a) + abs(dim1b-MC1b));
            // float error1 = abs(abs(dim2a-MC2a) + abs(dim2b-MC2b));
            // float error2 = abs(abs(dim3a-MC3a) + abs(dim3b-MC3b));
            // if (error0 < error1 && error0 < error2) {
            //     p1_x[current_index] = p1_x[bid];
            //     p1_y[current_index] = p1_y[bid];
            //     p1_z[current_index] = p1_z[bid];
            //     p0_x[current_index] = (p1_x[bid] + p0_x[bid]) / 2.0f;
            //     p0_y[current_index] = p0_y[bid];
            //     p0_z[current_index] = p0_z[bid];
            //     p1_x[bid] = (p1_x[bid] + p0_x[bid]) / 2.0f;
            // } else if (error1 < error0 && error1 < error2) {
            //     p1_x[current_index] = p1_x[bid];
            //     p1_y[current_index] = p1_y[bid];
            //     p1_z[current_index] = p1_z[bid];
            //     p0_y[current_index] = (p1_y[bid] + p0_y[bid]) / 2.0f;
            //     p0_x[current_index] = p0_x[bid];
            //     p0_z[current_index] = p0_z[bid];
            //     p1_y[bid] = (p1_y[bid] + p0_y[bid]) / 2.0f;
            // } else if (error2 < error0 && error2 < error1) {
            //     p1_x[current_index] = p1_x[bid];
            //     p1_y[current_index] = p1_y[bid];
            //     p1_z[current_index] = p1_z[bid];
            //     p0_z[current_index] = (p1_z[bid] + p0_z[bid]) / 2.0f;
            //     p0_x[current_index] = p0_x[bid];
            //     p0_y[current_index] = p0_y[bid];
            //     p1_z[bid] = (p1_z[bid] + p0_z[bid]) / 2.0f;
            // } else {
            //     float rnd_num = rnd(seed) * 3.0f;
            //     if (rnd_num < 1.0f) {
            //         p1_x[current_index] = p1_x[bid];
            //         p1_y[current_index] = p1_y[bid];
            //         p1_z[current_index] = p1_z[bid];
            //         p0_x[current_index] = (p1_x[bid] + p0_x[bid]) / 2.0f;
            //         p0_y[current_index] = p0_y[bid];
            //         p0_z[current_index] = p0_z[bid];
            //         p1_x[bid] = (p1_x[bid] + p0_x[bid]) / 2.0f;
            //     } else if (rnd_num < 2.0f) {
            //         p1_x[current_index] = p1_x[bid];
            //         p1_y[current_index] = p1_y[bid];
            //         p1_z[current_index] = p1_z[bid];
            //         p0_y[current_index] = (p1_y[bid] + p0_y[bid]) / 2.0f;
            //         p0_x[current_index] = p0_x[bid];
            //         p0_z[current_index] = p0_z[bid];
            //         p1_y[bid] = (p1_y[bid] + p0_y[bid]) / 2.0f;
            //     } else {
            //         p1_x[current_index] = p1_x[bid];
            //         p1_y[current_index] = p1_y[bid];
            //         p1_z[current_index] = p1_z[bid];
            //         p0_z[current_index] = (p1_z[bid] + p0_z[bid]) / 2.0f;
            //         p0_x[current_index] = p0_x[bid];
            //         p0_y[current_index] = p0_y[bid];
            //         p1_z[bid] = (p1_z[bid] + p0_z[bid]) / 2.0f;
            //     }
            // }


        }
    }
}

void psdr_cuda::init_tree(const float* cdf_x, const float* cdf_y, const float* cdf_z, int dimx, int dimy, int dimz, float* p0_x, float* p0_y, float* p0_z, float* p1_x, float* p1_y, float* p1_z) {
    int thread_size = 64;
    int block_size = (dimx*dimy*dimz-1) / thread_size + 1;
    build_tree<<<block_size,thread_size>>>( p0_x,  p0_y,  p0_z,  p1_x,  p1_y,  p1_z, cdf_x, dimx, dimy, dimz);
}

void psdr_cuda::generate_eval_point(int leaf_size, const float* p0_x, const float* p0_y, const float* p0_z, const float* p1_x, const float* p1_y, const float* p1_z, float* outx, float* outy, float* outz, int npass) {
    gen_eval_point<<<leaf_size, 36>>>(p0_x, p0_y, p0_z, p1_x, p1_y, p1_z, outx, outy, outz, npass);
}

int psdr_cuda::cut_grid(const float* eval_value, float* p0_x, float* p0_y, float* p0_z, float* p1_x, float* p1_y, float* p1_z,
                                                 float* eval_x, float* eval_y, float* eval_z, 
                                                 int fix_size, float thold, float wt1) {
    int thread_size = 512;
    int block_size = (fix_size-1) / thread_size + 1;

    int* global_index;
    cudaMalloc(&global_index, sizeof(int));
    cudaMemset(global_index, 0, sizeof(int));
    cut_tree<<<block_size, thread_size>>>(eval_value, p0_x,  p0_y,  p0_z,  p1_x,  p1_y, p1_z, eval_x, eval_y, eval_z,
                                            thold, wt1, fix_size, global_index);
    
    int app_size;
    cudaMemcpy(&app_size,global_index,sizeof(int),cudaMemcpyDeviceToHost);
    return app_size;
}
