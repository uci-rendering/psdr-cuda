#include <stdio.h>
#include <cuda.h>
#include <cublas.h>
#include <cuda/host/util.cuh>

#include <thrust/sort.h>

__global__ void cat_meow(){
    printf("Meow!!!\n");
}

void psdr_cuda::meow(int cat) {
    printf("Cats are meowing on GPU...\n");
    cat_meow<<<1,cat>>>();
    cudaDeviceSynchronize();
    printf("Cats are sleeping!\n");
}

__global__ void floatAdd_kernel(float* in1, float* in2, float* out)
{
    int i = threadIdx.x;
    out[i] = in1[i] + in2[i];
}

void psdr_cuda::float_add(float* in1, float* in2, float* out, int size) {
    floatAdd_kernel<<<1, size>>>(in1, in2, out);
}

__global__ void change_val(float* in1)
{
    int i = threadIdx.x;
    // in1[i] = 1.0f;
}

__global__ void test_val(float* in1, float* in2, float* in3, float* in4, float* in5, float* in6, int size)
{
   int i, j;
   float selected1, selected2, selected3, selected4, selected5, selected6;
   for (i = 1; i < size; i++){
      selected1 = in1[i];
      selected2 = in2[i];
      selected3 = in3[i];
      selected4 = in4[i];
      selected5 = in5[i];
      selected6 = in6[i];
      j = i - 1;
      while ((j >= 0) && (selected1+selected2+selected3 < in1[j]+in2[j]+in3[j])) {
         in1[j+1] = in1[j];
         in2[j+1] = in2[j];
         in3[j+1] = in3[j];
         in4[j+1] = in4[j];
         in5[j+1] = in5[j];
         in6[j+1] = in6[j];
         j--;
      }
      in1[j+1] = selected1;
      in2[j+1] = selected2;
      in3[j+1] = selected3;
      in4[j+1] = selected4;
      in5[j+1] = selected5;
      in6[j+1] = selected6;
   }
}

void psdr_cuda::edge_sort(cuda_edge ce, int size) {
	printf("test edge sort\n");
	printf("size=%d\n", size);
	test_val<<<1, 1>>>(ce.p0_x, ce.p0_y, ce.p0_z, ce.p1_x, ce.p1_y, ce.p1_z, size);
}

