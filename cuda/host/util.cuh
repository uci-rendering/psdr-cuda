namespace psdr_cuda {

	struct cuda_edge
	{
		float* p0_x;
		float* p0_y;
		float* p0_z;

		float* p1_x;
		float* p1_y;
		float* p1_z;
	};

	void meow(int cat);
	void float_add(float* in1, float* in2, float* out, int size);

	void edge_sort(cuda_edge ce, int size);
}