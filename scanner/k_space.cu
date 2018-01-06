#include"k_space.cuh"

kSpace::kSpace(const int dim_x, const int dim_y) :dim_x(dim_x), dim_y(dim_y){
	host_space = (cufftDoubleComplex**)malloc(sizeof(cufftDoubleComplex*)*dim_x);
	host_result = (cufftDoubleComplex**)malloc(sizeof(cufftDoubleComplex*)*dim_x);

	host_space[0] = (cufftDoubleComplex *)malloc(dim_y*dim_x*sizeof(cufftDoubleComplex));
	host_result[0] = (cufftDoubleComplex *)malloc(dim_y*dim_x*sizeof(cufftDoubleComplex));

	for (int i = 1; i < dim_x; i++){
		host_space[i] = host_space[i - 1] + dim_y;
		host_result[i] = host_result[i - 1] + dim_y;
	}

	host_space_pitch = dim_y * sizeof(cufftDoubleComplex);
	cudaMallocPitch(&dev_space, &dev_space_pitch, dim_y*sizeof(cufftDoubleComplex), dim_x);
}

size_t kSpace::index(int x, int y){
	return y + x*dim_y;
}

double kSpace::get_Mx(int x, int y){
	return host_space[x][y].x;
}

double kSpace::get_My(int x, int y){
	return host_space[x][y].y;
}

void kSpace::set_Mx(int x, int y, double val){
	host_space[x][y].x += val;
}
void kSpace::set_My(int x, int y, double val){
	host_space[x][y].y += val;
}
void kSpace::get_fft(){
	cudaMemcpy2D(dev_space,
		dev_space_pitch,
		host_space[0],
		host_space_pitch,
		dim_y*sizeof(cufftDoubleComplex),
		dim_x,
		cudaMemcpyHostToDevice);

	printf("Starting FFT Process\n");
	cufftHandle plan;
	cufftPlan2d(&plan, dim_x, dim_y, CUFFT_Z2Z);
	printf("Plan Built\n");
	cufftExecZ2Z(plan, dev_space, dev_space, CUFFT_INVERSE);
	cudaDeviceSynchronize();
	cudaMemcpy2D(host_result[0], host_space_pitch, dev_space, dev_space_pitch,
		dim_y * sizeof(cufftDoubleComplex), dim_x, cudaMemcpyDeviceToHost);

	recorder image_k_r("k_space_r");
	image_k_r.save_csv_complex<cufftDoubleComplex>(host_space, 'r', dim_x, dim_y);

	recorder image_k_i("k_space_i");
	image_k_i.save_csv_complex<cufftDoubleComplex>(host_space, 'i', dim_x, dim_y);

	recorder image_k("k_space");
	image_k.save_image_complex<cufftDoubleComplex>(host_space, 'm', dim_x, dim_y);

	recorder image("image");
	image.save_image_complex<cufftDoubleComplex>(host_result, 'm', dim_x, dim_y);
}

