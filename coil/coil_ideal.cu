/*
Ideal coil with perfect homogeneity.
*/

#include "coil_ideal.cuh"

__global__ void Coil_Ideal_GPU(Coil_Ideal** obj_ptr);

__device__ __host__ Coil_Ideal::Coil_Ideal(){
#ifndef __CUDA_ARCH__
	cudaMalloc(&dev_ptr, sizeof(Coil_Ideal**));
	Coil_Ideal_GPU << <1, 1 >> >(dev_ptr);
#endif
}

__device__ __host__ Vector3 Coil_Ideal::getField(const Sequence* sequence, Vector3 r, real time) const{
	return sequence->getG(r, time);
}

__device__ Vector3 Coil_Ideal::getField(const Sequence* sequence, Vector3 r, real time, real* G_tensor, real* RF_tensor, int thread) const{
	return sequence->getG(r, time, G_tensor, RF_tensor, thread);
}

__host__ const Coil** Coil_Ideal::devPointer() const{
	return (const Coil**) dev_ptr;
}

__global__ void Coil_Ideal_GPU(Coil_Ideal** obj_ptr){
	*obj_ptr = new Coil_Ideal;
}
