#include "sequence.cuh"

__global__ void make_tensors_kernel(Sequence** seq){
	const unsigned int tid = threadIdx.x + blockIdx.x*blockDim.x;
	//Ensures that we don't over-index the G_tensor if more threads used than time-steps required.
	if (tid < (*seq)->getSteps()){
		real element[9] = {0};
		real RF_element[3] = {0};

		(*seq)->get_G_tensor(tid, element, RF_element);

		for (int i = 0; i < 9; i++){
			(*seq)->G_tensor_t[9*tid+i] = element[i];
		}
		for (int i = 0; i < 3; i++){
			(*seq)->RF_tensor_t[3*tid+i] = RF_element[i];
		}
	}
}

__global__ void assign_tensors_ptr(Sequence** seq, real* G_tensor_t_ptr, real* RF_tensor_t_ptr){
	(*seq)->G_tensor_t = G_tensor_t_ptr;
	(*seq)->RF_tensor_t = RF_tensor_t_ptr;
}
