#ifndef _KERNELSETUP_H_
#define _KERNELSETUP_H_

#include <curand_kernel.h> //For curandState type.
#include "../master_def.h" //Number of threads/block are defined at compile time for kernels using a macro definition.

__global__ void setup_kernel( curandState * state, int seed );

template <class Params>
__global__ void setup_kernel(curandState * state, const Params *par )
{
    int tid = threadIdx.x + blockIdx.x*SIM_THREADS;

    curand_init( par->seed+tid, 0, 0, &state[tid] );
} 

__device__ inline int randomSign( curandState & localState ){
	return (curand( &localState ) % 2) ? 1 : -1;
}

__host__ inline int randomSign(){
	return ( rand() % 2) ? 1 : -1;
}

#endif
