#include "kernelSetup.cuh"

__global__ void setup_kernel( curandState * state, int seed )
{
    int tid = threadIdx.x + blockIdx.x*SIM_THREADS;
    curand_init( seed, tid, 0, &state[tid] );
}
