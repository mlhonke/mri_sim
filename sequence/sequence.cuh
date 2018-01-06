/*
Sequence: Superclass of all scanner sequences. Sequences are very diverse 
and therefore this superclass describes the basic methods required by a 
sequence to run successfully in this simulator.
Author: Michael Honke
Date: Oct. 14, 2016
*/

#ifndef SEQUENCE_CUH
#define SEQUENCE_CUH

#include<cuda_runtime.h>

#include "../master_def.h"
#include "../util/vector3.cuh"
#include "../util/misc.cuh"
#include "device_launch_parameters.h"
#include "../params/simuParams.cuh"
#include "pulses.cuh"
#include"../util/recorder.h"

class Sequence;

__global__ void make_tensors_kernel(Sequence** seq);
__global__ void assign_tensors_ptr(Sequence** seq, real* G_tensor_t_ptr, real* RF_tensor_t_ptr);

class Sequence{
protected:
	int n_sub_sequences;
	int readSteps;
	int readFactor;
	int steps;
	int num_pulses;
	int local_res_x;
	SimuParams* par;
	Pulse** pulse;
	real* G_tensor_t_devptr;
	real* RF_tensor_t_devptr;
	real* G_tensor_shared;
	real* RF_tensor_shared;
	real TR;
	int phase_enc_offset;
	int thread_assign_tensors;

	Sequence** dev_ptr;

	__device__ __host__ Sequence():thread_assign_tensors(ceil(G_SHARED_SIZE / (real) SIM_THREADS)){
		n_sub_sequences = 0;
	}

	__device__ __host__ Sequence(int n_sub_sequences)
		: Sequence(n_sub_sequences, 0, 0, 0, 0){
	}

	__device__ __host__ Sequence(int n_sub_sequences, int num_pulses, real phase_enc_offset, int _local_res_x, SimuParams* par)
		: n_sub_sequences(n_sub_sequences),
		  num_pulses(num_pulses),
		  phase_enc_offset(phase_enc_offset),par(par),
		  TR(par->TR),
		  thread_assign_tensors(ceil(G_SHARED_SIZE / (real) SIM_THREADS))
	{
		_local_res_x == 0 ? local_res_x = par->res_x : local_res_x = _local_res_x;//Check if a subsequence or not.
		steps = (TR/par->timestep)*local_res_x;
		pulse = new Pulse*[num_pulses];
		printf("Number of steps in sequence: %d\n", steps);

#if defined(ALLOC_G) && not defined (__CUDA_ARCH__)
		cudaMalloc((void**) &G_tensor_t_devptr, sizeof(real)*steps*9);
		cudaMalloc((void**) &RF_tensor_t_devptr, sizeof(real)*steps*3);

		//G_tensor_t = (real*) malloc(sizeof(real)*steps*9);
		//RF_tensor_t = (real*) malloc(sizeof(real)*steps*3);
		//printf("Pre-field allocation: %f MB\n", sizeof(real)*steps*12 / 1000000.0);
#endif
	}

	__host__ __device__ ~Sequence(void){
#if defined(ALLOC_G) && not defined (__CUDA_ARCH__)
		cudaFree(G_tensor_t_devptr);
		cudaFree(RF_tensor_t_devptr);
#endif
		free(pulse);
	}

	__host__ void make_tensors(){
		printf("Allocating tensors on device with %f MB of memory.\n", sizeof(real)*steps*12 / 1000000.0);
		assign_tensors_ptr<<<1,1>>>(dev_ptr, G_tensor_t_devptr, RF_tensor_t_devptr);
		int tensor_threads = 512;
		int tensor_blocks = ceil(steps/(real)tensor_threads);
		make_tensors_kernel<<<tensor_blocks, tensor_threads>>>(dev_ptr);
		printf("Done allocating tensors on device.\n");
		safe_cuda(cudaDeviceSynchronize());
	}

public:
	real* G_tensor_t;
	real* RF_tensor_t;

	virtual __host__ const Sequence* getSubSequences(int i) const = 0;
	virtual __device__ __host__ int get_k_start() const = 0;
	virtual __device__ __host__ int get_k_end() const = 0;
	virtual __device__ __host__ real getReadStart(real time) const = 0;
	virtual __device__ __host__ real getReadFinish(real time) const = 0;
	virtual __device__ __host__ Vector3 getK(int time) const = 0;

	virtual __device__ Vector3 getG(Vector3 r, real time, real* G_tensor_shared, real* RF_tensor_shared, int thread_id) const{
		//return getG(r, time);
		int timestep = time / par->timestep;
		int base_index = timestep * 9;
		int base_index_RF = timestep * 3;
		int index;

		if (timestep % G_SHARED_SIZE == 0){
			__syncthreads();
			for (int i = 0; i < thread_assign_tensors; i++){
				index = thread_id + SIM_THREADS*i;
				if ( index <= G_SHARED_SIZE){
					for (int j = 0; j < 9; j++){
						G_tensor_shared[index*9 + j] = G_tensor_t[base_index + index*9 + j];
					}
					for (int j = 0; j < 3; j++){
						RF_tensor_shared[index*3 + j] = RF_tensor_t[base_index_RF + index*3 + j];
					}
				}
			}

//			if (thread_id == 0){
//				for (int i = 0; i <= G_SHARED_SIZE; i++){
//					for (int j = 0; j < 9; j++){
//						G_tensor_shared[i*9 + j] = G_tensor_t[base_index + i*9 + j];
//					}
//					for (int j = 0; j < 3; j++){
//						RF_tensor_shared[i*3 + j] = RF_tensor_t[base_index_RF + i*3 + j];
//					}
//				}
//			}
		}

		base_index = (timestep % G_SHARED_SIZE) * 9;
		base_index_RF = (timestep % G_SHARED_SIZE) * 3;

		return Vector3(
				G_tensor_shared[base_index + 0]*r.x+
				G_tensor_shared[base_index + 1]*r.y+
				G_tensor_shared[base_index + 2]*r.z+RF_tensor_shared[base_index_RF + 0],
				G_tensor_shared[base_index + 3]*r.x+
				G_tensor_shared[base_index + 4]*r.y+
				G_tensor_shared[base_index + 5]*r.z+RF_tensor_shared[base_index_RF + 1],
				G_tensor_shared[base_index + 6]*r.x+
				G_tensor_shared[base_index + 7]*r.y+
				G_tensor_shared[base_index + 8]*r.z+RF_tensor_shared[base_index_RF + 2]
				);
	}

#ifdef ALLOC_G
	virtual __host__ __device__ Vector3 getG(Vector3 r, real time) const{
		int timestep = time / par->timestep;
		int base_index = timestep * 9;
		int base_index_RF = timestep * 3;

		return Vector3(
				G_tensor_t[base_index + 0]*r.x+
				G_tensor_t[base_index + 1]*r.y+
				G_tensor_t[base_index + 2]*r.z+RF_tensor_t[base_index_RF + 0],
				G_tensor_t[base_index + 3]*r.x+
				G_tensor_t[base_index + 4]*r.y+
				G_tensor_t[base_index + 5]*r.z+RF_tensor_t[base_index_RF + 1],
				G_tensor_t[base_index + 6]*r.x+
				G_tensor_t[base_index + 7]*r.y+
				G_tensor_t[base_index + 8]*r.z+RF_tensor_t[base_index_RF + 2]
				);
	}
#else
	//Should be called getB1 or getBgrad...
	virtual __device__ __host__ Vector3 getG(Vector3 r, real time) const{
		int seg = (int) (time / TR) + phase_enc_offset;
		real seg_time = time - seg*TR + phase_enc_offset*TR;
		Vector3 G(0.0,0.0,0.0);
		for (int i = 0; i<num_pulses; i++){
			if (pulse[i]->on(seg_time))
				G += pulse[i]->out(time, seg_time, r, seg);
		}
		return G;
	}
#endif

	virtual __host__ Vector3 getGCPU(Vector3 r, real time) const{
		int seg = (int) (time / TR) + phase_enc_offset;
		real seg_time = time - seg*TR + phase_enc_offset*TR;
		Vector3 G(0.0,0.0,0.0);
		for (int i = 0; i<num_pulses; i++){
			if (pulse[i]->on(seg_time))
				G += pulse[i]->out(time, seg_time, r, seg);
		}
		return G;
	}

	virtual __device__ __host__ void get_G_tensor(int step, real element[9], real RF_element[3]){
		real time = step * par->timestep;
		int seg = (int) (time / TR) + phase_enc_offset;
		real seg_time = time - seg*TR + phase_enc_offset*TR;
		Vector3 G(0.0,0.0,0.0);
		for (int i = 0; i<num_pulses; i++){
			if (pulse[i]->on(seg_time)){
				G = pulse[i]->out(time, seg_time, seg);
				if (pulse[i]->rDir == Vector3(1,0,0)){
					element[0] += G.x;
					element[3] += G.y;
					element[6] += G.z;
				}
				if (pulse[i]->rDir == Vector3(0,1,0)){
					element[1] += G.x;
					element[4] += G.y;
					element[7] += G.z;
				}
				if (pulse[i]->rDir == Vector3(0,0,1)){
					element[2] += G.x;
					element[5] += G.y;
					element[8] += G.z;
				}
				if (pulse[i]->rDir == Vector3(0,0,0)){
					RF_element[0] += G.x;
					RF_element[1] += G.y;
					RF_element[2] += G.z;
				}
			}
		}
	}

	virtual __host__ void save_pulse_diagram() const{
		recorder record_p("pulse");
		ofstream trial_p = record_p.setup_record_csv();
		Vector3 timing(0,0,0);

		for (int i = 0; i < steps; i++){
			timing = getGCPU(Vector3(1, 1, 1), i*par->timestep);
			trial_p << i << "," << timing.x << ',' << timing.y << ',' << timing.z << std::endl;
		}
	}

	virtual __device__ __host__ int getSteps() const{
		return steps;
	}
	virtual __device__ __host__ int getReadFactor() const{
		return readFactor;
	}
	virtual __device__ __host__ int getReadSteps() const{
		return readSteps;
	}
	virtual __host__ int getNSubSequences(void) const{
		return n_sub_sequences;
	}
	virtual __host__ Sequence** devPointer() const{
		return dev_ptr;
	}
};

#endif
