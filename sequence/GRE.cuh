#ifndef GRE_CU
#define GRE_CU

#include "sequence.cuh"
#include "pulses.cuh"

class GRE : public Sequence{
private:
	GRE* sub_seq;
	real G_max;
	real Gf;
	bool parallel = false;
	real readFinish;
	real readStart;
	real B0;
	real B1;
	real T_RF;
	int phase_steps;
	real T_read_duration;
	real ramp_time;
	real T_phase_dur;
	real TE;
	real delta_G;
	real BW;

public:
	__device__ __host__ GRE();
	__device__ __host__ GRE(SimuParams* par, int _phase_enc_offset = 0, int _local_res_x = 0);
	__device__ __host__ GRE(bool parallel, SimuParams* par);
	__device__ __host__ Vector3 getK(int readStep) const;
	__device__ __host__ int get_k_start() const;
	__device__ __host__ int get_k_end() const;
	__host__ const Sequence* getSubSequences(int i) const;
	__device__ __host__ real getReadStart(real time) const;
	__device__ __host__ real getReadFinish(real time) const;
};

#endif
