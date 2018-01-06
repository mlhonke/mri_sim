#ifndef PGSE_CU
#define PGSE_CU

#include "sequence.cuh"
#include "pulses.cuh"

class PGSE : public Sequence{
private:
	SimuParams* par;
	PGSE* sub_seq;
	real G_max;
	real Gf;
	bool parallel = false;
	real readFinish;
	real readStart;
	real TR;
	real B0;
	real B1;
	int phase_steps;
	real T_phase_start;
	real T_phase_end;
	real T_read_start;
	real T_read_end;
	real T_read_duration;
	real ramp_time;
	real T_phase_dur;
	real TE;
	real delta_G;
	int phase_enc_offset;
	int local_res_x;
	RFflip* excite;
	RFflip* inverse;
	PulseTrap* phase;
	PulseTrap* preread;
	PulseTrap* read;
	PulseTrap* spoil;

public:
	__device__ __host__ PGSE();
	__device__ __host__ PGSE(SimuParams* par, int _phase_enc_offset = 0, int _local_res_x = 0);
	__device__ __host__ PGSE(SimuParams* par, bool parallel);
	__device__ __host__ Vector3 getG(Vector3 r, real time) const;
	__device__ __host__ Vector3 getK(int readStep) const;
	__device__ __host__ int get_k_start() const;
	__device__ __host__ int get_k_end() const;
	__host__ const Sequence* getSubSequences(int i) const;
	__device__ __host__ real getReadStart(real time) const;
	__device__ __host__ real getReadFinish(real time) const;
};

#endif
