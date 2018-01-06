#include "PGSE.cuh"

__global__ void PGSE_GPU(Sequence** obj_ptr, SimuParams* par, int phase_enc_offset = 0, int _local_res_x = 0);

__device__ __host__ PGSE::PGSE() : Sequence(1){}

__device__ __host__ PGSE::PGSE(SimuParams* par, int _phase_enc_offset, int _local_res_x) : Sequence(1), par(par){
	_local_res_x == 0 ? local_res_x = par->res_x : local_res_x = _local_res_x;//Check if a subsequence or not.
	Gf = 0.000588;
	B0 = (par->B0).z;
	TR = par->TR;
	TE = par->TE;
	ramp_time = 0.05;
	phase_enc_offset = _phase_enc_offset;
	T_read_duration = (par->res_y)*(2.0 * PI * (1.0 / par->FOVy) * (1.0 / (GAMMA * Gf)));
	T_phase_dur = (T_read_duration)/2.0;
	readSteps = local_res_x * par->res_y;
	phase_steps = par->res_x - 1;
	steps = (TR/par->timestep)*local_res_x;
	G_max = (1.0 / par->FOVx) * (PI / GAMMA) * (par->res_x - 1)* (1.0 / (T_phase_dur - ramp_time)); //Max phase enc. gradient
	if (phase_steps == 0){
		delta_G = 0;
	}
	else {
		delta_G = (2.0 * G_max) / (par->res_x - 1);
	}
	readFactor = (int) ((T_read_duration) / (par->timestep*(par->res_y - 1)));

	excite = new RFflip(0, 0.01, PI / 2.0, par-> B0);
	phase = new PulseTrap(excite->end, T_phase_dur, G_max, Vector3(0,0,1), Vector3(0,1,0), ramp_time, delta_G);
	preread = new PulseTrap(excite->end, phase->duration, Gf, Vector3(0,0,1), Vector3(1,0,0), ramp_time, 0);
	inverse = new RFflip(preread->end, 0.01, 1.0*PI, par->B0);
	real iso_center_excite = (excite->start + excite->duration/2.0);
	real iso_center_inverse = (inverse->start + inverse->duration/2.0);
	real tau = (iso_center_inverse - iso_center_excite);
	read = new PulseTrap(excite->duration/2.0+2.0*tau - T_read_duration/2.0, T_read_duration, Gf, Vector3(0,0,1), Vector3(1,0,0), ramp_time, 0); //Read pulse
	spoil = new PulseTrap(read->end, TR - read->end, 0, Vector3(0,0,1), Vector3(1,0,0), ramp_time, 0);
#ifndef __CUDA_ARCH__
	safe_cuda(cudaMalloc(&dev_ptr, sizeof(PGSE**)));
	PGSE_GPU << <1, 1 >> >(dev_ptr, par->devPointer, _phase_enc_offset, _local_res_x);
#endif
}

__device__ __host__ PGSE::PGSE(SimuParams* par, bool parallel) : parallel(parallel){
	printf("Creating sub-sequences.\n");
	if (par->particle_concurrency){
		n_sub_sequences = par->num_streams;
		sub_seq = new PGSE[par->num_streams];
		for (int i = 0; i < par->num_streams; i++){
			sub_seq[i] = PGSE(par);
		}
	} else {
		n_sub_sequences = par->res_x;
		sub_seq = new PGSE[par->res_x];
		for (int i = 0; i < par->res_x; i++){
			sub_seq[i] = PGSE(par, i, 1);
		}
	}
}

__device__ __host__ Vector3 PGSE::getG(Vector3 r, real time) const{
	int seg = (int) (time / TR) + phase_enc_offset;
	real seg_time = time - seg*TR + TR*phase_enc_offset;
	if (seg_time <= excite->end){
		return excite->out(time, seg_time, r, seg);
	} else if (phase->start < seg_time && seg_time <= phase->end){
		return phase->out(time, seg_time, r, seg) + preread->out(time, seg_time, r, seg);
	} else if (inverse->start < seg_time && seg_time <= inverse->end){
		return inverse->out(time, seg_time, r, seg);
	} else if (read->start < seg_time && seg_time <= read->end){
		return read->out(time, seg_time, r, seg);
	} else if (spoil->start < seg_time && seg_time <= spoil->end) {
		return spoil->out(time, seg_time, r, seg);
	} else {
		return Vector3(0,0,0);
	}
}

__device__ __host__ Vector3 PGSE::getK(int readStep) const{
	int k_y = readStep % par->res_y;
	int k_x = readStep / par->res_y;
	return Vector3(k_x, k_y, 0);
}

__device__ __host__ int PGSE::get_k_start() const{
	return phase_enc_offset*par->res_y;
}
__device__ __host__ int PGSE::get_k_end() const{
	return (local_res_x+phase_enc_offset) * par->res_y;
}

__host__ const Sequence* PGSE::getSubSequences(int i) const{
	if (parallel)
		return &sub_seq[i];
	else
		return this;
}

__device__ __host__ real PGSE::getReadStart(real time) const{
	int seg = (int)(time / TR);
	return seg*TR + read->start;
}

__device__ __host__ real PGSE::getReadFinish(real time) const{
	int seg = (int)(time / TR);
	return seg*TR + read->end;
}

__global__ void PGSE_GPU(Sequence** obj_ptr, SimuParams* par, int _phase_enc_offset, int _local_res_x){
	*obj_ptr = new PGSE(par, _phase_enc_offset, _local_res_x);
}
