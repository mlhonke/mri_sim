#include "PGSE_D.cuh"

__global__ void PGSE_D_GPU(Sequence** obj_ptr, SimuParams* par, int phase_enc_offset = 0, int _local_res_x = 0);

__device__ __host__ PGSE_D::PGSE_D() : Sequence(1,0,0,0,0){}

__device__ __host__ PGSE_D::PGSE_D(SimuParams* par, int _phase_enc_offset, int _local_res_x) : Sequence(1, 12,_phase_enc_offset, _local_res_x, par){

	BW = 1/(10.0*par->timestep)*(1.0/2.0);
	Gf = 4.0*PI*BW*(1.0/par->FOVy)*(1.0/GAMMA);
	ramp_time = 0.005;
	phase_enc_offset = _phase_enc_offset;
	T_read_duration = par->res_y/(2*BW);
	T_phase_dur = T_read_duration/2.0;
	readSteps = local_res_x * par->res_y;
	phase_steps = par->res_x - 1.0;
	steps = (TR/par->timestep)*local_res_x;
	G_max = (1.0 / (par->FOVx)) * (PI / GAMMA) * (par->res_x-1.0)* (1.0 / (T_phase_dur)); //Max phase enc. gradient
	if (phase_steps == 0){
		delta_G = 0;
	}
	else {
		delta_G = (2.0 * G_max) / (par->res_x);
	}
	readFactor = (int) ((T_read_duration-2*ramp_time) / (par->timestep*(par->res_y - 1)));
	//readFactor = 10;

	Gd = 0.0;
	D_dur = 0.2;
	// Excitation RF
	pulse[0] = new RFflip(0, 0.1, PI / 1.85, par-> B0);
	// Phase encoding step
	pulse[1] = new PulseTrap(pulse[0]->end, T_phase_dur, G_max, Vector3(0,0,1), Vector3(0,1,0), ramp_time, delta_G);
	// Pre-read gradient
	pulse[2] = new PulseTrap(pulse[0]->end, pulse[1]->duration, Gf, Vector3(0,0,1), Vector3(1,0,0), ramp_time, 0);
	// Diffusion sensing gradient 1
	pulse[3] = new PulseTrap(pulse[2]->end, D_dur, Gd, Vector3(0,0,1), Vector3(0,0,1), ramp_time, 0);
	pulse[8] = new PulseTrap(pulse[2]->end, D_dur, 0.0, Vector3(0,0,1), Vector3(0,1,0), ramp_time, 0);
	pulse[9] = new PulseTrap(pulse[2]->end, D_dur, 0.0, Vector3(0,0,1), Vector3(1,0,0), ramp_time, 0);
	// 180 degree refocusing RF
	pulse[4] = new RFflip(pulse[3]->end, 0.1, 1.060*PI, par->B0);
	// Diffusion sensing gradient 2
	pulse[5] = new PulseTrap(pulse[4]->end, D_dur, Gd, Vector3(0,0,1), Vector3(0,0,1), ramp_time, 0);
	pulse[10] = new PulseTrap(pulse[4]->end, D_dur, 0.0, Vector3(0,0,1), Vector3(0,1,0), ramp_time, 0);
	pulse[11] = new PulseTrap(pulse[4]->end, D_dur, 0.0, Vector3(0,0,1), Vector3(1,0,0), ramp_time, 0);
	real iso_center_excite = (pulse[0]->start + pulse[0]->duration/2.0);
	real iso_center_inverse = (pulse[4]->start + pulse[4]->duration/2.0);
	real tau = (iso_center_inverse - iso_center_excite);
	// Read gradient
	pulse[6] = new PulseTrap(pulse[0]->duration/2.0+2.0*tau - T_read_duration/2.0, T_read_duration, Gf, Vector3(0,0,1), Vector3(1,0,0), ramp_time, 0);
	// Spoil gradient
	pulse[7] = new PulseTrap(pulse[6]->end, TR - pulse[6]->end, Gf, Vector3(0,0,1), Vector3(1,0,0), ramp_time, 0);

#ifndef __CUDA_ARCH__
	printf("Max phase encoding step: %f\n", G_max);
	printf("Max frequency encoding step: %f\n", Gf);
	safe_cuda(cudaMalloc(&dev_ptr, sizeof(PGSE_D**)));
	PGSE_D_GPU << <1, 1 >> >(dev_ptr, par->devPointer, _phase_enc_offset, _local_res_x);
	safe_cuda(cudaDeviceSynchronize());
#ifdef ALLOC_G
	make_tensors();
#endif
#endif
}

__device__ __host__ PGSE_D::PGSE_D(bool parallel, SimuParams* par) : Sequence(par->res_x, 12, 0, 0, par), parallel(parallel){
	printf("Creating sub-sequences.\n");
	sub_seq = new PGSE_D[par->res_x];
	for (int i = 0; i < par->res_x; i++){
		sub_seq[i] = PGSE_D(par, i, 1);
	}
}

__device__ __host__ Vector3 PGSE_D::getK(int readStep) const{
	int k_y = readStep % par->res_y;
	int k_x = readStep / par->res_y;
	return Vector3(k_x, k_y, 0);
}

__device__ __host__ int PGSE_D::get_k_start() const{
	return phase_enc_offset*par->res_y;
}
__device__ __host__ int PGSE_D::get_k_end() const{
	return (local_res_x+phase_enc_offset) * par->res_y;
}

__host__ const Sequence* PGSE_D::getSubSequences(int i) const{
	if (parallel)
		return &sub_seq[i];
	else
		return this;
}

__device__ __host__ real PGSE_D::getReadStart(real time) const{
	int seg = (int)(time / TR);
	return seg*TR + pulse[6]->start;
}

__device__ __host__ real PGSE_D::getReadFinish(real time) const{
	int seg = (int)(time / TR);
	return seg*TR + pulse[6]->end;
}

__global__ void PGSE_D_GPU(Sequence** obj_ptr, SimuParams* par, int _phase_enc_offset, int _local_res_x){
	*obj_ptr = new PGSE_D(par, _phase_enc_offset, _local_res_x);
}
