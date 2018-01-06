#include "GRE.cuh"
#include "pulses.cuh"

__global__ void GRE_GPU(Sequence** obj_ptr, SimuParams* par, int phase_enc_offset = 0, int _local_res_x = 0);

__device__ __host__ GRE::GRE():Sequence(){}

__device__ __host__ GRE::GRE(SimuParams* par, int _phase_enc_offset, int _local_res_x) : Sequence(1,5, phase_enc_offset, _local_res_x, par){
	BW = 1/(10*par->timestep)*(1.0/2.0);
	Gf = 4.0*PI*BW*(1.0/par->FOVy)*(1.0/GAMMA);
	ramp_time = 0.005;
	phase_enc_offset = _phase_enc_offset;
	T_read_duration = par->res_y/(2*BW);
	T_phase_dur = 0.5;
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
	readFactor = 10;
	pulse[0] = new RFflip(0.0, 0.1, PI/4.0, par->B0);
	pulse[1] = new PulseTrap(pulse[0]->end, T_phase_dur, G_max, Vector3(0,0,1), Vector3(0,1,0), ramp_time, delta_G);
	pulse[2] = new PulseTrap(pulse[1]->end, T_read_duration/2.0, -Gf, Vector3(0,0,1), Vector3(1,0,0), ramp_time, 0);
	pulse[3] = new PulseTrap(pulse[2]->end, T_read_duration, Gf, Vector3(0,0,1), Vector3(1,0,0), ramp_time, 0);
	pulse[4] = new PulseTrap(pulse[3]->end, T_phase_dur, -G_max, Vector3(0,0,1), Vector3(0,1,0), ramp_time, -delta_G);
#ifndef __CUDA_ARCH__
	printf("BW: %f\n", BW);
	printf("read-factor float %f\n", (T_read_duration / (par->timestep*(par->res_y - 1))));
	printf("Max phase encoding step: %f\n", G_max);
	printf("Max freq encoding step: %f\n", Gf);
	printf("Read factor: %d\n", readFactor);
	safe_cuda(cudaMalloc(&dev_ptr, sizeof(GRE**)));
	GRE_GPU << <1, 1 >> >(dev_ptr, par->devPointer, _phase_enc_offset, _local_res_x);
#ifdef ALLOC_G
	make_tensors();
#endif
#endif
}

__device__ __host__ GRE::GRE(bool parallel, SimuParams* par) : Sequence(par->res_x,5, 0, 0, par), parallel(parallel){
	printf("Creating sub-sequences.\n");
	sub_seq = new GRE[par->res_x];
	for (int i = 0; i < par->res_x; i++){
		sub_seq[i] = GRE(par, i, 1);
	}
}

__device__ __host__ Vector3 GRE::getK(int readStep) const{
	int k_y = readStep % par->res_y;
	int k_x = readStep / par->res_y;

	return Vector3(k_x, k_y, 0);
}

__device__ __host__ int GRE::get_k_start() const{
	return phase_enc_offset*par->res_y;
}

__device__ __host__ int GRE::get_k_end() const{
	return (local_res_x+phase_enc_offset) * par->res_y;
}

__host__ const Sequence* GRE::getSubSequences(int i) const{
	if (parallel)
		return &sub_seq[i];
	else
		return this;
}

__device__ __host__ real GRE::getReadStart(real time) const{
	int seg = (int)(time / TR);
	return seg*TR + pulse[3]->start;
}

__device__ __host__ real GRE::getReadFinish(real time) const{
	int seg = (int)(time / TR);
	return seg*TR + pulse[3]->end;
}

__global__ void GRE_GPU(Sequence** obj_ptr, SimuParams* par, int _phase_enc_offset, int _local_res_x){
	*obj_ptr = new GRE(par, _phase_enc_offset, _local_res_x);
}
