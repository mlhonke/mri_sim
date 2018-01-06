#include "pulses.cuh"

__device__ __host__ Pulse::Pulse(real start, real duration, real strength, Vector3 bDir, Vector3 rDir)
: start(start), duration(duration), strength(strength), bDir(bDir), rDir(rDir){
	end = start + duration;
#ifndef __CUDA_ARCH__
	printf("Pulse start: %f\nPulse end: %f\n================\n", start, end);
#endif
}

__device__ __host__ Pulse::Pulse(real start, real duration, Vector3 bDir, Vector3 rDir)
: Pulse(start, duration, 0, bDir, rDir){
}

__device__ __host__ Pulse::Pulse(real start, real duration)
: Pulse(start, duration, 0, Vector3(0,0,0), Vector3(0,0,0)){
}

__device__ __host__ Pulse::Pulse(){}

__device__ __host__ PulseTrap::PulseTrap(real start, real duration, real strength, Vector3 bDir, Vector3 rDir, real ramp, real delta)
:Pulse(start, duration, strength, bDir, rDir), ramp(ramp), delta(delta){}

__device__ __host__ PulseTrap::PulseTrap(){}

__device__ __host__ Vector3 PulseTrap::out(real time, real seg_time, Vector3 r, int seg) const{
	real local_strength = strength - delta * seg;
	real pulse_strength;

	real local_time = seg_time - start;
	if (local_time <= ramp)
		pulse_strength = (1 / ramp)*local_time * local_strength;
	else if (local_time > ramp && local_time <= (duration - ramp))
		pulse_strength = local_strength;
	else
		pulse_strength = (1 - (1 / ramp)*(local_time - duration + ramp))*local_strength;

	return Vector3(pulse_strength*bDir.x, pulse_strength*bDir.y, pulse_strength*bDir.z)*(r*rDir);
}

__device__ __host__ Vector3 PulseTrap::out(real time, real seg_time, int seg) const{
	real local_strength = strength - delta * seg;
	real pulse_strength;

	real local_time = seg_time - start;
	if (local_time <= ramp)
		pulse_strength = (1 / ramp)*local_time * local_strength;
	else if (local_time > ramp && local_time <= (duration - ramp))
		pulse_strength = local_strength;
	else
		pulse_strength = (1 - (1 / ramp)*(local_time - duration + ramp))*local_strength;

	return Vector3(pulse_strength*bDir.x, pulse_strength*bDir.y, pulse_strength*bDir.z);
}

__device__ __host__ RFflip::RFflip(real start, real duration, real angle, Vector3 B0)
:Pulse(start, duration), angle(angle){
	strength = angle / fabs(GAMMA*duration);
	w = B0.magnitude() * GAMMA;
}

__device__ __host__ RFflip::RFflip(){}

__device__ __host__ Vector3 RFflip::out(real time, real seg_time, Vector3 r, int seg) const{
	//time = seg_time - start;
	return Vector3(strength*cos(w*time), -strength*sin(w*time), 0);
}

__device__ __host__ Vector3 RFflip::out(real time, real seg_time, int seg) const{
	//time = seg_time - start;
	return Vector3(strength*cos(w*time), -strength*sin(w*time), 0);
}

__device__ __host__ real trap(real local_strength, real ramp, real start, real end, real local_time){
	if (local_time - start <= ramp)
		return (1 / ramp)*(local_time-start) * local_strength;
	else if (local_time - start > ramp && local_time <= (end - ramp))
		return local_strength;
	else
		return (1 - (1 / ramp)*(local_time - end + ramp))*local_strength;
}
