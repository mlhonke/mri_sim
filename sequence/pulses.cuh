#ifndef PULSES_CUH
#define PULSES_CUH

#include "../master_def.h"
#include "../util/vector3.cuh"

class Pulse{
public:
	real start;
	real end;
	real duration;
	real strength;
	Vector3 bDir;
	Vector3 rDir;
	__device__ __host__ Pulse(real start, real duration, real strength, Vector3 bDir, Vector3 rDir);
	__device__ __host__ Pulse(real start, real duration, Vector3 bDir, Vector3 rDir);
	__device__ __host__ Pulse(real start, real duration);
	__device__ __host__ Pulse();
	__device__ __host__ virtual Vector3 out(real time, real seg_time, Vector3 r, int seg) const = 0;
	__device__ __host__ virtual Vector3 out(real time, real seg_time, int seg) const = 0;
	__device__ __host__ inline virtual bool on(real time) const{
		if (time >= start && time < end)
			return true;
		else
			return false;
	}
};


class PulseTrap : public Pulse{
public:
	real ramp;
	real delta;
	__device__ __host__ PulseTrap(real start, real duration, real strength, Vector3 bDir, Vector3 rDir, real ramp, real delta);
	__device__ __host__ PulseTrap();
	__device__ __host__ Vector3 out(real time, real seg_time, int seg) const;
	__device__ __host__ Vector3 out(real time, real seg_time, Vector3 r, int seg) const;
};

class RFflip : public Pulse{
public:
	real angle;
	real w;
	__device__ __host__ RFflip();
	__device__ __host__ RFflip(real start, real duration, real angle, Vector3 B0);
	__device__ __host__ Vector3 out(real time, real seg_time, int seg) const;
	__device__ __host__ Vector3 out(real time, real seg_time, Vector3 r, int seg) const;
};

__device__ __host__ real trap(real local_strength, real ramp, real start, real end, real local_time);

#endif
