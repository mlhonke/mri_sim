#ifndef _BLOCHDIFF_H_
#define _BLOCHDIFF_H_

#include "../util/vector3.cuh"
#include <cuda.h>
#include <cuda_runtime.h>
#include "../coil/coil.cuh"

/*
Evaluation of Bloch derivative without relaxation
*/
__device__ __host__ Vector3 bloch(const Vector3 & M, const Vector3 & B);
/*
Evaluation of Bloch derivative with relaxation
*/
__device__ __host__ Vector3 bloch(const Vector3 & M, const Vector3 & B, const real T1, const real T2);
/*
Fourth Order Runge Kutta without Relaxation
*/
template <class FieldType>
__device__ __host__ void updateMagRK4(const Vector3 & r, real* Mx, real* My, real* Mz, const FieldType & B, real i, real h);
/*
Fourth Order Runge Kutta withRelaxation
*/
__device__ void updateMagRK4(const Vector3 &r, real* Mx, real* My, real* Mz, const Coil **C, Sequence **B, Vector3 B0, real T1, real T2, real i, real h, real* G_tensor, real* RF_tensor, int thread);
__device__ void updateMagRK4(const Vector3 &r, real* Mx, real* My, real* Mz, const Coil **C, Sequence **B, Vector3 B0, real T1, real T2, real i, real h);
__host__ void updateMagRK4(const Vector3 &r, real* Mx, real* My, real* Mz, const Coil *C, const Sequence *B, Vector3 B0, real T1, real T2, real i, real h);
__device__ __host__ void updateMag_noT2_rotate(real* Mx, real* My, real* Mz, const Vector3 & G, const Vector3 &B0, real timestep);
__device__ __host__ Vector3 updateMag_noT2_rotate(const Vector3 M, const Vector3 B, real timestep);
__host__ Vector3 updateMag_noT2_rotate_CPU_fast(Vector3 M, Vector3 B, real timestep);
__device__ __host__ Vector3 updateMag_noT2_Euler(Vector3 M, Vector3 B, real timestep);

template <class FieldType>
__device__  void updateMagDORP(const Vector3 & M, const Vector3 & r, const FieldType & B, real tn, real & h, real hmax, real hmin, Vector3 & Mnew){
	Vector3 k1 = bloch(M, B(r, tn))*h;
	Vector3 k2 = bloch((M + k1*.2), B(r, tn + .2*h))*h;
	Vector3 k3 = bloch((M + k1*(3.0 / 40.0) + k2*(9.0 / 40.0)), B(r, tn + .3*h))*h;
	Vector3 k4 = bloch((M + k1*(44.0 / 45.0) + k2*(-56.0 / 15.0) + k3*(32.0 / 9.0)), B(r, tn + .8*h))*h;
	Vector3 k5 = bloch((M + k1*(19372.0 / 6561.0) + k2*(-25360.0 / 2187.0) + k3*(64448.0 / 6561.0) + k4*(-212.0 / 729)), B(r, tn + (8.0 / 9.0)*h))*h;
	Vector3 k6 = bloch((M + k1*(9017.0 / 3168.0) + k2*(-355.0 / 33.0) + k3*(-46732.0 / 5247.0) + k4*(49.0 / 176.0) + k5*(-5103.0 / 18656.0)), B(r, tn + h))*h;
	Vector3 k7 = bloch((M + k1*(35.0 / 384.0) + k3*(500.0 / 1113.0) + k4*(125.0 / 192.0) + k5*(-2187.0 / 6784.0) + k6*(11.0 / 84.0)), B(r, tn + h))*h;
	Mnew = M + k1*(35.0 / 384.0) + k3*(500.0 / 1113.0) + k4*(125.0 / 192.0) + k5*(-2187.0 / 6784.0) + k6*(11.0 / 84.0);
	Vector3 Mnew5 = M + k1*(5179.0 / 57600.0) + k3*(7571.0 / 16695.0) + k4*(393.0 / 640.0) + k5*(-92097.0 / 339200.0) + k6*(187.0 / 2100.0) + k7*(1.0 / 40.0);
	Vector3 diff = Mnew5 - Mnew;
	Vector3 ratio;
	ratio.x = diff.x / (.000001 + .000001*max(Mnew5.x, Mnew.x));
	ratio.y = diff.y / (.000001 + .000001*max(Mnew5.y, Mnew.y));
	ratio.z = diff.z / (.000001 + .000001*max(Mnew5.z, Mnew.z));
	real err = sqrt(1.0 / 3.0*ratio.squareMagnitude());
	real hopt = h*(.98)*pow((1.0 / err), .2);
	if (hopt < hmax){ h = hmin; }
}

template <class FieldType>
__device__ void updateMag_DORP(const Vector3 & r, real* Mx, real* My, real* Mz, const FieldType & B, int i, real h){

	real accumTime = 0.0;
	real tk = i*h;
	real dt = h;
	real hmin = B.getStepTimeMin();

	Vector3 M(*Mx, *My, *Mz);
	Vector3 Mnew;

	do {

		updateMagDORP(M, r, B, tk, h, dt, hmin, Mnew);
		if (h < (dt - accumTime)){
			accumTime += h;
			updateMagDORP(M, r, B, tk, h, dt, hmin, Mnew);
		}

		else{
			accumTime += h;
		}

		tk = i*dt + accumTime;
		h = dt - accumTime;
		M = Mnew;

	} while (accumTime < dt);

	*Mx = M.x;
	*My = M.y;
	*Mz = M.z;

}

template <class FieldType>
__device__ void updateMag_rotate_noT2_multistep(const Vector3 & r, real* Mx, real* My, real* Mz, const FieldType & B, int i, real timestep){

	real accumTime = 0.0;
	real tk = i*timestep;
	real hmin = B.getStepTimeMin(tk);

	Vector3 M(*Mx, *My, *Mz);
	Vector3 Mnew;

	do {
		Mnew = updateMag_noT2_rotate(M, B(r, tk), hmin);
		accumTime += hmin;
		tk += hmin;
		M = Mnew;
	} while (accumTime < timestep);

	*Mx = M.x;
	*My = M.y;
	*Mz = M.z;

}

#endif
