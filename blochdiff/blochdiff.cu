#include "blochdiff.cuh"

/*
Evaluation of Bloch derivative without relaxation
*/
__device__ __host__ Vector3 bloch(const Vector3 & M, const Vector3 & B){
	return (M%B)*GAMMA;
}

/*
Evaluation of Bloch derivative with relaxation
*/
__device__ __host__ Vector3 bloch(const Vector3 & M, const Vector3 & B, const real T1, const real T2){
	return (M%B)*GAMMA - Vector3(M.x / T2, M.y / T2, (M.z - 1.0) / T1);
}

/*
Fourth Order Runge Kutta without Relaxation
*/
template <class FieldType>
__device__ __host__ void updateMagRK4(const Vector3 & r,
	real* Mx,
	real* My,
	real* Mz,
	const FieldType & B,
	real i,
	real h)
{
	real tn = i*h;
	Vector3 M(*Mx, *My, *Mz);
	Vector3 k1 = bloch(M, B(r, tn))*h;
	Vector3 k2 = bloch((M + k1*.5), B(r, tn + .5*h))*h;
	Vector3 k3 = bloch((M + k2*.5), B(r, tn + .5*h))*h;
	Vector3 k4 = bloch((M + k3), B(r, tn + h))*h;
	Vector3 finaldM = (k1 + k2*2.0 + k3*2.0 + k4)*(1.0 / 6.0);
	*Mx += finaldM.x;
	*My += finaldM.y;
	*Mz += finaldM.z;
}

/*
Fourth Order Runge Kutta with Relaxation
*/
__device__ void updateMagRK4(const Vector3 &r,
	real* Mx,
	real* My,
	real* Mz,
	const Coil **C,
	Sequence **B,
	Vector3 B0,
	real T1,
	real T2,
	real i,
	real h,
	real* G_tensor,
	real* RF_tensor,
	int thread)
{
	real tn = i*h;
	Vector3 M(*Mx, *My, *Mz);
	Vector3 k1 = bloch(M, (*C)->getField(*B, r, tn, G_tensor, RF_tensor, thread) + B0, T1, T2)*h;
	Vector3 k2 = bloch((M + k1*.5), (*C)->getField(*B, r, tn + .5*h, G_tensor, RF_tensor, thread) + B0, T1, T2)*h;
	Vector3 k3 = bloch((M + k2*.5), (*C)->getField(*B, r, tn + .5*h, G_tensor, RF_tensor, thread) + B0, T1, T2)*h;
	Vector3 k4 = bloch((M + k3), (*C)->getField(*B, r, tn + h, G_tensor, RF_tensor, thread) + B0, T1, T2)*h;
	Vector3 finaldM = (k1 + k2*2.0 + k3*2.0 + k4)*(1.0 / 6.0);
	*Mx += finaldM.x;
	*My += finaldM.y;
	*Mz += finaldM.z;
}

__device__ void updateMagRK4(const Vector3 &r,
	real* Mx,
	real* My,
	real* Mz,
	const Coil **C,
	Sequence **B,
	Vector3 B0,
	real T1,
	real T2,
	real i,
	real h)
{
	real tn = i*h;
	Vector3 M(*Mx, *My, *Mz);
	Vector3 k1 = bloch(M, ((*C)->getField(*B, r, tn) + B0), T1, T2)*h;
	Vector3 k2 = bloch((M + k1*.5), (*C)->getField(*B, r, tn + .5*h) + B0, T1, T2)*h;
	Vector3 k3 = bloch((M + k2*.5), (*C)->getField(*B, r, tn + .5*h) + B0, T1, T2)*h;
	Vector3 k4 = bloch((M + k3), (*C)->getField(*B, r, tn + h) + B0, T1, T2)*h;
	Vector3 finaldM = (k1 + k2*2.0 + k3*2.0 + k4)*(1.0 / 6.0);
	*Mx += finaldM.x;
	*My += finaldM.y;
	*Mz += finaldM.z;
}

//CPU Version, RK4 with relaxation.
__host__ void updateMagRK4(const Vector3 &r,
	real* Mx,
	real* My,
	real* Mz,
	const Coil *C,
	const Sequence *B,
	Vector3 B0,
	real T1,
	real T2,
	real i,
	real h)
{
	real tn = i*h;
	Vector3 M(*Mx, *My, *Mz);
	Vector3 k1 = bloch(M, ((C)->getField(B, r, tn) + B0), T1, T2)*h;
	Vector3 k2 = bloch((M + k1*.5), (C)->getField(B, r, tn + .5*h) + B0, T1, T2)*h;
	Vector3 k3 = bloch((M + k2*.5), (C)->getField(B, r, tn + .5*h) + B0, T1, T2)*h;
	Vector3 k4 = bloch((M + k3), (C)->getField(B, r, tn + h) + B0, T1, T2)*h;
	Vector3 finaldM = (k1 + k2*2.0 + k3*2.0 + k4)*(1.0 / 6.0);
	*Mx += finaldM.x;
	*My += finaldM.y;
	*Mz += finaldM.z;
}

__device__ __host__ void updateMag_noT2_rotate(real* Mx, real* My, real* Mz, const Vector3 & G, const Vector3 &B0, real timestep){
	Vector3 M(*Mx, *My, *Mz);
	Vector3 B = G + B0;
	real Bmag = B.magnitude();
	Vector3 C = (B*(1.0 / Bmag))*(M*(B*(1.0 / Bmag)));
	Vector3 u = M - C;
	Vector3 finalM = (C + u*cos(-GAMMA*Bmag*timestep) + ((B*(1.0 / Bmag)) % u)*sin(-GAMMA*Bmag*timestep));
	*Mx = finalM.x;
	*My = finalM.y;
	*Mz = finalM.z;
}

__device__ __host__ Vector3 updateMag_noT2_rotate(const Vector3 M, const Vector3 B, real timestep){
	real Bmag = B.magnitude();
	Vector3 C = (B*(1.0 / Bmag))*(M*(B*(1.0 / Bmag)));
	Vector3 u = M - C;
	return (C + u*cos(-GAMMA*Bmag*timestep) + ((B*(1.0 / Bmag)) % u)*sin(-GAMMA*Bmag*timestep));
}


__host__ Vector3 updateMag_noT2_rotate_CPU_fast(Vector3 M, Vector3 B, real timestep){
	real Bmag = B.magnitude();
	Vector3 C = (B*(1.0 / Bmag))*(M*(B*(1.0 / Bmag)));
	Vector3 u = M - C;
	return (C + u*cos(-GAMMA*Bmag*timestep) + ((B*(1.0 / Bmag)) % u)*sin(-GAMMA*Bmag*timestep));
}


__device__ __host__ Vector3 updateMag_noT2_Euler(Vector3 M, Vector3 B, real timestep){
	return (M%B)*timestep*GAMMA;
}
