/*
CylinderXY. Replaces old versions with a cylinder model that has a finite length.
Old CylinderXY changed to infCylinderXY to better reflect its model. Based on infCylinderXY
but with length parameter.
Added: July 18, 2016
Author: Michael Honke (based on original CylinderXY, now infCylinderXY, by Trevor Vincent).
*/

#ifndef _PRIMITIVE_CUH_
#define _PRIMITIVE_CUH_

#include <cuda_runtime.h>
#include <curand_kernel.h>
#include "../util/vector3.cuh"
#include "../util/misc.cuh"
#include "../util/compare.cuh"

class Primitive {
	/*
protected:
	Vector3 center;
	real D;
	real T2;
	real T1;
	real permeability;
	int region;
	real cylEPS;
	*/
protected:
	int num_particles;

public:
	
	virtual __host__ Primitive** devPointer() const = 0;
	virtual __device__ Vector3 unifRand(curandState localState) const = 0;
	virtual __host__ Vector3 unifRandCPU() const = 0;
	virtual __device__ __host__ bool inside(const Vector3 & r) const = 0;
	virtual __device__ __host__ bool inside(real x, real y, real z) const = 0;
	virtual __device__ __host__ bool intersect(const Vector3 & ri, const Vector3 & rf, real & v, Vector3 & n) const = 0;
	virtual __device__ __host__ Vector3 getNormal(const Vector3 &r) const = 0;
	virtual __device__ void randUnif(Vector3 & r, curandState & localState) const = 0;
	virtual __host__ void randUnif(Vector3 & r) const = 0;
	virtual __device__ __host__ Vector3 getCenter() const = 0;
	virtual __host__ void setCenter(Vector3 v) = 0;
	virtual __device__ __host__ int getRegion(const Vector3 & r) const = 0;
	virtual __device__ __host__  real getT2(const Vector3 & r) const = 0;
	virtual __device__ __host__ real getT2() const = 0;
	virtual __device__ __host__ real getT1() const = 0;
	virtual __device__ __host__ real getD(const Vector3 & r) const = 0;
	virtual __device__ __host__ real getD() const = 0;
	virtual __device__ __host__ real getPermeability() const = 0;
	virtual __host__ void setEPS(real _cylEPS) = 0;
	virtual __host__ void setRegion(int _region) = 0;
	virtual __device__ __host__ int getRegion() = 0;
	virtual __device__ __host__ int getNumParticles() const{
		return num_particles;
	}
	
};

#endif
