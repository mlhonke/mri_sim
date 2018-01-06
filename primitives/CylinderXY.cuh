/*
	CylinderXY. Augments old versions with a cylinder model that has a finite length.
	Old CylinderXY changed to infCylinderXY to better reflect its model. Based on infCylinderXY
	but with length parameter.
	Added: July 18, 2016
	Author: Michael Honke (based on original CylinderXY, now infCylinderXY, by Trevor Vincent).
	*/

#ifndef _CYLINDER_XY_H_
#define _CYLINDER_XY_H_

#include "primitive.cuh"

class Cylinder_XY : public Primitive{

private:
	real radius;
	real length;
	real left_end;
	real right_end;

	Vector3 center;
	real D;
	real T2;
	real T1;
	real permeability;
	int region;
	real cylEPS;

public:
	Cylinder_XY** dev_ptr;
	__device__ __host__ Cylinder_XY(Vector3 _center, real _radius, real _length, real _T2, real _T1, real _D, int _region, real _permeability, int _num_particles, real _eps = EPSILON);
	//__device__ __host__ ~Cylinder_XY(){};
	__host__ Primitive** devPointer() const;
	//this will need to be changed once we go to full generality
	__device__ Vector3 unifRand(curandState localState) const;
	__host__ Vector3 unifRandCPU() const;
	__device__ __host__ bool inside(const Vector3 & r) const;
	__device__ __host__ bool inside_on(const Vector3 & r) const;
	__device__ __host__ bool inside_on_side(const Vector3 & r) const;
	__device__ __host__ bool inside_on_end(const Vector3 & r) const;
	__device__ __host__ bool inside(real x, real y, real z) const;
	/*
	intersect: Determines if a particle has hit the side of the object.
	*/
	__device__ __host__ bool intersect(const Vector3 & ri, const Vector3 & rf, real & v, Vector3 & n) const;
	__device__ __host__ bool intersect_end(const Vector3 & ri, const Vector3 & rf, real & v, Vector3 & n, const Vector3 & n_dir, real end_z) const;
	__device__ __host__ bool intersect_side(const Vector3 & ri, const Vector3 & rf, real & v, Vector3 & n) const;
	//here r is a point on the surface
	//Needs to be updated for length parameter
	__device__ __host__ Vector3 getNormalSide(const Vector3 & r) const;
	__device__ __host__ Vector3 getNormalEnd(const Vector3 & r) const;
	__device__ __host__ Vector3 getNormal(const Vector3 & r) const;
	__device__ __host__ real getRadius() const;
	__device__ __host__ int getRegion(const Vector3 & r) const;
	__device__ __host__  real getT2(const Vector3 & r) const;
	__device__ __host__ real getT2() const;
	__device__ __host__ real getT1() const;
	__device__ __host__ real getD(const Vector3 & r) const;
	__device__ __host__ real getD() const;
	__device__ __host__ real getPermeability() const;
	__device__ __host__ Vector3 getCenter() const;
	__host__ void setCenter(Vector3 v);
	__host__ void setRadius(real _r);
	__host__ void setEPS(real _cylEPS);
	__host__ void setRegion(int _region);
	__device__ __host__ int getRegion();
	__device__ void randUnif(Vector3 & r, curandState & localState) const;
	__host__ void randUnif(Vector3 & r) const;
};

#endif
