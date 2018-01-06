/*
	CylinderXY. Augments old versions with a cylinder model that has a finite length.
	Old CylinderXY changed to infCylinderXY to better reflect its model. Based on infCylinderXY
	but with length parameter.
	Added: July 18, 2016
	Author: Michael Honke (based on original CylinderXY, now infCylinderXY, by Trevor Vincent).
	*/
#include "CylinderXY.cuh"

__global__ void Cylinder_XY_GPU(Cylinder_XY** obj_ptr, Vector3 _center, real _radius, real _length, real _T2, real _T1, real _D, int _region, real _permeability, int num_particles, real _eps);

__device__ __host__ Cylinder_XY::Cylinder_XY(Vector3 _center, real _radius, real _length, real _T2, real _T1, real _D, int _region, real _permeability, int _num_particles, real _eps){
#ifndef __CUDA_ARCH__
	cudaMalloc(&dev_ptr, sizeof(Cylinder_XY**));
	Cylinder_XY_GPU << <1, 1 >> >(dev_ptr, _center, _radius, _length, _T2, _T1, _D, _region, _permeability, num_particles, _eps);
#endif
	num_particles = _num_particles;
	center = _center;
	radius = _radius;
	length = _length;
	D = _D;
	T2 = _T2;
	T1 = _T1;
	region = _region;
	cylEPS = _eps;
	right_end = center.z + length / 2;
	left_end = center.z - length / 2;
}

//__device__ __host__ ~Cylinder_XY(){};

__host__ Primitive** Cylinder_XY::devPointer() const{
	return (Primitive**) dev_ptr;
}

//this will need to be changed once we go to full generality
__device__ Vector3 Cylinder_XY::unifRand(curandState localState) const{

	Vector3 r;
	do{
		r.x = radius*(2.0*curand_uniform(&localState) - 1) + center.x;
		r.y = radius*(2.0*curand_uniform(&localState) - 1) + center.y;
		r.z = length*(2.0*curand_uniform(&localState) - 1) / 2.0 + center.z;
	} while (!inside(r));

	return r;
}

__host__ Vector3 Cylinder_XY::unifRandCPU() const{

	Vector3 r;
	do{
		r.x = radius*(2.0*unifRandCPP() - 1) + center.x;
		r.y = radius*(2.0*unifRandCPP() - 1) + center.y;
		r.z = length*(2.0*unifRandCPP() - 1) / 2.0 + center.z;
	} while (!inside(r));

	return r;
}


__device__ __host__ bool Cylinder_XY::inside(const Vector3 & r) const{

	return ((r.x - center.x)*(r.x - center.x) + (r.y - center.y)*(r.y - center.y) < radius*radius) && (abs(r.z - center.z) < length / 2);
}

__device__ __host__ bool Cylinder_XY::inside_on(const Vector3 & r) const{
	const real rho = (r.x - center.x)*(r.x - center.x) + (r.y - center.y)*(r.y - center.y);
	const real cyc_z = abs(r.z - center.z);
	return ((rho < radius*radius || real_equal(rho, radius*radius, cylEPS)) && (cyc_z < length / 2 - EPSILON));
}

__device__ __host__ bool Cylinder_XY::inside_on_side(const Vector3 & r) const{
	const real rho = (r.x - center.x)*(r.x - center.x) + (r.y - center.y)*(r.y - center.y);
	const real cyc_z = abs(r.z - center.z);
	return (real_equal(rho, radius*radius, cylEPS) && (cyc_z < length / 2 - EPSILON));
}

__device__ __host__ bool Cylinder_XY::inside_on_end(const Vector3 & r) const{
	const real rho = (r.x - center.x)*(r.x - center.x) + (r.y - center.y)*(r.y - center.y);
	const real cyc_z = abs(r.z - center.z);
	return ((rho < radius*radius) && (real_equal(cyc_z, length/2, cylEPS)));
}

__device__ __host__ bool Cylinder_XY::inside(real x, real y, real z) const{

	return ((x - center.x)*(x - center.x) + (y - center.y)*(y - center.y) < radius*radius) && (abs(z - center.z) < length / 2);
}

/*
	intersect: Determines if a particle has hit the side of the object.
	*/
__device__ __host__ bool Cylinder_XY::intersect(const Vector3 & ri, const Vector3 & rf, real & v, Vector3 & n) const{
	//Particle is inside the cylinder and might hit the left end or particle is to the left of the cylinder and might hit the left end.
	if ((rf.z < left_end && inside_on(ri) == true) || (ri.z < left_end && rf.z > left_end)){
		return intersect_end(ri, rf, v, n, Vector3(0, 0, -1), left_end);
	} //Particle is inside the cylinder and might hit the right end or particle is to the right of the cylinder and might hit the right end.
	else if ((rf.z > right_end && inside_on(ri) == true) || (ri.z > right_end && rf.z < right_end)){
		return intersect_end(ri, rf, v, n, Vector3(0, 0, 1), right_end);
	} //Particle cannot possibly hit either end, but check if it hits the side.
	else {
		return intersect_side(ri, rf, v, n);
	}
}

__device__ __host__ bool Cylinder_XY::intersect_end(const Vector3 & ri, const Vector3 & rf, real & v, Vector3 & n, const Vector3 & n_dir, real end_z) const{
	Vector3 dr = rf - ri;
	//printf("End: %f\n", end_z);
	//printf("rf: %f %f %f\nri:%f %f %f\n", rf.x, rf.y, rf.z, ri.x, ri.y, ri.z);
	real line_param = (end_z - ri.z) / dr.z;
	//printf("Line_param = %f\n", line_param);
	real y_intersect = ri.y + dr.y * line_param;
	real x_intersect = ri.x + dr.x * line_param;
	//printf("y, x: %f %f\n", y_intersect, x_intersect);
	real rho = sqrt(pow(x_intersect - center.x, 2) + pow(y_intersect - center.y, 2));
	//printf("rho: %f\n", rho);

	//Particle is within the area of the cylinder end.
	if (rho <= radius){
		v = abs(end_z - ri.z) / abs(rf.z - ri.z);
		if (real_equal(v, 0.0, cylEPS)){//Then particle started on the wall and should move a full distance.
			return intersect_side(ri, rf, v, n);
		}
		Vector3 temp = (rf - ri)*v + ri;
		//printf("rfc: %f %f %f\n", temp.x, temp.y, temp.z);
		//printf("v: %f\n", v);
		n = n_dir;
		return true;
	} //Particle is outside the area of the cylinder end, but could possibly hit a cylinder side.
	else {
		return intersect_side(ri, rf, v, n);
	}
}

__device__ __host__ bool Cylinder_XY::intersect_side(const Vector3 & ri, const Vector3 & rf, real & v, Vector3 & n) const{

	Vector3 dr = rf - ri;
	real step_mag = dr.magnitude();

	// real a = dr.x*dr.x + dr.y*dr.y;
	// real b = 2*ri.x*dr.x + 2*ri.y*dr.y;
	// real c = ri.x*ri.x + ri.y*ri.y - radius*radius;

	real a = dr.x*dr.x + dr.y*dr.y;
	real b = 2.0*ri.x*dr.x - 2.0*dr.x*center.x + 2.0*ri.y*dr.y - 2.0*dr.y*center.y;
	real c = ri.x*ri.x + ri.y*ri.y - 2.0*ri.x*center.x - 2.0*ri.y*center.y + center.x*center.x + center.y*center.y - radius*radius;

	real q = -.5*(b + sgn(b)*sqrt(b*b - 4 * a*c));
	real root1 = q / a;
	real root2 = c / q;

	bool s1 = (root1 > 0.0 && root1 < 1.0 && b*b>4 * a*c && !real_equal(root1*step_mag, 0.0, cylEPS));
	bool s2 = (root2 > 0.0 && root2 < 1.0 && b*b>4 * a*c && !real_equal(root2*step_mag, 0.0, cylEPS));
	bool s3 = (fabs(root1) < fabs(root2));

	if ((s1 && s2 && s3) || (s1 && !s2)){
		v = root1;
		n = getNormalSide((rf - ri)*v + ri);
		return true;
	}

	else if ((s1 && s2 && !s3) || (s2 && !s1)){
		v = root2;
		n = getNormalSide((rf - ri)*v + ri);
		return true;
	}

	else {
		return false;
	}

}

//here r is a point on the surface
//Needs to be updated for length parameter
__device__ __host__ Vector3 Cylinder_XY::getNormal(const Vector3 & r) const{
	if (inside_on_side(r))
		return getNormalSide(r);
	else
		return getNormalEnd(r);
}

__device__ __host__ Vector3 Cylinder_XY::getNormalSide(const Vector3 & r) const{
	double n_x = r.x - center.x;
	double n_y = r.y - center.y;
	double mag = sqrt(n_x*n_x + n_y*n_y);
	return Vector3(n_x / mag, n_y / mag, 0.0);
}

__device__ __host__ Vector3 Cylinder_XY::getNormalEnd(const Vector3 & r) const{
	if (r.z < center.z)
		return Vector3(0.0,0.0,-1.0);
	else
		return Vector3(0.0,0.0,1.0);
}

__device__ __host__ real Cylinder_XY::getRadius() const{

	return radius;

}

__device__ __host__ int Cylinder_XY::getRegion(const Vector3 & r) const{
	return region;
}

__device__ __host__  real Cylinder_XY::getT2(const Vector3 & r) const{
	if (inside(r)){
		return T2;
	}
	return -1.0;
}

__device__ __host__ real Cylinder_XY::getT2() const{
	return T2;
}

__device__ __host__ real Cylinder_XY::getT1() const{
	return T1;
}

__device__ __host__ real Cylinder_XY::getD(const Vector3 & r) const{
	if (inside(r)){
		return D;
	}
	return -1.0;
}

__device__ __host__ real Cylinder_XY::getD() const{
	return D;
}

__device__ __host__ real Cylinder_XY::getPermeability() const{
	return permeability;
}

__device__ __host__ Vector3 Cylinder_XY::getCenter() const{
	return Vector3(center.x, center.y, center.z);
}

__host__ void Cylinder_XY::setCenter(Vector3 v){
	center.x = v.x;
	center.y = v.y;
	center.z = v.z;
}

__host__ void Cylinder_XY::setRadius(real _r){
	radius = _r;
}

__host__ void Cylinder_XY::setEPS(real _cylEPS){
	cylEPS = _cylEPS;
}

__host__ void Cylinder_XY::setRegion(int _region){
	region = _region;
}

__device__ __host__ int Cylinder_XY::getRegion(){
	return region;
}

__device__ void Cylinder_XY::randUnif(Vector3 & r, curandState & localState) const{
	do {
		r = Vector3((2.0*curand_uniform(&localState) - 1.0)*radius, (2.0*curand_uniform(&localState) - 1.0)*radius, (2.0*curand_uniform(&localState) - 1.0)*length / 2) + getCenter();
	} while (!inside(r));
}

__host__ void Cylinder_XY::randUnif(Vector3 & r) const{
	do {
		r = Vector3((2.0*unifRandCPP() - 1.0)*radius, (2.0*unifRandCPP() - 1.0)*radius, length*(2.0*unifRandCPP() - 1) / 2.0) + getCenter();
	} while (!inside(r));
}

__global__ void Cylinder_XY_GPU(Cylinder_XY** obj_ptr, Vector3 _center, real _radius, real _length, real _T2, real _T1, real _D, int _region, real _permeability, int num_particles, real _eps){
	if (threadIdx.x == 0 && blockIdx.x == 0)
		*obj_ptr = new Cylinder_XY(_center, _radius, _length, _T2, _T1, _D, _region, _permeability, num_particles, _eps);
}
