/*
homoCompPrimitive: 
Added: August 8, 2016
Author: Michael Honke
*/

#ifndef _STMB_PRIMITIVE_H_
#define _STMB_PRIMITIVE_H_

#include "..\dSim_master.cuh"

template <class basis>
class stmbPrimitive {

private:

	std::vector<basis> basis_set;
	std::vector<basis> basis_set_spawn;
	int n_spawn;

public:

	__device__ __host__ stmbPrimitive(){}

	__device__ __host__ stmbPrimitive(std::vector<basis> & _basis_set, std::vector<bool> & _basis_spawn){
		basis_set = _basis_set;
		int i = 0;

		for (std::vector<basis>::iterator cur_basis = _basis_spawn.begin(); cur_basis != _basis_spawn.end(); ++cur_basis, i++){
			if (*cur_basis == true){
				basis_set_spawn.push_back(basis_set[i]);
			}
		}
	}

	__device__ __host__ ~stmbPrimitive(){};

	__host__ void add_primitive(basis & new_primitive, bool spawn){
		basis_set.push_back(new_primitive);
		if (spawn == true){ 
			basis_set_spawn.push_back(new_primitive);
		}
	}

	__device__ Vector3 unifRand(curandState localState) const{
		int cur_basis = curand(&localState) % basis_set_spawn.size();
		return basis_set_spawn[cur_basis].unifRand(localState);
	}

	__host__ Vector3 unifRandCPU() const{
		int cur_basis = rand() % basis_set_spawn.size();
		return basis_set_spawn[cur_basis].unifRandCPU();
	}

	__device__ __host__ bool inside(const Vector3 & r) const{

		return ((r.x - center.x)*(r.x - center.x) + (r.y - center.y)*(r.y - center.y) < radius*radius) && (abs(r.z - center.z) < length / 2);
	}

	__device__ __host__ bool inside_on(const Vector3 & r) const{
		const real rho = (r.x - center.x)*(r.x - center.x) + (r.y - center.y)*(r.y - center.y);
		const real cyc_z = abs(r.z - center.z);
		return ((rho < radius*radius || real_equal(rho, radius*radius, cylEPS)) && (cyc_z < length / 2 - EPSILON));
	}

	__device__ __host__ bool inside(real x, real y, real z) const{

		return ((x - center.x)*(x - center.x) + (y - center.y)*(y - center.y) < radius*radius) && (abs(z - center.z) < length / 2);
	}

	/*
	intersect: Determines if a particle has hit the side of the object.
	*/
	__device__ __host__ bool intersect(const Vector3 & ri, const Vector3 & rf, real & v, Vector3 & n) const{
		//Particle is inside the cylinder and might hit the left end or particle is to the left of the cylinder and might hit the left end.
		if (rf.z < left_end && inside_on(ri) == true || ri.z < left_end && rf.z > left_end){
			return intersect_end(ri, rf, v, n, Vector3(0, 0, -1), left_end);
		} //Particle is inside the cylinder and might hit the right end or particle is to the right of the cylinder and might hit the right end.
		else if (rf.z > right_end && inside_on(ri) == true || ri.z > right_end && rf.z < right_end){
			return intersect_end(ri, rf, v, n, Vector3(0, 0, 1), right_end);
		} //Particle cannot possibly hit either end, but check if it hits the side.
		else {
			return intersect_side(ri, rf, v, n);
		}
	}

	//here r is a point on the surface
	//Needs to be updated for length parameter
	__device__ __host__ Vector3 getNormalSide(Vector3 & r) const{
		double n_x = r.x - center.x;
		double n_y = r.y - center.y;
		double mag = sqrt(n_x*n_x + n_y*n_y);
		return Vector3(n_x / mag, n_y / mag, 0.0);
	}

	__device__ __host__ real getRadius() const{

		return radius;

	}

	__device__ __host__ int getRegion(const Vector3 & r) const{
		return region;
	}

	__device__ __host__  real getT2(const Vector3 & r) const{
		if (inside(r)){
			return T2;
		}
		return -1.0;
	}

	__device__ __host__ real getT2() const{
		return T2;
	}

	__device__ __host__ real getT1() const{
		return T1;
	}

	__device__ __host__ real getD(const Vector3 & r) const{
		if (inside(r)){
			return D;
		}
		return -1.0;
	}

	__device__ __host__ real getD() const{
		return D;
	}

	__device__ __host__ real getPermeability() const{
		return permeability;
	}

	__device__ __host__ Vector3 getCenter() const{
		return Vector3(center.x, center.y, center.z);
	}

	__host__ void setCenter(Vector3 v){
		center.x = v.x;
		center.y = v.y;
		center.z = v.z;
	}

	__host__ void setRadius(real _r){
		radius = _r;
	}

	__host__ void setEPS(real _cylEPS){
		cylEPS = _cylEPS;
	}

	__host__ void setRegion(int _region){
		region = _region;
	}

	__host__ int getRegion(){
		return region;
	}

	__device__ void randUnif(Vector3 & r, curandState & localState) const{
		do {
			r = Vector3((2.0*curand_uniform(&localState) - 1.0)*radius, (2.0*curand_uniform(&localState) - 1.0)*radius, (2.0*curand_uniform(&localState) - 1.0)*length / 2) + getCenter();
		} while (!inside(r));
	}

	__host__ void randUnif(Vector3 & r) const{
		do {
			r = Vector3((2.0*unifRandCPP() - 1.0)*radius, (2.0*unifRandCPP() - 1.0)*radius, length*(2.0*unifRandCPP() - 1) / 2.0) + getCenter();
		} while (!inside(r));
	}

};

#endif