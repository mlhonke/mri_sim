#include "lattice.cuh"
#include "primitive.cuh"

__global__ void Lattice_GPU(Lattice** obj_ptr, real _a, real _b, real _c, real _T2, real _T1, real _D, int _basisSize){
	if (threadIdx.x == 0 && blockIdx.x == 0)
		*obj_ptr = new Lattice(_a, _b, _c, _T2, _T1, _D, _basisSize);
}

__device__ void Lattice::initializeInRegion( Primitive *** basis,
			     curandState & localState,
			     Vector3 & r,
				 int region
				 ) const{
	do {
		r.x = (curand_uniform(&localState)*2.0 - 1.0)*a;
		r.y = (curand_uniform(&localState)*2.0 - 1.0)*b;
		r.z = (curand_uniform(&localState)*2.0 - 1.0)*c;
	} while ( inRegion(basis,r) != region);
}

__device__ __host__ void Lattice::correctBoundary(Vector3 & r) const{
	if (!inLatticeCell(r)){
		lmod(r);
	}
}

__device__ __host__ void Lattice::lmod(Vector3 & r) const {
//    if( r.x > a) { r.x = fmod(r.x,a); }
//    if( r.x < 0.0) { r.x = fmod(r.x,a) + a; }
//
//    if( r.y > b) { r.y = fmod(r.y,b); }
//    if( r.y < 0.0) { r.y = fmod(r.y,b) + b; }
//
//    if( r.z > c) { r.z = fmod(r.z,c); }
//    if( r.z < 0.0) { r.z = fmod(r.z,c) + c; }
}

__device__ __host__ double Lattice::getT2(Primitive *** basis, Vector3 & r) const{
	Vector3 temp = r;
	lmod(temp);
  for (int i = 0; i < basisSize; i++){
		if((*basis[i])->inside(temp)){
			return (*basis[i])->getT2();
		}
  }
  return T2_lattice;
}

__device__ __host__ double Lattice::getT1(Primitive *** basis, Vector3 & r) const{
	Vector3 temp = r;
	lmod(temp);
  for (int i = 0; i < basisSize; i++){
		if((*basis[i])->inside(temp)){
			return (*basis[i])->getT1();
		}
  }
  return T1_lattice;
}

__device__ __host__ int Lattice::inRegion(Primitive *** basis, const Vector3 & r) const{
	int region;
	Vector3 temp = r;
	lmod(temp);
	for (int i = 0; i < basisSize; i++){
		region = (*basis[i])->getRegion(temp);
		if ((*basis[i])->inside(temp)){
			return region;
		}
	}
	return 0;
}
