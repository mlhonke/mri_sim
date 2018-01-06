#ifndef _LATTICE_H_
#define _LATTICE_H_

#include "primitive.cuh"

class Lattice;

__global__ void Lattice_GPU(Lattice** obj_ptr, real _a, real _b, real _c, real _T2, real _T1, real _D, int _basisSize);

class Lattice {

private:

  real a, b, c;
  real T2_lattice;
  real T1_lattice;
  real D_lattice;
  int basisSize;

public:
  Lattice** dev_ptr;

  __host__ __device__ Lattice(){}
  __host__ __device__ Lattice(real _a, real _b, real _c, real _T2, real _T1, real _D, int _basisSize){
#ifndef __CUDA_ARCH__
	cudaMalloc(&dev_ptr, sizeof(Lattice**));
	Lattice_GPU<<<1,1>>>(dev_ptr, _a, _b, _c, _T2, _T1, _D, _basisSize);
#endif
	T2_lattice = _T2;
	T1_lattice = _T1;
	D_lattice = _D;
	a=_a; b=_b; c=_c;
	basisSize = _basisSize;
  }

  __host__ Lattice** devPointer() const{
	  return dev_ptr;
  }

  __host__ __device__ real getA(){
	return a;
  }
  
  __host__ __device__ real getB(){
	return b;
  }
  
  __host__ __device__ real getC(){
	return c;
  }
  
  __device__ void initializeUniformly(Vector3 & r, curandState & localState) const
  {
    
	r.x = curand_uniform(&localState)*a;
	r.y = curand_uniform(&localState)*b;
	r.z = curand_uniform(&localState)*c;

  }
  
    __host__ void initializeUniformlyCPU(Vector3 & r) const
  {
    
    r.x = unifRandCPP()*a;
    r.y = unifRandCPP()*b;
    r.z = unifRandCPP()*c;

  }
  
	__device__ void initializeInRegion( Primitive *** basis, curandState & localState, Vector3 & r, int region) const;
	__device__ __host__ void correctBoundary(Vector3 & r) const;
	__device__ __host__ void lmod(Vector3 & r) const;
	__device__ __host__ double getT2(Primitive *** basis, Vector3 & r) const;
	__device__ __host__ double getT1(Primitive *** basis, Vector3 & r) const;

  __device__ __host__ double getD(Primitive *** basis, const Vector3 & r) const{

    double D;
	Vector3 temp = r;
	lmod(temp);
    for (int i = 0; i < basisSize; i++){
      if( D = (*basis[i])->getD(), D > 0 && (*basis[i])->inside(temp)){
	return D;
      }
    }

    return D_lattice;
  }
  
__device__ __host__ int inRegion(Primitive *** basis, const Vector3 & r) const;
  
  __device__ __host__ bool inLatticeCell (const Vector3 & r) const{
//    if ( r.x > a || r.x < 0.0 ||
//			  r.y > b || r.y < 0.0 ||
//					   r.z > c || r.z < 0.0 )
//      {
//	return false;
//      }
    return true;
  }

  __device__ __host__ bool intersectionCheckPermeable (Primitive*** basis,
					     const Vector3 & ri,
					     const Vector3 & rf,
					     real & v,
						 Vector3 & n,
						 real & permeability)const

  {

    //displacement vector

    Vector3 dr = rf - ri;
    
    //modded line 1

    Vector3 mod1_f = rf; lmod(mod1_f);
    Vector3 mod1_i = mod1_f - dr;
    
    //modded line 2
    Vector3 mod2_i = ri; lmod(mod2_i); 
    Vector3 mod2_f = mod2_i + dr;

    real vTemp = 0.0;
    real vBestEst = 10.0;

    for (int i = 0; i < basisSize; i++){

      if((*basis[i])->intersect(mod1_i, mod1_f, vTemp, n) && vTemp < vBestEst){
		vBestEst = vTemp;
		Vector3 intPoint = dr*vBestEst + mod1_i;
		n = (*basis[i])->getNormal( intPoint );
		permeability = (*basis[i])->getPermeability();
      }
      
      if((*basis[i])->intersect(mod2_i, mod2_f, vTemp, n) && vTemp < vBestEst){
		vBestEst = vTemp;
		Vector3 intPoint = dr*vBestEst + mod2_i;
		n = (*basis[i])->getNormal( intPoint );
		permeability = (*basis[i])->getPermeability();
      }

    } 

	v = vBestEst;
    return (vBestEst < 5.0);

  }
  
  __device__ __host__ bool intersectionCheckImpermeable(Primitive ***basis,
					     const Vector3 & ri,
					     const Vector3 & rf,
					     real & v,
						 Vector3 & n)const

  {

    //displacement vector

    Vector3 dr = rf - ri;
    
    //modded line 1

    Vector3 mod1_f = rf; lmod(mod1_f);
    Vector3 mod1_i = mod1_f - dr;
    
    //modded line 2
    Vector3 mod2_i = ri; lmod(mod2_i); 
    Vector3 mod2_f = mod2_i + dr;

    real vTemp = 0.0;
    real vBestEst = 10.0;

    for (int i = 0; i < basisSize; i++){

      if((*basis[i])->intersect(mod1_i, mod1_f, vTemp, n) && vTemp < vBestEst){
		vBestEst = vTemp;
		Vector3 intPoint = dr*vBestEst + mod1_i;
		n = (*basis[i])->getNormal( intPoint );
      }
      
      if((*basis[i])->intersect(mod2_i, mod2_f, vTemp, n) && vTemp < vBestEst){
		vBestEst = vTemp;
		Vector3 intPoint = dr*vBestEst + mod2_i;
		n = (*basis[i])->getNormal( intPoint );
      }

    } 

	v = vBestEst;
    return (vBestEst < 5.0);

  }
  
  __device__ __host__ int getBasisSize(){
  
	return basisSize;
  
  }
  

  __device__ __host__ void setBasisSize(int _basisSize){
	basisSize = _basisSize;
  }

};

#endif
