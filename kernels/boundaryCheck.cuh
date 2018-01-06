#ifndef _BOUNDARY_CHECK_H_
#define _BOUNDARY_CHECK_H_

#include <curand.h>
#include <curand_kernel.h>

class Lattice;
class Basis;
class Vector3;
class Primitive;

__device__ void boundaryNormal(Vector3 & ri, Vector3 & r, const real currentSpeed,  const Primitive** basis, const real timestep);
__device__ void boundaryCheckImpermeable(const Lattice ** lat, Primitive *** basis, const Vector3 & ri, Vector3 & r, const real currentSpeed,  const real timestep);
__device__ void boundaryCheckPermeable(const Lattice ** lat, Primitive *** basis, const Vector3 & ri, Vector3 & r, real accumtime, const real timestep, curandState & localState);

#endif
