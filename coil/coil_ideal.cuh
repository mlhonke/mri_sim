/*
Ideal coil with perfect homogeneity.
*/

#ifndef COIL_IDEAL_CU
#define COIL_IDEAL_CU

#include "coil.cuh"

class Coil_Ideal : public Coil{
public:
	Coil_Ideal** dev_ptr;

	__device__ __host__ Coil_Ideal();

	__device__ __host__ Vector3 getField(const Sequence* sequence, Vector3 r, real time) const;

	__device__ Vector3 getField(const Sequence* sequence, Vector3 r, real time, real* G_tensor, real* RF_tensor, int thread) const;

	__host__ const Coil** devPointer() const;
};

#endif
