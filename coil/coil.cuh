#ifndef COIL_CUH
#define COIL_CUH

#include "../util/vector3.cuh"
#include "../sequence/sequence.cuh"

class Sequence;

class Coil{
public:
	/*
	getG: Returns the gradient or B1 at a particular time and location in space. The sequence
	is mostly responsible for the gradient, but the coil is what produces the gradient.
	The gradient might be altered depending on the homogeneity of the coil.
	*/
	virtual __device__ __host__ Vector3 getField(const Sequence* sequence, Vector3 r, real time) const = 0;
	virtual __device__ Vector3 getField(const Sequence* sequence, Vector3 r, real time, real* G_tensor, real* RF_tensor, int thread) const = 0;
	virtual __host__ const Coil** devPointer() const = 0;
};

#endif
