#ifndef _COMPARE_H_
#define _COMPARE_H_

__host__ __device__ inline bool doub_equal(double a, double b) {
	if (fabs(a - b) < EPSILON) {return true;}
	else {return false;}
}

__host__ __device__ inline bool real_equal(real a, real b, real eps) {
	if (fabs(a - b) < eps) {return true;}
	else {return false;}
}

#endif
