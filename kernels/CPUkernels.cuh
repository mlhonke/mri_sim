#ifndef _CPU_KERNELS_
#define _CPU_KERNELS_

//#include "util/deviates.h"
#include "../params/simuParams.cuh"
#include "../blochdiff/blochdiff.cuh"
#include "../kernels/boundaryCheck.cuh"
#include "../coil/coil.cuh"
#include "../primitives/primitive.cuh"

void updateWalkersMagCPU(SimuParams *par, Primitive* basis, const Sequence* B,
		Coil* coil, int n_mags_track, real *signal_x, real *signal_y,
		real *signal_z);

//template<class Basis>
//__host__ void boundaryNormalCPU(Vector3 & ri, Vector3 & r,
//		const real currentSpeed, const Basis* basis, const real timestep) {
//
//	real v = 0.0;
//	real accumtime = 0.0;
//
//	while (basis->intersect(ri, r, v)) {
//
//		r = (r - ri) * v + ri;
//
//		//determine incident vector
//		const Vector3 I = (r - ri);
//		const real Imag = I.magnitude();
//
//		//calculate accumulated time
//		accumtime += Imag / currentSpeed;
//
//		//determine normal for reflection calculation
//		const Vector3 n = basis->getNormal(r);
//
//		//set last position to boundary position before rflection
//		ri = r;
//
//		//reflect, reflection travel time is the remaining time (timestep - accumtime)
//		r += (I - n * 2.0 * (n * I)) * (currentSpeed / Imag)
//				* (timestep - accumtime);
//	}
//
//}

#endif
