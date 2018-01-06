#ifndef _MAGACQUISITIONSTREAMCPU_H_
#define _MAGACQUISITIONSTREAMCPU_H_

#include <assert.h>
#include "../coil/coil.cuh"
#include "../sequence/sequence.cuh"
#include "../util/misc.cuh"
#include "../primitives/primitive.cuh"
#include "../primitives/lattice.cuh"
#include "magAcquisition.cuh"
#include "../params/simuParams.cuh"
#include "../kernels/CPUkernels.cuh"
#include "../util/recorder.h"

class ScanCPU{
private:
	//Post-simulation storage
	real* signal_x;
	real* signal_y;
	real* signal_z;
	//Simulation objects
	magAcquisition* acq;
	SimuParams* par;
	const Sequence* host_seq;
	Lattice* lattice;
	Primitive* basis;
	Coil* coil;
	//Properties
	int measurements;
	int number_of_particles;

public:
	__host__ ScanCPU();
	__host__ ScanCPU(
		magAcquisition* acq,
		SimuParams* par,
		const Sequence* host_seq,
		Primitive* basis,
		Coil* coil
		);
	__host__ void runScan();
	__host__ void saveScan();
};
#endif
