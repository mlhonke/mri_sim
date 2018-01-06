#ifndef _MAGACQUISITIONSTREAM_CUH_
#define _MAGACQUISITIONSTREAM_CUH_

#include <assert.h>
#include "../coil/coil.cuh"
#include "../sequence/sequence.cuh"
#include "../kernels/kernelMag.cuh"
#include "../kernels/kernelSetup.cuh"
#include "../util/pinnedVector.cuh"
#include "../util/cudaVector.cuh"
#include "../util/misc.cuh"
#include "../primitives/primitive.cuh"
#include "../primitives/lattice.cuh"
#include "magAcquisition.cuh"
#include "../params/simuParams.cuh"
#include "../util/recorder.h"

class Scan{
private:
	//Post-simulation storage
	cudaVector<real> signal_x;
	cudaVector<real> signal_y;
	cudaVector<real> signal_z;
	cudaVector<real> signal_x_total;
	cudaVector<real> signal_y_total;
	cudaVector<real> signal_z_total;
	cudaVector<real> dev_mx_track;
	cudaVector<real> dev_my_track;
	cudaVector<real> dev_mz_track;
	//Pre-simulation storage
	cudaVector<curandState> dev_states;
	cudaScalar<SimuParams> dev_par;
	//Simulation objects
	magAcquisition* acq;
	SimuParams* par;
	const Sequence* host_seq;
	const Lattice* lattice;
	const Primitive* basis;
	Primitive*** basis_dev_pointers;
	Primitive*** basis_dev_pointers_pointer;
	const Coil** coil;
	//Properties
	int devNum;
	const std::vector<int> *devSMP; //Number of streaming multiprocessors per device.
	cudaStream_t stream;
	int measurements;
	int number_of_particles;
	int num_blocks;

public:
	__host__ Scan();
	__host__ Scan(
		magAcquisition* acq,
		SimuParams* par,
		const Sequence* host_seq,
		const Primitive* basis,
		const Coil** coil,
		int devNum,
		const std::vector<int> devSMP,
		cudaStream_t stream
		);
	__host__ Scan(
		magAcquisition* acq,
		SimuParams* par,
		const Sequence* host_seq,
		Lattice* lattice,
		Primitive*** basis_dev_pointers,
		const Coil** coil,
		int devNum,
		const std::vector<int> devSMP,
		cudaStream_t stream
		);
	__host__ void runScan();
	__host__ void run_scan_lattice();
	__host__ void saveScan();
};
#endif
