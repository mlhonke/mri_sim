#include "magAcquisitionStreamCPU.cuh"

__host__ ScanCPU::ScanCPU(){}
__host__ ScanCPU::ScanCPU(
	magAcquisition* acq,
	SimuParams* par,
	const Sequence* host_seq,
	Primitive* basis,
	Coil* coil
	):
	acq(acq),
	par(par),
	host_seq(host_seq),
	basis(basis),
	coil(coil),
	measurements(par->measurements),
	number_of_particles(par->particles_per_stream)
{
	signal_x = new real[host_seq->getReadSteps()];
	signal_y = new real[host_seq->getReadSteps()];
	signal_z = new real[host_seq->getReadSteps()];
}

__host__ void ScanCPU::runScan(){
	updateWalkersMagCPU(
		par,
		basis,
		host_seq,
		coil,
		par->n_mags_track,
		signal_x,
		signal_y,
		signal_z
		);
}

__host__ void ScanCPU::saveScan(){
	for (int i = 0; i < host_seq->getReadSteps(); i++){
		acq->get_signal_x()[i] = (signal_x[i]);
		acq->get_signal_y()[i] = (signal_y[i]);
		acq->get_signal_z()[i] = (signal_z[i]);
	}
}
