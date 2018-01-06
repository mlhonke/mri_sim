#include "magAcquisitionStream.cuh"
#include "../kernels/kernelMagLattice.cuh"

__host__ Scan::Scan(){}
__host__ Scan::Scan(
	magAcquisition* acq,
	SimuParams* par,
	const Sequence* host_seq,
	const Primitive* basis,
	const Coil** coil,
	int devNum,
	const std::vector<int> devSMP,
	cudaStream_t stream
	):
	acq(acq),
	par(par),
	host_seq(host_seq),
	basis(basis),
	coil(coil),
	devNum(devNum),
	devSMP(&devSMP),
	stream(stream),
	measurements(par->measurements),
	number_of_particles(par->particles_per_stream)
{
	num_blocks = par->blocks;
	signal_x.malloc(host_seq->getReadSteps()*num_blocks, stream);
	signal_y.malloc(host_seq->getReadSteps()*num_blocks, stream);
	signal_z.malloc(host_seq->getReadSteps()*num_blocks, stream);
	signal_x_total.malloc(host_seq->getReadSteps(), stream);
	signal_y_total.malloc(host_seq->getReadSteps(), stream);
	signal_z_total.malloc(host_seq->getReadSteps(), stream);
	dev_states.malloc(number_of_particles, stream);
	dev_par.malloc(stream);
	dev_par = *par;
	if (par->n_mags_track > 0){
		dev_mx_track.malloc(measurements*par->steps*par->n_mags_track, stream);
		dev_my_track.malloc(measurements*par->steps*par->n_mags_track, stream);
		dev_mz_track.malloc(measurements*par->steps*par->n_mags_track, stream);
	}
	acq->set_tracked_particles(par->n_mags_track);
	safe_cuda(cudaGetLastError(), "Malloc");
	cudaStreamSynchronize(stream);
	dev_par.copyToDevice();
	safe_cuda(cudaGetLastError(), "Malloc2");
	setup_kernel <<< num_blocks, SIM_THREADS, 0, stream>>> (dev_states.getPointer(), par->getSeed());
	safe_cuda(cudaGetLastError(), "Setup");
}

__host__ Scan::Scan(
	magAcquisition* acq,
	SimuParams* par,
	const Sequence* host_seq,
	Lattice* lattice,
	Primitive*** basis_dev_pointers,
	const Coil** coil,
	int devNum,
	const std::vector<int> devSMP,
	cudaStream_t stream
	):
	acq(acq),
	par(par),
	host_seq(host_seq),
	lattice(lattice),
	basis_dev_pointers(basis_dev_pointers),
	coil(coil),
	devNum(devNum),
	devSMP(&devSMP),
	stream(stream),
	measurements(par->measurements),
	number_of_particles(par->particles_per_stream)
{
	size_t size_basis_dev_pointers = lattice->getBasisSize() * sizeof(Primitive**);
	cudaMalloc((void****)&basis_dev_pointers_pointer, size_basis_dev_pointers);
	cudaMemcpy(basis_dev_pointers_pointer, basis_dev_pointers, size_basis_dev_pointers, cudaMemcpyHostToDevice);
	num_blocks = par->blocks;
	signal_x.malloc(host_seq->getReadSteps()*num_blocks, stream);
	signal_y.malloc(host_seq->getReadSteps()*num_blocks, stream);
	signal_z.malloc(host_seq->getReadSteps()*num_blocks, stream);
	signal_x_total.malloc(host_seq->getReadSteps(), stream);
	signal_y_total.malloc(host_seq->getReadSteps(), stream);
	signal_z_total.malloc(host_seq->getReadSteps(), stream);
	dev_states.malloc(number_of_particles, stream);
	dev_par.malloc(stream);
	dev_par = *par;
	if (par->n_mags_track > 0){
		dev_mx_track.malloc(measurements*par->steps*par->n_mags_track, stream);
		dev_my_track.malloc(measurements*par->steps*par->n_mags_track, stream);
		dev_mz_track.malloc(measurements*par->steps*par->n_mags_track, stream);
	}
	acq->set_tracked_particles(par->n_mags_track);
	safe_cuda(cudaGetLastError(), "Malloc");
	cudaStreamSynchronize(stream);
	dev_par.copyToDevice();
	safe_cuda(cudaGetLastError(), "Malloc2");
	setup_kernel <<< num_blocks, SIM_THREADS, 0, stream>>> (dev_states.getPointer(), par->getSeed());
	safe_cuda(cudaGetLastError(), "Setup");
}

__host__ void Scan::runScan(){
	if (par->n_mags_track > 0){
		//printf("Pre walkers kernel\n");
		updateWalkersMag<true, false> << < num_blocks, SIM_THREADS, 0, stream >> > (
			dev_par.getPointer(),
			basis->devPointer(),
			host_seq->devPointer(),
			coil,
			dev_states.getPointer(),
			par->n_mags_track,
			dev_mx_track.getPointer(),
			dev_my_track.getPointer(),
			dev_mz_track.getPointer(),
			signal_x.getPointer(),
			signal_y.getPointer(),
			signal_z.getPointer()
			);
		//safe_cuda(cudaGetLastError(),"Walkers Track");
	}
	else {
		//printf("Non-tracking Walker's Kernel\n");
		updateWalkersMag<false, false> << < num_blocks, SIM_THREADS, 0, stream >> > (
			dev_par.getPointer(),
			basis->devPointer(),
			host_seq->devPointer(),
			coil,
			dev_states.getPointer(),
			par->n_mags_track,
			0,
			0,
			0,
			signal_x.getPointer(),
			signal_y.getPointer(),
			signal_z.getPointer()
			);
		//safe_cuda(cudaGetLastError(), "Walkers");
	}
}

__host__ void Scan::run_scan_lattice(){
	if (par->n_mags_track > 0){
		//printf("Pre walkers kernel\n");
		update_walkers_lattice_mag<true, false> << < num_blocks, SIM_THREADS, 0, stream >> > (
			dev_par.getPointer(),
			lattice->devPointer(),
			basis_dev_pointers_pointer,
			host_seq->devPointer(),
			coil,
			dev_states.getPointer(),
			par->n_mags_track,
			dev_mx_track.getPointer(),
			dev_my_track.getPointer(),
			dev_mz_track.getPointer(),
			signal_x.getPointer(),
			signal_y.getPointer(),
			signal_z.getPointer()
			);
		//safe_cuda(cudaGetLastError(),"Walkers Track");
	}
	else {
		//printf("Non-tracking Walker's Kernel\n");
		update_walkers_lattice_mag<false, false> << < num_blocks, SIM_THREADS, 0, stream >> > (
			dev_par.getPointer(),
			lattice->devPointer(),
			basis_dev_pointers_pointer,
			host_seq->devPointer(),
			coil,
			dev_states.getPointer(),
			par->n_mags_track,
			0,
			0,
			0,
			signal_x.getPointer(),
			signal_y.getPointer(),
			signal_z.getPointer()
			);
		//safe_cuda(cudaGetLastError(), "Walkers");
	}
}

__host__ void Scan::saveScan(){
	int threads_sum = 512;
	signal_x.sum(signal_x_total, threads_sum, NUM_SM, host_seq->getReadSteps(), num_blocks, stream);
	signal_y.sum(signal_y_total, threads_sum, NUM_SM, host_seq->getReadSteps(), num_blocks, stream);
	signal_z.sum(signal_z_total, threads_sum, NUM_SM, host_seq->getReadSteps(), num_blocks, stream);
	safe_cuda(cudaDeviceSynchronize());
	signal_x_total.copyFromDevice();
	signal_y_total.copyFromDevice();
	signal_z_total.copyFromDevice();
	safe_cuda(cudaDeviceSynchronize());
	signal_x_total.copyTo(acq->get_signal_x());
	signal_y_total.copyTo(acq->get_signal_y());
	signal_z_total.copyTo(acq->get_signal_z());
	if (par->n_mags_track > 0){
		dev_mx_track.copyFromDevice();
		dev_my_track.copyFromDevice();
		dev_mz_track.copyFromDevice();
		cudaStreamSynchronize(stream);
		dev_mx_track.copyTo(acq->getMxTracked());
		dev_my_track.copyTo(acq->getMyTracked());
		dev_mz_track.copyTo(acq->getMzTracked());
	}
}
