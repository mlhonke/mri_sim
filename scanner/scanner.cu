#include <chrono>
#include "scanner.cuh"
#include "../primitives/lattice.cuh"

Scanner::Scanner(	
	Sequence& sequence,
	Coil& coil,
	SimuParams& params) 
	: 
	sequence(&sequence),
	coil(&coil),
	params(&params)
	{
	scan_k = new kSpace(params.res_x, params.res_y);
	lattice_present = false;
}

Scanner::Scanner(
	Sequence& sequence,
	Coil& coil,
	SimuParams& params,
	Lattice& lattice)
	:
	sequence(&sequence),
	coil(&coil),
	params(&params),
	lattice(&lattice)
	{
	scan_k = new kSpace(params.res_x, params.res_y);
	lattice_present = true;
}

bool Scanner::scan(){
	bool scan_success;

	if (lattice_present){
		scan_success = scan_lattice();
	} else {
		scan_success = scan_single_basis();
	}

	return scan_success;
}

bool Scanner::scan_lattice(){
	printf("Building basis pointer array.\n");
	lattice->setBasisSize(primitives.size());
	printf("Detected %d basis objects\n", lattice->getBasisSize());
	Primitive** basis_dev_pointers[lattice->getBasisSize()];
	for (int i = 0; i < primitives.size(); i++){
		basis_dev_pointers[i] = primitives[i]->devPointer();
	}

	cudaStream_t streams[sequence->getNSubSequences()];
	printf("Starting scan.\n");
	printf("Number of sub sequences = %d.\n", sequence->getNSubSequences());
	//Scan* scans = new Scan[sequence->getNSubSequences()];
	Scan* scans[sequence->getNSubSequences()];

	for (int i = 0; i< sequence->getNSubSequences(); i++){
		acqs.push_back(new magAcquisition(params, sequence->getSubSequences(i)));
		cudaStreamCreate(&streams[i]);
		scans[i] = new Scan(acqs[i], params, sequence->getSubSequences(i), lattice, basis_dev_pointers, coil->devPointer(), 0, std::vector<int>(15), streams[i]);
	}

	cudaDeviceSynchronize();

	long start = std::chrono::duration_cast< std::chrono::milliseconds >(std::chrono::system_clock::now().time_since_epoch()).count();
	for (int i = 0; i < sequence->getNSubSequences(); i++){
		printf("Running scan %d/%d\n",i,sequence->getNSubSequences()-1);
		scans[i]->run_scan_lattice();
	}

	safe_cuda(cudaDeviceSynchronize(), "Post Scan\n");
	long end = std::chrono::duration_cast< std::chrono::milliseconds >(std::chrono::system_clock::now().time_since_epoch()).count();
	printf("Simulation Kernel Time: %ld\n", end-start);
	for (int i = 0; i < sequence->getNSubSequences(); i++){
		scans[i]->saveScan();
		cudaDeviceSynchronize();
		make_k_space(acqs[i], sequence->getSubSequences(i));
	}

	cudaDeviceSynchronize();

	scan_k->get_fft();

	return true;
}

bool Scanner::scan_single_basis(){
	cudaStream_t streams[sequence->getNSubSequences()];
	printf("Starting scan.\n");
	printf("Number of sub sequences = %d.\n", sequence->getNSubSequences());
	//Scan* scans = new Scan[sequence->getNSubSequences()];
	Scan* scans[sequence->getNSubSequences()];

	for (int i = 0; i< sequence->getNSubSequences(); i++){
		acqs.push_back(new magAcquisition(params, sequence->getSubSequences(i)));
		cudaStreamCreate(&streams[i]);
		scans[i] = new Scan(acqs[i], params, sequence->getSubSequences(i), primitives[0], coil->devPointer(), 0, std::vector<int>(15), streams[i]);
	}

	cudaDeviceSynchronize();
	long start = std::chrono::duration_cast< std::chrono::milliseconds >(std::chrono::system_clock::now().time_since_epoch()).count();
	for (int i = 0; i < sequence->getNSubSequences(); i++){
		printf("Running scan %d/%d\n",i,sequence->getNSubSequences()-1);
		scans[i]->runScan();
	}

	safe_cuda(cudaDeviceSynchronize(), "Post Scan\n");
	long end = std::chrono::duration_cast< std::chrono::milliseconds >(std::chrono::system_clock::now().time_since_epoch()).count();
	printf("Simulation Kernel Time: %ld\n", end-start);

	for (int i = 0; i < sequence->getNSubSequences(); i++){
		scans[i]->saveScan();
		cudaDeviceSynchronize();
		make_k_space(acqs[i], sequence->getSubSequences(i));
	}

	cudaDeviceSynchronize();

	scan_k->get_fft();

	return true;
}

bool Scanner::scanCPU(){
	printf("Starting scan.\n");
	printf("Number of sub sequences = %d.\n", sequence->getNSubSequences());
	//Scan* scans = new Scan[sequence->getNSubSequences()];
	ScanCPU* scans[sequence->getNSubSequences()];

	for (int i = 0; i< sequence->getNSubSequences(); i++){
		acqs.push_back(new magAcquisition(params, sequence->getSubSequences(i)));
		scans[i] = new ScanCPU(acqs[i], params, sequence->getSubSequences(i), primitives[0], coil);
	}

	long start = std::chrono::duration_cast< std::chrono::milliseconds >(std::chrono::system_clock::now().time_since_epoch()).count();
	for (int i = 0; i < sequence->getNSubSequences(); i++){
		printf("Running scan %d/%d\n",i,sequence->getNSubSequences()-1);
		scans[i]->runScan();
	}

	long end = std::chrono::duration_cast< std::chrono::milliseconds >(std::chrono::system_clock::now().time_since_epoch()).count();
	printf("Simulation Kernel Time: %ld\n", end-start);

	for (int i = 0; i < sequence->getNSubSequences(); i++){
		scans[i]->saveScan();
		make_k_space(acqs[i], sequence->getSubSequences(i));
	}

	scan_k->get_fft();

	return true;
}

bool Scanner::add_primitive(Primitive& new_primitive){
	primitives.push_back(&new_primitive);
	
	return true;
}

bool Scanner::make_k_space(magAcquisition *acq, const Sequence *seq){
	int kx;
	int ky;

	for (int i = seq->get_k_start(); i < seq->get_k_end(); i++){
		kx = seq->getK(i).x;
		ky = seq->getK(i).y;

		scan_k->set_Mx(kx, ky, (1-2*((kx+ky)&1))*acq->get_signal_x()[i%seq->getReadSteps()]);
		scan_k->set_My(kx, ky, (1-2*((kx+ky)&1))*acq->get_signal_y()[i%seq->getReadSteps()]);
	}

	return true;
}
