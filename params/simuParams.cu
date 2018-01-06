#include "simuParams.cuh"

__host__ SimuParams::SimuParams(){}

__host__ SimuParams::SimuParams(int _number_of_particles,
	int _particles_per_stream,
	real _TR,
	real _TE,
	real _timestep,
	int _n_mags_track,
	Vector3 m_initial,
	Vector3 _B0,
	int _res_x,
	int _res_y,
	real _FOVx,
	real _FOVy
	):seed_offset(0){

	number_of_particles = _number_of_particles;
	particles_per_stream = _particles_per_stream;
	if (number_of_particles != particles_per_stream)
		particle_concurrency = true;
	else
		particle_concurrency = false;
	num_streams = number_of_particles / particles_per_stream;
	TR = _TR;
	TE = _TE;
	measurements = 1;
	timestep = _timestep;
	n_mags_track = _n_mags_track;
	seed = time(NULL);
	mx_initial = m_initial.x;
	my_initial = m_initial.y;
	mz_initial = m_initial.z;
	time_end = timestep*steps;
	B0 = _B0;
	blocks = particles_per_stream / SIM_THREADS;
	res_x = _res_x;
	res_y = _res_y;
	FOVx = _FOVx;
	FOVy = _FOVy;

	//END PARAMETERS.
	steps = (TR / timestep) * res_x;
	printf("Number of steps in simulation: %d\n", steps);
	copyToDevice();
}

__host__ int SimuParams::getSeed(){
	return seed+(particles_per_stream*seed_offset++);
}

__host__ void SimuParams::copyToDevice(){
	cudaMalloc(&devPointer, sizeof(SimuParams));
	cudaMemcpy(devPointer, this, sizeof(SimuParams), cudaMemcpyHostToDevice);
}
