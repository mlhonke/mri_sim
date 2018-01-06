#ifndef _SIMUPARAMS_H_
#define _SIMUPARAMS_H_

#include <time.h>
#include "../master_def.h"
#include "../util/vector3.cuh"

class SimuParams {

public:

	int number_of_particles;
	int particles_per_stream;
	bool particle_concurrency;
	int num_streams;
	int steps;
	int measurements;
	int n_mags_track;
	unsigned long long seed;
	int blocks;
	int res_x;
	int res_y;
	int seed_offset;
	real TR;
	real TE;
	real timestep;
	real mx_initial;
	real my_initial;
	real mz_initial;
	real time_end;
	real FOVx;
	real FOVy;
	Vector3 B0;

	SimuParams* devPointer;

public:
	__host__ SimuParams();

	__host__ SimuParams(int _number_of_particles,
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
		);

	__host__ int getSeed();

	__host__ void copyToDevice();

};

#endif
