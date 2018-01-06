/*
kernelPhase.cuh is the main computation kernel for the high-field approximation simulation.
Responsibilities:
- Initialize the particle's position (r)
- Initialize the particle's phase at time zero.
- Move the particle randomly through space.
- Update the phase of the particle by calling the sequence (G) as it moves through space and time.
*/

#ifndef _KERNELPHASE_H_
#define _KERNELPHASE_H_

#include "util/deviates.h"
#include "params/simuParams.cuh"

template <class GradientType, class Basis>
__global__ void updateWalkersPhase(
	const SimuParamsPhase *par,
	const Basis* basis,
	const GradientType * G,
	curandState* globalState,
	real* phase
	)

{

	const unsigned int tid = threadIdx.x + blockIdx.x*blockDim.x;

	curandState localState = globalState[tid];

	real phi, theta;
	real speed = sqrt(6.0*basis->getD() / par->timestep);
	//printf("speed: %f\n", speed);
	real stepLength = speed*par->timestep;

	//Initialize the particle position at time zero in the basis.
	Vector3 r = basis[0].unifRand(localState);
	//printf("r: %f %f %f\n", r.x, r.y, r.z);

	//Initialize the particle phase at time zero.
	for (int i = 0; i < par->measurements; i++){
		//phase[tid + i*par->number_of_particles] = par->phase_initial;
		phase[tid + i*par->number_of_particles] = GAMMA*(G[i](r, 0.0, phase[tid + i*par->number_of_particles]))*par->timestep;
	}

	//Run the simulation (step by step, moving the particle).
	for (int i = 1; i < par->steps; i++){

		Vector3 ri = r;

#if defined LATTICE_PICKING

		r += Vector3(stepLength*randomSign(localState), stepLength*randomSign(localState), stepLength*randomSign(localState));

#else  //Pick a point on a sphere.
		phi = 2.0*PI*curand_uniform(&localState);
		theta = acos(2.0*curand_uniform(&localState) - 1);

		//r is updated with the most recent movement.
		r += Vector3(speed*par->timestep*sin(theta)*cos(phi), speed*par->timestep*sin(theta)*sin(phi), speed*par->timestep*cos(theta));
#endif

#if defined SPECULAR_REFLECTION

		boundaryNormal(ri, r, speed, basis, par->timestep);

#else

		if (!basis[0].inside(r)){
			r = ri;
		}

#endif
		//Need to update the particle phase in the new position. Do this for each sequence in the simulation.
		for (int j = 0; j < par->measurements; j++){
			//printf("%f, %f\n", i*par->timestep, phase[tid + j*par->number_of_particles]);
#ifdef USE_PERFECT_RF
			G[j].pulse(&phase[tid + j*par->number_of_particles], i, par->timestep);
#endif

			phase[tid + j*par->number_of_particles] += GAMMA*(G[j](r, par->timestep*i, phase[tid + j*par->number_of_particles]))*par->timestep;
		}

	}

	//globalState[tid] = localState; 

}

#endif
