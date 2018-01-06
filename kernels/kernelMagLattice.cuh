#ifndef _kernelMagLattice_H_
#define _kernelMagLattice_H_

#include "../master_def.h"
#include "../params/simuParams.cuh"
#include "../blochdiff/blochdiff.cuh"
#include "../kernels/boundaryCheck.cuh"
#include "../coil/coil.cuh"
#include "../primitives/primitive.cuh"
#include "../primitives/lattice.cuh"
#include <cub/cub.cuh>

template <bool trackingM, bool trackingX>
__global__ void update_walkers_lattice_mag(
	SimuParams *par,
	Lattice **lat,
	Primitive*** basis,
	Sequence** B,
	const Coil** coil,
	curandState* globalState,
	int n_mags_track,
	real *Mx_track,
	real *My_track,
	real *Mz_track,
	real *signal_x,
	real *signal_y,
	real *signal_z
	)
{
	typedef cub::BlockReduce<real, SIM_THREADS> BlockReduce;
	__shared__ typename BlockReduce::TempStorage temp_storage;

	const unsigned int tid = threadIdx.x + blockIdx.x*blockDim.x;
	curandState localState = globalState[tid];
	real Mx, My, Mz;
	real phi, theta, speed;
	Vector3 r, r_unbounded;
	real T1, T2, T2avg;

#ifdef INITIALIZE_IN_REGION
(*lat)->initializeInRegion(basis,localState,r, INITIALIZE_IN_REGION);
#else
(*lat)->initializeUniformly(r,localState);
#endif

	r_unbounded = r;
	speed = sqrt(6.0*(*lat)->getD(basis, r) / par->timestep);

	T1 = (*lat)->getT2(basis,r);
	T2 = (*lat)->getT1(basis,r);

	for (int i = 0; i < par->measurements; i++){

		Mx = par->mx_initial;
		My = par->my_initial;
		Mz = par->mz_initial;

#if defined RK4_RELAXATION
		updateMagRK4(r,
			&Mx,
			&My,
			&Mz,
			coil,
			B,
			par->B0,
			T1,
			T2,
			0.,
			par->timestep);
#elif defined DORP_NORELAXATION

		updateMag_DORP(r, &Mx[tid + i*par->number_of_particles],
			&My[tid + i*par->number_of_particles],
			&Mz[tid + i*par->number_of_particles],
			B[i],
			0.,
			par->timestep);

#elif defined ROTATE_MULTISTEP

		updateMag_rotate_noT2_multistep(r,
			&Mx[tid + i*par->number_of_particles],
			&My[tid + i*par->number_of_particles],
			&Mz[tid + i*par->number_of_particles],
			B[i],
			0.,
			par->timestep
			);

#elif defined RK4_NORELAXATION

		updateMagRK4(r, &Mx[tid + i*par->number_of_particles],
			&My[tid + i*par->number_of_particles],
			&Mz[tid + i*par->number_of_particles],
			B[i],
			0,
			par->timestep);

#else

		updateMag_noT2_rotate(&Mx,
			&My,
			&Mz,
			(*coil)->getField(*B, r, (real)0.0),
			par->B0,
			par->timestep);
#endif

	}

	if (trackingM){
		if (tid < n_mags_track){
			Mx_track[tid*(*B)->getSteps()] = Mx;
			My_track[tid*(*B)->getSteps()] = My;
			Mz_track[tid*(*B)->getSteps()] = Mz;
		}
	}

	if (trackingX){
		if (tid < n_mags_track){
			Mx_track[tid*(*B)->getSteps()] = r.x;
			My_track[tid*(*B)->getSteps()] = r.y;
			Mz_track[tid*(*B)->getSteps()] = r.z;
		}
	}

	int s = 0;
	for (int i = 1; i < (*B)->getSteps(); i++){
		real time = i * par->timestep;

		Vector3 ri = r;
		phi = 2.0*PI*curand_uniform(&localState);
		theta = acos(2.0*curand_uniform(&localState) - 1);
		r += Vector3(speed*par->timestep*sin(theta)*cos(phi), speed*par->timestep*sin(theta)*sin(phi), speed*par->timestep*cos(theta));

#if defined SPECULAR_REFLECTION
#if defined USE_PERMEABLE
	real accumtime = 0.0;
    boundaryCheckPermeable(lat, basis, ri, r, accumtime, par->timestep, localState);
#else
	boundaryCheckImpermeable(lat, basis, ri,r,speed, par->timestep);
#endif

#else // SIMPLE REJECTION

    if ( (*lat)->inRegion(basis,ri) != (*lat)->inRegion(basis,r) )
      {
		r = ri;
      }

#endif

#if defined USE_RELAXATION
#if defined USE_PERMEABLE
	real T2_i = lat->getT2(basis,ri);
	real T2_f = lat->getT2(basis,r);
	T2 = ((accumtime*T2_i) + ((par->timestep - accumtime)*T2_f))/par->timestep
	real T1_i = lat->getT1(basis,ri);
	real T1_f = lat->getT1(basis,r);
	T1 = ((accumtime*T1_i) + ((par->timestep - accumtime)*T1_f))/par->timestep
#else
	T2 = (*lat)->getT2(basis, ri);
	T1 = (*lat)->getT1(basis, ri);
#endif
#endif

	(*lat)->correctBoundary(r);
		/*
		par->timestep : The interval of time that the simulation uses. This will remain consistent throughout the simulation.
		*/
		for (int j = 0; j < par->measurements; j++){

#if defined USE_PERFECT_RF
			B[j].pulse(
				//The magnetizations in each direction for the current particle.
				&Mx[tid + j*par->number_of_particles],
				&My[tid + j*par->number_of_particles],
				&Mz[tid + j*par->number_of_particles],
				//The current simulation time
				i,
				par->timestep
				);
#endif

#if defined RK4_RELAXATION
			updateMagRK4(r,
				&Mx,
				&My,
				&Mz,
				coil,
				B,
				par->B0,
				T1,
				T2,
				i,
				par->timestep);

#elif defined DORP_NORELAXATION

			updateMag_DORP(
				r,
				&Mx[tid + j*par->number_of_particles],
				&My[tid + j*par->number_of_particles],
				&Mz[tid + j*par->number_of_particles],
				B[j],
				i,
				par->timestep
				);

#elif defined ROTATE_MULTISTEP

			updateMag_rotate_noT2_multistep(
				r,
				&Mx[tid + j*par->number_of_particles],
				&My[tid + j*par->number_of_particles],
				&Mz[tid + j*par->number_of_particles],
				B[j],
				i,
				par->timestep
				);

#else
			updateMag_noT2_rotate(&Mx,
				&My,
				&Mz,
				(*coil)->getField(*B, r, time),
				par->B0,
				par->timestep);

			// updateMagRK4( r, &Mx[tid + j*par->number_of_particles] ,
			// &My[tid + j*par->number_of_particles] ,
			// &Mz[tid + j*par->number_of_particles] ,
			// B[j],
			// i,
			// par->timestep );


#endif
			//printf ( " Kernel Mafter = %.6f %.6f %.6f \n",  Mx[tid + j*par->number_of_particles], My[tid + j*par->number_of_particles], Mz[tid + j*par->number_of_particles] );

		}
		if (trackingM){
			if (tid < n_mags_track){
				Mx_track[tid*(*B)->getSteps() + i] = Mx;
				My_track[tid*(*B)->getSteps() + i] = My;
				Mz_track[tid*(*B)->getSteps() + i] = Mz;
			}
		}

		if (trackingX){
			if (tid < n_mags_track){
				Mx_track[tid*(*B)->getSteps() + i] = r.x;
				My_track[tid*(*B)->getSteps() + i] = r.y;
				Mz_track[tid*(*B)->getSteps() + i] = r.z;
			}
		}

		//If timestep is a readstep:
		if (time >= (*B)->getReadStart(time) &&
			time < (*B)->getReadFinish(time) &&
			( (int)(i - (*B)->getReadStart(time)/par->timestep)) % (*B)->getReadFactor() == 0){
			//Wait for all threads in this block to be done.
			__syncthreads();

			//if (tid == 0)
				//printf("Time: %f and %d and %d and finish %f\n", time, (i - (int) ((*B)->getReadStart(time) / par->timestep)), (int)((*B)->getReadStart(time) / par->timestep), (*B)->getReadFinish(time));

			double w = (par->B0).z * GAMMA;
			real signal_x_block = BlockReduce(temp_storage).Sum(Mx*cos(w*time) - My*sin(w*time));
			real signal_y_block = BlockReduce(temp_storage).Sum(Mx*sin(w*time) + My*cos(w*time));
			real signal_z_block = BlockReduce(temp_storage).Sum(Mz);

			//Save the total signal from this block to global memory for later summation.
			if (threadIdx.x == 0){
				signal_x[s * par->blocks + blockIdx.x] = signal_x_block;
				signal_y[s * par->blocks + blockIdx.x] = signal_y_block;
				signal_z[s * par->blocks + blockIdx.x] = signal_z_block;
				s++;
			}
		}

	}
	//globalState[tid] = localState;
}

#endif
