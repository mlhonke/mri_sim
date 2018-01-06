#ifndef _kernelMag_H_
#define _kernelMag_H_

#include "../params/simuParams.cuh"
#include "../blochdiff/blochdiff.cuh"
#include "../kernels/boundaryCheck.cuh"
#include "../coil/coil.cuh"
#include "../primitives/primitive.cuh"
#include <cub/cub.cuh>

/*__global__ void updateWalkersMag(
	const SimuParams *par,
	Primitive** basis,
	Sequence** B,
	Coil** coil,
	curandState* globalState,
	real* Mx,
	real* My,
	real* Mz
	){
	updateWalkersMag<false,false>(par, basis, B, coil, globalState, Mx, My, Mz, 0, 0, 0, 0);
	}*/

template <bool trackingM, bool trackingX>
__global__ void updateWalkersMag(
	SimuParams *par,
	Primitive** basis,
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

	__shared__ real G_tensor[G_SHARED_SIZE*9 + 9];
	__shared__ real RF_tensor[G_SHARED_SIZE*3 + 3];
//	printf("Grad value %f\n", (*B)->getG(Vector3(1,1,1), 200*par->timestep, G_tensor, RF_tensor, threadIdx.x).magnitude());
//	printf("Grad value old %f\n", (*B)->getG(Vector3(1,1,1), 200*par->timestep).magnitude());
//
//	printf("Grad value %f\n", (*B)->getG(Vector3(1,1,1), 350*par->timestep, G_tensor, RF_tensor, threadIdx.x).magnitude());
//	printf("Grad value old %f\n", (*B)->getG(Vector3(1,1,1), 350*par->timestep).magnitude());

	real Mx;
	real My;
	real Mz;

	const unsigned int tid = threadIdx.x + blockIdx.x*blockDim.x;
	//printf("tid: %d \n", tid);
	curandState localState = globalState[tid];

	real phi, theta;
	real speed = sqrt(6.0*(*basis)->getD() / par->timestep);

	// printf ( " timestep = %.6f \n" , par->timestep); 
	// printf ( " measurements = %d \n" , par->measurements); 
	// printf ( " number_of_particles = %d \n" , par->number_of_particles); 
	// printf ( " steps = %d \n" , par->steps); 
	// printf ( " c = %.6f \n" , sub->getC() ); 

	Vector3 r;

#ifdef MASS_DEBUG
	if (tid < 64){
	real spacingx = (par->FOVx)/8;
	real spacingy = (par->FOVy)/8;
	r = Vector3(fmod(tid*spacingx, par->FOVx), ((int)((tid*spacingx)/ par->FOVy)) * spacingy, 0.0);
	printf("Init pos: %f %f\n", r.x, r.y);
	} else {
		r = (*basis)->unifRand(localState);
	}
#else
	r = (*basis)->unifRand(localState);
#endif

	real T1 = (*basis)->getT1();
	real T2 = (*basis)->getT2();

	for (int i = 0; i < par->measurements; i++){

		Mx = par->mx_initial;
		My = par->my_initial;
		Mz = par->mz_initial;

		//printf ( "0, %E, %E, %E \n",  Mx[tid + i*par->number_of_particles], My[tid + i*par->number_of_particles], Mz[tid + i*par->number_of_particles] );

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
			par->timestep,
			G_tensor,
			RF_tensor,
			threadIdx.x);
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
		//if (i*par->timestep > 5.08 && i*par->timestep < 5.5)
		//printf("%f, %E, %E, %E \n", i*par->timestep, Mx[tid], My[tid], Mz[tid]);
		r += Vector3(speed*par->timestep*sin(theta)*cos(phi), speed*par->timestep*sin(theta)*sin(phi), speed*par->timestep*cos(theta));
		//if (i == 1)
		//printf("%f, %f, %f\n", r.x, r.y, r.z);

#if defined SPECULAR_REFLECTION
		boundaryNormal(ri, r, speed, basis, par->timestep);
#elif defined NO_DIFF

#else
		if (!(*basis)->inside(r))
		{
			r = ri;
		}
#endif
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
				par->timestep,
				G_tensor,
				RF_tensor,
				threadIdx.x);

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
//			if (tid == 0){
//				printf("%f %f %f\n", (*B)->getG(r, time).x, (*B)->getG(r, time).y, (*B)->getG(r, time).z);
//				printf("%f %f %f\n", Mx, My, Mz);
//			}
		}
		if (trackingM){
			if (tid < n_mags_track){
				double w = (par->B0).z * GAMMA;
				Mx_track[tid*(*B)->getSteps() + i] = Mx*cos(w*time) - My*sin(w*time);
				My_track[tid*(*B)->getSteps() + i] = Mx*sin(w*time) + My*cos(w*time);
				Mz_track[tid*(*B)->getSteps() + i] = Mz;
//				Mx_track[tid*(*B)->getSteps() + i] = Mx;
//				My_track[tid*(*B)->getSteps() + i] = My;
//				Mz_track[tid*(*B)->getSteps() + i] = Mz;
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
