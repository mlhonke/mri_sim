#include "CPUkernels.cuh"
#include <iostream>

void updateWalkersMagCPU(SimuParams *par, Primitive* basis, const Sequence* B,
		Coil* coil, int n_mags_track, real *signal_x, real *signal_y,
		real *signal_z) {

	for (int tid = 0; tid < par->number_of_particles; tid++) {
		//printf("Current particle on CPU is: %d\n", tid);
		real Mx = 0;
		real My = 0;
		real Mz = 0;
		int s = 0;

		real phi, theta;
		real speed = sqrt(6.0 * basis->getD() / par->timestep);

		Vector3 r = basis->unifRandCPU();
		//printf("%f, %f, %f\n", r.x, r.y, r.z);

		for (int i = 0; i < par->measurements; i++) {
			Mx = par->mx_initial;
			My = par->my_initial;
			Mz = par->mz_initial;

			updateMagRK4(r, &Mx, &My, &Mz, coil, B, par->B0, basis->getT1(),
					basis->getT2(), i, par->timestep);
		}

		for (int i = 1; i < par->steps; i++) {
			real time = i * par->timestep;
			Vector3 ri = r;

			phi = 2.0 * PI * unifRandCPP();
			theta = acos(2.0 * unifRandCPP() - 1);

			r += Vector3(speed * par->timestep * sin(theta) * cos(phi),
					speed * par->timestep * sin(theta) * sin(phi),
					speed * par->timestep * cos(theta));

#if defined SPECULAR_REFLECTION
			boundaryNormalCPU(ri,r, speed, basis, par->timestep);
#else
			if (!basis->inside(r)) {
				r = ri;
			}
#endif
			for (int j = 0; j < par->measurements; j++) {
				updateMagRK4(r, &Mx, &My, &Mz, coil, B, par->B0, basis->getT1(),
						basis->getT2(), i, par->timestep);
			}
			if (time >= (B)->getReadStart(time)
					&& time < (B)->getReadFinish(time)
					&& ((int) (i - (B)->getReadStart(time) / par->timestep))
							% (B)->getReadFactor() == 0) {

				double w = (par->B0).z * GAMMA;

				//Save the total signal from this block to global memory for later summation.
				signal_x[s] += Mx * cos(w * time) - My * sin(w * time);
				signal_y[s] += Mx * sin(w * time) + My * cos(w * time);
				signal_z[s] += Mz;
				s++;
			}
		}
	}
}
