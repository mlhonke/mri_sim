#include "../util/vector3.cuh"
#include "../primitives/primitive.cuh"
#include "../primitives/lattice.cuh"

__device__ void boundaryNormal(Vector3 & ri, Vector3 & r, const real currentSpeed,  const Primitive** basis, const real timestep){

	real v = 0.0;
	real accumtime = 0.0;
	Vector3 n;
	//printf("r_i: %f %f\n", ri.x, ri.y);
	while ((*basis)->intersect(ri, r, v, n)){
				r = (r-ri)*v + ri;

				//printf("rbn: %f %f %f\n", r.x, r.y, r.z);
				//printf("r_mag: %f\n", sqrt(r.x*r.x + r.y*r.y));

				//determine incident vector
				Vector3 I = (r-ri);
				//printf("I: %f\n", I);
				real Imag = I.magnitude();

				//calculate accumulated time
				accumtime += Imag/currentSpeed;
				//printf("accumtime: %f\n", accumtime);

				//printf("n: %f %f\n", n.x, n.y);

				//set last position to boundary position before reflection
				ri = r;

				//reflect, reflection travel time is the remaining time (timestep - accumtime)
				r += (I - n*2.0*(n*I))*(currentSpeed / Imag)*(timestep - accumtime);

				//printf("r_mag_f: %f\n", sqrt(r.x*r.x + r.y*r.y));
				//printf("r_f_l: %f %f %f\n", r.x, r.y, r.z);
		}

	//printf("r_f: %f %f\n", r.x, r.y);

}

__device__ void boundaryCheckImpermeable(const Lattice ** lat, Primitive *** basis, const Vector3 & ri, Vector3 & r, const real currentSpeed,  const real timestep){
	real v = 0.0;
	real accumtime = 0.0;
	Vector3 n;
	Vector3 riCurrent = ri;

	while ((*lat)->intersectionCheckImpermeable(basis, riCurrent, r, v, n)){

		r = (r-riCurrent)*v + riCurrent;

		//determine incident vector
		Vector3 I = (r-riCurrent);
		real Imag = I.magnitude();

		//calculate accumulated time
		accumtime += Imag/currentSpeed;

		//set last position to boundary position before rflection
		riCurrent = r;

		//reflect, reflection travel time is the remaining time (timestep - accumtime)
		r += ( I - n*2.0*(n*I) )*(currentSpeed/Imag)*( timestep - accumtime );
	}
}

__device__ void boundaryCheckPermeable(const Lattice ** lat, Primitive *** basis, const Vector3 & ri, Vector3 & r, real accumtime, const real timestep, curandState & localState){

	real v = 0.0;
	Vector3 n;
	bool transmitted = false;
	real permeability;
	real currentSpeed = sqrt(6.0*(*lat)->getD(basis,ri)/timestep);
	Vector3 riCurrent = ri;

	while ((*lat)->intersectionCheckPermeable(basis, riCurrent, r, v, n, permeability) && transmitted == false){
		real transmissionSpeed = sqrt(6.0*(*lat)->getD(basis,r)/timestep);

		r = (r-riCurrent)*v + riCurrent;

		//determine incident vector
		Vector3 I = (r-riCurrent);
		real Imag = I.magnitude();

		//calculate accumulated time
		accumtime += Imag/currentSpeed;

		//reflect
		if (curand_uniform( &localState ) > permeability*4.0/currentSpeed){

			//set last position to boundary position before rflection
			riCurrent = r;

			//reflect, reflection travel time is the remaining time (timestep - accumtime)
			r += ( I - n*2.0*(n*I) )*(currentSpeed/Imag)*( timestep - accumtime );
		}
		else{

		r += (I/Imag)*transmissionSpeed*(timestep - accumtime);
		transmitted = true;
		}
	}
}
