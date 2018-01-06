/*
	RFpulse simulates the intended "ideal" effect of arbitrary angled pulses.
	Author: Michael Honke
	Date: June 23, 2016
*/

//#include "../dSim_master.cuh"

/*
	Pulse (B1) along the x direction, so rotate about the x-axis.
	Input: Pointers to Mx, My, Mz from the simulation. angle is in radians.
	Output: None, magnetization vectors are updated directly.
*/
__device__ __host__ void pulse_x(real* Mx, real* My, real* Mz, real angle){
	//Create a local scope x,y,z to store temporary results, since will have rotated My, before rotated Mz for example.
		real x = *Mx;
		real y = *My;
		real z = *Mz;

	/*
		Calculate and assign Mx, My, Mz new rotated values. Mx is unaffected.
		For numerical accuracy special cases PI and PI/2 pulses, which are 
		desired to be exact are implemented. For example sin(PI) != 0.000...
		due to PI being irrational.
	*/
		if (angle == PI){
			*My = -y;
			*Mz = -z;
		}
		else if ((angle <= PI / 2 + PI_THRESH) && (angle >= PI / 2 - PI_THRESH)){
			*My = -z;
			*Mz = y;
		}
		else {
			*My = cos(angle)*y - sin(angle)*z;
			*Mz = sin(angle)*y + cos(angle)*z;
		}
}
