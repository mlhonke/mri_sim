//Main simulator library.
#include "master_def.h"

//Specific coil, sequence... for this simulation.
#include <iostream>
#include "sequence/GRE.cuh"
#include "coil/coil_ideal.cuh"
#include "scanner/scanner.cuh"
#include "primitives/CylinderXY.cuh"
#include "params/simuParams.cuh"
#include "util/recorder.h"
#include "util/vector3.cuh"

int main(){
	//Simulation properties.
	SimuParams test_params(10240, //Number of particles.
		10240,					//Number of particles per stream.
		2.5,						//Sequence repeat time.
		0.5,						//Sequence echo time.
		0.001,						//Simulation timestep.
		0,							//Number of particles to track continual, individual magnetization.
		Vector3(0, 0, 1),			//Initial magnetization vector.
		Vector3(0, 0, 0.001),		//Main B0 field direction / strength.
		64,							//(vertical) resolution.
		64,							//(horizontal) resolution.
		5,							//(vertical) FOV.
		5							//(horizontal) FOV.
		);

	Coil_Ideal test_coil;
	GRE test_sequence(&test_params);
	test_sequence.save_pulse_diagram();
	Scanner test_scanner(test_sequence, test_coil, test_params);
	Cylinder_XY test_primitive(Vector3(0.0,0.0,0.0), 2.0, 2.0, 2.0, 4.0, 0, 0, 0, 10240);
	test_scanner.add_primitive(test_primitive);

	test_scanner.scan();

	//Post simulation commands.
	test_scanner.acqs[0]->save_signal("signal");
	test_scanner.acqs[0]->save_tracked("test_track");

//	Vector3 pos_vector;
//	recorder pos("positions");
//	ofstream pos_out = pos.setup_record_csv();
//	for (int i = 0; i < 1024; i++){
//		pos_vector = test_primitive.unifRandCPU();
//		pos_out << pos_vector.x << ',' << pos_vector.y << ',' << pos_vector.z << std::endl;
//	}

	cudaDeviceSynchronize();
	cudaDeviceReset();
	return 0;
}
