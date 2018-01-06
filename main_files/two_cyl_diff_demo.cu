//Main simulator library.
#include "master_def.h"

//Specific coil, sequence... for this simulation.
#include <iostream>
#include "sequence/PGSE_D.cuh"
#include "coil/coil_ideal.cuh"
#include "scanner/scanner.cuh"
#include "primitives/CylinderXY.cuh"
#include "params/simuParams.cuh"
#include "util/recorder.h"
#include "util/vector3.cuh"
#include "primitives/lattice.cuh"

int main(){
	//Simulation properties.
	SimuParams test_params(10240, //Number of particles.
		10240,						//Particles per stream/device.
		2.0,						//Sequence repeat time.
		0.5,						//Sequence echo time.
		0.001,						//Simulation timestep.
		1,							//Number of particles to track continual, individual magnetization.
		Vector3(0, 0, 1),			//Initial magnetization vector.
		Vector3(0, 0, 0.001),		//Main B0 field direction / strength.
		64,							//(vertical) resolution.
		64,							//(horizontal) resolution.
		4,							//(vertical) FOV.
		4							//(horizontal) FOV.
		);

	Coil_Ideal test_coil;
	PGSE_D test_sequence(&test_params);
	std::cout << "Done Sequence Initialization" << std::endl;
	Lattice test_lattice(3.0, 3.0, 0.5, 100.0, 100.0, 0.0, 2);
	Scanner test_scanner(test_sequence, test_coil, test_params, test_lattice);

	recorder record_p("pulse");
	ofstream trial_p = record_p.setup_record_csv();

	Cylinder_XY test_primitive(Vector3(-1, 0, 0), 0.5, 0.2, 0.5, 0.5, 0, 1, 0, 10240);
	test_scanner.add_primitive(test_primitive);
	Cylinder_XY test_primitive_2(Vector3(1, 0, 0), 0.5, 0.2, 0.5, 0.5, 0.0001, 1, 0, 10240);
	test_scanner.add_primitive(test_primitive_2);
	safe_cuda(cudaGetLastError());

	std::cout << test_sequence.getNSubSequences() << std::endl;

	test_scanner.scan();

	//Post simulation commands.
	test_scanner.acqs[0]->save_signal("signal");
	test_scanner.acqs[0]->save_tracked("test_track");

	return 0;
}
