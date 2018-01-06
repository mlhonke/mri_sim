//Main simulator macro definitions.
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
	SimuParams test_params(102400, //Number of particles.
		102400,						//Particles per stream/device.
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

	//Create base simulation objects.
	Coil_Ideal test_coil;
	PGSE_D test_sequence(&test_params);
	test_sequence.save_pulse_diagram();
	Lattice test_lattice(3.0, 3.0, 0.5, 100.0, 100.0, 0, 8);
	Scanner test_scanner(test_sequence, test_coil, test_params, test_lattice);

	//Declare shapes to be placed in scanner.
	Cylinder_XY test_primitive(Vector3(-1, -1, 0), 0.5, 0.2, 0.5, 1.0, 0, 1, 0, 10240);
	test_scanner.add_primitive(test_primitive);
	Cylinder_XY test_primitive_2(Vector3(1, -1, 0), 0.5, 0.2, 0.5, 1.0, 0, 1, 0, 10240);
	test_scanner.add_primitive(test_primitive_2);
	Cylinder_XY test_primitive_nose(Vector3(0, 0, 0), 0.5, 0.2, 0.5, 1.0, 0, 1, 0, 10240);
	test_scanner.add_primitive(test_primitive_nose);
	Cylinder_XY test_primitive_mouth(Vector3(-1.25, 1, 0), 0.25, 0.2, 0.5, 1.0, 0, 1, 0, 10240);
	test_scanner.add_primitive(test_primitive_mouth);
	Cylinder_XY test_primitive_mouth_2(Vector3(-0.60, 1.25, 0), 0.25, 0.2, 0.5, 1.0, 0.00001, 1, 0, 10240);
	test_scanner.add_primitive(test_primitive_mouth_2);
	Cylinder_XY test_primitive_mouth_3(Vector3(0, 1.4, 0), 0.25, 0.2, 0.5, 1.0, 0.00002, 1, 0, 10240);
	test_scanner.add_primitive(test_primitive_mouth_3);
	Cylinder_XY test_primitive_mouth_4(Vector3(0.60, 1.25, 0), 0.25, 0.2, 0.5, 1.0, 0.00003, 1, 0, 10240);
	test_scanner.add_primitive(test_primitive_mouth_4);
	Cylinder_XY test_primitive_mouth_5(Vector3(1.25, 1, 0), 0.25, 0.2, 0.5, 1.0, 0.00004, 1, 0, 10240);
	test_scanner.add_primitive(test_primitive_mouth_5);

	//Run the scan!
	test_scanner.scan();

	//Post simulation commands.
	test_scanner.acqs[0]->save_signal("signal");
	test_scanner.acqs[0]->save_tracked("test_track");

	return 0;
}
