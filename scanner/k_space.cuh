#ifndef K_SPACE_CU
#define K_SPACE_CU

#include <cufft.h>
#include "../util/recorder.h"
#include <fstream>

class kSpace{
public:
	const int dim_x;
	const int dim_y;
	cufftDoubleComplex *dev_space;
	cufftDoubleComplex **host_space;
	cufftDoubleComplex *dev_result;
	cufftDoubleComplex **host_result;
	size_t host_space_pitch;
	size_t dev_space_pitch;
	double* mag_image;

	kSpace(const int dim_x, const int dim_y);
	size_t index(int x, int y);
	double get_Mx(int x, int y);
	double get_My(int x, int y);
	void set_Mx(int x, int y, double val);
	void set_My(int x, int y, double val);
	void get_fft();
};

#endif
