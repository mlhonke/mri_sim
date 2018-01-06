#ifndef _MAG_ACQUISITION_H_
#define _MAG_ACQUISITION_H_

#include "../params/simuParams.cuh"
#include "../sequence/sequence.cuh"
#include "../util/recorder.h"
#include <vector>
#include <string>

class magAcquisition {

private:

	int numOfMeasurements;
	int numOfSteps;
	int n_mags_track = 0;
	int lastAllocMeasurement;
	int seed;
	int readSteps;
	int steps;

	std::vector<real> signal_x;
	std::vector<real> signal_y;
	std::vector<real> signal_z;
	std::vector<real> mx_tracked;
	std::vector<real> my_tracked;
	std::vector<real> mz_tracked;

	real ADC;
	real ADCmeansquare;
	int points;

public:

	magAcquisition(const SimuParams* params, const Sequence* _sequence){
		points = 0;
		readSteps = _sequence->getReadSteps();
		signal_x = std::vector<real>(readSteps);
		signal_y = std::vector<real>(readSteps);
		signal_z = std::vector<real>(readSteps);
		numOfMeasurements = params->measurements;
		seed = params->seed;
		numOfSteps = params->steps;
		lastAllocMeasurement = 0;
		steps = _sequence->getSteps();
	}

	//MH: Extra to track magnetization.
	void set_tracked_particles(int _n_mags_track){
		n_mags_track = _n_mags_track;
		mx_tracked = std::vector<real>(_n_mags_track*numOfSteps);
		my_tracked = std::vector<real>(_n_mags_track*numOfSteps);
		mz_tracked = std::vector<real>(_n_mags_track*numOfSteps);
	}

	int get_n_mags_track(){
		return n_mags_track;
	}

	std::vector<real> & get_signal_x(){
		return signal_x;
	}

	std::vector<real> & get_signal_y(){
		return signal_y;
	}

	std::vector<real> & get_signal_z(){
		return signal_z;
	}

	std::vector<real> & getMxTracked(){
		return mx_tracked;
	}

	std::vector<real> & getMyTracked(){
		return my_tracked;
	}

	std::vector<real> & getMzTracked(){
		return mz_tracked;
	}

	int & getSeed(){
		return seed;
	}

	void save_signal(string name){
		recorder record(name);
		ofstream trial = record.setup_record_csv();

		for (int i = 0; i < readSteps; i++){
			trial << i << "," << get_signal_x()[i] << "," << get_signal_y()[i] << "," << get_signal_z()[i] << std::endl;
		}
	}

	void save_tracked(string name){
		recorder record(name);
		ofstream trial = record.setup_record_csv();

		int skip_factor = 10;
		for (int i = 0; i < (steps/skip_factor)*n_mags_track; i++){
			trial << i << "," << getMxTracked()[i*skip_factor] << "," << getMyTracked()[i*skip_factor] << "," << getMzTracked()[i*skip_factor] << std::endl;
		}
	}
};

#endif
