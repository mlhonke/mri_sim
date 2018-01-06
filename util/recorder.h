/*
recorder.h Responsible for making record keeping of simulations convenient.

Author: Michael Honke
*/

#ifndef _RECORDER_H_
#define _RECORDER_H_
#define DIR_REC "/home/omega/Dev/dSimL/trials/"

using namespace std;

#include <string>
#include <fstream>
#include <sstream>
#include <ctime>
#include <cmath>
#include "../master_def.h"

class recorder{
public:
	ofstream setup_record_csv();
	ofstream setup_record_image();
	ofstream setup_record_csv(string custom_name);
	void txt_record(string file_name);
	template<typename T> void save_csv(T**data, int res_x, int res_y);
	template<typename T> void save_csv_complex(T**data, char mode, int res_x, int res_y);
	template<typename T> void save_image(T**data, int res_x, int res_y);
	template<typename T> void save_image_complex(T**data, char mode, int res_x, int res_y);
	recorder(string exp_name);
private:
	string _name;
	string _exp_name;
	string make_name(string ext);
};

template <typename T> void recorder::save_csv(T **data, int res_x, int res_y){
	ofstream output = recorder::setup_record_csv();

	for (int i = 0; i < res_x; i++){
		for (int j = 0; j < res_y; j++){
			output << data[i][j] << ",";
		}
		output << std::endl;
	}
}

template <typename T> void recorder::save_csv_complex(T **data, char mode, int res_x, int res_y){
	ofstream output = recorder::setup_record_csv();

	if (mode == 'r'){
		for (int i = 0; i < res_x; i++){
			for (int j = 0; j < res_y; j++){
				output << data[i][j].x << ",";
			}
			output << std::endl;
		}
	} else if (mode == 'i') {
		for (int i = 0; i < res_x; i++){
			for (int j = 0; j < res_y; j++){
				output << data[i][j].y << ",";
			}
			output << std::endl;
		}
	} else {
		for (int i = 0; i < res_x; i++){
			for (int j = 0; j < res_y; j++){
				output << std::sqrt( std::pow(data[i][j].x, 2) + std::pow(data[i][j].y, 2)) << ",";
			}
			output << std::endl;
		}
	}
}

template <typename T> void recorder::save_image(T **data, int res_x, int res_y){
	ofstream output = recorder::setup_record_image();

	for (int i = 0; i < res_x; i++){
		for (int j = 0; j < res_y; j++){
			output << data[i][j] << " ";
		}
		output << std::endl;
	}
}

template <typename T> void recorder::save_image_complex(T **data, char mode, int res_x, int res_y){
	ofstream output = recorder::setup_record_image();

	output << "P2" << std::endl;
	output << res_x << " " << res_y << std::endl;

	real max = 0;
	if (mode == 'r' || mode == 'i'){
		for (int i = 0; i < res_x; i++){
			for (int j = 0; j < res_y; j++){
				if (data[i][j].x > max)
					max = data[i][j].x;
				if (data[i][j].y > max)
					max = data[i][j].y;
			}
		}
	} else {
		for (int i = 0; i < res_x; i++){
			for (int j = 0; j < res_y; j++){
				if (std::sqrt( std::pow(data[i][j].x, 2) + std::pow(data[i][j].y, 2)) > max)
					max = std::sqrt( std::pow(data[i][j].x, 2) + std::pow(data[i][j].y, 2));
			}
		}
	}

	output << (int) 10240 << std::endl;
	real scale = max / 10240;

	if (mode == 'r'){
		for (int i = 0; i < res_x; i++){
			for (int j = 0; j < res_y; j++){
				output << (int) (data[i][j].x / scale) << " ";
			}
			output << std::endl;
		}
	} else if (mode == 'i') {
		for (int i = 0; i < res_x; i++){
			for (int j = 0; j < res_y; j++){
				output << (int) (data[i][j].y / scale) << " ";
			}
			output << std::endl;
		}
	} else {
		for (int i = 0; i < res_x; i++){
			for (int j = 0; j < res_y; j++){
				output << (int) (std::sqrt( std::pow(data[i][j].x, 2) + std::pow(data[i][j].y, 2)) / scale) << " ";
			}
			output << std::endl;
		}
	}
}

#endif
