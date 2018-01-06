#include "recorder.h"

recorder::recorder(string exp_name){
	_exp_name = exp_name;
}

string recorder::make_name(string ext){
	string name;
	time_t rawtime;
	struct tm * timeinfo;
	char buffer[80];
	time(&rawtime);
	timeinfo = localtime(&rawtime);
	strftime(buffer, 80, " %Y_%j_%H_%M_%S", timeinfo);
	string time_str(buffer);
	name = DIR_REC + _exp_name + time_str + ext;

	return name;
}

bool already_exists(string file_name){
	ifstream attempt(file_name);
	return attempt.good();
}

ofstream recorder::setup_record_csv(){
	string file_name = make_name(".csv");
	int t_dup = 0;

	while (already_exists(file_name)){
		t_dup++;
		file_name = make_name("_" + std::to_string(t_dup) + ".csv");
	}

	ofstream trial(file_name);

	return trial;
}

ofstream recorder::setup_record_image(){
	string file_name = make_name(".pgm");
	int t_dup = 0;

	while (already_exists(file_name)){
		t_dup++;
		file_name = make_name("_" + std::to_string(t_dup) + ".pgm");
	}

	ofstream trial(file_name);

	return trial;
}

ofstream recorder::setup_record_csv(string custom_name){
	string file_name = DIR_REC + custom_name + ".csv";
	int t_dup = 0;

	while (already_exists(file_name)){
		t_dup++;
		file_name = DIR_REC + custom_name + "_" + std::to_string(t_dup) + ".csv";
	}

	ofstream trial(file_name);

	return trial;
}

void recorder::txt_record(string file_name){
	std::ifstream infile(file_name, std::ifstream::binary);
	std::ofstream outfile(make_name(".txt"), std::ofstream::binary);

	// get size of file
	infile.seekg(0, infile.end);
	long size = infile.tellg();
	infile.seekg(0);

	// allocate memory for file content
	char* buffer = new char[size];

	// read content of infile
	infile.read(buffer, size);

	// write to outfile
	outfile.write(buffer, size);

	// release dynamically-allocated memory
	delete[] buffer;

	outfile.close();
	infile.close();
}
