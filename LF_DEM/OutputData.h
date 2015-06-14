//
//  OutputData.h
//  LF_DEM
//
//  Created by Ryohei Seto on 6/1/15.
//  Copyright (c) 2015 Ryohei Seto. All rights reserved.
//

#ifndef LF_DEM_OutputData_h
#define LF_DEM_OutputData_h

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <string>

using namespace std;

class OutputData {
private:
	int number_of_data;
	bool first_time;
	vector<string> output_data;
	vector<string> output_data_name;
public:
	OutputData():
	first_time(true) {}

	void init(int number_of_data_) {
		if (first_time) {
			number_of_data = number_of_data_;
			output_data.resize(number_of_data);
			output_data_name.resize(number_of_data);
			for (string &od : output_data) {
				od = "n";
			}
			for (string &odn : output_data_name) {
				odn = "blank";
			}
		}
	}
	
	template<typename T>
	void entryData(int num, string name, T value)
	{
		int index = num-1;
		ostringstream str_value;
		str_value << value;
		if (first_time) {
			if (output_data_name[index] != "blank") {
				cerr << "data["<< index << "] is redefined." << endl;
				exit(1);
			} else {
				output_data_name[index] = name;
			}
		}
		output_data[index] = str_value.str();
	}
	
	void exportFile(ofstream &fout_data)
	{
		if (first_time) {
			for (int i=0; i<number_of_data; i++) {
				fout_data << "#" << i+1 << ": ";
				fout_data << output_data_name[i];
				fout_data << endl;
			}
			first_time = false;
		}
		for (int i=0; i<number_of_data; i++) {
			fout_data << output_data[i] << ' ';
		}
		fout_data << endl;
	}
};
#endif
