//
//  OutputData.h
//  LF_DEM
//
//  Created by Ryohei Seto on 6/1/15.
//  Copyright (c) 2015 Ryohei Seto and Romain Mari. All rights reserved.
//

/**
 \class OutputData
 \brief Utility class to output data line-by-line, primitively dimension-aware.
 \author Ryohei Seto
 \author Romain Mari
 */

#ifndef LF_DEM_OutputData_h
#define LF_DEM_OutputData_h

#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <string>
#include <map>

class OutputData {
private:
	int number_of_data;
	bool first_time;
	std::string out_unit;
	std::vector < std::vector<std::string> > output_data;
	std::vector <std::string> output_data_name;
	std::vector <std::string> output_data_type;
	std::map <std::string, double> converter;
	std::ofstream fout;
	
	int getLineNumber() 
	{
		unsigned int line_nb = 0;
		for (const auto& col : output_data) {
			if (line_nb == 0 && col.size() > 0) {
				line_nb = col.size();
			}
			if (col.size() > 0 && col.size() != line_nb) {
				std::cerr << " Error: inconsistent output. Number of lines to output is heterogeneous." << std::endl;
				exit(1);
			} 
		}
		return line_nb;
	}

public:
	OutputData(): first_time(true) {}
	~OutputData() {
		fout.close();
	}
	void setFile (const std::string& fname,
				  const std::string& data_header)
	{
		fout.open(fname.c_str());
		fout << data_header;
	}
	
	void init(int number_of_data_, std::string output_unit) 
	{
		if (first_time) {
			out_unit = output_unit;
			number_of_data = number_of_data_;
			output_data.resize(number_of_data);
			output_data_name.resize(number_of_data);
			output_data_type.resize(number_of_data);
			for (auto& od : output_data) {
				od.clear();
			}
			for (std::string& odn : output_data_name) {
				odn = "blank";
			}
			for (std::string& odt : output_data_type) {
				odt = "none";
			}
		}
	}
	
	void setDimensionlessNumber(double dimensionless_number)
	// dimensionless_number = internal_force_unit/output_force_unit
	{
		if (dimensionless_number == 0) {
			std::cerr << "dimensionless_number (internal_force_unit/output_force_unit) = " << dimensionless_number << std::endl;
			dimensionless_number = 1; // @@@@ To be checked.
		}
		converter["none"] = 1;
		converter["viscosity"] = 6*M_PI;
		converter["stress"] = dimensionless_number;
		converter["time"] = 1/dimensionless_number;
		converter["rate"] = dimensionless_number;
		converter["velocity"] = dimensionless_number;
	}
	
	template<typename T>
	void entryData(int num,
				   std::string name,
				   std::string type,
				   T value)
	{
		int index = num-1;
		std::ostringstream str_value;
		str_value << converter[output_data_type[index]]*value;
		if (first_time) {
			output_data_name[index] = name;
			output_data_type[index] = type;
		}		
		output_data[index].push_back(str_value.str());
	}
	
	void writeToFile(std::string header) 
	{
		fout << header;
		writeToFile();
	}

	void writeToFile()
	{
		int line_nb = getLineNumber();
		if (first_time) {
			fout << "# data in " << out_unit << " units." << std::endl;	
			for (int i=0; i<number_of_data; i++) {
				fout << "#" << i+1 << ": ";
				fout << output_data_name[i];
				fout << std::endl;
			}
			first_time = false;
		}
		for (int i=0; i<line_nb; i++) {
			for (const auto& od : output_data) {
				if (!od.empty()) {
					fout << od[i] << " ";
				} else {
					fout << "n ";
				}
			}
			fout << std::endl;
		}
		for (auto& od : output_data) {
			od.clear();
		}
	}
};
#endif
