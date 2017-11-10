#include <string>
#include <iostream>
#include <fstream>
#include "System.h"
#include "States.h"

namespace State {

void outputStateBinary(std::string state_filename, const System &sys)
{
	/**
	 \brief Saves the current state of the simulation in a binary file.
	 
	 Depending on the type of simulation, we store the data differently, defined by
	 binary format version numbers, which is always the first data.
	 */
	std::ofstream state_export;
	state_export.open(state_filename.c_str(), std::ios::binary | std::ios::out);
	
	typedef std::underlying_type<StateFileFormat>::type format_type;
	format_type binary_format = 1;
	state_export.write((char*)&binary_format, sizeof(format_type));
	double strain = sys.get_cumulated_strain();
	state_export.write((char*)&strain, sizeof(double));
	double _time = sys.get_time();
	state_export.write((char*)&_time, sizeof(double));
	double _time_simu = sys.get_time_in_simulation_units();
	state_export.write((char*)&_time_simu, sizeof(double));
	state_export.close();
}

StateFileFormat getFileFormat(const std::string &filename)
{
	checkInFile(filename);
	std::ifstream input(filename.c_str(), std::ios::binary | std::ios::in);
	
	typedef std::underlying_type<StateFileFormat>::type format_type;
	format_type format_raw;
	input.read((char*)&format_raw, sizeof(format_type));
	return static_cast<StateFileFormat>(format_raw);
}

struct BasicCheckpoint readBasicCheckpoint(const std::string &filename)
{
	/**
	 \brief Saves the current state of the simulation in a binary file.
	 
	 Depending on the type of simulation, we store the data differently, defined by
	 binary format version numbers, which is always the first data.
	 */
	checkInFile(filename);
	std::ifstream input(filename.c_str(), std::ios::binary | std::ios::in);
	
	auto format = getFileFormat(filename);
	if (format != StateFileFormat::basic) {
		throw std::runtime_error("readBasicCheckpoint(): got incorrect binary format.");
	}
	BasicCheckpoint chkp;
	typedef std::underlying_type<StateFileFormat>::type format_type;
	format_type fmt;
	input.read((char*)&fmt, sizeof(format_type));
	input.read((char*)&chkp.clock.cumulated_strain, sizeof(decltype(chkp.clock.cumulated_strain)));
	input.read((char*)&chkp.clock.time_, sizeof(double));
	input.read((char*)&chkp.clock.time_in_simulation_units, sizeof(double));
	return chkp;
}

bool isZeroTimeChkp(BasicCheckpoint chkp){
	return chkp.clock.time_ == 0 \
	&& chkp.clock.time_in_simulation_units == 0 \
	&& chkp.clock.cumulated_strain == 0;
}
	
}
