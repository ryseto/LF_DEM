#include <string>
#include <iostream>
#include <fstream>
#include "States.h"
#include "CheckFile.h"

namespace State {	
	void outputStateBinary(std::string state_filename, double strain, double _time)
	{
		/**
		 \brief Saves the current state of the simulation in a binary file.
		 
		 Depending on the type of simulation, we store the data differently, defined by
		 binary format version numbers, which is always the first data.
		 */
		std::ofstream state_export;
		state_export.open(state_filename.c_str(), std::ios::binary | std::ios::out);
		
		typedef std::underlying_type<StateFileFormat>::type format_type;
		format_type binary_format = static_cast<format_type>(StateFileFormat::basic);
		state_export.write((char*)&binary_format, sizeof(format_type));
		state_export.write((char*)&strain, sizeof(double));
		state_export.write((char*)&_time, sizeof(double));
		state_export.close();
	}
	
	StateFileFormat getFileFormat(std::ifstream &input)
	{
		typedef std::underlying_type<StateFileFormat>::type format_type;
		format_type format_raw;
		input.read((char*)&format_raw, sizeof(format_type));
		return static_cast<StateFileFormat>(format_raw);
	}
	
	StateFileFormat getFileFormat(const std::string &filename)
	{
		checkInFile(filename);
		std::ifstream input(filename.c_str(), std::ios::binary | std::ios::in);
		
		return getFileFormat(input);
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
		
		auto format = getFileFormat(input);
		if (format != StateFileFormat::basic) {
			throw std::runtime_error("readBasicCheckpoint(): got incorrect binary format.");
		}
		BasicCheckpoint chkp;
		input.read((char*)&chkp.clock.cumulated_strain, sizeof(decltype(chkp.clock.cumulated_strain)));
		input.read((char*)&chkp.clock.time_, sizeof(double));
		return chkp;
	}
	
	bool isZeroTimeChkp(BasicCheckpoint chkp){
		return chkp.clock.time_ == 0 \
		&& chkp.clock.cumulated_strain == 0;
	}
}
