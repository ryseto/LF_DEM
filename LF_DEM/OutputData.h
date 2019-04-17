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
#include "DimensionalQty.h"

class OutputData {
private:
	bool first_time;
	bool restart_from_chkp;
	Dimensional::Unit out_unit;
	Dimensional::Unit internal_unit;
	Dimensional::UnitSystem units;
	std::map <std::string, std::vector<std::string> > output_data;
	std::vector <std::string> insert_order;
	std::vector <std::string> output_data_name;
	std::map <std::string, unsigned> output_data_width;
	std::ofstream fout;
	unsigned default_precision;

	unsigned getLineNumber();

	void initCol(std::string name, unsigned width);

public:
	OutputData(): first_time(true), restart_from_chkp(false), default_precision(6) {}
	~OutputData()
	{
		fout.close();
	}
	void setFile(const std::string& fname,
				 const std::string& data_header,
				 bool force_overwrite=false,
				 bool append=false);

	void setDefaultPrecision(unsigned precision)
	{
		default_precision = precision;
	}

	void setUnits(Dimensional::UnitSystem units_,
				  Dimensional::Unit output_unit);

	template<typename T>
	void entryData(std::string name,
				   Dimensional::Dimension dimension,
				   unsigned width,
				   T value,
				   unsigned precision=0);

	void writeFileHeader();
	void writeColsToFile();
	void writeToFile(std::string header);
	void writeToFile();
	template<typename T>
	T convertToOutput(Dimensional::Dimension dimension, T value);
};


template<typename T>
inline void OutputData::entryData(std::string name,
								  Dimensional::Dimension dimension,
								  unsigned width,
								  T value,
								  unsigned precision)
{
	if (first_time) {
		initCol(name, width);
	}
	unsigned output_precision;
	if (precision > 0) {
		output_precision = precision;
	} else {
		output_precision = default_precision;
	}
	std::ostringstream str_value;
	if (dimension != Dimensional::Dimension::none) {
		Dimensional::DimensionalQty<T> qty = {dimension, value, internal_unit};
		units.convertFromInternalUnit(qty, out_unit);
		str_value << std::setprecision(output_precision) << qty.value;
	} else {
		str_value << std::setprecision(output_precision) << value;
	}
	output_data[name].push_back(str_value.str());
}

template<typename T>
inline T OutputData::convertToOutput(Dimensional::Dimension dimension,
									 T value)
{
	if (dimension != Dimensional::Dimension::none) {
		Dimensional::DimensionalQty<T> qty = {dimension, value, internal_unit};
		units.convertFromInternalUnit(qty, out_unit);
		return qty.value;
	} else {
		return value;
	}
}

#endif
