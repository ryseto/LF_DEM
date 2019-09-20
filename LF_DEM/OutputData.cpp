#include "OutputData.h"
#include "CheckFile.h"

void OutputData::setFile(const std::string& fname,
						 const std::string& data_header,
						 bool force_overwrite,
						 bool append)
{
	if (!force_overwrite && !append) {
		std::ifstream file_test(fname.c_str());
		if (file_test.good()) {
			std::cerr << "The file '" << fname << "' already exists." << std::endl;
			std::cerr << "You need -f option to overwrite." << std::endl;
			exit(1);
		}
	}
	if (append) {
		checkInFile(fname);
		fout.open(fname.c_str(), std::ofstream::out | std::ofstream::app);
		if(fout.good()) {
 			std::cerr << "Successfully opening file '" << fname << "'"<< std::endl;
		}
		restart_from_chkp = true;
	} else {
		fout.open(fname.c_str(), std::ofstream::out);
		fout << data_header;
	}
}

unsigned OutputData::getLineNumber()
{
	unsigned line_nb = 0;
	for (const auto& data : output_data) {
		const auto &col = data.second;
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

void OutputData::initCol(std::string name,
						 unsigned width)
{
	if (output_data.find(name) != output_data.end()) {
		return;
	}
	std::vector <std::string> col;
	col.clear();
	output_data[name] = col;
	output_data_name.push_back(name);
	output_data_width[name] = width;
	insert_order.push_back(name);
}

void OutputData::writeFileHeader()
{
	fout << "# data in " << Dimensional::unit2suffix(out_unit) << " units." << std::endl;
	int i = 1;
	for (const auto &name : insert_order) {
		int width = output_data_width[name];
		if (width == 1) {
			fout << "#" << i << ": ";
		} else {
			fout << "#" << i << "-" << i+width-1 << ": ";
		}
		i += width;
		fout << name;
		fout << std::endl;
	}
	fout << std::endl;
}

void OutputData::writeColsToFile()
{
	int line_nb = getLineNumber();
	for (int i=0; i<line_nb; i++) {
		for (const auto& name : insert_order) {
			const auto &col = output_data[name];
			if (!col.empty()) {
				fout << col[i] << " ";
			} else {
				fout << "n ";
			}
		}
		fout << std::endl;
	}
	for (auto& od : output_data) {
		od.second.clear();
	}
}

void OutputData::writeToFile(std::string header)
{
	if (first_time) {
		if (!restart_from_chkp) {
			writeFileHeader();
		}
		first_time = false;
	}
	fout << header;
	writeColsToFile();
}

void OutputData::writeToFile()
{
	if (first_time) {
		if (!restart_from_chkp) {
			writeFileHeader();
		}
		first_time = false;
	}
	writeColsToFile();
}

void OutputData::setUnits(Dimensional::UnitSystem units_,
						  Dimensional::Unit output_unit)
{
	units = units_;
	internal_unit = units.getInternalUnit();
	out_unit = output_unit;
}
