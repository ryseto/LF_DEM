#ifndef __LF_DEM__CheckFile__
#define __LF_DEM__CheckFile__
#include <sstream>
#include <fstream>
#include <stdexcept>

inline void checkInFile(std::string filename)
{
	std::ifstream file_import(filename.c_str(), std::ios::in);
	if (!file_import) {
		std::ostringstream error_str;
		error_str  << " File '" << filename << "' not found." << std::endl;
		throw std::runtime_error(error_str.str());
	}
}

#endif
