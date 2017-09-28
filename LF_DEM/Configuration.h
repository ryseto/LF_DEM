#include <stdexcept>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <type_traits>
#include "global.h"
#include "Contact.h"
#include "vec3d.h"

#ifndef __LF_DEM__Configuration__
#define __LF_DEM__Configuration__


enum struct ConfFileFormat : int { // assign values as it is for output and otherwise may be compiler dependent.
	bin_format_base_old = 1,
	bin_format_base_new = 2,
	bin_format_fixed_vel = 3,
	txt_format_base_old = 4,
	txt_format_base_new = 5,
	txt_format_fixed_vel = 6,
	txt_format_circular_couette = 7
};

struct base_configuration {
	double lx, ly, lz;
	double volume_or_area_fraction;
	vec3d lees_edwards_disp;

	std::vector <vec3d> position;
	std::vector <double> radius;
	std::vector <double> angle;

	std::vector <struct contact_state> contact_states;
};

struct fixed_velo_configuration {
	double lx, ly, lz;
	double volume_or_area_fraction;
	vec3d lees_edwards_disp;

	std::vector <vec3d> position;
	std::vector <double> radius;
	std::vector <double> angle;
	std::vector<vec3d> fixed_velocities;

	std::vector <struct contact_state> contact_states;
};

struct circular_couette_configuration {
	int np_wall1;
	int np_wall2;
	double radius_in;
	double radius_out;

	double lx, ly, lz;
	double volume_or_area_fraction;
	vec3d lees_edwards_disp;

	std::vector <vec3d> position;
	std::vector <double> radius;
	std::vector <double> angle;

	std::vector <struct contact_state> contact_states;
};

inline std::vector <struct contact_state> readContactStatesBStream(std::istream &input, unsigned int np)
{
	int ncont;
	double dt_x, dt_y, dt_z, dr_x, dr_y, dr_z;
	std::vector <struct contact_state> cont_states;
	input.read((char*)&ncont, sizeof(unsigned int));
	std::iostream::pos_type file_pos = input.tellg();
	bool ushort_format = false;
	for (int i=0; i<ncont; i++) {
		unsigned int p0, p1;
		input.read((char*)&p0, sizeof(unsigned int));
		// hacky thing to guess if this is an old format with particle numbers as unsigned short
		if(p0>np){
			ushort_format = true;
			input.seekg(file_pos);
			break;
		}
		input.read((char*)&p1, sizeof(unsigned int));
		input.read((char*)&dt_x, sizeof(double));
		input.read((char*)&dt_y, sizeof(double));
		input.read((char*)&dt_z, sizeof(double));
		input.read((char*)&dr_x, sizeof(double));
		input.read((char*)&dr_y, sizeof(double));
		input.read((char*)&dr_z, sizeof(double));
		struct contact_state cs;
		cs.p0 = p0;
		cs.p1 = p1;
		cs.disp_tan = vec3d(dt_x, dt_y, dt_z);
		cs.disp_rolling = vec3d(dr_x, dr_y, dr_z);
		cont_states.push_back(cs);
	}
	if (ushort_format) {
		for (int i=0; i<ncont; i++) {
			unsigned short p0, p1;
			input.read((char*)&p0, sizeof(unsigned short));
			input.read((char*)&p1, sizeof(unsigned short));
			input.read((char*)&dt_x, sizeof(double));
			input.read((char*)&dt_y, sizeof(double));
			input.read((char*)&dt_z, sizeof(double));
			input.read((char*)&dr_x, sizeof(double));
			input.read((char*)&dr_y, sizeof(double));
			input.read((char*)&dr_z, sizeof(double));
			struct contact_state cs;
			cs.p0 = (int)p0;
			cs.p1 = (int)p1;
			cs.disp_tan = vec3d(dt_x, dt_y, dt_z);
			cs.disp_rolling = vec3d(dr_x, dr_y, dr_z);
			cont_states.push_back(cs);
		}
	}
	return cont_states;
}

inline std::pair<std::vector <vec3d>, std::vector <double>> readPositionsBStream(std::istream &input, unsigned int np)
{
	double x_, y_, z_, r_;
	std::vector <vec3d> position;
	std::vector <double> radius;
	for (unsigned int i=0; i<np; i++) {
		input.read((char*)&x_, sizeof(double));
		input.read((char*)&y_, sizeof(double));
		input.read((char*)&z_, sizeof(double));
		input.read((char*)&r_, sizeof(double));
		position.push_back(vec3d(x_,y_,z_));
		radius.push_back(r_);
	}
	return std::make_pair(position, radius);
}


inline std::string getMetaParameter(std::map<std::string,std::string> &meta_data,
                                    std::string &key)
{
	if (meta_data.find(key) != meta_data.end()) {
		return meta_data[key];
	} else {
		std::ostringstream error_str;
		error_str  << " Simulation:: parameter '" << key << "' not found in the header of the configuration file." << std::endl;
		throw std::runtime_error(error_str.str());
	}
}

inline std::string getMetaParameter(std::map<std::string,std::string> &meta_data,
                                    std::string &key,
                                    const std::string &default_val)
{
	if (meta_data.find(key) != meta_data.end()) {
		return meta_data[key];
	} else {
		return default_val;
	}
}

inline std::map<std::string,std::string> getConfMetaData(const std::string &line1, const std::string &line2)
{
	std::vector<std::string> l1_split = splitString(line1);
	std::vector<std::string> l2_split = splitString(line2);
	if (l1_split.size() != l2_split.size()) {
		throw std::runtime_error("Simulation:: Ill-formed header in the configuration file.\n");
	}
	std::map<std::string,std::string> meta_data;
	for (unsigned int i=1; i<l1_split.size(); i++) {
		meta_data[l1_split[i]] = l2_split[i];
	}
	return meta_data;
}

inline ConfFileFormat getBinaryConfigurationFileFormat(const std::string& filename_import_positions)
{
	checkInFile(filename_import_positions);
	std::ifstream file_import;
	file_import.open(filename_import_positions.c_str(), std::ios::binary | std::ios::in);
	typedef std::underlying_type<ConfFileFormat>::type format_type;
	format_type format_raw; // to read
	int switch_;
	file_import.read((char*)&switch_, sizeof(format_type));
	if (switch_ == -1) {
		file_import.read((char*)&format_raw, sizeof(format_type));
		return static_cast<ConfFileFormat>(format_raw);
	} else {
		return ConfFileFormat::bin_format_base_old; // may also be 1, but will be determined later
	}
}


inline ConfFileFormat getTxtConfigurationFileFormat(const std::string& filename_import_positions)
{
	checkInFile(filename_import_positions);
	std::ifstream file_import;
	file_import.open(filename_import_positions.c_str());

	std::string header_imported_configulation[2];
	getline(file_import, header_imported_configulation[0]);
	getline(file_import, header_imported_configulation[1]);

	auto meta_data = getConfMetaData(header_imported_configulation[0], header_imported_configulation[1]);
	if (meta_data.find("format") != meta_data.end()) {
		return static_cast<ConfFileFormat>(atoi(meta_data["format"].c_str()));
	} else if (meta_data.find("radius_in") != meta_data.end()) {
		return ConfFileFormat::txt_format_circular_couette;
	} else {
		return ConfFileFormat::txt_format_base_old;
	}
}

inline struct base_configuration readBinaryBaseConfiguration(const std::string& filename)
{
	checkInFile(filename);
	auto format = getBinaryConfigurationFileFormat(filename);
	if (format != ConfFileFormat::bin_format_base_old && format != ConfFileFormat::bin_format_base_new) {
		throw std::runtime_error("readBaseConfiguration(): got incorrect binary format.");
	}
	std::ifstream input(filename.c_str(), std::ios::binary | std::ios::in);
	struct base_configuration c;
	c.lees_edwards_disp.reset();

	int np, _switch;
	typedef std::underlying_type<ConfFileFormat>::type format_type;
	format_type fmt;
	if (format != ConfFileFormat::bin_format_base_old){
		input.read((char*)&_switch, sizeof(int));
		input.read((char*)&fmt, sizeof(format_type));
	}
	input.read((char*)&np, sizeof(int));
	input.read((char*)&c.volume_or_area_fraction, sizeof(double));
	input.read((char*)&c.lx, sizeof(double));
	input.read((char*)&c.ly, sizeof(double));
	input.read((char*)&c.lz, sizeof(double));
	input.read((char*)&c.lees_edwards_disp.x, sizeof(double));
	input.read((char*)&c.lees_edwards_disp.y, sizeof(double));

	std::tie(c.position, c.radius) = readPositionsBStream(input, np);
	if (c.ly == 0) { //2d
		c.angle.resize(c.position.size(), 0);
	}

	c.contact_states = readContactStatesBStream(input, np);
	return c;
}

inline struct fixed_velo_configuration readBinaryFixedVeloConfiguration(const std::string& filename)
{
	checkInFile(filename);
	auto format = getBinaryConfigurationFileFormat(filename);
	if (format != ConfFileFormat::bin_format_fixed_vel) {
		throw std::runtime_error("readBinaryFixedVeloConfiguration(): got incorrect binary format.");
	}
	std::ifstream input(filename.c_str(), std::ios::binary | std::ios::in);

	struct fixed_velo_configuration c;

	c.lees_edwards_disp.reset();
	int np, np_fixed;
	input.read((char*)&np, sizeof(int));
	input.read((char*)&np_fixed, sizeof(int));
	input.read((char*)&c.volume_or_area_fraction, sizeof(double));
	input.read((char*)&c.lx, sizeof(double));
	input.read((char*)&c.ly, sizeof(double));
	input.read((char*)&c.lz, sizeof(double));
	input.read((char*)&c.lees_edwards_disp.x, sizeof(double));
	input.read((char*)&c.lees_edwards_disp.y, sizeof(double));

	std::tie(c.position, c.radius) = readPositionsBStream(input, np);
	if (c.ly == 0) { //2d
		c.angle.resize(c.position.size(), 0);
	}

	double vx_, vy_, vz_;
	for (int i=0; i<np_fixed; i++) {
		input.read((char*)&vx_, sizeof(double));
		input.read((char*)&vy_, sizeof(double));
		input.read((char*)&vz_, sizeof(double));
		c.fixed_velocities.push_back(vec3d(vx_, vy_, vz_));
	}

	c.contact_states = readContactStatesBStream(input, np);
	return c;
}

template<typename T>
inline void setMetadataBase(std::istream &input, T &conf)
{

	std::string header_imported_configulation[2];
	getline(input, header_imported_configulation[0]);
	getline(input, header_imported_configulation[1]);

	auto meta_data = getConfMetaData(header_imported_configulation[0], header_imported_configulation[1]);
	std::string key, def;
	key = "lx";
	conf.lx = atof(getMetaParameter(meta_data, key).c_str());
	key = "ly";
	conf.ly = atof(getMetaParameter(meta_data, key).c_str());
	key = "lz";
	conf.lz = atof(getMetaParameter(meta_data, key).c_str());
	key = "vf";
	conf.volume_or_area_fraction = atof(getMetaParameter(meta_data, key).c_str());
	conf.lees_edwards_disp.reset();
	key = "dispx";
	def = "0";
	conf.lees_edwards_disp.x = atof(getMetaParameter(meta_data, key, def).c_str());
	key = "dispy";
	def = "0";
	conf.lees_edwards_disp.y = atof(getMetaParameter(meta_data, key, def).c_str());
}

inline struct base_configuration readTxtBaseConfiguration(const std::string& filename)
{
	/**
		\brief Import a text file base configuration.

		File format:
		# header

		x y z radius
		...
	 */
	checkInFile(filename);
	std::ifstream input(filename.c_str(), std::ios::in);

	struct base_configuration c;
	setMetadataBase(input, c);

	double x_, y_, z_, a_;
	while (input >> x_ >> y_ >> z_ >> a_) {
		c.position.push_back(vec3d(x_, y_, z_));
		c.radius.push_back(a_);
	}
	if (c.ly == 0) { //2d
		c.angle.resize(c.position.size(), 0);
	}
	return c;
}

inline struct fixed_velo_configuration readTxtFixedVeloConfiguration(const std::string& filename)
{
	/**
		\brief Import a text file configuration with imposed velocity particles.

		File format:
		# header

		x y z radius
		...
		x y z radius vx vy vz
		...
	 */
	// http://stackoverflow.com/questions/743191/how-to-parse-lines-with-differing-number-of-fields-in-c

	checkInFile(filename);
	std::ifstream input(filename.c_str(), std::ios::in);

	struct fixed_velo_configuration c;
	setMetadataBase(input, c);

	double x_, y_, z_, a_, vx_, vy_, vz_;
	std::string line;
	while (getline(input, line)) {
		std::istringstream is;
		is.str(line);
		if (!(is >> x_ >> y_ >> z_ >> a_ >> vx_ >> vy_ >> vz_)) {
			is.str(line);
			is >> x_ >> y_ >> z_ >> a_;
			c.position.push_back(vec3d(x_, y_, z_));
			c.radius.push_back(a_);
		} else {
			c.position.push_back(vec3d(x_, y_, z_));
			c.radius.push_back(a_);
			c.fixed_velocities.push_back(vec3d(vx_, vy_, vz_));
		}
	}
	if (c.ly == 0) { //2d
		c.angle.resize(c.position.size(), 0);
	}
	return c;
}

inline void setMetadataCircularCouette(std::istream &file_import, struct circular_couette_configuration &conf)
{
	std::string header_imported_configulation[2];
	getline(file_import, header_imported_configulation[0]);
	getline(file_import, header_imported_configulation[1]);

	auto meta_data = getConfMetaData(header_imported_configulation[0], header_imported_configulation[1]);
	std::string key, def;
	key = "lx";
	conf.lx = atof(getMetaParameter(meta_data, key).c_str());
	key = "ly";
	conf.ly = atof(getMetaParameter(meta_data, key).c_str());
	key = "lz";
	conf.lz = atof(getMetaParameter(meta_data, key).c_str());
	key = "vf";
	conf.volume_or_area_fraction = atof(getMetaParameter(meta_data, key).c_str());
	conf.lees_edwards_disp.reset();
	key = "dispx";
	def = "0";
	conf.lees_edwards_disp.x = atof(getMetaParameter(meta_data, key, def).c_str());
	key = "dispy";
	def = "0";
	conf.lees_edwards_disp.y = atof(getMetaParameter(meta_data, key, def).c_str());
	key = "np_wall1";
	def = "-1";
	conf.np_wall1 = atoi(getMetaParameter(meta_data, key, def).c_str());
	key = "np_wall2";
	def = "-1";
	conf.np_wall2 = atoi(getMetaParameter(meta_data, key, def).c_str());
	key = "radius_in";
	def = "0";
	conf.radius_in = atof(getMetaParameter(meta_data, key, def).c_str());
	key = "radius_out";
	def = "0";
	conf.radius_out = atof(getMetaParameter(meta_data, key, def).c_str());
}

inline struct circular_couette_configuration readTxtCircularCouetteConfiguration(const std::string& filename)
{
	/**
		\brief Import a text file circular Couette configuration.

		File format:
		# header

		x y z radius
		...
	 */

	checkInFile(filename);
 	std::ifstream input(filename.c_str(), std::ios::in);

	struct circular_couette_configuration c;
	setMetadataCircularCouette(input, c);

	double x_, y_, z_, a_;
	while (input >> x_ >> y_ >> z_ >> a_) {
		c.position.push_back(vec3d(x_, y_, z_));
		c.radius.push_back(a_);
	}
	if (c.ly == 0) { //2d
		c.angle.resize(c.position.size(), 0);
	}
	return c;
}

#endif /* defined(__LF_DEM__Configuration__) */
