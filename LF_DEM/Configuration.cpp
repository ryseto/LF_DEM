#include <stdexcept>
#include <string>
#include <iostream>
#include <sstream>
#include <map>
#include <set>
#include <algorithm>
#include "Contact.h"
#include "Configuration.h"
#include "System.h"
#include "LeesEdwards.h"
#include "CheckFile.h"
#include "Dimer_io.h"

inline void removeBlank(std::string& str)
{
	str.erase(std::remove_if(str.begin(), str.end(), (int(*)(int))isspace), str.end());
}

std::vector<std::string> splitString(const std::string& str)
{
	std::string stripped_str = str;
	std::size_t brk = 0;

	std::vector<std::string> elements;

	while (stripped_str.length() > 0) {
		brk = stripped_str.find(" ");
		std::string first_part;
		if (brk < std::string::npos) {
			brk += 1;
			first_part = stripped_str.substr(0, brk);
			stripped_str = stripped_str.substr(brk, std::string::npos);
		} else {
			first_part = stripped_str.substr(0, brk);
			stripped_str = "";
		}
		removeBlank(first_part);
		if ( first_part.length()>0 ) {
			elements.push_back(first_part);
		}
	}
	return elements;
}

std::pair<std::vector <vec3d>, std::vector <double>> readPositionsBStream(std::istream &input, unsigned int np)
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

std::string getMetaParameter(std::map<std::string,std::string> &meta_data, std::string &key)
{
	if (meta_data.find(key) != meta_data.end()) {
		return meta_data[key];
	} else {
		std::ostringstream error_str;
		error_str  << " Simulation:: parameter '" << key << "' not found in the header of the configuration file." << std::endl;
		throw std::runtime_error(error_str.str());
	}
}

std::string getMetaParameter(std::map<std::string,std::string> &meta_data,
							 std::string &key,
							 const std::string &default_val)
{
	if (meta_data.find(key) != meta_data.end()) {
		return meta_data[key];
	} else {
		return default_val;
	}
}

std::map<std::string,std::string> getConfMetaData(const std::string &line1, const std::string &line2)
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

ConfFileFormat getBinaryConfigurationFileFormat(const std::string& filename_import_positions)
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

ConfFileFormat getTxtConfigurationFileFormat(const std::string& filename_import_positions)
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
	} else if (meta_data.find("z_bot") != meta_data.end()) {
		return ConfFileFormat::txt_format_fixed_vel;
	} else if (meta_data.find("np_wall1") != meta_data.end()) {
		return ConfFileFormat::txt_format_fixed_vel;
	} else {
		return ConfFileFormat::txt_format_base_old;
	}
}

template<typename T>
void fillBaseConfiguration(T &vel_conf,  const struct base_configuration &base) 
{
	vel_conf.lx = base.lx;
	vel_conf.ly = base.ly;
	vel_conf.lz = base.lz;
	vel_conf.volume_or_area_fraction = base.volume_or_area_fraction;
	vel_conf.position = base.position;
	vel_conf.radius = base.radius;
	vel_conf.angle = base.angle;
	vel_conf.contact_states = base.contact_states;
}

struct base_shear_configuration readBinaryBaseShearConfiguration_old(const std::string& filename)
{
	checkInFile(filename);
	auto format = getBinaryConfigurationFileFormat(filename);
	if (format != ConfFileFormat::bin_format_base_old && format != ConfFileFormat::bin_format_base_new) {
		throw std::runtime_error("readBinaryBaseShearConfiguration(): got incorrect binary format.");
	}
	std::ifstream input(filename.c_str(), std::ios::binary | std::ios::in);
	struct base_shear_configuration c;
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

	c.contact_states = Interactions::Contact_ios::readStatesBStream(input, np);
	return c;
}

struct base_shear_configuration readBinaryBaseShearConfiguration(const std::string &filename)
{
	checkInFile(filename);
	auto format = getBinaryConfigurationFileFormat(filename);
	if (format == ConfFileFormat::bin_format_base_old && format != ConfFileFormat::bin_format_base_new) {
		return readBinaryBaseShearConfiguration_old(filename);
	}
	std::set<ConfFileFormat> allowed_formats = \
	{
		ConfFileFormat::bin_format_base_shear,
	};
	if (!allowed_formats.count(format)) {
		throw std::runtime_error("readBinaryBaseShearConfiguration(): got incorrect binary format.");
	}
	std::ifstream input(filename.c_str(), std::ios::binary | std::ios::in);
	struct base_shear_configuration c;

	int _switch;
	typedef std::underlying_type<ConfFileFormat>::type format_type;
	format_type fmt;
	input.read((char*)&_switch, sizeof(int));
	input.read((char*)&fmt, sizeof(format_type));

	auto base = readBinaryBaseConfiguration(input);
	fillBaseConfiguration<struct base_shear_configuration>(c, base);

	c.lees_edwards_disp.reset();
	input.read((char*)&c.lees_edwards_disp.x, sizeof(double));
	input.read((char*)&c.lees_edwards_disp.y, sizeof(double));

	return c;
}

struct base_configuration readBinaryBaseConfiguration(std::ifstream &input)
{
	// When modifying this function, don't forget to make the equivalent modification for the writeBinaryBaseConfiguration function.
	struct base_configuration c;
	unsigned np;
	input.read((char*)&np, sizeof(unsigned));
	input.read((char*)&c.volume_or_area_fraction, sizeof(double));
	input.read((char*)&c.lx, sizeof(double));
	input.read((char*)&c.ly, sizeof(double));
	input.read((char*)&c.lz, sizeof(double));

	double x_, y_, z_, r_, a_;
	for (unsigned i=0; i<np; i++) {
		input.read((char*)&x_, sizeof(double));
		input.read((char*)&y_, sizeof(double));
		input.read((char*)&z_, sizeof(double));
		input.read((char*)&r_, sizeof(double));
		c.position.push_back(vec3d(x_,y_,z_));
		c.radius.push_back(r_);
		if (c.ly == 0) {
			input.read((char*)&a_, sizeof(double)); // angle
			c.angle.push_back(a_);
		}
	}
	c.contact_states = Interactions::Contact_ios::readStatesBStream(input, np);
	return c;
}

void writeBinaryBaseConfiguration(std::ofstream &conf_export, const struct base_configuration &conf) 
{
	// When modifying this function, don't forget to make the equivalent modification for the readBinaryBaseConfiguration function.
	std::vector<std::vector<double>> pos(conf.position.size());
	unsigned dims = 4;
	if (conf.ly == 0) {
		dims = 5;
	}
	for (unsigned i=0; i<conf.position.size(); i++) {
		pos[i].resize(dims);
		pos[i][0] = conf.position[i].x;
		pos[i][1] = conf.position[i].y;
		pos[i][2] = conf.position[i].z;
		pos[i][3] = conf.radius[i];
		if (conf.ly == 0) {
			pos[i][4] = conf.angle[i];
		}
	}
	unsigned np = conf.position.size();
	conf_export.write((char*)&np, sizeof(unsigned));
	conf_export.write((char*)&conf.volume_or_area_fraction, sizeof(double));
	conf_export.write((char*)&conf.lx, sizeof(double));
	conf_export.write((char*)&conf.ly, sizeof(double));
	conf_export.write((char*)&conf.lz, sizeof(double));
	for (unsigned i=0; i<conf.position.size(); i++) {
		conf_export.write((char*)&pos[i][0], dims*sizeof(double));
	}
	Interactions::Contact_ios::writeStatesBStream(conf_export, conf.contact_states);
}

std::pair<struct base_shear_configuration, std::vector<Interactions::Dimer::DimerState>> readBinaryDimerConfiguration(const std::string &filename)
{

	checkInFile(filename);
	auto format = getBinaryConfigurationFileFormat(filename);
	std::set<ConfFileFormat> allowed_formats = \
	{
		ConfFileFormat::bin_dimers,
	};
	if (!allowed_formats.count(format)) {
		throw std::runtime_error("readBinaryDimerConfiguration(): got incorrect binary format.");
	}
	std::ifstream input(filename.c_str(), std::ios::binary | std::ios::in);
	struct base_shear_configuration c;

	int _switch;
	typedef std::underlying_type<ConfFileFormat>::type format_type;
	format_type fmt;
	input.read((char*)&_switch, sizeof(int));
	input.read((char*)&fmt, sizeof(format_type));

	auto base = readBinaryBaseConfiguration(input);
	fillBaseConfiguration<struct base_shear_configuration>(c, base);
	auto dimers = Interactions::Dimer::io::readStatesBStream(input);

	c.lees_edwards_disp.reset();
	input.read((char*)&c.lees_edwards_disp.x, sizeof(double));
	input.read((char*)&c.lees_edwards_disp.y, sizeof(double));

	return std::make_pair(c, dimers);
}



std::vector<vec3d> readBinaryFixedVelocities(std::ifstream &input)
{
	std::vector<vec3d> fixed_velocities;
	unsigned np_fixed;
	input.read((char*)&np_fixed, sizeof(unsigned));
	double vx_, vy_, vz_;
	for (unsigned i=0; i<np_fixed; i++) {
		input.read((char*)&vx_, sizeof(double));
		input.read((char*)&vy_, sizeof(double));
		input.read((char*)&vz_, sizeof(double));
		fixed_velocities.push_back(vec3d(vx_, vy_, vz_));
	}
	return fixed_velocities;
}

void writeBinaryFixedVelocities(std::ofstream &conf_export, const std::vector<vec3d> &fixed_velocities)
{
	unsigned np_fixed = fixed_velocities.size();
	conf_export.write((char*)&np_fixed, sizeof(unsigned));
	std::vector<std::vector<double> > vel(np_fixed);
	for (unsigned i=0; i<np_fixed; i++) {
		vel[i].resize(3);
		vel[i][0] = fixed_velocities[i].x;
		vel[i][1] = fixed_velocities[i].y;
		vel[i][2] = fixed_velocities[i].z;
	}
	for (unsigned i=0; i<np_fixed; i++) {
		conf_export.write((char*)&vel[i][0], 3*sizeof(double));
	}
}

struct delayed_adhesion_configuration readBinaryDelayedAdhesionConfiguration(std::string filename)
{
	checkInFile(filename);
	auto format = getBinaryConfigurationFileFormat(filename);
	std::set<ConfFileFormat> allowed_formats = \
	{
		ConfFileFormat::bin_delayed_adhesion,
	};
	if (!allowed_formats.count(format)) {
		throw std::runtime_error("readBinaryDelayedAdhesionConfiguration(): got incorrect binary format.");
	}
	std::ifstream input(filename.c_str(), std::ios::binary | std::ios::in);
	struct delayed_adhesion_configuration c;

	int _switch;
	typedef std::underlying_type<ConfFileFormat>::type format_type;
	format_type fmt;
	input.read((char*)&_switch, sizeof(int));
	input.read((char*)&fmt, sizeof(format_type));

	c.base = readBinaryBaseConfiguration(input);
	c.adhesion_states = Interactions::TActAdhesion::readStatesBStream(input);
	c.lees_edwards_disp.reset();
	input.read((char*)&c.lees_edwards_disp.x, sizeof(double));
	input.read((char*)&c.lees_edwards_disp.y, sizeof(double));
	
	return c;
}

struct fixed_velo_configuration readBinaryFixedVeloConfigurationOld(const std::string& filename)
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

	c.contact_states = Interactions::Contact_ios::readStatesBStream(input, np);
	return c;
}

struct fixed_velo_configuration readBinaryFixedVeloConfiguration(const std::string& filename)
{
	checkInFile(filename);
	auto format = getBinaryConfigurationFileFormat(filename);
	if (format == ConfFileFormat::bin_format_fixed_vel) {
		return readBinaryFixedVeloConfigurationOld(filename);
	}

	if (format != ConfFileFormat::bin_format_fixed_vel_shear) {
		throw std::runtime_error("readBinaryFixedVeloConfiguration(): got incorrect binary format.");
	}
	std::ifstream input(filename.c_str(), std::ios::binary | std::ios::in);

	struct fixed_velo_configuration c;

	int _switch;
	typedef std::underlying_type<ConfFileFormat>::type format_type;
	format_type fmt;
	input.read((char*)&_switch, sizeof(int));
	input.read((char*)&fmt, sizeof(format_type));

	auto base = readBinaryBaseConfiguration(input);
	fillBaseConfiguration<struct fixed_velo_configuration>(c, base);
	c.fixed_velocities = readBinaryFixedVelocities(input);

	c.lees_edwards_disp.reset();
	input.read((char*)&c.lees_edwards_disp.x, sizeof(double));
	input.read((char*)&c.lees_edwards_disp.y, sizeof(double));

	return c;
}

void writeBinaryHeader(std::ofstream &conf_export, ConfFileFormat format)
{
	int conf_switch = -1; // older formats did not have labels, -1 signs for a labeled binary
	typedef std::underlying_type<ConfFileFormat>::type format_type;
	conf_export.write((char*)&conf_switch, sizeof(int));
	conf_export.write((char*)&format, sizeof(format_type));
}

void outputBinaryConfiguration(const System &sys, 
							   std::string conf_filename,
							   ConfFileFormat format)
{
	/**
	 \brief Saves the current configuration of the system in a binary file.
	 */
	std::set<ConfFileFormat> allowed_formats = \
	{
		ConfFileFormat::bin_delayed_adhesion,
		ConfFileFormat::bin_format_fixed_vel_shear,
		ConfFileFormat::bin_format_base_shear,
		ConfFileFormat::bin_dimers
	};
	if (!allowed_formats.count(format)) {
		throw std::runtime_error("outputBinaryConfiguration(): got incorrect binary format.");
	}

	std::ofstream conf_export;
	conf_export.open(conf_filename.c_str(), std::ios::binary | std::ios::out);
	writeBinaryHeader(conf_export, format);

	writeBinaryBaseConfiguration(conf_export, sys.getBaseConfiguration());
	if (format == ConfFileFormat::bin_format_fixed_vel_shear) {
		if (sys.res_solver->velo_assignor) {
			writeBinaryFixedVelocities(conf_export, sys.res_solver->velo_assignor->getFixedVel().vel);
		}
	} else if (format == ConfFileFormat::bin_delayed_adhesion) {
		Interactions::TActAdhesion::writeStatesBStream(conf_export,
										 			   *(sys.interaction));
	} else if (format == ConfFileFormat::bin_dimers) {
		Interactions::Dimer::io::writeStatesBStream(conf_export, *(sys.dimer_manager));
	}
	auto shear_disp = sys.lees->getShearDisp();
	conf_export.write((char*)&(shear_disp.x), sizeof(double));
	conf_export.write((char*)&(shear_disp.y), sizeof(double));

	conf_export.close();
}

struct base_shear_configuration readTxtBaseConfiguration(const std::string& filename)
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

	struct base_shear_configuration c;
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

void setMetadataFixedVelo(std::istream &file_import, struct fixed_velo_configuration &conf)
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
	key = "np_wall1";
	def = "-1";
	conf.np_wall1 = atof(getMetaParameter(meta_data, key, def).c_str());
	key = "np_wall2";
	def = "-1";
	conf.np_wall2 = atof(getMetaParameter(meta_data, key, def).c_str());
	key = "z_bot";
	key = "0";
	conf.z_bot = atof(getMetaParameter(meta_data, key, def).c_str());
	key = "z_top";
	key = "0";
	conf.z_top = atof(getMetaParameter(meta_data, key, def).c_str());
}

struct fixed_velo_configuration readTxtFixedVeloConfiguration(const std::string& filename)
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
	setMetadataFixedVelo(input, c);

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
			//std::cerr << x_ << ' ' << y_ << ' ' << z_ << std::endl;
		} else {
			c.position.push_back(vec3d(x_, y_, z_));
			c.radius.push_back(a_);
			c.fixed_velocities.push_back(vec3d(vx_, vy_, vz_));
			//std::cerr << vx_ << ' ' <<  vy_ << ' ' <<  vz_ << std::endl;
		}
	}
	if (c.ly == 0) { //2d
		c.angle.resize(c.position.size(), 0);
	}
	return c;
}

void setMetadataCircularCouette(std::istream &file_import, struct circular_couette_configuration &conf)
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

struct circular_couette_configuration readTxtCircularCouetteConfiguration(const std::string& filename)
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
