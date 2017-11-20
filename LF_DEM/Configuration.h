#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <type_traits>
#include "vec3d.h"
#include "Contact.h"
#include "TimeActivatedAdhesion_io.h"
#include "TimeActivatedAdhesion.h"

#ifndef __LF_DEM__Configuration__
#define __LF_DEM__Configuration__


enum struct ConfFileFormat : int { // assign values as it is for output and otherwise may be compiler dependent.
	// Deprecated formats are not output anymore by LF_DEM, but should still be readable
	bin_format_base_old = 1,    		// deprecated
	bin_format_base_new = 2,			// deprecated
	bin_format_fixed_vel = 3,			// deprecated
	txt_format_base_old = 4,
	txt_format_base_new = 5,
	txt_format_fixed_vel = 6,
	txt_format_circular_couette = 7,
	bin_format_base_shear = 8,
	bin_format_fixed_vel_shear = 9,
	bin_delayed_adhesion = 10,
};


struct base_configuration {
	double lx, ly, lz;
	double volume_or_area_fraction;

	std::vector <vec3d> position;
	std::vector <double> radius;
	std::vector <double> angle;

	std::vector <struct contact_state> contact_states;
};

struct base_shear_configuration {
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

struct delayed_adhesion_configuration {
	struct base_configuration base;
	vec3d lees_edwards_disp;
	std::vector <struct TActAdhesion::State> adhesion_states;
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

ConfFileFormat getBinaryConfigurationFileFormat(const std::string& filename_import_positions);

ConfFileFormat getTxtConfigurationFileFormat(const std::string& filename_import_positions);

struct base_shear_configuration readBinaryBaseShearConfiguration(const std::string& filename);

struct base_configuration readBinaryBaseConfiguration(std::ifstream &input);

struct delayed_adhesion_configuration readBinaryDelayedAdhesionConfiguration(std::string filename);

struct fixed_velo_configuration readBinaryFixedVeloConfiguration(const std::string& filename);

struct base_shear_configuration readTxtBaseConfiguration(const std::string& filename);

struct fixed_velo_configuration readTxtFixedVeloConfiguration(const std::string& filename);

struct circular_couette_configuration readTxtCircularCouetteConfiguration(const std::string& filename);

#endif /* defined(__LF_DEM__Configuration__) */
