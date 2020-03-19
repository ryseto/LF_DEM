#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <type_traits>
#include "vec3d.h"
#include "Contact.h"
#include "TimeActivatedAdhesion_io.h"
#include "TimeActivatedAdhesion.h"
#include "ActivatedAdhesionState.h"
#include "DimerState.h"

#ifndef __LF_DEM__Configuration__
#define __LF_DEM__Configuration__

class System;

enum struct ConfFileFormat : int { // assign values as it is for output and otherwise may be compiler dependent.
	// Deprecated formats are not output anymore by LF_DEM, but should still be readable
	bin_format_base_old = 1,			// deprecated
	bin_format_base_new = 2,			// deprecated
	bin_format_fixed_vel = 3,			// deprecated
	txt_format_base_old = 4,
	txt_format_base_new = 5,
	txt_format_fixed_vel = 6,
	txt_format_circular_couette = 7,
	bin_format_base_shear = 8,
	bin_format_fixed_vel_shear = 9,
	bin_delayed_adhesion = 10,
	bin_dimers = 11,
	bin_activated_adhesion = 12
};

struct base_configuration {
	double lx, ly, lz;
	double volume_or_area_fraction;
	
	std::vector <vec3d> position;
	std::vector <double> radius;
	std::vector <double> angle;

	std::vector <struct Interactions::contact_state> contact_states;
};

struct base_shear_configuration {
	double lx, ly, lz;
	double volume_or_area_fraction;
	vec3d lees_edwards_disp;

	std::vector <vec3d> position;
	std::vector <double> radius;
	std::vector <double> angle;

	std::vector <struct Interactions::contact_state> contact_states;
};

struct fixed_velo_configuration {
	double lx, ly, lz;
	double volume_or_area_fraction;
	vec3d lees_edwards_disp;
	int np_wall1;
	int np_wall2;
	double z_bot;
	double z_top;
	std::vector <vec3d> position;
	std::vector <double> radius;
	std::vector <double> angle;
	std::vector<vec3d> fixed_velocities;

	std::vector <struct Interactions::contact_state> contact_states;
};

struct delayed_adhesion_configuration {
	struct base_configuration base;
	vec3d lees_edwards_disp;
	std::vector <struct Interactions::TActAdhesion::State> adhesion_states;
};

struct activated_adhesion_configuration {
	struct base_configuration base;
	vec3d lees_edwards_disp;
	std::vector <struct Interactions::ActAdhesion::State> adhesion_states;
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

	std::vector <struct Interactions::contact_state> contact_states;
};

ConfFileFormat getBinaryConfigurationFileFormat(const std::string& filename_import_positions);

ConfFileFormat getTxtConfigurationFileFormat(const std::string& filename_import_positions);

struct base_shear_configuration readBinaryBaseShearConfiguration(const std::string& filename);

struct base_configuration readBinaryBaseConfiguration(std::ifstream &input);

struct delayed_adhesion_configuration readBinaryDelayedAdhesionConfiguration(std::string filename);

struct activated_adhesion_configuration readBinaryActivatedAdhesionConfiguration(std::string filename);

struct fixed_velo_configuration readBinaryFixedVeloConfiguration(const std::string& filename);

std::pair<struct base_shear_configuration, std::vector<Interactions::Dimer::DimerState>> readBinaryDimerConfiguration(const std::string &filename);

struct base_shear_configuration readTxtBaseConfiguration(const std::string& filename);

struct fixed_velo_configuration readTxtFixedVeloConfiguration(const std::string& filename);

struct circular_couette_configuration readTxtCircularCouetteConfiguration(const std::string& filename);

void outputBinaryConfiguration(const System &sys, 
							   std::string conf_filename,
							   ConfFileFormat format);
#endif /* defined(__LF_DEM__Configuration__) */
