#include <stdexcept>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include "global.h"
#include "Contact.h"
#include "vec3d.h"

#ifndef __LF_DEM__Configuration__
#define __LF_DEM__Configuration__

struct base_configuration {
  int np;
  double lx, ly, lz;
  vec3d initial_lees_edwards_disp;

  std::vector <vec3d> initial_position;
  std::vector <double> radius;

  std::vector <struct contact_state> contact_states;
};

struct fixed_velo_configuration {
  int np;
  int np_fixed;
  double lx, ly, lz;
  vec3d initial_lees_edwards_disp;

  std::vector <vec3d> initial_position;
  std::vector <double> radius;
  std::vector<vec3d> fixed_velocities;

  std::vector <struct contact_state> contact_states;
};

struct circular_couette_configuration {
  int np;
  int np_fixed;
  int np_wall1;
	int np_wall2;
	double radius_in;
	double radius_out;

  double lx, ly, lz;
  vec3d initial_lees_edwards_disp;

  std::vector <vec3d> initial_position;
  std::vector <double> radius;
  std::vector<vec3d> fixed_velocities;

  std::vector <struct contact_state> contact_states;
};

inline std::vector <struct contact_state> readContactStatesBStream(std::istream &input, unsigned int np) {
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
    if((int)p0>np){
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

inline std::pair<std::vector <vec3d>, std::vector <double>> readPositionsBStream(std::istream &input, unsigned int np) {
  double x_, y_, z_, r_;
  std::vector <vec3d> position;
  std::vector <double> radius;
  for (int i=0; i<np; i++) {
    input.read((char*)&x_, sizeof(double));
    input.read((char*)&y_, sizeof(double));
    input.read((char*)&z_, sizeof(double));
    input.read((char*)&r_, sizeof(double));
    position.push_back(vec3d(x_,y_,z_));
    radius.push_back(r_);
  }
  return std::make_pair(position, radius);
}


inline int getBinaryConfigurationFileFormat(const std::string& filename_import_positions) {
  std::ifstream file_import;
  file_import.open(filename_import_positions.c_str(), std::ios::binary | std::ios::in);
  if (!file_import) {
    std::ostringstream error_str;
    error_str  << " Position file '" << filename_import_positions << "' not found." << std::endl;
    throw std::runtime_error(error_str.str());
  }
  int binary_format_version;
  int switch_;
  int np_fixed;
  double lx, ly, lz;
  file_import.read((char*)&switch_, sizeof(int));
  if (switch_ == -1) {
    file_import.read((char*)&binary_format_version, sizeof(int));
  } else {
    binary_format_version = BIN_FORMAT_BASE_NEW; // may also be 1, but will be determined later
  }
  file_import.close();

  return binary_format_version;
}

inline struct base_configuration readBaseConfiguration(std::istream &input) {
  struct base_configuration c;
  double volume_or_area_fraction; // now unused
  c.initial_lees_edwards_disp.reset();
  input.read((char*)&c.np, sizeof(int));

  input.read((char*)&volume_or_area_fraction, sizeof(double));
  input.read((char*)&c.lx, sizeof(double));
  input.read((char*)&c.ly, sizeof(double));
  input.read((char*)&c.lz, sizeof(double));
  input.read((char*)&c.initial_lees_edwards_disp.x, sizeof(double));
  input.read((char*)&c.initial_lees_edwards_disp.y, sizeof(double));

  std::tie(c.initial_position, c.radius) = readPositionsBStream(input, c.np);
  c.contact_states = readContactStatesBStream(input, c.np);
  return c;
}

inline struct fixed_velo_configuration readFixedVeloConfiguration(std::istream &input) {
  struct fixed_velo_configuration c;
  double volume_or_area_fraction; // now unused

  c.initial_lees_edwards_disp.reset();
  int switch_, version_format_;
  input.read((char*)&switch_, sizeof(int));
  input.read((char*)&version_format_, sizeof(int));
  if(switch_ != -1 || version_format_ != 3) {
    throw std::runtime_error("readFixedVeloConfiguration(): got incorrect binary format.");
  }

  input.read((char*)&c.np, sizeof(int));
  input.read((char*)&c.np_fixed, sizeof(int));
  input.read((char*)&volume_or_area_fraction, sizeof(double));
  input.read((char*)&c.lx, sizeof(double));
  input.read((char*)&c.ly, sizeof(double));
  input.read((char*)&c.lz, sizeof(double));
  input.read((char*)&c.initial_lees_edwards_disp.x, sizeof(double));
  input.read((char*)&c.initial_lees_edwards_disp.y, sizeof(double));

  std::tie(c.initial_position, c.radius) = readPositionsBStream(input, c.np);

  double vx_, vy_, vz_;
  for (int i=0; i<c.np_fixed; i++) {
    input.read((char*)&vx_, sizeof(double));
    input.read((char*)&vy_, sizeof(double));
    input.read((char*)&vz_, sizeof(double));
    c.fixed_velocities.push_back(vec3d(vx_, vy_, vz_));
  }

  c.contact_states = readContactStatesBStream(input, c.np);
  return c;
}


#endif /* defined(__LF_DEM__Configuration__) */
