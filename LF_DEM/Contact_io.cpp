#include "Contact.h"

namespace Contact_ios {

template<typename T> 
struct contact_state readStateBStream(std::istream &input)
{
	T p0, p1;
	double dt_x, dt_y, dt_z, dr_x, dr_y, dr_z;
	input.read((char*)&p0, sizeof(T));
	input.read((char*)&p1, sizeof(T));
	input.read((char*)&dt_x, sizeof(decltype(dt_x)));
	input.read((char*)&dt_y, sizeof(decltype(dt_y)));
	input.read((char*)&dt_z, sizeof(decltype(dt_z)));
	input.read((char*)&dr_x, sizeof(decltype(dr_x)));
	input.read((char*)&dr_y, sizeof(decltype(dr_y)));
	input.read((char*)&dr_z, sizeof(decltype(dr_z)));
	struct contact_state cs;
	cs.p0 = (int)p0;
	cs.p1 = (int)p1;
	cs.disp_tan = vec3d(dt_x, dt_y, dt_z);
	cs.disp_rolling = vec3d(dr_x, dr_y, dr_z);
	return cs;
}

std::vector <struct contact_state> readStatesBStream(std::istream &input, unsigned np)
{
	unsigned ncont;
	input.read((char*)&ncont, sizeof(decltype(ncont)));

	// hacky thing to guess if this is an old format with particle numbers as unsigned short
	bool ushort_format = false;
	unsigned p0;
	std::iostream::pos_type file_pos = input.tellg();
	input.read((char*)&p0, sizeof(decltype(p0)));
	if(p0>np){
		ushort_format = true;
	}
	input.seekg(file_pos);

	std::vector <struct contact_state> cont_states;
	if (ushort_format) {
		for (unsigned i=0; i<ncont; i++) {
			cont_states.push_back(readStateBStream<unsigned short>(input));
		}
	} else {
		for (unsigned i=0; i<ncont; i++) {
			cont_states.push_back(readStateBStream<unsigned>(input));
		}
	}
	return cont_states;
}


void writeStatesBStream(std::ostream &conf_export, const std::vector <struct contact_state> &cs)
{
	unsigned ncont = cs.size();
	conf_export.write((char*)&ncont, sizeof(unsigned int));
	for (int i=0; i<ncont; i++) {
		conf_export.write((char*)&(cs[i].p0), sizeof(unsigned int));
		conf_export.write((char*)&(cs[i].p1), sizeof(unsigned int));
		conf_export.write((char*)&(cs[i].disp_tan.x), sizeof(double));
		conf_export.write((char*)&(cs[i].disp_tan.y), sizeof(double));
		conf_export.write((char*)&(cs[i].disp_tan.z), sizeof(double));
		conf_export.write((char*)&(cs[i].disp_rolling.x), sizeof(double));
		conf_export.write((char*)&(cs[i].disp_rolling.y), sizeof(double));
		conf_export.write((char*)&(cs[i].disp_rolling.z), sizeof(double));
	}
}
} // namespace Contact_io
