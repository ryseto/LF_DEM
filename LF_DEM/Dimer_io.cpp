#include "DimerState.h"
#include "DimerManager.h"
#include "Dimer_io.h"
#include "CheckFile.h"
#include "vec3d.h"

namespace Interactions {

namespace Dimer {

namespace io 
{

std::vector<UnloadedDimerState> readTxtDimer(const std::string& filename)
{
	/**
		\brief Import a text file base configuration.

		File format:
		# header

		i j relaxed_length
		...
	 */
	checkInFile(filename);
	std::ifstream input(filename.c_str(), std::ios::in);

	std::vector<UnloadedDimerState> c;
	UnloadedDimerState st;
	while (input >> st.p0 >> st.p1 >> st.relaxed_length) {		
		c.push_back(st);
	}
	return c;
}

std::vector <struct DimerState> readStatesBStream(std::istream &input)
{
	unsigned ndimer;
	input.read((char*)&ndimer, sizeof(decltype(ndimer)));
	std::vector <struct DimerState> states;
	for (unsigned i=0; i<ndimer; i++) {
		unsigned p0, p1;
		double dt_x, dt_y, dt_z, dr_x, dr_y, dr_z, relaxed_length;
		input.read((char*)&p0, sizeof(unsigned));
		input.read((char*)&p1, sizeof(unsigned));
		input.read((char*)&relaxed_length, sizeof(decltype(dt_x)));
		input.read((char*)&dt_x, sizeof(decltype(dt_x)));
		input.read((char*)&dt_y, sizeof(decltype(dt_y)));
		input.read((char*)&dt_z, sizeof(decltype(dt_z)));
		input.read((char*)&dr_x, sizeof(decltype(dr_x)));
		input.read((char*)&dr_y, sizeof(decltype(dr_y)));
		input.read((char*)&dr_z, sizeof(decltype(dr_z)));
		struct DimerState s;
		s.p0 = (int)p0;
		s.p1 = (int)p1;
		s.sliding_st.stretch = vec3d(dt_x, dt_y, dt_z);
		s.sliding_st.relaxed_length = relaxed_length;
		s.rotation_st.stretch = vec3d(dr_x, dr_y, dr_z);
		s.rotation_st.relaxed_length = 0;
		states.push_back(s);
	}
	return states;
}


void writeStatesBStream(std::ostream &conf_export, const DimerManager &dm) {
	unsigned ndimer = dm.size();
	conf_export.write((char*)&ndimer, sizeof(unsigned int));
	for (const auto &dimer: dm) {
		auto ds = dimer->getState();
		conf_export.write((char*)&(ds.p0), sizeof(unsigned int));
		conf_export.write((char*)&(ds.p1), sizeof(unsigned int));
		conf_export.write((char*)&(ds.sliding_st.relaxed_length), sizeof(double));
		conf_export.write((char*)&(ds.sliding_st.stretch.x), sizeof(double));
		conf_export.write((char*)&(ds.sliding_st.stretch.y), sizeof(double));
		conf_export.write((char*)&(ds.sliding_st.stretch.z), sizeof(double));
		conf_export.write((char*)&(ds.rotation_st.stretch.x), sizeof(double));
		conf_export.write((char*)&(ds.rotation_st.stretch.y), sizeof(double));
		conf_export.write((char*)&(ds.rotation_st.stretch.z), sizeof(double));
	}
}

} //namespace io

} // namespace Dimer

} // namespace Interactions
