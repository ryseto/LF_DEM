#include <string>

#ifndef __LF_DEM__States__
#define __LF_DEM__States__

class System;

namespace State {

enum class StateFileFormat : int {
	basic = 1
};

struct Clock {
	double time_; ///< time elapsed since beginning of the time evolution.
	double time_in_simulation_units; ///< time elapsed since beginning of the time evolution. \b note: this is measured in Simulation (output) units, not in internal System units.
	double cumulated_strain;
};

struct BasicCheckpoint {
	struct Clock clock;
};

static BasicCheckpoint zero_time_basicchkp = {{0, 0, 0}};
void outputStateBinary(std::string state_filename, const System &sys);
struct BasicCheckpoint readBasicCheckpoint(const std::string &filename);
bool isZeroTimeChkp(BasicCheckpoint chkp);
}
#endif // #ifndef __LF_DEM__States__
