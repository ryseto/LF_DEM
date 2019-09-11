#include <string>

#ifndef __LF_DEM__States__
#define __LF_DEM__States__

namespace State {
	
enum class StateFileFormat : int 
{
	basic = 1,
};

struct Clock {
	double time_; ///< time elapsed since beginning of the time evolution.
	double cumulated_strain;
};

struct BasicCheckpoint {
	struct Clock clock;
};

void outputStateBinary(std::string state_filename, double strain, double _time);
struct BasicCheckpoint readBasicCheckpoint(const std::string &filename);
static BasicCheckpoint zero_time_basicchkp = {{0, 0}};
bool isZeroTimeChkp(BasicCheckpoint chkp);
}
#endif // #ifndef __LF_DEM__States__
