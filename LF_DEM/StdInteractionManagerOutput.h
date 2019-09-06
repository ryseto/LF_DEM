#ifndef __LF_DEM__StdInteractionManagerOutput__
#define __LF_DEM__StdInteractionManagerOutput__

namespace Parameters {
struct ParameterSet;
}

class OutputData;
class ParticleVelocity;

namespace Interactions {
class StdInteractionManager;

class StdInteractionManagerOutput {
public:
	// StdInteractionManagerOutput();
	void output(const StdInteractionManager &manager, ParticleVelocity *vel, const struct Parameters::ParameterSet &p, OutputData &outdata_int);
};

}

#endif
