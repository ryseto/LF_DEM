#ifndef __LF_DEM__InteractionManagerOutput__
#define __LF_DEM__InteractionManagerOutput__

class OutputData;
class ParticleVelocity;
namespace Parameters {
class ParameterSet;
}

namespace Interactions {
class StdInteractionManager;
void output(const StdInteractionManager &manager, ParticleVelocity *vel, const struct Parameters::ParameterSet &p, OutputData &outdata_int);

namespace Dimer {
class DimerManager;
void output(const DimerManager &manager, ParticleVelocity *vel, const struct Parameters::ParameterSet &p, OutputData &outdata_int);
}

}

#endif
