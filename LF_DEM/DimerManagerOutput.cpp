#include "DimerManager.h"
#include "ParameterSet.h"
#include "OutputData.h"
#include "PairwiseConfig.h"


namespace Interactions {
namespace Dimer {
void output(const DimerManager &manager, ParticleVelocity *vel, const struct Parameters::ParameterSet &p, OutputData &outdata_int)
{
	struct PairVelocity pvel;
	manager.pdist->setVelocityState(vel);
	for (const auto &inter: manager.interactions) {
		unsigned int i, j;
		std::tie(i, j) = inter->get_par_num();
		manager.pdist->getVelocities(i, j, pvel);

		outdata_int.entryData("particle 1 label", Dimensional::Dimension::none, 1, i);
		outdata_int.entryData("particle 2 label", Dimensional::Dimension::none, 1, j);
		outdata_int.entryData("normal vector, oriented from particle 1 to particle 2", \
							  Dimensional::Dimension::none, 3, inter->getUnitSeparation());
		outdata_int.entryData("Separation", \
							  Dimensional::Dimension::none, 1,  inter->getSeparation().norm());
		outdata_int.entryData("Relaxed separation", \
							  Dimensional::Dimension::none, 1,  inter->getRelaxedLength());
	
		outdata_int.entryData("total force", Dimensional::Dimension::Force, 3, std::get<0>(inter->getForceTorque(pvel)));
		outdata_int.entryData("dashpot force", Dimensional::Dimension::Force, 3, std::get<0>(inter->getForceTorqueDashpot(pvel)));
		outdata_int.entryData("spring force", Dimensional::Dimension::Force, 3, std::get<0>(inter->getForceTorqueSpring()));
	}
}
}
} // namespace Interactions