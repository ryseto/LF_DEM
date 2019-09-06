#include "StdInteractionManagerOutput.h"
#include "StdInteractionManager.h"
#include "ParameterSet.h"
#include "OutputData.h"
#include "PairwiseConfig.h"
#include "ParticleConfig.h"


namespace Interactions {

// StdInteractionManagerOutput::StdInteractionManagerOutput(std::shared_ptr<StdInteractionManager> interactions) :
// interaction(interactions)
// {

// }

void StdInteractionManagerOutput::output(const StdInteractionManager &manager, ParticleVelocity *vel, const struct Parameters::ParameterSet &p, OutputData &outdata_int)
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
		outdata_int.entryData("dimensionless gap = s-2, s = 2r/(a1+a2)", \
							  Dimensional::Dimension::none, 1,  inter->getReducedGap());
		
		/* [NOTE]
		 * Lubrication forces are reference values
		 * in the Brownian case. The force balancing
		 * velocities are recalculated without
		 * including the Brownian forces.
		 * It seems there is no better way to visualize
		 * the lubrication forces.
		 */

		if (Interactions::has_lubrication(p.lub)) {
			double lub_npart = 0;
			vec3d lub_tpart;
			if (inter->lubrication) {
				auto force = inter->lubrication->getTotalForce(pvel, 0.5*(manager.Einf->E[i] + manager.Einf->E[j]));
				lub_npart = - dot(force, inter->getUnitSeparation());
				lub_tpart = force + lub_npart*inter->getUnitSeparation();
			}
			outdata_int.entryData("normal part of the lubrication force (positive for compression)", Dimensional::Dimension::Force, 1, lub_npart);
			outdata_int.entryData("tangential part of the lubrication force", Dimensional::Dimension::Force, 3, lub_tpart);
		}
		/*
		 * Contact forces are the sums of spring forces and dashpot forces.
		 * (It can be negative even repulsive contact force).
		 */
		int cs = 0;
		double c_npart = 0;
		vec3d c_tpart;
		// Sym2Tensor stress_contact;
		if (inter->contact) {
			cs = inter->contact->getFrictionState();
			c_npart = - inter->contact->getNormalForceValue(pvel);
			c_tpart = inter->contact->getTangentialForce(pvel);
			// stress_contact = inter->contact->getContactStressXF();
		}
		outdata_int.entryData("contact state "
							  "(0 = no contact, "
							  "1 = frictionless contact, "
							  "2 = non-sliding frictional, "
							  "3 = sliding frictional)",
							  Dimensional::Dimension::none, 1, cs);
		outdata_int.entryData("norm of the normal part of the contact force", Dimensional::Dimension::Force, 1, c_npart);
		outdata_int.entryData("tangential part of the contact force", Dimensional::Dimension::Force, 3, c_tpart);
		// outdata_int.entryData("Viscosity contribution of contact xF", Dimensional::Dimension::Stress, 1, \
		// 					  doubledot(stress_contact, sys.imposed_flow->sym_grad_u/sr)/sr);

		if (Interactions::has_repulsion(p.repulsion)) {
			double r_npart = 0;
			if (inter->repulsion) {
				r_npart = inter->repulsion->getForceNorm();
			}
			outdata_int.entryData("norm of the normal repulsive force", Dimensional::Dimension::Force, 1, r_npart);
		}

		if (Interactions::has_delayed_adhesion(p.TA_adhesion)) {
			double a_npart = 0;
			double tratio = 0;
			if (inter->delayed_adhesion) {
				a_npart = inter->delayed_adhesion->getForceNorm();
				tratio = inter->delayed_adhesion->ratioUptimeToActivation();
			}
			outdata_int.entryData("norm of the normal adhesion force", Dimensional::Dimension::Force, 1, a_npart);
			outdata_int.entryData("adhesion ratio uptime to activation time", Dimensional::Dimension::none, 1, tratio);
		}
	}
}
}