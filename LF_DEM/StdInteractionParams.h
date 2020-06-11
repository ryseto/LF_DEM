#ifndef __LF_DEM__StdInteractionParams__
#define __LF_DEM__StdInteractionParams__
#include <memory>
#include "ContactParams.h"
#include "LubricationParams.h"
#include "RepulsiveForceParams.h"
#include "VanDerWaalsParams.h"
#include "ActivatedAdhesion_Params.h"
#include "AgeingContactParams.h"


namespace Interactions
{

struct StdInteractionParams {
	std::unique_ptr<ContactParams> contp; 
	std::unique_ptr<Lub::LubParams> lubp;
	std::unique_ptr<RepulsiveForceParams> repp;
	std::unique_ptr<vanDerWaalsForceParams> vdwp;
	std::unique_ptr<ActAdhesion::Params> actadhp;
	std::unique_ptr<AgeingContactParams> ageing_contp;

	void add(const ContactParams &p) {
		contp = std::unique_ptr<ContactParams>(new ContactParams (p));
	}
	void add(const Lub::LubParams &p) {
		lubp = std::unique_ptr<Lub::LubParams>(new Lub::LubParams (p));
	}
	void add(const RepulsiveForceParams &p) {
		repp = std::unique_ptr<RepulsiveForceParams>(new RepulsiveForceParams (p));
	}
	void add(const vanDerWaalsForceParams &p) {
		vdwp = std::unique_ptr<vanDerWaalsForceParams>(new vanDerWaalsForceParams (p));
	}
	void add(const ActAdhesion::Params &p) {
		actadhp = std::unique_ptr<ActAdhesion::Params>(new ActAdhesion::Params (p));
	}
	void add(const AgeingContactParams &p) {
		ageing_contp = std::unique_ptr<AgeingContactParams>(new AgeingContactParams (p));
	}
};

}

#endif