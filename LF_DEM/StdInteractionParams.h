#ifndef __LF_DEM__StdInteractionParams__
#define __LF_DEM__StdInteractionParams__
#include <memory>
#include "ContactParams.h"
#include "LubricationParams.h"
#include "RepulsiveForceParams.h"
#include "VanDerWaalsParams.h"

namespace Interactions
{

struct StdInteractionParams {
	std::unique_ptr<ContactParams> contp; 
	std::unique_ptr<Lub::LubParams> lubp;
	std::unique_ptr<RepulsiveForceParams> repp;
	std::unique_ptr<vanDerWaalsForceParams> vdwp;

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
};

}

#endif