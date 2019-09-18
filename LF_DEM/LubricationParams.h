#ifndef __LF_DEM__LubricationParams__
#define __LF_DEM__LubricationParams__

#include <exception>

namespace Interactions
{

namespace Lub
{

struct LubParams {
	/*
		 * Leading term of lubrication force is 1/reduced_gap,
		 * with reduced_gap the gap
		 * reduced_gap = 2r/(a0+a1) - 2.
		 * we set a cutoff for the lubrication interaction,
		 * such that the lub term is proportional to:
		 *
		 * 1/(reduced_gap+regularization_length)
		 * when reduced_gap > 0.
		 */
	double max_gap;					///< Lubrication range (in interparticle gap distance) [0.5]
	double regularization_length;   ///< Lubrication regularization length ("roughness length") [1e-3]
	std::string model;              ///< Lubrication type. "none": no lubrication, "normal": 1/xi lubrication (only squeeze mode), "tangential": "normal" plus log(1/xi) terms (shear and pump modes) ["tangential"]
	bool smooth;
};

inline void setupLubricationParameters(LubParams &p)
{
	std::string indent = "  setupLubricationParameters::\t";
	if (p.max_gap < 0) {
		throw std::runtime_error(indent+"lub_max_gap<0 is forbidden.");
	}
	if (p.model != "normal" &&
		p.model != "none" &&
		p.model != "tangential") {
		throw std::runtime_error(indent+"unknown lubrication_model "+p.model+"\n");
	}
}

} // namespace Lub

inline bool has_lubrication(const Lub::LubParams &p)
{
	return p.model != "none";
}

 
} // namespace Interactions

#endif