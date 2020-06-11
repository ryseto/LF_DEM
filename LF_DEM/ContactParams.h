#ifndef __LF_DEM__ContactParams__
#define __LF_DEM__ContactParams__

namespace Interactions {

enum class FrictionModel : unsigned {
	frictionless = 0,
	Coulomb = 1,
	criticalload = 2,
	criticalload_mu_inf = 3,
	infinity = 4,
	ft_max = 5,
	Coulomb_max = 6,
	ageing_Coulomb = 7
};

struct ContactParams {
	double kn;						///< Particle stiffness: normal spring constant [0h]
	double kt;						///< Particle stiffness: tangential spring constant [0kn]
	double kr;						///< Particle stiffness: rolling spring constant [0kn]
	FrictionModel friction_model;   ///< Friction model from the list above  [1]
	double relaxation_time;			///< Relaxation time (normal) of the contact model. Sets the normal dashpot. If <0, use normal lubrication at contact as normal dashpot. [1e-3input_unit]
	double relaxation_time_tan;		///< Relaxation time (tangential) of the contact model. Sets the tangential dashpot. If <0, use tangential lubrication at contact as tangential dashpot. [-1input_unit]
	double critical_load;			///< Amplitude of the critical load [0]
	double ft_max;					///< max tangential force in friction_model = 5 [0 guarranted_unit]
	double mu_static;				///< friction coefficient (static) [1]
	double mu_dynamic;				///< friction coefficient (dynamic). If -1, mu_dynamic = mu_static [-1]
	double mu_rolling;				///< friction coefficient (rolling) [0]
	double adhesion;				///< Amplitude of the adhesion [0 guarranted_unit]
};

inline bool has_friction(FrictionModel fm)
{
	return fm != FrictionModel::frictionless;
}

inline bool has_rolling_friction(const ContactParams &p)
{
	return p.mu_rolling != 0;
}

inline bool has_adhesion(const ContactParams &p)
{
	return p.adhesion > 0;
}

inline bool has_dashpot(const ContactParams &p)
{
	return p.relaxation_time != 0 || p.relaxation_time_tan != 0;
}


inline void setupContactParameters(ContactParams &p)
{
	std::string indent = "  setupContactParameters::\t";
	if (p.friction_model == FrictionModel::frictionless) {
		std::cout << indent+"friction model: no friction" << std::endl;
		p.mu_static = 0;
	} else if (p.friction_model == FrictionModel::Coulomb) {
		std::cout << indent+"friction model: Coulomb" << std::endl;
	} else if (p.friction_model == FrictionModel::criticalload || p.friction_model == FrictionModel::criticalload_mu_inf) {
		std::cout << indent+"friction model: Coulomb + Critical Load" << std::endl;
	} else if (p.friction_model == FrictionModel::infinity) {
		std::cout << indent+"friction model: Static friction (mu = infinity)" << std::endl;
	} else if (p.friction_model == FrictionModel::ft_max) {
		std::cout << indent+"friction_model: Max tangential force" << std::endl;
	} else if (p.friction_model == FrictionModel::Coulomb_max) {
		std::cout << indent+"friction_model: Coulomb law + Max tangential force" << std::endl;
	} else {
		throw std::runtime_error(indent+"Error: unknown friction model\n");
	}
	if (p.mu_dynamic < 0) {
		p.mu_dynamic = p.mu_static;
	}
}

} // namespace Interaction

#endif