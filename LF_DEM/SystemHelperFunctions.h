#ifndef __LF_DEM__SystemHF__
#define __LF_DEM__SystemHF__

inline double evaluateMaxAngVelocity(const System & sys)
{
	double _max_ang_velocity = 0;
	for (unsigned i = 0; i < sys.get_np(); i++) {
		vec3d na_ang_velocity_tmp = sys.na_ang_velocity[i];
		if (na_ang_velocity_tmp.norm() > _max_ang_velocity) {
			_max_ang_velocity = na_ang_velocity_tmp.norm();
		}
	}
	return _max_ang_velocity;
}

inline double evaluateMinGap(const System & sys)
{
	double _min_reduced_gap = sys.p.interaction_range;
	unsigned int p0, p1;
	for (const auto &inter: sys.interaction) {
		std::tie(p0, p1) = inter.get_par_num();
		if ((int)p0 < sys.np_mobile && // exclude fixed-fixed
			inter.get_reduced_gap() < _min_reduced_gap) {
			_min_reduced_gap = inter.get_reduced_gap();
		}
	}
	return _min_reduced_gap;
}

inline double evaluateAvgContactGap(const System & sys)
{
	double _avg_reduced_gap = 0;
	int nb=0;
	for (const auto &inter: sys.interaction) {
		if (inter.contact.is_active()) {
			_avg_reduced_gap += inter.get_reduced_gap();
			nb++;
		}
	}
	_avg_reduced_gap /= nb;
	return _avg_reduced_gap;
}

inline double evaluateMaxContactGap(const System & sys)
{
	double _max_contact_gap = 0;
	for (const auto &inter: sys.interaction) {
		if (inter.contact.is_active() &&
			inter.get_reduced_gap() > _max_contact_gap) {
			_max_contact_gap = inter.get_reduced_gap();
		}
	}
	return _max_contact_gap;
}

inline double evaluateMaxVelocity(const System &sys)
{
	double sq_max_velocity = 0;
	for (unsigned i = 0; i < sys.get_np(); i++) {
		vec3d na_velocity_tmp = sys.na_velocity[i];
		if (na_velocity_tmp.sq_norm() > sq_max_velocity) {
			sq_max_velocity = na_velocity_tmp.sq_norm();
		}
	}
	return sqrt(sq_max_velocity);
}

double evaluateMaxInterNormalVelocity(const System & sys);
double evaluateMaxContactSlidingVelocity(const System & sys);
double evaluateMaxContactRollingVelocity(const System & sys);
double evaluateMaxNAVelocityComponent(const System & sys, std::string component);

inline double evaluateMaxDispTan(const System &sys)
{
	double _max_disp_tan = 0;
	for (const auto &inter: sys.interaction) {
		if (inter.contact.disp_tan.norm() > _max_disp_tan) {
			_max_disp_tan = inter.contact.disp_tan.norm();
		}
	}
	return _max_disp_tan;
}

inline double evaluateMaxDispRolling(const System &sys)
{
	double _max_disp_rolling = 0;
	for (const auto &inter: sys.interaction) {
		if (inter.contact.disp_rolling.norm() > _max_disp_rolling) {
			_max_disp_rolling = inter.contact.disp_rolling.norm();
		}
	}
	return _max_disp_rolling;
}

inline std::pair<unsigned int, unsigned int> countNumberOfContact(const System &sys)
{
	unsigned int contact_nb = 0;
	unsigned int fric_contact_nb = 0;
	for (const auto &inter: sys.interaction) {
		if (inter.contact.is_active()) {
			contact_nb ++;
			if (inter.contact.is_frictional()) {
				fric_contact_nb ++;
			}
		}
	}
	return std::make_pair(contact_nb, fric_contact_nb);
}

inline std::pair<unsigned, double> getTAAdhesionActivityStatistics(const System &sys)
{
	assert(sys.delayed_adhesion);
	unsigned active_nb = 0;
	unsigned active_dormant_nb = 0;
	for (const auto &inter: sys.interaction) {
		auto state = inter.delayed_adhesion->getState();
		if (state.activity != TActAdhesion::Activity::inactive) {
			active_dormant_nb++;
			if (state.activity == TActAdhesion::Activity::active) {
				active_nb++;
			}
		}
	}
	double active_ratio = 0;
	if (active_dormant_nb>0) {
		active_ratio = (double)(active_nb)/active_dormant_nb;
	}
	return std::make_pair(active_nb, active_ratio);
}

inline double getPotentialEnergy(const System &sys)
{
	double total_energy = 0;
	for (const auto &inter: sys.interaction) {
		if (inter.contact.is_active()){
			total_energy += inter.contact.calcEnergy();
		}
		if (sys.repulsiveforce) {
			total_energy += inter.repulsion.calcEnergy();
		}
	}
	return total_energy;
}

inline void isInContact(const System &sys, std::vector<int> &isincontact)
{
	isincontact.resize(sys.get_np(), 0);
	for (unsigned i=0; i<sys.get_np(); i ++) {
		for (auto& inter : sys.interaction_list[i]) {
			if (inter->contact.is_active()) {
				isincontact[i] = 1;
			}
		}
	}
 }

#endif /* defined(__LF_DEM__SystemHF__) */
