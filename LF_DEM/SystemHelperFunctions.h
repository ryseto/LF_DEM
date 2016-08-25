#ifndef __LF_DEM__SystemHF__
#define __LF_DEM__SystemHF__

inline double evaluateMaxAngVelocity(const System & sys)
{
	double _max_ang_velocity = 0;
	for (int i = 0; i < sys.get_np(); i++) {
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
	for (int k=0; k<sys.nb_interaction; k++) {
		if (sys.interaction[k].is_active()) {
			std::tie(p0, p1) = sys.interaction[k].get_par_num();
			if ((int)p0 < sys.np_mobile && // exclude fixed-fixed
				sys.interaction[k].get_reduced_gap() < _min_reduced_gap) {
				_min_reduced_gap = sys.interaction[k].get_reduced_gap();
			}
		}
	}
	return _min_reduced_gap;
}

inline double evaluateAvgContactGap(const System & sys)
{
	double _avg_reduced_gap = 0;
	int nb=0;
	for (int k=0; k<sys.nb_interaction; k++) {
		if (sys.interaction[k].is_active() &&
			sys.interaction[k].contact.is_active()) {
			_avg_reduced_gap += sys.interaction[k].get_reduced_gap();
			nb++;
		}
	}
	_avg_reduced_gap /= nb;
	return _avg_reduced_gap;
}

inline double evaluateMaxContactGap(const System & sys)
{
	double _max_contact_gap = 0;
	for (int k= 0; k<sys.nb_interaction; k++) {
		if (sys.interaction[k].is_active() &&
			sys.interaction[k].contact.is_active() &&
			sys.interaction[k].get_reduced_gap() > _max_contact_gap) {
			_max_contact_gap = sys.interaction[k].get_reduced_gap();
		}
	}
	return _max_contact_gap;
}

inline double evaluateMaxVelocity(const System &sys)
{
	double sq_max_velocity = 0;
	for (int i = 0; i < sys.get_np(); i++) {
		vec3d na_velocity_tmp = sys.na_velocity[i];
		if (na_velocity_tmp.sq_norm() > sq_max_velocity) {
			sq_max_velocity = na_velocity_tmp.sq_norm();
		}
	}
	return sqrt(sq_max_velocity);
}

inline double evaluateMaxDispTan(const System &sys)
{
	double _max_disp_tan = 0;
	for (int k= 0; k<sys.nb_interaction; k++) {
		if (sys.interaction[k].is_active() &&
			sys.interaction[k].contact.disp_tan.norm() > _max_disp_tan) {
			_max_disp_tan = sys.interaction[k].contact.disp_tan.norm();
		}
	}
	return _max_disp_tan;
}

inline double evaluateMaxDispRolling(const System &sys)
{
	double _max_disp_rolling = 0;
	for (int k= 0; k<sys.nb_interaction; k++) {
		if (sys.interaction[k].is_active() &&
			sys.interaction[k].contact.disp_rolling.norm() > _max_disp_rolling) {
			_max_disp_rolling = sys.interaction[k].contact.disp_rolling.norm();
		}
	}
	return _max_disp_rolling;
}

inline double evaluateMaxFcNormal(const System &sys)
{
	double max_fc_normal_ = 0;
	for (int k=0; k<sys.nb_interaction; k++) {
		if (sys.interaction[k].is_active() &&
			sys.interaction[k].contact.is_active() &&
			sys.interaction[k].contact.getNormalForce().norm() > max_fc_normal_) {
				max_fc_normal_ = sys.interaction[k].contact.getNormalForce().norm();
		}
	}
	return max_fc_normal_;
}

inline std::pair<unsigned int, unsigned int> countNumberOfContact(const System &sys)
{
	unsigned int contact_nb = 0;
	unsigned int fric_contact_nb = 0;
	for (int k=0; k<sys.nb_interaction; k++) {
		if (sys.interaction[k].is_active() &&
			sys.interaction[k].contact.is_active()) {
			contact_nb ++;
			if (sys.interaction[k].contact.is_frictional()) {
				fric_contact_nb ++;
			}
		}
	}
	return std::make_pair(contact_nb, fric_contact_nb);
}

inline double getPotentialEnergy(const System &sys)
{
	double total_energy = 0;
	for (int k=0; k<sys.nb_interaction; k++) {
		if (sys.interaction[k].is_active()) {
			if (sys.interaction[k].contact.is_active()){
				total_energy += sys.interaction[k].contact.calcEnergy();
			}
			if (sys.repulsiveforce) {
				total_energy += sys.interaction[k].repulsion.calcEnergy();
			}
		}
	}
	return total_energy;
}
#endif /* defined(__LF_DEM__SystemHF__) */
