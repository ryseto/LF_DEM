//
//  SystemHelperFunctions.cpp
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 11/14/12.
//  Copyright (c) 2012-2015 Ryohei Seto and Romain Mari. All rights reserved.
//

#include "System.h"
#include <stdexcept>

using namespace std;


double evaluateMaxInterNormalVelocity(const System & sys)
{
	double max_normal_velocity = 0;
	for (const auto &inter: sys.interaction) {
		double normal_velocity = std::abs(inter.getNormalVelocity());
		if (normal_velocity > max_normal_velocity) {
			max_normal_velocity = normal_velocity;
		}
	}
	return max_normal_velocity;
}

double evaluateMaxContactSlidingVelocity(const System & sys)
{
	if (sys.friction) {
		double sq_max_sliding_velocity = 0;
		for (const auto &inter: sys.interaction) {
			if (inter.contact.is_active()) {
				double sq_sliding_velocity = inter.contact.getSlidingVelocity().sq_norm();
				if (sq_sliding_velocity > sq_max_sliding_velocity) {
					sq_max_sliding_velocity = sq_sliding_velocity;
				}
			}
		}
		return std::sqrt(sq_max_sliding_velocity);
	} else {
		return 0;
	}
}

double evaluateMaxContactRollingVelocity(const System & sys)
{
	if (sys.rolling_friction) {
		double sq_max_rolling_velocity = 0;
		for (const auto &inter: sys.interaction) {
			if (inter.contact.is_active()) {
				double sq_rolling_velocity = inter.contact.getRollingVelocity().sq_norm();
				if (sq_rolling_velocity > sq_max_rolling_velocity) {
					sq_max_rolling_velocity = sq_rolling_velocity;
				}
			}
		}
		return std::sqrt(sq_max_rolling_velocity);
	} else {
		return 0;
	}
}

double evaluateMaxNAVelocityComponent(const System & sys, std::string component)
{
	double sq_max_na_velocity = 0;
	for (unsigned i=0; i<sys.get_np(); i++) {
		auto sq_na_velocity = sys.na_velo_components.at(component).vel[i].sq_norm();
		if (sq_na_velocity > sq_max_na_velocity) {
				sq_max_na_velocity = sq_na_velocity;
		}
	}
	return std::sqrt(sq_max_na_velocity);
}
