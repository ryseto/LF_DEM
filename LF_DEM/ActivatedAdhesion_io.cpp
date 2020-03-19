//
//  TimeActivatedAdhesion_io.cpp
//  LF_DEM
//
//  Created by Romain Mari on 01/11/17.
//  Copyright (c) 2017 Ryohei Seto and Romain Mari. All rights reserved.
//


#include <assert.h>
#include "ActivatedAdhesion_io.h"
#include "StdInteractionManager.h"

namespace Interactions {

namespace ActAdhesion {

std::vector <struct State> readStatesBStream(std::istream &input)
{
	unsigned ncont;
	typedef std::underlying_type<Activity>::type activity_type;
	activity_type activity;
	State state;
	std::vector <struct State> state_vec;

	input.read((char*)&ncont, sizeof(unsigned));
	for (unsigned i=0; i<ncont; i++) {
		input.read((char*)&state.p0, sizeof(unsigned));
		input.read((char*)&state.p1, sizeof(unsigned));
		input.read((char*)&activity, sizeof(activity_type));

		state.activity = static_cast<Activity>(activity);
		state_vec.push_back(state);
	}
	return state_vec;
}

void writeStatesBStream(std::ostream &output,
						const Interactions::StdInteractionManager &interactions)
{
	std::vector<State> adhesion_states;
	for (auto &inter: interactions) {
		if (inter->act_adhesion) {
			adhesion_states.push_back(inter->act_adhesion->getState());
		}
	}

	unsigned ncont = adhesion_states.size();
	output.write((char*)&ncont, sizeof(unsigned));
	typedef std::underlying_type<Activity>::type activity_type;
	
	for (auto &state: adhesion_states) {
		activity_type activity = static_cast<activity_type>(state.activity);
		output.write((char*)&state.p0, sizeof(unsigned));
		output.write((char*)&state.p1, sizeof(unsigned));
		output.write((char*)&activity, sizeof(activity_type));
	}
}

void setupInteractions(Interactions::StdInteractionManager &interactions, 
					   const std::vector <struct State> &adhesion_states)
{
	/**
	 \brief Set a list of adhesions with their state variables.
	 
	 Used to restart the simulation from a given state.
	 */
	for (const auto& state : adhesion_states) {
		bool existing_interaction = false;
		for (auto &inter: interactions) {
			unsigned int p0, p1;
			std::tie(p0, p1) = inter->get_par_num();
			if (p0 == state.p0 && p1 == state.p1) {
				inter->act_adhesion->setState(state);
				existing_interaction = true;
			}
		}
		assert(existing_interaction); // otherwise you probably messed up the initialisation
	}
}

} //namespace TActAdhesion

}