//
//  TimeActivatedAdhesion_io.cpp
//  LF_DEM
//
//  Created by Romain Mari on 01/11/17.
//  Copyright (c) 2017 Ryohei Seto and Romain Mari. All rights reserved.
//


#include <assert.h>
#include "TimeActivatedAdhesion_io.h"


namespace TActAdhesion {


std::vector <struct State> readStatesBStream(std::istream &input)
{
	unsigned ncont;
	Activity activity;
	State state;	
	std::vector <struct State> state_vec;

	input.read((char*)&ncont, sizeof(unsigned));
	for (unsigned i=0; i<ncont; i++) {
		unsigned p0, p1;
		input.read((char*)&p0, sizeof(unsigned));
		input.read((char*)&p1, sizeof(unsigned));
		input.read((char*)&activity, sizeof(unsigned));
		input.read((char*)&state.uptime, sizeof(double));

		state.activity = activity;
		state_vec.push_back(state);
	}
	return state_vec;
}

void writeStatesBStream(std::ostream &output,
			    				const std::vector<Interaction> &interactions)
{
	std::vector<State> adhesion_states;
	for (auto &inter: interactions) {
		adhesion_states.push_back(inter.delayed_adhesion->getState());
	}

	unsigned ncont = adhesion_states.size();
	output.write((char*)&ncont, sizeof(unsigned));
	for (auto &state: adhesion_states) {
		output.write((char*)&state.p0, sizeof(unsigned));
		output.write((char*)&state.p1, sizeof(unsigned));
		output.write((char*)&state.activity, sizeof(unsigned));
		output.write((char*)&state.uptime, sizeof(double));
	}
}

void setupInteractions(std::vector<Interaction> &interactions, 
					          const std::vector <struct State> &adhesion_states, 
					          double time_now)
{
	/**
		\brief Set a list of adhesions with their state variables.

		Used to restart the simulation from a given state.
	 */
	for (const auto& state : adhesion_states) {
		bool existing_interaction = false;
		for (auto &inter: interactions) {
			unsigned int p0, p1;
			std::tie(p0, p1) = inter.get_par_num();
			if (p0 == state.p0 && p1 == state.p1) {
				inter.delayed_adhesion->setState(state, time_now);
				existing_interaction = true;
			}
		}
		assert(existing_interaction); // otherwise you probably messed up the initialisation
	}
}

} //namespace TActAdhesion