//
//  TimeActivatedAdhesion.h
//  LF_DEM
//
//  Created by Romain Mari on 01/11/17.
//  Copyright (c) 2017 Ryohei Seto and Romain Mari. All rights reserved.
//

/**
 \class TimeActivatedAdhesion
 \brief Adhesive force activated only after a given contact time.
 \author Romain Mari
 */


#ifndef __LF_DEM__TimeActivatedAdhesion_IO__
#define __LF_DEM__TimeActivatedAdhesion_IO__
#include "Interaction.h"
#include "TimeActivatedAdhesion.h"

namespace TActAdhesion {


inline std::vector <struct State> readStatesBStream(std::istream &input)
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

inline void writeStatesBStream(std::ostream &output,
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

} //namespace TActAdhesion
#endif //ifndef __LF_DEM__TimeActivatedAdhesion_IO__
