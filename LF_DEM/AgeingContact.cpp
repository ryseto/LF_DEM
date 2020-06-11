//
//
//  Created by Romain Mari on 9/06/20.
//  Copyright (c) 2020 Romain Mari. All rights reserved.
//
#include "AgeingContact.h"

namespace Interactions 
{

AgeingContact::AgeingContact(PairwiseInteraction* interaction_, 
							 const ContactParams &p, 
							 const AgeingContactParams &ageing_params, 
							 double norm_dashpot_coeff, 
							 double tan_dashpot_coeff) :
Contact(interaction_, p, norm_dashpot_coeff, tan_dashpot_coeff),
ageing(ageing_params),
state({0})
{}

void AgeingContact::incrementDisplacements(double dt, const struct PairVelocity &vel)
{
	state.contact_time += dt;
	mu_static = mu_dynamic + ageing.amplitude*std::log(state.contact_time/ageing.timescale);
	Contact::incrementDisplacements(dt, vel);
}

struct AgeingState AgeingContact::getAgeingState() const
{
	return state;
}

void AgeingContact::setState(const struct AgeingState& as)
{
	state = as;
}

void AgeingContact::saveState()
{
	prev_state = state;	
}

void AgeingContact::restoreState()
{
	state = prev_state;
}

} // namespace Interactions