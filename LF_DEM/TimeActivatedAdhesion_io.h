//
//  TimeActivatedAdhesion_io.h
//  LF_DEM
//
//  Created by Romain Mari on 01/11/17.
//  Copyright (c) 2017 Ryohei Seto and Romain Mari. All rights reserved.
//


#ifndef __LF_DEM__TimeActivatedAdhesion_IO__
#define __LF_DEM__TimeActivatedAdhesion_IO__

#include "Interaction.h"
#include "TimeActivatedAdhesion.h"

namespace TActAdhesion {

std::vector <struct State> readStatesBStream(std::istream &input);

void writeStatesBStream(std::ostream &output,
			    				const std::vector<Interaction> &interactions);

void setupInteractions(std::vector<Interaction> &interactions, 
					          const std::vector <struct State> &adhesion_states, 
					          double time_now);

} //namespace TActAdhesion
#endif //ifndef __LF_DEM__TimeActivatedAdhesion_IO__
