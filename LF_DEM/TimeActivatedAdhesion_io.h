//
//  TimeActivatedAdhesion_io.h
//  LF_DEM
//
//  Created by Romain Mari on 01/11/17.
//  Copyright (c) 2017 Ryohei Seto and Romain Mari. All rights reserved.
//


#ifndef __LF_DEM__TimeActivatedAdhesion_IO__
#define __LF_DEM__TimeActivatedAdhesion_IO__

#include "TimeActivatedAdhesion.h"

namespace Interactions {
class StdInteractionManager;
namespace TActAdhesion {

std::vector <struct State> readStatesBStream(std::istream &input);

void writeStatesBStream(std::ostream &output,
						const Interactions::StdInteractionManager &interactions);

void setupInteractions(Interactions::StdInteractionManager &interactions, 
					   const std::vector <struct State> &adhesion_states,
					   double time_now);

} //namespace TActAdhesion
} //namespace Interactions

#endif //ifndef __LF_DEM__TimeActivatedAdhesion_IO__
