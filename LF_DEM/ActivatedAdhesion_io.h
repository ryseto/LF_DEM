//
//  TimeActivatedAdhesion_io.h
//  LF_DEM
//
//  Created by Romain Mari on 01/11/17.
//  Copyright (c) 2017 Ryohei Seto and Romain Mari. All rights reserved.
//


#ifndef __LF_DEM__ActivatedAdhesion_IO__
#define __LF_DEM__ActivatedAdhesion_IO__

#include "ActivatedAdhesion.h"

namespace Interactions {
class StdInteractionManager;
namespace ActAdhesion {

std::vector <struct State> readStatesBStream(std::istream &input);

void writeStatesBStream(std::ostream &output,
						const Interactions::StdInteractionManager &interactions);

void setupInteractions(Interactions::StdInteractionManager &interactions, 
					   const std::vector <struct State> &adhesion_states);

} //namespace ActAdhesion
} //namespace Interactions

#endif //ifndef __LF_DEM__ActivatedAdhesion_IO__
