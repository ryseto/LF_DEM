//
//  Cluster.cpp
//  LF_DEM
//
//  Created by Ryohei Seto on 7/9/13.
//  Copyright (c) 2013 Ryohei Seto. All rights reserved.
//

#include "Cluster.h"

Cluster::Cluster(int i, vec3d pos_){
	member.push_back(i);
	pos.push_back(pos_);
};