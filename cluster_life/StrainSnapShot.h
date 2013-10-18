//
//  StrainSnapShot.h
//  LF_DEM
//
//  Created by Ryohei Seto on 10/10/13.
//  Copyright (c) 2013 Ryohei Seto. All rights reserved.
//

#ifndef __LF_DEM__StrainSnapShot__
#define __LF_DEM__StrainSnapShot__

#include <iostream>
#include "Cluster.h"
#include <vector>


class StrainSnapShot{
private:
	
protected:
public:
	StrainSnapShot(){;}
	~StrainSnapShot(){;}
	vector <Cluster*> cluster;
	int num_cluster;
};


#endif /* defined(__LF_DEM__StrainSnapShot__) */
