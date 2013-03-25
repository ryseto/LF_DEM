//
//  common.h
//  LF_DEM
//
//  Created by Ryohei Seto on 1/27/13.
//  Copyright (c) 2013 Ryohei Seto. All rights reserved.
//

#ifndef LF_DEM_common_h
#define LF_DEM_common_h

struct stresslet {
	/* element[]= {S_xx, S_xy, S_xz, S_yz, S_yy}
	 */
	double elm[5];
};

#endif
