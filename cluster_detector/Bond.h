//
//  Bond.h
//  LF_DEM
//
//  Created by Ryohei Seto on 7/9/13.
//  Copyright (c) 2013 Ryohei Seto. All rights reserved.
//

#ifndef __LF_DEM__Bond__
#define __LF_DEM__Bond__

#include <iostream>
using namespace std;

class Bond{
private:


protected:
public:
	unsigned int par_num[2];
	int contact_state;
	void set(int i, int j, int contact_state_){
		contact_state = contact_state_;
		if (i < j){
			par_num[0] = i;
			par_num[1] = j;
		} else {
			par_num[0] = j;
			par_num[1] = i;
		}
	}
	
	int next(int i){
		if (contact_state == 1){
			if (i == par_num[0]){
				return par_num[1];
			} else if(i == par_num[1]){
				return par_num[0];
			}
		}
		return -1;
	}
	
};
#endif /* defined(__LF_DEM__Bond__) */

