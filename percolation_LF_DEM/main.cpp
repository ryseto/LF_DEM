//
//  main.cpp
//  percolation_LF_DEM
//
//  Created by Ryohei Seto on 7/24/13.
//  Copyright (c) 2013 Ryohei Seto. All rights reserved.
//

#include <iostream>
#include "System.h"

using namespace std;

int main(int argc, const char * argv[])
{
	System sys;
	string file_par = argv[1];
	string file_int = "int" + file_par.substr(3,-1);
	cerr << file_par << endl;
	cerr << file_int << endl;
	sys.setImportFile_par(file_par);
	sys.setImportFile_int(file_int);
	while(true){
		sys.importConfiguration();
		sys.max_cluster_size = 0;
		for (int i=0; i< sys.np; i++){
			sys.check_percolation(i);
		}
		cerr << sys.strain << ' ' << (1.*sys.max_cluster_size)/sys.np << endl;
		if (sys.strain >= 40){
			break;
		}
	}

    return 0;
}

