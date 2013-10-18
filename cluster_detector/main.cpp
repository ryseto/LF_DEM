//
//  main.cpp
//  cluster_detector
//
//  Created by Ryohei Seto on 7/9/13.
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
	string file_yap = "cls_y_" + file_par.substr(3,-1);
	string file_data = "cls_d_" + file_par.substr(3,-1);
	string file_cluster = "cluster" + file_par.substr(3,-1);
	cerr << file_par << endl;
	cerr << file_int << endl;
	sys.setImportFile_par(file_par);
	sys.setImportFile_int(file_int);
	sys.openYapFile(file_yap);
	sys.openDataFile(file_data);
	sys.openClusterFile(file_cluster);
	while(true){
		sys.clear();
		sys.importConfiguration();
		sys.buildCluster();
		sys.output_yap();
		sys.output_clusters();
		if (sys.fin_p.eof()){
			break;
		}
	}
	

	
    return 0;
}

