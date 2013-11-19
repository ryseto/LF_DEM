//
//  main.cpp
//  cluster_life
//
//  Created by Ryohei Seto on 10/10/13.
//  Copyright (c) 2013 Ryohei Seto. All rights reserved.
//

#include <iostream>
#include <iomanip>
#include <fstream>
#include "Cluster.h"
#include "StrainSnapShot.h"
using namespace std;
int main(int argc, const char * argv[])
{
	ifstream fin;
	fin.open(argv[1]);
	char buf[1000];
	
	vector <Cluster*> cluster;
	vector <StrainSnapShot*> sss;
	int num_cluster = 0;
	cluster.push_back(new Cluster());
	std::string delimiter = " ";
	while(fin.getline(buf, 1000)) {
		if (buf[0] == NULL) {
			exit(1);
		} else {
			string line = buf;
			size_t pos = 0;
			string token;
			while ((pos = line.find(delimiter)) != std::string::npos) {
				token = line.substr(0, pos);
				cluster[num_cluster]->add(stoi(token));
				line.erase(0, pos + delimiter.length());
			}
			cerr << "n = " << num_cluster << endl;
			cluster[num_cluster]->printmember();
			cerr << "---" << endl;
			num_cluster++;
			cluster.push_back(new Cluster());
		}
		
	}
	// insert code here...
	std::cout << "Hello, World!\n";
    return 0;
}







