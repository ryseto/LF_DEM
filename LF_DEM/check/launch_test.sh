#!/bin/sh

~/recherche/codes/LF_DEM/LF_DEM/LF_DEM -p 10 -k knkt_mu1.dat -n confs/D3N500VF0.45Bidi1.4_0.5Cubic_1_colloid_determine_knkt_0.04_p10.bin colloid_determine_knkt_0.04.txt
~/recherche/codes/LF_DEM/LF_DEM/LF_DEM -r 0.2 confs/D3N500VF0.45Bidi1.4_0.5Cubic_1_repulsive_r0.2.bin -k knkt_mu1.dat -n repulsive.txt
~/recherche/codes/LF_DEM/LF_DEM/LF_DEM -s 5 -k knkt_mu1.dat confs/D3N500VF0.45Bidi1.4_0.5Cubic_1_repulsive_s5.bin repulsive.txt
~/recherche/codes/LF_DEM/LF_DEM/LF_DEM -p 100 -k knkt_mu1_peclet_ov0.02.txt -n confs/D3N500VF0.45Bidi1.4_0.5Cubic_1_colloid_p100.bin colloid.txt
~/recherche/codes/LF_DEM/LF_DEM/LF_DEM -p 0.1 -k knkt_mu1_peclet_ov0.02.txt -n confs/D3N500VF0.45Bidi1.4_0.5Cubic_1_colloid_p0.1.bin colloid.txt


diff 
