//
//  BoxSet.h
//  LF_DEM
//
//  Copyright (c) 2013-2015 Ryohei Seto and Romain Mari. All rights reserved.
//

/**
 \class BoxSet
 \brief Set of Box objects making a partition of the simulation box.
 \author Ryohei Seto
 \author Romain Mari
 */

#ifndef __LF_DEM__BoxSet__
#define __LF_DEM__BoxSet__
#include <set>
#include <unordered_set>
#include <vector>
#include "vec3d.h"
#include "Box.h"

#define DELETE(x) if(x){delete [] x; x = NULL;}
class System;

class BoxSet{
private:
	double box_xsize;
	double box_ysize;
	double box_zsize;
	int x_box_nb;
	int y_box_nb;
	int z_box_nb;
	int box_nb;
	bool _is_boxed;
	std::set <Box*> Boxes;
	std::set <Box*> BulkBoxes;
	std::set <Box*> TopBoxes;
	std::set <Box*> BottomBoxes;
	std::set <Box*> TopBottomBoxes;
	std::vector <Box*> box_labels;
	System* sys;
	int amax, bmax, cmax; // amax = min( x_box_nb, 3), bmax = min( y_box_nb, 3), cmax = min( z_box_nb, 3)
	Box* WhichBox(vec3d*);
	Box* WhichBox_periodize(vec3d);
	void updateNeighbors();
	// init methods
	void allocateBoxes();
	void positionBoxes();
	void assignNeighbors();
	void assignNeighborsBulk();
	void assignNeighborsTop();
	void assignNeighborsBottom();
	void assignNeighborsTopBottom();
	Box** boxMap;
	std::vector<vec3d> top_probing_positions;
	std::vector<vec3d> bottom_probing_positions;
public:
	BoxSet(){;}
	~BoxSet();
	void init(double interaction_dist, System *sys_);
	/*****
	 update()
	 
	 To be called at each time step.
	 It updates the neighborhood relations betwenn boxes.
	 Those relations change at each time step for boxes on top or bottom
	 *****/
	void update();
	/*****
	 is_boxed()
	 
	 Can be called before calling an other method of BoxSet.
	 If false, than calls to other method may usually be avoided.
	 They can be performed anyway though, and they are normally safe (and useless)
	 
	 is_boxed() tells if the boxing is effective.
	 If the system size is small, the neighborhood of a box may contain
	 the whole system. If this happens, the boxing consists of only one box (as it
	 is useless, if not ill defined, to do something else), and is_boxed() returns false.
	 
	 *****/
	bool is_boxed();
	/*****
	 WhichBox(vec3d pos)
	 
	 returns a pointer on the box containg position pos
	 *****/
	Box* WhichBox(vec3d);
	/*****
	 box(int i)
	 boxes particles i
	 should be called after moving particle i
	 *****/
	void box(int i);
	/*****
	 neighborhood_begin(int i) and neighborhood_end(int i)
	 gives iterators to beginning and ending point of the container including
	 all particles in the box containing particle i and in the adjacent boxes.
	 *****/
	std::vector<int>::iterator neighborhood_begin(int i);
	std::vector<int>::iterator neighborhood_end(int i);
	std::vector <int>& neighborhood(int i);
	void printBoxNetwork();
	void printBoxContainers();
	void printNeighborhoodContainers();
	void printBoxMap();
};
#endif /* defined(__LF_DEM__BoxSet__) */
