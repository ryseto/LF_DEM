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
#include <iostream>
#include <iomanip>
#include <fstream>
#include "vec3d.h"
#include "Box.h"

#define DELETE(x) if(x){delete [] x; x = NULL;}
class System;

class BoxSet{
private:
	double box_xsize;
	double box_ysize;
	double box_zsize;
	std::size_t x_box_nb;
	std::size_t y_box_nb;
	std::size_t z_box_nb;
	std::size_t yz_box_nb;
	std::size_t box_nb;
	bool _is_boxed;
	std::set <Box*> Boxes;
	std::set <Box*> BulkBoxes;
	std::set <Box*> TopBoxes;
	std::set <Box*> BottomBoxes;
	std::set <Box*> TopBottomBoxes;
	std::vector <Box*> box_labels;
	System* sys;
	int nb_movingneighbor_part;

	/*****
	 WhichBox(vec3d pos)
	 returns a pointer on the box containg position pos
	 *****/
	Box* whichBox(const vec3d&);
	Box* whichBox(unsigned int box_label);

	void updateNeighbors();
	void updateNeighborsExtFlow();
	
	// init methods
	void allocateBoxes();
	void positionBoxes();
	void assignNeighbors();
	void assignNeighborsBulk();
	void assignNeighborsTop();
	void assignNeighborsBottom();
	void assignNeighborsTopBottom();
	
	void assignNeighborsStaticExtFlow();
	void assignNeighobrsDynamicExtFlow(Box& bx, const vec3d &pos);

	Box** boxMap;
	std::vector<vec3d> top_probing_positions;
	std::vector<vec3d> bottom_probing_positions;
	
	unsigned int whichBoxLabel(const vec3d&);
	unsigned int labelRight(int unsigned label_) {
		return label_+yz_box_nb;
	}
	unsigned int labelLeft(int unsigned label_) {
		return label_-yz_box_nb;
	}
	unsigned int labelUp(int unsigned label_) {
		return label_+1;
	}
	unsigned int labelDown(int unsigned label_) {
		return label_-1;
	}
	std::vector<int> next_labels;
	std::vector < std::vector<int> > next_of_periodicimage;
	vec3d origin_ext_flow;
	std::vector<vec3d> box_corners;

public:
	BoxSet(){;}
	~BoxSet();
	void init(double interaction_dist, System *sys_);
	void initExtFlow(double interaction_dist, const vec3d origin, System *sys_);
	/*****
	 update()

	 To be called at each time step.
	 It updates the neighborhood relations betwenn boxes.
	 Those relations change at each time step for boxes on top or bottom
	 *****/
	void update();
	void updateExtFlow();
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
	const std::vector <int>& neighborhood(int i);
	void printBoxNetwork();
	void printBoxContainers();
	void printNeighborhoodContainers();
	void printBoxMap();
	void yaplotBox(std::ofstream &fout_boxing);
	void checkNeighobrType(std::set<unsigned int> &box_labels);
	int boxType(int i)
	{
		return boxMap[i]->type;
	}
	int boxNType(int i)
	{
		return boxMap[i]->type_neighborhood;
	}
	vec3d periodicDiffShift(int i, int j);

};
#endif /* defined(__LF_DEM__BoxSet__) */
