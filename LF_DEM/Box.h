//
//  Box.h
//  LF_DEM
//
//  Copyright (c) 2013-2015 Ryohei Seto and Romain Mari. All rights reserved.
//

/**
 \class Box
 \brief Box object holding identities of the particles in a subset of the whole suspension, to be used in a BoxSet object.
 \author Ryohei Seto
 \author Romain Mari
 */

#ifndef __LF_DEM__Box__
#define __LF_DEM__Box__
#define __vector_container__
#ifndef __vector_container__
#define __flist_container__
#endif
#include "vec3d.h"
#include <set>
#include <vector>

class Box{
private:
	std::set <Box*> _neighbors;
	std::set <Box*> _moving_neighbors;
	std::set <int> container;
	std::vector <int> neighborhood_container;
	vec3d position;

public:
	Box(){};
	~Box();

	void addStaticNeighbor(Box* neigh_box);
	void addMovingNeighbor(Box* neigh_box);
	void resetMovingNeighbors();
	const std::set <Box*> & getNeighborBox(){return _neighbors;}

	vec3d getPosition() const {return position;}
	void setPosition(vec3d pos) {position = pos;}
	// void is_top(bool it) {_is_top = it;};
	// void is_bottom(bool ib) {_is_bottom = ib;};
	// bool is_top() const {return _is_top;};
	// bool is_bottom() const {return _is_bottom;};

	void add(int);
	void remove(int);

	const std::set <int> & getContainer() const {return container;}
	const std::vector <int> & getNeighborhoodContainer() const {return neighborhood_container;};
	void buildNeighborhoodContainer();
};

#endif /* defined(__LF_DEM__Box__) */
