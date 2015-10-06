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
	bool _is_bottom;
	bool _is_top;
	
public:
	Box();
	~Box();
	vec3d position;
	std::set <int> container;
	void addStaticNeighbor(Box* neigh_box);
	void addMovingNeighbor(Box* neigh_box);
	void reset_moving_neighbors();
	std::set <Box*> neighbors()
	{
		return _neighbors;
	}
	
	void is_top(bool);
	void is_bottom(bool);
	bool is_top();
	bool is_bottom();
	void add(int);
	void remove(int);
	
	std::set<int>::iterator begin()
	{
		return container.begin();
	}
	
	std::set<int>::iterator end()
	{
		return container.end();
	}
	
	std::vector<int>::iterator neighborhood_begin()
	{
		return neighborhood_container.begin();
	}
	
	std::vector<int>::iterator neighborhood_end()
	{
		return neighborhood_container.end();
	}
	
	std::vector <int> neighborhood_container;
	
	size_t container_size()
	{
		return container.size();
	}
	
	void build_neighborhood_container();
};

#endif /* defined(__LF_DEM__Box__) */
