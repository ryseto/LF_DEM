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

namespace Boxing
{

class Box{
private:
	std::vector <Box*> static_neighbors;
	std::vector <Box*> moving_neighbors;
	std::set <int> container;
	std::vector <int> neighborhood_container;
	vec3d position;

public:
	Box()
	{
		type = 0;
	};
	~Box();

	void addStaticNeighbor(Box* neigh_box);
	void addMovingNeighbor(Box* neigh_box);
	void resetMovingNeighbors();
	const std::vector <Box*> & getNeighborBox()
	{
		return static_neighbors;
	}
	vec3d getPosition() const {return position;}
	double getYPosition() const {return position.y;}
	void setPosition(vec3d pos) {position = pos;}
	void add(int);
	void remove(int);
	const std::set <int> & getContainer() const {return container;}
	const std::vector <int> & getNeighborhoodContainer() const {return neighborhood_container;}
	void buildNeighborhoodContainer();
	void resetTypes()
	{
		type = 0;
		type_neighborhood = 0;
	}
	/********************************************************************************
	 * `type' of cell is determined according to the location.
	 * 1:bulk
	 * 2:left
	 * 3:right
	 * 4:bottom
	 * 8:top
	 * 6 =2+4:left bottom
	 * 7 =3+4:right bottom
	 * 10=2+8:left top
	 * 11=3+8:right top
	 *  10 - 8 --11
	 *  |         |
	 *  2    1    3
	 *  |         |
	 *  6 -- 4 ---7
	 * If all self + 8 neighboring cells are 1, `type_neighborhood' is also 1.
	 * If one of 2/3/8/4 is included, `type_neighborhood' is 2/3/8/4.
	 * if two of 2/3/8/4 is included, `type_neighborhood' is the sum of them,
	 * i.e., 2 and 8 is included in self + 8 neighboring cells, the value is 10.
	 *****************************************************************************/
	int type;
	int type_neighborhood;
};

} // namespace Boxing

#endif /* defined(__LF_DEM__Box__) */
