#include <string>
#include "SimpleShearBoxSet.h"
#include "LeesEdwards.h"

namespace Boxing {

SimpleShearBoxSet::SimpleShearBoxSet(double interaction_dist, unsigned np, const BC::LeesEdwardsBC &lees)
{
	std::string indent = "  BoxSet::\t";
	std::cout << indent << "Setting up Cell List System ... ";
	boxMap.resize(np);
	for (unsigned i=0; i<np; i++) {
		boxMap[i] = NULL;
	}
	auto sysbox = lees.getContainer();
	x_box_nb = (unsigned int)(sysbox.lx/interaction_dist);
	y_box_nb = (unsigned int)(sysbox.ly/interaction_dist);
	z_box_nb = (unsigned int)(sysbox.lz/interaction_dist);
	if (x_box_nb == 0) {
		x_box_nb = 1;
	}
	if (y_box_nb == 0) {
		y_box_nb = 1;
	}
	if (z_box_nb == 0) {
		z_box_nb = 1;
	}
	if (x_box_nb < 4 && y_box_nb < 4 && z_box_nb < 4) { // boxing useless: a neighborhood is the whole system
		x_box_nb = 1;
		y_box_nb = 1;
		z_box_nb = 1;
		yz_box_nb = 1;

		_is_boxed = false;
		box_xsize = sysbox.lx;
		box_ysize = sysbox.ly;
		box_zsize = sysbox.lz;
		box_nb = 1;
		auto it = Boxes.insert(new Box());
		Box* const b = (*it.first);
		b->setPosition({0, 0, 0});
		TopBottomBoxes.insert(b);
		box_labels.push_back(b);
		b->type = 1;
		
	} else {
		yz_box_nb = y_box_nb*z_box_nb;
		_is_boxed = true;
		box_xsize = sysbox.lx/x_box_nb;
		box_ysize = sysbox.ly/y_box_nb;
		box_zsize = sysbox.lz/z_box_nb;
		int m1p1[] = {-1, 1};
		for (int a : m1p1) {
			for (int b : m1p1) {
				auto far_corner = 1.4999999*vec3d(a*box_xsize, b*box_ysize, box_zsize);
				top_probing_positions.push_back(far_corner);
				top_probing_positions.push_back(far_corner-vec3d(a*box_xsize, 0, 0));
				top_probing_positions.push_back(far_corner-vec3d(a*box_xsize, b*box_ysize, 0));
				top_probing_positions.push_back(far_corner-vec3d(0, b*box_ysize, 0));
				far_corner = 1.4999999*vec3d(a*box_xsize, b*box_ysize, -box_zsize);
				bottom_probing_positions.push_back(far_corner);
				bottom_probing_positions.push_back(far_corner-vec3d(a*box_xsize, 0, 0));
				bottom_probing_positions.push_back(far_corner-vec3d(a*box_xsize, b*box_ysize, 0));
				bottom_probing_positions.push_back(far_corner-vec3d(0, b*box_ysize, 0));
			}
		}
		box_nb = x_box_nb*yz_box_nb;
		allocateBoxes();
		// give them their position
		positionBoxes();
		// tell them their neighbors
		assignNeighbors(lees);
	}
	std::cout << " [ok]" << std::endl;
}

void SimpleShearBoxSet::positionBoxes()
{
	// position boxes
	auto it = Boxes.begin();
	for (unsigned int ix=0; ix<x_box_nb; ix++) {
		for (unsigned int iy=0; iy<y_box_nb; iy++) {
			for (unsigned int iz=0; iz<z_box_nb; iz++) {
				Box* const bx = (*it);
				bx->type = 1; ///@@@@@@@@@@@@
				bx->setPosition({box_xsize*(ix+0.5), box_ysize*(iy+0.5), box_zsize*(iz+0.5)}); // the center of the box
				unsigned int label = ix*yz_box_nb+iy*z_box_nb+iz;
				box_labels[label] = bx;
				if (iz == 0 && iz < z_box_nb-1) {// bottom box
					BottomBoxes.insert(bx);
				}
				if (iz == z_box_nb-1 && iz > 0) {//top box
					TopBoxes.insert(bx);
				}
				if (iz == 0 && iz == z_box_nb-1) {// bottom box
					TopBottomBoxes.insert(bx);
				}
				if (iz > 0 && iz < z_box_nb-1) {// bulk box
					BulkBoxes.insert(bx);
				}
				it++;
			}
		}
	}
}

void SimpleShearBoxSet::assignNeighborsBulk(const BC::LeesEdwardsBC &lees)
{
	for (auto& bx : BulkBoxes) {
		auto pos = bx->getPosition();
		vec3d delta;
		int m10p1[] = {-1, 0, 1};
		for (const auto& a : m10p1) {
			delta.x = a*box_xsize;
			for (const auto& b : m10p1) {
				delta.y = b*box_ysize;
				for (const auto& c : m10p1) {
					delta.z = c*box_zsize;
					bx->addStaticNeighbor(whichBox(lees.periodized(pos+delta)));
				}
			}
		}
	}
}

void SimpleShearBoxSet::assignNeighborsBottom(const BC::LeesEdwardsBC &lees)
{
	for (auto& bx : BottomBoxes) {
		auto pos = bx->getPosition();
		vec3d delta;
		// boxes  at same level and above first: these are fixed once and for all in the simulation
		int m10p1[] = {-1, 0, 1};
		int p10[] = {0, 1};
		for (const auto& a : m10p1) {
			delta.x = a*box_xsize;
			for (const auto& b : m10p1) {
				delta.y = b*box_ysize;
				for (const auto& c : p10) {
					delta.z = c*box_zsize;
					bx->addStaticNeighbor(whichBox(lees.periodized(pos+delta)));
				}
			}
		}
		for (const auto& delta_prob : bottom_probing_positions) {
			bx->addMovingNeighbor(whichBox(lees.periodized(pos+delta_prob)));
		}
	}
}

void SimpleShearBoxSet::assignNeighborsTop(const BC::LeesEdwardsBC &lees)
{
	for (auto& bx : TopBoxes) {
		auto pos = bx->getPosition();
		vec3d delta;
		// boxes  at same level and bottom first: these are fixed once and for all in the simulation
		int m10p1[] = {-1, 0, 1};
		int m10[] = {-1, 0};
		for (const auto & a : m10p1) {
			delta.x = a*box_xsize;
			for (const auto& b : m10p1) {
				delta.y = b*box_ysize;
				for (const auto& c : m10) {
					delta.z = c*box_zsize;
					bx->addStaticNeighbor(whichBox(lees.periodized(pos+delta)));
				}
			}
		}
		for (const auto& delta_prob : top_probing_positions) {
			bx->addMovingNeighbor(whichBox(lees.periodized(pos+delta_prob)));
		}
	}
}

void SimpleShearBoxSet::assignNeighborsTopBottom(const BC::LeesEdwardsBC &lees)
{
	for (auto& bx : TopBottomBoxes) {
		auto pos = bx->getPosition();
		vec3d delta;
		// boxes at same level first: these are fixed once and for all in the simulation
		int m10p1[] = {-1, 0, 1};
		for (const auto& a : m10p1) {
			delta.x = a*box_xsize;
			for (const auto& b : m10p1) {
				delta.y = b*box_ysize;
				delta.z = 0;
				bx->addStaticNeighbor(whichBox(lees.periodized(pos+delta)));
			}
		}
		for (const auto& delta_prob : top_probing_positions) {
			bx->addMovingNeighbor(whichBox(lees.periodized(pos+delta_prob)));
		}
		for (const auto& delta_prob : bottom_probing_positions) {
			bx->addMovingNeighbor(whichBox(lees.periodized(pos+delta_prob)));
		}
	}
}


void SimpleShearBoxSet::assignNeighbors(const BC::LeesEdwardsBC &lees)
{
	// bulk boxes
	assignNeighborsBulk(lees);
	// bottom boxes
	assignNeighborsBottom(lees);
	// top boxes
	assignNeighborsTop(lees);
	// top/bottom boxes
	assignNeighborsTopBottom(lees);
}

/*****
 UpdateNeighbors()
 
 At each time step, we need to check if neighborhood on top and bottom boxes have changed.
 Bulk boxes do not need to be updated.
 *****/
void SimpleShearBoxSet::updateNeighbors(const BC::LeesEdwardsBC &lees)
{
	/**
	 \brief Update the neighbors of top and bottom boxes have changed.
	 
	 To be called when the boundary conditions have changed.
	 **/
	
	for (auto& bx : TopBoxes) {
		bx->resetMovingNeighbors();
		auto pos = bx->getPosition();
		for (const auto& delta_prob : top_probing_positions) {
			bx->addMovingNeighbor(whichBox(lees.periodized(pos+delta_prob)));
		}
	}
	
	for (auto& bx : BottomBoxes) {
		bx->resetMovingNeighbors();
		auto pos = bx->getPosition();
		for (const auto& delta_prob : bottom_probing_positions) {
			bx->addMovingNeighbor(whichBox(lees.periodized(pos+delta_prob)));
		}
	}
	
	for (auto& bx : TopBottomBoxes) {
		bx->resetMovingNeighbors();
		auto pos = bx->getPosition();
		for (const auto& delta_prob : top_probing_positions) {
			bx->addMovingNeighbor(whichBox(lees.periodized(pos+delta_prob)));
		}
		for (const auto& delta_prob : bottom_probing_positions) {
			bx->addMovingNeighbor(whichBox(lees.periodized(pos+delta_prob)));
		}
	}
}


//public methods
void SimpleShearBoxSet::update(const BC::LeesEdwardsBC &lees)
{
	if (is_boxed()) {
		updateNeighbors(lees);
	}
	for (const auto& bx : Boxes) {
		bx->buildNeighborhoodContainer();
	}
}

void SimpleShearBoxSet::box(int i, vec3d position)
{
	Box* b = whichBox(position);
	addToBox(i, b);
}

}