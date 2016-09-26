#include <stdexcept>
#include <sstream>
#include "BoxSet.h"
#include "System.h"
using namespace std;

void BoxSet::init(double interaction_dist, System* sys_)
{
	string indent = "  BoxSet::\t";
	cout << indent << "Setting up Cell List System ... ";
	sys = sys_;
	boxMap = new Box* [sys->get_np()];
	for (int i=0; i<sys->get_np(); i++) {
		boxMap[i] = NULL;
	}
	double xratio = sys->get_lx()/interaction_dist;
	double yratio = sys->get_ly()/interaction_dist;
	double zratio = sys->get_lz()/interaction_dist;
	x_box_nb = (unsigned int)xratio;
	y_box_nb = (unsigned int)yratio;
	z_box_nb = (unsigned int)zratio;
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
		_is_boxed = false;
		box_xsize = sys->get_lx();
		box_ysize = sys->get_ly();
		box_zsize = sys->get_lz();
		box_nb = 1;

		auto it = Boxes.insert(new Box());
		Box* const b = (*it.first);
		b->setPosition({0, 0, 0});
		b->is_bottom(true);
		b->is_top(true);
		TopBottomBoxes.insert(b);
		box_labels.push_back(b);
	} else {
		_is_boxed = true;
		box_xsize = sys->get_lx()/x_box_nb;
		box_ysize = sys->get_ly()/y_box_nb;
		box_zsize = sys->get_lz()/z_box_nb;
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
		box_nb = x_box_nb*y_box_nb*z_box_nb;
		allocateBoxes();
		// give them their position
		positionBoxes();
		// tell them their neighbors
		assignNeighbors();
	}
	cout << " [ok]" << endl;
}

void BoxSet::allocateBoxes()
{
	for (unsigned int i=0; i<box_nb; i++) {
		Boxes.insert(new Box());
	}
	box_labels.resize(box_nb);
}

void BoxSet::positionBoxes()
{
	// position boxes
	auto it = Boxes.begin();

	for (unsigned int ix=0; ix<x_box_nb; ix++) {
		for (unsigned int iy=0; iy<y_box_nb; iy++) {
			for (unsigned int iz=0; iz<z_box_nb; iz++) {
				Box* const bx = (*it);
				bx->setPosition({box_xsize*(ix+0.5), box_ysize*(iy+0.5), box_zsize*(iz+0.5)}); // the center of the box
				int label = ix*y_box_nb*z_box_nb+iy*z_box_nb+iz;
				box_labels[label] = bx;
				if (iz == 0 && iz < z_box_nb-1) {// bottom box
					bx->is_bottom(true);
					BottomBoxes.insert(bx);
				}
				if (iz == z_box_nb-1 && iz > 0) {//top box
					bx->is_top(true);
					TopBoxes.insert(bx);
				}
				if (iz == 0 && iz == z_box_nb-1) {// bottom box
					bx->is_bottom(true);
					bx->is_top(true);
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


void BoxSet::assignNeighborsBulk()
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
					bx->addStaticNeighbor(whichBox(sys->periodized(pos+delta)));
				}
			}
		}
	}
}

void BoxSet::assignNeighborsBottom()
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
					bx->addStaticNeighbor(whichBox(sys->periodized(pos+delta)));
				}
			}
		}
		for (const auto& delta_prob : bottom_probing_positions) {
			bx->addMovingNeighbor(whichBox(sys->periodized(pos+delta_prob)));
		}
	}
}

void BoxSet::assignNeighborsTop()
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
					bx->addStaticNeighbor(whichBox(sys->periodized(pos+delta)));
				}
			}
		}
		for (const auto& delta_prob : top_probing_positions) {
			bx->addMovingNeighbor(whichBox(sys->periodized(pos+delta_prob)));
		}
	}
}

void BoxSet::assignNeighborsTopBottom()
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
				bx->addStaticNeighbor(whichBox(sys->periodized(pos+delta)));
			}
		}

		for (const auto& delta_prob : top_probing_positions) {
			bx->addMovingNeighbor(whichBox(sys->periodized(pos+delta_prob)));
		}
		for (const auto& delta_prob : bottom_probing_positions) {
			bx->addMovingNeighbor(whichBox(sys->periodized(pos+delta_prob)));
		}
	}
}

void BoxSet::assignNeighbors()
{
	// bulk boxes
	assignNeighborsBulk();
	// bottom boxes
	assignNeighborsBottom();
	// top boxes
	assignNeighborsTop();
	// top/bottom boxes
	assignNeighborsTopBottom();
}

BoxSet::~BoxSet()
{
	Boxes.clear();
	BulkBoxes.clear();
	TopBoxes.clear();
	BottomBoxes.clear();
	TopBottomBoxes.clear();
	box_labels.clear();
	DELETE(boxMap);
}

/*****
 UpdateNeighbors()

 At each time step, we need to check if neighborhood on top and bottom boxes have changed.
 Bulk boxes do not need to be updated.
 *****/
void BoxSet::updateNeighbors()
{
	/**
	 \brief Update the neighbors of top and bottom boxes have changed.

		To be called when the boundary conditions have changed.
	 **/

	for (auto& bx : TopBoxes) {
		bx->resetMovingNeighbors();
		auto pos = bx->getPosition();
		for (const auto& delta_prob : top_probing_positions) {
			bx->addMovingNeighbor(whichBox(sys->periodized(pos+delta_prob)));
		}
	}

	for (auto& bx : BottomBoxes) {
		bx->resetMovingNeighbors();
		auto pos = bx->getPosition();
		for (const auto& delta_prob : bottom_probing_positions) {
			bx->addMovingNeighbor(whichBox(sys->periodized(pos+delta_prob)));
		}
	}

	for (auto& bx : TopBottomBoxes) {
		bx->resetMovingNeighbors();
		auto pos = bx->getPosition();
		for (const auto& delta_prob : top_probing_positions) {
			bx->addMovingNeighbor(whichBox(sys->periodized(pos+delta_prob)));
		}
		for (const auto& delta_prob : bottom_probing_positions) {
			bx->addMovingNeighbor(whichBox(sys->periodized(pos+delta_prob)));
		}
	}
}

//public methods
void BoxSet::update()
{
	if (is_boxed()) {
		updateNeighbors();
	}
	for (const auto& bx : Boxes) {
		bx->buildNeighborhoodContainer();
	}
}

bool BoxSet::is_boxed()
{
	return _is_boxed;
}

Box* BoxSet::whichBox(const vec3d &pos)
{
	unsigned int ix = (unsigned int)(pos.x/box_xsize);
	unsigned int iy;
	if (sys->twodimension) {
		iy = 0;
	} else {
		iy = (unsigned int)(pos.y/box_ysize);
	}
	unsigned int iz = (unsigned int)(pos.z/box_zsize);
	unsigned int label = ix*y_box_nb*z_box_nb+iy*z_box_nb+iz;
	if (label > box_labels.size()-1) {
		ostringstream error_str;
		error_str  << " BoxSet: trying to box position out of boundaries \"" << pos	<< "\"" << endl;
		throw runtime_error(error_str.str());
	}
	return box_labels[label];
}

void BoxSet::box(int i)
{
	Box* b = whichBox(sys->position[i]);
	if (b != boxMap[i]) {
		b->add(i);
		if (boxMap[i] != NULL) {
			boxMap[i]->remove(i);
		}
		boxMap[i] = b;
	}
}

const vector<int>& BoxSet::neighborhood(int i){
	return (boxMap[i])->getNeighborhoodContainer();
}

void BoxSet::printBoxNetwork()
{
	for (const auto& bx : Boxes) {
		const auto& neighbors = bx->getNeighborBox();
		for (const auto& neighbor_box : neighbors) {
			cerr << " "  << neighbors.size() << " " << bx->getPosition() << " ";
			cerr << neighbor_box->getPosition() << " " << bx->is_top() << " ";
			cerr << bx->is_bottom() << endl;
		}
	}
}

void BoxSet::printBoxContainers()
{
	for (const auto& bx : Boxes) {
		for (const auto& j : bx->getContainer()) {
			cerr << bx->getPosition() << " " << j << endl;
		}
	}
}

void BoxSet::printNeighborhoodContainers()
{
	for (const auto& bx : Boxes) {
		for (const auto& j : bx->getNeighborhoodContainer()) {
			cerr << bx->getPosition() << " " << j << endl;
		}
	}
}

void BoxSet::printBoxMap()
{
	for (int i=0; i<sys->get_np(); i++) {
		cerr << i << " " << boxMap[i]->getPosition() << endl;
	}
}
