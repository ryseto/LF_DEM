#include <stdexcept>
#include <sstream>
#include "BoxSet.h"
using namespace std;

namespace Boxing {

void BoxSet::allocateBoxes()
{
	for (unsigned int i=0; i<box_nb; i++) {
		Boxes.insert(new Box());
	}
	box_labels.resize(box_nb);
}

bool BoxSet::is_boxed() const
{
	return _is_boxed;
}

Box* BoxSet::whichBox(unsigned int box_label) const
{
	if (box_label >= box_labels.size()) {
		ostringstream error_str;
		throw runtime_error(error_str.str());
	}
	return box_labels[box_label];
}

unsigned int BoxSet::whichBoxLabel(const vec3d &pos) const
{
	unsigned int ix = (unsigned int)(pos.x/box_xsize);
	unsigned int iy;
	if (box_ysize == 0) {
		iy = 0;
	} else {
		iy = (unsigned int)(pos.y/box_ysize);
	}
	unsigned int iz = (unsigned int)(pos.z/box_zsize);
	unsigned int label = ix*yz_box_nb+iy*z_box_nb+iz;
	if (label >= box_labels.size()) {
		// We first make sure it's not a rounding issue
		if (std::nextafter(ix*box_xsize, (ix+1u)*box_xsize) > pos.x) {
			ix--;
		}
		if (std::nextafter(iy*box_ysize, (iy+1u)*box_ysize) > pos.y) {
			iy--;
		}
		if (std::nextafter(iz*box_zsize, (iz+1u)*box_zsize) > pos.z) {
			iz--;
		}
		label = ix*yz_box_nb+iy*z_box_nb+iz;
		if (label >= box_labels.size()) {
			ostringstream error_str;
			error_str << " BoxSet: trying to box position out of boundaries \"" << pos	<< "\"" << endl;
			error_str << " pos " << pos << endl;
			error_str << " box size " << x_box_nb*box_xsize << " " << y_box_nb*box_ysize << " " << z_box_nb*box_zsize << endl;
			throw runtime_error(error_str.str());
		}
	}
	return label;
}

Box* BoxSet::whichBox(const vec3d &pos) const
{
	return box_labels[whichBoxLabel(pos)];
}

void BoxSet::addToBox(int i, Box* b)
{
	if (b != boxMap[i]) {
		b->add(i);
		if (boxMap[i] != NULL) {
			boxMap[i]->remove(i);
		}
		boxMap[i] = b;
	}
}

const vector<int>& BoxSet::neighborhood(int i) const
{
	return (boxMap[i])->getNeighborhoodContainer();
}

void BoxSet::printBoxNetwork() const
{
	for (const auto& bx : Boxes) {
		const auto& neighbors = bx->getNeighborBox();
		for (const auto& neighbor_box : neighbors) {
			cerr << " "  << neighbors.size() << " " << bx->getPosition() << " ";
			cerr << neighbor_box->getPosition() << endl;
		}
	}
}

void BoxSet::printBoxContainers() const
{
	for (const auto& bx : Boxes) {
		for (const auto& j : bx->getContainer()) {
			cerr << bx->getPosition() << " " << j << endl;
		}
	}
}

void BoxSet::printNeighborhoodContainers() const
{
	for (const auto& bx : Boxes) {
		for (const auto& j : bx->getNeighborhoodContainer()) {
			cerr << bx->getPosition() << " " << j << endl;
		}
	}
}

void BoxSet::printBoxMap() const
{
	for (unsigned i=0; i<boxMap.size(); i++) {
		cerr << i << " " << boxMap[i]->getPosition() << endl;
	}
}

}