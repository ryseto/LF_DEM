#include "Box.h"
using namespace std;

Box::Box(){
	_is_top=false;
	_is_bottom=false;
	// #ifdef __vector_container__
	//   container = new vector <int>;
	// #endif
	// #ifdef __flist_container__
	//   container = new forward_list <int>;
	// #endif
	//container = new set <int>;
}

Box::~Box(){
	if (_neigh_nb > 0) {
		delete [] _neighbors;
	}
	if (_moving_neigh_nb > 0) {
		delete [] _moving_neighbors;
		delete [] _probing_positions;
	}
	//delete container;
}

/* reset every moving neighbor:
 * important for top/bottom boxes, as the number of moving neighbors
 * can be smaller than _moving_neigh_nb
 */
void Box::reset_moving_neighbors()
{
	for (int label=0; label<_moving_neigh_nb; label++) {
		_neighbors[label+_still_neigh_nb] = NULL;
		_moving_neighbors[label] = NULL;
	}
}

void Box::neigh_nb(int n, int moving_n)
{
	_neigh_nb = n;
	_moving_neigh_nb = moving_n;
	_still_neigh_nb = n-moving_n;
	if (_neigh_nb > 0) {
		_neighbors = new Box* [_neigh_nb];
	}
	if (_moving_neigh_nb > 0) {
		_moving_neighbors = new Box* [_moving_neigh_nb];
		_probing_positions = new vec3d [_moving_neigh_nb];
		_probe_nb = _moving_neigh_nb;
	}
	reset_moving_neighbors();
}

bool Box::can_be_added(int label, Box* neigh_box)
{
	if (neigh_box == this) {
		return false;
	}
	for (int i=0; i<label; i++) {
		if (neigh_box == _neighbors[i]) {
			return false;
		}
	}
	return true;
}

bool Box::neighbor(int label, Box* neigh_box)
{
	// we have to check that we don't list two times the same neighbor (for top AND bottom boxes)
	// or don't the the box as its own neighbor (for every box)
	if (can_be_added(label, neigh_box)) {
		_neighbors[label] = neigh_box;
		if (label >= _still_neigh_nb) {
			_moving_neighbors[label-_still_neigh_nb] = neigh_box;
		}
		return true;
	} else {
		return false;
	}
}

void Box::probing_positions(int label, const vec3d &pos)
{
	//  cout << label-_still_neigh_nb << endl;
	_probing_positions[label-_still_neigh_nb] = pos;
}

bool Box::moving_neighbor(int moving_label, Box* neigh_box)
{
	bool success = neighbor(moving_label+_still_neigh_nb, neigh_box);
	return success;
}

void Box::is_top(bool it)
{
	_is_top = it;
}

void Box::is_bottom(bool ib)
{
	_is_bottom = ib;
}

bool Box::is_top()
{
	return _is_top;
}

bool Box::is_bottom()
{
	return _is_bottom;
}

// void
// Box::add(int i){
// #ifdef __vector_container__
//   container->push_back(i);
// #endif
// #ifdef __flist_container__
//   container->push_front(i);
// #endif
// }

// void
// Box::remove(int i){
// #ifdef __vector_container__
//   container->push_back(i);
// #endif
// #ifdef __flist_container__
//   container->push_front(i);
// #endif
// }

void Box::add(int i)
{
	container.insert(i);
}

void Box::remove(int i)
{
	container.erase(i);
}

void Box::build_neighborhood_container()
{
	neighborhood_container.clear();
	int size = (int)container_size();
	for (int i=0; i<neigh_nb(); i++) {
		size += _neighbors[i]->container_size();
	}
	neighborhood_container.resize(size);
	int j = 0;
	// own box
	for (set<int>::iterator it = begin(); it != end(); it++) {
		neighborhood_container[j] = *it;
		j++;
	}
	// neighboring boxes
	for (int i=0; i<neigh_nb(); i++) {
		for(set<int>::iterator it = (_neighbors[i])->begin();
			it != (_neighbors[i])->end(); it++) {
			neighborhood_container[j] = *it;
			j++;
		}
	}
}
