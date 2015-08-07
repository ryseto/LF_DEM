#include "Box.h"
using namespace std;

Box::Box(){
	_is_top=false;
	_is_bottom=false;
}

Box::~Box(){
	_neighbors.clear();
	_moving_neighbors.clear();
}

/* reset every moving neighbor:
 * important for top/bottom boxes, as the number of moving neighbors
 * can be smaller than _moving_neigh_nb
 */
void Box::reset_moving_neighbors()
{
	for (const auto& box : _moving_neighbors) {
		_neighbors.erase(box);
	}
	_moving_neighbors.clear();
}

void Box::addStaticNeighbor(Box* neigh_box)
{
	if (neigh_box == this) {
		return;
	}
	_neighbors.insert(neigh_box);
}

void Box::addMovingNeighbor(Box* neigh_box)
{
	if (neigh_box == this) {
		return;
	}
	_neighbors.insert(neigh_box);
	_moving_neighbors.insert(neigh_box);
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
	size_t size = container_size();
	for (const auto& box : _neighbors) {
		size += box->container_size();
	}
	neighborhood_container.resize(size);
	int j = 0;
	// own box
	for (const int &k : container) {
		neighborhood_container[j] = k;
		j++;
	}
	// neighboring boxes
	for (const auto& box : _neighbors) {
		for (const int &k : box->container) {
			neighborhood_container[j] = k;
			j++;
		}
	}
}
