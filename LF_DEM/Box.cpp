#include "Box.h"
using namespace std;

Box::~Box()
{
	_neighbors.clear();
	_moving_neighbors.clear();
}

/* reset every moving neighbor:
 * important for top/bottom boxes, as the number of moving neighbors
 * can be smaller than _moving_neigh_nb
 */
void Box::resetMovingNeighbors()
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

void Box::add(int i)
{
	container.insert(i);
}

void Box::remove(int i)
{
	container.erase(i);
}

void Box::buildNeighborhoodContainer()
{
	neighborhood_container.clear();
	size_t size = container.size();
	for (const auto& box : _neighbors) {
		size += box->getContainer().size();
	}
	neighborhood_container.resize(size);
	int j = 0;
	// own box
	for (const int& k : container) {
		neighborhood_container[j] = k;
		j++;
	}
	// neighboring boxes
	for (const auto& box : _neighbors) {
		for (const int& k : box->container) {
			neighborhood_container[j] = k;
			j++;
		}
	}
}
