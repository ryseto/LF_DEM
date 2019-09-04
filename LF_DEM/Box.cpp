#include "Box.h"

namespace Boxing {

Box::~Box()
{
	static_neighbors.clear();
	moving_neighbors.clear();
}

/* reset every moving neighbor:
 * important for top/bottom boxes, as the number of moving neighbors
 * can be smaller than _moving_neigh_nb
 */
void Box::resetMovingNeighbors()
{
	moving_neighbors.clear();
}

void Box::addStaticNeighbor(Box* neigh_box)
{
	if (neigh_box == this) {
		return;
	}
	static_neighbors.push_back(neigh_box);
}

void Box::addMovingNeighbor(Box* neigh_box)
{
	if (neigh_box == this) {
		return;
	}
	moving_neighbors.push_back(neigh_box);
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
	if (type == 0) {
		return;
	}
	neighborhood_container.clear();
	size_t size = container.size();
	for (const auto& box : static_neighbors) {
		if (box->type != 0) {
			size += box->getContainer().size();
		}
	}
	for (const auto& box : moving_neighbors) {
		if (box->type != 0) {
			size += box->getContainer().size();
		}
	}
	neighborhood_container.resize(size);
	int j = 0;
	// own box
	for (const int& k : container) {
		neighborhood_container[j] = k;
		j++;
	}
	// static neighboring boxes
	for (const auto& box : static_neighbors) {
		if (box->type != 0) {
			for (const int& k : box->container) {
				neighborhood_container[j] = k;
				j++;
			}
		}
	}
	// moving neighboring boxes
	for (const auto& box : moving_neighbors) {
		if (box->type != 0) {
			for (const int& k : box->container) {
				neighborhood_container[j] = k;
				j++;
			}
		}
	}
}

} // namespace Boxing