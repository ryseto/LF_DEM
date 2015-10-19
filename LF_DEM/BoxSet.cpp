#include "BoxSet.h"
#include "System.h"
using namespace std;

void BoxSet::init(double interaction_dist, System *sys_)
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
	x_box_nb = (int)xratio;
	y_box_nb = (int)yratio;
	z_box_nb = (int)zratio;
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
		b->position.reset();
		b->is_bottom(true);
		b->is_top(true);
		TopBottomBoxes.insert(b);
		box_labels.push_back(b);
	} else {
		_is_boxed = true;
		box_xsize = sys->get_lx()/x_box_nb;
		box_ysize = sys->get_ly()/y_box_nb;
		box_zsize = sys->get_lz()/z_box_nb;
		int m1p1 [2] = {-1,1};
		for (int a : m1p1 ) {
			for (int b : m1p1 ) {
				vec3d far_corner = 1.4999999*vec3d(a*box_xsize,b*box_ysize,box_zsize);
				top_probing_positions.push_back(far_corner);
				top_probing_positions.push_back(far_corner-vec3d(a*box_xsize,0,0));
				top_probing_positions.push_back(far_corner-vec3d(a*box_xsize,b*box_ysize,0));
				top_probing_positions.push_back(far_corner-vec3d(0,b*box_ysize,0));

				far_corner = 1.4999999*vec3d(a*box_xsize,b*box_ysize,-box_zsize);
				bottom_probing_positions.push_back(far_corner);
				bottom_probing_positions.push_back(far_corner-vec3d(a*box_xsize,0,0));
				bottom_probing_positions.push_back(far_corner-vec3d(a*box_xsize,b*box_ysize,0));
				bottom_probing_positions.push_back(far_corner-vec3d(0,b*box_ysize,0));
			}
		}

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
	box_nb = x_box_nb*y_box_nb*z_box_nb;

	for (int i=0; i<box_nb;i++) {
		Boxes.insert(new Box());
	}
	box_labels.resize(box_nb);
}

void BoxSet::positionBoxes()
{
	if (x_box_nb > 3) {
		amax = 3;
	} else {
		amax = x_box_nb;
	}
	if (y_box_nb > 3) {
		bmax = 3;
	} else {
		bmax = y_box_nb;
	}
	if (z_box_nb > 3) {
		cmax = 3;
	} else {
		cmax = z_box_nb;
	}

	// position boxes
	auto it = Boxes.begin();

	for (int ix=0; ix<x_box_nb; ix++) {
		for (int iy=0; iy<y_box_nb; iy++) {
			for (int iz=0; iz<z_box_nb; iz++) {
				Box* const bx = (*it);
				bx->position = vec3d(box_xsize*(ix+0.5), box_ysize*(iy+0.5), box_zsize*(iz+0.5)); // the center of the box
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
	for (auto & bx : BulkBoxes) {
		vec3d pos = bx->position;
		vec3d delta;
		int m10p1 [3] = {-1,0,1};
		for (const auto& a : m10p1) {
			delta.x = a*box_xsize;
			for (const auto& b : m10p1) {
				delta.y = b*box_ysize;
				for (const auto& c : m10p1) {
					delta.z = c*box_zsize;
					bx->addStaticNeighbor(WhichBox(pos+delta));
				}
			}
		}
	}
}

void BoxSet::assignNeighborsBottom()
{
	for (auto & bx : BottomBoxes) {
		vec3d pos = bx->position;
		vec3d delta;

		// boxes  at same level and above first: these are fixed once and for all in the simulation
		int m10p1 [3] = {-1,0,1};
		int p10 [2] = {0,1};
		for (const auto & a : m10p1) {
			delta.x = a*box_xsize;
			for (const auto& b : m10p1) {
				delta.y = b*box_ysize;
				for (const auto& c : p10) {
					delta.z = c*box_zsize;
					bx->addStaticNeighbor(WhichBox(pos+delta));
				}
			}
		}

		for (const auto& delta_prob : bottom_probing_positions){
			bx->addMovingNeighbor(WhichBox(pos+delta_prob));
		}
	}
}

void BoxSet::assignNeighborsTop()
{
	for (auto & bx : TopBoxes) {
		vec3d pos = bx->position;
		vec3d delta;

		// boxes  at same level and bottom first: these are fixed once and for all in the simulation
		int m10p1 [3] = {-1,0,1};
		int m10 [2] = {-1,0};
		for (const auto & a : m10p1) {
			delta.x = a*box_xsize;
			for (const auto& b : m10p1) {
				delta.y = b*box_ysize;
				for (const auto& c : m10) {
					delta.z = c*box_zsize;
					bx->addStaticNeighbor(WhichBox(pos+delta));
				}
			}
		}

		for (const auto& delta_prob : top_probing_positions){
			bx->addMovingNeighbor(WhichBox(pos+delta_prob));
		}
	}
}


void BoxSet::assignNeighborsTopBottom()
{
	for (auto & bx : TopBottomBoxes) {
		vec3d pos = bx->position;
		vec3d delta;

		// boxes at same level first: these are fixed once and for all in the simulation
		int m10p1 [3] = {-1,0,1};
		for (const auto & a : m10p1) {
			delta.x = a*box_xsize;
			for (const auto& b : m10p1) {
				delta.y = b*box_ysize;
				delta.z = 0;
				bx->addStaticNeighbor(WhichBox(pos+delta));
			}
		}

		for (const auto& delta_prob : top_probing_positions){
			bx->addMovingNeighbor(WhichBox(pos+delta_prob));
		}
		for (const auto& delta_prob : bottom_probing_positions){
			bx->addMovingNeighbor(WhichBox(pos+delta_prob));
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
BoxSet::~BoxSet(){
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

	for (auto & bx : TopBoxes) {
		bx->reset_moving_neighbors();
		vec3d pos = bx->position;
		vec3d delta;

		for (const auto& delta_prob : top_probing_positions){
			bx->addMovingNeighbor(WhichBox(pos+delta_prob));
		}
	}

	for (auto & bx : BottomBoxes) {
		bx->reset_moving_neighbors();
		vec3d pos = bx->position;
		vec3d delta;

		for (const auto& delta_prob : bottom_probing_positions){
			bx->addMovingNeighbor(WhichBox(pos+delta_prob));
		}
	}

	for (auto & bx : TopBottomBoxes) {
		bx->reset_moving_neighbors();
		vec3d pos = bx->position;
		vec3d delta;

		for (const auto& delta_prob : top_probing_positions){
			bx->addMovingNeighbor(WhichBox(pos+delta_prob));
		}
		for (const auto& delta_prob : bottom_probing_positions){
			bx->addMovingNeighbor(WhichBox(pos+delta_prob));
		}
	}
}

//public methods
void BoxSet::update()
{
	if (is_boxed()) {
		updateNeighbors();
	}
	for (const auto & bx : Boxes) {
		bx->build_neighborhood_container();
	}

	// if(sys->get_shear_strain()>0.237){
	// 	printBoxNetwork();exit(1);
	// }
	//	printBoxContainers(); exit(1);
	//	printNeighborhoodContainers(); exit(1);
	// printBoxMap(); exit(1);
}

bool BoxSet::is_boxed()
{
	return _is_boxed;
}

Box* BoxSet::WhichBox(vec3d pos)
{
	return WhichBox(&pos);
}

Box* BoxSet::WhichBox(vec3d *pos)
{
	sys->periodize(*pos);
	int ix = (int)(pos->x/box_xsize);
	int iy;
	if (sys->twodimension) {
		iy = 0;
	} else {
		iy = (int)(pos->y/box_ysize);
	}
	int iz = (int)(pos->z/box_zsize);
	int label = ix*y_box_nb*z_box_nb+iy*z_box_nb+iz;

	return box_labels[label];
}

void BoxSet::box(int i)
{
	Box *b = WhichBox(sys->position[i]);
	if (b != boxMap[i]) {
		b->add(i);
		if (boxMap[i] != NULL) {
			boxMap[i]->remove(i);
		}
		boxMap[i] = b;
	}
}

vector <int> & BoxSet::neighborhood(int i){
	return (boxMap[i])->neighborhood_container;
}

vector<int>::iterator BoxSet::neighborhood_begin(int i)
{
	return (boxMap[i])->neighborhood_begin();
}

vector<int>::iterator BoxSet::neighborhood_end(int i)
{
	return (boxMap[i])->neighborhood_end();
}

void BoxSet::printBoxNetwork()
{
	for (const auto & bx : Boxes) {
		const auto & neighbors = bx->neighbors();
		for (const auto & neighbor_box : neighbors) {
			cerr << " "  << neighbors.size() << " " << bx->position << " ";
			cerr << neighbor_box->position << " " << bx->is_top() << " ";
			cerr << bx->is_bottom() << endl;
		}
	}
}

void BoxSet::printBoxContainers()
{
	for (const auto & bx : Boxes) {
		for(const auto& j : bx->container){
			cerr << bx->position << " " << j << endl;
		}
	}
}

void BoxSet::printNeighborhoodContainers()
{
	for (const auto & bx : Boxes) {
		for(const auto& j : bx->neighborhood_container){
			cerr << bx->position << " " << j << endl;
		}
	}
}

void BoxSet::printBoxMap()
{
	for (int i=0; i<sys->get_np(); i++) {
		cerr << i << " " << boxMap[i]->position << endl;
	}
}
