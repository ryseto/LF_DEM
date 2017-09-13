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
	origin_ext_flow.reset();
	x_box_nb = (unsigned int)(sys->get_lx()/interaction_dist);
	y_box_nb = (unsigned int)(sys->get_ly()/interaction_dist);
	z_box_nb = (unsigned int)(sys->get_lz()/interaction_dist);
	if (x_box_nb == 0) {
		x_box_nb = 1;
	}
	if (y_box_nb == 0) {
		y_box_nb = 1;
	}
	if (z_box_nb == 0) {
		z_box_nb = 1;
	}
	yz_box_nb = y_box_nb*z_box_nb;
	if (x_box_nb < 4 && y_box_nb < 4 && z_box_nb < 4) { // boxing useless: a neighborhood is the whole system
		_is_boxed = false;
		box_xsize = sys->get_lx();
		box_ysize = sys->get_ly();
		box_zsize = sys->get_lz();
		box_nb = 1;
		auto it = Boxes.insert(new Box());
		Box* const b = (*it.first);
		b->setPosition({0, 0, 0});
		TopBottomBoxes.insert(b);
		box_labels.push_back(b);
		b->type = 1;
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
		box_nb = x_box_nb*yz_box_nb;
		allocateBoxes();
		// give them their position
		positionBoxes();
		// tell them their neighbors
		assignNeighbors();
	}
	cout << " [ok]" << endl;
}

void BoxSet::initExtFlow(double interaction_dist,
						 const vec3d origin,
						 System *sys_)
{
	string indent = "  BoxSet::\t";
	cout << indent << "Setting up Cell List System ... \n";
	sys = sys_;
	boxMap = new Box* [sys->get_np()];
	for (int i=0; i<sys->get_np(); i++) {
		boxMap[i] = NULL;
	}
	origin_ext_flow = origin;
	x_box_nb = (unsigned int)(sys->get_lx_ext_flow()/interaction_dist);
	y_box_nb = (unsigned int)(sys->get_ly_ext_flow()/interaction_dist);
	z_box_nb = (unsigned int)(sys->get_lz_ext_flow()/interaction_dist);
	if (x_box_nb == 0) {
		x_box_nb = 1;
	}
	if (y_box_nb <= 3) {
		y_box_nb = 1;
	}
	if (z_box_nb == 0) {
		z_box_nb = 1;
	}
	yz_box_nb = y_box_nb*z_box_nb;
	cout << indent << "nb : " << x_box_nb << ' ' << y_box_nb << ' ' << z_box_nb << endl;
	_is_boxed = true;
	box_xsize = sys->get_lx_ext_flow()/x_box_nb;
	box_ysize = sys->get_ly_ext_flow()/y_box_nb;
	box_zsize = sys->get_lz_ext_flow()/z_box_nb;
	box_nb = x_box_nb*yz_box_nb;
	allocateBoxes();
	// give them their position
	positionBoxes();
	// tell them their neighbors
	assignNeighborsStaticExtFlow();
	next_labels.resize(8);
	next_labels[0] = (int)(-yz_box_nb-1);
	next_labels[1] = -1;
	next_labels[2] = (int)( yz_box_nb-1);
	next_labels[3] = (int)(-yz_box_nb);
	next_labels[4] = (int)( yz_box_nb);
	next_labels[5] = (int)(-yz_box_nb+1);
	next_labels[6] = 1;
	next_labels[7] = (int)( yz_box_nb+1);
	box_corners.resize(4);
	box_corners[0] = vec3d(-0.5*box_xsize, 0, -0.5*box_zsize);
	box_corners[1] = vec3d( 0.5*box_xsize, 0, -0.5*box_zsize);
	box_corners[2] = vec3d( 0.5*box_xsize, 0,  0.5*box_zsize);
	box_corners[3] = vec3d(-0.5*box_xsize, 0,  0.5*box_zsize);

	next_of_periodicimage.resize(4);
	if (sys->twodimension) {
		nb_movingneighbor_part = 3;
	} else {
		nb_movingneighbor_part = 11;
	}
	// left bottom
	next_of_periodicimage[0].resize(nb_movingneighbor_part);
	next_of_periodicimage[0][0] = -1;
	next_of_periodicimage[0][1] = (int)(-yz_box_nb-1);
	next_of_periodicimage[0][2] = (int)(-yz_box_nb);
	// right bottom
	next_of_periodicimage[1].resize(nb_movingneighbor_part);
	next_of_periodicimage[1][0] = (int)(yz_box_nb);
	next_of_periodicimage[1][1] = (int)(yz_box_nb-1);
	next_of_periodicimage[1][2] = (int)(-1);
	// right top
	next_of_periodicimage[2].resize(nb_movingneighbor_part);
	next_of_periodicimage[2][0] = 1;
	next_of_periodicimage[2][1] = (int)(yz_box_nb+1);
	next_of_periodicimage[2][2] = (int)(yz_box_nb);
	// left top
	next_of_periodicimage[3].resize(nb_movingneighbor_part);
	next_of_periodicimage[3][0] = (int)(-yz_box_nb);
	next_of_periodicimage[3][1] = (int)(-yz_box_nb+1);
	next_of_periodicimage[3][2] = 1;
	if (sys->twodimension == false) {
		for (int j=0; j<4; j++) {
			next_of_periodicimage[j][3] =  (int)( z_box_nb);
			next_of_periodicimage[j][7] =  (int)(-z_box_nb);
			for (int k=0; k<3; k++) {
				next_of_periodicimage[j][4+k] = next_of_periodicimage[j][k]+(int)z_box_nb;// 4, 5, 6
				next_of_periodicimage[j][8+k] = next_of_periodicimage[j][k]-(int)z_box_nb;// 8, 9, 10
			}
		}
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
	if (!sys->ext_flow) {
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
	} else {
		auto it = Boxes.begin();
		for (unsigned int ix=0; ix<x_box_nb; ix++) {
			for (unsigned int iy=0; iy<y_box_nb; iy++) {
				for (unsigned int iz=0; iz<z_box_nb; iz++) {
					Box* const bx = (*it);
					// The origin of extensional flow is (0, 0, 0)
					bx->setPosition({box_xsize*(ix+0.5)-origin_ext_flow.x,
						box_ysize*(iy+0.5), box_zsize*(iz+0.5)-origin_ext_flow.z}); // the center of the box
					unsigned int label = ix*yz_box_nb+iy*z_box_nb+iz;
					box_labels[label] = bx;
					it++;
				}
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

void BoxSet::assignNeighborsStaticExtFlow()
{
	vec3d delta;
	int bmax = sys->twodimension ? 0 : 1;
	for (auto& bx : Boxes) {
		auto pos = bx->getPosition()+origin_ext_flow;
		for (int b=-bmax; b<=bmax; b++) {
			delta.y = b*box_ysize;
			for (int a=-1; a<=1; a++) {
				delta.x = a*box_xsize;
				for (int c=-1; c<=1; c++) {
					if ((a == 0 && b == 0 && c == 0) == false) {
						delta.z = c*box_zsize;
						vec3d pos_next = pos+delta;
						if (pos_next.x >= 0
							&& pos_next.z >= 0
							&& pos_next.x < sys->get_lx_ext_flow()
							&& pos_next.z < sys->get_lz_ext_flow()) {
							if (pos_next.y < 0) {
								pos_next.y += sys->get_ly_ext_flow();
							} else if (pos_next.y > sys->get_ly_ext_flow()) {
								pos_next.y -= sys->get_ly_ext_flow();
							}
							bx->addStaticNeighbor(whichBox(pos_next));
						}
					}
				}
			}
		}
	}
}

void BoxSet::assignNeighobrsDynamicExtFlow(Box& bx, const vec3d &pos)
{
	// This is for extensional flow
	// We can optimize this part more.
	//	int bmax = sys->twodimension ? 0 : 1; //
	int iy = 0;
	if (sys->twodimension == false) {
		iy = (unsigned int)(pos.y/box_ysize);
	}
	for (int j=0; j<4; j++) {
		vec3d pos_corner = pos+box_corners[j];
		if (pos_corner.x >= 0 && pos_corner.z >= 0
			&& pos_corner.x < sys->get_lx_ext_flow()
			&& pos_corner.z < sys->get_lz_ext_flow()) {
			unsigned int bl = whichBoxLabel(pos_corner);
			Box* bx_next = box_labels[bl];
			if (bx_next->type != 0) {
				bx.addMovingNeighbor(bx_next);
			}
			for (int k=0; k<nb_movingneighbor_part; k++) {
				int bl_next_pi = bl+next_of_periodicimage[j][k]; // k = 0,1,2
				if (k >= 3) {
					if (iy == y_box_nb-1) {
						if (k < 7) {
							// next_of_periodicimage[j][5+k] = next_of_periodicimage[j][k]+z_box_nb; // k = 3,4,5,6
							bl_next_pi -= yz_box_nb;
						}
					} else if (iy == 0) {
						if (k >= 7) {
							// next_of_periodicimage[j][8+k] = next_of_periodicimage[j][k]-z_box_nb;  // k=7, 8, 9, 10
							bl_next_pi += yz_box_nb;
						}
					}
				}
				if (bl_next_pi < box_nb) {
					bx_next = box_labels[bl_next_pi];
					if (bx_next->type != 0) {
						bx.addMovingNeighbor(bx_next);
					}
				}
			}
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

void BoxSet::updateNeighborsExtFlow()
{
	vec3d pos_pd;
	vec3d pos;
	for (auto& bx : Boxes) {
		if (bx->type != 0 && bx->type_neighborhood >= 2) {
			bx->resetMovingNeighbors();
			pos = bx->getPosition()+origin_ext_flow;
			if (bx->type_neighborhood == 2) {          // left
				pos_pd = pos + sys->box_axis1;         //--->right
				assignNeighobrsDynamicExtFlow(*bx, pos_pd);
			} else if (bx->type_neighborhood == 3) {   // right
				pos_pd = pos - sys->box_axis1;         //--->left
				assignNeighobrsDynamicExtFlow(*bx, pos_pd);
			} else if (bx->type_neighborhood == 4) {   // bot
				pos_pd = pos + sys->box_axis2;         //--->top
				assignNeighobrsDynamicExtFlow(*bx, pos_pd);
			} else if (bx->type_neighborhood == 8) {   // top
				pos_pd = pos - sys->box_axis2;         //--->bottom
				assignNeighobrsDynamicExtFlow(*bx, pos_pd);
			} else if (bx->type_neighborhood == 6) {   // left bottom
				pos_pd = pos + sys->box_axis1;         //--->right
				assignNeighobrsDynamicExtFlow(*bx, pos_pd);
				pos_pd = pos + sys->box_axis2;         //--->top
				assignNeighobrsDynamicExtFlow(*bx, pos_pd);
				pos_pd = pos + sys->box_diagonal_6_11; //--->right top
				assignNeighobrsDynamicExtFlow(*bx, pos_pd);
			} else if (bx->type_neighborhood == 7) {   // right bottom
				pos_pd = pos - sys->box_axis1;         //--->left
				assignNeighobrsDynamicExtFlow(*bx, pos_pd);
				pos_pd = pos + sys->box_axis2;         //--->top
				assignNeighobrsDynamicExtFlow(*bx, pos_pd);
				pos_pd = pos + sys->box_diagonal_7_10; //--->left top
				assignNeighobrsDynamicExtFlow(*bx, pos_pd);
			} else if (bx->type_neighborhood == 10) {  // left top
				pos_pd = pos + sys->box_axis1;         //--->right
				assignNeighobrsDynamicExtFlow(*bx, pos_pd);
				pos_pd = pos - sys->box_axis2;         //--->bottom
				assignNeighobrsDynamicExtFlow(*bx, pos_pd);
				pos_pd = pos - sys->box_diagonal_7_10; //--->right bot
				assignNeighobrsDynamicExtFlow(*bx, pos_pd);
			} else if (bx->type_neighborhood == 11) {  // right top
				pos_pd = pos - sys->box_axis1;         //--->left
				assignNeighobrsDynamicExtFlow(*bx, pos_pd);
				pos_pd = pos - sys->box_axis2;         //--->botom
				assignNeighobrsDynamicExtFlow(*bx, pos_pd);
				pos_pd = pos - sys->box_diagonal_6_11; //--->left bot
				assignNeighobrsDynamicExtFlow(*bx, pos_pd);
			}
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

void BoxSet::updateExtFlow()
{
	double l1 = sys->box_axis1.norm();
	double l2 = sys->box_axis2.norm();
	double delta = 0.02;
	vec3d delta_n1 = delta*sys->box_axis1/l1;
	vec3d delta_n2 = delta*sys->box_axis2/l2;
	for (auto& bx : Boxes) {
		bx->resetTypes();
	}
	set <unsigned int> active_box;
	set <unsigned int> left_edge;
	set <unsigned int> bottom_edge;
	/* Top edge */
	vec3d pos = origin_ext_flow+sys->box_axis2;
	int jmax = l1/delta;
	for (int j=0; j<2*jmax; j++) {
		unsigned int box_label = whichBoxLabel(pos);
		if (box_labels[box_label]->type == 0) {
			if (j <= jmax ) {
				active_box.insert(box_label);
			}
			box_labels[box_label]->type +=8; // Top
		}
		pos += delta_n1;
		if (pos.x >= sys->get_lx_ext_flow()) {
			break;
		}
	}
	/* Bottom edge */
	pos = origin_ext_flow+sys->box_axis1;
	for (int j=0; j<2*jmax; j++) {
		unsigned int box_label = whichBoxLabel(pos);
		if (box_labels[box_label]->type == 0) {
			if (j <= jmax ) {
				active_box.insert(box_label);
				bottom_edge.insert(box_label);
			}
			box_labels[box_label]->type += 4; // bottom
		}
		pos -= delta_n1;
		if (pos.x <= 0) {
			break;
		}
	}
	/* Left edge */
	pos = origin_ext_flow+sys->box_axis2;
	jmax = l2/delta;
	for (int j=0; j<2*jmax; j++) {
		unsigned int box_label = whichBoxLabel(pos);
		if (box_labels[box_label]->type == 0
			|| box_labels[box_label]->type == 8
			|| box_labels[box_label]->type == 4) {
			box_labels[box_label]->type += 2; // left
			if (j <= jmax ) {
				left_edge.insert(box_label);
				active_box.insert(box_label);
			}
		}
		pos -= delta_n2;
		if (pos.z <= 0) {
			break;
		}
	}
	/* Right edge */
	pos = origin_ext_flow+sys->box_axis1;
	for (int j=0; j<2*jmax; j++) {
		unsigned int box_label = whichBoxLabel(pos);
		if (box_labels[box_label]->type == 0
			|| box_labels[box_label]->type == 4
			|| box_labels[box_label]->type == 8) {
			if (j <= jmax ) {
				active_box.insert(box_label);
			}
			box_labels[box_label]->type += 3; // right
		}
		pos += delta_n2;
		if (pos.z >= sys->get_lz_ext_flow()) {
			break;
		}
	}
	/* Fill Bulk part from left edge */
	for (auto bl : left_edge) {
		if (box_labels[bl]->type == 2 || box_labels[bl]->type == 6) {
			unsigned int bl_right = labelRight(bl);
			while (box_labels[bl_right]->type == 0) {
				box_labels[bl_right]->type = 1;
				active_box.insert(bl_right);
				bl_right = labelRight(bl_right);
			}
		}
	}
	/* Fill Bulk part from bottom */
	for (auto bl : bottom_edge) {
		if (box_labels[bl]->type == 4 || box_labels[bl]->type == 6) {
			unsigned int bl_up = labelUp(bl);
			while (box_labels[bl_up]->type == 0) {
				box_labels[bl_up]->type = 1;
				active_box.insert(bl_up);
				bl_up = labelUp(bl_up);
			}
		}
	}
	checkNeighobrType(active_box);
	if (!sys->twodimension) {
		for (auto bl: active_box) {
			unsigned int bl_y = bl+z_box_nb;
			unsigned int bl_y_max = bl+yz_box_nb;
			while (bl_y < bl_y_max) {
				//cerr << box_labels[bl_y]->getYPosition() << ' ' ;
				box_labels[bl_y]->type = box_labels[bl]->type;
				box_labels[bl_y]->type_neighborhood = box_labels[bl]->type_neighborhood;
				bl_y += z_box_nb;
			}
			//cerr << endl;
		}
	}
	updateNeighborsExtFlow();
	for (const auto& bx : Boxes) {
		bx->buildNeighborhoodContainer();
	}
}

void BoxSet::checkNeighobrType(set <unsigned int> &box_labels_)
{
	for (auto bl: box_labels_) {
		bool check_bulk = true;
		int check_cooner = 0; // 6 7 10 11
		int check_sides = 0;  // 2: left 3: right
		int check_topbot = 0; // 4: bot 8: top
		for (int j=0; j<8; j++) {
			int bln = bl + next_labels[j];
			if (box_labels[bln]->type == 6 || box_labels[bln]->type == 7
				|| box_labels[bln]->type == 10 || box_labels[bln]->type == 11) {
				check_bulk = false;
				check_cooner = box_labels[bln]->type;
				break;
			} else {
				if (box_labels[bln]->type == 2 || box_labels[bln]->type == 3) {
					check_bulk = false;
					check_sides = box_labels[bln]->type;
				}
				if (box_labels[bln]->type == 4 || box_labels[bln]->type == 8) {
					check_bulk = false;
					check_topbot = box_labels[bln]->type;
				}
			}
		}
		if (check_bulk) {
			box_labels[bl]->type_neighborhood = 1;
		} else {
			if (check_cooner != 0) {
				box_labels[bl]->type_neighborhood = check_cooner;
			} else {
				//    10 - 8 --11
				//    |         |
				//    2    1    3
				//    |         |
				//    6 -- 4 ---7
				// 	box_labels[bl]->type_neighborhood = check_sides + check_topbot
				if (check_sides == 2        && check_topbot == 4) {
					box_labels[bl]->type_neighborhood = 6; //left bottom
				} else if (check_sides == 2 && check_topbot == 8) {
					box_labels[bl]->type_neighborhood = 10; // left top
				} else if (check_sides == 3 && check_topbot == 4) {
					box_labels[bl]->type_neighborhood = 7; // right bot
				} else if (check_sides == 3 && check_topbot == 8) {
					box_labels[bl]->type_neighborhood = 11; //right top
				} else {
					if (check_topbot == 0) {
						box_labels[bl]->type_neighborhood = check_sides;
					} else if (check_sides == 0) {
						box_labels[bl]->type_neighborhood = check_topbot;
					} else {
						cerr << "not expected to be here" << endl;
						exit(1);
					}
				}
			}
		}
	}
}

vec3d BoxSet::periodicDiffShift(int i, int j)
{
	int type_i = boxMap[i]->type_neighborhood;
	int type_j = boxMap[j]->type_neighborhood;
	//    10 - 8 --11
	//    |         |
	//    2    1    3
	//    |         |
	//    6 -- 4 ---7
	vec3d pd_shift(0,0,0);
	if (type_i == 1 || type_j == 1 || type_i == type_j) {
		// When one of them is in bulk, they interact within the main frame.
		// When both of them are in the same region, they interact within the main frame.
		return pd_shift;
	}
	if (type_i == 2) { // left
		if (type_j == 3 || type_j== 7 || type_j == 11) { // including right
			pd_shift = -sys->box_axis1;
		}
	} else if (type_i == 3) { // right
		if (type_j== 2 || type_j== 6 || type_j == 10) { // including left
			pd_shift = sys->box_axis1;
		}
	} else if (type_i == 4) { // bottom
		if (type_j == 8 || type_j== 10 || type_j == 11) { // including top
			pd_shift = -sys->box_axis2;
			return pd_shift;
		}
	} else if (type_i == 8) { // top
		if (type_j == 4 || type_j == 6 || type_j == 7) { // including bot
			pd_shift = sys->box_axis2;
		}
	} else if (type_i == 6) { // left bottom
		if (type_j == 3 || type_j == 7) { // right and right bot
			pd_shift = -sys->box_axis1;
		} else if (type_j == 8 || type_j == 10) { // top or top left
			pd_shift = -sys->box_axis2;
		} else if (type_j == 11) { //right top
			pd_shift = -sys->box_diagonal_6_11;
		}
	} else if (type_i == 7) { // right bottom
		if (type_j == 2 || type_j == 6) { // left and left bot
			pd_shift = sys->box_axis1;
		} else if (type_j == 8 || type_j == 11) {//top and top right
			pd_shift = -sys->box_axis2;
		} else if (type_j == 10) { // left top
			pd_shift = -sys->box_diagonal_7_10;
		}
	} else if (type_i == 10) { // left top
		if (type_j == 3 || type_j == 11) { //right and right top
			pd_shift = -sys->box_axis1;
		} else if (type_j == 4 || type_j == 6) { //bottom and bottom left
			pd_shift = sys->box_axis2;
		} else if (type_j == 7) { // right bot
			pd_shift = sys->box_diagonal_7_10;
		}
	} else if (type_i == 11) { // right top
		if (type_j == 2 || type_j == 10) { //left and left top
			pd_shift = sys->box_axis1;
		} else if (type_j == 4 || type_j == 7) { //bot and bot right
			pd_shift = sys->box_axis2;
		} else if (type_j == 6) { // left bot
			pd_shift = sys->box_diagonal_6_11;
		}
	}
	return pd_shift;
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
	unsigned int label = ix*yz_box_nb+iy*z_box_nb+iz;
	if (label >= box_labels.size()) {
		ostringstream error_str;
		error_str  << " BoxSet: trying to box position out of boundaries \"" << pos	<< "\"" << endl;
		throw runtime_error(error_str.str());
	}
	return box_labels[label];
}

Box* BoxSet::whichBox(unsigned int box_label)
{
	if (box_label >= box_labels.size()) {
		ostringstream error_str;
		throw runtime_error(error_str.str());
	}
	return box_labels[box_label];
}

unsigned int BoxSet::whichBoxLabel(const vec3d &pos)
{
	unsigned int ix = (unsigned int)(pos.x/box_xsize);
	unsigned int iy;
	if (sys->twodimension) {
		iy = 0;
	} else {
		iy = (unsigned int)(pos.y/box_ysize);
	}
	unsigned int iz = (unsigned int)(pos.z/box_zsize);
	unsigned int label = ix*yz_box_nb+iy*z_box_nb+iz;
	if (label >= box_labels.size()) {
		ostringstream error_str;
		error_str << " BoxSet: trying to box position out of boundaries \"" << pos	<< "\"" << endl;
		error_str << pos - origin_ext_flow << endl;
		throw runtime_error(error_str.str());
	}
	return label;
}

void BoxSet::box(int i)
{
//	Box* b = whichBox(sys->position[i]);
	Box* b = whichBox(sys->position[i]+origin_ext_flow);
	if (b != boxMap[i]) {
		b->add(i);
		if (boxMap[i] != NULL) {
			boxMap[i]->remove(i);
		}
		boxMap[i] = b;
	}
}

const vector<int>& BoxSet::neighborhood(int i)
{
	return (boxMap[i])->getNeighborhoodContainer();
}

void BoxSet::printBoxNetwork()
{
	for (const auto& bx : Boxes) {
		const auto& neighbors = bx->getNeighborBox();
		for (const auto& neighbor_box : neighbors) {
			cerr << " "  << neighbors.size() << " " << bx->getPosition() << " ";
			cerr << neighbor_box->getPosition() << endl;
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

void BoxSet::yaplotBox(std::ofstream &fout_boxing)
{
	vec3d dx(0.5*box_xsize, 0, 0);
	vec3d dz(0, 0, 0.5*box_zsize);
	vec3d dy(0, 0.01, 0);
	fout_boxing << "y 3\n";
	for (const auto& bx : Boxes) {
		//if (bx->is_active()) {
		vec3d center = bx->getPosition();
		vec3d left_bottom = center - dx - dz + dy;
		vec3d right_bottom = center + dx - dz + dy;
		vec3d right_top = center + dx + dz + dy;
		vec3d left_top = center - dx + dz + dy;
		int box_color = bx->type_neighborhood;
		//int box_color = bx->type;
		if (box_color == 0) {
			fout_boxing << "@ 0" << endl;
		} else if (box_color == 1) {
			fout_boxing << "@ 2" << endl;
		} else if (box_color == 2 || box_color == 3) {
			fout_boxing << "@ 3" << endl;
		} else if (box_color == 4 || box_color == 8) {
			fout_boxing << "@ 5" << endl;
		}  else if (box_color == 6) {
			fout_boxing << "@ 6" << endl;
		} else if (box_color == 7) {
			fout_boxing << "@ 6" << endl;
		} else if (box_color == 10) {
			fout_boxing << "@ 6" << endl;
		} else if (box_color == 11) {
			fout_boxing << "@ 6" << endl;
		} else {
			cerr << "unexpected bx type" << endl;
			cerr << box_color << endl;
			exit(1);
		}
		if (sys->twodimension) {
			fout_boxing << "p 4 "<< left_bottom;
			fout_boxing << ' ' << right_bottom;
			fout_boxing << ' ' << right_top;
			fout_boxing << ' ' << left_top;
			fout_boxing << endl;
		} else {
			if (box_color != 0) {
				fout_boxing << "l "<< left_bottom << ' ' << right_bottom << endl;
				fout_boxing << "l " << right_bottom << ' '<< right_top << endl;
				fout_boxing << "l " << right_top << ' '<< left_top << endl;
				fout_boxing << "l " << left_top << ' '<< left_bottom << endl;
			}
		}
	}
	if (0) {
		fout_boxing << "r 1\n";
		for (const auto& bx : Boxes) {
			if (bx->type != 0 && bx->type_neighborhood >= 1) {
				if (bx->type_neighborhood == 1) {
					fout_boxing << "@ 2" << endl;
				} else if (bx->type_neighborhood == 2 || bx->type_neighborhood == 3) {
					fout_boxing << "@ 3" << endl;
				} else if (bx->type_neighborhood == 4 || bx->type_neighborhood == 8) {
					fout_boxing << "@ 4" << endl;
				} else {
					fout_boxing << "@ 0" << endl;
				}
				fout_boxing << "c " <<  bx->getPosition() << endl;
			}
		}
	}
}
