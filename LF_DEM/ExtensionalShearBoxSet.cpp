#include <string>
#include "ExtensionalShearBoxSet.h"

namespace Boxing {

ExtensionalShearBoxSet::ExtensionalShearBoxSet(double interaction_dist,
											   const vec3d origin,
											   unsigned np, 
											   Geometry::box3d sysbox,
											   bool twodimension_)
{
	std::string indent = "  BoxSet::\t";
	std::cout << indent << "Setting up Cell List System ... \n";
	boxMap.resize(np);
	for (unsigned i=0; i<np; i++) {
		boxMap[i] = NULL;
	}
	origin_ext_flow = origin;
	twodimension = twodimension_;
	x_box_nb = (unsigned int)(sysbox.lx/interaction_dist);
	y_box_nb = (unsigned int)(sysbox.ly/interaction_dist);
	z_box_nb = (unsigned int)(sysbox.lz/interaction_dist);
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
	std::cout << indent << "nb : " << x_box_nb << ' ' << y_box_nb << ' ' << z_box_nb << std::endl;
	_is_boxed = true;
	box_xsize = sysbox.lx/x_box_nb;
	box_ysize = sysbox.ly/y_box_nb;
	box_zsize = sysbox.lz/z_box_nb;
	box_nb = x_box_nb*yz_box_nb;
	allocateBoxes();
	// give them their position
	positionBoxes();
	// tell them their neighbors
	assignNeighborsStaticExtFlow(sysbox);
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
	if (twodimension) {
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
	if (twodimension == false) {
		for (int j=0; j<4; j++) {
			next_of_periodicimage[j][3] =  (int)( z_box_nb);
			next_of_periodicimage[j][7] =  (int)(-z_box_nb);
			for (int k=0; k<3; k++) {
				next_of_periodicimage[j][4+k] = next_of_periodicimage[j][k]+(int)z_box_nb;// 4, 5, 6
				next_of_periodicimage[j][8+k] = next_of_periodicimage[j][k]-(int)z_box_nb;// 8, 9, 10
			}
		}
	}
	std::cout << " [ok]" << std::endl;
}

void ExtensionalShearBoxSet::positionBoxes()
{
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

void ExtensionalShearBoxSet::assignNeighborsStaticExtFlow(Geometry::box3d sysbox)
{
	vec3d delta;
	int bmax = twodimension ? 0 : 1;
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
							&& pos_next.x < sysbox.lx
							&& pos_next.z < sysbox.lz) {
							if (pos_next.y < 0) {
								pos_next.y += sysbox.ly;
							} else if (pos_next.y > sysbox.ly) {
								pos_next.y -= sysbox.ly;
							}
							bx->addStaticNeighbor(whichBox(pos_next));
						}
					}
				}
			}
		}
	}
}

void ExtensionalShearBoxSet::assignNeighobrsDynamicExtFlow(Box& bx, const vec3d &pos, const Geometry::box3d &sysbox)
{
	// This is for extensional flow
	// We can optimize this part more.
	//	int bmax = sys->twodimension ? 0 : 1; //
	int iy = 0;
	if (twodimension == false) {
		iy = (unsigned int)(pos.y/box_ysize);
	}
	for (int j=0; j<4; j++) {
		vec3d pos_corner = pos+box_corners[j];
		if (pos_corner.x >= 0 && pos_corner.z >= 0
			&& pos_corner.x < sysbox.lx
			&& pos_corner.z < sysbox.lz) {
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


void ExtensionalShearBoxSet::updateNeighborsExtFlow(const struct Geometry::box3d &sysbox,
													const struct BC::ExtFlowAxis &ax)
{
	vec3d pos_pd;
	vec3d pos;
	for (auto& bx : Boxes) {
		if (bx->type != 0 && bx->type_neighborhood >= 2) {
			bx->resetMovingNeighbors();
			pos = bx->getPosition()+origin_ext_flow;
			if (bx->type_neighborhood == 2) {          // left
				pos_pd = pos + ax.box_axis1;         //--->right
				assignNeighobrsDynamicExtFlow(*bx, pos_pd, sysbox);
			} else if (bx->type_neighborhood == 3) {   // right
				pos_pd = pos - ax.box_axis1;         //--->left
				assignNeighobrsDynamicExtFlow(*bx, pos_pd, sysbox);
			} else if (bx->type_neighborhood == 4) {   // bot
				pos_pd = pos + ax.box_axis2;         //--->top
				assignNeighobrsDynamicExtFlow(*bx, pos_pd, sysbox);
			} else if (bx->type_neighborhood == 8) {   // top
				pos_pd = pos - ax.box_axis2;         //--->bottom
				assignNeighobrsDynamicExtFlow(*bx, pos_pd, sysbox);
			} else if (bx->type_neighborhood == 6) {   // left bottom
				pos_pd = pos + ax.box_axis1;         //--->right
				assignNeighobrsDynamicExtFlow(*bx, pos_pd, sysbox);
				pos_pd = pos + ax.box_axis2;         //--->top
				assignNeighobrsDynamicExtFlow(*bx, pos_pd, sysbox);
				pos_pd = pos + ax.box_diagonal_6_11; //--->right top
				assignNeighobrsDynamicExtFlow(*bx, pos_pd, sysbox);
			} else if (bx->type_neighborhood == 7) {   // right bottom
				pos_pd = pos - ax.box_axis1;         //--->left
				assignNeighobrsDynamicExtFlow(*bx, pos_pd, sysbox);
				pos_pd = pos + ax.box_axis2;         //--->top
				assignNeighobrsDynamicExtFlow(*bx, pos_pd, sysbox);
				pos_pd = pos + ax.box_diagonal_7_10; //--->left top
				assignNeighobrsDynamicExtFlow(*bx, pos_pd, sysbox);
			} else if (bx->type_neighborhood == 10) {  // left top
				pos_pd = pos + ax.box_axis1;         //--->right
				assignNeighobrsDynamicExtFlow(*bx, pos_pd, sysbox);
				pos_pd = pos - ax.box_axis2;         //--->bottom
				assignNeighobrsDynamicExtFlow(*bx, pos_pd, sysbox);
				pos_pd = pos - ax.box_diagonal_7_10; //--->right bot
				assignNeighobrsDynamicExtFlow(*bx, pos_pd, sysbox);
			} else if (bx->type_neighborhood == 11) {  // right top
				pos_pd = pos - ax.box_axis1;         //--->left
				assignNeighobrsDynamicExtFlow(*bx, pos_pd, sysbox);
				pos_pd = pos - ax.box_axis2;         //--->botom
				assignNeighobrsDynamicExtFlow(*bx, pos_pd, sysbox);
				pos_pd = pos - ax.box_diagonal_6_11; //--->left bot
				assignNeighobrsDynamicExtFlow(*bx, pos_pd, sysbox);
			}
		}
	}
}

void ExtensionalShearBoxSet::updateExtFlow(const struct Geometry::box3d &sysbox,
										   const struct BC::ExtFlowAxis &ax)
{
	double l1 = ax.box_axis1.norm();
	double l2 = ax.box_axis2.norm();
	double delta = 0.02;
	vec3d delta_n1 = delta*ax.box_axis1/l1;
	vec3d delta_n2 = delta*ax.box_axis2/l2;
	for (auto& bx : Boxes) {
		bx->resetTypes();
	}
	std::set <unsigned int> active_box;
	std::set <unsigned int> left_edge;
	std::set <unsigned int> bottom_edge;
	/* Top edge */
	vec3d pos = origin_ext_flow+ax.box_axis2;
	int jmax = (int)(l1/delta);
	for (int j=0; j<2*jmax; j++) {
		unsigned int box_label = whichBoxLabel(pos);
		if (box_labels[box_label]->type == 0) {
			if (j <= jmax ) {
				active_box.insert(box_label);
			}
			box_labels[box_label]->type +=8; // Top
		}
		pos += delta_n1;
		if (pos.x >= sysbox.lx) {
			break;
		}
	}
	/* Bottom edge */
	pos = origin_ext_flow+ax.box_axis1;
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
	pos = origin_ext_flow+ax.box_axis2;
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
	pos = origin_ext_flow+ax.box_axis1;
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
		if (pos.z >= sysbox.lz) {
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
	if (!twodimension) {
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
	updateNeighborsExtFlow(sysbox, ax);
	for (const auto& bx : Boxes) {
		bx->buildNeighborhoodContainer();
	}
}

void ExtensionalShearBoxSet::checkNeighobrType(std::set <unsigned int> &box_labels_)
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
						std::cerr << "not expected to be here" << std::endl;
						exit(1);
					}
				}
			}
		}
	}
}

vec3d ExtensionalShearBoxSet::periodicDiffShift(int i, int j, const struct BC::ExtFlowAxis &ax)
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
			pd_shift = -ax.box_axis1;
		}
	} else if (type_i == 3) { // right
		if (type_j== 2 || type_j== 6 || type_j == 10) { // including left
			pd_shift = ax.box_axis1;
		}
	} else if (type_i == 4) { // bottom
		if (type_j == 8 || type_j== 10 || type_j == 11) { // including top
			pd_shift = -ax.box_axis2;
			return pd_shift;
		}
	} else if (type_i == 8) { // top
		if (type_j == 4 || type_j == 6 || type_j == 7) { // including bot
			pd_shift = ax.box_axis2;
		}
	} else if (type_i == 6) { // left bottom
		if (type_j == 3 || type_j == 7) { // right and right bot
			pd_shift = -ax.box_axis1;
		} else if (type_j == 8 || type_j == 10) { // top or top left
			pd_shift = -ax.box_axis2;
		} else if (type_j == 11) { //right top
			pd_shift = -ax.box_diagonal_6_11;
		}
	} else if (type_i == 7) { // right bottom
		if (type_j == 2 || type_j == 6) { // left and left bot
			pd_shift = ax.box_axis1;
		} else if (type_j == 8 || type_j == 11) {//top and top right
			pd_shift = -ax.box_axis2;
		} else if (type_j == 10) { // left top
			pd_shift = -ax.box_diagonal_7_10;
		}
	} else if (type_i == 10) { // left top
		if (type_j == 3 || type_j == 11) { //right and right top
			pd_shift = -ax.box_axis1;
		} else if (type_j == 4 || type_j == 6) { //bottom and bottom left
			pd_shift = ax.box_axis2;
		} else if (type_j == 7) { // right bot
			pd_shift = ax.box_diagonal_7_10;
		}
	} else if (type_i == 11) { // right top
		if (type_j == 2 || type_j == 10) { //left and left top
			pd_shift = ax.box_axis1;
		} else if (type_j == 4 || type_j == 7) { //bot and bot right
			pd_shift = ax.box_axis2;
		} else if (type_j == 6) { // left bot
			pd_shift = ax.box_diagonal_6_11;
		}
	}
	return pd_shift;
}

void ExtensionalShearBoxSet::box(int i, vec3d position)
{
	//	Box* b = whichBox(position);
	Box* b = whichBox(position+origin_ext_flow);
	addToBox(i, b);
}

void ExtensionalShearBoxSet::yaplotBox(std::ofstream &fout_boxing)
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
			fout_boxing << "@ 0" << std::endl;
		} else if (box_color == 1) {
			fout_boxing << "@ 2" << std::endl;
		} else if (box_color == 2 || box_color == 3) {
			fout_boxing << "@ 3" << std::endl;
		} else if (box_color == 4 || box_color == 8) {
			fout_boxing << "@ 5" << std::endl;
		}  else if (box_color == 6) {
			fout_boxing << "@ 6" << std::endl;
		} else if (box_color == 7) {
			fout_boxing << "@ 6" << std::endl;
		} else if (box_color == 10) {
			fout_boxing << "@ 6" << std::endl;
		} else if (box_color == 11) {
			fout_boxing << "@ 6" << std::endl;
		} else {
			std::cerr << "unexpected bx type" << std::endl;
			std::cerr << box_color << std::endl;
			exit(1);
		}
		if (twodimension) {
			fout_boxing << "p 4 "<< left_bottom;
			fout_boxing << ' ' << right_bottom;
			fout_boxing << ' ' << right_top;
			fout_boxing << ' ' << left_top;
			fout_boxing << std::endl;
		} else {
			if (box_color != 0) {
				fout_boxing << "l "<< left_bottom << ' ' << right_bottom << std::endl;
				fout_boxing << "l " << right_bottom << ' '<< right_top << std::endl;
				fout_boxing << "l " << right_top << ' '<< left_top << std::endl;
				fout_boxing << "l " << left_top << ' '<< left_bottom << std::endl;
			}
		}
	}
	if (0) {
		fout_boxing << "r 1\n";
		for (const auto& bx : Boxes) {
			if (bx->type != 0 && bx->type_neighborhood >= 1) {
				if (bx->type_neighborhood == 1) {
					fout_boxing << "@ 2" << std::endl;
				} else if (bx->type_neighborhood == 2 || bx->type_neighborhood == 3) {
					fout_boxing << "@ 3" << std::endl;
				} else if (bx->type_neighborhood == 4 || bx->type_neighborhood == 8) {
					fout_boxing << "@ 4" << std::endl;
				} else {
					fout_boxing << "@ 0" << std::endl;
				}
				fout_boxing << "c " <<  bx->getPosition() << std::endl;
			}
		}
	}
}

}