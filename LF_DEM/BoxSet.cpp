#include "BoxSet.h"
#include "System.h"
using namespace std;

void
BoxSet::init(double interaction_dist, System *sys_){
	sys = sys_;
	boxMap = new Box* [sys->np()];
	for (int i=0; i<sys->np(); i++) {
		boxMap[i] = NULL;
	}
	double xratio = sys->lx()/interaction_dist;
	double yratio = sys->ly()/interaction_dist;
	double zratio = sys->lz()/interaction_dist;
	
	x_box_nb = (int)xratio;
	y_box_nb = (int)yratio;
	z_box_nb = (int)zratio;
	
	if (x_box_nb == 0) {
		x_box_nb = 1;
		box_xsize = sys->lx();
	}
	if (y_box_nb == 0) {
		y_box_nb = 1;
		box_ysize = sys->ly();
	}
	if (z_box_nb == 0) {
		z_box_nb = 1;
		box_zsize = sys->lz();
	}
	
	if (x_box_nb < 4 && y_box_nb < 4 && z_box_nb < 4) { // boxing useless: a neighborhood is the whole system
		cerr << "boxing useless: a neighborhood is the whole system" << endl;
		_is_boxed = false;
		box_xsize = sys->lx();
		box_ysize = sys->ly();
		box_zsize = sys->lz();
		box_xsize_half = 0.5*box_xsize;
		box_ysize_half = 0.5*box_ysize;
		box_zsize_half = 0.5*box_zsize;
		box_nb = 1;
		top_box_nb = 0;
		bottom_box_nb = 0;
		topbottom_box_nb = 1;
		bulk_box_nb = 0;
		Boxes = new Box* [box_nb];
		Boxes[0] = new Box;
		Boxes[0]->position.x = 0;
		Boxes[0]->position.y = 0;
		Boxes[0]->position.z = 0;
		Boxes[0]->is_bottom(true);
		Boxes[0]->is_top(true);
		TopBottomBoxes = new Box* [topbottom_box_nb];
		TopBottomBoxes[0] = Boxes[0];
		(Boxes[0])->neigh_nb(0);
	} else {
		_is_boxed = true;
		box_xsize = sys->lx()/x_box_nb;
		box_ysize = sys->ly()/y_box_nb;
		box_zsize = sys->lz()/z_box_nb;
		box_xsize_half = 0.5*box_xsize;
		box_ysize_half = 0.5*box_ysize;
		box_zsize_half = 0.5*box_zsize;
		allocateBoxes();
		// give them their position
		positionBoxes();
		// tell them their neighbors
		assignNeighbors();
	}
	cerr << box_nb << " " << interaction_dist << " " << box_xsize << endl;
}

void
BoxSet::allocateBoxes(){
	box_nb = x_box_nb*y_box_nb*z_box_nb;
	top_box_nb = x_box_nb*y_box_nb;
	bottom_box_nb = top_box_nb;
	bulk_box_nb = box_nb-top_box_nb-bottom_box_nb;
	topbottom_box_nb = 0;
	if (bulk_box_nb < 0) { // there is only one layer in the z direction ( ie BottomBoxes == TopBoxes )
		bulk_box_nb = 0;
		topbottom_box_nb = top_box_nb;
		top_box_nb = 0;
		bottom_box_nb = 0;
	}
	cerr << endl << "Boxer allocating :" << endl;
	cerr << box_nb << " boxes" << endl;
	cerr << top_box_nb << " top_boxes" << endl;
	cerr << bottom_box_nb << " bottom_boxes" << endl;
	cerr << bulk_box_nb << " bulk_boxes" << endl;
	cerr << topbottom_box_nb << " topbottom_boxes" << endl;
	Boxes = new Box* [box_nb];
	if (top_box_nb > 0) {
		TopBoxes = new Box* [top_box_nb];
	}
	if (bottom_box_nb > 0) {
		BottomBoxes = new Box* [bottom_box_nb];
	}
	if (topbottom_box_nb > 0) {
		TopBottomBoxes = new Box* [topbottom_box_nb];
	}
	if (bulk_box_nb > 0) {
		BulkBoxes = new Box* [bulk_box_nb];
	}
}

void
BoxSet::positionBoxes(){
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
	
	int neigh_nb = amax*bmax*cmax-1;
	double x,y,z;
	int label = 0;
	int toplabel = 0;
	int bottomlabel = 0;
	int topbottomlabel = 0;
	int bulklabel = 0;
	// position boxes
	for (int ix=0; ix<x_box_nb; ix++) {
		x = box_xsize*ix;
		for (int iy=0; iy<y_box_nb; iy++) {
			y = box_ysize*iy;
			for (int iz=0; iz<z_box_nb; iz++) {
				int extra_neigh=0;
				z = box_zsize*iz;
				Boxes[label] = new Box;
				Boxes[label]->position.x = x;
				Boxes[label]->position.y = y;
				Boxes[label]->position.z = z;
				if (iz == 0 && iz < z_box_nb-1) {// bottom box
					Boxes[label]->is_bottom(true);
					BottomBoxes[bottomlabel] = Boxes[label];
					if (x_box_nb > 3) {
						extra_neigh += bmax;// extra neighbors for bottom layer
					}
					bottomlabel++;
					(Boxes[label])->neigh_nb(neigh_nb+extra_neigh, amax*bmax+extra_neigh);
				}
				
				if (iz == z_box_nb-1 && iz > 0) {//top box
					Boxes[label]->is_top(true);
					TopBoxes[toplabel] = Boxes[label];
					if (x_box_nb > 3) {
						extra_neigh += bmax; // idem for top
					}
					toplabel++;
					(Boxes[label])->neigh_nb(neigh_nb+extra_neigh, amax*bmax+extra_neigh);
				}
				
				if (iz == 0 && iz == z_box_nb-1) {// bottom box
					Boxes[label]->is_bottom(true);
					Boxes[label]->is_top(true);
					TopBottomBoxes[topbottomlabel] = Boxes[label];
					if (x_box_nb > 3) {
						extra_neigh += 2*bmax;// top and bottom extra neighbors
					}
					topbottomlabel++;
					(Boxes[label])->neigh_nb(neigh_nb+extra_neigh, amax*bmax+extra_neigh);
				}
				if (iz > 0 && iz < z_box_nb-1) {// bulk box
					BulkBoxes[bulklabel] = Boxes[label];
					bulklabel++;
					(Boxes[label])->neigh_nb(neigh_nb);
				}
				label++;
			}
		}
	}
}

void
BoxSet::assignNeighborsBulk(){
	for (int i=0; i<bulk_box_nb; i++) {
		vec3d pos = BulkBoxes[i]->position;
		vec3d delta;
		pos.x += box_xsize_half;
		pos.y += box_ysize_half;
		pos.z += box_zsize_half;
		int label=0;
		for (int a=0; a<amax; a++) {
			delta.x = (a-1)*box_xsize;
			for (int b=0; b<bmax; b++) {
				delta.y = (b-1)*box_ysize;
				for (int c=0; c<cmax; c++) {
					delta.z = (c-1)*box_zsize;
					bool successful_add = (BulkBoxes[i])->neighbor(label, WhichBox(pos+delta));
					if (successful_add) {
						label++;
					}
				}
			}
		}
	}
}

void
BoxSet::assignNeighborsBottom(){
	for (int i=0; i<bottom_box_nb; i++) {
		vec3d pos = BottomBoxes[i]->position;
		vec3d delta;
		pos.x += box_xsize_half;
		pos.y += box_ysize_half;
		pos.z += box_zsize_half;
		int label=0;
		// boxes  at same level and above first: these are fixed once and for all in the simulation
		for (int a=0; a<amax; a++) {
			delta.x = (a-1)*box_xsize;
			for (int b=0; b<bmax; b++) {
				delta.y = (b-1)*box_ysize; // same level and above
				for (int c=0; c<2; c++) {
					delta.z = c*box_zsize; // same level and above
					bool successful_add = (BottomBoxes[i])->neighbor(label, WhichBox(pos+delta));
					if (successful_add) {
						label++;
					}
				}
			}
		}
		// identities of boxes below are added at the end, and they will be updated at each time step
		pos.x -= 0.499999999*box_xsize;
		delta.z = -box_zsize;  // below
		for (int a=0; a<amax+1; a++) {  // one more line of boxes
			delta.x = (a-1)*box_xsize;
			for (int b=0; b<bmax; b++) {
				delta.y = (b-1)*box_ysize;
				(BottomBoxes[i])->probing_positions(label, pos+delta);
				bool successful_add = (BottomBoxes[i])->neighbor(label, WhichBox(pos+delta));
				if (successful_add) {
					label++;
				}
			}
		}
	}
}

void
BoxSet::assignNeighborsTop(){
	for (int i=0; i<top_box_nb; i++) {
		vec3d pos = TopBoxes[i]->position;
		vec3d delta;
		pos.x += box_xsize_half;
		pos.y += box_ysize_half;
		pos.z += box_zsize_half;
		int label = 0;
		// boxes at same level and below first: these are fixed once and for all in the simulation
		for (int a=0; a<amax; a++) {
			delta.x = (a-1)*box_xsize;
			for (int b=0; b<bmax; b++) {
				delta.y = (b-1)*box_ysize;  // below and same level
				for (int c=-1; c<1; c++) {
					delta.z = c*box_zsize; // below and same level
					bool successful_add = (TopBoxes[i])->neighbor(label, WhichBox(pos+delta));
					if (successful_add) {
						label++;
					}
				}
			}
		}
		// identities of boxes above are added at the end, and they will be updated at each time step
		pos.x -= 0.499999999*box_xsize;
		delta.z = +1.*box_zsize;  // above
		for (int a=0; a<amax+1; a++) {  // one more line of boxes
			delta.x = (a-1)*box_xsize;
			for (int b=0; b<bmax; b++) {
				delta.y = (b-1)*box_ysize;
				(TopBoxes[i])->probing_positions(label, pos+delta);
				bool successful_add = (TopBoxes[i])->neighbor(label, WhichBox(pos+delta));
				if (successful_add) {
					label++;
				}
			}
		}
	}
}

void
BoxSet::assignNeighborsTopBottom(){
	for (int i=0; i<topbottom_box_nb; i++) {
		vec3d pos = TopBottomBoxes[i]->position;
		vec3d delta;
		pos.x += box_xsize_half;
		pos.y += box_ysize_half;
		pos.z += box_zsize_half;
		int label = 0;
		// boxes at same level first: these are fixed once and for all in the simulation
		delta.z = 0;  // same level
		for (int a=0; a<amax; a++) {
			delta.x = (a-1)*box_xsize;
			for (int b=0; b<bmax; b++) {
				delta.y = (b-1)*box_ysize;
				bool successful_add = (TopBottomBoxes[i])->neighbor(label, WhichBox(pos+delta));
				if (successful_add) {
					label++;
				}
			}
		}
		// identities of boxes above and below are added, and they will be updated at each time step
		pos.x -= 0.499999999*box_xsize;
		delta.z = box_zsize;  // above
		for (int a=0; a<amax+1; a++) {  // one more line of boxes
			delta.x = (a-1)*box_xsize;
			for (int b=0; b<bmax; b++) {
				delta.y = (b-1)*box_ysize;
				(TopBottomBoxes[i])->probing_positions(label, pos+delta);
				bool successful_add = (TopBottomBoxes[i])->neighbor(label, WhichBox(pos+delta));
				if (successful_add) {
					label++;
				}
			}
		}
		delta.z = -box_zsize;  // below
		for (int a=0; a<amax+1; a++) {  // one more line of boxes
			delta.x = (a-1)*box_xsize;
			for (int b=0; b<bmax; b++) {
				delta.y = (b-1)*box_ysize;
				(TopBottomBoxes[i])->probing_positions(label, pos+delta);
				bool successful_add = (TopBottomBoxes[i])->neighbor(label, WhichBox(pos+delta));
				if (successful_add) {
					label++;
				}
			}
		}
	}
}

void
BoxSet::assignNeighbors(){
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
	for (int i=0; i<box_nb; i++) {
		delete Boxes[i];
	}
	DELETE(Boxes);
	if (bulk_box_nb > 0) {
		DELETE(BulkBoxes);
	}
	if (_is_boxed) {
		DELETE(TopBoxes);
		DELETE(BottomBoxes);
	} else {
		DELETE(TopBottomBoxes);
	}
	DELETE(boxMap);
}

void
BoxSet::updateNeighbors(Box* b){
	b->reset_moving_neighbors();
	vec3d *probes = b->probing_positions();
	int probes_nb = b->probe_nb();
	int moving_label = 0;
	for (int i=0; i < probes_nb; i++) {
		//	  cout <<  i << " " << b->position.x << " " << b->position.y<< " " << b->position.z << " "<< (WhichBox(probes[i]))->position.x << " " << (WhichBox(probes[i]))->position.y<< " " << (WhichBox(probes[i]))->position.z << endl;
		bool successful_add = b->moving_neighbor(moving_label, WhichBox(probes[i]));
		if (successful_add) {
			moving_label++;
		}
	}
}

/*****
 UpdateNeighbors()
 
 At each time step, we need to check if neighborhood on top and bottom boxes have changed.
 Bulk boxes do not need to be updated.
 *****/
void
BoxSet::updateNeighbors(){
	// top boxes
	for(int i=0; i<top_box_nb; i++) {
		updateNeighbors(TopBoxes[i]);
	}
	// bottom boxes
	for(int i=0; i<bottom_box_nb; i++) {
		updateNeighbors(BottomBoxes[i]);
	}
	// topbottom boxes
	for(int i=0; i<topbottom_box_nb; i++) {
		updateNeighbors(TopBottomBoxes[i]);
	}
}

//public methods
void
BoxSet::update(){
	if (is_boxed()) {
		updateNeighbors();
	}
	for(int i=0; i<box_nb; i++) {
		Boxes[i]->build_neighborhood_container();
	}
}

bool
BoxSet::is_boxed(){
	return _is_boxed;
}

Box*
BoxSet::WhichBox(vec3d pos){
	return WhichBox(&pos);
}

Box*
BoxSet::WhichBox(vec3d *pos){
	sys->periodize(*pos);
	int ix = (int)(pos->x/box_xsize);
	int iy;
	if (sys->dimension == 2) {
		iy = 0;
	} else {
		iy = (int)(pos->y/box_ysize);
	}
	int iz = (int)(pos->z/box_zsize);
	int label = ix*y_box_nb*z_box_nb+iy*z_box_nb+iz;
	// cout << " Which Box : " << label << endl;
	// cout << x_box_nb << " " << y_box_nb << " " << z_box_nb<< " "<< ix  <<" "<< iy << " " << iz << endl;
	// cout <<pos->x << " " << pos->y << " " << pos->z << endl;
	// cout << (Boxes[label])->position.x<< " " << (Boxes[label])->position.y<< " " << (Boxes[label])->position.z<< endl;
	return Boxes[label];
}

void
BoxSet::box(int i){
	Box *b = WhichBox(sys->position[i]);
	if (b != boxMap[i]) {
		b->add(i);
		if (boxMap[i] != NULL) {
			boxMap[i]->remove(i);
		}
		boxMap[i] = b;
	}
}

vector<int>::iterator
BoxSet::neighborhood_begin(int i){
	return (boxMap[i])->neighborhood_begin();
}

vector<int>::iterator
BoxSet::neighborhood_end(int i){
	return (boxMap[i])->neighborhood_end();
}

void
BoxSet::printBoxNetwork(){
	for (int i=0; i<box_nb; i++) {
		Box** neighbors = (Boxes[i])->neighbors();
		int nneigh = (Boxes[i])->neigh_nb();
		for (int j = 0; j< nneigh;j++) {
			//		  if(!(Boxes[i])->is_top() && !(Boxes[i])->is_bottom()){
			if ((Boxes[i])->is_top()) {
				cerr << i << " " <<box_xsize << " " <<nneigh << " " << (Boxes[i])->position.x << " ";
				cerr << (Boxes[i])->position.y <<  " " << (neighbors[j])->position.x << " ";
				cerr << (neighbors[j])->position.y << " " << (Boxes[i])->is_top() << " ";
				cerr << (Boxes[i])->is_bottom() << endl;
			}
		}
	}
}
