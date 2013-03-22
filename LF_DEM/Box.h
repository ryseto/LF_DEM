#ifndef __LF_DEM__Box__
#define __LF_DEM__Box__
#define __vector_container__
#ifndef __vector_container__
#define __flist_container__
#endif
#include "vec3d.h"
#include <set>
#include <vector>
using namespace std;

class Box{
private:
	Box **_neighbors;
	Box **_moving_neighbors;
	int _neigh_nb;
	int _moving_neigh_nb;
	int _still_neigh_nb;
	bool _is_bottom;
	bool _is_top;
	vec3d *_probing_positions;// if box is at top/bottom, these positions allow to figure out which boxes are above/below you
	int _probe_nb;
	bool can_be_added(int, Box*);
	
	/* #ifdef __vector_container__ */
	/*   vector <int> *container; */
	/* #endif */
	/* #ifdef __flist_container__ */
	/*   forward_list <int> *container; */
	/* #endif */
	/*   map <int,int> key; */
	vector <int> neighborhood_container;
protected:
	
public:
	Box();
	~Box();
	vec3d position;
	set <int> *container;
	void neigh_nb(int n, int moving_n=0);
	bool neighbor(int label, Box* neigh_box);
	bool moving_neighbor(int moving_label, Box* neigh_box);
	void reset_moving_neighbors();
	Box** neighbors(){
		return _neighbors;
	}
	int neigh_nb(){
		return _neigh_nb;
	}
	int probe_nb(){
		return _probe_nb;
	}
	vec3d* probing_positions(){
		return _probing_positions;
	}
	void probing_positions(int label, const vec3d &pos);
	void is_top(bool);
	void is_bottom(bool);
	bool is_top();
	bool is_bottom();
	void add(int);
	void remove(int);
	set<int>::iterator begin(){
		return container->begin();
	}
	set<int>::iterator end(){
		return container->end();
	}
	vector<int>::iterator neighborhood_begin(){
		return neighborhood_container.begin();
	}
	vector<int>::iterator neighborhood_end(){
		return neighborhood_container.end();
	}
	size_t container_size(){
		return container->size();
	}
	void build_neighborhood_container();
};

#endif /* defined(__LF_DEM__Box__) */
