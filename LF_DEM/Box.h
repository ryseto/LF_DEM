#ifndef __LF_DEM__Box__
#define __LF_DEM__Box__

#include "vec3d.h"

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

 protected:
  
 public:

  vec3d position;
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
  void probing_positions(int label, vec3d pos);

  void is_top(bool);
  void is_bottom(bool);
  bool is_top();
  bool is_bottom();

  Box();
  ~Box();

};

#endif /* defined(__LF_DEM__Box__) */
