#ifndef __LF_DEM__MatrixBlocks__
#define __LF_DEM__MatrixBlocks__
#include <array>
#include "vec3d.h"

struct ODBlock {
  /**
  An off-diagonal block of the resistance matrix R_FU computing the lubrication forces/torques
  from the non-affine velocities/angular velocities: F = R_FU.U
  General structure
  =================
  The off-diagonal blocks of R_FU are 6by6 matrices with terms coming
  from A, B, \tilde B and C matrices as defined in Jeffrey and Onishi 1984.
  They are made of 3by3 subblocks which, loosely speaking, have different origins as:
  | A   \tildeB |
  |	B      C    |
  Symmetries
  ==========
  The "A" and "C" subblocks are symmetric,
	and "B" and "\tildeB" sublocks are antisymmetric (see Jeffrey and Onishi 1984).
  As a consequence, each ODBlock has 24 independent elements.
  The stored values are organized in columns (internally as 6 std::array's, one for each column in the block),
  and correspond to the following elements:
     | col0[0]  .        .        0        .        .       |
  "A"| col0[1]  col1[0]  .        col3[0]  0        .       |
     | col0[2]  col1[1]  col2[0]  col3[1]  col4[0]  0       |
     | 0        .        .        col3[2]  .        .       |
  "B"| col0[3]  0        .        col3[3]  col4[1]  .       |
     | col0[4]  col1[2]  0        col3[4]  col4[2]  col5[0] |
                                           "C"

  The full block can be reconstructed as:
     | col0[0]  col0[1]  col0[2]  0       -col3[0] -col3[1] |
  "A"| col0[1]  col1[0]  col1[1]  col3[0]  0       -col4[0] |
     | col0[2]  col1[1]  col2[0]  col3[1]  col4[0]  0       |
     | 0       -col0[3] -col0[4]  col3[2]  col3[3]  col3[4] |
  "B"| col0[3]  0       -col1[2]  col3[3]  col4[1]  col4[2] |
     | col0[4]  col1[2]  0        col3[4]  col4[2]  col5[0] |
                                           "C"
	*/
	std::array<double, 5> col0;
	std::array<double, 3> col1;
	std::array<double, 1> col2;
	std::array<double, 5> col3;
	std::array<double, 3> col4;
	std::array<double, 1> col5;
};

inline void resetODBlock(struct ODBlock &b)
{
	b.col0.fill(0);
	b.col1.fill(0);
	b.col2.fill(0);
	b.col3.fill(0);
	b.col4.fill(0);
	b.col5.fill(0);
}

struct DBlock {
	/**
	A diagonal block of the resistance matrix R_FU computing the lubrication forces/torques
	from the non-affine velocities/angular velocities: F = R_FU.U

	General structure
	=================
	The diagonal blocks of R_FU are 6by6 matrices with terms coming
	from A, B, \tilde B and C matrices as defined in Jeffrey and Onishi 1984.
	They are made of 3by3 subblocks which, loosely speaking, have different origins as:
	| A   \tildeB |
	|	B      C    |

	Symmetries
	==========
	The blocks are symmetric, and their "B" and "\tildeB" sublocks are antisymmetric (see Jeffrey and Onishi 1984).
 	As a consequence, each DBlock has 18 independent elements (6 diagonal terms and 12 off diagonal).

	The stored values are organized in columns (internally as 6 std::array's, one for each column in the block),
	and correspond to the following elements:
    | col0[0]  .        .        .        .        .       |
 "A"| col0[1]  col1[0]  .        .        .        .       |
    | col0[2]  col1[1]  col2[0]  .        .        .       |
    | 0        .        .        col3[0]  .        .       |
 "B"| col0[3]  0        .        col3[1]  col4[0]  .       |
    | col0[4]  col1[2]  0        col3[2]  col4[1]  col5[0] |
                                          "C"

   The full block can be reconstructed as:
    | col0[0]  col0[1]  col0[2]  0       -col0[3] -col0[4] |
 "A"| col0[1]  col1[0]  col1[1]  col0[3]  0       -col1[2] |
    | col0[2]  col1[1]  col2[0]  col0[4]  col1[2]  0       |
    | 0       -col0[3] -col0[4]  col3[0]  col3[1]  col3[2] |
 "B"| col0[3]  0       -col1[2]  col3[1]  col4[0]  col4[1] |
    | col0[4]  col1[2]  0        col3[2]  col4[1]  col5[0] |
                                          "C"

*/

	std::array<double, 5> col0;
	std::array<double, 3> col1;
	std::array<double, 1> col2;
	std::array<double, 3> col3;
	std::array<double, 2> col4;
	std::array<double, 1> col5;

	inline DBlock& operator += (const DBlock &b)
	{
		for (unsigned int i=0; i<col0.size(); i++) {
			col0[i] += b.col0[i];
		}
		for (unsigned int i=0; i<col1.size(); i++) {
			col1[i] += b.col1[i];
		}
		col2[0] += b.col2[0];
		for (unsigned int i=0; i<col3.size(); i++) {
			col3[i] += b.col3[i];
		}
		for (unsigned int i=0; i<col4.size(); i++) {
			col4[i] += b.col4[i];
		}
		col5[0] += b.col5[0];
		return *this;
	}
};

inline void resetDBlock(struct DBlock& b)
{
	b.col0.fill(0);
	b.col1.fill(0);
	b.col2.fill(0);
	b.col3.fill(0);
	b.col4.fill(0);
	b.col5.fill(0);
}

namespace ResistanceBlocks {

enum class Label {
	ForceVel,       // top left, symmetric
	ForceAngVel,	// top right, antisymmetric
	TorqueVel,		// bottom left, antisymmetric
	TorqueAngVel	// bottom right, symmetric
};

class DBlockBuilder
{
public:
	DBlockBuilder();
	void addIdentity(Label loc, double prefactor);
	void addDyadic(Label loc, vec3d v, double prefactor);
	void addComplementaryDyadic(Label loc, vec3d v, double prefactor);
	void addVectorProduct(Label loc, vec3d v);

	DBlock block;
private:

};

class ODBlockBuilder
{
public:
	ODBlockBuilder();
	void addIdentity(Label loc, double prefactor);
	void addDyadic(Label loc, vec3d v, double prefactor);
	void addComplementaryDyadic(Label loc, vec3d v, double prefactor);
	void addVectorProduct(Label loc, vec3d v);

	ODBlock block;
private:

};
}

#endif
