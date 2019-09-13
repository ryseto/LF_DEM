#include <stdexcept>
#include "MatrixBlocks.h"

namespace ResistanceBlocks {

DBlockBuilder::DBlockBuilder()
{
	resetDBlock(block);
}


void DBlockBuilder::addIdentity(Label loc, double prefactor)
{
	if (loc == Label::ForceVel) {
		block.col0[0] += prefactor;
		block.col1[0] += prefactor;
		block.col2[0] += prefactor;
	} else if (loc == Label::TorqueAngVel) {
		block.col3[0] += prefactor;
		block.col4[0] += prefactor;
		block.col5[0] += prefactor;
	} else {
		throw std::runtime_error(" DBlockBuilder:: Impossible to add identity to an asymmetric subblock.");
	}
}

void DBlockBuilder::addDyadic(Label loc, vec3d v, double prefactor)
{
	if (loc == Label::ForceVel) {
		block.col0[0] += prefactor*v.x*v.x;
		block.col0[1] += prefactor*v.x*v.y;
		block.col0[2] += prefactor*v.x*v.z;
		block.col1[0] += prefactor*v.y*v.y;
		block.col1[1] += prefactor*v.y*v.z;
		block.col2[0] += prefactor*v.z*v.z;
	} else if (loc == Label::TorqueAngVel) {
		block.col3[0] += prefactor*v.x*v.x;
		block.col3[1] += prefactor*v.x*v.y;
		block.col3[2] += prefactor*v.x*v.z;
		block.col4[0] += prefactor*v.y*v.y;
		block.col4[1] += prefactor*v.y*v.z;
		block.col5[0] += prefactor*v.z*v.z;
	} else {
		throw std::runtime_error(" DBlockBuilder:: Impossible to add projection to an asymmetric subblock.");
	}
}

void DBlockBuilder::addComplementaryDyadic(Label loc, vec3d v, double prefactor)
{
	// Unnormalized projection on v
	// 
	// Writing u = v/|v|
	// returns prefactor*|v|^2*(1-uu)
	addIdentity(loc, prefactor*v.sq_norm());
	addDyadic(loc, v, -prefactor); 
}

void DBlockBuilder::addVectorProduct(Label loc, vec3d v)
{
	if (loc == Label::ForceAngVel) {
		throw std::runtime_error(" DBlockBuilder:: Cannot specify ForceAngVel subblock, please specify TorqueVel subblock.");
	}
	if (loc == Label::TorqueVel) {
		block.col0[3] += v.z;
		block.col0[4] += -v.y;
		block.col1[2] += v.x;
	} else {
		throw std::runtime_error(" DBlockBuilder:: Impossible to add vector product to a symmetric subblock.");
	}
}

ODBlockBuilder::ODBlockBuilder()
{
	resetODBlock(block);
}


void ODBlockBuilder::addIdentity(Label loc, double prefactor)
{
	if (loc == Label::ForceVel) {
		block.col0[0] += prefactor;
		block.col1[0] += prefactor;
		block.col2[0] += prefactor;
	} else if (loc == Label::TorqueAngVel) {
		block.col3[2] += prefactor;
		block.col4[1] += prefactor;
		block.col5[0] += prefactor;
	} else {
		throw std::runtime_error(" ODBlockBuilder:: Impossible to add identity to an asymmetric subblock.");
	}
}

void ODBlockBuilder::addDyadic(Label loc, vec3d v, double prefactor)
{
	if (loc == Label::ForceVel) {
		block.col0[0] += prefactor*v.x*v.x;
		block.col0[1] += prefactor*v.x*v.y;
		block.col0[2] += prefactor*v.x*v.z;
		block.col1[0] += prefactor*v.y*v.y;
		block.col1[1] += prefactor*v.y*v.z;
		block.col2[0] += prefactor*v.z*v.z;
	} else if (loc == Label::TorqueAngVel) {
		block.col3[2] += prefactor*v.x*v.x;
		block.col3[3] += prefactor*v.x*v.y;
		block.col3[4] += prefactor*v.x*v.z;
		block.col4[1] += prefactor*v.y*v.y;
		block.col4[2] += prefactor*v.y*v.z;
		block.col5[0] += prefactor*v.z*v.z;
	} else {
		throw std::runtime_error(" ODBlockBuilder:: Impossible to add projection to an asymmetric subblock.");
	}
}

void ODBlockBuilder::addComplementaryDyadic(Label loc, vec3d v, double prefactor)
{
	// prefactor*( (v.v)Id - vv )
	addIdentity(loc, prefactor*v.sq_norm());
	addDyadic(loc, v, -prefactor); 
}

void ODBlockBuilder::addVectorProduct(Label loc, vec3d v)
{
	if (loc == Label::TorqueVel) {
		block.col0[3] += v.z;
		block.col0[4] += -v.y;
		block.col1[2] += v.x;
	} else if (loc == Label::ForceAngVel) {
		block.col3[0] += v.z;
		block.col3[1] += -v.y;
		block.col4[0] += v.x;
	} else {
		throw std::runtime_error(" ODBlockBuilder:: Impossible to add vector product to a symmetric subblock.");
	}
}

} // namespace ResistanceBlocks