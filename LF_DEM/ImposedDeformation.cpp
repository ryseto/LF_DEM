#include <tuple>
#include "ImposedDeformation.h"

namespace Geometry {

void ImposedDeformation::setShape(Sym2Tensor sym, vec3d antisym)
{
	sym_shape = sym;
	antisym_shape = antisym;
	shape = sym2FullMatrix(sym) + vec2AntiSymFullMatrix(antisym_shape);

	grad_u = shear_rate*shape;
	sym_grad_u = shear_rate*sym_shape;
	antisym_grad_u = shear_rate*antisym_shape;
	printShape();
}

void ImposedDeformation::setShape(const matrix &shape_)
{
	shape = shape_;
	std::tie(sym_shape, antisym_shape) = symmetryDecomposition(shape);

	grad_u = shear_rate*shape;
	sym_grad_u = shear_rate*sym_shape;
	antisym_grad_u = shear_rate*antisym_shape;
	printShape();
}

void ImposedDeformation::setRate(double rate, bool verbose)
{
	shear_rate = rate;
	if (verbose) {
		std::cout << "\tImposedDeformation:: Setting deformation rate "<< getRate() << std::endl;
	}
	grad_u = shear_rate*shape;
	sym_grad_u = shear_rate*sym_shape;
	antisym_grad_u = shear_rate*antisym_shape;
}

void ImposedDeformation::printShape()
{
	std::cout << "\tImposedDeformation:: Flow shape: " << std::endl;
	for (unsigned i=0; i<3; i++) {
		std::cout << "\t\t\t";
		for (unsigned j=0; j<3; j++) {
			std::cout << " " << shape.get(i, j);
		}
		std::cout << std::endl;
	}
}
}
