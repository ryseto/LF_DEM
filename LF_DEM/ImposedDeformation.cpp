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
}

void ImposedDeformation::setShape(const matrix &shape_)
{
	shape = shape_;
	std::tie(sym_shape, antisym_shape) = symmetryDecomposition(shape);

	grad_u = shear_rate*shape;
	sym_grad_u = shear_rate*sym_shape;
	antisym_grad_u = shear_rate*antisym_shape;
}

void ImposedDeformation::setRate(double rate)
{
	shear_rate = rate;
	std::cout << "\tImposedDeformation : Setting deformation rate "<< getRate() << std::endl;
	grad_u = shear_rate*shape;
	sym_grad_u = shear_rate*sym_shape;
	antisym_grad_u = shear_rate*antisym_shape;
}

}
