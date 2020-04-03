#ifndef LF_DEM_ImposedDeformation_h
#define LF_DEM_ImposedDeformation_h

#include "vec3d.h"
#include "Sym2Tensor.h"
#include "Matrix.h"

namespace Geometry {

class ImposedDeformation {
public:
	ImposedDeformation() : shear_rate(0) {};

	void setShape(Sym2Tensor sym, vec3d antisym);
	void setShape(const matrix &shape_);
	void setRate(double rate, bool verbose=false);
	bool zero_shear() {return shear_rate == 0;};

	const matrix& getShape() {return shape;};
	const Sym2Tensor& getSymShape()  {return sym_shape;};
	const vec3d& getAntiSymShape()  {return antisym_shape;};

	const matrix& getGradU() {return grad_u;};
	const Sym2Tensor& getSymGradU()  {return sym_grad_u;};
	const vec3d& getAntiSymGradU()  {return antisym_grad_u;};

	double getRate() {return shear_rate;};

private:
	matrix shape;
	Sym2Tensor sym_shape; // E/shear_rate: "shape" of the flow
	vec3d antisym_shape;  // omega/shear_rate: "shape" of the flow

	matrix grad_u;
	Sym2Tensor sym_grad_u;
	vec3d antisym_grad_u;

	double shear_rate;

	void printShape();
};

}
#endif