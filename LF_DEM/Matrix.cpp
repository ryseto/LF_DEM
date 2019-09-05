#include "Matrix.h"
#include "Sym2Tensor.h"

std::pair<Sym2Tensor, vec3d> symmetryDecomposition(const matrix &m)
{
	Sym2Tensor sympart;
	sympart.setSymmetrize(m);

	matrix antisympart = 0.5*(m + (-1.)*m.transpose());
	vec3d omega = {-antisympart.get(1, 2), antisympart.get(0, 2), -antisympart.get(0, 1)};

	return std::make_pair(sympart, omega);
}

matrix vec2AntiSymFullMatrix(vec3d omega) 
{
	return matrix( 0,      -omega.z, omega.y,
				   omega.z, 0,      -omega.x,
				  -omega.y, omega.x, 0       );
}

matrix sym2FullMatrix(Sym2Tensor sym)
{
	return matrix(sym.elm[0], sym.elm[1], sym.elm[2],
				  sym.elm[1], sym.elm[4], sym.elm[3],
				  sym.elm[2], sym.elm[3], sym.elm[5]);
}

