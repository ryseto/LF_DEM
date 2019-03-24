//
//  SolventFlow.cpp
//  LF_DEM
//
//  Created by Ryohei Seto on 2019/03/16.
//  Copyright © 2019 Ryohei Seto. All rights reserved.
//

#include "SolventFlow.h"
#include "System.h"
#include <vector>

SolventFlow::SolventFlow()
{
	psolver = new Eigen::SimplicialLDLT <SpMat>;
}

SolventFlow::~SolventFlow()
{
	delete [] psolver;
}

void SolventFlow::init(System* sys_)
{
	sys = sys_;
	nx = 20;
	nz = 20;
	n = nx*nz;
	dx = sys->get_lx()/nx;
	dz = sys->get_lz()/nz; // @@ need to be changed.
	pressure.resize(n);
	u_x.resize(n);
	u_z.resize(n);
	u_sol_x.resize(n);
	u_sol_z.resize(n);
	u_sol_ast_x.resize(n);
	u_sol_ast_z.resize(n);
	u_particle_x.resize(n);
	u_particle_z.resize(n);
	div_u_sol_ast.resize(n);
	phi.resize(n);
	pos.resize(n);
	for (int iz=0; iz<nz; iz++) {
		for (int ix=0; ix<nx; ix++) {
			pos[ix+nx*iz].set(dx*ix+dx/2, 0, dz*iz+dz/2);
		}
	}
	for (int k = 0; k < n; k++) {
		u_sol_x[k] = 0;
		u_sol_z[k] = 0;
	}
}

void SolventFlow::initPoissonSolver()
{
	std::vector< std::vector<double> > lmat;
	lmat.resize(n);
	for (int l=0; l<n; l++){
		lmat[l].resize(n);
	}
	for (int l=0; l<n; l++){
		for (int m=0; m<n; m++){
			lmat[l][m] = 0;
		}
	}
	
	double exz = 2*(1/(dx*dx)+1/(dz*dz));
	double ex = 1/(dx*dx);
	double ez = 1/(dz*dz);
	for (int zi = 0; zi < nz; zi++){
		for (int xi = 0; xi < nx; xi++){
			// b(xi, zi)
			int k = q(xi, zi);
			lmat[k][k] = -exz;  //p(xi,zi) --> b(xi,zi)
			lmat[k][q(xi+1, zi  )] =   ex;  //p(xi+1,zi) --> b(xi,zi)
			lmat[k][q(xi-1, zi  )] =   ex;  //p(xi-1,zi) --> b(xi,zi)
			lmat[k][q(xi  , zi+1)] =   ez;  //p(xi,zi+1) --> b(xi,zi)
			lmat[k][q(xi  , zi-1)] =   ez;  //p(xi,zi-1) --> b(xi,zi)
		}
	}
	std::vector<T> t_lmat;            // list of non-zeros coefficients
	for (int l=0; l<n; l++) {
		for (int k=0; k<n; k++) {
			if (lmat[l][k] != 0) {
				// std::cerr << k << ' ' << l << ' ' << lmat[l][k] << std::endl;
				t_lmat.push_back(T(k, l, lmat[l][k]));
			}
		}
	}
	std::cerr << t_lmat.size() << std::endl;
	lap_mat.resize(n, n);
	lap_mat.setFromTriplets(t_lmat.begin(), t_lmat.end());
	b.resize(n);
	x.resize(n);
	psolver->analyzePattern(lap_mat);   // for this step the numerical values of A are not used
	psolver->factorize(lap_mat);
	if (psolver->info() != Eigen::Success) {
		std::cerr << "decomposition failed" << std::endl;
		return;
	}
}

void SolventFlow::calcMeshVelocity()
{
	for (int k=0; k<n; k++) {
		u_x[k] = 0;
		u_z[k] = 0;
		phi[k] = 0;
	}
	double x, z, radius;
	double ux, uz;
	double ave_range = 1.5;
	double sq_ave_range = ave_range*ave_range;
	std::vector <int> mesh_id;
	std::vector <double> ux_values;
	std::vector <double> uz_values;
	std::vector <double> phi_values;
	for (int i=0; i < sys->get_np(); i++) {
		x = sys->position[i].x;
		z = sys->position[i].z;
		ux = sys->na_velocity[i].x;
		uz = sys->na_velocity[i].z;
		radius = sys->radius[i];
		int ix = x /dx;
		int iz = z /dz;
		mesh_id.clear();
		ux_values.clear();
		uz_values.clear();
		phi_values.clear();
		double total_weight_ux = 0;
		double total_weight_uz = 0;
		double total_weight_phi = 0;
		double particle_volume = M_PI*sys->radius[i]*sys->radius[i];
		for (int l=-1; l<=2; l++){
			for (int m=-1; m<=2; m++){
				int iix = ix+l;
				int iiz = iz+m;
				double xo = x;
				double zo = z;
				if (iix < 0) {
					iix += nx;
					xo += sys->get_lx();
				} else if (iix >= nx) {
					iix -= nx;
					xo -= sys->get_lx();
				}
				if (iiz < 0) {
					iiz += nz;
					zo += sys->get_lz();
				} else if (iiz >= nz) {
					iiz -= nz;
					zo -= sys->get_lz();
				}
				int k = iix + iiz*nx;
				double x_p = pos[k].x;
				double z_p = pos[k].z;
				double x_u = x_p - dx/2;
				double z_u = z_p - dz/2;
				double xx = (xo-x_p)*(xo-x_p);
				double zz = (zo-z_p)*(zo-z_p);
				double dist_sq_phi = xx+zz;
				double dist_sq_ux =  (xo-x_u)*(xo-x_u) + zz;
				double dist_sq_uz =  xx + (zo-z_u)*(zo-z_u);
				double w_ux = exp(-dist_sq_ux/sq_ave_range);
				double w_uz = exp(-dist_sq_uz/sq_ave_range);
				double w_phi = exp(-dist_sq_phi/sq_ave_range);
				total_weight_ux += w_ux;
				total_weight_uz += w_uz;
				total_weight_phi += w_phi;
				mesh_id.push_back(k);
				double phi = particle_volume/(dx*dz);
				ux_values.push_back(ux*radius*w_ux);
				uz_values.push_back(uz*radius*w_uz);
				phi_values.push_back(phi*w_phi);
			}
		}
		for (int l=0; l<mesh_id.size(); l++) {
			int k = mesh_id[l];
			u_x[k] += ux_values[l]/total_weight_ux;
			u_z[k] += uz_values[l]/total_weight_uz;
			phi[k] += phi_values[l]/total_weight_phi;
		}
	}
}

void SolventFlow::update(double pressure_difference_)
{
	static std::ofstream fout_tmp("debug.dat");
	pressure_difference = pressure_difference_;
	calcMeshVelocity();
	double rho = 0.1;
	d_tau = sys->dt/rho;
	for (int k=0; k<n; k++) {
		double coeff = d_tau*phi[k];
		u_sol_ast_x[k] = u_sol_x[k] + coeff*u_x[k];
		u_sol_ast_z[k] = u_sol_z[k] + coeff*u_z[k];
	}
	calcVelocityDivergence();
	solvePressure();
	updateSolventFlow();
	fout_tmp << std::endl;
}

void SolventFlow::calcVelocityDivergence()
{
	double dux_dx, duz_dz;
	for (int j=0; j<nz; j++){
		for (int i=0; i<nx; i++){
			int ip1 = i+1;
			if (ip1 == nx) {
				ip1 = 0;
			}
			int jp1 = j+1;
			if (jp1 == nz) {
				jp1 = 0;
			}
			int k = i + j*nx;
			dux_dx = (u_sol_ast_x[ip1 + j*nx] - u_sol_ast_x[k])/dx;
			duz_dz = (u_sol_ast_z[i + jp1*nx] - u_sol_ast_z[k])/dz;
			div_u_sol_ast[k] = dux_dx + duz_dz;
		}
	}
}

void SolventFlow::solvePressure()
{
	double delta_P_dxdx = pressure_difference/(dx*dx);
	for (int j=0; j<nz; j++){
		for (int i=0; i<nx; i++){
			int k = i+nx*j;
			b(k) = (1/d_tau)*div_u_sol_ast[k];
			if (i == 0) {
				b(k) += -delta_P_dxdx;
			} else if (i == nx-1) {
				b(k) += delta_P_dxdx;
			}
		}
	}
	x = psolver->solve(b);
	if (psolver->info() != Eigen::Success) {
		// solving failed
		std::cerr << "solving failed" << std::endl;
		return;
	}
	double mean_p = x.mean();
	// std::cerr << "mean_p " << mean_p << std::endl;
	for (int i=0; i<n; i++) {
		pressure[i] = x(i)-mean_p;
	}
}

void SolventFlow::updateSolventFlow()
{
	// -grad P = u - vp
	// u = vp -grad P
	int im1, jm1;
	double pd;
	for (int j=0; j<nz; j++){
		for (int i=0; i<nx; i++){
			int k = q(i,j);
			im1 = i-1;
			pd = 0;
			if (im1 == -1) {
				im1 = nx-1;
				pd = pressure_difference;
			}
			jm1 = j-1;
			if (jm1 == -1) {
				jm1 = nz-1;
			}
			u_sol_x[k] = u_sol_ast_x[k] - d_tau*(pressure[k] - (pressure[q(im1,j)]+pd))/dx;
			u_sol_z[k] = u_sol_ast_z[k] - d_tau*(pressure[k] -  pressure[q(i,jm1)])/dz;
//			std::cerr << k << ' ' << q(im1,j) << ' ' << q(i, jm1) << ' ' << pressure[k] - pressure[q(im1,j)] << ' ' << pressure[k] - pressure[q(i,jm1)] << std::endl;
		}
	}
}

vec3d SolventFlow::localFlow(const vec3d &p)
{
	int i = p.x/dx;
	int j = p.z/dz;
	int k = q(i, j);
	double x_diff = (p.x-(pos[k].x-dx/2))/dx;
	double z_diff = (p.z-(pos[k].z-dz/2))/dz;
	double vx = u_sol_x[k] + (u_sol_x[q(i+1,j)]-u_sol_x[k])*x_diff;
	double vz = u_sol_z[k] + (u_sol_z[q(i,j+1)]-u_sol_z[k])*z_diff;
	return vec3d(vx, 0, vz);
}

double SolventFlow::meanVelocity()
{
	double mean_velocity = 0;
	for (int k=0; k<n; k++) {
		mean_velocity += u_sol_x[k] ;
	}
	return mean_velocity/n;
}

void SolventFlow::velocityProfile(std::ofstream &fout_fp)
{
	for (int j=0; j<nz; j++){
		double total_ux_sol = 0;
		double total_ux_na = 0;
		for (int i=0; i<nx; i++){
			int k = q(i,j);
			total_ux_sol += u_sol_x[q(i,j)];
			total_ux_na += u_x[q(i,j)];
		}
		double mean_ux_sol = total_ux_sol/nx;
		double mean_ux_na = total_ux_na/nx;
		
		fout_fp << pos[q(0,j)].z << ' ' << mean_ux_sol << ' ' << mean_ux_na+mean_ux_sol << std::endl;
	}
	fout_fp << std::endl;
	
	
	
}


void SolventFlow::outputYaplot(std::ofstream &fout_flow)
{
	static bool first = true;
	if (first) {
		first = false;
	} else {
		fout_flow << std::endl;
	}
	fout_flow << "y 1" << std::endl;
	fout_flow << "@ 0" << std::endl;

	for (int i=0; i < sys->get_np(); i++) {
		double x = sys->position[i].x-sys->get_lx()/2;
		double z = sys->position[i].z-sys->get_lz()/2;
		double a = sys->radius[i];
		fout_flow << "r " << a << std::endl;
		fout_flow << "c " << x  << ' ' << 0 << ' ' << z << std::endl;
	}
	fout_flow << "y 2" << std::endl;
	fout_flow << "@ 3" << std::endl;
	fout_flow << "r " << 0.1 << std::endl;
	double vfactor = 1;
	vec3d o_shift(sys->get_lx()/2, 0, sys->get_lz()/2);
	if (0) {
		for (int j=0; j<nz; j++){
			for (int i=0; i<nx; i++){
				vec3d po = pos[i + nx*j] - o_shift;
				//	double vx = vfactor*u_particle_x[i + nx*j];
				//	double vz = vfactor*u_particle_z[i + nx*j];
				//			double vx = vfactor*u_sol_x[i + nx*j];
				//		double vz = vfactor*u_sol_z[i + nx*j];
				//double vx = vfactor*u_sol_x[i + nx*j];
				//double vz = vfactor*u_sol_z[i + nx*j];
				
				double vx = vfactor*u_x[i + nx*j];
				double vz = vfactor*u_z[i + nx*j];
				
				fout_flow << "s ";
				fout_flow << po.x - dx/2      << ' ' << -0.02 << ' ' << po.z << ' ';
				fout_flow << po.x - dx/2 + vx << ' ' << -0.02 << ' ' << po.z << std::endl;
				fout_flow << "s ";
				fout_flow << po.x << ' ' << -0.02 << ' ' << po.z - dz/2 << ' ';
				fout_flow << po.x << ' ' << -0.02 << ' ' << po.z - dz/2 + vz << std::endl;
			}
		}
	}
	fout_flow << "y 2" << std::endl;
	fout_flow << "@ 3" << std::endl;
	fout_flow << "r " << 0.1 << std::endl;
	for (int j=0; j<nz; j++){
		for (int i=0; i<nx; i++){
			vec3d po = pos[i + nx*j] - o_shift;
			int k = q(i,j);
			double vx = vfactor*(u_sol_x[k] +u_sol_x[q(i+1,j)]);
			double vz = vfactor*(u_sol_z[k] +u_sol_z[q(i,j+1)]);
			fout_flow << "s ";
			fout_flow << po.x      << ' ' << -0.02 << ' ' << po.z << ' ';
			fout_flow << po.x + vx << ' ' << -0.02 << ' ' << po.z + vz << std::endl;
		}
	}
	fout_flow << "y 3" << std::endl;
	fout_flow << "@ 4" << std::endl;
	fout_flow << "r " << 0.1 << std::endl;
	for (int j=0; j<nz; j++){
		for (int i=0; i<nx; i++){
			vec3d po = pos[i + nx*j] - o_shift;
			int k = q(i,j);
			double vx = vfactor*(u_x[k] +u_x[q(i+1,j)]);
			double vz = vfactor*(u_z[k] +u_z[q(i,j+1)]);
			fout_flow << "s ";
			fout_flow << po.x      << ' ' << -0.02 << ' ' << po.z << ' ';
			fout_flow << po.x + vx << ' ' << -0.02 << ' ' << po.z + vz << std::endl;
		}
	}

	
	
	if (1) {
		fout_flow << "y 5" << std::endl;
		for (int j=0; j<nz; j++){
			for (int i=0; i<nx; i++){
				vec3d po = pos[i + nx*j] - o_shift;
				double s = phi[i + nx*j];
				//double s = 100*div_u_sol_ast[i+nx*j];
				if (s > 0 ) {
					fout_flow << "@ 4" << std::endl;
				} else {
					fout_flow << "@ 6" << std::endl;
				}
				fout_flow << "r " << abs(s) << std::endl;
				fout_flow << "c ";
				fout_flow << po.x << ' ' << -0.02 << ' ' << po.z << std::endl;
			}
		}
	}
	
	fout_flow << "y 3" << std::endl;
	fout_flow << "@ 5" << std::endl;
	fout_flow << "r " << 0.05 << std::endl;
	for (int i=0; i<nx; i++){
		double x = i*dx-o_shift.x;
		fout_flow << "s ";
		fout_flow << x << ' ' << -0.02 << ' ' <<              -o_shift.z << ' ';
		fout_flow << x << ' ' << -0.02 << ' ' << sys->get_lz()-o_shift.z << std::endl;
	}
	for (int i=0; i<nz; i++){
		double z =  i*dz-o_shift.z;
		fout_flow << "s ";
		fout_flow <<              -o_shift.x << ' ' << -0.02 << ' ' << z << ' ';
		fout_flow << sys->get_lx()-o_shift.x << ' ' << -0.02 << ' ' << z << std::endl;
	}

	fout_flow << "y 4" << std::endl;
	
	double p_max = 0;
	for (int i=0; i< n; i++) {
		if (p_max < pressure[i]) {
			p_max = pressure[i];
		}
	}
	for (int j=0; j<nz; j++){
		for (int i=0; i<nx; i++){
			//			double s =div_u_particle[i+nx*j];
			double s = pressure[i+nx*j]/pressure_difference;
			if (s > 0 ) {
				fout_flow << "@ 4" << std::endl;
			} else {
				fout_flow << "@ 6" << std::endl;
			}
			vec3d po = pos[i + nx*j] - o_shift;
			fout_flow << "r "<< abs(s) << std::endl;
			fout_flow << "c " << po.x << " -0.03 " << po.z << std::endl;
		}
	}
}
