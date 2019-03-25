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
	pressure_difference = 0;
	sys = sys_;
	nx = 20;
	nz = 20;
	n = nx*nz;
	dx = sys->get_lx()/nx;
	dz = sys->get_lz()/nz;
	smooth_length = dx/2;
	sq_smooth_length = smooth_length*smooth_length;
	pressure.resize(n);
	u_diff_x.resize(n);
	u_diff_z.resize(n);
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
	for (int j = 0; j < nz; j++){
		for (int i = 0; i < nx; i++){
			// b(xi, zi)
			int k = i+j*nx;
			lmat[k][k] = -exz;  //p(xi,zi) --> b(xi,zi)
			lmat[k][meshNb(i+1, j)] = ex;  //p(xi+1,zi) --> b(xi,zi)
			lmat[k][meshNb(i-1, j)] = ex;  //p(xi-1,zi) --> b(xi,zi)
			lmat[k][meshNb(i, j+1)] = ez;  //p(xi,zi+1) --> b(xi,zi)
			lmat[k][meshNb(i, j-1)] = ez;  //p(xi,zi-1) --> b(xi,zi)
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

double SolventFlow::weightFunc(double r_sq)
{
	return exp(-r_sq/sq_smooth_length);
}

void SolventFlow::particleVelocityDiffToMesh()
{
	double cell_area = dx*dz;
	std::vector <int> mesh_nb;
	std::vector <double> udx_values;
	std::vector <double> udz_values;
	std::vector <double> phi_values;
	for (int k=0; k<n; k++) {
		u_diff_x[k] = 0;
		u_diff_z[k] = 0;
		phi[k] = 0;
	}
	for (int i=0; i < sys->get_np(); i++) {
		double x = sys->position[i].x;
		double z = sys->position[i].z;
		double ud_x = sys->na_velocity[i].x;
		double ud_z = sys->na_velocity[i].z;
		double radius = sys->radius[i];
		int ix = x /dx;
		int iz = z /dz;
		mesh_nb.clear();
		udx_values.clear();
		udz_values.clear();
		phi_values.clear();
		double total_weight_udx = 0;
		double total_weight_udz = 0;
		double total_weight_phi = 0;
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
				double w_udx  = weightFunc((xo-x_u)*(xo-x_u) + zz);
				double w_udz  = weightFunc(xx + (zo-z_u)*(zo-z_u));
				double w_phi = weightFunc(xx+zz);
				total_weight_udx += w_udx;
				total_weight_udz += w_udz;
				total_weight_phi += w_phi;
				mesh_nb.push_back(k);
				udx_values.push_back(ud_x*w_udx);
				udz_values.push_back(ud_z*w_udz);
				phi_values.push_back(w_phi);
			}
		}
		double particle_volume = M_PI*sys->radius[i]*sys->radius[i];
		for (int m=0; m<mesh_nb.size(); m++) {
			int k = mesh_nb[m];
			u_diff_x[k] += udx_values[m]/total_weight_udx;
			u_diff_z[k] += udz_values[m]/total_weight_udz;
			phi[k] += (particle_volume*phi_values[m]/cell_area)/total_weight_phi;
		}
	}
	if (false) {
		static std::ofstream fout_tmp("hoge.yap");
		fout_tmp << "y 1" << std::endl;
		fout_tmp << "@ 0" << std::endl;
		for (int i=0; i < sys->get_np(); i++) {
			double x = sys->position[i].x-sys->get_lx()/2;
			double z = sys->position[i].z-sys->get_lz()/2;
			double a = sys->radius[i];
			fout_tmp << "r " << a << std::endl;
			fout_tmp << "c " << x  << ' ' << 0 << ' ' << z << std::endl;
		}
		fout_tmp << "y 2" << std::endl;
		fout_tmp << "@ 3" << std::endl;
		double total_area = 0;
		for (int k=0; k<n; k++) {
			double x_p = pos[k].x;
			double z_p = pos[k].z;
			double r = sqrt(cell_area*phi[k]/M_PI);
			fout_tmp << "r " << r << std::endl;
			fout_tmp << "c " << x_p-sys->get_lx()/2  << ' ' << 0 << ' ' << z_p-sys->get_lz()/2 << std::endl;
			total_area += r*r*M_PI;
		}
		std::cerr << "phi = " << total_area/(sys->get_lx()*sys->get_lz()) << std::endl;
		fout_tmp << std::endl;
	}
}

void SolventFlow::update(double pressure_difference_)
{
	static std::ofstream fout_tmp("debug.dat");
	double flux = calcFlux();
	double target_flux = 1;
	//pressure_difference = 10*pressure_difference_;
	double pd_increment = 1e-5;
	if (flux < target_flux) {
		pressure_difference += pd_increment;
	} else {
		pressure_difference -= pd_increment;
	}
	double rho = 0.01;
	d_tau = sys->dt/rho;
	particleVelocityDiffToMesh();
	predictorStep();
	calcVelocityDivergence();
	solvePressure();
	correctorStep();
	fout_tmp << std::endl;
}

void SolventFlow::predictorStep()
{
	double phi_max = 0.84;
	double inv_res_max = (1-phi_max)*(1-phi_max)/phi_max;
	double coeff = 0.1*inv_res_max;
	double viscosity = 1;
	double dx2 = dx*dx;
	double dz2 = dz*dz;
	for (int j=0; j<nz; j++) {
		for (int i=0; i<nx; i++) {
			int k = i + nx*j;
			double porefrac = 1-phi[k];
			double res_coeff = coeff*phi[k]/(porefrac*porefrac);
			u_sol_ast_x[k] = u_sol_x[k] + res_coeff*d_tau*u_diff_x[k];
			//- viscous_stabiliser*u_sol_x[k]; // u_x = u_p - u_s
			if (sys->p.boundary_conditions == 1) {
				/* periodic boundary condtions in x directions.
				 * Non-slip boundy conditions at z = 0 and z=Lz.
				 *
				 */
				if (j == 0) {
					u_sol_ast_z[k] = 0;
				} else {
					u_sol_ast_z[k] = u_sol_z[k] + res_coeff*d_tau*u_diff_z[k];
					//- viscous_stabiliser*u_sol_z[k];
				}
			} else {
				/* periodic boundary condtions in x and z directions.
				 */
				u_sol_ast_z[k] = u_sol_z[k] + res_coeff*d_tau*u_diff_z[k];
				//- viscous_stabiliser*u_sol_z[k];
			}
			if (true) {
				int ip1 = i+1;
				if (ip1 == nx) {
					ip1 = 0;
				}
				int im1 = i-1;
				if (im1 == -1) {
					im1 = 0;
				}
				int jp1 = j+1;
				if (jp1 == nz) {
					jp1 = 0;
				}
				int jm1 = j-1;
				if (jm1 == -1) {
					jm1 = 0;
				}
				int ir = ip1+j*nx; //right
				int il = im1+j*nx; //left
				int ju = i+jp1*nx; //up
				int jd = i+jm1*nx; //down
				u_sol_ast_x[k] += viscosity*d_tau*((u_sol_x[ir]-2*u_sol_x[k]+u_sol_x[il])/(dx2)+(u_sol_x[ju]-2*u_sol_x[k]+u_sol_x[jd])/(dz2));
				u_sol_ast_z[k] += viscosity*d_tau*((u_sol_z[ir]-2*u_sol_z[k]+u_sol_z[il])/(dx2)+(u_sol_z[ju]-2*u_sol_z[k]+u_sol_z[jd])/(dz2));
			}
		}
	}
}

void SolventFlow::calcVelocityDivergence()
{
	double dux_dx, duz_dz;
	for (int j=0; j<nz; j++) {
		for (int i=0; i<nx; i++) {
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

void SolventFlow::correctorStep()
{
	for (int j=0; j<nz; j++){
		for (int i=0; i<nx; i++){
			int k = i+nx*j;
			int im1 = i-1;
			double pd = 0;
			if (im1 == -1) {
				im1 = nx-1;
				pd = pressure_difference;
			}
			int jm1 = j-1;
			if (jm1 == -1) {
				jm1 = nz-1;
			}
			u_sol_x[k] = u_sol_ast_x[k] - d_tau*(pressure[k] - (pressure[im1+j*nx]+pd))/dx;
			if (sys->p.boundary_conditions == 1) {
				/* periodic boundary condtions in x directions.
				 * Non-slip boundy conditions at z = 0 and z=Lz.
				 *
				 */
				if (j == 0) {
					u_sol_z[k] = 0;
				} else {
					u_sol_z[k] = u_sol_ast_z[k] - d_tau*(pressure[k] - pressure[i+jm1*nx])/dz;
				}
			} else {
				/* periodic boundary condtions in x and z directions.
				 */
				u_sol_z[k] = u_sol_ast_z[k] - d_tau*(pressure[k] - pressure[i+jm1*nx])/dz;
			}
		}
	}
}

vec3d SolventFlow::localFlow(const vec3d &p)
{
	int i = p.x/dx;
	int j = p.z/dz;
	int k = i+j*nx;
	int ip1 = i+1;
	if (ip1 == nx) {
		ip1 = 0;
	}
	int jp1 = j+1;
	if (jp1 == nz) {
		jp1 = 0;
	}
	double x_diff = (p.x-(pos[k].x-dx/2))/dx;
	double z_diff = (p.z-(pos[k].z-dz/2))/dz;
	double vx = u_sol_x[k] + (u_sol_x[ip1+j*nx]-u_sol_x[k])*x_diff;
	double vz = u_sol_z[k] + (u_sol_z[i+jp1*nx]-u_sol_z[k])*z_diff;
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

double SolventFlow::calcFlux()
{
	if (true) {
		int i=(int)(0.5*nx);
		double flux_total = 0;
		for (int j=0; j<nz; j++) {
			int k = i + nx*j;
			flux_total += u_sol_x[k]*dz;
		}
		return flux_total/sys->get_lz();
	} else {
		std::vector <double> flux(nx);
		double flux_min = 99999;
		double flux_max = 0;
		for (int i=0; i< nx; i++) {
			double flux_total = 0;
			for (int j=0; j<nz; j++) {
				int k = i + nx*j;
				flux_total += u_sol_x[k]*dz;
			}
			flux[i] = flux_total/sys->get_lz();
			if (flux[i] < flux_min) {
				flux_min = flux[i];
			}
			if (flux[i] > flux_max) {
				flux_max = flux[i];
			}
		}
		double flux_sum = 0;
		for (auto fl : flux) {
			flux_sum += fl;
		}
		double flux_mean = flux_sum/nx;
		std::cerr << (flux_max  - flux_min)/flux_mean << std::endl;
		return flux_mean;
	}
}

void SolventFlow::velocityProfile(std::ofstream &fout_fp)
{
	for (int j=0; j<nz; j++) {
		double total_ux_sol = 0;
		double total_ux_na = 0;
		for (int i=0; i<nx; i++){
			int k = i+nx*j;
			total_ux_sol += u_sol_x[k];
			total_ux_na += u_diff_x[k];
		}
		double mean_ux_sol = total_ux_sol/nx;
		double mean_ux_na = total_ux_na/nx;
		// na = up - us
		// na + us = up - us + us = up
		fout_fp << pos[meshNb(0,j)].z << ' ' << mean_ux_sol << ' ' << mean_ux_na+mean_ux_sol << std::endl;
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
	for (int i=0; i < sys->np_mobile; i++) {
		double x = sys->position[i].x-sys->get_lx()/2;
		double z = sys->position[i].z-sys->get_lz()/2;
		double a = sys->radius[i];
		fout_flow << "r " << a << std::endl;
		fout_flow << "c " << x  << ' ' << 0 << ' ' << z << std::endl;
	}
	fout_flow << "y 6" << std::endl;
	fout_flow << "@ 2" << std::endl;
	for (int i=sys->np_mobile; i < sys->get_np(); i++) {
		double x = sys->position[i].x-sys->get_lx()/2;
		double z = sys->position[i].z-sys->get_lz()/2;
		double a = sys->radius[i];
		fout_flow << "r " << a << std::endl;
		fout_flow << "c " << x  << ' ' << 0 << ' ' << z << std::endl;
		if (sys->position[i].z == 0) {
			fout_flow << "c " << x  << ' ' << 0 << ' ' << z+sys->get_lz() << std::endl;
		}
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
				
				double vx = vfactor*u_diff_x[i + nx*j];
				double vz = vfactor*u_diff_z[i + nx*j];
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
			int k = i+j*nx;
			double vx = vfactor*(u_sol_x[k] +u_sol_x[meshNb(i+1,j)]);
			double vz = vfactor*(u_sol_z[k] +u_sol_z[meshNb(i,j+1)]);
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
			int k = meshNb(i,j);
			double udx = vfactor*(u_diff_x[k] +u_diff_x[meshNb(i+1,j)]);
			double udz = vfactor*(u_diff_z[k] +u_diff_z[meshNb(i,j+1)]);
			fout_flow << "s ";
			fout_flow << po.x      << ' ' << -0.02 << ' ' << po.z << ' ';
			fout_flow << po.x + udx << ' ' << -0.02 << ' ' << po.z + udz << std::endl;
		}
	}
	
	if (1) {
		fout_flow << "y 5" << std::endl;
		double cell_area = dx*dz;
		for (int j=0; j<nz; j++){
			for (int i=0; i<nx; i++){
				vec3d po = pos[i + nx*j] - o_shift;
				double r = sqrt(cell_area*phi[i + nx*j]/M_PI);
				//double s = 100*div_u_sol_ast[i+nx*j];
				if (r > 0 ) {
					fout_flow << "@ 4" << std::endl;
				} else {
					fout_flow << "@ 6" << std::endl;
				}
				fout_flow << "r " << abs(r) << std::endl;
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

int SolventFlow::meshNb(int xi, int zi)
{
	if (xi >= nx) {
		xi -= nx;
	} else if (xi <= -1) {
		xi += nx;
	}
	if (zi >= nz) {
		zi -= nz ;
	} else if (zi <= -1) {
		zi += nz;
	}
	return xi+zi*nx;
}

