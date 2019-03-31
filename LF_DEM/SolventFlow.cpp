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

SolventFlow::SolventFlow():
settling(false),
channel_flow(false)
{
	psolver = new Eigen::SimplicialLDLT <SpMat>;
}

SolventFlow::~SolventFlow()
{
	delete [] psolver;
}

void SolventFlow::init(System* sys_, std::string simulation_type)
{
	sys = sys_;
	if (simulation_type == "settling") {
		settling = true;
	} else if (simulation_type == "channel flow") {
		channel_flow = true;
	} else {
		std::ostringstream error_str;
		error_str << "Incorrect simulation type\n";
		throw std::runtime_error(error_str.str());
	}
	
	if (settling) {
		std::cerr << "settling simulation" << std::endl;
		pressure_difference = 0;
		sys->body_force = true;
	}
	if (channel_flow) {
		pressure_difference = 0;
		target_flux = 1;
		sys->body_force = false;
	}
	viscosity = 1;
	nx = sys->p.mesh_nx;
	nz = sys->p.mesh_nz;
	n = nx*nz;
	dx = sys->get_lx()/nx;
	dz = sys->get_lz()/nz;
	cell_area = dx*dz;
	numerical_Re = 0.01;
	smooth_length = dx/2;
	sq_smooth_length = smooth_length*smooth_length;
	pressure.resize(n, 0);
	u_diff_x.resize(n, 0);
	u_sol_x.resize(n, 0);
	u_sol_ast_x.resize(n, 0);
	div_u_sol_ast.resize(n, 0);
	phi_ux.resize(n,0);
	pos.resize(n);
	if (sys->p.boundary_conditions == 0) {
		jmax_uz = nz;
		u_diff_z.resize(n, 0);
		u_sol_z.resize(n, 0);
		u_sol_ast_z.resize(n, 0);
		phi_uz.resize(n, 0);
		omega.resize(n, 0);
		strain_rate_xx.resize(n, 0);
		strain_rate_xz.resize(n, 0);
		strain_rate_zz.resize(n, 0);
	} else {
		ux_top = 0;
		ux_bot = 0;
		jmax_uz = nz+1;
		u_diff_z.resize(nx*(nz+1), 0);
		u_sol_z.resize(nx*(nz+1), 0);
		u_sol_ast_z.resize(nx*(nz+1), 0);
		phi_uz.resize(nx*(nz+1), 0);
		omega.resize(nx*(nz+1), 0);
		strain_rate_xx.resize(nx*(nz+1), 0);
		strain_rate_xz.resize(nx*(nz+1), 0);
		strain_rate_zz.resize(nx*(nz+1), 0);
	}
	for (int j=0; j<nz; j++) {
		for (int i=0; i<nx; i++) {
			pos[i+nx*j].set(dx*i+dx/2, 0, dz*j+dz/2);
		}
	}
	initPoissonSolver();
}

void SolventFlow::initPoissonSolver()
{
	std::vector< std::vector<double> > lmat;
	lmat.resize(n);
	for (int l=0; l<n; l++) {
		lmat[l].resize(n, 0);
	}
	//	double exz = 2*(1/(dx*dx)+1/(dz*dz));
	double dx2i = 1/(dx*dx);
	double dz2i = 1/(dz*dz);
	for (int j = 0; j < nz; j++) {
		for (int i = 0; i < nx; i++) {
			int k = i+j*nx;
			if (sys->p.boundary_conditions == 0) {
				// b(xi, zi)
				lmat[k][k]              += -2*(dx2i+dz2i);  //p(xi,zi) --> b(xi,zi)
				lmat[k][meshNb(i+1, j)] += dx2i;  //p(xi+1,zi) --> b(xi,zi)
				lmat[k][meshNb(i-1, j)] += dx2i;  //p(xi-1,zi) --> b(xi,zi)
				lmat[k][meshNb(i, j+1)] += dz2i;  //p(xi,zi+1) --> b(xi,zi)
				lmat[k][meshNb(i, j-1)] += dz2i;  //p(xi,zi-1) --> b(xi,zi)
			} else if (sys->p.boundary_conditions == 1) {
				if (j == 0) {
					// dp/dz = 0 at the bottom wall
					lmat[k][k]              += -2*dx2i-dz2i;
					lmat[k][meshNb(i+1, j)] += dx2i;  //p(xi+1,zi) --> b(xi,zi)
					lmat[k][meshNb(i-1, j)] += dx2i;  //p(xi-1,zi) --> b(xi,zi)
					lmat[k][meshNb(i, j+1)] += dz2i;  //p(xi,zi+1) --> b(xi,zi)
					// lmat[k][meshNb(i, j-1)] += 0;
				} else if (j == nz-1) {
					// dp/dz = 0 at the top wall
					lmat[k][k]              += -2*dx2i-dz2i;
					lmat[k][meshNb(i+1, j)] += dx2i;  //p(xi+1,zi) --> b(xi,zi)
					lmat[k][meshNb(i-1, j)] += dx2i;  //p(xi-1,zi) --> b(xi,zi)
					// lmat[k][meshNb(i, j+1)] += 0;
					lmat[k][meshNb(i, j-1)] += dz2i;  //j = nz-1, k = i+(nz-1)*nx =
				} else {
					lmat[k][k]              += -2*dx2i-2*dz2i;  //p(xi,zi) --> b(xi,zi)
					lmat[k][meshNb(i+1, j)] += dx2i;  //p(xi+1,zi) --> b(xi,zi)
					lmat[k][meshNb(i-1, j)] += dx2i;  //p(xi-1,zi) --> b(xi,zi)
					lmat[k][meshNb(i, j+1)] += dz2i; // j = [1, nz-2]
					lmat[k][meshNb(i, j-1)] += dz2i; // j = [1, nz-2]
				}
			}
		}
	}
	Eigen::initParallel();
	//Eigen::setNbThreads(4);
	// std::cerr << "number of threads = " << Eigen::nbThreads() << std::endl;
	//	std::fill(lmat[0].begin(), lmat[0].end(), 0);
	//lmat[0][0] = 1;
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
	std::fill(u_diff_x.begin(), u_diff_x.end(), 0);
	std::fill(u_diff_z.begin(), u_diff_z.end(), 0);
	std::fill(phi_ux.begin(), phi_ux.end(), 0);
	std::fill(phi_uz.begin(), phi_uz.end(), 0);
	for (int i=0; i < sys->get_np(); i++) {
		double x = sys->position[i].x;
		double z = sys->position[i].z;
		double ud_x = sys->na_velocity[i].x;
		double ud_z = sys->na_velocity[i].z;
		double radius = sys->radius[i];
		double particle_volume = M_PI*radius*radius;

		int ix = x/dx;
		int iz = z/dz;
		mesh_nb.clear();
		//phi_mesh_nb.clear();
		udx_values.clear();
		udz_values.clear();
		phi_ux_values.clear();
		phi_uz_values.clear();
		double total_weight_udx = 0;
		double total_weight_udz = 0;
		for (int m=-1; m<=2; m++) {
			for (int l=-1; l<=2; l++) {
				int iix = ix+l;
				int iiz = iz+m;
				double xo = x;
				double zo = z;
				if (iix <= -1) {
					iix += nx;
					xo += sys->get_lx();
				} else if (iix >= nx) {
					iix -= nx;
					xo -= sys->get_lx();
				}
				if (sys->p.boundary_conditions == 0) {
					if (iiz <= -1) {
						iiz += nz;
						zo += sys->get_lz();
					} else if (iiz >= nz) {
						iiz -= nz;
						zo -= sys->get_lz();
					}
				}
				if (iiz >= 0 && iiz < jmax_uz) {
					int k = iix + iiz*nx;
					mesh_nb.push_back(k);
					double xx = (xo-pos[k].x)*(xo-pos[k].x);
					double zz = (zo-pos[k].z)*(zo-pos[k].z);
					double xu = iix*dx;
					double zu = iiz*dz;
					double w_udx = (m != 2 ? weightFunc((xo-xu)*(xo-xu) + zz) : 0);
					total_weight_udx += w_udx;
					udx_values.push_back(ud_x*w_udx);
					phi_ux_values.push_back(w_udx);
					double w_udz = (l != 2 ? weightFunc(xx + (zo-zu)*(zo-zu)) : 0);
					total_weight_udz += w_udz;
					udz_values.push_back(ud_z*w_udz);
					phi_uz_values.push_back(w_udz);
				}
			}
		}
		for (int m=0; m< mesh_nb.size(); m++) {
			int k = mesh_nb[m];
			if (k < n) {
				u_diff_x[k] += udx_values[m]/total_weight_udx;
				phi_ux[k] += (particle_volume*phi_ux_values[m]/cell_area)/total_weight_udx;
			}
			u_diff_z[k] += udz_values[m]/total_weight_udz;
			phi_uz[k] += (particle_volume*phi_uz_values[m]/cell_area)/total_weight_udz;
		}
	}
}

void SolventFlow::update(double pressure_difference_)
{
	static std::ofstream fout_tmp("debug.dat");
	double flux = calcFlux();
	if (channel_flow) {
		//pressure_difference = 10*pressure_difference_;
		double pd_increment = 1e-4;
		if (flux < target_flux) {
			pressure_difference += pd_increment;
		} else {
			pressure_difference -= pd_increment;
		}
	}
	d_tau = sys->dt/numerical_Re;
	particleVelocityDiffToMesh();
	predictorStep();
	calcVelocityDivergence();
	solvePressure();
	correctorStep();
	calcOmega();
	fout_tmp << std::endl;
}

double SolventFlow::porousResistance(double area_fraction)
{
	//static const double phi_max = 1.2;
	/*
	 static const double phi_st = 0.8;
	 static const double inv_res_max = (phi_max-phi_st)/phi_st;
	 return inv_res_max*volume_fraction/(phi_max-volume_fraction);
	 */
	return (27.0/8)*area_fraction;
}

void SolventFlow::predictorStep()
{
	double dx2 = dx*dx;
	double dz2 = dz*dz;
	for (int j=0; j<nz; j++) {
		int jp1, jm1;
		if (sys->p.boundary_conditions == 0) {
			jp1 = (j == nz-1 ? 0 : j+1);
			jm1 = (j == 0 ? nz-1 : j-1);
		} else {
			jp1 = j+1;
			jm1 = j-1;
		}
		for (int i=0; i<nx; i++) {
			int k = i + nx*j;
			int ip1 = (i == nx-1 ? 0 : i+1);
			int im1 = (i == 0 ? nx-1 : i-1);
			int ir = ip1 + j  *nx; //right
			int il = im1 + j  *nx; //left
			int ju = i   + jp1*nx; //up
			int jd = i   + jm1*nx; //down
			double res_coeff_ux = porousResistance(phi_ux[k]);
			double res_coeff_uz = porousResistance(phi_uz[k]);
			if (sys->p.boundary_conditions == 0) {
				/* periodic boundary condtions in x and z directions.
				 */
				double dd_ux = (u_sol_x[ir]-2*u_sol_x[k]+u_sol_x[il])/dx2+(u_sol_x[ju]-2*u_sol_x[k]+u_sol_x[jd])/dz2;
				double dd_uz = (u_sol_z[ir]-2*u_sol_z[k]+u_sol_z[il])/dx2+(u_sol_z[ju]-2*u_sol_z[k]+u_sol_z[jd])/dz2;
				double fx = viscosity*dd_ux + res_coeff_ux*u_diff_x[k];
				double fz = viscosity*dd_uz + res_coeff_uz*u_diff_z[k];
				u_sol_ast_x[k] = u_sol_x[k] + d_tau*fx;
				u_sol_ast_z[k] = u_sol_z[k] + d_tau*fz;
			} else {
				/* periodic boundary condtions in x directions.
				 * Non-slip boundy conditions at z = 0 and z=Lz.
				 */
				if (j == 0) {
					// bottom
					double ux_bot_m1 = 2*ux_bot - u_sol_x[k]; // (ux(i,-1)+ux(i,0))/2 = ux_bot
					double dd_ux = (u_sol_x[ir]-2*u_sol_x[k]+u_sol_x[il])/dx2+(u_sol_x[ju]-2*u_sol_x[k]+ux_bot_m1)/dz2;
					double fx = viscosity*dd_ux + res_coeff_ux*u_diff_x[k];
					u_sol_ast_x[k] = u_sol_x[k] + d_tau*fx;
					u_sol_ast_z[k] = 0; // Just on the bottom wall
				} else {
					double dd_ux, dd_uz;
					if (j == nz-1) {
						// top cell (The upper segments are top wall)
						/* (ux_top_p1 + u_sol_x[k])/2 = ux_top
						 * --->
						 * ux_top_p1 = 2*ux_top - u_sol_x[k];
						 */
						double ux_top_p1 = 2*ux_top - u_sol_x[k];
						dd_ux = (u_sol_x[ir]-2*u_sol_x[k]+u_sol_x[il])/dx2+(ux_top_p1-2*u_sol_x[k]+u_sol_x[jd])/dz2;
					} else {
						dd_ux = (u_sol_x[ir]-2*u_sol_x[k]+u_sol_x[il])/dx2+(u_sol_x[ju]-2*u_sol_x[k]+u_sol_x[jd])/dz2;
					}
					
					dd_uz = (u_sol_z[ir]-2*u_sol_z[k]+u_sol_z[il])/dx2+(u_sol_z[ju]-2*u_sol_z[k]+u_sol_z[jd])/dz2;
					double fx = viscosity*dd_ux + res_coeff_ux*u_diff_x[k];
					double fz = viscosity*dd_uz + res_coeff_uz*u_diff_z[k];
					u_sol_ast_x[k] = u_sol_x[k] + d_tau*fx;
					u_sol_ast_z[k] = u_sol_z[k] + d_tau*fz;
				}
			}
		}
	}
	if (sys->p.boundary_conditions == 1) {
		int j = nz;
		for (int i=0; i<nx; i++) {
			int k = i + nx*j;
			u_sol_ast_z[k] = 0;
		}
	}
}

void SolventFlow::calcVelocityDivergence()
{
	// This works for both boundary conditions.
	// u_sol_ast_z[i, j = 0] = 0 (top and bottom)
	if (sys->p.boundary_conditions == 0) {
		for (int j=0; j<nz; j++) {
			int jp1 = (j == nz-1 ? 0 : j+1);
			for (int i=0; i<nx; i++) {
				int ip1 = (i == nx-1 ? 0 : i+1);
				int k = i + j*nx;
				double dux_dx = (u_sol_ast_x[ip1 + j*nx] - u_sol_ast_x[k])/dx;
				double duz_dz = (u_sol_ast_z[i + jp1*nx] - u_sol_ast_z[k])/dz;
				div_u_sol_ast[k] = dux_dx + duz_dz;
			}
		}
	} else {
		for (int j=0; j<nz; j++) {
			int jp1 = j+1;
			for (int i=0; i<nx; i++) {
				int ip1 = (i == nx-1 ? 0 : i+1);
				int k = i + j*nx;
				double dux_dx = (u_sol_ast_x[ip1 + j*nx] - u_sol_ast_x[k])/dx;
				double duz_dz = (u_sol_ast_z[i + jp1*nx] - u_sol_ast_z[k])/dz;
				div_u_sol_ast[k] = dux_dx + duz_dz;
			}
		}
	}
}

void SolventFlow::solvePressure()
{
	double delta_P_dxdx = pressure_difference/(dx*dx);
	for (int k=0; k<n; k++) {
		b(k) = (1/d_tau)*div_u_sol_ast[k];
	}
	if (pressure_difference != 0) {
		for (int j=0; j<nz; j++) {
			int nx_j = nx*j;
			// i = 0 --> k = nx*j
			b(nx_j) += -delta_P_dxdx;
			// i = nx-1 ---> k = nx-1+nx*j
			b(nx-1+nx_j) += delta_P_dxdx;
		}
	}
	x = psolver->solve(b);
	if (psolver->info() != Eigen::Success) {
		// solving failed
		std::cerr << "solving failed" << std::endl;
		return;
	}
	double mean_p = x.mean();
	for (int i=0; i<n; i++) {
		pressure[i] = x(i)-mean_p;
	}
}

void SolventFlow::correctorStep()
{
	if (sys->p.boundary_conditions == 0) {
		/* periodic boundary condtions in x and z directions.
		 */
		for (int i=0; i<nx; i++) {
			int im1 = (i == 0 ? nx-1 : i-1);
			double pd = (i == 0 ? pressure_difference : 0);
			for (int j=0; j<nz; j++) {
				int k = i+nx*j;
				int jm1 = (j == 0 ? nz-1 : j-1);
				u_sol_x[k] = u_sol_ast_x[k] - d_tau*(pressure[k] - (pressure[im1+j*nx]+pd))/dx;
				u_sol_z[k] = u_sol_ast_z[k] - d_tau*(pressure[k] - pressure[i+jm1*nx])/dz;
			}
		}
	} else {
		/* periodic boundary condtions in x directions.
		 * Non-slip boundy conditions at z = 0 and z=Lz.
		 * (uz[i, j = nz] = 0,)
		 */
		for (int i=0; i<nx; i++) {
			int im1 = (i == 0 ? nx-1 : i-1);
			double pd = (i == 0 ? pressure_difference : 0);
			{// bottom
				u_sol_x[i] = u_sol_ast_x[i] - d_tau*(pressure[i] - (pressure[im1]+pd))/dx;
				u_sol_z[i] = 0;
			}
			for (int j=1; j<nz; j++) {
				int k = i+nx*j;
				u_sol_x[k] = u_sol_ast_x[k] - d_tau*(pressure[k] - (pressure[im1+j*nx]+pd))/dx;
				u_sol_z[k] = u_sol_ast_z[k] - d_tau*(pressure[k] - pressure[i+(j-1)*nx])/dz;
			}
			{// top
				int j=nz;
				u_sol_z[i+j*nx] = 0;
			}
		}
	}
}

void SolventFlow::calcOmega()
{
	if (sys->p.boundary_conditions == 0) {
		for (int j=0; j< nz; j++) {
			int jp1 = (j == nz-1 ? 0 : j+1);
			int jm1 = (j == 0 ? nz-1 : j-1);
			for (int i=0; i<nx; i++) {
				int k =i+nx*j;
				int im1 = (i == 0 ? nx-1 : i-1);
				int ip1 = (i == nx-1 ? 0 : i+1);
				double d_ux_d_x = (u_sol_x[ip1+j*nx]-u_sol_x[k])/dx;
				double d_uz_d_z = (u_sol_x[i+jp1*nx]-u_sol_x[k])/dz;
				double d_ux_d_z = (u_sol_x[k]-u_sol_x[i+jm1*nx])/dz;
				double d_uz_d_x = (u_sol_z[k]-u_sol_z[im1+nx*j])/dx;
				omega[k] = -(d_uz_d_x-d_ux_d_z)/2; // left_bottom
				strain_rate_xx[k] = d_ux_d_x; // center
				strain_rate_xz[k] = (d_ux_d_z + d_uz_d_x)/2; // left_bottom
				strain_rate_zz[k] = d_uz_d_z; // center
			}
		}
	} else {
		for (int j=0; j<=nz; j++) {
			int jp1 = j+1;
			int jm1 = j-1;
			for (int i=0; i<nx; i++) {
				int k =i+nx*j;
				int ip1 = (i == nx-1 ? 0 : i+1);
				if (j == 0) {
					// (u[j] + u[j-1])/2 = ux_bot
					// u[j-1] = 2*ux_bot - u[j];
					// u[j] - u[j-1] = u[j] -(2*ux_bot - u[j]) = 2*u[j] -2*ux_bot
					double d_ux_d_z = 2*(u_sol_x[k]-ux_bot)/dz; // left_bottom
					double d_uz_d_x = 0; // left_bottom (uz[j=0]= 0)
					double d_ux_d_x = (u_sol_x[ip1]-u_sol_x[k])/dx; //center
					double d_uz_d_z = (u_sol_z[i+nx]-u_sol_z[k])/dz; //center
					omega[k] = -(d_uz_d_x-d_ux_d_z)/2;
					strain_rate_xx[k] = d_ux_d_x; // center
					strain_rate_xz[k] = (d_ux_d_z + d_uz_d_x)/2; // left_bottom
					strain_rate_zz[k] = d_uz_d_z; // center
				} else if (j == nz) {
					// (u[j]+u[j-1])/2 = ux_top
					// u[j] = 2ux_top - u[j-1]
					// u[j]-u[j-1] = 2ux_top - u[j-1] - u[j-1]
					//
					int k =i+nx*j;
					int im1 = (i == 0 ? nx-1 : i-1);
					double d_ux_d_x = (-u_sol_x[i+jm1*nx]+u_sol_x[im1+jm1*nx])/dx;
					double d_uz_d_z = 0;
					double d_ux_d_z = 2*(ux_top-u_sol_x[i+jm1*nx])/dz; // left_bottom
					double d_uz_d_x = 0; // left_bottom
					omega[k] = -(d_uz_d_x-d_ux_d_z)/2; // left_bottom
					strain_rate_xx[k] = d_ux_d_x; // center
					strain_rate_xz[k] = (d_ux_d_z + d_uz_d_x)/2; // left_bottom
					strain_rate_zz[k] = d_uz_d_z; // center
				} else {
					int k =i+nx*j;
					int im1 = (i == 0 ? nx-1 : i-1);
					int ip1 = (i == nx-1 ? 0 : i+1);
					double d_ux_d_x = (u_sol_x[ip1+j*nx]-u_sol_x[k])/dx;
					double d_uz_d_z = (u_sol_z[i+jp1*nx]-u_sol_z[k])/dz;
					double d_ux_d_z = (u_sol_x[k]-u_sol_x[i+jm1*nx])/dz;
					double d_uz_d_x = (u_sol_z[k]-u_sol_z[im1+nx*j])/dx;
					omega[k] = -(d_uz_d_x-d_ux_d_z)/2; // left_bottom
					strain_rate_xx[k] = d_ux_d_x; // center
					strain_rate_xz[k] = (d_ux_d_z + d_uz_d_x)/2; // left_bottom
					strain_rate_zz[k] = d_uz_d_z; // center
				}
			}
		}
	}
}

void SolventFlow::localFlow(const vec3d &p,
							vec3d &u_local,
							vec3d &omega_local,
							std::vector<double> &e_local)
{
	double x_dx = p.x/dx;
	double z_dz = p.z/dz;
	int i = (int)x_dx;
	int j = (int)z_dz;
	int k = i+j*nx;
	int ip1 = (i == nx-1 ? 0 : i+1);
	double x_diff = x_dx - i;
	double z_diff = z_dz - j;
	u_local.x = u_sol_x[k] + (u_sol_x[ip1+j*nx]-u_sol_x[k])*x_diff;
	double sxz = strain_rate_xz[k];
	if (sys->p.boundary_conditions == 0) {
		int jp1 = (j == nz-1 ? 0 : j+1);
		u_local.z = u_sol_z[k] + (u_sol_z[i+jp1*nx]-u_sol_z[k])*z_diff;
		omega_local.y = omega[k] + (omega[ip1+j*nx]-omega[k])*x_diff + (omega[i+jp1*nx]-omega[k])*z_diff;
		e_local[1] = sxz + (strain_rate_xz[ip1+j*nx]-sxz)*x_diff + (strain_rate_xz[i+jp1*nx]-sxz)*z_diff;
	} else {
		int jp1 = j+1;
		u_local.z = u_sol_z[k] + (u_sol_z[i+jp1*nx]-u_sol_z[k])*z_diff;
		omega_local.y = omega[k] + (omega[ip1+j*nx]-omega[k])*x_diff + (omega[i+jp1*nx]-omega[k])*z_diff;
		e_local[1] = sxz + (strain_rate_xz[ip1+j*nx]-sxz)*x_diff + (strain_rate_xz[i+jp1*nx]-sxz)*z_diff;
	}
	e_local[0] = strain_rate_xx[k];
	e_local[2] = strain_rate_zz[k];
	
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
	
	double vfactor = 1;
	vec3d o_shift(sys->get_lx()/2, 0, sys->get_lz()/2);
	if (0) {
		fout_flow << "y 2" << std::endl;
		fout_flow << "@ 3" << std::endl;
		fout_flow << "r " << 0.1 << std::endl;
		for (int j=0; j<=nz; j++){
			for (int i=0; i<nx; i++){
				vec3d po = pos[i + nx*j] - o_shift;
				//	double vx = vfactor*u_particle_x[i + nx*j];
				//	double vz = vfactor*u_particle_z[i + nx*j];
				//			double vx = vfactor*u_sol_x[i + nx*j];
				//		double vz = vfactor*u_sol_z[i + nx*j];
				//double vx = vfactor*u_sol_x[i + nx*j];
				//double vz = vfactor*u_sol_z[i + nx*j];
				if (j < nz) {
					double vx = vfactor*u_sol_x[i + nx*j];
					fout_flow << "s ";
					fout_flow << po.x - dx/2      << ' ' << -0.02 << ' ' << po.z << ' ';
					fout_flow << po.x - dx/2 + vx << ' ' << -0.02 << ' ' << po.z << std::endl;
				}
				double vz = vfactor*u_sol_z[i + nx*j];
				fout_flow << "s ";
				fout_flow << po.x << ' ' << -0.02 << ' ' << po.z - dz/2 << ' ';
				fout_flow << po.x << ' ' << -0.02 << ' ' << po.z - dz/2 + vz << std::endl;
			}
		}
		fout_flow << "y 3" << std::endl;
		fout_flow << "@ 4" << std::endl;
		fout_flow << "r " << 0.1 << std::endl;
		for (int j=0; j<nz; j++){
			for (int i=0; i<nx; i++){
				vec3d po = pos[i + nx*j] - o_shift;
				int k = meshNb(i,j);
				if (j < nz) {
					double vx = vfactor*u_diff_x[k];
					fout_flow << "s ";
					fout_flow << po.x - dx/2      << ' ' << -0.02 << ' ' << po.z << ' ';
					fout_flow << po.x - dx/2 + vx << ' ' << -0.02 << ' ' << po.z << std::endl;
				}
				double vz = vfactor*u_diff_z[k];
				fout_flow << "s ";
				fout_flow << po.x << ' ' << -0.02 << ' ' << po.z - dz/2 << ' ';
				fout_flow << po.x << ' ' << -0.02 << ' ' << po.z - dz/2 + vz << std::endl;
			}
		}
		
	} else {
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
	}
	
	if (0) {
		fout_flow << "y 5" << std::endl;
//		double cell_area = dx*dz;
		fout_flow << "@ 10" << std::endl;
		
		for (int j=0; j<nz; j++){
			for (int i=0; i<nx; i++){
				vec3d po = pos[i + nx*j] - o_shift;
				// double r = sqrt(cell_area*phi[i + nx*j]/M_PI);
				//double r = porousResistance(phi[i + nx*j]);
				//double s = 100*div_u_sol_ast[i+nx*j];
				//				if (r > 0 ) {
				//					fout_flow << "@ 4" << std::endl;
				//				} else {
				//					fout_flow << "@ 6" << std::endl;
				//				}
				//fout_flow << "r " << abs(r) << std::endl;
				fout_flow << "c ";
				fout_flow << po.x << ' ' << -0.02 << ' ' << po.z << std::endl;
			}
		}
	}
	
	fout_flow << "y 9" << std::endl;
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
	double pr_max = 0;
	for (auto pr : pressure) {
		if (abs(pr) > pr_max) {
			pr_max = abs(pr);
		}
	}
	for (int j=0; j<nz; j++){
		for (int i=0; i<nx; i++){
			//			double s =div_u_particle[i+nx*j];
			double s = pressure[i+nx*j]/pr_max;
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

