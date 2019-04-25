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
sedimentation(false),
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
	length_scale = 10;
	conv_factor = length_scale; // for d = 2 (square if d = 3)
	six_pi = 6*M_PI;
	if (simulation_type == "sedimentation") {
		sedimentation = true;
	} else if (simulation_type == "channel flow") {
		channel_flow = true;
	} else {
		std::ostringstream error_str;
		error_str << "Incorrect simulation type\n";
		throw std::runtime_error(error_str.str());
	}
	average_pressure_x.setRelaxationTime(sys->p.sflow_pcontrol_rtime);
	if (sedimentation) {
		std::cerr << "sedimentation simulation" << std::endl;
		pressure_difference_x = 0;
		sys->body_force = true;
	}
	if (channel_flow) {
		pressure_difference_x = 0;
		sys->body_force = false;
	}
	nx = sys->p.sflow_nx;
	dx = sys->get_lx()/nx;
	if (sys->p.sflow_nz == -1) {
		std::cerr << std::setprecision(16) << (sys->get_lz()/dx) << std::endl;
		nz = std::lround(sys->get_lz()/dx);
		std::cerr << "nz = " << nz << std::endl;
	} else {
		nz = sys->p.sflow_nz;
	}
	dz = sys->get_lz()/nz;
	n = nx*nz;
	if (abs(dx-dz) > 1e-8) {
		std::ostringstream error_str;
		error_str << "dx = " << dx << "  dz = " << dz << "\n";
		error_str << "dx != dz. Modify sflow_nx or sflow_nz\n";
		throw std::runtime_error(error_str.str());
	}
	cell_area = dx*dz;
	smooth_length = dx/2;
	sq_smooth_length = smooth_length*smooth_length;
	pressure.resize(n, 0);
	u_diff_x.resize(n, 0);
	u_sol_x.resize(n, 0);
	u_x.resize(n, 0);
	u_sol_ast_x.resize(n, 0);
	div_u_sol_ast.resize(n, 0);
	gr_phi_Ud_phi_div_Ud.resize(n, 0);
	phi_ux.resize(n,0);
	
	if (sys->p.sflow_boundary_conditions == 0) {
		pos.resize(n);
		jmax_uz = nz;
		u_diff_z.resize(n, 0);
		u_z.resize(n, 0);
		u_sol_z.resize(n, 0);
		u_sol_ast_z.resize(n, 0);
		phi_uz.resize(n, 0);
		omega.resize(n, 0);
		strain_rate_xx.resize(n, 0);
		strain_rate_xz.resize(n, 0);
		strain_rate_zz.resize(n, 0);
	} else {
		pos.resize(nx*(nz+1));
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
	for (int j=0; j<jmax_uz; j++) {
		for (int i=0; i<nx; i++) {
			pos[i+nx*j].set(dx*i+dx/2, 0, dz*j+dz/2);
		}
	}
	initPoissonSolver();
	
	double particle_area = 0;
	for (int i=0; i<sys->np_mobile; i++) {
		particle_area += M_PI*sys->radius[i]*sys->radius[i];
	}
	average_area_fraction = particle_area/sys->get_lx()/sys->get_lz();
	std::cerr << " area_fraction = " << average_area_fraction << std::endl;
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
			if (sys->p.sflow_boundary_conditions == 0) {
				// b(xi, zi)
				lmat[k][k]              += -2*(dx2i+dz2i);  //p(xi,zi) --> b(xi,zi)
				lmat[k][meshNb(i+1, j)] += dx2i;  //p(xi+1,zi) --> b(xi,zi)
				lmat[k][meshNb(i-1, j)] += dx2i;  //p(xi-1,zi) --> b(xi,zi)
				lmat[k][meshNb(i, j+1)] += dz2i;  //p(xi,zi+1) --> b(xi,zi)
				lmat[k][meshNb(i, j-1)] += dz2i;  //p(xi,zi-1) --> b(xi,zi)
			} else if (sys->p.sflow_boundary_conditions == 1) {
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
	rhs_vector.resize(n);
	pressure_vector.resize(n);
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
	// set zero for mesh values.
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
				if (sys->p.sflow_boundary_conditions == 0) {
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
					if (m != 2) {
						/* ux data is left-middle of the cell labbeled as (i, j).
						 * Thus, the particle date is distributed
						 * from i-1 to i+2 along x      --> l = -1...2
						 * but from j-1 to j+1 along z. --> m = -1...1
						 */
						double zo_z = zo-pos[k].z;
						double xu = iix*dx;
						double xo_xu = xo-xu;
						double w_udx = weightFunc(xo_xu*xo_xu + zo_z*zo_z);
						total_weight_udx += w_udx;
						udx_values.push_back(ud_x*w_udx);
						phi_ux_values.push_back(w_udx);
					} else {
						/* mesh_nb includes this point.
						 * Therefore, 0 need to be given.
						 */
						udx_values.push_back(0);
						phi_ux_values.push_back(0);
					}
					if (l != 2) {
						/* uz data is bottom-center of the cell labbeled as (i, j).
						 * Thus, the particle date is distributed
						 * from i-1 to i+1 along x      --> l = -1...1
						 * but from j-1 to j+2 along z. --> m = -1...2
						 */
						double xo_x = xo-pos[k].x;
						double zu = iiz*dz;
						double zo_zu = zo-zu;
						double w_udz = weightFunc(xo_x*xo_x + zo_zu*zo_zu);
						total_weight_udz += w_udz;
						udz_values.push_back(ud_z*w_udz);
						phi_uz_values.push_back(w_udz);
					} else {
						/* mesh_nb includes this point.
						 * Therefore, 0 need to be given.
						 */
						udz_values.push_back(0);
						phi_uz_values.push_back(0);
					}
				}
			}
		}
		/* The distributed data are normalized by the total weights.
		 * The normalized data are set to the mesh points.
		 */
		for (int m=0; m< mesh_nb.size(); m++) {
			int k = mesh_nb[m];
			if (k < n) {
				u_diff_x[k] += udx_values[m]/(total_weight_udx);
				phi_ux[k] += (particle_volume*phi_ux_values[m]/cell_area)/total_weight_udx;
			}
			u_diff_z[k] += udz_values[m]/(total_weight_udz);
			phi_uz[k] += (particle_volume*phi_uz_values[m]/cell_area)/total_weight_udz;
		}
	}
	for (int k=0; k<n; k++) {
		u_diff_x[k] *= 1/conv_factor;
		u_diff_z[k] *= 1/conv_factor;
	}
}

void SolventFlow::pressureController()
{
	u_ave = calcAverageU();
	average_pressure_x.update(pressure_difference_x, sys->get_time()/sys->p.sflow_re);
	if (channel_flow) {
		//pressure_difference = 10*pressure_difference_;
		if (u_ave.x < sys->p.sflow_target_flux) {
			pressure_difference_x += sys->p.sflow_pcontrol_increment;
		} else {
			pressure_difference_x -= sys->p.sflow_pcontrol_increment;
		}
	} else {
		//double pd_increment = sys->p.sflow_pcontrol_increment;
		double diff_x = u_ave.x-target_flux;
		//if (u_sol_ave.x < target_flux) {
		pressure_difference_x += -diff_x*sys->p.sflow_pcontrol_increment*sys->dt;
		//		} else {
		//		pressure_difference_x -= diff_x*sys->p.sflow_pcontrol_increment*sys->dt;
		//}
		pressure_difference_x += - sys->p.sflow_pcontrol_damper*(pressure_difference_x - average_pressure_x.get())*sys->dt;
		// pressure_difference_x -= (pressure_difference_x-average_pressure.get());
	}
}

void SolventFlow::update(double pressure_difference_)
{
	//static std::ofstream fout_tmp("debug.dat");
	d_tau = sys->dt/sys->p.sflow_re;
	particleVelocityDiffToMesh();
	predictorStep();
	calcVelocityDivergence();
	solvePressure();
	correctorStep();
	for (int k=0; k<n; k++) {
		u_x[k] = conv_factor*(u_sol_x[k] + phi_ux[k]*u_diff_x[k]);
		u_z[k] = conv_factor*(u_sol_z[k] + phi_uz[k]*u_diff_z[k]);
	}
	calcVelocityGradients();
}

double SolventFlow::porousResistance(double area_fraction)
{
	/*
	 * 27/8=3.375
	 */
	//double porosity = (1-area_fraction)/(1-average_area_fraction);
	double porosity = (1-area_fraction);
	double poro_factor = sys->p.sflow_Darcy_coeff;
	for (int k=0; k<sys->p.sflow_Darcy_power; k++) {
		poro_factor *= porosity;
	}
	return 3.375*area_fraction/poro_factor;
}

void SolventFlow::predictorStep()
{
	double dx2 = dx*dx;
	double dz2 = dz*dz;
	double scalefactor = length_scale*length_scale;
	for (int j=0; j<nz; j++) {
		int jp1, jm1;
		if (sys->p.sflow_boundary_conditions == 0) {
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
			double res_coeff_ux = scalefactor*porousResistance(phi_ux[k]);
			double res_coeff_uz = scalefactor*porousResistance(phi_uz[k]);
			if (sys->p.sflow_boundary_conditions == 0) {
				/* periodic boundary condtions in x and z directions.
				 */
				double dd_ux = (u_sol_x[ir]-2*u_sol_x[k]+u_sol_x[il])/dx2+(u_sol_x[ju]-2*u_sol_x[k]+u_sol_x[jd])/dz2;
				double dd_uz = (u_sol_z[ir]-2*u_sol_z[k]+u_sol_z[il])/dx2+(u_sol_z[ju]-2*u_sol_z[k]+u_sol_z[jd])/dz2;
				double ux_duxdx = u_sol_x[k]*(u_sol_x[ir]-u_sol_x[il])/(2*dx);
				double uz_duxdz = 0.25*(u_sol_z[k]+u_sol_z[il]+u_sol_z[ju]+u_sol_z[im1+jp1*nx])*(u_sol_x[ju]-u_sol_x[jd])/(2*dz);
				double uz_duzdz = u_sol_z[k]*(u_sol_z[ju]-u_sol_z[jd])/(2*dz);
				double uz_duzdx = 0.25*(u_sol_x[k]+u_sol_x[ir]+u_sol_x[jd]+u_sol_x[ip1+jm1*nx])*(u_sol_z[ir]-u_sol_z[il])/(2*dx);
				double fx = dd_ux - ux_duxdx - uz_duxdz + res_coeff_ux*u_diff_x[k];
				double fz = dd_uz - uz_duzdz - uz_duzdx + res_coeff_uz*u_diff_z[k];
				u_sol_ast_x[k] = u_sol_x[k] + d_tau*fx;
				u_sol_ast_z[k] = u_sol_z[k] + d_tau*fz;
				// - sys->p.sf_zfriction*u_sol_ave.z);
				/* sf_zfriction*u_sol_z[k] term is added
				 * to stabilize view center along z direction.
				 * When periodic boundary conditions are used, nothing bounds the z velocity.
				 * We may control by pressure difference along z-direction.
				 * But, the expected average value is zero and the overall velocity along z direction
				 * is not important. Thus, we introduce the friction to v = 0 frame.
				 */
			} else {
				/* periodic boundary condtions in x directions.
				 * Non-slip boundy conditions at z = 0 and z=Lz.
				 */
				if (j == 0) {
					// bottom
					double ux_bot_m1 = 2*ux_bot - u_sol_x[k]; // (ux(i,-1)+ux(i,0))/2 = ux_bot
					double dd_ux = (u_sol_x[ir]-2*u_sol_x[k]+u_sol_x[il])/dx2+(u_sol_x[ju]-2*u_sol_x[k]+ux_bot_m1)/dz2;
					double fx = dd_ux + res_coeff_ux*u_diff_x[k];
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
					double fx = dd_ux + res_coeff_ux*u_diff_x[k];
					double fz = dd_uz + res_coeff_uz*u_diff_z[k];
					u_sol_ast_x[k] = u_sol_x[k] + d_tau*fx;
					u_sol_ast_z[k] = u_sol_z[k] + d_tau*fz;
				}
			}
		}
	}
	if (sys->p.sflow_boundary_conditions == 1) {
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
	if (sys->p.sflow_boundary_conditions == 0) {
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
	// Grad phi . (U-u_s)
	if (sys->p.sflow_boundary_conditions == 0) {
		for (int j=0; j<nz; j++) {
			int jp1 = (j == nz-1 ? 0 : j+1);
			for (int i=0; i<nx; i++) {
				int ip1 = (i == nx-1 ? 0 : i+1);
				int k = i + j*nx;
				int k_r = ip1+j*nx; // right
				int k_u = i+jp1*nx; // up
				double dphix = (phi_ux[k_r] - phi_ux[k])/dx;
				double dphiz = (phi_uz[k_u] - phi_uz[k])/dz;
				double u_diff_x_c = (u_diff_x[k_r]+u_diff_x[k])/2;
				double u_diff_z_c = (u_diff_z[k_u]+u_diff_z[k])/2;
				double d_udiffx_x = (u_diff_x[k_r]-u_diff_x[k])/dx;
				double d_udiffz_z = (u_diff_z[k_u]-u_diff_z[k])/dz;
				double phi_c = (phi_ux[k_r] + phi_ux[k] + phi_uz[k] + phi_uz[k])/4;
				gr_phi_Ud_phi_div_Ud[k] = dphix*u_diff_x_c + dphiz*u_diff_z_c + phi_c*(d_udiffx_x+d_udiffz_z);
				//gr_phi_Ud_phi_div_Ud[k] = 0;
			}
		}
	} else {
		// @@@
	}
}

void SolventFlow::solvePressure()
{
	for (int k=0; k<n; k++) {
		rhs_vector(k) = (div_u_sol_ast[k]+gr_phi_Ud_phi_div_Ud[k])/(six_pi*d_tau);
	}
	if (pressure_difference_x != 0) {
		double delta_P_dxdx = pressure_difference_x/(dx*dx);
		for (int j=0; j<nz; j++) {
			int j_nx = j*nx;
			// i = 0 --> k = j*nx
			rhs_vector(j_nx) += -delta_P_dxdx;
			// i = nx-1 ---> k = nx-1+j*nx
			rhs_vector(nx-1+j_nx) += delta_P_dxdx;
		}
	}
	pressure_vector = psolver->solve(rhs_vector);
	if (psolver->info() != Eigen::Success) {
		// solving failed
		std::cerr << "solving failed" << std::endl;
		return;
	}
	double mean_p = pressure_vector.mean();
	for (int i=0; i<n; i++) {
		pressure[i] = pressure_vector(i)-mean_p;
	}
}

void SolventFlow::correctorStep()
{
	if (sys->p.sflow_boundary_conditions == 0) {
		/* periodic boundary condtions in x and z directions.
		 */
		double sixpi_d_tau_dx = six_pi*d_tau/dx;
		double sixpi_d_tau_dz = six_pi*d_tau/dz;
		for (int i=0; i<nx; i++) {
			int im1 = (i == 0 ? nx-1 : i-1);
			double pd_x = (i == 0 ? pressure_difference_x : 0);
			for (int j=0; j<nz; j++) {
				int k = i+nx*j;
				int jm1 = (j == 0 ? nz-1 : j-1);
				u_sol_x[k] = u_sol_ast_x[k] - sixpi_d_tau_dx*(pressure[k] - (pressure[im1+j*nx]+pd_x));
				u_sol_z[k] = u_sol_ast_z[k] - sixpi_d_tau_dz*(pressure[k] - pressure[i+jm1*nx]);
			}
		}
	} else {
		/* periodic boundary condtions in x directions.
		 * Non-slip boundy conditions at z = 0 and z=Lz.
		 * (uz[i, j = nz] = 0,)
		 */
		for (int i=0; i<nx; i++) {
			int im1 = (i == 0 ? nx-1 : i-1);
			double pd_x = (i == 0 ? pressure_difference_x : 0);
			{// bottom
				u_sol_x[i] = u_sol_ast_x[i] - d_tau*(pressure[i] - (pressure[im1]+pd_x))/dx;
				u_sol_z[i] = 0;
			}
			for (int j=1; j<nz; j++) {
				int k = i+nx*j;
				u_sol_x[k] = u_sol_ast_x[k] - six_pi*d_tau*(pressure[k] - (pressure[im1+j*nx]+pd_x))/dx;
				u_sol_z[k] = u_sol_ast_z[k] - six_pi*d_tau*(pressure[k] - pressure[i+(j-1)*nx])/dz;
			}
			{// top
				int j=nz;
				u_sol_z[i+j*nx] = 0;
			}
		}
	}
}

void SolventFlow::calcVelocityGradients()
{
	if (sys->p.sflow_boundary_conditions == 0) {
		for (int j=0; j<nz; j++) {
			int jm1 = (j == 0 ? nz-1 : j-1);
			int jp1 = (j == nz-1 ? 0 : j+1);
			for (int i=0; i<nx; i++) {
				int k = i+nx*j;
				int im1 = (i == 0 ? nx-1 : i-1);
				int ip1 = (i == nx-1 ? 0 : i+1);
				strain_rate_xx[k] = (u_x[ip1+j*nx]-u_x[k]       )/dx; // center d_ux_d_x
				strain_rate_zz[k] = (u_z[i+jp1*nx]-u_z[k]       )/dz; // center d_uz_d_z
				double d_uz_d_x   = (u_z[k]       -u_z[im1+j*nx])/dx; // left_bottom
				double d_ux_d_z   = (u_x[k]       -u_x[i+jm1*nx])/dz; // left_bottom
				omega[k]          = (d_ux_d_z - d_uz_d_x)/2; // left_bottom
				strain_rate_xz[k] = (d_ux_d_z + d_uz_d_x)/2; // left_bottom
			}
		}
	} else {
		for (int j=0; j<=nz; j++) {
			int jp1 = j+1;
			int jm1 = j-1;
			for (int i=0; i<nx; i++) {
				int k = i+nx*j;
				int ip1 = (i == nx-1 ? 0 : i+1);
				double d_ux_d_z, d_uz_d_x, d_ux_d_x, d_uz_d_z;
				if (j == 0) {
					// (u[j] + u[j-1])/2 = ux_bot
					// u[j-1] = 2*ux_bot - u[j];
					// u[j] - u[j-1] = u[j] -(2*ux_bot - u[j]) = 2*u[j] -2*ux_bot
					d_ux_d_z = 2*(u_sol_x[k]-ux_bot)/dz; // left_bottom
					d_uz_d_x = 0; // left_bottom (uz[j=0]= 0)
					d_ux_d_x = (u_sol_x[ip1]-u_sol_x[k])/dx; //center
					d_uz_d_z = (u_sol_z[i+nx]-u_sol_z[k])/dz; //center
				} else if (j == nz) {
					// (u[j]+u[j-1])/2 = ux_top
					// u[j] = 2ux_top - u[j-1]
					// u[j]-u[j-1] = 2ux_top - u[j-1] - u[j-1]
					int im1 = (i == 0 ? nx-1 : i-1);
					d_ux_d_x = (-u_sol_x[i+jm1*nx]+u_sol_x[im1+jm1*nx])/dx;
					d_uz_d_z = 0;
					d_ux_d_z = 2*(ux_top-u_sol_x[i+jm1*nx])/dz; // left_bottom
					d_uz_d_x = 0; // left_bottom
				} else {
					int k =i+nx*j;
					int im1 = (i == 0 ? nx-1 : i-1);
					int ip1 = (i == nx-1 ? 0 : i+1);
					d_ux_d_x = (u_sol_x[ip1+j*nx]-u_sol_x[k])/dx;
					d_uz_d_z = (u_sol_z[i+jp1*nx]-u_sol_z[k])/dz;
					d_ux_d_z = (u_sol_x[k]-u_sol_x[i+jm1*nx])/dz;
					d_uz_d_x = (u_sol_z[k]-u_sol_z[im1+nx*j])/dx;
				}
				omega[k] = -(d_uz_d_x-d_ux_d_z)/2;
				strain_rate_xx[k] = d_ux_d_x; // center
				strain_rate_xz[k] = (d_ux_d_z + d_uz_d_x)/2; // left_bottom
				strain_rate_zz[k] = d_uz_d_z; // center
			}
		}
	}
}

void SolventFlow::localFlow(const vec3d &p,
							vec3d &u_local,
							vec3d &omega_local,
							std::vector<double> &e_local)
{
	/* ux-data  (xo,      zo+dz/2)
	 * uz-data  (xo+dx/2, zo)
	 * omega-data  (xo, zo)
	 */
	double x_dx = p.x/dx;
	double z_dz = p.z/dz;
	int i = (int)x_dx;
	int j = (int)z_dz;
	int ip1 = (i == nx-1 ? 0 : i+1);
	int jp1 = (j == nz-1 ? 0 : j+1);
	int jnx = j*nx;
	int jp1nx = jp1*nx;
	/* x_diff = (p.x - i*dx)/dx
	 *        = p.x/dx - i
	 *        = x_dx - i
	 */
	double x_diff = x_dx - i;
	double z_diff = z_dz - j;
	int kk[4];
	/********************************************************
	 * ux interpolation
	 * (0, dy/2)
	 ********************************************************/
	double z_diffS;
	if (z_diff < 0.5) {
		int jm1 = (j == 0 ? nz-1 : j-1);
		z_diffS = z_diff + 0.5;
		kk[0] = i   + jm1*nx;
		kk[1] = ip1 + jm1*nx;
		kk[2] = i   + jnx;
		kk[3] = ip1 + jnx;
	} else {
		z_diffS = z_diff - 0.5;
		kk[0] = i   + jnx;
		kk[1] = ip1 + jnx;
		kk[2] = i   + jp1nx;
		kk[3] = ip1 + jp1nx;
	}
	u_local.x = u_x[kk[0]];
	u_local.x += (u_x[kk[1]]-u_x[kk[0]])*x_diff;
	u_local.x += (u_x[kk[2]]-u_x[kk[0]])*z_diffS;
	u_local.x += (u_x[kk[3]]-u_x[kk[2]]-u_x[kk[1]]+u_x[kk[0]])*x_diff*z_diffS;
	/********************************************************
	 * uz interpolation
	 * (dx/2, 0)
	 ********************************************************/
	double x_diffS;
	if (x_diff < 0.5) {
		int im1 = (i == 0 ? nx-1 : i-1);
		x_diffS = x_diff + 0.5;
		kk[0] = im1 + jnx;
		kk[1] = i   + jnx;
		kk[2] = im1 + jp1nx;
		kk[3] = i   + jp1nx;
	} else {
		x_diffS = x_diff - 0.5;
		kk[0] = i   + jnx;
		kk[1] = ip1 + jnx;
		kk[2] = i   + jp1nx;
		kk[3] = ip1 + jp1nx;
	}
	u_local.z = u_z[kk[0]];
	u_local.z += (u_z[kk[1]]-u_z[kk[0]])*x_diffS;
	u_local.z += (u_z[kk[2]]-u_z[kk[0]])*z_diff;
	u_local.z += (u_z[kk[3]]-u_z[kk[2]]-u_z[kk[1]]+u_z[kk[0]])*x_diffS*z_diff;
	/********************************************************
	 * omega interpolation
	 * (0, 0)
	 ********************************************************/
	kk[0] = i   + jnx;
	kk[1] = ip1 + jnx; // +dx
	kk[2] = i   + jp1nx; // +dz
	kk[3] = ip1 + jp1nx; // +dx+dz
	omega_local.y = omega[kk[0]];
	omega_local.y += (omega[kk[1]]-omega[kk[0]])*x_diff;
	omega_local.y += (omega[kk[2]]-omega[kk[0]])*z_diff;
	omega_local.y += (omega[kk[3]]-omega[kk[2]]-omega[kk[1]]+omega[kk[0]])*x_diff*z_diff;
	e_local[1] = strain_rate_xz[kk[0]];
	e_local[1] += (strain_rate_xz[kk[1]]-strain_rate_xz[kk[0]])*x_diff;
	e_local[1] += (strain_rate_xz[kk[2]]-strain_rate_xz[kk[0]])*z_diff;
	e_local[1] += (strain_rate_xz[kk[3]]-strain_rate_xz[kk[2]]-strain_rate_xz[kk[1]]+strain_rate_xz[kk[0]])*x_diff*z_diff;
	/***************************************************************/
	// (dx/2, dy/2)
	double xs_dx = x_dx+0.5; // (p.x+0.5*dx)/dx
	double zs_dz = z_dz+0.5; // (p.z+0.5*dz)/dz
	i = (int)xs_dx;
	j = (int)zs_dz;
	x_diff = xs_dx - i;
	z_diff = zs_dz - j;
	i--; // -0.5 = +0.5 - 1
	j--; // -0.5 = +0.5 - 1
	if (i == -1) {i += nx;}
	if (j == -1) {j += nz;}
	ip1 = (i == nx-1 ? 0 : i+1);
	jp1 = (j == nz-1 ? 0 : j+1);
	kk[0] = i   + jnx;
	kk[1] = ip1 + jnx; // +dx
	kk[2] = i   + jp1nx; // +dz
	kk[3] = ip1 + jp1nx; // +dx+dz
	e_local[0] = strain_rate_xx[kk[0]];
	e_local[0] += (strain_rate_xx[kk[1]]-strain_rate_xx[kk[0]])*x_diff;
	e_local[0] += (strain_rate_xx[kk[2]]-strain_rate_xx[kk[0]])*z_diff;
	e_local[0] += (strain_rate_xx[kk[3]]-strain_rate_xx[kk[2]]-strain_rate_xx[kk[1]]+strain_rate_xx[kk[0]])*x_diff*z_diff;

	e_local[2] = strain_rate_zz[kk[0]];
	e_local[2] += (strain_rate_zz[kk[1]]-strain_rate_zz[kk[0]])*x_diff;
	e_local[2] += (strain_rate_zz[kk[2]]-strain_rate_zz[kk[0]])*z_diff;
	e_local[2] += (strain_rate_zz[kk[3]]-strain_rate_zz[kk[2]]-strain_rate_zz[kk[1]]+strain_rate_zz[kk[0]])*x_diff*z_diff;
}

double SolventFlow::flowFiledDissipation()
{
	double energy_dissipation = 0;
	for (int j=0; j<nz; j++) {
		for (int i=0; i<nx; i++) {
			int k = i+j*nx;
			energy_dissipation += strain_rate_xx[k]*strain_rate_xx[k];
			energy_dissipation += 2*strain_rate_xz[k]*strain_rate_xz[k];
			energy_dissipation += strain_rate_zz[k]*strain_rate_zz[k];
		}
	}
	return energy_dissipation/sys->get_lx()/sys->get_lz();
}

double SolventFlow::particleDissipation()
{
	double energy_dissipation = 0;
	vec3d u, omega;
	for (int i=0; i<sys->np_mobile; i++) {
		energy_dissipation += 0.5*(sys->na_velocity[i].sq_norm()+ sys->na_ang_velocity[i].sq_norm());
	}
	return energy_dissipation/sys->get_lx()/sys->get_lz();
}

double SolventFlow::meanVelocity()
{
	double mean_velocity = 0;
	for (int k=0; k<n; k++) {
		mean_velocity += u_sol_x[k] ;
	}
	return mean_velocity/n;
}

vec3d SolventFlow::calcAverageU()
{
	int i = std::lround(0.5*nx);
	double u_x_total = 0;
	for (int j=0; j<nz; j++) {
		int k = i + nx*j;
		u_x_total += u_x[k];
	}
	int j = std::lround(0.5*nz);
	double u_z_total = 0;
	for (int i=0; i<nx; i++) {
		int k = i + nx*j;
		u_z_total += u_z[k];
	}
	return vec3d(u_x_total/nz, 0 , u_z_total/nx);
}

void SolventFlow::velocityProfile(std::ofstream &fout_fp)
{
	Eigen::VectorXd ux(nz);
	Eigen::VectorXd ux_na(nz);
	Eigen::VectorXd phi(nz);

	for (int j=0; j<nz; j++) {
		double total_ux = 0;
		double total_ux_na = 0;
		double total_phi = 0;
		for (int i=0; i<nx; i++){
			int k = i+nx*j;
			total_ux    += u_x[k];
			total_ux_na += conv_factor*u_diff_x[k];
			total_phi += phi_ux[k];
		}
		ux(j) = total_ux/nx;
		ux_na(j) = total_ux_na/nx;
		// na = up - us
		phi(j) = total_phi/nx;
		// na + us = up - us + us = up
	}
	double mean_ux = ux.mean();
	/*
	 * F_bf = up - u
	 * ux_na is the distributed up-u on mesh points.
	 * The total value (<ux_na>*[#mesh points] must be N*F_bf/
	 * std::cerr << "<up-u> = " << ux_na.mean()*nx*nz/sys->np_mobile << std::endl;
	 *
	 */
	std::cerr << "<up-u> = " << ux_na.mean()*nx*nz/sys->np_mobile << std::endl;
	for (int j=0; j<nz; j++) {
		double ux_out = ux(j) - mean_ux;
		fout_fp << pos[meshNb(0,j)].z << ' '; // 1: z
		fout_fp << ux_out << ' '; // 2: total velocity profile
		fout_fp << ux_na(j) + ux_out << ' '; // 3: u_particle
		fout_fp << phi(j); // 4: volume fraction
		fout_fp << std::endl;
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
				double vx = vfactor*(u_sol_x[k] +u_sol_x[meshNb(i+1,j)])/2;
				double vz = vfactor*(u_sol_z[k] +u_sol_z[meshNb(i,j+1)])/2;
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
				double udx = vfactor*(u_diff_x[k] +u_diff_x[meshNb(i+1,j)])/2;
				double udz = vfactor*(u_diff_z[k] +u_diff_z[meshNb(i,j+1)])/2;
				fout_flow << "s ";
				fout_flow << po.x      << ' ' << -0.02 << ' ' << po.z << ' ';
				fout_flow << po.x + udx << ' ' << -0.02 << ' ' << po.z + udz << std::endl;
			}
		}
	}
	
	if (1) {
		fout_flow << "y 5" << std::endl;
//		double cell_area = dx*dz;
		fout_flow << "@ 10" << std::endl;
		for (int j=0; j<nz; j++){
			for (int i=0; i<nx; i++){
				vec3d po = pos[i + nx*j] - o_shift;
				// double r = sqrt(cell_area*phi[i + nx*j]/M_PI);
				int ip1 = (i == nx-1 ? 0 : i+1);
				int jp1 = (j == nz-1 ? 0 : j+1);
				double r0 = porousResistance(phi_ux[i + nx*j]);
				double r1 = porousResistance(phi_ux[ip1 + nx*j]);
				double r2 = porousResistance(phi_uz[i + nx*j]);
				double r3 = porousResistance(phi_uz[i + nx*jp1]);
//				double r = (r0 + r1 + r2 +r3)/4;
				double r = 	0.01*gr_phi_Ud_phi_div_Ud[i + nx*j];
				//double s = 100*div_u_sol_ast[i+nx*j];
				//if (r > 0 ) {
				//					fout_flow << "@ 4" << std::endl;
				//				} else {
				//					fout_flow << "@ 6" << std::endl;
				//				}
				
				fout_flow << "r " << 0.1*abs(r) << std::endl;
				fout_flow << "c ";
				fout_flow << po.x << ' ' << -0.02 << ' ' << po.z << std::endl;
			}
		}
	}
	
	if (0) {
		fout_flow << "y 5" << std::endl;
		//		double cell_area = dx*dz;
		fout_flow << "@ 5" << std::endl;
		for (int j=0; j<nz; j++){
			for (int i=0; i<nx; i++){
				int k = i + j*nx;
				vec3d po = pos[k] - o_shift;
				double r = strain_rate_xz[k];
				fout_flow << "r " << abs(r) << std::endl;
				fout_flow << "c ";
				fout_flow << po.x -dx/2<< ' ' << -0.02 << ' ' << po.z - dz/2 << std::endl;
			}
		}
	}
	
	if (0) {
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

