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
using namespace std;

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
//	Eigen::setNbThreads(4);
	Eigen::initParallel();
	sys = sys_;
	re_num = sys->p.sflow_ReNum_p;
	//	if (sys->twodimension) {
	//		length_scale = sqrt(sys->p.sflow_ReNum/sys->p.sflow_ReNum_p);
	//	} else {
	//		//length_scale = pow(sys->p.sflow_ReNum/sys->p.sflow_ReNum_p, 0.3333);
	//		exit(1);
	//	}
	//std::cerr << "Ro/a = " << length_scale << std::endl;
	/* flow unit --> particle dynamics unit
	 * u in PD unit = (u in SF unit) * (R0/a)^{d-1}
	 */
	//conv_factor = length_scale; // for d = 2 (square if d = 3)
	//alpha = re_num/length_scale;
	six_pi = 6*M_PI;
	if (simulation_type == "sedimentation") {
		sedimentation = true;
		pressure_grad_x = 0;
	} else if (simulation_type == "channel flow") {
		channel_flow = true;
		pressure_grad_x = sys->pressure_drop;
		pressure_difference_x = pressure_grad_x*sys->get_lx();
		std::cerr << "pressure_grad_x = " << pressure_grad_x << std::endl;
	} else if (simulation_type == "simple shear") {
		simple_shear = true;
		std::ostringstream error_str;
		error_str << "Solvent flow algoritm for simple shear is not implemented yet. \n";
		error_str << "Lees--Edwards boundary condtions is complicated.\n";
		throw std::runtime_error(error_str.str());
	} else {
		std::ostringstream error_str;
		error_str << "Incorrect simulation type\n";
		throw std::runtime_error(error_str.str());
	}
	if (sedimentation) {
		average_pressure_x.setRelaxationTime(sys->p.sflow_pcontrol_rtime);
		std::cerr << "sedimentation simulation" << std::endl;
	}

	nx = (int)(sys->get_lx()/sys->p.sflow_dx);
	dx = sys->get_lx()/nx;
	std::cerr << sys->p.sflow_dx << " --> " << dx << std::endl;
	if (sys->get_lx() == sys->get_lz()) {
		dz = dx;
		nz = nx;
	} else {
		std::cerr << std::setprecision(16) << (sys->get_lz()/dx) << std::endl;
		nz = std::lround(sys->get_lz()/dx);
		std::cerr << "nz = " << nz << std::endl;
		dz = sys->get_lz()/nz;
	}
	n = nx*nz;
	if (sys->p.sflow_boundary_conditions == 0) {
		nz_adpt = nz;
	} else {
		nz_adpt = nz+1;
		ux_top = 0;
		ux_bot = 0;
	}
	n_adpt = nx*nz_adpt;
//	if (abs(dx-dz) > 1e-8) {
//		std::ostringstream error_str;
//		error_str << "dx = " << dx << "  dz = " << dz << "\n";
//		error_str << "dx != dz. Modify sflow_nx or sflow_nz\n";
//		throw std::runtime_error(error_str.str());
//	}
	cell_area = dx*dz;
	system_volume = sys->get_lx()*sys->get_lz()*2;
	smooth_length = sys->p.sflow_smooth_length;
	sq_smooth_length = smooth_length*smooth_length;
	pressure.resize(n, 0);
	Urel_x.resize(n, 0);
	Urel_z.resize(n_adpt, 0);
	u_x.resize(n, 0);
	u_z.resize(n_adpt, 0);
	u_sol_x.resize(n, 0);
	u_sol_z.resize(n_adpt, 0);
	u_sol_ast_x.resize(n, 0);
	u_sol_ast_z.resize(n_adpt, 0);
	div_u_sol_ast.resize(n, 0);
	div_u.resize(n, 0);
	phi_ux.resize(n,0);
	phi_uz.resize(n_adpt, 0);
	Urel_phi_x.resize(n, 0); // Urel * phi / (1-phi)
	Urel_phi_z.resize(n_adpt, 0); // Urel * phi / (1-phi)
	pos.resize(n_adpt);
	omega.resize(n_adpt, 0);
	strain_rate_xx.resize(n, 0);
	strain_rate_xz.resize(n_adpt, 0);
	strain_rate_zz.resize(n, 0);
	
	for (int j=0; j<nz_adpt; j++) {
		for (int i=0; i<nx; i++) {
			pos[i+nx*j].set(dx*i+dx/2, 0, dz*j+dz/2);
		}
	}
	initPoissonSolver();
	double particle_area = 0;
	for (int i=0; i<sys->np_mobile; i++) {
		particle_area += M_PI*sys->conf->radius[i]*sys->conf->radius[i];
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
			lmat[k][meshNb(i+1, j)] += dx2i;  //p(xi+1,zi) --> b(xi,zi)
			lmat[k][meshNb(i-1, j)] += dx2i;  //p(xi-1,zi) --> b(xi,zi)
			if (sys->p.sflow_boundary_conditions == 0) {
				lmat[k][k]              += -2*(dx2i+dz2i);  //p(xi,zi) --> b(xi,zi)
				lmat[k][meshNb(i, j+1)] += dz2i;  //p(xi,zi+1) --> b(xi,zi)
				lmat[k][meshNb(i, j-1)] += dz2i;  //p(xi,zi-1) --> b(xi,zi
			} else {
				if (j == 0) {
					// dp/dz = 0 at the bottom wall --> p(i, j-1) = p(i, j)
					lmat[k][k]              += -(2*dx2i+dz2i);  //p(xi,zi) --> b(xi,zi)
					lmat[k][meshNb(i, j+1)] += dz2i;  //p(xi,zi+1) --> b(xi,zi)
				} else if (j == nz-1) {
					// dp/dz = 0 at the top wall --> p(i, j+1) = p(i, j)
					lmat[k][k]              += -(2*dx2i+dz2i);  //p(xi,zi) --> b(xi,zi)
					lmat[k][meshNb(i, j-1)] += dz2i;  //j = nz-1, k = i+(nz-1)*nx =
				} else {
					lmat[k][k]              += -2*(dx2i+dz2i);  //p(xi,zi) --> b(xi,zi)
					lmat[k][meshNb(i, j+1)] += dz2i;  //p(xi,zi+1) --> b(xi,zi)
					lmat[k][meshNb(i, j-1)] += dz2i;  //p
				}
			}
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
	//	std::cerr << t_lmat.size() << std::endl;
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
	std::fill(Urel_x.begin(), Urel_x.end(), 0);
	std::fill(Urel_z.begin(), Urel_z.end(), 0);
	std::fill(phi_ux.begin(), phi_ux.end(), 0);
	std::fill(phi_uz.begin(), phi_uz.end(), 0);
	for (int i=0; i < sys->get_np(); i++) {
		double x = sys->conf->position[i].x;
		double z = sys->conf->position[i].z;
		double Uus_x = sys->na_velocity.vel[i].x; // This is U - u (note: not U - u_s)
		double Uus_z = sys->na_velocity.vel[i].z; // This is U - u
		double radius = sys->conf->radius[i];
		double particle_volume = M_PI*radius*radius;
//		if (i >= sys->np_mobile) {
//			particle_volume *= 1;
//		}
		int ix = (int)(x/dx);
		int iz = (int)(z/dz);
		mesh_nb_x.clear();
		mesh_nb_z.clear();
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
				if (iiz >= 0) {
					int k = iix + iiz*nx;
					if (iiz < nz && m != 2) {
						mesh_nb_x.push_back(k);
						/* ux data is left-middle of the cell labbeled as (i, j).
						 * Thus, the particle date is distributed
						 * from i-1 to i+2 along x      --> l = -1...2
						 * but from j-1 to j+1 along z. --> m = -1...1
						 */
						double zo_z = zo-pos[k].z; // pos_z = j*dz + dz/2
						double xu = iix*dx;
						double xo_xu = xo-xu;
						double w_udx = weightFunc(xo_xu*xo_xu + zo_z*zo_z);
						total_weight_udx += w_udx;
						udx_values.push_back(Uus_x*w_udx);
						phi_ux_values.push_back(w_udx);
					}
					if (iiz < nz_adpt && l != 2) {
						mesh_nb_z.push_back(k);
						/* uz data is bottom-center of the cell labbeled as (i, j).
						 * Thus, the particle date is distributed
						 * from i-1 to i+1 along x      --> l = -1...1
						 * but from j-1 to j+2 along z. --> m = -1...2
						 */
						double xo_x = xo-pos[k].x; // pos_x = i*dx + dd/2
						double zu = iiz*dz;
						double zo_zu = zo-zu;
						double w_udz = weightFunc(xo_x*xo_x + zo_zu*zo_zu);
						total_weight_udz += w_udz;
						udz_values.push_back(Uus_z*w_udz);
						phi_uz_values.push_back(w_udz);
					}
				}
			}
		}
		/* The distributed data are normalized by the total weights.
		 * The normalized data are set to the mesh points.
		 */
		for (int m=0; m< mesh_nb_x.size(); m++) {
			int k = mesh_nb_x[m];
			Urel_x[k] += udx_values[m]/total_weight_udx;
			phi_ux[k] += (particle_volume*phi_ux_values[m]/cell_area)/total_weight_udx;
		}
		for (int m=0; m< mesh_nb_z.size(); m++) {
			int k = mesh_nb_z[m];
			Urel_z[k] += udz_values[m]/total_weight_udz;
			phi_uz[k] += (particle_volume*phi_uz_values[m]/cell_area)/total_weight_udz;
		}
	}
	// Here, from U-u to U - us by the factor 1/(1-phi)
	for (int k=0; k<n; k++) {
		Urel_phi_x[k] = Urel_x[k]*phi_ux[k]/(1-phi_ux[k]);
	}
	for (int k=0; k<n_adpt; k++) {
		Urel_phi_z[k] = Urel_z[k]*phi_uz[k]/(1-phi_uz[k]);
	}
	/**** Doctor codes ****/
	//doctor_phi();
}

void SolventFlow::doctor_phi()
{
	/* Averages of local volume fractions are equal the global value */
	double phi_x_total = 0;
	double phi_z_total = 0;
	for (int k=0; k<n; k++) {
		phi_x_total += phi_ux[k];
		phi_z_total += phi_uz[k];
	}
	std::cerr << phi_x_total /n << ' ' << phi_z_total /n << std::endl;
}

void SolventFlow::pressureController()
{
	u_ave = calcAverageU(); // fluid unit
	average_pressure_x.update(pressure_grad_x, sys->get_time());
	if (channel_flow) {
		//
	} else if (sedimentation) {
		double diff_x = u_ave.x-target_flux;
		pressure_grad_x += -diff_x*sys->p.sflow_pcontrol_increment*sys->dt;
		pressure_grad_x += -sys->p.sflow_pcontrol_damper*(pressure_grad_x - average_pressure_x.get())*sys->dt;
	}
}

void SolventFlow::update()
{
	//static std::ofstream fout_tmp("debug.dat");
	particleVelocityDiffToMesh();
	predictorStep();
	calcVelocityDivergence();
	solvePressure();
	correctorStep();
	for (int k=0; k<n; k++) {
		/* u = suspension total velocity */
		u_x[k] = u_sol_x[k] + Urel_phi_x[k];// particle unit
		u_z[k] = u_sol_z[k] + Urel_phi_z[k];// particle unit
	}
	if (sys->p.sflow_boundary_conditions == 1) {
		for (int k=n; k<n_adpt; k++) {
			u_z[k] = u_sol_z[k] + Urel_phi_z[k];// particle unit
		}
	}
	if (0) {
		calcSuspensionVelocityDivergence();
	}
	//	std::cerr << "p = " << pressure[0] << ' ' << pos[0].x << ' ' <<pos[0].z << std::endl;
	calcVelocityGradients();
	return;
}

double SolventFlow::porousResistance(double phi)
{
	/*
	 * 27/8=3.375
	 */
	//double porosity = (1-area_fraction)/(1-average_area_fraction);
	double porosity = 1-phi;
	double poro_factor = 1;
	for (int k=0; k<sys->p.sflow_Darcy_power; k++) {
		poro_factor *= porosity;
	}
	return 3.375*phi*sys->p.sflow_Darcy_coeff/poro_factor;
}

void SolventFlow::predictorStep()
{
	double dx2 = dx*dx;
	double dz2 = dz*dz;
	double dd_ux, dd_uz, ux_duxdx, uz_duxdz, uz_duzdz, ux_duzdx;
	for (int j=0; j<nz_adpt; j++) {
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
			double res_coeff_ux = porousResistance(phi_ux[k]);
			double res_coeff_uz = porousResistance(phi_uz[k]);
			if (sys->p.sflow_boundary_conditions == 0) {
				/* periodic boundary condtions in x and z directions.
				 */
				dd_ux = (u_sol_x[ir]-2*u_sol_x[k]+u_sol_x[il])/dx2+(u_sol_x[ju]-2*u_sol_x[k]+u_sol_x[jd])/dz2;
				dd_uz = (u_sol_z[ir]-2*u_sol_z[k]+u_sol_z[il])/dx2+(u_sol_z[ju]-2*u_sol_z[k]+u_sol_z[jd])/dz2;
				ux_duxdx = u_sol_x[k]*(u_sol_x[ir]-u_sol_x[il])/(2*dx);
				uz_duxdz = 0.25*(u_sol_z[k]+u_sol_z[il]+u_sol_z[ju]+u_sol_z[im1+jp1*nx])*(u_sol_x[ju]-u_sol_x[jd])/(2*dz);
				uz_duzdz = u_sol_z[k]*(u_sol_z[ju]-u_sol_z[jd])/(2*dz);
				ux_duzdx = 0.25*(u_sol_x[k]+u_sol_x[ir]+u_sol_x[jd]+u_sol_x[ip1+jm1*nx])*(u_sol_z[ir]-u_sol_z[il])/(2*dx);
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
				double u_sol_x_k;
				double u_sol_x_ir;
				double u_sol_x_il;
				double u_sol_x_jd;
				double u_sol_z_jd;
				double u_sol_x_ju;
				double u_sol_z_ju;
				double u_sol_z_ilju;
				double u_sol_x_irjd;
				if (j == 0) {
					u_sol_x_jd   = 2*ux_bot - u_sol_x[k]; // (ux(i,-1)+ux(i,0))/2 = ux_bot
					u_sol_x_irjd = 2*ux_bot - u_sol_x[ip1];
					u_sol_z_jd   = 0;
				} else {
					u_sol_x_jd   = u_sol_x[jd];
					u_sol_x_irjd = u_sol_x[ip1+jm1*nx];
					u_sol_z_jd   = u_sol_z[jd];
				}
				if (j == nz-1) {
					u_sol_x_ju = 2*ux_top - u_sol_x[k];
				} else {
					u_sol_x_ju = u_sol_x[ju];
				}
				if (j == nz) {
					u_sol_x_k  = 2*ux_top - u_sol_x[jd];
					u_sol_x_ir = 2*ux_top - u_sol_x[ip1 + jm1*nx];
					u_sol_x_il = 2*ux_top - u_sol_x[im1 + jm1*nx];
					u_sol_z_ju = 0;
					u_sol_z_ilju = 0;
 				} else {
					u_sol_x_k  = u_sol_x[k];
					u_sol_x_ir = u_sol_x[ir];
					u_sol_x_il = u_sol_x[il];
					u_sol_z_ju = u_sol_z[ju];
					u_sol_z_ilju = u_sol_z[im1+jp1*nx];
				}
				dd_ux = (u_sol_x_ir -2*u_sol_x_k +u_sol_x_il) /dx2+(u_sol_x_ju-2*u_sol_x_k +u_sol_x_jd)/dz2;
				dd_uz = (u_sol_z[ir]-2*u_sol_z[k]+u_sol_z[il])/dx2+(u_sol_z_ju-2*u_sol_z[k]+u_sol_z_jd)/dz2;
				ux_duxdx = u_sol_x_k*(u_sol_x_ir-u_sol_x_il)/(2*dx);
				uz_duxdz = 0.25*(u_sol_z[k]+u_sol_z[il]+u_sol_z_ju+u_sol_z_ilju)*(u_sol_x_ju-u_sol_x_jd)/(2*dz);

				ux_duzdx = 0.25*(u_sol_x_k+u_sol_x_ir+u_sol_x_jd+u_sol_x_irjd)*(u_sol_z[ir]-u_sol_z[il])/(2*dx);
				uz_duzdz = u_sol_z[k]*(u_sol_z_ju-u_sol_z_jd)/(2*dz);

			}
			//std::cerr << Urel_x[k] << ' ' << Urel_z[k] << std::endl;
			if (j != nz) {
				double fx = dd_ux + res_coeff_ux*Urel_x[k];
				u_sol_ast_x[k] = u_sol_x[k] + sys->dt*(fx/re_num - (ux_duxdx + uz_duxdz));
				//u_sol_ast_x[k] = u_sol_x[k] + sys->dt*fx/re_num;
			}
			double fz = dd_uz + res_coeff_uz*Urel_z[k];
			u_sol_ast_z[k] = u_sol_z[k] + sys->dt*(fz/re_num - (ux_duzdx + uz_duzdz));
			//u_sol_ast_z[k] = u_sol_z[k] + sys->dt*fz/re_num;
			
			//				u_sol_ast_x[k] = u_sol_x[k] + sys->dt*(fx/re_num);
			//				u_sol_ast_z[k] = u_sol_z[k] + sys->dt*(fz/re_num);
		}
	}

}

void SolventFlow::calcVelocityDivergence()
{
	/*
	 *
	 \nabla^2 p = \frac{Re }{6 \pi \Delta t}
	 \left( \nabla \cdot \boldsymbol{u}_s^{\ast}
	 + \frac{1}{(1-\phi)^2} \nabla \phi \cdot (U-u)
	 + \frac{\phi}{1-\phi} \nabla \cdot (U-u)  \right)
	 *
	 */
	
	// This works for both boundary conditions.
	// u_sol_ast_z[i, j = 0] = 0 (top and bottom)
	if (sys->p.sflow_boundary_conditions == 0) {
		for (int j=0; j<nz; j++) {
			int j_nx = j*nx;
			int jp1_nx = (j == nz-1 ? 0 : (j+1)*nx);
			for (int i=0; i<nx; i++) {
				int ip1 = (i == nx-1 ? 0 : i+1);
				int k   = i   + j_nx;
				int k_r = ip1 + j_nx; // right
				int k_u = i   + jp1_nx; // up
				double dux_dx = ((u_sol_ast_x[k_r]+Urel_phi_x[k_r]) - (u_sol_ast_x[k]+Urel_phi_x[k]))/dx;
				double duz_dz = ((u_sol_ast_z[k_u]+Urel_phi_z[k_u]) - (u_sol_ast_z[k]+Urel_phi_z[k]))/dz;
				div_u_sol_ast[k] = dux_dx + duz_dz; // center
			}
		}
	} else {
		for (int j=0; j<nz; j++) {
			int j_nx = j*nx;
			int jp1_nx = (j+1)*nx;
			for (int i=0; i<nx; i++) {
				int ip1 = (i != nx-1 ? i+1 : 0);
				int k   = i   + j_nx;
				int k_r = ip1 + j_nx; // right
				int k_u = i   + jp1_nx; // up
				double dux_dx = ((u_sol_ast_x[k_r]+Urel_phi_x[k_r]) - (u_sol_ast_x[k]+Urel_phi_x[k]))/dx;
				double duz_dz = ((u_sol_ast_z[k_u]+Urel_phi_z[k_u]) - (u_sol_ast_z[k]+Urel_phi_z[k]))/dz;
				//double dux_dx = (u_sol_ast_x[k_r] - u_sol_ast_x[k])/dx;
				//double duz_dz = (u_sol_ast_z[k_u] - u_sol_ast_z[k])/dz;
				div_u_sol_ast[k] = dux_dx + duz_dz; // center
			}
		}
	}
}

void SolventFlow::calcSuspensionVelocityDivergence()
{
	/*
	 *
	 \nabla^2 p = \frac{Re }{6 \pi \Delta t}
	 \left( \nabla \cdot \boldsymbol{u}_s^{\ast}
	 + \frac{1}{(1-\phi)^2} \nabla \phi \cdot (U-u)
	 + \frac{\phi}{1-\phi} \nabla \cdot (U-u)  \right)
	 *
	 */
	
	// This works for both boundary conditions.
	// u_sol_ast_z[i, j = 0] = 0 (top and bottom)
	if (sys->p.sflow_boundary_conditions == 0) {
		for (int j=0; j<nz; j++) {
			int j_nx = j*nx;
			int jp1_nx = (j != nz-1 ? (j+1)*nx : 0);
			for (int i=0; i<nx; i++) {
				int ip1 = (i == nx-1 ? 0 : i+1);
				int k   = i   + j_nx;
				int k_r = ip1 + j_nx; // right
				int k_u = i   + jp1_nx; // up
				double dux_dx = (u_x[k_r] - u_x[k])/dx;
				double duz_dz = (u_z[k_u] - u_z[k])/dz;
				div_u[k] = dux_dx + duz_dz; // center
			}
		}
	} else {
		for (int j=0; j<nz; j++) {
			int j_nx = j*nx;
			int jp1_nx = (j+1)*nx;
			for (int i=0; i<nx; i++) {
			int ip1 = (i == nx-1 ? 0 : i+1);
				int k   = i   + j_nx;
				int k_r = ip1 + j_nx; // right
				int k_u = i   + jp1_nx; // up
				double dux_dx = (u_x[k_r] - u_x[k])/dx;
				double duz_dz = (u_z[k_u] - u_z[k])/dz;
				div_u[k] = dux_dx + duz_dz; // center
			}
		}
	}
}

void SolventFlow::solvePressure()
{
	double rhs_coeff = re_num/(six_pi*sys->dt);
	for (int k=0; k<n; k++) {
		rhs_vector(k) = div_u_sol_ast[k]*rhs_coeff;
	}

	if (pressure_grad_x != 0) {
		double delta_P_dxdx = pressure_difference_x/(dx*dx);
		for (int j=0; j<nz; j++) {
			int j_nx = j*nx;
			// i = 0 --> k = j*nx
			rhs_vector(j_nx) += -delta_P_dxdx;
			// i = nx-1 ---> k = nx-1+j*nx
			rhs_vector(nx-1+j_nx) += delta_P_dxdx;
		}
	}
	double rhs_mean = rhs_vector.mean();
	for (int i=0; i<n; i++) {
		rhs_vector(i) -= rhs_mean;
	}
	pressure_vector = psolver->solve(rhs_vector);
	if (psolver->info() != Eigen::Success) {
		// solving failed
		std::cerr << "solving failed" << std::endl;
		return;
	}
	double mean_p = pressure_vector.mean();
	for (int i=0; i<n; i++) {
		pressure[i] = pressure_vector(i) - mean_p;
		//pressure[i] = pressure_vector(i);
	}
}

void SolventFlow::correctorStep()
{
	double sixpi_dt_Re = six_pi*sys->dt/re_num;
	double sixpi_dt_Re_dx = sixpi_dt_Re/dx;
	double sixpi_dt_Re_dz = sixpi_dt_Re/dz;
	if (sys->p.sflow_boundary_conditions == 0) {
		/* periodic boundary condtions in x and z directions.
		 */
		for (int i=0; i<nx; i++) {
			int im1 = (i == 0 ? nx-1 : i-1);
			double pd_x = (i == 0 ? pressure_difference_x : 0);
			for (int j=0; j<nz; j++) {
				int k = i+nx*j;
				int jm1 = (j == 0 ? nz-1 : j-1);
				u_sol_x[k] = u_sol_ast_x[k] - sixpi_dt_Re_dx*(pressure[k] - (pressure[im1+j*nx]+pd_x));
				u_sol_z[k] = u_sol_ast_z[k] - sixpi_dt_Re_dz*(pressure[k] - pressure[i+jm1*nx]);
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
			// bottom  j = 0
			u_sol_x[i] = u_sol_ast_x[i] - sixpi_dt_Re_dx*(pressure[i] - (pressure[im1]+pd_x));
			u_sol_z[i] = 0;

			for (int j=1; j<nz; j++) {
				int k = i+nx*j;
				u_sol_x[k] = u_sol_ast_x[k] - sixpi_dt_Re_dx*(pressure[k] - (pressure[im1+j*nx]+pd_x));
				u_sol_z[k] = u_sol_ast_z[k] - sixpi_dt_Re_dz*(pressure[k] - pressure[i+(j-1)*nx]);
			}
			// top
			// j = nz;
			u_sol_z[i+nx*nz] = 0;
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
				strain_rate_xx[k] = (u_x[ip1+j*nx]-u_x[k])/dx; // center d_ux_d_x; // sxx - (sxx+szz)/2 = sxx/2-szz/2
				strain_rate_zz[k] = (u_z[i+jp1*nx]-u_z[k])/dz; // szz - (sxx+szz)/2 = szz/2-dxx/2
				double d_uz_d_x   = (u_z[k]-u_z[im1+j*nx])/dx; // left_bottom
				double d_ux_d_z   = (u_x[k]-u_x[i+jm1*nx])/dz; // left_bottom
				omega[k]          = (d_ux_d_z - d_uz_d_x)/2; // left_bottom
				strain_rate_xz[k] = (d_ux_d_z + d_uz_d_x)/2; // left_bottom
			}
		}
	} else {
		for (int j=0; j<nz; j++) {
			int jp1_nx = (j+1)*nx;
			for (int i=0; i<nx; i++) {
				int k = i+nx*j;
				int ip1 = (i == nx-1 ? 0 : i+1);
				double d_ux_d_x = (u_x[ip1+j*nx]-u_x[k])/dx;
				double d_uz_d_z = (u_z[i+jp1_nx]-u_z[k])/dz;
				strain_rate_xx[k] = d_ux_d_x; // center
				strain_rate_zz[k] = d_uz_d_z; // center
			}
		}

		for (int j=0; j<nz_adpt; j++) {
			int jm1_nx = (j-1)*nx;
			for (int i=0; i<nx; i++) {
				int k = i+nx*j;
				double d_ux_d_z, d_uz_d_x;
				int im1 = (i == 0 ? nx-1 : i-1);
				if (j == 0) {
					d_ux_d_z = 2*(u_x[k]-ux_bot)/dz;
					d_uz_d_x = 0;
				} if (j == nz) {
					d_ux_d_z = 2*(ux_top-u_x[i+jm1_nx])/dz;
					d_uz_d_x = 0;
				} else {
					d_ux_d_z = (u_x[k]-u_x[i  +jm1_nx])/dz;
					d_uz_d_x = (u_z[k]-u_z[im1+j  *nx])/dx;
				}
				omega[k]          = (d_ux_d_z - d_uz_d_x)/2; // left_bottom
				strain_rate_xz[k] = (d_ux_d_z + d_uz_d_x)/2; // left_bottom

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
	int jp1;
	int jm1;
	if (sys->p.sflow_boundary_conditions == 0) {
		jp1 = (j == nz-1 ? 0 : j+1);
		jm1 = (j == 0 ? nz-1 : j-1);
	} else {
		jp1 = j+1;
		jm1 = j-1;
	}
	/* x_diff = (p.x - i*dx)/dx
	 *        = p.x/dx - i
	 *        = x_dx - i
	 */
	double x_diff = x_dx - i; // p.x/dx - xo/dx
	double z_diff = z_dz - j;
	int kk[4];
	/********************************************************
	 * ux interpolation
	 * (0, dy/2)
	 ********************************************************/
	double z_diffS;
	double ux0, ux1, ux2, ux3;
	if (z_diff < 0.5) {
		z_diffS = z_diff + 0.5;
		kk[2] = i   + j*nx;
		kk[3] = ip1 + j*nx;
		ux2 = u_x[kk[2]];
		ux3 = u_x[kk[3]];
		if (jm1 != -1) {
			kk[0] = i   + jm1*nx;
			kk[1] = ip1 + jm1*nx;
			ux0 = u_x[kk[0]];
			ux1 = u_x[kk[1]];
		} else {
			ux0 = 2*ux_bot - ux2; // (ux(i,-1)+ux(i,0))/2 = ux_bot
			ux1 = 2*ux_bot - ux3; // (ux(i,-1)+ux(i,0))/2 = ux_bot
		}
	} else {
		z_diffS = z_diff - 0.5;
		kk[0] = i   + j*nx;
		kk[1] = ip1 + j*nx;
		ux0 = u_x[kk[0]];
		ux1 = u_x[kk[1]];
		if (jp1 != nz) {
			kk[2] = i   + jp1*nx;
			kk[3] = ip1 + jp1*nx;
			ux2 = u_x[kk[2]];
			ux3 = u_x[kk[3]];
		} else {
			ux2 = 2*ux_top - ux0;
			ux3 = 2*ux_top - ux1;
		}
	}
	u_local.x = ux0 + (ux1-ux0)*x_diff + (ux2-ux0)*z_diffS + (ux3-ux2-ux1+ux0)*x_diff*z_diffS;

	/********************************************************
	 * uz interpolation
	 * (dx/2, 0)
	 ********************************************************/
	double x_diffS;
	if (x_diff < 0.5) {
		int im1 = (i == 0 ? nx-1 : i-1);
		x_diffS = x_diff + 0.5;
		kk[0] = im1 + j*nx;
		kk[1] = i   + j*nx;
		kk[2] = im1 + jp1*nx;
		kk[3] = i   + jp1*nx;
	} else {
		x_diffS = x_diff - 0.5;
		kk[0] = i   + j*nx;
		kk[1] = ip1 + j*nx;
		kk[2] = i   + jp1*nx;
		kk[3] = ip1 + jp1*nx;
	}
	u_local.z = u_z[kk[0]];
	u_local.z += (u_z[kk[1]]-u_z[kk[0]])*x_diffS;
	u_local.z += (u_z[kk[2]]-u_z[kk[0]])*z_diff;
	u_local.z += (u_z[kk[3]]-u_z[kk[2]]-u_z[kk[1]]+u_z[kk[0]])*x_diffS*z_diff;
	/********************************************************
	 * omega interpolation
	 * (0, 0)
	 ********************************************************/
	kk[0] = i   + j  *nx;
	kk[1] = ip1 + j  *nx; // +dx
	kk[2] = i   + jp1*nx; // +dz
	kk[3] = ip1 + jp1*nx; // +dx+dz
	omega_local.y = omega[kk[0]];
	omega_local.y += (omega[kk[1]]-omega[kk[0]])*x_diff;
	omega_local.y += (omega[kk[2]]-omega[kk[0]])*z_diff;
	omega_local.y += (omega[kk[3]]-omega[kk[2]]-omega[kk[1]]+omega[kk[0]])*x_diff*z_diff;
	e_local[1] =   strain_rate_xz[kk[0]];
	e_local[1] += (strain_rate_xz[kk[1]]-strain_rate_xz[kk[0]])*x_diff;
	e_local[1] += (strain_rate_xz[kk[2]]-strain_rate_xz[kk[0]])*z_diff;
	e_local[1] += (strain_rate_xz[kk[3]]-strain_rate_xz[kk[2]]-strain_rate_xz[kk[1]]+strain_rate_xz[kk[0]])*x_diff*z_diff;
	/***************************************************************/
	/* (dx/2, dy/2)
	 * Get integers i and j from a shifted origin (-dx/2, -dz/2)
	 * instead of (dx/2, dz/2).
	 * In this way, we can convert all points to positive integers.
	 * Then, we can shift the integer values (-1, -1) for right indexes.
	 */
	
	double xs_dx = x_dx+0.5; // (p.x+0.5*dx)/dx     x_dx = p.x/dx;
	double zs_dz = z_dz+0.5; // (p.z+0.5*dz)/dz
	i = (int)xs_dx;
	j = (int)zs_dz;
	x_diff = xs_dx - i;
	z_diff = zs_dz - j;
	if (i == 0) {
		i = nx-1;
	} else {
		i--; // -0.5 = +0.5 - 1
	}
	j--; // -0.5 = +0.5 - 1
	ip1 = (i == nx-1 ? 0 : i+1);
	double sx0, sx1, sx2, sx3;
	double sz0, sz1, sz2, sz3;
	if (sys->p.sflow_boundary_conditions == 0) {
		if (j == -1) {j += nz;}
		jp1 = (j == nz-1 ? 0 : j+1);
		kk[0] = i   + j*nx;
		kk[1] = ip1 + j*nx; // +dx
		kk[2] = i   + jp1*nx; // +dz
		kk[3] = ip1 + jp1*nx; // +dx+dz
		//
		sx0 = strain_rate_xx[kk[0]];
		sx1 = strain_rate_xx[kk[1]];
		sx2 = strain_rate_xx[kk[2]];
		sx3 = strain_rate_xx[kk[3]];
		//
		sz0 = strain_rate_zz[kk[0]];
		sz1 = strain_rate_zz[kk[1]];
		sz2 = strain_rate_zz[kk[2]];
		sz3 = strain_rate_zz[kk[3]];
	} else {
		if (j == -1) {
			//double d_ux_d_x = (u_x[ip1+j*nx]-u_x[k])/dx;
			//double d_uz_d_z = (u_z[i+jp1_nx]-u_z[k])/dz;
			//strain_rate_xx[k] = d_ux_d_x; // center
			//strain_rate_zz[k] = d_uz_d_z; // center
			//	kk[0] = i   + j*nx;
			//	kk[1] = ip1 + j*nx; // +dx
			kk[2] = i; // +dz     (j+1 = 0);
			kk[3] = ip1; // +dx+dz
			sx2 = strain_rate_xx[kk[2]];
			sx3 = strain_rate_xx[kk[3]];
			sx0 = -sx2;
			sx1 = -sx3;
			// Dzz (-1) = 0;
			sz0 = 0;
			sz1 = 0;
			sz2 = strain_rate_zz[kk[2]];
			sz3 = strain_rate_zz[kk[3]];
		} else if (j == nz -1) {
			kk[0] = i   + j*nx;
			kk[1] = ip1 + j*nx; // +dx
			//	kk[2] = i   + jp1*nx; // +dz
			//	kk[3] = ip1 + jp1*nx; // +dx+dz
			sx0 = strain_rate_xx[kk[0]];
			sx1 = strain_rate_xx[kk[1]];
			sx2 = -sx0;
			sx3 = -sx1;
			// Dzz (nz) = 0;
			sz0 = strain_rate_zz[kk[0]];
			sz1 = strain_rate_zz[kk[1]];
			sz2 = 0;
			sz3 = 0;
		} else {
			jp1 = j+1;
			kk[0] = i   + j*nx;
			kk[1] = ip1 + j*nx; // +dx
			kk[2] = i   + jp1*nx; // +dz
			kk[3] = ip1 + jp1*nx; // +dx+dz
			//
			sx0 = strain_rate_xx[kk[0]];
			sx1 = strain_rate_xx[kk[1]];
			sx2 = strain_rate_xx[kk[2]];
			sx3 = strain_rate_xx[kk[3]];
			//
			sz0 = strain_rate_zz[kk[0]];
			sz1 = strain_rate_zz[kk[1]];
			sz2 = strain_rate_zz[kk[2]];
			sz3 = strain_rate_zz[kk[3]];
		}
	}
	double sxx = sx0 + (sx1-sx0)*x_diff + (sx2-sx0)*z_diff + (sx3-sx2-sx1+sx0)*x_diff*z_diff;
	double szz = sz0 + (sz1-sz0)*x_diff + (sz2-sz0)*z_diff + (sz3-sz2-sz1+sz0)*x_diff*z_diff;
	double s = sxx-szz;
	//std::cerr << sxx + szz << std::endl;
	//e_local[0] = sxx;
	//e_local[2] = szz;
	e_local[0] = 0.5*s;
	e_local[2] = -0.5*s;
}

double SolventFlow::flowFiledDissipation()
{
	double energy_dissipation = 0;
	for (int k = 0; k<n; k++) {
		energy_dissipation +=   strain_rate_xx[k]*strain_rate_xx[k];
		energy_dissipation += 2*strain_rate_xz[k]*strain_rate_xz[k];
		energy_dissipation +=   strain_rate_zz[k]*strain_rate_zz[k];
	}
	return energy_dissipation/system_volume; // particle dynamic unit
}

double SolventFlow::particleDissipation()
{
	/*
	 * lubrication force x velocity need to be calculate.
	 */
	double energy_dissipation = 0;
	for (int i=0; i<sys->np_mobile; i++) {
		energy_dissipation += 0.5*sys->na_velocity.vel[i].sq_norm(); // W = F*U  (lu
	}
	return energy_dissipation/system_volume;
}

vec3d SolventFlow::calcAverageU()
{
	// particle unit.
	int i = std::lround(0.5*nx);
	double u_x_total = 0;
	for (int j=0; j<nz; j++) {
		int k = i + nx*j;
		u_x_total += u_x[k];//  + phi_ux[k]*Urel_x[k];
	}

	i = 0;
	double u_x_total0 = 0;
	for (int j=0; j<nz; j++) {
		int k = i + nx*j;
		u_x_total0 += u_x[k];
		//u_x_total0 += u_x[k] + phi_ux[k]*Urel_x[k];
	}
	if (0) {
		double err =  (u_x_total - u_x_total0)/u_x_total;
		//std::cerr << u_x_total << ' ' << u_x_total0 << std::endl;
		if (abs(err) > 0.001) {
			std::cerr << "err = " << err  << std::endl;
		}
	}
	
	int j = std::lround(0.5*nz);
	double u_z_total = 0;
	for (int i=0; i<nx; i++) {
		int k = i + nx*j;
		u_z_total += u_z[k];// + phi_uz[k]*Urel_z[k];
	}
	return vec3d(u_x_total/nz, 0 , u_z_total/nx);
}

void SolventFlow::velocityProfile(std::ofstream &fout_fp)
{
	// particle unit.
	Eigen::VectorXd ux(nz);
	Eigen::VectorXd ux_na(nz);
	Eigen::VectorXd phi(nz);
	for (int j=0; j<nz; j++) {
		double total_ux = 0;
		double total_ux_na = 0;
		double total_phi = 0;
		for (int i=0; i<nx; i++) {
			int k = i+nx*j;
			total_ux    += u_x[k];//         n + phi_ux[k]*Urel_x[k];
			total_ux_na += Urel_x[k];
			total_phi   += phi_ux[k];
		}
		ux(j) = total_ux/nx;
		ux_na(j) = total_ux_na/nx;
		// na = up - us
		phi(j) = total_phi/nx;
		// na + us = up - us + us = up
	}
	/*
	 * F_bf = up - u
	 * ux_na is the distributed up-u on mesh points.
	 * The total value (<ux_na>*[#mesh points] must be N*F_bf/
	 * std::cerr << "<up-u> = " << ux_na.mean()*nx*nz/sys->np_mobile << std::endl;
	 *
	 */
	std::cerr << "<up-us> = " << ux_na.mean()*nx*nz/sys->np_mobile << std::endl;
	for (int j=0; j<nz; j++) {
		fout_fp << pos[meshNb(0,j)].z << ' '; // 1: z
		fout_fp << ux(j) << ' ';              // 2: total velocity profile
		fout_fp << ux_na(j) + ux(j) << ' ';   // 3: u_particle
		fout_fp << phi(j);                    // 4: volume fraction
		fout_fp << std::endl;
	}
	fout_fp << std::endl;
		
}

void SolventFlow::outputYaplot(std::ofstream &fout_flow)
{
	double vfactor = 0.1;
	static bool first = true;
	if (first) {
		first = false;
	} else {
		fout_flow << std::endl;
	}
	vec3d pos0(sys->get_lx()/2,0, sys->get_lz()/2);
	fout_flow << "y 1" << std::endl;
	fout_flow << "@ 0" << std::endl;
	for (int i=0; i < sys->np_mobile; i++) {
		double x = sys->conf->position[i].x-pos0.x;
		double z = sys->conf->position[i].z-pos0.z;
		double a = sys->conf->radius[i];
		fout_flow << "r " << a << std::endl;
		fout_flow << "c " << x  << ' ' << 0 << ' ' << z << std::endl;
	}
//	fout_flow << "y 6" << std::endl;
//	fout_flow << "@ 2" << std::endl;
//	for (int i=sys->np_mobile; i < sys->get_np(); i++) {
//		double x = sys->conf->position[i].x-sys->get_lx()/2;
//		double z = sys->conf->position[i].z-sys->get_lz()/2;
//		double a = sys->conf->radius[i];
//		fout_flow << "r " << a << std::endl;
//		fout_flow << "c " << x  << ' ' << 0 << ' ' << z << std::endl;
//		if (sys->conf->position[i].z == 0) {
//			fout_flow << "c " << x  << ' ' << 0 << ' ' << z+sys->get_lz() << std::endl;
//		}
//	}
	
	fout_flow << "y 8" << endl;
	fout_flow << "@ 4" << endl;
	fout_flow << "r 0.3" << endl;
	for (const auto &inter: *(sys->interaction)) {
		unsigned int p0 = inter->get_p0();
		unsigned int p1 = inter->get_p1();
		double gap = inter->getReducedGap();
		if (gap < -0.1 ) {
			vec3d nvec = inter->getUnitSeparation()/inter->getUnitSeparation().norm();
			cerr << gap << ' ' << p0 << ' ' << p1 << ' ' << sys->np_mobile << endl;
			nvec.cerr();
			fout_flow << "s " << sys->conf->position[p0].x     -pos0.x << " 0 " << sys->conf->position[p0].z     -pos0.z;
			fout_flow << " "  << sys->conf->position[p0].x+nvec.x-pos0.x << " 0 " << sys->conf->position[p0].z+nvec.z-pos0.z << endl;

			fout_flow << "s " << sys->conf->position[p1].x     -pos0.x << " 0 " << sys->conf->position[p1].z         -pos0.z;
			fout_flow << " "  << sys->conf->position[p1].x-nvec.x-pos0.x << " 0 " << sys->conf->position[p1].z-nvec.z-pos0.z << endl;

		}
	}
	
	fout_flow << "y 7" << std::endl;
	fout_flow << "@ 3" << std::endl;
	for (int i=0; i < sys->np_mobile; i++) {
		double x = sys->conf->position[i].x-sys->get_lx()/2;
		double z = sys->conf->position[i].z-sys->get_lz()/2;
	//	double a = sys->conf->radius[i];
		vec3d v = vfactor*sys->vel_bg.vel[i];
		fout_flow << "l " << x  << ' ' << 0 << ' ' << z ;
		fout_flow << ' ' << x + v.x << ' ' << 0 << ' ' << z + v.z << std::endl;
	}
	
	fout_flow << "y 6" << std::endl;
	fout_flow << "@ 2" << std::endl;
	for (int i=sys->np_mobile; i < sys->get_np(); i++) {
		double x = sys->conf->position[i].x-sys->get_lx()/2;
		double z = sys->conf->position[i].z-sys->get_lz()/2;
		double a = sys->conf->radius[i];
		fout_flow << "r " << a << std::endl;
		fout_flow << "c " << x  << ' ' << 0 << ' ' << z << std::endl;
		if (sys->conf->position[i].z == 0) {
			fout_flow << "c " << x  << ' ' << 0 << ' ' << z+sys->get_lz() << std::endl;
		}
	}
	

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
					double vx = vfactor*Urel_x[k];
					fout_flow << "s ";
					fout_flow << po.x - dx/2      << ' ' << -0.02 << ' ' << po.z << ' ';
					fout_flow << po.x - dx/2 + vx << ' ' << -0.02 << ' ' << po.z << std::endl;
				}
				double vz = vfactor*Urel_z[k];
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
				//double vx = vfactor*(u_sol_x[k] +u_sol_x[meshNb(i+1,j)])/2;
				//double vz = vfactor*(u_sol_z[k] +u_sol_z[meshNb(i,j+1)])/2;
				
				double vx = vfactor*(u_x[k] +u_x[meshNb(i+1,j)])/2;
				double vz = vfactor*(u_z[k] +u_z[meshNb(i,j+1)])/2;

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
				double udx = vfactor*(Urel_x[k] +Urel_x[meshNb(i+1,j)])/2;
				double udz = vfactor*(Urel_z[k] +Urel_z[meshNb(i,j+1)])/2;
				fout_flow << "s ";
				fout_flow << po.x      << ' ' << -0.02 << ' ' << po.z << ' ';
				fout_flow << po.x + udx << ' ' << -0.02 << ' ' << po.z + udz << std::endl;
			}
		}
	}
	
	if (0) {
//		double cell_area = dx*dz;
		fout_flow << "@ 10" << std::endl;
		fout_flow << "y 5" << std::endl;
		double factor = 0.05;
		for (int j=0; j<nz_adpt; j++){
			for (int i=0; i<nx; i++){
				int k = i + nx*j;
				vec3d po = pos[k] - o_shift;
				if (j < nz) {
					//double r0 = phi_ux2[i   + nx*j];
					double r0 = factor*porousResistance(phi_ux[k]);

					fout_flow << "r " << abs(r0) << std::endl;
					fout_flow << "c ";
					fout_flow << po.x - dx/2 << ' ' << -0.02 << ' ' << po.z << std::endl;
				}
				//				double r2 = phi_uz2[k];
				double r2 = factor*porousResistance(phi_uz[k]);
				fout_flow << "r " << abs(r2) << std::endl;
				fout_flow << "c ";
				fout_flow << po.x   << ' ' << -0.02 << ' ' << po.z -dz/2<< std::endl;
				// double r = sqrt(cell_area*phi[i + nx*j]/M_PI);
//				int ip1 = (i == nx-1 ? 0 : i+1);
//				int jp1 = (j == nz-1 ? 0 : j+1);
//				int jp1 = j+1;

				
//				double r0 = porousResistance(phi_ux[i + nx*j]);
//				double r1 = porousResistance(phi_ux[ip1 + nx*j]);
//				double r2 = porousResistance(phi_uz[i + nx*j]);
//				double r3 = porousResistance(phi_uz[i + nx*jp1]);

			//	double r1 = phi_ux[ip1 + nx*j];

//				double r3 = phi_uz[i + nx*jp1];
				//double r = (r0 + r1 + r2 +r3)/4;
				//double r = 	gr_phi_Ud_phi_div_Ud[i + nx*j];
				//double s = 100*div_u_1_ast[i+nx*j];
				//				if (r > 0 ) {
				//					fout_flow << "@ 4" << std::endl;
				//				} else {
				//					fout_flow << "@ 6" << std::endl;
				//				}
	//			fout_flow << "r " << abs(r1) << std::endl;
		//		fout_flow << "c ";
//				fout_flow << po.x + dx/2 << ' ' << -0.02 << ' ' << po.z << std::endl;
//				fout_flow << "r " << abs(r3) << std::endl;
//				fout_flow << "c ";
	//			fout_flow << po.x   << ' ' << -0.02 << ' ' << po.z +dz/2<< std::endl;
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
			double z = i*dz-o_shift.z;
			fout_flow << "s ";
			fout_flow <<              -o_shift.x << ' ' << -0.02 << ' ' << z << ' ';
			fout_flow << sys->get_lx()-o_shift.x << ' ' << -0.02 << ' ' << z << std::endl;
		}
	}
		
	fout_flow << "y 4" << std::endl;
//	double p_max = 0;
//	for (int i=0; i< n; i++) {
//		if (p_max < pressure[i]) {
//			p_max = pressure[i];
//		}
//	}
//	double pr_max = 0;
//	for (auto pr : pressure) {
//		if (abs(pr) > pr_max) {
//			pr_max = abs(pr);
//		}
//	}

	for (int j=0; j<nz; j++){
		for (int i=0; i<nx; i++){
			//			double s =div_u_particle[i+nx*j];
			double s = pressure[i+nx*j]/pressure_difference_x;

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
	double pressure_total = 0;
	for (int k = 0; k<n; k++) {
		pressure_total += pressure[k];
	}
	cerr << "<p> =" << pressure_total / n << endl;
	
	if (0) {
		fout_flow << "y 9" << std::endl;
		for (int j=0; j<nz; j++){
			for (int i=0; i<nx; i++){
				//			double s =div_u_particle[i+nx*j];
				double s = 10*div_u_sol_ast[i+nx*j];
				if (s > 0 ) {
					fout_flow << "@ 4" << std::endl;
				} else {
					fout_flow << "@ 6" << std::endl;
				}
				//			cerr << s << endl;
				vec3d po = pos[i + nx*j] - o_shift;
				fout_flow << "r "<< abs(s) << std::endl;
				fout_flow << "c " << po.x << " -0.03 " << po.z << std::endl;
			}
		}
		
		
		fout_flow << "y 8" << std::endl;
		for (int j=0; j<nz; j++){
			for (int i=0; i<nx; i++){
				//			double s =div_u_particle[i+nx*j];
				double s = 10*div_u[i+nx*j];
				if (s > 0 ) {
					fout_flow << "@ 4" << std::endl;
				} else {
					fout_flow << "@ 6" << std::endl;
				}
				//			cerr << s << endl;
				vec3d po = pos[i + nx*j] - o_shift;
				fout_flow << "r "<< abs(s) << std::endl;
				fout_flow << "c " << po.x << " -0.03 " << po.z << std::endl;
			}
		}
	}
	
//	fout_flow << "y 9" << std::endl;
//	for (int j=0; j<nz_adpt; j++){
//		for (int i=0; i<nx; i++){
//			//			double s =div_u_particle[i+nx*j];
//			double s = omega[i+nx*j];
//			if (s > 0 ) {
//				fout_flow << "@ 4" << std::endl;
//			} else {
//				fout_flow << "@ 6" << std::endl;
//			}
//			vec3d po = pos[i + nx*j] - o_shift;
//			fout_flow << "r "<< abs(s) << std::endl;
//			fout_flow << "c " << po.x -dx/2<< " -0.03 " << po.z -dz/2<< std::endl;
//		}
//	}
}

int SolventFlow::meshNb(int xi, int zi)
{
	if (xi >= nx) {
		xi -= nx;
	} else if (xi <= -1) {
		xi += nx;
	}
	if (sys->p.sflow_boundary_conditions == 0) {
		if (zi >= nz) {
			zi -= nz ;
		} else if (zi <= -1) {
			zi += nz;
		}
	}
	int k = xi+zi*nx;
	if (k < 0 || k > n_adpt) {
		cerr << "k = " << k << ' ' << xi << ' ' << zi<<  endl;
		exit(1);
	}
	return k;
}

double SolventFlow::get_pressure_grad_x()
{
	/* flow unit --> particle dynamics unit
	 * grad p in PD unit = (grad p in SF unit) * (a/R0)^{3-d}
	 */
	return pressure_grad_x;
}
