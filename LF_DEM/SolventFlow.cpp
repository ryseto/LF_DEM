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
	nx = 30;
	nz = 30;
	n = nx*nz;
	dx = sys->get_lx()/nx;
	dz = sys->get_lz()/nz;
	pressure.resize(n);
	u_solvent_x.resize(n);
	u_solvent_z.resize(n);
	u_particle_x.resize(n);
	u_particle_z.resize(n);
	div_u_particle.resize(n);

	pos.resize(n);
	for (int iz=0; iz<nz; iz++) {
		for (int ix=0; ix<nx; ix++) {
			pos[ix+nx*iz].set(dx*ix+dx/2, 0, dz*iz+dz/2);
		}
	}
}

void SolventFlow::initPoissonSolver()
{
	std::vector< std::vector<double> > lmat;
	lmat.resize(n);
	for (int l=0; l<n; l++){
		lmat[l].resize(n);
	}
	for (auto l : lmat) {
		for (auto elm : l) {
			elm = 0;
		}
	}
	double exz = 2*(1/(dx*dx)+1/(dz*dz));
	double ex = 1/(dx*dx);
	double ez = 1/(dz*dz);

	for (int zi = 0; zi < nz; zi++){
		for (int xi = 0; xi < nx; xi++){
			// b(xi, zi)
			int xip1 = xi+1;
			if (xip1 == nx) {
				xip1 = 0;
			}
			int xim1 = xi-1;
			if (xim1 == -1) {
				xim1 = nx-1;
			}
			int zip1 = zi+1;
			if (zip1 == nz) {
				zip1 = 0;
			}
			int zim1 = zi-1;
			if (zim1 == -1) {
				zim1 = nz-1;
			}
			lmat[q(xi, zi)][q(xi,   zi  )] = -exz;  //p(xi,zi) --> b(xi,zi)
			lmat[q(xi, zi)][q(xip1, zi  )] =   ex;  //p(xi+1,zi) --> b(xi,zi)
			lmat[q(xi, zi)][q(xim1, zi  )] =   ex;  //p(xi-1,zi) --> b(xi,zi)
			lmat[q(xi, zi)][q(xi  , zip1)] =   ez;  //p(xi,zi+1) --> b(xi,zi)
			lmat[q(xi, zi)][q(xi  , zim1)] =   ez;  //p(xi,zi-1) --> b(xi,zi)
		}
	}
	if (true) {
		int cnt_nz = 0;
		for (int l=0; l<n; l++) {
			for (int k=0; k<n; k++) {
				if ( lmat[l][k] != 0) {
					cnt_nz ++;
				}
			}
		}
		std::cerr << "cnt_nz=  " << cnt_nz << std::endl;
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
	std::ofstream mat_out("mat.dat");
	for (int l=0; l<n; l++) {
		for (int k=0; k<n; k++) {
			int i = k%nx;
			int j = (k - i)/nx;
			mat_out << lmat[l][k] << ' ';
			
		}
		mat_out << std::endl;
	}
	mat_out.close();
}

void SolventFlow::updateParticleVelocity()
{
	for (int k=0; k<n; k++) {
		u_particle_x[k] = 0;
		u_particle_z[k] = 0;
	}
	double x, z;
	double vx, vz;

	double ave_range = 1.5;
	double sq_ave_range = ave_range*ave_range;

	int i_range = 4*ave_range/dx;
	for (int i=0; i < sys->get_np(); i++) {
		x = sys->position[i].x;
		z = sys->position[i].z;
		vx = sys->velocity[i].x;
		vz = sys->velocity[i].z;
		int ix = x /dx;
		int iz = z /dz;
		for (int l=-i_range; l<=i_range; l++){
			for (int m=-i_range; m<=i_range; m++){
				int iix = ix+l;
				int iiz = iz+m;
				double xo = x;
				if (iiz >= 0 && iiz < nz) {
					if (iix < 0) {
						iix += nx;
						xo += sys->get_lx();
					} else if (iix >= nx) {
						iix -= nx;
						xo -= sys->get_lx();
					}
					double x_p = pos[iix +nx*iiz].x;
					double x_u = x_p -dx/2;
					double z_p = pos[iix +nx*iiz].z;
					double z_u = z_p -dz/2;
					double dist_sq_ux = (xo-x_u)*(xo-x_u) + (z-z_p)*(z-z_p);
					double dist_sq_uz = (xo-x_p)*(xo-x_p) + (z-z_u)*(z-z_u);
					double weight_ux = exp(-dist_sq_ux/sq_ave_range);
					double weight_uz = exp(-dist_sq_uz/sq_ave_range);
					u_particle_x[iix + nx*iiz] += weight_ux*vx;
					u_particle_z[iix + nx*iiz] += weight_uz*vz;
				}
			}
		}
	}
	calcVelocityDivergence();
	solvePressure();
}

void SolventFlow::calcVelocityDivergence()
{
	double dux_dx, duz_dz;
	for (int j=0; j<nz; j++){
		for (int i=0; i<nx; i++){
			int i_next = i+1;
			if (i_next == nx) {
				i_next = 0;
			}
			int j_next = j+1;
			if (j_next == nz) {
				j_next = 0;
			}
			dux_dx = (u_particle_x[i_next + j*nx] - u_particle_x[i + j*nx])/dx;
			duz_dz = (u_particle_z[i + j_next*nx] - u_particle_z[i + j*nx])/dz;
			div_u_particle[i+nx*j] = dux_dx + duz_dz;
		}
	}
}

void SolventFlow::solvePressure()
{
	double delta_P_dxdx = 1;
	for (int j=0; j<nz; j++){
		for (int i=0; i<nx; i++){
			//b.setConstant(i, random());
			b(i+nx*j) = div_u_particle[i+nx*j];
			if (i == 0) {
				b(i+nx*j) += -delta_P_dxdx;
			} else if (i == nx-1) {
				b(i+nx*j) += delta_P_dxdx;
			}
		}
	}
	x = psolver->solve(b);
	if (psolver->info() != Eigen::Success) {
		// solving failed
		std::cerr << "solving failed" << std::endl;
		return;
	}
	for (int i=0; i < n; i++) {
		pressure[i] = x(i);
	}
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
	fout_flow << "r 1" << std::endl;
	for (int i=0; i < sys->get_np(); i++) {
		double x = sys->position[i].x-sys->get_lx()/2;
		double z = sys->position[i].z-sys->get_lz()/2;
		fout_flow << "c " << x  << ' ' << 0 << ' ' << z << std::endl;
	}
	fout_flow << "y 2" << std::endl;
	fout_flow << "@ 3" << std::endl;
	fout_flow << "r " << 0.1 << std::endl;
	double vfactor = 0.2;
	vec3d o_shift(sys->get_lx()/2, 0, sys->get_lz()/2);
	for (int i=0; i<nx; i++){
		for (int j=0; j<nz; j++){
			vec3d po = pos[i + nx*j] - o_shift;
			double vx = vfactor*u_particle_x[i + nx*j];
			double vz = vfactor*u_particle_z[i + nx*j];
			fout_flow << "s ";
			fout_flow << po.x - dx/2      << ' ' << -0.02 << ' ' << po.z << ' ';
			fout_flow << po.x - dx/2 + vx << ' ' << -0.02 << ' ' << po.z << std::endl;
			fout_flow << "s ";
			fout_flow << po.x << ' ' << -0.02 << ' ' << po.z - dz/2 << ' ';
			fout_flow << po.x << ' ' << -0.02 << ' ' << po.z - dz/2 + vz << std::endl;

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
	std::cerr << p_max << std::endl;
	for (int j=0; j<nz; j++){
		for (int i=0; i<nx; i++){
			//			double s =div_u_particle[i+nx*j];
			double s = pressure[i+nx*j]/10;
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
