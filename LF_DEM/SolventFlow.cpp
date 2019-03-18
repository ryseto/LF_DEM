//
//  SolventFlow.cpp
//  LF_DEM
//
//  Created by Ryohei Seto on 2019/03/16.
//  Copyright © 2019 Ryohei Seto. All rights reserved.
//

#include "SolventFlow.h"
#include "System.h"
using namespace std;

SolventFlow::SolventFlow()
{
	cerr << "." << endl;
	
}

SolventFlow::~SolventFlow()
{
//	DELETE(Pressure)
//	DELETE(u_solvent)
//	DELETE(u_particle)
}

void SolventFlow::init(System* sys_)
{
	sys = sys_;
	nx = 30;
	nz = 30;
	n = nx*nz;
	dx = sys->get_lx()/nx;
	dz = sys->get_lz()/nz;
	
//	Pressure = new double [nx*nz];
//	u_solvent = new vec3d [nx*nz];
//	u_particle = new vec3d [nx*nz];
	Pressure.resize(n);
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
}


void SolventFlow::calcVelocityDivergence()
{
	double dux_dx, duz_dz;
	for (int i=0; i<nx; i++){
		for (int j=0; j<nz; j++){
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

void SolventFlow::outputYaplot(std::ofstream &fout_flow)
{
	static bool first = true;
	if (first) {
		first = false;
	} else {
		fout_flow << endl;
	}
	fout_flow << "y 1" << endl;
	fout_flow << "@ 0" << endl;
	fout_flow << "r 1" << endl;
	for (int i=0; i < sys->get_np(); i++) {
		double x = sys->position[i].x-sys->get_lx()/2;
		double z = sys->position[i].z-sys->get_lz()/2;
		fout_flow << "c " << x  << ' ' << 0 << ' ' << z << endl;
	}
	fout_flow << "y 2" << endl;
	fout_flow << "@ 3" << endl;
	fout_flow << "r " << 0.1 << endl;
	double vfactor = 0.2;
	vec3d o_shift(sys->get_lx()/2, 0, sys->get_lz()/2);
	for (int i=0; i<nx; i++){
		for (int j=0; j<nz; j++){
			vec3d po = pos[i + nx*j] - o_shift;
			double vx = vfactor*u_particle_x[i + nx*j];
			double vz = vfactor*u_particle_z[i + nx*j];
			fout_flow << "s ";
			fout_flow << po.x - dx/2      << ' ' << -0.02 << ' ' << po.z << ' ';
			fout_flow << po.x - dx/2 + vx << ' ' << -0.02 << ' ' << po.z << endl;
			fout_flow << "s ";
			fout_flow << po.x << ' ' << -0.02 << ' ' << po.z - dz/2 << ' ';
			fout_flow << po.x << ' ' << -0.02 << ' ' << po.z - dz/2 + vz << endl;

		}
	}
	fout_flow << "y 3" << endl;
	fout_flow << "@ 5" << endl;
	fout_flow << "r " << 0.05 << endl;
	for (int i=0; i<nx; i++){
		double x = i*dx-o_shift.x;
		fout_flow << "s ";
		fout_flow << x << ' ' << -0.02 << ' ' <<              -o_shift.z << ' ';
		fout_flow << x << ' ' << -0.02 << ' ' << sys->get_lz()-o_shift.z << endl;
	}
	for (int i=0; i<nz; i++){
		double z =  i*dz-o_shift.z;
		fout_flow << "s ";
		fout_flow <<              -o_shift.x << ' ' << -0.02 << ' ' << z << ' ';
		fout_flow << sys->get_lx()-o_shift.x << ' ' << -0.02 << ' ' << z << endl;
	}

	fout_flow << "y 4" << endl;

	for (int i=0; i<nx; i++){
		for (int j=0; j<nz; j++){
			if (div_u_particle[i+nx*j] > 0 ) {
				fout_flow << "@ 4" << endl;
			} else {
				fout_flow << "@ 6" << endl;
			}
			vec3d po = pos[i + nx*j] - o_shift;
			fout_flow << "r "<< 0.4*abs(div_u_particle[i+nx*j]) << endl;
			fout_flow << "c " << po.x << " -0.03 " << po.z << endl;
		}
	}
}
