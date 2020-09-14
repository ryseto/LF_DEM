//
//  MagneticInteractionParams.h
//  LF_DEM
//
//  Created by Zhiyuan Zhao on 2020/9/1.
//  Copyright © 2020 Ryohei Seto. All rights reserved.
//

#ifndef __LF_DEM__MagneticInteractionParams__
#define __LF_DEM__MagneticInteractionParams__

namespace Interactions
{

struct MagneticInteractionParams
{
    double beta_dipoledipole;                                       // Intensity of magnetic dipole-dipole interaction: \miu_0*m^2/(\pi*a^3*k_B*T)
//    bool janus_magnetic_particle;
//    int dipole_orient;
//    double dipole_shift;
};

inline bool has_magnetic_int(struct MagneticInteractionParams p)
{
    return p.beta_dipoledipole > 0;
}

//inline bool setMagneticParameters(MagneticInteractionParams &p)
//{
//    string indent = "  setupMagneticForce::\t";
//    if (p.magnetic_field_type == 1) {
//        cout << indent+"static homogenous magnetic field" << endl;
//    } else if (p.magnetic_field_type == 2) {
//        cout << indent+"static magnetic field gradient" << endl;
//    } else if (p.magnetic_field_type == 3) {
//        cout << indent+"alternating homogeneous magnetic field" << endl;
//    } else if (p.magnetic_field_type == 4) {
//        cout << indent+"rotating homogeneous magnetic field" << endl;
//    } else {
//        throw runtime_error(indent+"Error: unknown magnetic field type\n");
//    }
//}

}

#endif /* __LF_DEM__MagneticInteractionParams__ */
