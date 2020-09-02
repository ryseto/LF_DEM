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
    int magnetic_field_type;                            // Type of applied external magnetic field
    double langevin_parameter;                          // Ratio of magnetic energy to thermal energy: \miu_0 * m * H_0 / (k_B * T)
    double magnetic_field_freq;                         // Dimensionless angular frequency of magnetic field: 2 * \pi * f / Dr
    double beta_DipoleDipole;                           // Intensity of magnetic dipole-dipole interaction: \miu_0 * m^2 /(\pi * a^3 * k_B * T)
};

inline bool has_magnetic_int(struct MagneticInteractionParams p)
{
    return p.langevin_parameter > 0;
}

}

#endif /* __LF_DEM__MagneticInteractionParams__ */
