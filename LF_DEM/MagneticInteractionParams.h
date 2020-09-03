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
    double beta_dipoledipole;                           // Intensity of magnetic dipole-dipole interaction: \miu_0 * m^2 /(\pi * a^3 * k_B * T)
};

inline bool has_magnetic_int(struct MagneticInteractionParams p)
{
    return p.beta_dipoledipole > 0;
}

}

#endif /* __LF_DEM__MagneticInteractionParams__ */
