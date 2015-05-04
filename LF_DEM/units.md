Unit scales
===========

Of course all those scales are not independent (choosing
the time scale should be enough), but I here list them all for
convenience, as they are used in the code.

WARNING: because of a current inconsistency, multiplying a dimensionless
viscosity by a dimesionless shear rate does *not* give a stress in the
correct stress units, that is, if $\eta\dot\gamma= S$
then $\tilde{\eta}\tilde{\dot\gamma} = 6\pi\tilde{S}$, where
$\tilde{X}=X/X_0$ is the dimensionless $X$.

Some scales are common to all the scales sets:
-   viscosity: suspending fluid viscosity $\eta_0$
-   resistance ($\boldsymbol{R}_{\mathrm{FU}}$): $6\pi\eta_0 a$
-   resistance ($\boldsymbol{R}_{\mathrm{FE}}$): $6\pi\eta_0 a^2$
-   resistance ($\boldsymbol{R}_{\mathrm{SU}}$): $6\pi\eta_0/a$

That for example means that in the code the Stokes drag is always $\tilde{\boldsymbol{F}}_{\mathrm{SD}} = - (\tilde{\boldsymbol{U}}-\tilde{\boldsymbol{U}}_{\infty})$, as the force-velocity part of the diagonal of $\tilde{\boldsymbol{R}}_{\mathrm{FU}}$ is 1.

Hydrodynamic scales
-------------------

-   length: particle radius $a$
-   time: inverse shear rate $1/\dot\gamma$
-   velocity: $a\dot\gamma$
-   force: $6\pi\eta_0 \dot\gamma a^2$
-   stress: $6\pi\eta_0 \dot\gamma$
-   temperature/energy: $6\pi\eta_0 \dot\gamma a^3$

For example, with these units:
- the Brownian force is $\tilde{\boldsymbol{F}}_{\mathrm{B}}
  = \sqrt{\frac{2}{\Delta\tilde{t}{\mathrm{Pe}}}}
  \tilde{\boldsymbol{L}}^{T} \cdot \boldsymbol{\psi}$
- the stress coming from Einstein viscosity is $(1+\frac{5}{2}\phi)/(6\pi)$

Brownian scales
---------------

-   length: particle radius $a$
-   time: diffusion time: $6\pi\eta_0a^3/kT$
-   velocity: $kT/6\pi\eta_0a^2$
-   force: $kT/a$
-   stress: $kT/a^3$
-   temperature/energy: $kT$

For example, with these units:
- the Brownian force is $\tilde{\boldsymbol{F}}_{\mathrm{B}}
  = \sqrt{\frac{2}{\Delta\tilde{t}}}
  \tilde{\boldsymbol{L}}^{T} \cdot \boldsymbol{\psi}$
- the stress coming from Einstein viscosity is $(1+\frac{5}{2}\phi)\mathrm{Pe}/(6\pi)$

Repulsive scales
----------------

Repulsive force $F_R(h)=F_R^\ast \ e^{-h/\lambda_D}$
-   length: particle radius $a$
-   time: $6\pi\eta_0a^2/F_R^{\ast}$
-   velocity: $F_R^{\ast}/6\pi\eta_0a$
-   force: $F_R^{\ast}$
-   stress: $F_R^{\ast}/a^2$
-   temperature/energy: $F_R^{\ast}a$
