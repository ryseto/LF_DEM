Unit scales
===========

Of course all those scales are not independent (choosing the time scale
should be enough), but I here list them all for convenience, as they are
used in the code.

WARNING: because of a current inconsistency, multiplying a dimensionless
viscosity by a dimesionless shear rate does *not* give a stress in the
correct stress units, that is, if :math:`\eta\dot\gamma= S` then
:math:`\tilde{\eta}\tilde{\dot\gamma} = 6\pi\tilde{S}`, where
:math:`\tilde{X}=X/X_0` is the dimensionless :math:`X`.

Some scales are common to all the scales sets: 

- length: particle radius :math:`a`
- viscosity: suspending fluid viscosity :math:`\eta_0`
- resistance (:math:`\b{R}_{\mathrm{FU}}`): :math:`6\pi\eta_0 a`
- resistance (:math:`\b{R}_{\mathrm{FE}}`): :math:`6\pi\eta_0 a^2`
- resistance (:math:`\b{R}_{\mathrm{SU}}`): :math:`6\pi\eta_0/a`

That for example means that in the code the Stokes drag is always
:math:`\tilde{\b{F}}_{\mathrm{SD}} = - (\tilde{\b{U}}-\tilde{\b{U}}_{\infty})`,
as the force-velocity part of the diagonal of
:math:`\tilde{\b{R}}_{\mathrm{FU}}` is 1.

Hydrodynamic scales
-------------------

-  time: inverse shear rate :math:`1/\dot\gamma`
-  velocity: :math:`a\dot\gamma`
-  force: :math:`6\pi\eta_0 \dot\gamma a^2`
-  stress: :math:`6\pi\eta_0 \dot\gamma`
-  temperature/energy: :math:`6\pi\eta_0 \dot\gamma a^3`

For example, with these units:

- the Brownian force is :math:`\tilde{\b{F}}_{\mathrm{B}}  = \sqrt{\frac{2}{\Delta\tilde{t}{\mathrm{Pe}}}}  \tilde{\b{L}}^{T} \cdot \b{\psi}`
- the stress coming from Einstein viscosity is :math:`(1+\frac{5}{2}\phi)/(6\pi)`

Brownian scales
---------------

-  time: diffusion time: :math:`6\pi\eta_0a^3/kT`
-  velocity: :math:`kT/6\pi\eta_0a^2`
-  force: :math:`kT/a`
-  stress: :math:`kT/a^3`
-  temperature/energy: :math:`kT`

For example, with these units:

- the Brownian force is :math:`\tilde{\b{F}}_{\mathrm{B}}  = \sqrt{\frac{2}{\Delta\tilde{t}}}  \tilde{\b{L}}^{T} \cdot \b{\psi}`
- the stress coming from Einstein viscosity is :math:`(1+\frac{5}{2}\phi)\mathrm{Pe}/(6\pi)`

Repulsive scales
----------------

Repulsive force :math:`F_R(h)=F_R^\ast \ e^{-h/\lambda_D}`

- time: :math:`6\pi\eta_0a^2/F_R^{\ast}`
- velocity: :math:`F_R^{\ast}/6\pi\eta_0a`
- force: :math:`F_R^{\ast}`
- stress: :math:`F_R^{\ast}/a^2`
- temperature/energy: :math:`F_R^{\ast}a`
