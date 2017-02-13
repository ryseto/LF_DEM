#Notes on the methods


Force balance:
\[
-R_{FU}^{lub}(U-U_{\infty}) - R_{FU}^{dash}U + R_{FE}^{lub}:E_{\infty} + F(X) = 0
\]
Calling $R_{FU} = R_{FU}^{lub}+R_{FU}^{dash}$, this is:
\[
-R_{FU}(U-U_{\infty}) - R_{FU}^{dash}U_{\infty} + R_{FE}^{lub}:E_{\infty} + F(X) = 0
\]

Stress:
\[
 S = R_{SU}^{lub}(U-U_{\infty})
 + R_{SE}^{lub}:E_{\infty}
 + x(F(X)-R_{FU}^{dash}U)
\]

This is solved in different ways depending on what is controlled (strain rate, shear stress, pressure, etc).

## Rate controlled

Force balance is solved as
\[
U-U_{\infty} =
R_{FU}^{-1}\left[ - R_{FU}^{dash}U_{\infty} + R_{FE}^{lub}:E_{\infty} + F(X) \right]
\]

We know $E_{\infty}$ and $U_{\infty}$, so we are done. From $U-U_{\infty}$, we can compute the stress.

## Stress controlled

We have to sort out the velocity pieces proportional to $\dot\gamma$ from those independent from it by decomposing $U=U_0 + \dot\gamma U_1$

Force balance is solved in two separate parts:
\[
  U_0 = R_{FU}^{-1}F(X)
\]
\[
  U_1 = U_{\infty}/\dot\gamma
  + R_{FU}^{-1}\left[ - R_{FU}^{dash}U_{\infty}/\dot\gamma + R_{FE}^{lub}:E_{\infty}/\dot\gamma \right]
\]

Note that we do not know $\dot\gamma$ but we do know both
$E_{\infty}/\dot\gamma$ and $U_{\infty}/\dot\gamma$,
so we can obtain explicitely $U_0$ and $U_1$ at this stage.

Now we can rewrite the stress as:
\[
S = R_{SU}^{lub}U_0 + x(F(X)-R_{FU}^{dash}U_0)
  + \dot\gamma \left[ R_{SU}^{lub}(U_1-U_{\infty}/\dot\gamma)
  + R_{SE}^{lub}:E_{\infty}/\dot\gamma
  - xR_{FU}^{dash}U_1
  \right]\\
  \equiv S_0 + \dot\gamma S_1
\]
Note that, like $U_0$ and $U_1$, $S_0$ and $S_1$
are completely independent of $\dot\gamma$,
and are computed explicitely.

To control the shear stress, we need to look only
at $\sigma = S:E_{\infty}/\dot\gamma$,
in which case the above relation is an affine equation for $\dot\gamma$.
It is trivial to solve:
\[
  \dot\gamma =
  \left(
    \sigma - S_0:E_{\infty}/\dot\gamma
  \right)\left(S_1:E_{\infty}/\dot\gamma\right)^{-1}
\]

After getting $\dot\gamma$, we can go back to the velocity
using $U=U_0 + \dot\gamma U_1$.
