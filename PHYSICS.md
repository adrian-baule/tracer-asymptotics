# Physics and Equations

## Setup

A passive tracer is held fixed at the origin. A single AOUP starts at position
y = (b cos theta, b sin theta) and evolves freely. The tracer displacement increment
due to this single swimmer is (to first order in tracer mobility mu):

    A(y) = integral_0^T  f_x( X(u), n(u) )  du

where f_x is the x-projection of the hydrodynamic dipole force. The full tracer
displacement in a suspension of density rho is a Poisson sum over independent
such increments.

## AOUP dynamics

    dX/dt = V
    tau_c dV/dt = -V + sqrt(2 D_A) xi(t)

xi(t) is 2D Gaussian white noise, zero mean, unit variance per component.

Stationary velocity distribution: V ~ N(0, D_A/tau_c) per component.
Typical active speed: v_A = sqrt(D_A / tau_c).
Diffusion coefficient at long times: D_eff = D_A * tau_c.

Euler-Maruyama discretisation (timestep dt):

    X(t+dt) = X(t) + V(t) * dt
    V(t+dt) = V(t) * (1 - dt/tau_c) + sqrt(2 D_A / tau_c) * sqrt(dt) * N(0,1)

where N(0,1) is a 2D standard normal vector.

## Hydrodynamic dipole force

    F_hyd(x, n) = (p / |x|^3) * ( 3*(n.x)^2/|x|^2 - 1 ) * x

where:
  x = X(u)          vector from tracer (origin) to swimmer
  n = V(u)/|V(u)|   swimmer orientation unit vector
  p                 dipole strength (signed; p < 0 for pushers, p > 0 for pullers)
  |x|               Euclidean norm of x

The x-component used in A(y):

    f_x(x, n) = (p / |x|^3) * ( 3*(n.x)^2/|x|^2 - 1 ) * x_1

where x_1 is the x-coordinate of x.

Hard-core regularisation: if |x| < b_min, set f_x = 0.

Validity of first-order approximation: mu * |p| / b_min^2 << v_A.

## Importance sampling

True starting distribution in 2D (uniform spatial density):
    p_true(b) proportional to b,   b in [b_min, b_max]
    theta ~ Uniform[0, 2*pi)

Biased distribution (log-uniform in b, exponent gamma = -2):
    p_bias(b) proportional to b^(1+gamma) = 1/b

Sampling: b = b_min * (b_max/b_min)^U,  U ~ Uniform(0,1)

Importance weight:
    w(b) = p_true(b) / p_bias(b) proportional to b^2

Normalise: w_i <- w_i / mean(w) so that mean(w) = 1.

Weighted estimator for observable O:
    <O> = sum_i w_i O_i / sum_i w_i

Effective sample size:
    N_eff = (sum w_i)^2 / sum(w_i^2)

## Theoretical prediction

In the long-time limit T -> infinity, the Poisson-summed displacement distribution
converges to the bilateral Mittag-Leffler distribution M_{1/4}:

    p(y, t) = (1 / sqrt(C h(1/t))) * M_{1/4}( y / sqrt(C h(1/t)) )

where for d=2:
    h(s) = log(1/s) / (4 pi D_eff)    (Brownian motion result; AOUP generalisation TBD)

and C is the non-universal constant:
    C = integral dx integral dy  u(x,y) f(x) f(y)

with u(x,y) the subleading term in the Laplace-space propagator decomposition.

The moments of M_{1/4} are:
    integral z^n M_{1/4}(z) dz = (1 + (-1)^n)/2 * n! / Gamma(n/4 + 1)

In particular: <z^2> = 2 / Gamma(3/4), <z^4> = 24 / Gamma(2) = 24.

## Symmetry check

By symmetry of the dipole force under theta -> theta + pi and the uniform
distribution over theta, the mean displacement satisfies:
    <A> = 0

This is a basic sanity check for the simulation.

## Key dimensionless groups

    tau_c / (b_min^2 / D_A)   : ratio of persistence time to diffusion time across b_min
    p / (v_A b_min^2)         : dimensionless force strength at hard core
    b_max / b_min             : dynamic range of impact parameters
    T / tau_c                 : simulation duration in units of persistence time
    T * D_eff / b_max^2       : number of diffusive crossings of outer boundary
