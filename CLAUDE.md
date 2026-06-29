# AOUP Tracer Simulation

Simulation of passive tracer displacement statistics in a dilute active suspension,
via single-scattering increments from individual swimmer trajectories.

## Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

Produces five binaries in `build/`:

| Binary | Process | Key parameters |
|--------|---------|---------------|
| `tracer_aoup`  | Active Ornstein-Uhlenbeck | `--tau_c`, `--D_A` |
| `tracer_abp`   | Active Brownian Particle  | `--v_A`, `--D_r` |
| `tracer_rtp`   | Run-and-Tumble Particle   | `--v_A`, `--omega` |
| `tracer_levy1` | Levy walk, beta > 1 (diffusive)      | `--v_A`, `--beta` |
| `tracer_levy2` | Levy walk, beta < 1 (superdiffusive) | `--v_A`, `--beta`, `--tau_0` |
| `tracer_bmshort` | 2D Brownian motion, Gaussian force | `--D_bm`, `--V_0`, `--sigma` |
| `tracer_bmlong`  | 2D Brownian motion, Coulomb force, A1+A2 diagnostic | `--D_bm`, `--sigma`, `--b_min`, `--mu` |

## Swimmer processes

- **AOUP**: velocity relaxes via OU process; `v_A = sqrt(D_A/tau_c)`.
- **ABP**: fixed speed `v_A`, orientation diffuses with rate `D_r`; persistence length `l_p = v_A/D_r`.
- **RTP**: fixed speed `v_A`, tumbles at Poisson rate `omega`; persistence length `l_p = v_A/omega`.
- **levy1**: fixed speed `v_A`, run times drawn from `P(tau) ~ tau^{-(1+beta)}` with `beta > 1`
  (finite mean run time, diffusive long-time behaviour). Minimum run time = `dt`.
- **bmshort**: 2D Brownian motion swimmer (`dX = sqrt(2*D_bm)*dW`), radial Gaussian force
  `f(x) = (V_0/sigma^2)*exp(-|x|^2/(2*sigma^2))*x1`. Used for troubleshooting: force is isotropic,
  no singularity, theory fully worked out. Limit distribution = M_0 (bilateral exponential),
  CF = 1/(1+k^2). Key check: `<A^2>(T) ~ C*log(T)` (logarithmic growth in T). The log scaling
  arises because 2D BM is null-recurrent: the occupation time of the Gaussian force region
  (radius ~ sigma) starting from impact parameter b scales as `(sigma^2/D)*log(DT/b^2)`, so
  the variance of the scattering increment grows as `log(T)` rather than `T`.
  Extra output: `mu2_bmshort.csv` with weighted `<A(t)^2>` at integer times t=1,...,T.
  Validate with `python scripts/check_C.py` (fits C*log(t) to large-t data).
- **bmlong**: 2D BM swimmer (`dX = sqrt(2*D_bm)*dW`), 2D Coulomb force
  `F(x) = sigma*x/|x|^3`. Back-reaction diagnostic: accumulates both
  `A1x(t) = mu * int_0^t F_x(Y_s) ds` (first-order / adiabatic) and
  `A2x(t) = -mu^2 * int_0^t [grad F(Y_u)^T Avec(u)]_x du` (second-order Picard correction).
  Outputs `data/var_bmlong.csv` with `t, VarA1, VarA2, meanA1, meanA2, ratio`.
  Key check: `ratio = Var(A2)/Var(A1)` vs `log(t)`.
  Flat => adiabatic truncation benign; linear in `log(t)` => logarithmic back-reaction promotion.
  Theory: `Var(A1) ~ t*log(t)` (marginal Coulomb in 2D); `Var(A2) ~ t*(log t)^a`, `1 <= a <= 3`.
  `b_max = sqrt(4*D_bm*T)` (auto-computed; scales with sqrt(T)).
  Validate with `python scripts/check_backreaction.py`.
  Note: A2 is O(mu^2) and noisier; use large N_traj for clean large-t statistics.

- **levy2**: fixed speed `v_A`, run times drawn from `P(tau > t) = (tau_0/t)^beta` with `0 < beta < 1`
  (divergent mean run time, superdiffusive long-time behaviour). `tau_0` is an explicit physical
  parameter (minimum run time, independent of `dt`); requires `tau_0 >= 10*dt`.
  Expected to exhibit `M_{beta/2}` Mittag-Leffler limit distribution with `beta/2 < 1/2`.

## Run

```bash
mkdir -p data
./build/tracer_aoup
./build/tracer_levy2 --beta 0.5 --tau_0 1.0
```

Output per run: `data/output_{process}.csv` and `data/config_{process}.txt`.

## Tests

```bash
ctest --test-dir build -V
```

## Cluster (Apocrita)

```bash
bash scripts/submit_apocrita.sh
```
