# Tracer Asymptotics

Simulation of passive tracer displacement statistics in a dilute active suspension,
via single-scattering increments from individual swimmer trajectories. All dynamics
are in **d = 2**.

## Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

Or a single binary directly (what the cluster runs use):

```bash
g++ -std=c++17 -O3 -march=native -fopenmp -DPROCESS_RTP_EXACT \
    -Iinclude src/main.cpp src/simulation.cpp src/force.cpp src/sampling.cpp \
    -o tracer_rtp_exact
```

One binary per process, selected by a compile-time `-DPROCESS_*` flag. Nine binaries:

| Binary | Swimmer | Force | Key parameters |
|--------|---------|-------|----------------|
| `tracer_aoup`   | Active Ornstein-Uhlenbeck | dipole | `--tau_c`, `--D_A` |
| `tracer_abp`    | Active Brownian Particle  | dipole | `--v_A`, `--D_r` |
| `tracer_rtp`    | Run-and-Tumble Particle   | dipole | `--v_A`, `--omega` |
| `tracer_levy1`  | Levy walk, beta > 1 | dipole | `--v_A`, `--beta` |
| `tracer_levy2`  | Levy walk, beta < 1 | dipole | `--v_A`, `--beta`, `--tau_0` |
| `tracer_bmshort`| 2D Brownian motion | Gaussian | `--D_bm`, `--V_0`, `--sigma` |
| `tracer_bmlong` | 2D Brownian motion | Coulomb  | `--D_bm`, `--sigma`, `--b_min`, `--mu` |
| `tracer_bmshort_exact` | 2D Brownian motion | Gaussian | `--D_bm`, `--V_0`, `--sigma`, `--mu_list` |
| `tracer_rtp_exact`     | Run-and-Tumble     | dipole   | `--v_A`, `--omega`, `--p`, `--mu_list` |

## Two tracer models

**Adiabatic (`A1`)** — the tracer is held at the origin; the swimmer `Y` evolves
freely and the increment is `A1(t) = mu * int_0^t F(Y_s) ds`. Used by the seven
non-`_exact` binaries.

**Exact (`_exact` binaries)** — the tracer is a state variable obeying
`dX/dt = mu * F(Y - X)`, with the force evaluated at the *true* separation. The
swimmer is unaffected by the tracer, so several `mu` values are integrated
simultaneously on one swimmer path (`--mu_list`, at most `MAX_MU = 8`), which is
far cheaper than separate runs and removes trajectory noise from `mu` comparisons.
Each exact run also accumulates `A1` in the same pass as a reference.

The exact and adiabatic results agree as `mu -> 0`: `X(mu)/mu -> A1`, verified to
5 digits at `mu = 1e-4`. This is a **weak-coupling (Picard) expansion in mu**, not
Kubo linear response — see "Coupling scale" below.

## Swimmer processes

- **AOUP**: velocity relaxes via OU process; `v_A = sqrt(D_A/tau_c)`.
- **ABP**: fixed speed `v_A`, orientation diffuses with rate `D_r`; `l_p = v_A/D_r`.
- **RTP**: fixed speed `v_A`, tumbles at Poisson rate `omega`; `l_p = v_A/omega`.
  2D effective diffusion `D_eff = v_A^2/(2*omega)`.
- **levy1**: `P(tau) ~ tau^{-(1+beta)}`, `beta > 1` (finite mean run time,
  diffusive). Minimum run time = `dt`.
- **levy2**: `P(tau > t) = (tau_0/t)^beta`, `0 < beta < 1` (divergent mean run
  time, superdiffusive). `tau_0` is physical and independent of `dt`; requires
  `tau_0 >= 10*dt` (the code warns below `100*dt`). Expected `M_{beta/2}`
  Mittag-Leffler limit with `beta/2 < 1/2`.
- **bmshort**: `dY = sqrt(2*D_bm)*dW`, radial Gaussian force
  `F(x) = (V_0/sigma^2)*exp(-|x|^2/(2*sigma^2))*x`. Isotropic, bounded, no
  singularity — the troubleshooting case.
- **bmlong**: same BM swimmer, 2D Coulomb force `F(x) = sigma*x/|x|^3`
  (`sigma` is the strength here, not a range). Marginal in 2D.

## Forces

Set by `ForceParams::type` in `include/force.h`; `force_x` and `force_2d` dispatch.

- `DIPOLE`: `F = (p/r^3)*(3*(n.xhat)^2 - 1) * x`, zero inside `r < b_min`.
  Orientation-dependent, `|F| ~ p/r^2`, sign-changing angular factor.
- `GAUSSIAN`: `F = (V_0/sigma^2)*exp(-r^2/(2*sigma^2)) * x`. Orientation ignored.
- `COULOMB`: `F = sigma*x/r^3`, with `r` clamped to `b_min`.

## Second-order Picard correction (A2)

The non-exact binaries also accumulate
`A2x(t) = -mu^2 * int_0^t [grad F(Y_u)^T Avec(u)]_x du`.

**Time ordering is critical**: `A2` must be updated with `Avec` *before* `Avec` is
advanced, so the inner integral excludes the current instant (`s < u`). Verified
against a brute-force double loop to 1e-10.

For the dipole processes the gradient is obtained by **central finite differences**
of `force_x` (`fd_grad_Fx`, step `--fd_step`, default 1e-4), holding orientation
fixed — so one code path serves all five dipole processes with no hand-derived
tensors. Points straddling `r = b_min` are flagged and the rate printed.
bmshort and bmlong use analytic gradients instead.

## Sampling

Impact parameter `b` is log-uniform on `[b_min, b_max]` with weight `w = b^2`
(normalised to mean 1), representing a uniform swimmer density in a disc of
radius `b_max`. Diagnostics `N_eff/N_traj`, weighted mean/Var/Kurt go to stdout.

`b_max` is **auto-computed** (sentinel `-1`) for the BM and `rtp_exact` processes
as `2*sqrt(4*D_eff*T)` — `D_eff = D_bm` for BM, `v_A^2/(2*omega)` for RTP — so the
sampling disc tracks the swimmer's diffusive spread. Override with `--b_max`.
Other processes default to `b_max = 50`.

## Coupling scale (dipole only)

The dipole force diverges as `p/r^2`, so the tracer's excursion during one
encounter at impact parameter `b` is `~ mu*p/(b*v_A)`. Requiring it to stay small
compared to `b` gives

    mu* = v_A * b_min^2 / p

Below `mu*` the exact result is a smooth correction to adiabatic; at and above it
the tracer can be entrained into the force-free core (`r < b_min`) and the
response becomes **non-monotonic in mu**. `rtp_exact` therefore defaults to
`b_min = 1.0` (not 0.5), which puts `mu* = 1` and keeps `dt = 1e-3` adequate near
the core. The bounded Gaussian force has no such threshold.

## Output

Files land in `--output_dir` (default `data/`, created if absent).

| File | Written by | Contents |
|------|-----------|----------|
| `samples_{process}.csv` | all except `rtp_exact` | `b,A1,A2,w` per trajectory (exact: `b,X_mu*,A1,w`) |
| `var_{process}.csv` | the 7 adiabatic processes | `t,VarA1,VarA2,meanA1,meanA2,ratio,kurtA1` |
| `var_exact_bmshort.csv` | `bmshort_exact` | `t,VarX_{mu},meanX_{mu},...,VarA1,meanA1` |
| `var_exact_rtp.csv` | `rtp_exact` | same layout |
| `config_{process}.txt` | all | parameters actually used |

All `var_*` files are checkpointed at **integer times** `t = 1..floor(T)`, so
`T = 1e5` gives 100k rows (~25 MB).

`rtp_exact` deliberately writes **no samples file** — at `T = 1e5` it is
impractically large and the variance statistics carry everything needed.

## Results so far

Both exact systems were run at `T = 1e5`, `N_traj = 1e5`,
`mu_list = 0.01,0.03,0.1,0.3,1`.

- **bmshort**: `Var ~ C*t` (normal diffusion). Local slope `dlnVar/dlnt = 1`.
- **rtp+dipole**: `Var ~ C*t*ln(t)`. Local slope tracks `1 + 1/ln(t)`
  parameter-free for `mu <= 0.3` — the dipole is marginal in 2D, like Coulomb.
- **Back-reaction enters at first order in mu**, not `mu^2`:
  `C(mu) = C_0*(1 + a*mu)` with `a = 0.48` (bmshort), `0.59` (rtp). The
  `A1`-`A2` cross-correlation does not vanish.
- The fitted intercept `C_0` reproduces the independently-accumulated adiabatic
  value to 0.01% (bmshort) and 0.4% (rtp) — a good end-to-end check.
- Linearity holds to `mu ~ 0.5`. Above it bmshort saturates *below* the linear
  extrapolation; rtp overshoots *above* it (3.0x vs 1.6x at `mu = 1`).
- **`rtp_exact` at `mu = 1` is not converged**: its curve still declines at
  `T = 1e5` and the plateau shifts 14% with the averaging window (vs ~5% for
  `mu < 1`). Exclude it from fits or run longer.

Plateaus are extracted as the mean of `Var/(mu^2 * t)` (or `/(mu^2 * t*ln t)`)
over `t` in `[1e4, 1e5]`; quote the spread across several windows as the error,
since within-window scatter understates it (the points are correlated).

## Numerical caveats

- **Timestep**: the dipole force is discontinuous at `r = b_min`, so the Euler
  sum converges slowly there. Absolute `Var(A1)` still shifts ~20% between
  `dt = 1e-3` and `6.25e-5`. Ratios and `t`-scaling are robust; treat absolute
  prefactors as unconverged.
- **Heavy tails**: dipole and small-`b_min` Gaussian increments have large
  kurtosis (30-60), so variance estimates are dominated by rare close
  encounters. `bmshort` (`b_min = 0.01`) is visibly noisier at late `t` than
  `rtp_exact` (`b_min = 1`).
- Weights span `[b_min^2, b_max^2]` before normalisation; very small `w` just
  means a close encounter, not a bug.

## Run

```bash
./build/tracer_aoup
./build/tracer_levy2 --beta 0.5 --tau_0 1.0
./build/tracer_rtp_exact --N_traj 100000 --T 10000 --mu_list 0.1,0.3,0.6,1.0
./build/tracer_bmshort_exact --N_traj 100000 --T 1000 --mu_list 0.03,0.1,0.3,1.0
```

Cost scales as `N_traj * T/dt`. The dipole force has no `exp` call and is cheap:
`N_traj = 1e5`, `T = 1e4` (1e7 steps/trajectory) runs in ~10 min on 48 cores.

Full CLI: `--tau_c --D_A --v_A --D_r --omega --beta --tau_0 --D_bm --V_0 --sigma
--mu --mu_list --fd_step --dt --p --b_min --b_max --gamma --T --N_traj --seed
--output_dir --output`

## Tests

```bash
ctest --test-dir build -V
```

`test_force`, `test_sampling`, `test_swimmers`, `test_bmlong`, `test_fd_a2`.
The last covers the FD gradient against analytic Coulomb, `A2` against a
brute-force double loop (dipole and Coulomb), and ensemble symmetry
`<A1> = <A2> = 0`.

## Cluster (Apocrita)

```bash
bash scripts/submit_apocrita.sh
```

Submits one SLURM job per process, with per-process `N_traj`, `T` and walltime.

## Known stale items

- `scripts/check_C.py` reads `data/mu2_bmshort.csv`, which is **no longer
  produced** — bmshort now writes `var_bmshort.csv` (different columns). The
  script needs updating before it will run.
- `src/aoup.cpp`, `include/aoup.h`, `tests/test_aoup.cpp` are **not in the
  build** and appear to be superseded by the `swimmer_factory` design.
- `config_bmshort_exact.txt` does not record `mu_vals` (only the rtp config
  does); read the `mu` values off the CSV header instead.
