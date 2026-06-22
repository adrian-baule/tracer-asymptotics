# Parameter Choices and Physical Reasoning

## Suggested default parameter set

| Parameter | Default | Units       | Notes                                      |
|-----------|---------|-------------|--------------------------------------------|
| tau_c     | 1.0     | time        | sets the time unit                         |
| D_A       | 1.0     | length^2/time | v_A = sqrt(D_A/tau_c) = 1.0             |
| p         | 1.0     | length^2    | dipole strength; try p = -1.0 for pushers  |
| b_min     | 0.5     | length      | hard-core cutoff                           |
| b_max     | 50.0    | length      | outer cutoff                               |
| gamma     | -2.0    | -           | importance sampling exponent               |
| dt        | 1e-3    | time        | timestep = 1e-3 tau_c                      |
| T         | 1e3     | time        | total time per trajectory = 1e3 tau_c      |
| N_traj    | 1e5     | -           | number of trajectories                     |
| seed      | 42      | -           | RNG seed for reproducibility               |

## Reasoning for each parameter

### tau_c = 1.0, D_A = 1.0
Sets natural units. The active speed is v_A = 1, the diffusion coefficient is
D_eff = D_A * tau_c = 1. All lengths are in units of v_A * tau_c (persistence length).

### b_min = 0.5
Must satisfy mu * |p| / b_min^2 << v_A. With mu = 1, p = 1, v_A = 1:
    1.0 / 0.25 = 4.0
This is not deeply in the valid regime. To be safer use b_min = 1.0 or 2.0 for
initial tests, then decrease to check sensitivity.

### b_max = 50.0
Contributions from b >> sqrt(D_eff * T) = sqrt(1000) ~ 31.6 are negligible
since the AOUP cannot diffuse that far and return in time T. b_max = 50 is
slightly above this scale. Can be reduced to 30 for speed without much loss.

### gamma = -2.0
Gives p_bias(b) proportional to 1/b (log-uniform). This ensures equal statistical
weight per decade of b. Since the force scales as 1/b^2 and the interaction time
scales as b^2/D_eff, the contribution to A scales as p/b^2 * b^2/D_eff = p/D_eff,
roughly independent of b. Log-uniform sampling is therefore close to optimal.
Monitor N_eff/N_traj: if < 10%, try gamma = -1.5.

### dt = 1e-3
Must satisfy:
  (a) dt << tau_c = 1.0              (resolve AOUP dynamics)        -> dt << 1
  (b) dt << b_min^2 / D_A = 0.25    (resolve close encounters)     -> dt << 0.25
  (c) dt << b_min / v_A = 0.5       (resolve transit across b_min) -> dt << 0.5
dt = 1e-3 satisfies all three comfortably.

### T = 1e3
Must be long enough for the AOUP to return to the vicinity of the origin many times.
The return time to within b_min ~ 0.5 from a starting distance b scales as b^2/D_eff.
For b_max = 50, this is 2500 tau_c -- longer than T = 1000.
This means trajectories starting near b_max will not have returned by time T.
However, their contribution to A is small (they start far away). Check convergence
by running T = 100, 1000, 10000 and comparing distributions.

### N_traj = 1e5
Gives O(10^5) weighted samples. With N_eff/N_traj ~ 10-30% (expected for gamma = -2),
effective sample size is O(10^4). This should resolve the distribution body but not
the far tails. Increase to 1e6 for publication-quality results.

## Convergence checks to run

1. Vary T in {100, 1000, 10000}: distribution should stabilise.
2. Vary b_min in {0.1, 0.5, 1.0, 2.0}: check sensitivity (non-universal constant C
   depends on b_min, but the shape M_{1/4} should not).
3. Vary b_max in {20, 50, 100}: distribution should be insensitive if b_max >> sqrt(D_eff * T).
4. Vary gamma in {0, -1, -2, -3}: monitor N_eff and distribution consistency.
5. Vary N_traj in {1e4, 1e5, 1e6}: check statistical convergence.

## Estimated runtime

With T = 1e3, dt = 1e-3: 1e6 steps per trajectory.
With N_traj = 1e5: 1e11 total steps.
At ~1e8 steps/second/thread in optimised C++: ~1000 seconds single-threaded.
With 8 OpenMP threads: ~125 seconds.

For initial testing use T = 1e2 (1e5 steps/trajectory):
With N_traj = 1e5 and 8 threads: ~1 second. Good for debugging.
