# AOUP Tracer Displacement Statistics

Simulation of the single-scattering displacement distribution of a passive tracer
interacting with an active Ornstein-Uhlenbeck particle (AOUP) via a hydrodynamic
dipole force in d=2.

## Physics background

The project tests the universal displacement statistics derived in:

- Baule (2026), "Comment on Universal dynamics of a passive particle driven by Brownian motion"
- Baule (2023), "Universal Poisson statistics of a passive tracer diffusing in dilute active suspensions", PNAS 120, e2308226120

To first order in swimmer density rho and tracer mobility mu, the tracer displacement
is a Poisson sum of independent single-scattering increments

    A(y) = integral_0^T  f_x( X(u), n(u) )  du

where X(u) is a free AOUP trajectory starting at y, n(u) = V(u)/|V(u)| is the
instantaneous orientation, and f_x is the x-projection of the hydrodynamic dipole force.

In the long-time limit, the rescaled displacement distribution is predicted to converge
to the bilateral Mittag-Leffler distribution M_{1/4}, universal across swimmer models
and force profiles (subject to a short-range decay condition).

## Repository structure

```
aoup_tracer/
├── README.md
├── CMakeLists.txt
├── PHYSICS.md                  # Detailed physics and equations
├── src/
│   ├── main.cpp                # Entry point, parameter parsing
│   ├── simulation.cpp          # Main simulation loop
│   ├── aoup.cpp                # AOUP integrator (Euler-Maruyama)
│   ├── force.cpp               # Hydrodynamic dipole force
│   └── sampling.cpp            # Importance sampling utilities
├── include/
│   ├── simulation.h
│   ├── aoup.h
│   ├── force.h
│   └── sampling.h
├── scripts/
│   ├── plot_distribution.py    # Plot weighted histogram vs theory
│   ├── plot_convergence.py     # Check convergence with T
│   └── mittag_leffler.py       # Reference M_{1/4} distribution
├── tests/
│   ├── test_force.cpp          # Unit tests for force evaluation
│   ├── test_aoup.cpp           # Unit tests for AOUP integrator
│   └── test_sampling.cpp       # Unit tests for importance weights
├── docs/
│   ├── requirements.txt        # Full simulation requirements (plain text)
│   └── parameters.md           # Parameter choices and physical reasoning
└── data/                       # Output directory (gitignored)
```

## Building

Requires: C++17 compiler, OpenMP

```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j4
```

## Running

```bash
./aoup_tracer --tau_c 1.0 --D_A 1.0 --p 1.0 \
              --b_min 0.5 --b_max 50.0 \
              --dt 1e-3 --T 1e3 \
              --N_traj 100000 \
              --gamma -2.0 \
              --output data/output.csv \
              --seed 42
```

See `docs/parameters.md` for guidance on parameter choices.

## Output

The simulation writes a CSV file with columns:

    b, A, w

where b is the impact parameter, A is the scattering increment, and w is the
normalised importance weight. Post-processing scripts in `scripts/` compute
weighted histograms and compare against the theoretical prediction.

## Validation

Run the convergence check first:

```bash
python scripts/plot_convergence.py --data_dir data/
```

This compares distributions obtained with T = 100, 1000, 10000 to check
that the limit distribution has been reached.
