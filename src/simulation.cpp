#include "simulation.h"
#include "aoup.h"
#include "force.h"
#include "sampling.h"
#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>
#include <omp.h>

void run_simulation(const SimParams& params) {
    int N = params.N_traj;
    int n_steps = static_cast<int>(params.T / params.aoup.dt);

    // Precompute constants that are the same for every trajectory / thread
    const double dt        = params.aoup.dt;
    const double decay     = 1.0 - dt / params.aoup.tau_c;
    const double noise_std = std::sqrt(2.0 * params.aoup.D_A / params.aoup.tau_c * dt);
    const double p         = params.force.p;
    const double b_min2    = params.force.b_min * params.force.b_min;
    const double b_min_sq  = b_min2;  // alias for clarity

    std::vector<double> A_vals(N), b_vals(N), w_vals(N);

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        std::mt19937_64 rng(params.seed + tid * 1000007ULL);
        std::normal_distribution<double> nd(0.0, 1.0);  // one object per thread, reused every step
        std::uniform_real_distribution<double> ud(0.0, 1.0);

        // Precompute sampling constants once per thread
        const double log_ratio  = std::log(params.sampling.b_max / params.sampling.b_min);
        const double sigma_v    = std::sqrt(params.aoup.D_A / params.aoup.tau_c);

        #pragma omp for schedule(dynamic, 64)
        for (int i = 0; i < N; ++i) {
            // Sample starting position
            double b     = params.sampling.b_min * std::exp(ud(rng) * log_ratio); // log-uniform (gamma=-2)
            double theta = ud(rng) * (2.0 * M_PI);
            double px    = b * std::cos(theta);
            double py    = b * std::sin(theta);

            // Initial AOUP state
            double x  = px, y  = py;
            double vx = sigma_v * nd(rng), vy = sigma_v * nd(rng);

            double A = 0.0;

            for (int s = 0; s < n_steps; ++s) {
                // Orientation (unit vector in velocity direction)
                double vmag = std::sqrt(vx*vx + vy*vy);
                double nx, ny;
                if (vmag > 1e-15) { nx = vx / vmag; ny = vy / vmag; }
                else              { nx = 1.0; ny = 0.0; }

                // Hydrodynamic dipole force x-component with hard-core cutoff
                double r2 = x*x + y*y;
                if (r2 >= b_min_sq) {
                    double r     = std::sqrt(r2);
                    double r3    = r2 * r;
                    double ndotx = nx*x + ny*y;
                    double cos2  = (ndotx * ndotx) / r2;
                    A += (p / r3) * (3.0 * cos2 - 1.0) * x * dt;
                }

                // Euler-Maruyama step
                x  += vx * dt;
                y  += vy * dt;
                vx  = vx * decay + noise_std * nd(rng);
                vy  = vy * decay + noise_std * nd(rng);
            }

            b_vals[i] = b;
            A_vals[i] = A;
            w_vals[i] = b * b;  // w(b) ~ b^(-gamma) = b^2 for gamma=-2
        }
    }

    // Normalise weights so mean(w) = 1
    double sum_w = 0.0;
    for (int i = 0; i < N; ++i) sum_w += w_vals[i];
    double inv_mean_w = N / sum_w;
    for (int i = 0; i < N; ++i) w_vals[i] *= inv_mean_w;

    // Weighted diagnostics
    double N_eff = effective_sample_size(w_vals.data(), N);

    double wsum = 0.0, mean_A = 0.0;
    for (int i = 0; i < N; ++i) { wsum += w_vals[i]; mean_A += w_vals[i] * A_vals[i]; }
    mean_A /= wsum;

    double var_A = 0.0, kurt_A = 0.0;
    for (int i = 0; i < N; ++i) {
        double d = A_vals[i] - mean_A;
        double d2 = d * d;
        var_A  += w_vals[i] * d2;
        kurt_A += w_vals[i] * d2 * d2;
    }
    var_A  /= wsum;
    kurt_A /= wsum * var_A * var_A;  // excess kurtosis denominator

    std::cout << "N_eff / N_traj    = " << N_eff / N << "\n";
    std::cout << "Weighted <A>      = " << mean_A   << " (should be ~0)\n";
    std::cout << "Weighted Var(A)   = " << var_A    << "\n";
    std::cout << "Weighted Kurt(A)  = " << kurt_A   << "\n";

    // Write config file alongside the data file
    {
        std::string config_path = params.output.substr(0, params.output.rfind('/') + 1) + "config.txt";
        std::ofstream cfg(config_path);
        cfg << "tau_c   = " << params.aoup.tau_c       << "\n"
            << "D_A     = " << params.aoup.D_A         << "\n"
            << "dt      = " << params.aoup.dt          << "\n"
            << "p       = " << params.force.p          << "\n"
            << "b_min   = " << params.force.b_min      << "\n"
            << "b_max   = " << params.sampling.b_max   << "\n"
            << "gamma   = " << params.sampling.gamma   << "\n"
            << "T       = " << params.T                << "\n"
            << "N_traj  = " << params.N_traj           << "\n"
            << "seed    = " << params.seed             << "\n";
        std::cout << "Config written to " << config_path << "\n";
    }

    // Write raw data
    std::ofstream out(params.output);
    out << "b,A,w\n";
    out.precision(15);
    for (int i = 0; i < N; ++i)
        out << b_vals[i] << "," << A_vals[i] << "," << w_vals[i] << "\n";
    std::cout << "Output written to " << params.output << "\n";
}
