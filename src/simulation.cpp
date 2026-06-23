#include "simulation.h"
#include "force.h"
#include "sampling.h"
#include "swimmer_factory.h"
#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>
#include <omp.h>

void run_simulation(const SimParams& params) {
    int N = params.N_traj;
    int n_steps = static_cast<int>(params.T / params.process.dt);

    const double dt      = params.process.dt;
    const double p       = params.force.p;
    const double b_min_sq = params.force.b_min * params.force.b_min;
    const double log_ratio = std::log(params.sampling.b_max / params.sampling.b_min);

    std::vector<double> A_vals(N), b_vals(N), w_vals(N);

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        std::mt19937_64 rng(params.seed + tid * 1000007ULL);
        std::uniform_real_distribution<double> ud(0.0, 1.0);

        auto swimmer = make_swimmer(params.process);

        #pragma omp for schedule(dynamic, 64)
        for (int i = 0; i < N; ++i) {
            double b     = params.sampling.b_min * std::exp(ud(rng) * log_ratio);
            double theta = ud(rng) * (2.0 * M_PI);
            double x0    = b * std::cos(theta);
            double y0    = b * std::sin(theta);

            SwimmerState state = swimmer->init(x0, y0, rng);

            double A = 0.0;
            for (int s = 0; s < n_steps; ++s) {
                auto n  = swimmer->orientation(state);
                double x = state.x, y = state.y;
                double r2 = x*x + y*y;
                if (r2 >= b_min_sq) {
                    double r     = std::sqrt(r2);
                    double r3    = r2 * r;
                    double ndotx = n[0]*x + n[1]*y;
                    double cos2  = (ndotx * ndotx) / r2;
                    A += (p / r3) * (3.0 * cos2 - 1.0) * x * dt;
                }
                swimmer->step(state, rng);
            }

            b_vals[i] = b;
            A_vals[i] = A;
            w_vals[i] = b * b;  // w(b) ~ b^2 for gamma=-2 log-uniform sampling
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
        double d  = A_vals[i] - mean_A;
        double d2 = d * d;
        var_A  += w_vals[i] * d2;
        kurt_A += w_vals[i] * d2 * d2;
    }
    var_A  /= wsum;
    kurt_A /= wsum * var_A * var_A;

    std::cout << "N_eff / N_traj    = " << N_eff / N << "\n";
    std::cout << "Weighted <A>      = " << mean_A   << " (should be ~0)\n";
    std::cout << "Weighted Var(A)   = " << var_A    << "\n";
    std::cout << "Weighted Kurt(A)  = " << kurt_A   << "\n";

    // Write config
    {
        std::string config_path = params.output.substr(0, params.output.rfind('/') + 1) + "config_" + process_name() + ".txt";
        std::ofstream cfg(config_path);
        cfg << "process = " << process_name()          << "\n"
            << "tau_c   = " << params.process.tau_c    << "\n"
            << "D_A     = " << params.process.D_A      << "\n"
            << "v_A     = " << params.process.v_A      << "\n"
            << "D_r     = " << params.process.D_r      << "\n"
            << "omega   = " << params.process.omega     << "\n"
            << "beta    = " << params.process.beta      << "\n"
            << "tau_0   = " << params.process.tau_0     << "\n"
            << "dt      = " << params.process.dt        << "\n"
            << "p       = " << params.force.p           << "\n"
            << "b_min   = " << params.force.b_min       << "\n"
            << "b_max   = " << params.sampling.b_max    << "\n"
            << "gamma   = " << params.sampling.gamma    << "\n"
            << "T       = " << params.T                 << "\n"
            << "N_traj  = " << params.N_traj            << "\n"
            << "seed    = " << params.seed              << "\n";
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
