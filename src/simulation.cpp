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

    std::vector<double> A_vals(N), b_vals(N), w_vals(N);

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        std::mt19937_64 rng(params.seed + tid * 1000007ULL);

        #pragma omp for schedule(dynamic, 64)
        for (int i = 0; i < N; ++i) {
            double b     = sample_b(params.sampling, rng);
            double theta = sample_theta(rng);
            double x0    = b * std::cos(theta);
            double y0    = b * std::sin(theta);

            AOUPState state = aoup_init(x0, y0, params.aoup, rng);

            double A = 0.0;
            for (int s = 0; s < n_steps; ++s) {
                auto n = aoup_orientation(state);
                A += force_x(state.x, state.y, n[0], n[1], params.force) * params.aoup.dt;
                aoup_step(state, params.aoup, rng);
            }

            b_vals[i] = b;
            A_vals[i] = A;
            w_vals[i] = importance_weight(b, params.sampling);
        }
    }

    // Normalise weights
    double sum_w = 0.0;
    for (int i = 0; i < N; ++i) sum_w += w_vals[i];
    for (int i = 0; i < N; ++i) w_vals[i] /= (sum_w / N); // mean(w) = 1

    // Diagnostics
    double N_eff = effective_sample_size(w_vals.data(), N);
    double mean_A = 0.0, var_A = 0.0;
    double sum_w2 = 0.0;
    for (int i = 0; i < N; ++i) {
        mean_A += w_vals[i] * A_vals[i];
        sum_w2 += w_vals[i];
    }
    mean_A /= sum_w2;
    for (int i = 0; i < N; ++i)
        var_A += w_vals[i] * (A_vals[i] - mean_A) * (A_vals[i] - mean_A);
    var_A /= sum_w2;

    std::cout << "N_eff / N_traj = " << N_eff / N << "\n";
    std::cout << "Weighted mean(A) = " << mean_A << " (should be ~0)\n";
    std::cout << "Weighted var(A)  = " << var_A << "\n";

    // Write output
    std::ofstream out(params.output);
    out << "b,A,w\n";
    for (int i = 0; i < N; ++i)
        out << b_vals[i] << "," << A_vals[i] << "," << w_vals[i] << "\n";
    std::cout << "Output written to " << params.output << "\n";
}
