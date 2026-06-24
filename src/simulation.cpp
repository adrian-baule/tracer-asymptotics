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
    int N       = params.N_traj;
    int n_steps = static_cast<int>(params.T / params.process.dt);
    const double dt = params.process.dt;
    const double log_ratio = std::log(params.sampling.b_max / params.sampling.b_min);

#if defined(PROCESS_BMSHORT)
    // Gaussian force: precompute V_0/sigma^2 and 1/(2*sigma^2)
    const double inv_s2     = 1.0 / (params.force.sigma * params.force.sigma);
    const double force_pref = params.force.V_0 * inv_s2;
    const double half_inv_s2 = 0.5 * inv_s2;
    const int T_int = static_cast<int>(params.T);  // integer-time checkpoints

    // Thread-local accumulators for mu2, reduced after parallel region
    int n_threads_used = 1;
    #pragma omp parallel
    { n_threads_used = omp_get_num_threads(); }
    std::vector<std::vector<double>> tl_wA2(n_threads_used, std::vector<double>(T_int, 0.0));
    std::vector<std::vector<double>> tl_w  (n_threads_used, std::vector<double>(T_int, 0.0));
#else
    // Dipole force: precompute constants
    const double p         = params.force.p;
    const double b_min_sq  = params.force.b_min * params.force.b_min;
#endif

    std::vector<double> A_vals(N), b_vals(N), w_vals(N);

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        std::mt19937_64 rng(params.seed + tid * 1000007ULL);
        std::uniform_real_distribution<double> ud(0.0, 1.0);

        auto swimmer = make_swimmer(params.process);

#if defined(PROCESS_BMSHORT)
        std::vector<double> A_snap(T_int, 0.0);
#endif

        #pragma omp for schedule(dynamic, 64)
        for (int i = 0; i < N; ++i) {
            double b     = params.sampling.b_min * std::exp(ud(rng) * log_ratio);
            double theta = ud(rng) * (2.0 * M_PI);
            double x0    = b * std::cos(theta);
            double y0    = b * std::sin(theta);

            SwimmerState state = swimmer->init(x0, y0, rng);

            double A = 0.0;

#if defined(PROCESS_BMSHORT)
            // Accumulate A and record at integer-time checkpoints
            int next_snap = 0;  // next checkpoint index (t = next_snap + 1)
            for (int s = 0; s < n_steps; ++s) {
                double x = state.x, y = state.y;
                double r2 = x*x + y*y;
                A += force_pref * std::exp(-half_inv_s2 * r2) * x * dt;
                swimmer->step(state, rng);
                // Record A whenever time crosses an integer: t = (s+1)*dt
                double t_now = (s + 1) * dt;
                while (next_snap < T_int && (next_snap + 1) <= t_now + 1e-12) {
                    A_snap[next_snap] = A;
                    ++next_snap;
                }
            }
            // Fill any remaining checkpoints (if T is not exactly an integer)
            while (next_snap < T_int) { A_snap[next_snap] = A; ++next_snap; }

            double w = b * b;
            for (int k = 0; k < T_int; ++k) {
                tl_wA2[tid][k] += w * A_snap[k] * A_snap[k];
                tl_w  [tid][k] += w;
            }
#else
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
#endif

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

#if defined(PROCESS_BMSHORT)
    // Reduce thread-local mu2 accumulators and write mu2_bmshort.csv
    {
        std::string dir = params.output.substr(0, params.output.rfind('/') + 1);
        std::string mu2_path = dir + "mu2_bmshort.csv";
        std::ofstream mu2_out(mu2_path);
        mu2_out << "t,mu2\n";
        mu2_out.precision(15);
        for (int k = 0; k < T_int; ++k) {
            double sw2 = 0.0, sw = 0.0;
            for (int t = 0; t < n_threads_used; ++t) {
                sw2 += tl_wA2[t][k];
                sw  += tl_w  [t][k];
            }
            mu2_out << (k + 1) << "," << sw2 / sw << "\n";
        }
        std::cout << "mu2 written to " << mu2_path << "\n";
    }
#endif

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
            << "D_bm    = " << params.process.D_bm      << "\n"
            << "V_0     = " << params.process.V_0       << "\n"
            << "sigma   = " << params.process.sigma      << "\n"
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
