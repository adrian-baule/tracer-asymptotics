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
    const double inv_s2      = 1.0 / (params.force.sigma * params.force.sigma);
    const double force_pref  = params.force.V_0 * inv_s2;
    const double half_inv_s2 = 0.5 * inv_s2;
    const int T_int = static_cast<int>(params.T);

    int n_threads_used = 1;
    #pragma omp parallel
    { n_threads_used = omp_get_num_threads(); }
    // Accumulators: w, A1 (mean/var/kurt), A2 (mean/var), matching bmlong column order
    std::vector<std::vector<double>> tl_w    (n_threads_used, std::vector<double>(T_int, 0.0));
    std::vector<std::vector<double>> tl_wA1  (n_threads_used, std::vector<double>(T_int, 0.0));
    std::vector<std::vector<double>> tl_wA1sq(n_threads_used, std::vector<double>(T_int, 0.0));
    std::vector<std::vector<double>> tl_wA1_4(n_threads_used, std::vector<double>(T_int, 0.0));
    std::vector<std::vector<double>> tl_wA2  (n_threads_used, std::vector<double>(T_int, 0.0));
    std::vector<std::vector<double>> tl_wA2sq(n_threads_used, std::vector<double>(T_int, 0.0));

#elif defined(PROCESS_BMLONG)
    // Coulomb force parameters
    const double sigma_c = params.force.sigma;  // Coulomb strength
    const double b_min_c = params.force.b_min;
    const double mu      = params.mu;
    const int T_int = static_cast<int>(params.T);

    int n_threads_used = 1;
    #pragma omp parallel
    { n_threads_used = omp_get_num_threads(); }
    // Accumulators per thread: w, A1, A1sq, A1_4 (kurtosis), A2, A2sq
    std::vector<std::vector<double>> tl_w    (n_threads_used, std::vector<double>(T_int, 0.0));
    std::vector<std::vector<double>> tl_wA1  (n_threads_used, std::vector<double>(T_int, 0.0));
    std::vector<std::vector<double>> tl_wA1sq(n_threads_used, std::vector<double>(T_int, 0.0));
    std::vector<std::vector<double>> tl_wA1_4(n_threads_used, std::vector<double>(T_int, 0.0));
    std::vector<std::vector<double>> tl_wA2  (n_threads_used, std::vector<double>(T_int, 0.0));
    std::vector<std::vector<double>> tl_wA2sq(n_threads_used, std::vector<double>(T_int, 0.0));
    long long total_clamps = 0;

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
        std::vector<double> A1_snap(T_int, 0.0);
        std::vector<double> A2_snap(T_int, 0.0);
#elif defined(PROCESS_BMLONG)
        std::vector<double> A1_snap(T_int, 0.0);
        std::vector<double> A2_snap(T_int, 0.0);
        long long local_clamps = 0;
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
            // A1x = int_0^t F_x(Y_s) ds  (Gaussian force, mu=1)
            // A2x = -int_0^t [grad F(Y_u)^T Avec(u)]_x du  (Picard correction, mu=1)
            // Gaussian force gradient:
            //   d_x F_x = force_pref * exp(-r^2/2s^2) * (1 - x1^2/s^2)
            //   d_y F_x = force_pref * exp(-r^2/2s^2) * (-x1*x2/s^2)
            // ORDERING: A2 updated with pre-update Avec (s < u time ordering).
            double Avec_x = 0.0, Avec_y = 0.0, A2x = 0.0;
            int next_snap = 0;
            for (int s = 0; s < n_steps; ++s) {
                double x = state.x, y = state.y;
                double r2  = x*x + y*y;
                double fac = force_pref * std::exp(-half_inv_s2 * r2);
                double Fx  = fac * x;
                double Fy  = fac * y;
                double dxFx = fac * (1.0 - x*x * inv_s2);
                double dyFx = fac * (-x*y * inv_s2);

                // Step 1: update A2 using current Avec (s < u ordering)
                A2x -= (dxFx * Avec_x + dyFx * Avec_y) * dt;

                // Step 2: advance Avec
                Avec_x += Fx * dt;
                Avec_y += Fy * dt;

                // Step 3: advance swimmer
                swimmer->step(state, rng);

                double t_now = (s + 1) * dt;
                while (next_snap < T_int && (next_snap + 1) <= t_now + 1e-12) {
                    A1_snap[next_snap] = Avec_x;   // mu=1 implicit
                    A2_snap[next_snap] = A2x;
                    ++next_snap;
                }
            }
            while (next_snap < T_int) {
                A1_snap[next_snap] = Avec_x;
                A2_snap[next_snap] = A2x;
                ++next_snap;
            }
            A = Avec_x;  // final A1 for main output

            double w = b * b;
            for (int k = 0; k < T_int; ++k) {
                double a1   = A1_snap[k];
                double a1sq = a1 * a1;
                tl_w    [tid][k] += w;
                tl_wA1  [tid][k] += w * a1;
                tl_wA1sq[tid][k] += w * a1sq;
                tl_wA1_4[tid][k] += w * a1sq * a1sq;
                tl_wA2  [tid][k] += w * A2_snap[k];
                tl_wA2sq[tid][k] += w * A2_snap[k] * A2_snap[k];
            }

#elif defined(PROCESS_BMLONG)
            // A1x(t) = mu * Avec_x(t)  where Avec_x = int_0^t Fx ds
            // A2x(t) = -mu^2 * int_0^t [grad F(Y_u)^T Avec(u)]_x du
            //
            // ORDERING: A2 is updated with the CURRENT Avec (memory up to s < u),
            // then Avec is advanced. This respects the s < u time-ordering of the
            // nested integral (inner integral excludes the current instant).
            double Avec_x = 0.0, Avec_y = 0.0, A2x = 0.0;
            int next_snap = 0;
            for (int s = 0; s < n_steps; ++s) {
                auto ce = coulomb_eval(state.x, state.y, sigma_c, b_min_c);
                if (ce.clamped) ++local_clamps;

                // Step 1: update A2 using Avec BEFORE advancing Avec (s < u ordering)
                double gx = ce.dxFx * Avec_x + ce.dyFx * Avec_y;
                A2x -= mu * mu * gx * dt;

                // Step 2: advance Avec (memory for future steps)
                Avec_x += ce.Fx * dt;
                Avec_y += ce.Fy * dt;

                // Step 3: advance swimmer
                swimmer->step(state, rng);

                // Checkpoint: record A1x = mu*Avec_x and A2x at integer times
                double t_now = (s + 1) * dt;
                while (next_snap < T_int && (next_snap + 1) <= t_now + 1e-12) {
                    A1_snap[next_snap] = mu * Avec_x;
                    A2_snap[next_snap] = A2x;
                    ++next_snap;
                }
            }
            while (next_snap < T_int) {
                A1_snap[next_snap] = mu * Avec_x;
                A2_snap[next_snap] = A2x;
                ++next_snap;
            }
            A = mu * Avec_x;  // final A1x for main output

            double w = b * b;
            for (int k = 0; k < T_int; ++k) {
                double a1   = A1_snap[k];
                double a1sq = a1 * a1;
                tl_w    [tid][k] += w;
                tl_wA1  [tid][k] += w * a1;
                tl_wA1sq[tid][k] += w * a1sq;
                tl_wA1_4[tid][k] += w * a1sq * a1sq;
                tl_wA2  [tid][k] += w * A2_snap[k];
                tl_wA2sq[tid][k] += w * A2_snap[k] * A2_snap[k];
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
            w_vals[i] = b * b;
        }

#if defined(PROCESS_BMLONG)
        #pragma omp atomic
        total_clamps += local_clamps;
#endif
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

    // Helper: reduce thread-local accumulators for one checkpoint k
    // Returns {sw, swA1, swA1sq, swA1_4} (bmshort) — defined inline below per process

#if defined(PROCESS_BMSHORT)
    // Write var_bmshort.csv: t, VarA1, VarA2, meanA1, meanA2, ratio, kurtA1
    // Column order matches var_bmlong.csv.
    {
        std::string dir = params.output.substr(0, params.output.rfind('/') + 1);
        std::string out_path = dir + "var_bmshort.csv";
        std::ofstream fout(out_path);
        if (!fout) { std::cerr << "Error: cannot open " << out_path << "\n"; return; }
        fout << "t,VarA1,VarA2,meanA1,meanA2,ratio,kurtA1\n";
        fout.precision(15);
        for (int k = 0; k < T_int; ++k) {
            double sw = 0.0, swA1 = 0.0, swA1sq = 0.0, swA1_4 = 0.0, swA2 = 0.0, swA2sq = 0.0;
            for (int t = 0; t < n_threads_used; ++t) {
                sw     += tl_w    [t][k];
                swA1   += tl_wA1  [t][k];
                swA1sq += tl_wA1sq[t][k];
                swA1_4 += tl_wA1_4[t][k];
                swA2   += tl_wA2  [t][k];
                swA2sq += tl_wA2sq[t][k];
            }
            double mA1  = swA1   / sw;
            double mA2  = swA2   / sw;
            double vA1  = swA1sq / sw - mA1 * mA1;
            double vA2  = swA2sq / sw - mA2 * mA2;
            double ratio = (vA1 > 0.0) ? vA2 / vA1 : 0.0;
            double kA1   = (vA1 > 0.0) ? (swA1_4 / sw) / (vA1 * vA1) : 0.0;
            fout << (k + 1) << "," << vA1 << "," << vA2 << ","
                 << mA1 << "," << mA2 << "," << ratio << "," << kA1 << "\n";
        }
        std::cout << "var_bmshort written to " << out_path << "\n";
    }

#elif defined(PROCESS_BMLONG)
    // Write var_bmlong.csv: t, VarA1, VarA2, meanA1, meanA2, ratio, kurtA1
    {
        std::string dir = params.output.substr(0, params.output.rfind('/') + 1);
        std::string out_path = dir + "var_bmlong.csv";
        std::ofstream fout(out_path);
        if (!fout) { std::cerr << "Error: cannot open " << out_path << "\n"; return; }
        fout << "t,VarA1,VarA2,meanA1,meanA2,ratio,kurtA1\n";
        fout.precision(15);
        double final_ratio = 0.0;
        for (int k = 0; k < T_int; ++k) {
            double sw = 0.0, swA1 = 0.0, swA1sq = 0.0, swA1_4 = 0.0, swA2 = 0.0, swA2sq = 0.0;
            for (int t = 0; t < n_threads_used; ++t) {
                sw     += tl_w    [t][k];
                swA1   += tl_wA1  [t][k];
                swA1sq += tl_wA1sq[t][k];
                swA1_4 += tl_wA1_4[t][k];
                swA2   += tl_wA2  [t][k];
                swA2sq += tl_wA2sq[t][k];
            }
            double mA1   = swA1   / sw;
            double mA2   = swA2   / sw;
            double vA1   = swA1sq / sw - mA1 * mA1;
            double vA2   = swA2sq / sw - mA2 * mA2;
            double ratio = (vA1 > 0.0) ? vA2 / vA1 : 0.0;
            double kA1   = (vA1 > 0.0) ? (swA1_4 / sw) / (vA1 * vA1) : 0.0;
            fout << (k + 1) << "," << vA1 << "," << vA2 << ","
                 << mA1 << "," << mA2 << "," << ratio << "," << kA1 << "\n";
            if (k == T_int - 1) final_ratio = ratio;
        }
        double clamp_rate = static_cast<double>(total_clamps) / (static_cast<double>(N) * n_steps);
        std::cout << "var_bmlong written to " << out_path << "\n";
        std::cout << "b_max used        = " << params.sampling.b_max
                  << "  (= 2*sqrt(4*D_bm*T) = " << 2.0 * std::sqrt(4.0 * params.process.D_bm * params.T) << ")\n";
        std::cout << "mu                = " << params.mu << "\n";
        std::cout << "Clamp rate        = " << clamp_rate << "  (fraction of steps where r < b_min)\n";
        std::cout << "Var(A2)/Var(A1) at t=T = " << final_ratio << "\n";
    }
#endif

    // Write config (only parameters relevant for this process and force)
    {
        std::string config_path = params.output.substr(0, params.output.rfind('/') + 1) + "config_" + process_name() + ".txt";
        std::ofstream cfg(config_path);
        cfg << "process = " << process_name() << "\n";

#if defined(PROCESS_AOUP)
        cfg << "tau_c   = " << params.process.tau_c << "\n"
            << "D_A     = " << params.process.D_A   << "\n"
            << "p       = " << params.force.p       << "\n"
            << "b_min   = " << params.force.b_min   << "\n";
#elif defined(PROCESS_ABP)
        cfg << "v_A     = " << params.process.v_A   << "\n"
            << "D_r     = " << params.process.D_r   << "\n"
            << "p       = " << params.force.p       << "\n"
            << "b_min   = " << params.force.b_min   << "\n";
#elif defined(PROCESS_RTP)
        cfg << "v_A     = " << params.process.v_A   << "\n"
            << "omega   = " << params.process.omega  << "\n"
            << "p       = " << params.force.p       << "\n"
            << "b_min   = " << params.force.b_min   << "\n";
#elif defined(PROCESS_LEVY1)
        cfg << "v_A     = " << params.process.v_A   << "\n"
            << "beta    = " << params.process.beta   << "\n"
            << "p       = " << params.force.p       << "\n"
            << "b_min   = " << params.force.b_min   << "\n";
#elif defined(PROCESS_LEVY2)
        cfg << "v_A     = " << params.process.v_A   << "\n"
            << "beta    = " << params.process.beta   << "\n"
            << "tau_0   = " << params.process.tau_0  << "\n"
            << "p       = " << params.force.p       << "\n"
            << "b_min   = " << params.force.b_min   << "\n";
#elif defined(PROCESS_BMSHORT)
        cfg << "D_bm    = " << params.process.D_bm   << "\n"
            << "V_0     = " << params.process.V_0    << "\n"
            << "sigma   = " << params.process.sigma   << "\n"
            << "b_min   = " << params.force.b_min    << "\n";
#elif defined(PROCESS_BMLONG)
        cfg << "D_bm    = " << params.process.D_bm   << "\n"
            << "sigma   = " << params.force.sigma     << "\n"
            << "b_min   = " << params.force.b_min     << "\n"
            << "mu      = " << params.mu              << "\n";
#endif
        cfg << "b_max   = " << params.sampling.b_max << "\n"
            << "gamma   = " << params.sampling.gamma << "\n"
            << "dt      = " << params.process.dt     << "\n"
            << "T       = " << params.T              << "\n"
            << "N_traj  = " << params.N_traj         << "\n"
            << "seed    = " << params.seed           << "\n";
        std::cout << "Config written to " << config_path << "\n";
    }

    // Write raw data
    std::ofstream out(params.output);
    if (!out) { std::cerr << "Error: cannot open " << params.output << "\n"; return; }
    out << "b,A,w\n";
    out.precision(15);
    for (int i = 0; i < N; ++i)
        out << b_vals[i] << "," << A_vals[i] << "," << w_vals[i] << "\n";
    std::cout << "Output written to " << params.output << "\n";
}
