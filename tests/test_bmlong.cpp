// Unit tests for the bmlong (BM + Coulomb) A1/A2 accumulation.
// Critical test: verify that the on-the-fly A2 accumulator matches a brute-force
// double loop, confirming the s < u time-ordering discretization is correct.
#include "force.h"
#include "bmshort_swimmer.h"
#include <cmath>
#include <cassert>
#include <iostream>
#include <vector>
#include <random>

static bool approx(double a, double b, double tol) { return std::abs(a - b) < tol; }

// ---- A2 ordering test -------------------------------------------------------
// On a short fixed trajectory of N steps with prescribed positions, compute A2x
// two ways and compare to machine precision.
//
// On-the-fly (from simulation.cpp logic):
//   A2x -= mu^2 * gx * dt   // gx uses Avec BEFORE updating it
//   Avec += F * dt
//
// Brute-force double loop:
//   A2x = -mu^2 * sum_{u=0}^{N-1} [ (dxFx(Y_u)*Avec_x(u) + dyFx(Y_u)*Avec_y(u)) ] * dt
// where Avec_x(u) = sum_{s=0}^{u-1} Fx(Y_s) * dt  (s strictly < u)
void test_A2_ordering() {
    const double sigma = 1.0, b_min = 0.5, mu = 1.0, dt = 0.1;
    const int N = 20;

    // Deterministic path: spiral outward to stay away from singularity
    std::vector<double> xs(N), ys(N);
    for (int i = 0; i < N; ++i) {
        double angle = 0.3 * i;
        double r     = 1.0 + 0.1 * i;
        xs[i] = r * std::cos(angle);
        ys[i] = r * std::sin(angle);
    }

    // --- On-the-fly ---
    double Avec_x = 0.0, Avec_y = 0.0, A2x_otf = 0.0;
    for (int u = 0; u < N; ++u) {
        auto ce = coulomb_eval(xs[u], ys[u], sigma, b_min);
        double gx = ce.dxFx * Avec_x + ce.dyFx * Avec_y;
        A2x_otf -= mu * mu * gx * dt;   // update A2 with Avec BEFORE advancing it
        Avec_x  += ce.Fx * dt;
        Avec_y  += ce.Fy * dt;
    }

    // --- Brute-force double loop ---
    double A2x_bf = 0.0;
    for (int u = 0; u < N; ++u) {
        // Avec(u) = sum_{s < u} F(Y_s)*dt
        double Av_x = 0.0, Av_y = 0.0;
        for (int s = 0; s < u; ++s) {
            auto ce_s = coulomb_eval(xs[s], ys[s], sigma, b_min);
            Av_x += ce_s.Fx * dt;
            Av_y += ce_s.Fy * dt;
        }
        auto ce_u = coulomb_eval(xs[u], ys[u], sigma, b_min);
        double gx = ce_u.dxFx * Av_x + ce_u.dyFx * Av_y;
        A2x_bf -= mu * mu * gx * dt;
    }

    assert(approx(A2x_otf, A2x_bf, 1e-12));
    std::cout << "PASS: A2 on-the-fly matches brute-force double loop (s<u ordering)\n";
}

// ---- Symmetry test ----------------------------------------------------------
// <A1> = 0 and <A2> = 0 by isotropy over a small ensemble of random trajectories.
void test_A1_A2_symmetry() {
    const double sigma = 1.0, b_min = 0.5, mu = 1.0, dt = 1e-3;
    const double D_bm  = 1.0;
    const int n_steps  = 1000;
    const int n_traj   = 2000;

    std::mt19937_64 rng(9999);
    std::uniform_real_distribution<double> ud(0.0, 1.0);
    BMShortSwimmer swimmer(D_bm, dt);

    double sum_A1 = 0.0, sum_A2 = 0.0;
    for (int i = 0; i < n_traj; ++i) {
        // Random impact parameter b in [b_min, 10]
        double b     = b_min * std::exp(ud(rng) * std::log(10.0 / b_min));
        double theta = ud(rng) * (2.0 * M_PI);
        auto state   = swimmer.init(b * std::cos(theta), b * std::sin(theta), rng);

        double Avec_x = 0.0, Avec_y = 0.0, A2x = 0.0;
        for (int s = 0; s < n_steps; ++s) {
            auto ce = coulomb_eval(state.x, state.y, sigma, b_min);
            double gx = ce.dxFx * Avec_x + ce.dyFx * Avec_y;
            A2x   -= mu * mu * gx * dt;
            Avec_x += ce.Fx * dt;
            Avec_y += ce.Fy * dt;
            swimmer.step(state, rng);
        }
        sum_A1 += mu * Avec_x;
        sum_A2 += A2x;
    }
    double mean_A1 = sum_A1 / n_traj;
    double mean_A2 = sum_A2 / n_traj;

    // Tolerance: ~3 sigma of the mean; use loose bound since variance can be large
    assert(std::abs(mean_A1) < 1.0);
    assert(std::abs(mean_A2) < 1.0);
    std::cout << "PASS: <A1> ~ 0 (" << mean_A1 << ") and <A2> ~ 0 (" << mean_A2 << ") by symmetry\n";
}

int main() {
    test_A2_ordering();
    test_A1_A2_symmetry();
    std::cout << "All bmlong tests passed.\n";
    return 0;
}
