// Tests for the FD-based A2 accumulator (process-agnostic finite-difference gradient).
//
// T1: FD gradient of Coulomb force agrees with analytic gradient to ~1e-4 relative.
// T2: A2 on-the-fly accumulator matches brute-force double loop (s<u ordering),
//     for BOTH Coulomb (no orientation) and Dipole (oriented) forces.
// T3: FD-based A2 for Coulomb agrees with analytic-gradient A2 (coulomb_eval) to FD accuracy.
// T4: Symmetry: <A1> ~ 0 and <A2> ~ 0 for a small dipole ensemble.
#include "force.h"
#include "bmshort_swimmer.h"
#include "abp_swimmer.h"
#include <cmath>
#include <cassert>
#include <iostream>
#include <vector>
#include <random>

static bool approx(double a, double b, double tol) { return std::abs(a - b) < tol; }

static bool rel_approx(double a, double b, double tol) {
    double den = 0.5 * (std::abs(a) + std::abs(b));
    if (den < 1e-30) return std::abs(a - b) < tol;
    return std::abs(a - b) / den < tol;
}

// ---- T1: FD gradient vs analytic Coulomb gradient ---------------------------
void test_fd_vs_analytic_coulomb() {
    const double sigma = 1.5, b_min = 0.5, h = 1e-4;
    ForceParams p;
    p.type  = ForceType::COULOMB;
    p.sigma = sigma;
    p.b_min = b_min;

    // Several non-clamped points
    std::vector<std::pair<double,double>> pts = {
        {1.0, 0.5}, {2.0, -1.0}, {0.8, 1.5}, {-1.2, 0.9}, {3.0, -2.0}
    };
    for (auto [x1, x2] : pts) {
        auto fg  = fd_grad_Fx(x1, x2, 0.0, 0.0, p, h);
        auto ce  = coulomb_eval(x1, x2, sigma, b_min);
        assert(rel_approx(fg.dxFx, ce.dxFx, 1e-4));
        assert(rel_approx(fg.dyFx, ce.dyFx, 1e-4));
    }
    std::cout << "PASS T1: FD gradient agrees with analytic Coulomb gradient to 1e-4 relative\n";
}

// Helper: compute A2 via on-the-fly accumulator on a prescribed path.
// Returns A2x at the end.
static double a2_otf(const std::vector<double>& xs, const std::vector<double>& ys,
                      const std::vector<double>& n1s, const std::vector<double>& n2s,
                      const ForceParams& p, double mu, double dt, double h) {
    double Avec_x = 0.0, Avec_y = 0.0, A2x = 0.0;
    int N = static_cast<int>(xs.size());
    for (int u = 0; u < N; ++u) {
        auto fv = force_2d(xs[u], ys[u], n1s[u], n2s[u], p);
        auto fg = fd_grad_Fx(xs[u], ys[u], n1s[u], n2s[u], p, h);
        double gx = fg.dxFx * Avec_x + fg.dyFx * Avec_y;
        A2x   -= mu * mu * gx * dt;
        Avec_x += fv[0] * dt;
        Avec_y += fv[1] * dt;
    }
    return A2x;
}

// Helper: compute A2 via brute-force double loop.
static double a2_bf(const std::vector<double>& xs, const std::vector<double>& ys,
                     const std::vector<double>& n1s, const std::vector<double>& n2s,
                     const ForceParams& p, double mu, double dt, double h) {
    int N = static_cast<int>(xs.size());
    double A2x = 0.0;
    for (int u = 0; u < N; ++u) {
        // Avec(u) = sum_{s < u} F(Y_s) * dt
        double Av_x = 0.0, Av_y = 0.0;
        for (int s = 0; s < u; ++s) {
            auto fv_s = force_2d(xs[s], ys[s], n1s[s], n2s[s], p);
            Av_x += fv_s[0] * dt;
            Av_y += fv_s[1] * dt;
        }
        auto fg_u = fd_grad_Fx(xs[u], ys[u], n1s[u], n2s[u], p, h);
        double gx = fg_u.dxFx * Av_x + fg_u.dyFx * Av_y;
        A2x -= mu * mu * gx * dt;
    }
    return A2x;
}

// ---- T2a: A2 ordering for Coulomb (no orientation) --------------------------
void test_a2_ordering_coulomb() {
    const double sigma = 1.0, b_min = 0.5, mu = 1.0, dt = 0.1, h = 1e-4;
    const int N = 25;
    ForceParams p;
    p.type  = ForceType::COULOMB;
    p.sigma = sigma;
    p.b_min = b_min;

    std::vector<double> xs(N), ys(N), n1s(N, 0.0), n2s(N, 0.0);
    for (int i = 0; i < N; ++i) {
        double angle = 0.3 * i;
        double r     = 1.0 + 0.1 * i;
        xs[i] = r * std::cos(angle);
        ys[i] = r * std::sin(angle);
    }

    double otf = a2_otf(xs, ys, n1s, n2s, p, mu, dt, h);
    double bf  = a2_bf (xs, ys, n1s, n2s, p, mu, dt, h);
    assert(approx(otf, bf, 1e-10));
    std::cout << "PASS T2a: A2 on-the-fly matches brute-force (Coulomb, s<u ordering)\n";
}

// ---- T2b: A2 ordering for Dipole (oriented force) ---------------------------
void test_a2_ordering_dipole() {
    const double p_str = 1.0, b_min = 0.5, mu = 1.0, dt = 0.1, h = 1e-4;
    const int N = 25;
    ForceParams p;
    p.type  = ForceType::DIPOLE;
    p.p     = p_str;
    p.b_min = b_min;

    std::vector<double> xs(N), ys(N), n1s(N), n2s(N);
    for (int i = 0; i < N; ++i) {
        double angle = 0.3 * i;
        double r     = 1.0 + 0.1 * i;
        xs[i] = r * std::cos(angle);
        ys[i] = r * std::sin(angle);
        // Rotate orientation separately to exercise non-trivial n dependence
        double phi = 0.5 * i;
        n1s[i] = std::cos(phi);
        n2s[i] = std::sin(phi);
    }

    double otf = a2_otf(xs, ys, n1s, n2s, p, mu, dt, h);
    double bf  = a2_bf (xs, ys, n1s, n2s, p, mu, dt, h);
    assert(approx(otf, bf, 1e-10));
    std::cout << "PASS T2b: A2 on-the-fly matches brute-force (Dipole+orientation, s<u ordering)\n";
}

// ---- T3: FD-based Coulomb A2 agrees with analytic-gradient Coulomb A2 -------
void test_fd_coulomb_vs_analytic_a2() {
    const double sigma = 1.0, b_min = 0.5, mu = 1.0, dt = 0.05, h = 1e-4;
    const int N = 30;
    ForceParams p;
    p.type  = ForceType::COULOMB;
    p.sigma = sigma;
    p.b_min = b_min;

    std::vector<double> xs(N), ys(N);
    for (int i = 0; i < N; ++i) {
        double angle = 0.2 * i;
        double r     = 1.2 + 0.05 * i;
        xs[i] = r * std::cos(angle);
        ys[i] = r * std::sin(angle);
    }

    // FD-based A2
    std::vector<double> n1s(N, 0.0), n2s(N, 0.0);
    double A2_fd = a2_otf(xs, ys, n1s, n2s, p, mu, dt, h);

    // Analytic-gradient A2 (using coulomb_eval)
    double Avec_x = 0.0, Avec_y = 0.0, A2_an = 0.0;
    for (int u = 0; u < N; ++u) {
        auto ce = coulomb_eval(xs[u], ys[u], sigma, b_min);
        double gx = ce.dxFx * Avec_x + ce.dyFx * Avec_y;
        A2_an  -= mu * mu * gx * dt;
        Avec_x += ce.Fx * dt;
        Avec_y += ce.Fy * dt;
    }

    // Should agree to FD accuracy (~h^2 * scale)
    assert(rel_approx(A2_fd, A2_an, 1e-3));
    std::cout << "PASS T3: FD-based Coulomb A2 agrees with analytic-gradient A2 (rel err < 1e-3)\n";
}

// ---- T4: Symmetry test (dipole, small ensemble) -----------------------------
void test_symmetry_dipole() {
    const double p_str = 1.0, b_min = 0.5, mu = 1.0, dt = 1e-3, h = 1e-4;
    const double v_A = 1.0, D_r = 1.0;
    const int n_steps = 500, n_traj = 3000;
    ForceParams fp;
    fp.type  = ForceType::DIPOLE;
    fp.p     = p_str;
    fp.b_min = b_min;

    std::mt19937_64 rng(12345);
    std::uniform_real_distribution<double> ud(0.0, 1.0);
    ABPSwimmer swimmer(v_A, D_r, dt);

    double sum_A1 = 0.0, sum_A2 = 0.0;
    for (int i = 0; i < n_traj; ++i) {
        double b     = b_min * std::exp(ud(rng) * std::log(10.0 / b_min));
        double theta = ud(rng) * (2.0 * M_PI);
        auto state   = swimmer.init(b * std::cos(theta), b * std::sin(theta), rng);

        double Avec_x = 0.0, Avec_y = 0.0, A2x = 0.0;
        for (int s = 0; s < n_steps; ++s) {
            auto nv = swimmer.orientation(state);
            auto fv = force_2d(state.x, state.y, nv[0], nv[1], fp);
            auto fg = fd_grad_Fx(state.x, state.y, nv[0], nv[1], fp, h);
            A2x   -= mu * mu * (fg.dxFx * Avec_x + fg.dyFx * Avec_y) * dt;
            Avec_x += fv[0] * dt;
            Avec_y += fv[1] * dt;
            swimmer.step(state, rng);
        }
        sum_A1 += mu * Avec_x;
        sum_A2 += A2x;
    }
    double mean_A1 = sum_A1 / n_traj;
    double mean_A2 = sum_A2 / n_traj;
    assert(std::abs(mean_A1) < 1.0);
    assert(std::abs(mean_A2) < 1.0);
    std::cout << "PASS T4: Dipole <A1>~0 (" << mean_A1 << ") <A2>~0 (" << mean_A2 << ") by symmetry\n";
}

int main() {
    test_fd_vs_analytic_coulomb();
    test_a2_ordering_coulomb();
    test_a2_ordering_dipole();
    test_fd_coulomb_vs_analytic_a2();
    test_symmetry_dipole();
    std::cout << "All fd_a2 tests passed.\n";
    return 0;
}
