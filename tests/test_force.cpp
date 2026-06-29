// Basic unit tests for force evaluation
#include "force.h"
#include <cassert>
#include <cmath>
#include <iostream>

static ForceParams dipole(double p_val, double b_min) {
    ForceParams p; p.type = ForceType::DIPOLE; p.p = p_val; p.b_min = b_min; return p;
}

static ForceParams gaussian(double V_0, double sigma) {
    ForceParams p; p.type = ForceType::GAUSSIAN; p.V_0 = V_0; p.sigma = sigma; return p;
}

void test_hard_core() {
    auto p = dipole(1.0, 0.5);
    double f = force_x(0.1, 0.0, 1.0, 0.0, p);
    assert(f == 0.0);
    std::cout << "PASS: hard core cutoff\n";
}

void test_symmetry() {
    auto p = dipole(1.0, 0.1);
    double f1 = force_x( 2.0, 1.0, 1.0, 0.0, p);
    double f2 = force_x(-2.0,-1.0, 1.0, 0.0, p);
    assert(std::abs(f1 + f2) < 1e-12);
    std::cout << "PASS: antisymmetry under x -> -x\n";
}

void test_known_value() {
    // On-axis: x=(r,0), n=(1,0): cos^2=1, f_x = (p/r^3)*(3-1)*r = 2p/r^2
    auto p = dipole(1.0, 0.1);
    double r = 3.0;
    double f = force_x(r, 0.0, 1.0, 0.0, p);
    double expected = 2.0 * p.p / (r * r);
    assert(std::abs(f - expected) < 1e-12);
    std::cout << "PASS: on-axis value\n";
}

void test_perpendicular() {
    // Perpendicular: x=(r,0), n=(0,1): cos^2=0, f_x = (p/r^3)*(-1)*r = -p/r^2
    auto p = dipole(1.0, 0.1);
    double r = 3.0;
    double f = force_x(r, 0.0, 0.0, 1.0, p);
    double expected = -p.p / (r * r);
    assert(std::abs(f - expected) < 1e-12);
    std::cout << "PASS: perpendicular value\n";
}

void test_gaussian_symmetry() {
    // f_x(x1,x2) = -f_x(-x1,-x2)
    auto p = gaussian(1.0, 1.0);
    double f1 = force_x( 1.5, 0.7, 0.0, 0.0, p);
    double f2 = force_x(-1.5,-0.7, 0.0, 0.0, p);
    assert(std::abs(f1 + f2) < 1e-14);
    std::cout << "PASS: Gaussian force antisymmetry\n";
}

void test_gaussian_maximum() {
    // f_x(x1,0) is maximised at x1 = sigma (peak of x*exp(-x^2/(2s^2)))
    auto p = gaussian(1.0, 1.5);
    double sigma = p.sigma;
    double f_peak     = force_x(sigma,        0.0, 0.0, 0.0, p);
    double f_below    = force_x(sigma * 0.9,  0.0, 0.0, 0.0, p);
    double f_above    = force_x(sigma * 1.1,  0.0, 0.0, 0.0, p);
    assert(f_peak > f_below && f_peak > f_above);
    std::cout << "PASS: Gaussian force maximum at r = sigma\n";
}

void test_gaussian_decay() {
    // Force at 5*sigma should be negligible vs at sigma
    auto p = gaussian(1.0, 1.0);
    double f_sigma = force_x(p.sigma,       0.0, 0.0, 0.0, p);
    double f_far   = force_x(5.0 * p.sigma, 0.0, 0.0, 0.0, p);
    assert(std::abs(f_far / f_sigma) < 1e-6);
    std::cout << "PASS: Gaussian force negligible at 5*sigma\n";
}

void test_gaussian_orientation_invariant() {
    // Gaussian force must not depend on n
    auto p = gaussian(1.0, 1.0);
    double f1 = force_x(1.0, 0.5, 1.0, 0.0, p);
    double f2 = force_x(1.0, 0.5, 0.0, 1.0, p);
    double f3 = force_x(1.0, 0.5, 0.7, 0.7, p);
    assert(std::abs(f1 - f2) < 1e-14 && std::abs(f1 - f3) < 1e-14);
    std::cout << "PASS: Gaussian force independent of orientation\n";
}

void test_coulomb_x_value() {
    // f_x = sigma * x1 / r^3; on-axis (x2=0): f_x = sigma / x1^2
    double sigma = 2.0, b_min = 0.1, x1 = 3.0;
    auto e = coulomb_eval(x1, 0.0, sigma, b_min);
    double expected = sigma / (x1 * x1);
    assert(std::abs(e.Fx - expected) < 1e-12);
    std::cout << "PASS: Coulomb f_x magnitude\n";
}

void test_coulomb_antisymmetry() {
    // f_x(-x1,-x2) = -f_x(x1,x2)
    double sigma = 1.5, b_min = 0.1;
    auto e1 = coulomb_eval( 2.0,  1.0, sigma, b_min);
    auto e2 = coulomb_eval(-2.0, -1.0, sigma, b_min);
    assert(std::abs(e1.Fx + e2.Fx) < 1e-12);
    std::cout << "PASS: Coulomb antisymmetry\n";
}

void test_coulomb_gradient_fd() {
    // Finite-difference check of dxFx and dyFx at (x1,x2) = (2.0, 1.5)
    double sigma = 1.0, b_min = 0.01, h = 1e-6;
    double x1 = 2.0, x2 = 1.5;
    auto ep = coulomb_eval(x1 + h, x2, sigma, b_min);
    auto em = coulomb_eval(x1 - h, x2, sigma, b_min);
    double fd_dxFx = (ep.Fx - em.Fx) / (2.0 * h);
    auto e  = coulomb_eval(x1, x2, sigma, b_min);
    assert(std::abs(fd_dxFx - e.dxFx) < 1e-5);

    auto ep2 = coulomb_eval(x1, x2 + h, sigma, b_min);
    auto em2 = coulomb_eval(x1, x2 - h, sigma, b_min);
    double fd_dyFx = (ep2.Fx - em2.Fx) / (2.0 * h);
    assert(std::abs(fd_dyFx - e.dyFx) < 1e-5);
    std::cout << "PASS: Coulomb gradient finite-difference check\n";
}

void test_coulomb_clamp() {
    // At r < b_min the clamp flag is set and force is evaluated at b_min
    double sigma = 1.0, b_min = 0.5;
    auto e = coulomb_eval(0.1, 0.0, sigma, b_min);
    assert(e.clamped);
    // After clamping to (b_min, 0): Fx = sigma / b_min^2
    assert(std::abs(e.Fx - sigma / (b_min * b_min)) < 1e-12);
    std::cout << "PASS: Coulomb clamp at r < b_min\n";
}

int main() {
    test_hard_core();
    test_symmetry();
    test_known_value();
    test_perpendicular();
    test_gaussian_symmetry();
    test_gaussian_maximum();
    test_gaussian_decay();
    test_gaussian_orientation_invariant();
    test_coulomb_x_value();
    test_coulomb_antisymmetry();
    test_coulomb_gradient_fd();
    test_coulomb_clamp();
    std::cout << "All force tests passed.\n";
    return 0;
}
