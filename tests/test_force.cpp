// Basic unit tests for force evaluation
#include "force.h"
#include <cassert>
#include <cmath>
#include <iostream>

void test_hard_core() {
    ForceParams p = {1.0, 0.5};
    double f = force_x(0.1, 0.0, 1.0, 0.0, p);
    assert(f == 0.0);
    std::cout << "PASS: hard core cutoff\n";
}

void test_symmetry() {
    // f_x should be antisymmetric: f_x(x, n) = -f_x(-x, n) when n is fixed
    ForceParams p = {1.0, 0.1};
    double f1 = force_x( 2.0, 1.0, 1.0, 0.0, p);
    double f2 = force_x(-2.0,-1.0, 1.0, 0.0, p);
    assert(std::abs(f1 + f2) < 1e-12);
    std::cout << "PASS: antisymmetry under x -> -x\n";
}

void test_known_value() {
    // On-axis: x = (r, 0), n = (1, 0): cos^2 = 1, f_x = (p/r^3)*(3-1)*r = 2p/r^2
    ForceParams p = {1.0, 0.1};
    double r = 3.0;
    double f = force_x(r, 0.0, 1.0, 0.0, p);
    double expected = 2.0 * p.p / (r * r);
    assert(std::abs(f - expected) < 1e-12);
    std::cout << "PASS: on-axis value\n";
}

void test_perpendicular() {
    // Perpendicular: x = (r, 0), n = (0, 1): cos^2 = 0, f_x = (p/r^3)*(-1)*r = -p/r^2
    ForceParams p = {1.0, 0.1};
    double r = 3.0;
    double f = force_x(r, 0.0, 0.0, 1.0, p);
    double expected = -p.p / (r * r);
    assert(std::abs(f - expected) < 1e-12);
    std::cout << "PASS: perpendicular value\n";
}

int main() {
    test_hard_core();
    test_symmetry();
    test_known_value();
    test_perpendicular();
    std::cout << "All force tests passed.\n";
    return 0;
}
