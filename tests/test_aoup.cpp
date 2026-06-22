// Basic unit tests for AOUP integrator
#include "aoup.h"
#include <cassert>
#include <cmath>
#include <iostream>

void test_stationary_variance() {
    // After many steps, velocity variance should be D_A / tau_c per component
    AOUPParams p = {1.0, 1.0, 1e-3};
    std::mt19937_64 rng(12345);
    AOUPState state = aoup_init(0.0, 0.0, p, rng);
    int N = 100000;
    double sum_vx2 = 0.0;
    for (int i = 0; i < N; ++i) {
        aoup_step(state, p, rng);
        sum_vx2 += state.vx * state.vx;
    }
    double var_vx = sum_vx2 / N;
    double expected = p.D_A / p.tau_c;
    assert(std::abs(var_vx - expected) / expected < 0.05);
    std::cout << "PASS: stationary velocity variance = " << var_vx << " (expected " << expected << ")\n";
}

void test_orientation_unit() {
    AOUPParams p = {1.0, 1.0, 1e-3};
    std::mt19937_64 rng(42);
    AOUPState state = aoup_init(0.0, 0.0, p, rng);
    for (int i = 0; i < 1000; ++i) {
        aoup_step(state, p, rng);
        auto n = aoup_orientation(state);
        double mag = std::sqrt(n[0]*n[0] + n[1]*n[1]);
        assert(std::abs(mag - 1.0) < 1e-12);
    }
    std::cout << "PASS: orientation is always a unit vector\n";
}

int main() {
    test_stationary_variance();
    test_orientation_unit();
    std::cout << "All AOUP tests passed.\n";
    return 0;
}
