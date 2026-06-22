// Basic unit tests for importance sampling
#include "sampling.h"
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

void test_log_uniform() {
    SamplingParams p = {0.5, 50.0, -2.0};
    std::mt19937_64 rng(99);
    int N = 100000;
    int n_below = 0;
    double log_bmin = std::log(p.b_min);
    double log_bmax = std::log(p.b_max);
    double log_mid  = 0.5 * (log_bmin + log_bmax);
    for (int i = 0; i < N; ++i) {
        double b = sample_b(p, rng);
        assert(b >= p.b_min && b <= p.b_max);
        if (std::log(b) < log_mid) n_below++;
    }
    // Should be ~50% below geometric midpoint
    double frac = (double)n_below / N;
    assert(std::abs(frac - 0.5) < 0.01);
    std::cout << "PASS: log-uniform sampling, fraction below midpoint = " << frac << "\n";
}

void test_weight_normalisation() {
    SamplingParams p = {0.5, 50.0, -2.0};
    std::mt19937_64 rng(7);
    int N = 10000;
    std::vector<double> w(N);
    double sum = 0.0;
    for (int i = 0; i < N; ++i) {
        double b = sample_b(p, rng);
        w[i] = importance_weight(b, p);
        sum += w[i];
    }
    // After normalisation mean should be 1
    for (int i = 0; i < N; ++i) w[i] /= (sum / N);
    double mean = 0.0;
    for (int i = 0; i < N; ++i) mean += w[i];
    mean /= N;
    assert(std::abs(mean - 1.0) < 1e-10);
    std::cout << "PASS: normalised weights have mean = " << mean << "\n";
}

int main() {
    test_log_uniform();
    test_weight_normalisation();
    std::cout << "All sampling tests passed.\n";
    return 0;
}
