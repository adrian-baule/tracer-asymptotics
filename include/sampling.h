#pragma once
#include <random>

struct SamplingParams {
    double b_min;   // minimum impact parameter
    double b_max;   // maximum impact parameter
    double gamma;   // importance sampling exponent: p_bias(b) ~ b^(1+gamma)
};

// Sample impact parameter b from biased distribution p_bias(b) ~ b^(1+gamma)
// For gamma = -2: log-uniform sampling, b = b_min * (b_max/b_min)^U
double sample_b(const SamplingParams& params, std::mt19937_64& rng);

// Sample angle theta uniformly in [0, 2*pi)
double sample_theta(std::mt19937_64& rng);

// Importance weight w(b) = p_true(b) / p_bias(b), unnormalised
// p_true(b) ~ b (uniform 2D density), p_bias(b) ~ b^(1+gamma)
// w(b) ~ b^(-gamma)
double importance_weight(double b, const SamplingParams& params);

// Compute effective sample size N_eff = (sum w)^2 / sum(w^2)
double effective_sample_size(const double* weights, int N);
