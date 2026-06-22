#include "sampling.h"
#include <cmath>

double sample_b(const SamplingParams& params, std::mt19937_64& rng) {
    std::uniform_real_distribution<double> ud(0.0, 1.0);
    double u = ud(rng);
    if (std::abs(params.gamma + 2.0) < 1e-10) {
        return params.b_min * std::pow(params.b_max / params.b_min, u);
    }
    double alpha  = 2.0 + params.gamma;
    double bmin_a = std::pow(params.b_min, alpha);
    double bmax_a = std::pow(params.b_max, alpha);
    return std::pow(bmin_a + u * (bmax_a - bmin_a), 1.0 / alpha);
}

double sample_theta(std::mt19937_64& rng) {
    std::uniform_real_distribution<double> ud(0.0, 2.0 * M_PI);
    return ud(rng);
}

double importance_weight(double b, const SamplingParams& params) {
    return std::pow(b, -params.gamma);
}

double effective_sample_size(const double* weights, int N) {
    double sum_w = 0.0, sum_w2 = 0.0;
    for (int i = 0; i < N; ++i) { sum_w += weights[i]; sum_w2 += weights[i]*weights[i]; }
    return (sum_w * sum_w) / sum_w2;
}
