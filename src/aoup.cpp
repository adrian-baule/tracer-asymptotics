#include "aoup.h"
#include <cmath>

AOUPState aoup_init(double x0, double y0, const AOUPParams& params, std::mt19937_64& rng) {
    double sigma_v = std::sqrt(params.D_A / params.tau_c);
    std::normal_distribution<double> nd(0.0, sigma_v);
    return {x0, y0, nd(rng), nd(rng)};
}

void aoup_step(AOUPState& state, const AOUPParams& params, std::mt19937_64& rng) {
    double dt = params.dt;
    double decay = 1.0 - dt / params.tau_c;
    double noise_std = std::sqrt(2.0 * params.D_A / params.tau_c * dt);
    std::normal_distribution<double> nd(0.0, 1.0);
    state.x  += state.vx * dt;
    state.y  += state.vy * dt;
    state.vx  = state.vx * decay + noise_std * nd(rng);
    state.vy  = state.vy * decay + noise_std * nd(rng);
}

std::array<double, 2> aoup_orientation(const AOUPState& state) {
    double mag = std::sqrt(state.vx * state.vx + state.vy * state.vy);
    if (mag < 1e-15) return {1.0, 0.0};
    return {state.vx / mag, state.vy / mag};
}
