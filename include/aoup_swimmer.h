#pragma once
#include "swimmer.h"
#include <cmath>

class AOUPSwimmer : public Swimmer {
public:
    AOUPSwimmer(double tau_c, double D_A, double dt)
        : tau_c_(tau_c), D_A_(D_A), dt_(dt),
          decay_(1.0 - dt / tau_c),
          noise_std_(std::sqrt(2.0 * D_A / tau_c * dt)),
          sigma_v_(std::sqrt(D_A / tau_c)) {}

    SwimmerState init(double x0, double y0, std::mt19937_64& rng) const override {
        std::normal_distribution<double> nd(0.0, sigma_v_);
        return {x0, y0, nd(rng), nd(rng), 0.0};
    }

    void step(SwimmerState& s, std::mt19937_64& rng) const override {
        std::normal_distribution<double> nd(0.0, 1.0);
        s.x  += s.vx * dt_;
        s.y  += s.vy * dt_;
        s.vx  = s.vx * decay_ + noise_std_ * nd(rng);
        s.vy  = s.vy * decay_ + noise_std_ * nd(rng);
    }

    std::array<double, 2> orientation(const SwimmerState& s) const override {
        double mag = std::sqrt(s.vx * s.vx + s.vy * s.vy);
        if (mag < 1e-15) return {1.0, 0.0};
        return {s.vx / mag, s.vy / mag};
    }

private:
    double tau_c_, D_A_, dt_, decay_, noise_std_, sigma_v_;
};
