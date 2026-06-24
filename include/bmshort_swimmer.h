#pragma once
#include "swimmer.h"
#include <cmath>
#include <stdexcept>

class BMShortSwimmer : public Swimmer {
public:
    BMShortSwimmer(double D, double dt) : D_(D), dt_(dt), noise_std_(std::sqrt(2.0 * D * dt)) {
        if (D <= 0.0)
            throw std::invalid_argument("BMShortSwimmer requires D > 0");
    }

    SwimmerState init(double x0, double y0, std::mt19937_64&) const override {
        return {x0, y0, 0.0, 0.0, 0.0};  // vx=vy=0 unused
    }

    void step(SwimmerState& s, std::mt19937_64& rng) const override {
        std::normal_distribution<double> nd(0.0, 1.0);
        s.x += noise_std_ * nd(rng);
        s.y += noise_std_ * nd(rng);
    }

    // Gaussian force is isotropic; orientation is unused — return dummy value
    std::array<double, 2> orientation(const SwimmerState&) const override {
        return {1.0, 0.0};
    }

private:
    double D_, dt_, noise_std_;
};
