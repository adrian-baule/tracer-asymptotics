#pragma once
#include "swimmer.h"
#include <cmath>

class ABPSwimmer : public Swimmer {
public:
    ABPSwimmer(double v_A, double D_r, double dt)
        : v_A_(v_A), D_r_(D_r), dt_(dt),
          angle_std_(std::sqrt(2.0 * D_r * dt)) {}

    SwimmerState init(double x0, double y0, std::mt19937_64& rng) const override {
        std::uniform_real_distribution<double> ud(0.0, 2.0 * M_PI);
        double theta = ud(rng);
        return {x0, y0, v_A_ * std::cos(theta), v_A_ * std::sin(theta), 0.0};
    }

    void step(SwimmerState& s, std::mt19937_64& rng) const override {
        std::normal_distribution<double> nd(0.0, 1.0);
        double theta = std::atan2(s.vy, s.vx) + angle_std_ * nd(rng);
        s.x  += s.vx * dt_;
        s.y  += s.vy * dt_;
        s.vx  = v_A_ * std::cos(theta);
        s.vy  = v_A_ * std::sin(theta);
    }

    std::array<double, 2> orientation(const SwimmerState& s) const override {
        return {s.vx / v_A_, s.vy / v_A_};
    }

private:
    double v_A_, D_r_, dt_, angle_std_;
};
