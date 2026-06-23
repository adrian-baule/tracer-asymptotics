#pragma once
#include "swimmer.h"
#include <cmath>

class LevySwimmer : public Swimmer {
public:
    LevySwimmer(double v_A, double beta, double dt)
        : v_A_(v_A), beta_(beta), dt_(dt), inv_beta_(1.0 / beta) {}

    SwimmerState init(double x0, double y0, std::mt19937_64& rng) const override {
        std::uniform_real_distribution<double> ud(0.0, 1.0);
        double theta       = ud(rng) * (2.0 * M_PI);
        double remaining   = draw_run_time(ud(rng));
        return {x0, y0, v_A_ * std::cos(theta), v_A_ * std::sin(theta), remaining};
    }

    void step(SwimmerState& s, std::mt19937_64& rng) const override {
        if (s.aux <= 0.0) {
            std::uniform_real_distribution<double> ud(0.0, 1.0);
            double theta = ud(rng) * (2.0 * M_PI);
            s.vx  = v_A_ * std::cos(theta);
            s.vy  = v_A_ * std::sin(theta);
            s.aux = draw_run_time(ud(rng));
        }
        s.x   += s.vx * dt_;
        s.y   += s.vy * dt_;
        s.aux -= dt_;
    }

    std::array<double, 2> orientation(const SwimmerState& s) const override {
        return {s.vx / v_A_, s.vy / v_A_};
    }

private:
    double v_A_, beta_, dt_, inv_beta_;

    // P(tau) ~ tau^{-(1+beta)}: tau = dt * (U^{-1/beta} - 1)
    double draw_run_time(double u) const {
        return dt_ * (std::pow(u, -inv_beta_) - 1.0);
    }
};
