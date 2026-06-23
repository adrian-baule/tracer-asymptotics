#pragma once
#include "swimmer.h"
#include <cmath>

class RTPSwimmer : public Swimmer {
public:
    RTPSwimmer(double v_A, double omega, double dt)
        : v_A_(v_A), omega_(omega), dt_(dt), tumble_prob_(omega * dt) {}

    SwimmerState init(double x0, double y0, std::mt19937_64& rng) const override {
        std::uniform_real_distribution<double> ud(0.0, 2.0 * M_PI);
        double theta = ud(rng);
        return {x0, y0, v_A_ * std::cos(theta), v_A_ * std::sin(theta), 0.0};
    }

    void step(SwimmerState& s, std::mt19937_64& rng) const override {
        std::uniform_real_distribution<double> ud(0.0, 1.0);
        if (ud(rng) < tumble_prob_) {
            double theta = ud(rng) * (2.0 * M_PI);
            s.vx = v_A_ * std::cos(theta);
            s.vy = v_A_ * std::sin(theta);
        }
        s.x += s.vx * dt_;
        s.y += s.vy * dt_;
    }

    std::array<double, 2> orientation(const SwimmerState& s) const override {
        return {s.vx / v_A_, s.vy / v_A_};
    }

private:
    double v_A_, omega_, dt_, tumble_prob_;
};
