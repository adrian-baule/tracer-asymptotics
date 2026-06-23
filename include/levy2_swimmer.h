#pragma once
#include "swimmer.h"
#include <cmath>
#include <iostream>
#include <stdexcept>

class Levy2Swimmer : public Swimmer {
public:
    Levy2Swimmer(double v_A, double beta, double dt, double tau_0)
        : v_A_(v_A), beta_(beta), dt_(dt), tau_0_(tau_0), inv_beta_(1.0 / beta) {
        if (beta <= 0.0 || beta >= 1.0)
            throw std::invalid_argument("Levy2Swimmer requires 0 < beta < 1");
        if (tau_0 < 10.0 * dt)
            throw std::invalid_argument("Levy2Swimmer requires tau_0 >= 10*dt");
        if (tau_0 < 100.0 * dt)
            std::cerr << "Warning: tau_0/dt = " << tau_0 / dt
                      << ". Recommend tau_0 >> dt for levy2 (beta<1). "
                      << "The minimum run time tau_0 should be a physical parameter, "
                      << "not set by dt.\n";
    }

    SwimmerState init(double x0, double y0, std::mt19937_64& rng) const override {
        std::uniform_real_distribution<double> ud(0.0, 1.0);
        double theta = ud(rng) * (2.0 * M_PI);
        double aux   = draw_run_time(ud(rng));
        return {x0, y0, v_A_ * std::cos(theta), v_A_ * std::sin(theta), aux};
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
    double v_A_, beta_, dt_, tau_0_, inv_beta_;

    // P(tau > t) = (tau_0/t)^beta for t >= tau_0
    double draw_run_time(double u) const {
        return tau_0_ * std::pow(u, -inv_beta_);
    }
};
