#include "swimmer_factory.h"
#include "aoup_swimmer.h"
#include "abp_swimmer.h"
#include "rtp_swimmer.h"
#include "levy_swimmer.h"
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>

static bool approx(double a, double b, double tol) { return std::abs(a - b) < tol; }

static void check(bool cond, const char* msg) {
    if (!cond) { std::cerr << "FAIL: " << msg << "\n"; std::exit(1); }
    std::cout << "PASS: " << msg << "\n";
}

int main() {
    std::mt19937_64 rng(12345);
    const int N = 100000;
    const double dt = 1e-3;

    // ---- AOUP: stationary velocity variance = D_A/tau_c per component ----
    {
        AOUPSwimmer sw(1.0, 1.0, dt);
        double var = 0.0;
        for (int i = 0; i < N; ++i) {
            auto s = sw.init(0.0, 0.0, rng);
            var += s.vx * s.vx + s.vy * s.vy;
        }
        var /= (2.0 * N);  // two components, D_A/tau_c = 1
        check(approx(var, 1.0, 0.05), "AOUP stationary velocity variance = D_A/tau_c");
    }

    // ---- ABP: |velocity| = v_A at all times ----
    {
        ABPSwimmer sw(1.0, 1.0, dt);
        auto state = sw.init(0.0, 0.0, rng);
        bool ok = true;
        for (int i = 0; i < 10000; ++i) {
            double speed = std::sqrt(state.vx*state.vx + state.vy*state.vy);
            if (!approx(speed, 1.0, 1e-12)) { ok = false; break; }
            sw.step(state, rng);
        }
        check(ok, "ABP |velocity| = v_A at all times");
    }

    // ---- RTP: |velocity| = v_A at all times ----
    {
        RTPSwimmer sw(1.0, 1.0, dt);
        auto state = sw.init(0.0, 0.0, rng);
        bool ok = true;
        for (int i = 0; i < 10000; ++i) {
            double speed = std::sqrt(state.vx*state.vx + state.vy*state.vy);
            if (!approx(speed, 1.0, 1e-12)) { ok = false; break; }
            sw.step(state, rng);
        }
        check(ok, "RTP |velocity| = v_A at all times");
    }

    // ---- RTP: tumble rate matches omega ----
    {
        const double omega = 1.0;
        RTPSwimmer sw(1.0, omega, dt);
        auto state = sw.init(0.0, 0.0, rng);
        int tumbles = 0;
        double prev_theta = std::atan2(state.vy, state.vx);
        const int steps = 1000000;
        for (int i = 0; i < steps; ++i) {
            sw.step(state, rng);
            double theta = std::atan2(state.vy, state.vx);
            if (std::abs(theta - prev_theta) > 1e-10) ++tumbles;
            prev_theta = theta;
        }
        double measured_rate = (double)tumbles / (steps * dt);
        check(approx(measured_rate, omega, 0.05), "RTP tumble rate matches omega");
    }

    // ---- Levy: run lengths follow power-law (CCDF slope ~ -beta) ----
    {
        const double beta = 1.5;
        LevySwimmer sw(1.0, beta, dt);
        // collect run lengths by simulating many runs
        std::vector<double> run_lengths;
        run_lengths.reserve(5000);
        auto state = sw.init(0.0, 0.0, rng);
        double run = state.aux;
        for (int i = 0; i < 5000000 && run_lengths.size() < 5000; ++i) {
            double prev_vx = state.vx;
            sw.step(state, rng);
            if (std::abs(state.vx - prev_vx) > 1e-10) {
                if (run > 0) run_lengths.push_back(run);
                run = state.aux;
            }
        }
        // Estimate slope of log CCDF at upper tail
        std::sort(run_lengths.begin(), run_lengths.end());
        int n = (int)run_lengths.size();
        // use top 20% of runs for slope estimate
        int lo = 4 * n / 5;
        double log_tau1 = std::log(run_lengths[lo]);
        double log_tau2 = std::log(run_lengths[n - 1]);
        double log_ccdf1 = std::log((double)(n - lo) / n);
        double log_ccdf2 = std::log(1.0 / n);
        double slope = (log_ccdf2 - log_ccdf1) / (log_tau2 - log_tau1);
        check(approx(slope, -beta, 0.3), "Levy run lengths follow power-law P(tau)~tau^{-(1+beta)}");
    }

    // ---- All processes: orientation is a unit vector ----
    {
        auto test_unit = [&](Swimmer& sw, const char* name) {
            auto state = sw.init(1.0, 0.5, rng);
            bool ok = true;
            for (int i = 0; i < 1000; ++i) {
                auto n = sw.orientation(state);
                double mag = std::sqrt(n[0]*n[0] + n[1]*n[1]);
                if (!approx(mag, 1.0, 1e-12)) { ok = false; break; }
                sw.step(state, rng);
            }
            std::string msg = std::string(name) + " orientation is always a unit vector";
            check(ok, msg.c_str());
        };
        AOUPSwimmer aoup(1.0, 1.0, dt); test_unit(aoup, "AOUP");
        ABPSwimmer  abp(1.0, 1.0, dt);  test_unit(abp,  "ABP");
        RTPSwimmer  rtp(1.0, 1.0, dt);  test_unit(rtp,  "RTP");
        LevySwimmer levy(1.0, 1.5, dt); test_unit(levy, "Levy");
    }

    // ---- All processes: mean displacement ~ 0 (isotropy) ----
    {
        auto test_isotropy = [&](Swimmer& sw, const char* name) {
            double mx = 0.0, my = 0.0;
            const int steps = 10000;
            for (int traj = 0; traj < 200; ++traj) {
                auto state = sw.init(0.0, 0.0, rng);
                for (int i = 0; i < steps; ++i) sw.step(state, rng);
                mx += state.x;
                my += state.y;
            }
            mx /= 200; my /= 200;
            // rough check: mean displacement much less than sqrt(D*T) ~ O(100)
            bool ok = std::abs(mx) < 20.0 && std::abs(my) < 20.0;
            std::string msg = std::string(name) + " mean displacement ~ 0 (isotropy)";
            check(ok, msg.c_str());
        };
        AOUPSwimmer aoup(1.0, 1.0, dt); test_isotropy(aoup, "AOUP");
        ABPSwimmer  abp(1.0, 1.0, dt);  test_isotropy(abp,  "ABP");
        RTPSwimmer  rtp(1.0, 1.0, dt);  test_isotropy(rtp,  "RTP");
        LevySwimmer levy(1.0, 1.5, dt); test_isotropy(levy, "Levy");
    }

    std::cout << "All swimmer tests passed.\n";
    return 0;
}
