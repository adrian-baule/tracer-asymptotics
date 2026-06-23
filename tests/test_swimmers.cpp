#include "swimmer_factory.h"
#include "aoup_swimmer.h"
#include "abp_swimmer.h"
#include "rtp_swimmer.h"
#include "levy1_swimmer.h"
#include "levy2_swimmer.h"
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>

static bool approx(double a, double b, double tol) { return std::abs(a - b) < tol; }

static void check(bool cond, const char* msg) {
    if (!cond) { std::cerr << "FAIL: " << msg << "\n"; std::exit(1); }
    std::cout << "PASS: " << msg << "\n";
}

// Sample run times directly from init() — avoids waiting for long runs to complete
static std::vector<double> sample_run_times(Swimmer& sw, std::mt19937_64& rng, int n) {
    std::vector<double> times(n);
    for (int i = 0; i < n; ++i)
        times[i] = sw.init(0.0, 0.0, rng).aux;
    return times;
}

// Hill estimator for the tail index alpha where P(X > x) ~ x^{-alpha}
// Returns -alpha (negative, matching CCDF slope convention).
// Uses the top k order statistics; k ~ sqrt(n) is a common choice.
static double hill_slope(std::vector<double>& samples, int k = 0) {
    std::sort(samples.begin(), samples.end());
    int n = (int)samples.size();
    if (k <= 0) k = (int)std::sqrt((double)n);
    if (k >= n) k = n - 1;
    // X_{(1)} >= X_{(2)} >= ... so in sorted-ascending order X_{(i)} = samples[n-i]
    double threshold = samples[n - k - 1];
    double sum = 0.0;
    for (int i = 0; i < k; ++i)
        sum += std::log(samples[n - 1 - i] / threshold);
    double alpha_hat = k / sum;
    return -alpha_hat;  // CCDF slope = -alpha
}

int main() {
    std::mt19937_64 rng(12345);
    const int N  = 100000;
    const double dt = 1e-3;

    // ---- AOUP: stationary velocity variance = D_A/tau_c per component ----
    {
        AOUPSwimmer sw(1.0, 1.0, dt);
        double var = 0.0;
        for (int i = 0; i < N; ++i) {
            auto s = sw.init(0.0, 0.0, rng);
            var += s.vx * s.vx + s.vy * s.vy;
        }
        var /= (2.0 * N);
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

    // ---- Levy1 (beta=1.5): |velocity| = v_A at all times ----
    {
        Levy1Swimmer sw(1.0, 1.5, dt);
        auto state = sw.init(0.0, 0.0, rng);
        bool ok = true;
        for (int i = 0; i < 10000; ++i) {
            double speed = std::sqrt(state.vx*state.vx + state.vy*state.vy);
            if (!approx(speed, 1.0, 1e-12)) { ok = false; break; }
            sw.step(state, rng);
        }
        check(ok, "Levy1 |velocity| = v_A at all times");
    }

    // ---- Levy1 (beta=1.5): run times follow P(tau)~tau^{-(1+beta)} ----
    {
        const double beta = 1.5;
        Levy1Swimmer sw(1.0, beta, dt);
        auto runs = sample_run_times(sw, rng, 50000);
        double slope = hill_slope(runs);
        check(approx(slope, -beta, 0.15),
              "Levy1 run times: CCDF slope ~ -beta (beta=1.5)");
    }

    // ---- Levy1: construction throws for beta <= 1 ----
    {
        bool threw = false;
        try { Levy1Swimmer sw(1.0, 0.8, dt); }
        catch (const std::invalid_argument&) { threw = true; }
        check(threw, "Levy1 throws std::invalid_argument for beta <= 1");
    }

    // ---- Levy2 (beta=0.5): |velocity| = v_A at all times ----
    {
        Levy2Swimmer sw(1.0, 0.5, dt, 1.0);
        auto state = sw.init(0.0, 0.0, rng);
        bool ok = true;
        for (int i = 0; i < 10000; ++i) {
            double speed = std::sqrt(state.vx*state.vx + state.vy*state.vy);
            if (!approx(speed, 1.0, 1e-12)) { ok = false; break; }
            sw.step(state, rng);
        }
        check(ok, "Levy2 |velocity| = v_A at all times");
    }

    // ---- Levy2 (beta=0.5, tau_0=1.0): run times follow P(tau>t)=(tau_0/t)^beta ----
    {
        const double beta = 0.5, tau_0 = 1.0;
        Levy2Swimmer sw(1.0, beta, dt, tau_0);
        auto runs = sample_run_times(sw, rng, 50000);
        double slope = hill_slope(runs);
        check(approx(slope, -beta, 0.1),
              "Levy2 run times: CCDF slope ~ -beta (beta=0.5)");
    }

    // ---- Levy2: construction throws for beta >= 1 ----
    {
        bool threw = false;
        try { Levy2Swimmer sw(1.0, 1.5, dt, 1.0); }
        catch (const std::invalid_argument&) { threw = true; }
        check(threw, "Levy2 throws std::invalid_argument for beta >= 1");
    }

    // ---- Levy2: construction throws for tau_0 < 10*dt ----
    {
        bool threw = false;
        try { Levy2Swimmer sw(1.0, 0.5, dt, dt); }  // tau_0 = dt < 10*dt
        catch (const std::invalid_argument&) { threw = true; }
        check(threw, "Levy2 throws std::invalid_argument for tau_0 < 10*dt");
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
        AOUPSwimmer  aoup(1.0, 1.0, dt);       test_unit(aoup,  "AOUP");
        ABPSwimmer   abp(1.0, 1.0, dt);        test_unit(abp,   "ABP");
        RTPSwimmer   rtp(1.0, 1.0, dt);        test_unit(rtp,   "RTP");
        Levy1Swimmer lv1(1.0, 1.5, dt);        test_unit(lv1,   "Levy1");
        Levy2Swimmer lv2(1.0, 0.5, dt, 0.1);  test_unit(lv2,   "Levy2");
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
            bool ok = std::abs(mx) < 20.0 && std::abs(my) < 20.0;
            std::string msg = std::string(name) + " mean displacement ~ 0 (isotropy)";
            check(ok, msg.c_str());
        };
        AOUPSwimmer  aoup(1.0, 1.0, dt);       test_isotropy(aoup,  "AOUP");
        ABPSwimmer   abp(1.0, 1.0, dt);        test_isotropy(abp,   "ABP");
        RTPSwimmer   rtp(1.0, 1.0, dt);        test_isotropy(rtp,   "RTP");
        Levy1Swimmer lv1(1.0, 1.5, dt);        test_isotropy(lv1,   "Levy1");
        Levy2Swimmer lv2(1.0, 0.5, dt, 0.1);  test_isotropy(lv2,   "Levy2");
    }

    std::cout << "All swimmer tests passed.\n";
    return 0;
}
