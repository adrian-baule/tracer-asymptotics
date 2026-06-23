#pragma once
#include <random>
#include <array>

struct SwimmerState {
    double x, y;
    double vx, vy;
    double aux;  // Levy walk: remaining run time; ignored by other processes
};

class Swimmer {
public:
    virtual ~Swimmer() = default;
    virtual SwimmerState init(double x0, double y0, std::mt19937_64& rng) const = 0;
    virtual void step(SwimmerState& state, std::mt19937_64& rng) const = 0;
    virtual std::array<double, 2> orientation(const SwimmerState& state) const = 0;
};
