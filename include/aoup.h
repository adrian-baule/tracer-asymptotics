#pragma once
#include <random>
#include <array>

// State of a single AOUP: position and velocity in 2D
struct AOUPState {
    double x, y;   // position
    double vx, vy; // velocity
};

// Parameters for AOUP dynamics
struct AOUPParams {
    double tau_c;  // persistence time
    double D_A;    // active diffusivity
    double dt;     // timestep
};

// Initialise AOUP at position (x0, y0) with velocity drawn from stationary distribution
AOUPState aoup_init(double x0, double y0, const AOUPParams& params, std::mt19937_64& rng);

// Advance AOUP by one timestep (Euler-Maruyama)
void aoup_step(AOUPState& state, const AOUPParams& params, std::mt19937_64& rng);

// Return orientation unit vector n = V/|V| (undefined if |V|=0, handled by caller)
std::array<double, 2> aoup_orientation(const AOUPState& state);
