#pragma once
#include <array>

// Parameters for the hydrodynamic dipole force
struct ForceParams {
    double p;      // dipole strength (p < 0: pusher, p > 0: puller)
    double b_min;  // hard-core cutoff: force set to zero if |x| < b_min
};

// Hydrodynamic dipole force x-component
// x: vector from tracer (at origin) to swimmer
// n: swimmer orientation unit vector
// Returns f_x = (p/|x|^3) * (3*(n.x)^2/|x|^2 - 1) * x_1
// Returns 0.0 if |x| < b_min
double force_x(double x1, double x2, double n1, double n2, const ForceParams& params);

// Full 2D force vector (for diagnostics)
std::array<double, 2> force_2d(double x1, double x2, double n1, double n2, const ForceParams& params);
