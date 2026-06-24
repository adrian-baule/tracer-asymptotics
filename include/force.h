#pragma once
#include <array>

enum class ForceType { DIPOLE, GAUSSIAN };

struct ForceParams {
    ForceType type  = ForceType::DIPOLE;
    // Dipole parameters
    double p        = 1.0;
    double b_min    = 0.5;
    // Gaussian parameters
    double V_0      = 1.0;
    double sigma    = 1.0;
};

// Force x-component dispatches on ForceParams::type.
// For DIPOLE: f_x = (p/|x|^3) * (3*(n.x)^2/|x|^2 - 1) * x1, zero if |x| < b_min.
// For GAUSSIAN: f_x = (V_0/sigma^2) * exp(-|x|^2/(2*sigma^2)) * x1; n ignored.
double force_x(double x1, double x2, double n1, double n2, const ForceParams& params);

// Full 2D force vector (for diagnostics)
std::array<double, 2> force_2d(double x1, double x2, double n1, double n2, const ForceParams& params);
