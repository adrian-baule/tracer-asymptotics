#pragma once
#include <array>

enum class ForceType { DIPOLE, GAUSSIAN, COULOMB };

struct ForceParams {
    ForceType type  = ForceType::DIPOLE;
    // Dipole parameters
    double p        = 1.0;
    double b_min    = 0.5;
    // Gaussian parameters (bmshort)
    double V_0      = 1.0;
    double sigma    = 1.0;   // Gaussian range (bmshort) or Coulomb strength (bmlong)
};

// Force x-component dispatches on ForceParams::type.
// For DIPOLE:   f_x = (p/r^3)*(3*(n.x)^2/r^2-1)*x1, zero if r < b_min.
// For GAUSSIAN: f_x = (V_0/sigma^2)*exp(-r^2/(2*sigma^2))*x1; n ignored.
// For COULOMB:  f_x = sigma*x1/r^3, clamped to r >= b_min.
double force_x(double x1, double x2, double n1, double n2, const ForceParams& params);

// Full 2D force vector (for diagnostics / A1 accumulation)
std::array<double, 2> force_2d(double x1, double x2, double n1, double n2, const ForceParams& params);

// Coulomb evaluation bundle: force vector + x-gradient components needed for A2.
// r is clamped to b_min before evaluation; clamped=true if the clamp was active.
struct CoulombEval {
    double Fx, Fy;   // force vector components
    double dxFx;     // d/dx1 F_x = sigma*(1/r^3 - 3*x1^2/r^5)
    double dyFx;     // d/dx2 F_x = sigma*(-3*x1*x2/r^5)
    bool   clamped;  // true if r < b_min was enforced
};
CoulombEval coulomb_eval(double x1, double x2, double sigma, double b_min);
