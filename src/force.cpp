#include "force.h"
#include <cmath>

static double dipole_x(double x1, double x2, double n1, double n2, const ForceParams& p) {
    double r2 = x1*x1 + x2*x2;
    double r  = std::sqrt(r2);
    if (r < p.b_min) return 0.0;
    double r3    = r2 * r;
    double ndotx = n1*x1 + n2*x2;
    double cos2  = (ndotx * ndotx) / r2;
    return (p.p / r3) * (3.0 * cos2 - 1.0) * x1;
}

static double gaussian_x(double x1, double x2, const ForceParams& p) {
    double r2     = x1*x1 + x2*x2;
    double inv_s2 = 1.0 / (p.sigma * p.sigma);
    return (p.V_0 * inv_s2) * std::exp(-0.5 * r2 * inv_s2) * x1;
}

CoulombEval coulomb_eval(double x1, double x2, double sigma, double b_min) {
    double r2      = x1*x1 + x2*x2;
    double r       = std::sqrt(r2);
    bool   clamped = (r < b_min);
    if (clamped) {
        // Clamp: treat swimmer as sitting at radius b_min along (x1,x2) direction
        // (or at (b_min,0) if exactly at origin)
        double scale = (r > 0.0) ? b_min / r : 1.0;
        x1 *= scale; x2 *= scale;
        r   = b_min;  r2 = b_min * b_min;
    }
    double r3   = r2 * r;
    double r5   = r3 * r2;
    double Fx   = sigma * x1 / r3;
    double Fy   = sigma * x2 / r3;
    double dxFx = sigma * (1.0/r3 - 3.0*x1*x1/r5);
    double dyFx = sigma * (-3.0*x1*x2/r5);
    return {Fx, Fy, dxFx, dyFx, clamped};
}

FDGrad fd_grad_Fx(double x1, double x2, double n1, double n2,
                   const ForceParams& p, double h) {
    // Central differences: d_x F_x and d_y F_x, orientation n fixed.
    // near_core guard: flag if any shifted point enters the hard core (r < b_min).
    double bm = p.b_min;
    auto nc = [bm](double a, double b) {
        return std::sqrt(a*a + b*b) < bm;
    };
    bool near_core = nc(x1+h, x2) || nc(x1-h, x2) || nc(x1, x2+h) || nc(x1, x2-h);
    double dxFx = (force_x(x1+h, x2, n1, n2, p) - force_x(x1-h, x2, n1, n2, p)) / (2.0*h);
    double dyFx = (force_x(x1, x2+h, n1, n2, p) - force_x(x1, x2-h, n1, n2, p)) / (2.0*h);
    return {dxFx, dyFx, near_core};
}

double force_x(double x1, double x2, double n1, double n2, const ForceParams& params) {
    switch (params.type) {
        case ForceType::DIPOLE:   return dipole_x(x1, x2, n1, n2, params);
        case ForceType::GAUSSIAN: return gaussian_x(x1, x2, params);
        case ForceType::COULOMB:  return coulomb_eval(x1, x2, params.sigma, params.b_min).Fx;
    }
    return 0.0;
}

std::array<double, 2> force_2d(double x1, double x2, double n1, double n2, const ForceParams& params) {
    if (params.type == ForceType::GAUSSIAN) {
        double r2     = x1*x1 + x2*x2;
        double inv_s2 = 1.0 / (params.sigma * params.sigma);
        double fac    = (params.V_0 * inv_s2) * std::exp(-0.5 * r2 * inv_s2);
        return {fac * x1, fac * x2};
    }
    if (params.type == ForceType::COULOMB) {
        auto e = coulomb_eval(x1, x2, params.sigma, params.b_min);
        return {e.Fx, e.Fy};
    }
    // DIPOLE
    double r2 = x1*x1 + x2*x2;
    double r  = std::sqrt(r2);
    if (r < params.b_min) return {0.0, 0.0};
    double r3    = r2 * r;
    double ndotx = n1*x1 + n2*x2;
    double cos2  = (ndotx * ndotx) / r2;
    double prefac = (params.p / r3) * (3.0 * cos2 - 1.0);
    return {prefac * x1, prefac * x2};
}
