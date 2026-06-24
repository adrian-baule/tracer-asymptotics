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

double force_x(double x1, double x2, double n1, double n2, const ForceParams& params) {
    switch (params.type) {
        case ForceType::DIPOLE:   return dipole_x(x1, x2, n1, n2, params);
        case ForceType::GAUSSIAN: return gaussian_x(x1, x2, params);
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
    double r2 = x1*x1 + x2*x2;
    double r  = std::sqrt(r2);
    if (r < params.b_min) return {0.0, 0.0};
    double r3    = r2 * r;
    double ndotx = n1*x1 + n2*x2;
    double cos2  = (ndotx * ndotx) / r2;
    double prefac = (params.p / r3) * (3.0 * cos2 - 1.0);
    return {prefac * x1, prefac * x2};
}
