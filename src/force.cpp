#include "force.h"
#include <cmath>

double force_x(double x1, double x2, double n1, double n2, const ForceParams& params) {
    double r2 = x1*x1 + x2*x2;
    double r  = std::sqrt(r2);
    if (r < params.b_min) return 0.0;
    double r3    = r2 * r;
    double ndotx = n1*x1 + n2*x2;
    double cos2  = (ndotx * ndotx) / r2;
    return (params.p / r3) * (3.0 * cos2 - 1.0) * x1;
}

std::array<double, 2> force_2d(double x1, double x2, double n1, double n2, const ForceParams& params) {
    double r2 = x1*x1 + x2*x2;
    double r  = std::sqrt(r2);
    if (r < params.b_min) return {0.0, 0.0};
    double r3    = r2 * r;
    double ndotx = n1*x1 + n2*x2;
    double cos2  = (ndotx * ndotx) / r2;
    double prefac = (params.p / r3) * (3.0 * cos2 - 1.0);
    return {prefac * x1, prefac * x2};
}
