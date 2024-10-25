// ODE113.h
#ifndef ODE113_H
#define ODE113_H

#include <vector>

// ODE113 function declaration
std::vector<std::vector<double>> ODE113(
    std::vector<double>(*f)(double, const std::vector<double>&, double),
    const std::vector<double>& time_points,
    const std::vector<double>& y0,
    double tol,
    double hmax = 1.0 ,
    double hmin = 1e-6,
    double mu = 398600
);

#endif // ODE113_H
