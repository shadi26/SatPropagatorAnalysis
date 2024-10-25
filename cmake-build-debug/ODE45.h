// ODE45.h
#ifndef ODE45_H
#define ODE45_H

#include <vector>

// ODE45 function declaration
std::vector<std::vector<double>> ode45(
    std::vector<double>(*func)(double, const std::vector<double>&, double),
    const std::vector<double>& t_span, const std::vector<double>& y0, double mu, double rtol , double atol );

#endif // ODE45_H
