#ifndef ODE113_H
#define ODE113_H

#include <vector>

// Structure to hold the result of ODE113
struct ODE113Result {
    std::vector<double> time;               // Time points
    std::vector<std::vector<double>> values; // Corresponding solution values
};

// ODE113 function declaration
ODE113Result ODE113(
    std::vector<double> (*ode)(double, const std::vector<double>&, double), // ODE function
    const std::vector<double>& time_points,  // Time span or specific points
    const std::vector<double>& y0,          // Initial conditions
    double mu,                               // Gravitational parameter or other constant
    double rel_tol = 1e-9,                   // Relative tolerance
    double abs_tol = 1e-9,                   // Absolute tolerance
    double hmax = 1.0,                       // Maximum step size
    double hmin = 1e-6                       // Minimum step size
);

#endif // ODE113_H
