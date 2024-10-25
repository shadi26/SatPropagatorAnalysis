#include <iostream>
#include <vector>
#include <cmath>
#include <tuple>
#include <algorithm>
#include <iomanip>
#include "CommonFunctions.h"
#include "ODE45.h"
// Helper function for vector operations (operator overloading provided at the top of the code)
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif


// Function to perform a single RK45 step
std::tuple<double, std::vector<double>, double> rk45_step(
    std::vector<double>(*func)(double, const std::vector<double>&, double),
    double t, const std::vector<double>& y, double h, double rtol, double atol, double mu) {

    // Coefficients from the provided Butcher tableau
    std::vector<double> a = {0, 1.0 / 5.0, 3.0 / 10.0, 4.0 / 5.0, 8.0 / 9.0, 1.0, 1.0};
    std::vector<std::vector<double>> b = {
        {0},
        {1.0 / 5.0},
        {3.0 / 40.0, 9.0 / 40.0},
        {44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0},
        {19372.0 / 6561.0, -25360.0 / 2187.0, 64448.0 / 6561.0, -212.0 / 729.0},
        {9017.0 / 3168.0, -355.0 / 33.0, 46732.0 / 5247.0, 49.0 / 176.0, -5103.0 / 18656.0},
        {35.0 / 384.0, 0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0}
    };

    std::vector<double> c4 = {35.0 / 384.0, 0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0, 0};
    std::vector<double> c5 = {5179.0 / 57600.0, 0, 7571.0 / 16695.0, 393.0 / 640.0, -92097.0 / 339200.0, 187.0 / 2100.0, 1.0 / 40.0};

    // Calculate the k values using the overloaded operators
std::vector<double> k1 = h * func(t, y, mu);
std::vector<double> k2 = h * func(t + a[1] * h, y + b[1][0] * k1, mu);
std::vector<double> k3 = h * func(t + a[2] * h, y + b[2][0] * k1 + b[2][1] * k2, mu);
std::vector<double> k4 = h * func(t + a[3] * h, y + b[3][0] * k1 + b[3][1] * k2 + b[3][2] * k3, mu);
std::vector<double> k5 = h * func(t + a[4] * h, y + b[4][0] * k1 + b[4][1] * k2 + b[4][2] * k3 + b[4][3] * k4, mu);
std::vector<double> k6 = h * func(t + a[5] * h, y + b[5][0] * k1 + b[5][1] * k2 + b[5][2] * k3 + b[5][3] * k4 + b[5][4] * k5, mu);
std::vector<double> k7 = h * func(t + a[6] * h, y + b[6][0] * k1 + b[6][1] * k2 + b[6][2] * k3 + b[6][3] * k4 + b[6][4] * k5 + b[6][5] * k6, mu);


    // Compute the 4th and 5th order solutions
    std::vector<double> y4 = y;
    std::vector<double> y5 = y;
    for (size_t i = 0; i < y.size(); ++i) {
        y4[i] += c4[0] * k1[i] + c4[2] * k3[i] + c4[3] * k4[i] + c4[4] * k5[i] + c4[5] * k6[i];
        y5[i] += c5[0] * k1[i] + c5[2] * k3[i] + c5[3] * k4[i] + c5[4] * k5[i] + c5[5] * k6[i] + c5[6] * k7[i];
    }

    // Estimate the error
    double error = 0.0;
    for (size_t i = 0; i < y.size(); ++i) {
        double err_i = std::abs(y5[i] - y4[i]) / (atol + rtol * std::max(std::abs(y4[i]), std::abs(y5[i])));
        if (err_i > error) {
            error = err_i;
        }
    }

    // Adjust the step size
    double h_new = h * std::min(2.0, std::max(0.1, 0.9 / std::pow(error, 0.2)));
    return {t + h, y5, h_new};
}

// ODE45 solver using adaptive step-size control with the RK45 method
std::vector<std::vector<double>> ode45(
    std::vector<double>(*func)(double, const std::vector<double>&, double),
    const std::vector<double>& t_span, const std::vector<double>& y0, double mu = 398600, double rtol = 1e-6, double atol = 1e-6) {

    std::vector<double> tout = {t_span[0]};
    std::vector<std::vector<double>> yout = {y0};

    double t = t_span[0];
    std::vector<double> y = y0;

    for (size_t i = 1; i < t_span.size(); ++i) {
        double h = t_span[i] - t_span[i - 1];  // Step size
        while (t < t_span[i]) {
            auto [t_next, y_next, h_new] = rk45_step(func, t, y, h, rtol, atol, mu);
            t = t_next;
            y = y_next;
            h = h_new;
        }

        tout.push_back(t);
        yout.push_back(y);
    }

    return yout;
}



/*
int main() {
    // Define initial conditions and time span
    std::vector<double> r0 = {-3314.58372, 2809.02481, -5400.46222};  // Initial position vector
    std::vector<double> v0 = {-3.42701791, -6.62341508, -1.34238849};  // Initial velocity vector
    std::vector<double> y0 = {r0[0], r0[1], r0[2], v0[0], v0[1], v0[2]};  // Initial state vector

    double t_start = 0;
    double t_end = 3600;
    int n_points = 32;

    // Generate Gauss-Lobatto points for the given time range
    std::vector<double> gauss_lobatto_tspan = gaussLobattoPoints(n_points, t_start, t_end);

    // Define satellite motion function
    auto func = [](double t, const std::vector<double>& y, double mu) -> std::vector<double> {
        return satellite_motion(t, y, mu);
    };

    // Solve ODE using ode45 method
    auto result = ode45(func, gauss_lobatto_tspan, y0);

    // Print results
    std::cout << "Results (Time | Position X, Y, Z | Velocity X, Y, Z):\n";
    std::cout << std::setw(8) << "Time" << std::setw(12) << "PosX" << std::setw(12) << "PosY"
              << std::setw(12) << "PosZ" << std::setw(12) << "VelX" << std::setw(12) << "VelY"
              << std::setw(12) << "VelZ" << "\n";
    std::cout << std::string(80, '-') << "\n";

    for (size_t i = 0; i < result.size(); ++i) {
        std::cout << std::setw(8) << gauss_lobatto_tspan[i] << std::setw(12) << result[i][0] << std::setw(12) << result[i][1]
                  << std::setw(12) << result[i][2] << std::setw(12) << result[i][3] << std::setw(12) << result[i][4]
                  << std::setw(12) << result[i][5] << "\n";
    }

    return 0;
}
*/
