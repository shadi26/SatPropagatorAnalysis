// RK4.cpp
#include <vector>
#include <iostream>
#include <cmath>
#include <iomanip>
#include "CommonFunctions.h"
#include "RK4.h"

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif


// 4th order Runge-Kutta method (RK4)
std::vector<std::vector<double>> rk4(
    std::vector<double>(*odefun)(double, const std::vector<double>&),
    const std::vector<double>& t_gauss_lobatto, 
    const std::vector<double>& y0
) {
    std::vector<std::vector<double>> yout(t_gauss_lobatto.size(), std::vector<double>(y0.size()));
    std::vector<double> y = y0;
    yout[0] = y;

    // Perform the integration using RK4
    for (size_t i = 1; i < t_gauss_lobatto.size(); ++i) {
        double h = t_gauss_lobatto[i] - t_gauss_lobatto[i - 1];

        std::vector<double> k1 = odefun(t_gauss_lobatto[i - 1], y);
        std::vector<double> k2 = odefun(t_gauss_lobatto[i - 1] + h / 2.0, vectorAdd(y, scalarMultiply(k1, h / 2.0)));
        std::vector<double> k3 = odefun(t_gauss_lobatto[i - 1] + h / 2.0, vectorAdd(y, scalarMultiply(k2, h / 2.0)));
        std::vector<double> k4 = odefun(t_gauss_lobatto[i - 1] + h, vectorAdd(y, scalarMultiply(k3, h)));

        for (size_t j = 0; j < y.size(); ++j) {
            y[j] += h * (k1[j] + 2.0 * k2[j] + 2.0 * k3[j] + k4[j]) / 6.0;
        }
        yout[i] = y;
    }

    return yout;
}

/*
int main() {
    // Define the initial conditions
    std::vector<double> r0 = {-3829.29, 5677.86, -1385.16};  // Initial position
    std::vector<double> v0 = {-1.69535, -0.63752, 7.33375};  // Initial velocity
    std::vector<double> y0 = {r0[0], r0[1], r0[2], v0[0], v0[1], v0[2]};  // Initial state vector

    double TSPAN_start = 0, TSPAN_end = 3600;
    int n_points = 20;  // Number of Gauss-Lobatto points
    std::vector<double> gauss_lobatto_tspan = gaussLobattoPoints(n_points, TSPAN_start, TSPAN_end);

    // Define the satellite motion function
    auto func = [](double t, const std::vector<double>& y, double mu) -> std::vector<double> {
        return satelliteMotion(t, y);
    };

    // Run RK8 integration
    std::vector<std::vector<double>> results = RK8(func, gauss_lobatto_tspan, y0);

    // Print results in organized format
    std::cout << "Results (Time | Position X, Y, Z | Velocity X, Y, Z):\n";
    std::cout << std::setw(8) << "Time" << std::setw(12) << "PosX" << std::setw(12) << "PosY"
              << std::setw(12) << "PosZ" << std::setw(12) << "VelX" << std::setw(12) << "VelY"
              << std::setw(12) << "VelZ" << "\n";
    std::cout << std::string(80, '-') << "\n";

    for (size_t i = 0; i < results.size(); ++i) {
        std::cout << std::setw(8) << gauss_lobatto_tspan[i] << std::setw(12) << results[i][0] << std::setw(12) << results[i][1]
                  << std::setw(12) << results[i][2] << std::setw(12) << results[i][3] << std::setw(12) << results[i][4]
                  << std::setw(12) << results[i][5] << "\n";
    }

    return 0;
}
*/





