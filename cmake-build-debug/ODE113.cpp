#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <tuple>
#include <iomanip>
#include "CommonFunctions.h"

// Helper function for vector operations (operator overloading provided at the top of the code)
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif



// Function to perform a single step of the classic 4th-order Runge-Kutta method
std::vector<double> rk4_step(std::vector<double>(*f)(double, const std::vector<double>&, double), double t, const std::vector<double>& y, double h, double mu) {
    std::vector<double> k1 = h * f(t, y, mu);
    std::vector<double> k2 = h * f(t + 0.5 * h, y + 0.5 * k1, mu);
    std::vector<double> k3 = h * f(t + 0.5 * h, y + 0.5 * k2, mu);
    std::vector<double> k4 = h * f(t + h, y + k3, mu);
    return y + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
}

// Adams-Bashforth-Moulton method (ODE113) using fixed time points to solve ODEs
std::vector<std::vector<double>> ODE113(
    std::vector<double>(*f)(double, const std::vector<double>&, double),
    const std::vector<double>& time_points,
    const std::vector<double>& y0,
    double tol,
    double hmax ,
    double hmin ,
    double mu
){
    std::vector<double> t_values_list = {time_points[0]};
    std::vector<std::vector<double>> y_values_list = {y0};

    // Initial conditions for the first few steps using RK4
    for (int i = 0; i < 3 && i < time_points.size() - 1; ++i) {
        double t = time_points[i];
        double h = time_points[i + 1] - t;
        std::vector<double> y_next = rk4_step(f, t, y_values_list.back(), h, mu);
        t_values_list.push_back(time_points[i + 1]);
        y_values_list.push_back(y_next);
    }

    // Adams-Bashforth-Moulton method
    for (int i = 3; i < time_points.size() - 1; ++i) {
        double t = time_points[i];
        double h = time_points[i + 1] - t;

        // Adams-Bashforth predictor
        std::vector<double> y_pred = y_values_list.back() + h / 24.0 * (
            55.0 * f(time_points[i], y_values_list.back(), mu)
            - 59.0 * f(time_points[i - 1], y_values_list[i - 1], mu)
            + 37.0 * f(time_points[i - 2], y_values_list[i - 2], mu)
            - 9.0 * f(time_points[i - 3], y_values_list[i - 3], mu)
        );

        // Adams-Moulton corrector
        std::vector<double> y_correct = y_values_list.back() + h / 24.0 * (
            9.0 * f(time_points[i + 1], y_pred, mu)
            + 19.0 * f(time_points[i], y_values_list[i], mu)
            - 5.0 * f(time_points[i - 1], y_values_list[i - 1], mu)
            + 1.0 * f(time_points[i - 2], y_values_list[i - 2], mu)
        );

        t_values_list.push_back(time_points[i + 1]);
        y_values_list.push_back(y_correct);
    }

    return y_values_list;
}

/*
int main() {
// Initial conditions
std::vector<double> r0 = {-3829.29, 5677.86, -1385.16};  // Initial position vector
std::vector<double> v0 = {-1.69535, -0.63752, 7.33375};  // Initial velocity vector
std::vector<double> y0(6);
std::copy(r0.begin(), r0.end(), y0.begin());
std::copy(v0.begin(), v0.end(), y0.begin() + 3);

// Define the time range and the number of Gauss-Lobatto points
double t_start = 0;
double t_end = 3600;  // For example, a time span from 0 to 3600 seconds
int n_points = 20;    // Number of Gauss-Lobatto points

// Generate Gauss-Lobatto points within the time range
std::vector<double> t_gauss_lobatto = gaussLobattoPoints(n_points, t_start, t_end);

// Define the function for satellite motion
auto func = [](double t, const std::vector<double>& y) -> std::vector<double> {
return satelliteMotion(t, y);
};

// Perform integration using ODE113 with Gauss-Lobatto points
std::vector<std::vector<double>> results = ODE113(satellite_motion, t_gauss_lobatto, y0, 1e-6);


// Print results in an organized format
std::cout << std::setw(8) << "Time" << std::setw(12) << "PosX" << std::setw(12) << "PosY"
<< std::setw(12) << "PosZ" << std::setw(12) << "VelX" << std::setw(12) << "VelY"
<< std::setw(12) << "VelZ" << std::endl;
std::cout << std::string(80, '-') << std::endl;

for (size_t i = 0; i < results.size(); ++i) {
std::cout << std::fixed << std::setprecision(2);
std::cout << std::setw(8) << t_gauss_lobatto[i] << std::setw(12) << results[i][0] << std::setw(12) << results[i][1]
<< std::setw(12) << results[i][2] << std::setw(12) << results[i][3] << std::setw(12) << results[i][4]
<< std::setw(12) << results[i][5] << std::endl;
}

return 0;
}
*/


