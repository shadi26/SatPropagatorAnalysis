#include <iostream>
#include <fstream>
#include <chrono>
#include <cmath>
#include <vector>
#include <tuple>
#include <iomanip>
#include <algorithm>
#include <string>
#include <sstream>
#include <ctime>
#include <filesystem>
#include "CommonFunctions.h"
#include "RK4.h"
#include "RK8.h"
#include "ODE45.h"
#include "ODE78.h"
#include "ODE113.h"
#include "coefficients78.h"
#include <numeric>

// Define constants
const double EARTH_RADIUS_KM = 6378.137;  // Earth's radius in kilometers
const int NUM_SEGMENTS = 16;
const int NUM_GAUSS_LOBATTO_POINTS = 32;  // Number of points per segment, can be adjusted


// Function to save results to CSV file
void save_results_to_csv(const std::string& satellite_name, const std::string& algorithm_name,
                         const std::vector<double>& time_points, const std::vector<double>& position_differences) {
    std::string directory = "results";
    std::filesystem::create_directory(directory);

    std::string filename = directory + "/" + satellite_name + "_" + algorithm_name + "_results.csv";
    std::ofstream file(filename);

    if (file.is_open()) {
        file << "Time (s), Error (km)\n";
        for (size_t i = 0; i < time_points.size(); ++i) {
            file << time_points[i] << "," << position_differences[i] << "\n";
        }
        std::cout << "Results saved to " << filename << std::endl;
    } else {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }

    file.close();
}

// Function to compare RK4 algorithm and print execution time for each time point
void compare_RK4_algorithm(const std::string& algorithm_name,
                           std::vector<std::vector<double>> (*algorithm_func)(
                               std::vector<double>(*)(double, const std::vector<double>&),
                               const std::vector<double>&,
                               const std::vector<double>&
                           ),
                           const std::vector<double>& r0, const std::vector<double>& v0, const std::vector<double>& time_points) {

    std::cout << "Running " << algorithm_name << " algorithm...\n";

    // Concatenate initial position and velocity
    std::vector<double> y0 = r0;
    y0.insert(y0.end(), v0.begin(), v0.end());

    std::vector<double> execution_times;

    // Loop over each time point to measure execution time for each step
    for (size_t i = 0; i < time_points.size(); ++i) {
        // Define the segment of time_points up to the current time point
        std::vector<double> segment_time_points(time_points.begin(), time_points.begin() + i + 1);

        // Start time for each segment
        auto start_time = std::chrono::high_resolution_clock::now();

        // Execute the algorithm on the current segment
        std::vector<std::vector<double>> results = algorithm_func(satelliteMotion, segment_time_points, y0);

        // End time for each segment
        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end_time - start_time;
        double execution_time = elapsed.count();

        // Store the execution time for the current time point
        execution_times.push_back(execution_time);
    }

    // Print out the execution time for each time point
    std::cout << "Time (s) | Execution Time (s)\n";
    for (size_t i = 0; i < time_points.size(); ++i) {
        std::cout << time_points[i] << " | " << execution_times[i] << " seconds\n";
    }

    // Save results to CSV
    save_results_to_csv("Satellite", algorithm_name, time_points, execution_times);

    // Output total execution time
    double total_execution_time = std::accumulate(execution_times.begin(), execution_times.end(), 0.0);
    std::cout << algorithm_name << " total execution time: " << total_execution_time << " seconds\n";
}


#include <chrono>  // For measuring execution time

void compare_RK8_algorithm(
    const std::string& algorithm_name,
    std::vector<std::vector<double>> (*algorithm_func)(
        std::vector<double>(*)(double, const std::vector<double>&, double),
        const std::vector<double>&,
        const std::vector<double>&,
        double),
    const std::vector<double>& r0, const std::vector<double>& v0,
    const std::vector<double>& time_points, double h = 0.1) {

    std::cout << "Running " << algorithm_name << " algorithm...\n";

    // Concatenate initial position and velocity
    std::vector<double> y0 = r0;
    y0.insert(y0.end(), v0.begin(), v0.end());

    // Vector to store execution time for each time point
    std::vector<double> execution_times;

    // Execute the RK8 algorithm and record execution time at each time point
    for (const auto& t : time_points) {
        auto start_time = std::chrono::high_resolution_clock::now();

        // Execute the RK8 algorithm up to the current time point
        std::vector<std::vector<double>> results = algorithm_func(satellite_motion, {0, t}, y0, h);

        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end_time - start_time;
        execution_times.push_back(elapsed.count());
    }

    // Display or save the execution times with time points
    std::cout << "Execution times for " << algorithm_name << " algorithm:\n";
    std::cout << "Time(s) | Execution Time(s)\n";
    for (size_t i = 0; i < time_points.size(); ++i) {
        std::cout << time_points[i] << "   " << execution_times[i] << "\n";
    }

    // Save results to CSV with execution times instead of position differences
    save_results_to_csv("Satellite", algorithm_name, time_points, execution_times);
}


void compare_ODE45_algorithm(
    const std::string& algorithm_name,
    std::vector<std::vector<double>> (*algorithm_func)(
        std::vector<double>(*)(double, const std::vector<double>&, double),
        const std::vector<double>&,
        const std::vector<double>&,
        double, double, double),
    const std::vector<double>& r0, const std::vector<double>& v0,
    const std::vector<double>& time_points, double mu = 398600, double rtol = 1e-6, double atol = 1e-6) {

    std::cout << "Running " << algorithm_name << " algorithm...\n";

    // Concatenate initial position and velocity into y0
    std::vector<double> y0 = r0;
    y0.insert(y0.end(), v0.begin(), v0.end());

    // Vector to store execution time for each time point
    std::vector<double> execution_times;

    // Execute the ODE45 algorithm and record execution time at each time point
    for (const auto& t : time_points) {
        auto start_time = std::chrono::high_resolution_clock::now();

        // Execute the ODE45 algorithm up to the current time point
        std::vector<std::vector<double>> results = algorithm_func(satellite_motion, {0, t}, y0, mu, rtol, atol);

        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end_time - start_time;
        execution_times.push_back(elapsed.count());
    }


    // Display or save the execution times with time points
    std::cout << "Execution times for " << algorithm_name << " algorithm:\n";
    std::cout << "Time(s) | Execution Time(s)\n";
    for (size_t i = 0; i < time_points.size(); ++i) {
        std::cout << time_points[i] << "   " << execution_times[i] << "\n";
    }

    // Save results to CSV with execution times instead of position differences
    save_results_to_csv("Satellite", algorithm_name, time_points, execution_times);
}

void compare_ODE78_algorithm(const std::string& algorithm_name,
                             const std::vector<double>& r0,
                             const std::vector<double>& v0,
                             const std::vector<double>& time_points,
                             const std::vector<double>& b,
                             const std::vector<double>& bh,
                             const std::vector<double>& c,
                             const std::map<int, std::map<int, double>>& a,
                             double rtol = 1e-3, double atol = 1e-6) {

    std::cout << "Running " << algorithm_name << " algorithm...\n";

    // Concatenate initial position and velocity into a single vector
    std::vector<double> y0 = r0;
    y0.insert(y0.end(), v0.begin(), v0.end());

    // Validate size of y0
    if (y0.empty()) {
        std::cerr << "Error: Initial state vector y0 is empty.\n";
        return;
    }

    // Vector to store execution times for each time point
    std::vector<double> execution_times;

    // Define the ODE function inline for satellite motion (example model)
    auto ode_func = [](double t, const std::vector<double>& y) -> std::vector<double> {
        if (y.size() < 6) {  // Adjust based on expected state vector size (e.g., 6 for position + velocity)
            std::cerr << "Error: State vector y has insufficient size.\n";
            return std::vector<double>(y.size(), 0.0); // Return zero vector as fallback
        }
        std::vector<double> dydt(y.size(), 0.0);
        // Example of satellite motion dynamics placeholder (replace with actual equations)
        dydt[0] = y[3];  // Example: dx/dt = vx
        dydt[1] = y[4];  // Example: dy/dt = vy
        dydt[2] = y[5];  // Example: dz/dt = vz
        // Define derivatives for velocity components (dynamics equations required)
        return dydt;
    };

    // Execute the ODE78 algorithm for each time point
    for (const auto& t : time_points) {
        auto start_time = std::chrono::high_resolution_clock::now();

        // Execute the ODE78 algorithm up to the current time point
        std::vector<std::vector<double>> results = ode78(ode_func, {0, t}, y0, b, bh, c, a, rtol, atol);

        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end_time - start_time;
        execution_times.push_back(elapsed.count());
    }

    // Display or save the execution times with time points
    std::cout << "Execution times for " << algorithm_name << " algorithm:\n";
    std::cout << "Time(s) | Execution Time(s)\n";
    for (size_t i = 0; i < time_points.size(); ++i) {
        std::cout << time_points[i] << "   " << execution_times[i] << "\n";
    }

    // Save execution times to CSV instead of position differences
    save_results_to_csv("Satellite", algorithm_name, time_points, execution_times);
}

// Function to compare ODE113 algorithm execution and print execution time
void compare_ODE113_algorithm(
    const std::string& algorithm_name,
    std::vector<double>(*ode_func)(double, const std::vector<double>&, double),
    const std::vector<double>& r0,
    const std::vector<double>& v0,
    const std::vector<double>& time_points,
    double tol = 1e-6,
    double hmax = 1.0,
    double hmin = 1e-6,
    double mu = 398600
) {
    std::cout << "Running " << algorithm_name << " algorithm...\n";

    // Concatenate initial position and velocity into y0
    std::vector<double> y0 = r0;
    y0.insert(y0.end(), v0.begin(), v0.end());

    // Vector to store execution times for each time point
    std::vector<double> execution_times;

    // Execute the ODE113 algorithm for each time point
    for (const auto& t : time_points) {
        auto start_time = std::chrono::high_resolution_clock::now();

        // Execute the ODE113 algorithm up to the current time point
        std::vector<std::vector<double>> results = ODE113(ode_func, {0, t}, y0, tol, hmax, hmin, mu);

        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end_time - start_time;
        execution_times.push_back(elapsed.count());
    }

    // Display or save the execution times with time points
    std::cout << "Execution times for " << algorithm_name << " algorithm:\n";
    std::cout << "Time(s) | Execution Time(s)\n";
    for (size_t i = 0; i < time_points.size(); ++i) {
        std::cout << time_points[i] << "   " << execution_times[i] << "\n";
    }

    // Save execution times to CSV
    save_results_to_csv("Satellite", algorithm_name, time_points, execution_times);
}


// Function definitions for each algorithm (rk4, rk8, ode45, ode78, ODE113)
// Assume these functions are implemented in other files as per your setup
extern std::vector<std::vector<double>>  rk4(
    std::vector<double>(*odefun)(double, const std::vector<double>&),
    const std::vector<double>& t_gauss_lobatto,
    const std::vector<double>& y0
);
extern std::vector<std::vector<double>> RK8(
    std::vector<double>(*f)(double, const std::vector<double>&, double),
    const std::vector<double>& t_gauss_lobatto, const std::vector<double>& Y0, double h );

extern std::vector<std::vector<double>> ode78(
    std::function<std::vector<double>(double, const std::vector<double>&)> ode_func,
    const std::vector<double>& t_span,
    const std::vector<double>& y0,
    const std::vector<double>& b,  // 8th-order weights
    const std::vector<double>& bh, // 7th-order weights
    const std::vector<double>& c,  // Time nodes (Gauss-Lobatto points)
    const std::map<int, std::map<int, double>>& a,  // Coupling coefficients (Butcher tableau)
    double rtol ,
    double atol
);

// Function to compute and print total execution time using Gauss-Lobatto points in segments
void compare_RK4_algorithm_with_segments(const std::string& algorithm_name,
                                     std::vector<std::vector<double>> (*algorithm_func)(
                                         std::vector<double>(*)(double, const std::vector<double>&),
                                         const std::vector<double>&,
                                         const std::vector<double>&),
                                     const std::vector<double>& r0, const std::vector<double>& v0,
                                     double total_time) {
    double segment_duration = total_time / NUM_SEGMENTS;
    std::vector<double> y0 = r0;
    y0.insert(y0.end(), v0.begin(), v0.end());

    double total_execution_time = 0.0;

    for (int i = 0; i < NUM_SEGMENTS; ++i) {
        double t_start = i * segment_duration;
        double t_end = (i + 1) * segment_duration;

        // Generate Gauss-Lobatto points for the current segment
        std::vector<double> segment_time_points = gaussLobattoPoints(NUM_GAUSS_LOBATTO_POINTS, t_start, t_end);

        auto start_time = std::chrono::high_resolution_clock::now();

        // Execute algorithm on segment
        std::vector<std::vector<double>> results = algorithm_func(satelliteMotion, segment_time_points, y0);

        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end_time - start_time;

        // Accumulate total execution time
        total_execution_time += elapsed.count();

        // Update initial conditions for next segment based on last result
        y0 = results.back();
    }

    // Print the total execution time
    std::cout << algorithm_name << " total execution time: " << total_execution_time << " seconds\n";

    // Save total execution time to CSV
    std::string directory = "results";
    std::filesystem::create_directory(directory);
    std::string filename = directory + "/" + algorithm_name + "_total_execution_time.csv";
    std::ofstream file(filename);

    if (file.is_open()) {
        file << "Algorithm,Total Execution Time (s)\n";
        file << algorithm_name << "," << total_execution_time << "\n";
        std::cout << "Total execution time saved to " << filename << std::endl;
    } else {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }

    // Print the final position
    std::cout << "Final position after " << total_time << " seconds: "
              << y0[0] << ", " << y0[1] << ", " << y0[2] << "\n";

    // Save total execution time and final position to CSV
    file << "Algorithm,Total Execution Time (s),Final Position\n";
    file << algorithm_name << "," << total_execution_time << ","
         << y0[0] << "," << y0[1] << "," << y0[2] << "\n";

    file.close();
}

// Function to compute and print total execution time for RK8 algorithm using Gauss-Lobatto points in segments
void compare_RK8_algorithm_with_segments(const std::string& algorithm_name,
                                         std::vector<std::vector<double>> (*algorithm_func)(
                                             std::vector<double>(*)(double, const std::vector<double>&, double),
                                             const std::vector<double>&,
                                             const std::vector<double>&,
                                             double),
                                         const std::vector<double>& r0, const std::vector<double>& v0,
                                         double total_time, double h = 0.1) {
    double segment_duration = total_time / NUM_SEGMENTS;
    std::vector<double> y0 = r0;
    y0.insert(y0.end(), v0.begin(), v0.end());

    double total_execution_time = 0.0;

    for (int i = 0; i < NUM_SEGMENTS; ++i) {
        double t_start = i * segment_duration;
        double t_end = (i + 1) * segment_duration;

        // Generate Gauss-Lobatto points for the current segment
        std::vector<double> segment_time_points = gaussLobattoPoints(NUM_GAUSS_LOBATTO_POINTS, t_start, t_end);

        auto start_time = std::chrono::high_resolution_clock::now();

        // Execute RK8 algorithm on segment with step size h
        std::vector<std::vector<double>> results = algorithm_func(satellite_motion, segment_time_points, y0, h);

        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end_time - start_time;

        // Accumulate total execution time
        total_execution_time += elapsed.count();

        // Update initial conditions for next segment based on last result
        y0 = results.back();
    }

    // Print the total execution time
    std::cout << algorithm_name << " total execution time: " << total_execution_time << " seconds\n";

    // Save total execution time to CSV
    std::string directory = "results";
    std::filesystem::create_directory(directory);
    std::string filename = directory + "/" + algorithm_name + "_total_execution_time.csv";
    std::ofstream file(filename);

    if (file.is_open()) {
        file << "Algorithm,Total Execution Time (s)\n";
        file << algorithm_name << "," << total_execution_time << "\n";
        std::cout << "Total execution time saved to " << filename << std::endl;
    } else {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }

    // Print the final position
    std::cout << "Final position after " << total_time << " seconds: "
              << y0[0] << ", " << y0[1] << ", " << y0[2] << "\n";

    // Save total execution time and final position to CSV
    file << "Algorithm,Total Execution Time (s),Final Position\n";
    file << algorithm_name << "," << total_execution_time << ","
         << y0[0] << "," << y0[1] << "," << y0[2] << "\n";
    file.close();
}

// Function to compute and print total execution time using Gauss-Lobatto points in segments for ODE45
void compare_ODE45_algorithm_with_segments(const std::string& algorithm_name,
                                           std::vector<std::vector<double>> (*algorithm_func)(
                                               std::vector<double>(*)(double, const std::vector<double>&, double),
                                               const std::vector<double>&,
                                               const std::vector<double>&,
                                               double, double, double),
                                           const std::vector<double>& r0, const std::vector<double>& v0,
                                           double total_time, double mu = 398600, double rtol = 1e-6, double atol = 1e-6) {
    double segment_duration = total_time / NUM_SEGMENTS;
    std::vector<double> y0 = r0;
    y0.insert(y0.end(), v0.begin(), v0.end());

    double total_execution_time = 0.0;

    for (int i = 0; i < NUM_SEGMENTS; ++i) {
        double t_start = i * segment_duration;
        double t_end = (i + 1) * segment_duration;

        // Generate Gauss-Lobatto points for the current segment
        std::vector<double> segment_time_points = gaussLobattoPoints(NUM_GAUSS_LOBATTO_POINTS, t_start, t_end);

        auto start_time = std::chrono::high_resolution_clock::now();

        // Execute the ODE45 algorithm on the current segment
        std::vector<std::vector<double>> results = algorithm_func(satellite_motion, segment_time_points, y0, mu, rtol, atol);

        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end_time - start_time;

        // Accumulate total execution time
        total_execution_time += elapsed.count();

        // Update initial conditions for the next segment based on the last result
        y0 = results.back();
    }

    // Print the total execution time
    std::cout << algorithm_name << " total execution time: " << total_execution_time << " seconds\n";

    // Save total execution time to CSV
    std::string directory = "results";
    std::filesystem::create_directory(directory);
    std::string filename = directory + "/" + algorithm_name + "_total_execution_time.csv";
    std::ofstream file(filename);

    if (file.is_open()) {
        file << "Algorithm,Total Execution Time (s)\n";
        file << algorithm_name << "," << total_execution_time << "\n";
        std::cout << "Total execution time saved to " << filename << std::endl;
    } else {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }

    // Print the final position
    std::cout << "Final position after " << total_time << " seconds: "
              << y0[0] << ", " << y0[1] << ", " << y0[2] << "\n";

    // Save total execution time and final position to CSV
    file << "Algorithm,Total Execution Time (s),Final Position\n";
    file << algorithm_name << "," << total_execution_time << ","
         << y0[0] << "," << y0[1] << "," << y0[2] << "\n";

    file.close();
}

void compare_ODE78_algorithm_with_segments(double total_time1) {
    if (!std::isnan(total_time1)) {  // Check if total_time1 is not NaN
        double total_time = total_time1;
        // Continue with your function logic
    } else {
        std::cerr << "Error: total_time1 is NaN." << std::endl;
    }
    double total_time = total_time1;
    int num_segments = 16;
    double segment_time = total_time / num_segments;
    int n_points = 32;  // Number of Gauss-Lobatto points
    double mu = 398600; // Gravitational parameter
    double tol = 1e-16;
    double rtol = tol / 1000;
    double hmax = 0.01;
    double hmin = 1e-10;

    // Initial position and velocity vectors
    std::vector<double> r0 = {-19946.988367, 11684.423264, 43511.217135};  // Initial position vector
    std::vector<double> v0 = {-0.803367, -1.762325, 0.200044};  // Initial velocity vector
    std::vector<double> y0 = r0;
    y0.insert(y0.end(), v0.begin(), v0.end());

    std::vector<double> time_points;
    std::vector<double> position_differences;

    // Convert c from map<int, double> to vector<double>
    std::vector<double> c_vector(c.size());
    for (const auto& pair : c) {
        c_vector[pair.first - 1] = pair.second;  // Convert map to vector, adjusting index
    }

    auto func = [&mu](double t, const std::vector<double>& y) { return satellite_motion(t, y, mu); };

    // Record the start time
    auto start = std::chrono::high_resolution_clock::now();

    for (int segment = 0; segment < num_segments; ++segment) {
        double t_start = segment * segment_time;
        double t_end = t_start + segment_time;
        auto gauss_lobatto_tspan = gaussLobattoPoints(n_points, t_start, t_end);

        // Run ODE78 with Gauss-Lobatto points for the current segment
        auto segment_result = ode78(func, gauss_lobatto_tspan, y0, b, bh, c_vector, a, tol, rtol);

        // Calculate SGP4 position comparison at each time step
        for (const auto& result : segment_result) {
            double t_point = result[0];
            std::vector<double> current_position(result.begin() + 1, result.begin() + 4);  // Extract position

            time_points.push_back(t_point);
        }

        // Update y0 for the next segment with the last state of this segment
        y0[0] = segment_result.back()[1];
        y0[1] = segment_result.back()[2];
        y0[2] = segment_result.back()[3];
        y0[3] = segment_result.back()[4];
        y0[4] = segment_result.back()[5];
        y0[5] = segment_result.back()[6];
    }

    // Record the end time
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;


    std::cout << "ODE78 Execution time: " << elapsed.count() << " seconds\n";
    std::cout << "ODE78 Final Position: X=" << y0[0] << " Y=" << y0[1] << " Z=" << y0[2] << std::endl;
}

// Function to compute and print total execution time using Gauss-Lobatto points in segments for ODE113
void compare_ODE113_algorithm_with_segments(const std::string& algorithm_name,
                                            std::vector<std::vector<double>> (*algorithm_func)(
                                                std::vector<double>(*)(double, const std::vector<double>&, double),
                                                const std::vector<double>&,
                                                const std::vector<double>&,
                                                double, double, double, double),
                                            const std::vector<double>& r0, const std::vector<double>& v0,
                                            double total_time,
                                            double tol, double hmax, double hmin, double mu) {
    double segment_duration = total_time / NUM_SEGMENTS;
    std::vector<double> y0 = r0;
    y0.insert(y0.end(), v0.begin(), v0.end());

    double total_execution_time = 0.0;

    for (int i = 0; i < NUM_SEGMENTS; ++i) {
        double t_start = i * segment_duration;
        double t_end = (i + 1) * segment_duration;

        // Generate Gauss-Lobatto points for the current segment
        std::vector<double> segment_time_points = gaussLobattoPoints(NUM_GAUSS_LOBATTO_POINTS, t_start, t_end);

        auto start_time = std::chrono::high_resolution_clock::now();

        // Execute the ODE113 algorithm on the current segment
        std::vector<std::vector<double>> results = algorithm_func(satellite_motion, segment_time_points, y0, tol, hmax, hmin, mu);

        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end_time - start_time;

        // Accumulate total execution time
        total_execution_time += elapsed.count();

        // Update initial conditions for the next segment based on the last result
        y0 = results.back();
    }

    // Print the total execution time
    std::cout << algorithm_name << " total execution time: " << total_execution_time << " seconds\n";

    // Save total execution time to CSV
    std::string directory = "results";
    std::filesystem::create_directory(directory);
    std::string filename = directory + "/" + algorithm_name + "_total_execution_time.csv";
    std::ofstream file(filename);

    if (file.is_open()) {
        file << "Algorithm,Total Execution Time (s)\n";
        file << algorithm_name << "," << total_execution_time << "\n";
        std::cout << "Total execution time saved to " << filename << std::endl;
    } else {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }

    // Print the final position
    std::cout << "Final position after " << total_time << " seconds: "
              << y0[0] << ", " << y0[1] << ", " << y0[2] << "\n";

    // Save total execution time and final position to CSV
    file << "Algorithm,Total Execution Time (s),Final Position\n";
    file << algorithm_name << "," << total_execution_time << ","
         << y0[0] << "," << y0[1] << "," << y0[2] << "\n";

    file.close();
}

int main() {
    // Initial conditions
    std::vector<double> r0 = {-19946.988367, 11684.423264, 43511.217135};  // Initial position vector
    std::vector<double> v0 = {-0.803367, -1.762325, 0.200044};  // Initial velocity vector

    // Define parameters for ODE113
    double tol = 1e-6;   // Tolerance
    double hmax = 1.0;   // Maximum step size
    double hmin = 1e-6;  // Minimum step size
    double mu = 398600;  // Gravitational parameter for Earth
    double total_time = 1000000.0;  // Total simulation time in seconds

    // Define the number of Gauss-Lobatto points you want
    const int num_points = 32;  // Adjust as needed for the desired resolution

    // Generate Gauss-Lobatto points from 0 to 3600
    std::vector<double> time_points = gaussLobattoPoints(num_points, 0.0, 3600.0);

    // Compare different algorithms
    compare_RK4_algorithm("RK4", rk4, r0, v0, time_points);

    // Call the specialized function for RK8
    compare_RK8_algorithm("RK8", RK8, r0, v0, time_points, 0.1);

    // Call the specialized function for ODE45
    compare_ODE45_algorithm("ODE45", ode45, r0, v0, time_points, 398600, 1e-6, 1e-6);

    // Run and compare the ODE45 algorithm
    // Convert c_map (std::map<int, double>) to c (std::vector<double>)
    std::vector<double> c_new;
    for (const auto& pair : c) {
        c_new.push_back(pair.second);  // Push only the value part of each map entry into the vector
    }

    // Run and compare the ODE78 algorithm
    compare_ODE78_algorithm("ODE78", r0, v0, time_points, b, bh, c_new, a, 1e-3, 1e-6);

    // Run and compare the ODE113 algorithm
    compare_ODE113_algorithm("ODE113", satellite_motion, r0, v0, time_points, tol, hmax, hmin, mu);


    // Compare RK4
    compare_RK4_algorithm_with_segments("RK4", rk4, r0, v0, 1000000.0);

    // Step size for RK8
    double step_size = 0.1;  // This is the h parameter used in the RK8 algorithm

    // Run the RK8 algorithm in segments with Gauss-Lobatto points and print the total execution time
    compare_RK8_algorithm_with_segments("RK8", RK8, r0, v0, 1000000.0, step_size);


    // Define constants for ODE45
    const double TOTAL_TIME = 1000000.0;  // Total simulation time in seconds
    const double MU = 398600;          // Gravitational parameter for Earth
    const double RTOL = 1e-6;          // Relative tolerance
    const double ATOL = 1e-6;          // Absolute tolerance

    // Run the ODE45 algorithm in segments with Gauss-Lobatto points and print the total execution time
    compare_ODE45_algorithm_with_segments("ODE45", ode45, r0, v0, TOTAL_TIME, MU, RTOL, ATOL);

    // Relative and absolute tolerances for ODE78
    double rtol = 1e-3;
    double atol = 1e-6;


    // Run the comparison for the ODE78 algorithm using segments with Gauss-Lobatto points
    compare_ODE78_algorithm_with_segments(total_time);


    // Call the comparison function for ODE113
    compare_ODE113_algorithm_with_segments("ODE113", ODE113, r0, v0, total_time, tol, hmax, hmin, mu);




    return 0;
}
