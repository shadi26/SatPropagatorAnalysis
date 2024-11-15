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

// Function to save final positions to a CSV file
void save_final_positions_to_csv(const std::string& satellite_name, const std::string& algorithm_name,
                                 const std::vector<double>& time_points,
                                 const std::vector<std::vector<double>>& positions) {
    std::string directory = "results";
    std::filesystem::create_directory(directory);

    std::string filename = directory + "/" + satellite_name + "_" + algorithm_name + "_final_positions.csv";
    std::ofstream file(filename);

    if (file.is_open()) {
        file << "Time (s), X (km), Y (km), Z (km)\n";
        for (size_t i = 0; i < time_points.size(); ++i) {
            file << time_points[i] << "," << positions[i][0] << "," << positions[i][1] << "," << positions[i][2] << "\n";
        }
        std::cout << "Final positions saved to " << filename << std::endl;
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
    std::vector<std::vector<double>> final_positions;

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

        final_positions.push_back({results.back()[0], results.back()[1], results.back()[2]});
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
    save_final_positions_to_csv("Satellite", algorithm_name, time_points, final_positions);

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
    std::vector<std::vector<double>> final_positions;

    // Execute the RK8 algorithm and record execution time at each time point
    for (const auto& t : time_points) {
        auto start_time = std::chrono::high_resolution_clock::now();

        // Execute the RK8 algorithm up to the current time point
        std::vector<std::vector<double>> results = algorithm_func(satellite_motion, {0, t}, y0, h);

        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end_time - start_time;
        execution_times.push_back(elapsed.count());
        final_positions.push_back({results.back()[0], results.back()[1], results.back()[2]});
    }

    // Display or save the execution times with time points
    std::cout << "Execution times for " << algorithm_name << " algorithm:\n";
    std::cout << "Time(s) | Execution Time(s)\n";
    for (size_t i = 0; i < time_points.size(); ++i) {
        std::cout << time_points[i] << "   " << execution_times[i] << "\n";
    }

    // Save results to CSV with execution times instead of position differences
    save_results_to_csv("Satellite", algorithm_name, time_points, execution_times);
    save_final_positions_to_csv("Satellite", algorithm_name, time_points, final_positions);
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
    std::vector<std::vector<double>> final_positions;

    // Execute the ODE45 algorithm and record execution time at each time point
    for (const auto& t : time_points) {
        auto start_time = std::chrono::high_resolution_clock::now();

        // Execute the ODE45 algorithm up to the current time point
        std::vector<std::vector<double>> results = algorithm_func(satellite_motion, {0, t}, y0, mu, rtol, atol);

        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end_time - start_time;
        execution_times.push_back(elapsed.count());
        final_positions.push_back({results.back()[0], results.back()[1], results.back()[2]});
    }


    // Display or save the execution times with time points
    std::cout << "Execution times for " << algorithm_name << " algorithm:\n";
    std::cout << "Time(s) | Execution Time(s)\n";
    for (size_t i = 0; i < time_points.size(); ++i) {
        std::cout << time_points[i] << "   " << execution_times[i] << "\n";
    }

    // Save results to CSV with execution times instead of position differences
    save_results_to_csv("Satellite", algorithm_name, time_points, execution_times);
    save_final_positions_to_csv("Satellite", algorithm_name, time_points, final_positions);
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
    std::vector<std::vector<double>> final_positions;

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
        final_positions.push_back({results.back()[0], results.back()[1], results.back()[2]});
    }

    // Display or save the execution times with time points
    std::cout << "Execution times for " << algorithm_name << " algorithm:\n";
    std::cout << "Time(s) | Execution Time(s)\n";
    for (size_t i = 0; i < time_points.size(); ++i) {
        std::cout << time_points[i] << "   " << execution_times[i] << "\n";
    }

    // Save execution times to CSV instead of position differences
    save_results_to_csv("Satellite", algorithm_name, time_points, execution_times);
    save_final_positions_to_csv("Satellite", algorithm_name, time_points, final_positions);
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
    std::vector<std::vector<double>> final_positions;

    // Execute the ODE113 algorithm for each time point
    for (const auto& t : time_points) {
        auto start_time = std::chrono::high_resolution_clock::now();

        // Execute the ODE113 algorithm up to the current time point
        std::vector<std::vector<double>> results = ODE113(ode_func, {0, t}, y0, tol, hmax, hmin, mu);

        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end_time - start_time;
        execution_times.push_back(elapsed.count());
        final_positions.push_back({results.back()[0], results.back()[1], results.back()[2]});
    }

    // Display or save the execution times with time points
    std::cout << "Execution times for " << algorithm_name << " algorithm:\n";
    std::cout << "Time(s) | Execution Time(s)\n";
    for (size_t i = 0; i < time_points.size(); ++i) {
        std::cout << time_points[i] << "   " << execution_times[i] << "\n";
    }

    // Save execution times to CSV
    save_results_to_csv("Satellite", algorithm_name, time_points, execution_times);
    save_final_positions_to_csv("Satellite", algorithm_name, time_points, final_positions);
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

void compare_RK4_algorithm_with_segments(const std::string& algorithm_name,
                                         std::vector<std::vector<double>> (*algorithm_func)(
                                                 std::vector<double>(*)(double, const std::vector<double>&),
                                                 const std::vector<double>&,
                                                 const std::vector<double>&),
                                         const std::vector<double>& r0, const std::vector<double>& v0,
                                         double total_time) {
    double segment_duration = total_time / NUM_SEGMENTS; // Duration of each segment
    std::vector<double> y0 = r0; // Initial position
    y0.insert(y0.end(), v0.begin(), v0.end()); // Append initial velocity to form state vector

    double total_execution_time = 0.0;

    // Vectors to store execution times, final positions, and time points
    std::vector<double> execution_times;
    std::vector<std::vector<double>> final_positions;
    std::vector<double> time_points;

    for (int segment = 0; segment < NUM_SEGMENTS; ++segment) {
        double t_start = segment * segment_duration;
        double t_end = (segment + 1) * segment_duration;

        // Generate Gauss-Lobatto points for the current segment
        std::vector<double> segment_time_points = gaussLobattoPoints(NUM_GAUSS_LOBATTO_POINTS, t_start, t_end);

        // Start measuring execution time for this segment
        auto start_time = std::chrono::high_resolution_clock::now();

        // Execute the algorithm for all Gauss-Lobatto points in this segment
        std::vector<std::vector<double>> results = algorithm_func(satelliteMotion, segment_time_points, y0);

        // Stop measuring execution time for this segment
        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end_time - start_time;

        // Accumulate total execution time
        double execution_time = elapsed.count();
        total_execution_time += execution_time;

        // Iterate through each result and store positions and times
        for (size_t i = 0; i < segment_time_points.size(); ++i) {
            time_points.push_back(segment_time_points[i]); // Store time point
            final_positions.push_back({results[i][0], results[i][1], results[i][2]}); // Store position
        }

        // Update y0 to the last state of this segment
        y0 = results.back();

        // Store the execution time for the entire segment
        execution_times.push_back(execution_time);
    }

    // Save execution times to a CSV file (only for segments)
    save_results_to_csv("Satellite", algorithm_name, time_points, execution_times);

    // Save final positions to a CSV file
    save_final_positions_to_csv("Satellite", algorithm_name, time_points, final_positions);

    // Print total execution time
    std::cout << algorithm_name << " total execution time: " << total_execution_time << " seconds\n";

    // Print the final position after the last segment
    std::cout << "Final position after " << total_time << " seconds: "
              << y0[0] << ", " << y0[1] << ", " << y0[2] << "\n";
}





void compare_RK8_algorithm_with_segments(const std::string& algorithm_name,
                                         std::vector<std::vector<double>> (*algorithm_func)(
                                                 std::vector<double>(*)(double, const std::vector<double>&, double),
                                                 const std::vector<double>&,
                                                 const std::vector<double>&,
                                                 double),
                                         const std::vector<double>& r0, const std::vector<double>& v0,
                                         double total_time, double h = 0.1) {
    double segment_duration = total_time / NUM_SEGMENTS; // Duration of each segment
    std::vector<double> y0 = r0; // Initial position
    y0.insert(y0.end(), v0.begin(), v0.end()); // Append initial velocity to form state vector

    double total_execution_time = 0.0;

    // Vectors to store execution times, final positions, and time points
    std::vector<double> execution_times;
    std::vector<std::vector<double>> final_positions;
    std::vector<double> time_points;

    for (int segment = 0; segment < NUM_SEGMENTS; ++segment) {
        double t_start = segment * segment_duration;
        double t_end = (segment + 1) * segment_duration;

        // Generate Gauss-Lobatto points for the current segment
        std::vector<double> segment_time_points = gaussLobattoPoints(NUM_GAUSS_LOBATTO_POINTS, t_start, t_end);

        auto start_time = std::chrono::high_resolution_clock::now();

        // Execute the RK8 algorithm for all Gauss-Lobatto points in this segment
        std::vector<std::vector<double>> results = algorithm_func(satellite_motion, segment_time_points, y0, h);

        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end_time - start_time;

        // Accumulate total execution time for the segment
        double execution_time = elapsed.count();
        total_execution_time += execution_time;

        // Iterate through each result and store positions and times
        for (size_t i = 0; i < segment_time_points.size(); ++i) {
            time_points.push_back(segment_time_points[i]); // Store time point
            final_positions.push_back({results[i][0], results[i][1], results[i][2]}); // Store position
        }

        // Update y0 to the last state of this segment
        y0 = results.back();

        // Store execution time for the entire segment
        execution_times.push_back(execution_time);
    }

    // Save execution times to a CSV file (segment-level granularity)
    save_results_to_csv("Satellite", algorithm_name, time_points, execution_times);

    // Save final positions to a CSV file (time point-level granularity)
    save_final_positions_to_csv("Satellite", algorithm_name, time_points, final_positions);

    // Print total execution time
    std::cout << algorithm_name << " total execution time: " << total_execution_time << " seconds\n";

    // Print the final position after the last segment
    std::cout << "Final position after " << total_time << " seconds: "
              << y0[0] << ", " << y0[1] << ", " << y0[2] << "\n";
}


void compare_ODE45_algorithm_with_segments(const std::string& algorithm_name,
                                           std::vector<std::vector<double>> (*algorithm_func)(
                                                   std::vector<double>(*)(double, const std::vector<double>&, double),
                                                   const std::vector<double>&,
                                                   const std::vector<double>&,
                                                   double, double, double),
                                           const std::vector<double>& r0, const std::vector<double>& v0,
                                           double total_time, double mu = 398600, double rtol = 1e-6, double atol = 1e-6) {
    double segment_duration = total_time / NUM_SEGMENTS; // Duration of each segment
    std::vector<double> y0 = r0; // Initial position
    y0.insert(y0.end(), v0.begin(), v0.end()); // Append initial velocity to form state vector

    double total_execution_time = 0.0;

    // Vectors to store execution times, final positions, and time points
    std::vector<double> execution_times;
    std::vector<std::vector<double>> final_positions;
    std::vector<double> time_points;

    for (int segment = 0; segment < NUM_SEGMENTS; ++segment) {
        double t_start = segment * segment_duration;
        double t_end = (segment + 1) * segment_duration;

        // Generate Gauss-Lobatto points for the current segment
        std::vector<double> segment_time_points = gaussLobattoPoints(NUM_GAUSS_LOBATTO_POINTS, t_start, t_end);

        auto start_time = std::chrono::high_resolution_clock::now();

        // Execute the ODE45 algorithm for all Gauss-Lobatto points in this segment
        std::vector<std::vector<double>> results = algorithm_func(satellite_motion, segment_time_points, y0, mu, rtol, atol);

        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end_time - start_time;

        // Accumulate total execution time for the segment
        double execution_time = elapsed.count();
        total_execution_time += execution_time;

        // Iterate through each result and store positions and times
        for (size_t i = 0; i < segment_time_points.size(); ++i) {
            time_points.push_back(segment_time_points[i]); // Store time point
            final_positions.push_back({results[i][0], results[i][1], results[i][2]}); // Store position
        }

        // Update y0 to the last state of this segment
        y0 = results.back();

        // Store execution time for the entire segment
        execution_times.push_back(execution_time);
    }

    // Save execution times to a CSV file (segment-level granularity)
    save_results_to_csv("Satellite", algorithm_name, time_points, execution_times);

    // Save final positions to a CSV file (time point-level granularity)
    save_final_positions_to_csv("Satellite", algorithm_name, time_points, final_positions);

    // Print total execution time
    std::cout << algorithm_name << " total execution time: " << total_execution_time << " seconds\n";

    // Print the final position after the last segment
    std::cout << "Final position after " << total_time << " seconds: "
              << y0[0] << ", " << y0[1] << ", " << y0[2] << "\n";
}


void compare_ODE78_algorithm_with_segments(double total_time1) {
    if (std::isnan(total_time1)) {  // Check if total_time1 is NaN
        std::cerr << "Error: total_time1 is NaN." << std::endl;
        return;
    }

    double total_time = total_time1;
    int num_segments = 16;                // Number of segments
    double segment_time = total_time / num_segments;
    int n_points = 32;                    // Number of Gauss-Lobatto points
    double mu = 398600;                   // Gravitational parameter
    double tol = 1e-16;                   // Tolerance
    double rtol = tol / 1000;
    double hmax = 0.01;
    double hmin = 1e-10;

    // Initial position and velocity vectors
    std::vector<double> r0 = {-19946.988367, 11684.423264, 43511.217135};  // Initial position
    std::vector<double> v0 = {-0.803367, -1.762325, 0.200044};             // Initial velocity
    std::vector<double> y0 = r0;
    y0.insert(y0.end(), v0.begin(), v0.end());

    // Vectors to store results
    std::vector<double> time_points;
    std::vector<std::vector<double>> final_positions;
    std::vector<double> execution_times;

    // Convert `c` from map<int, double> to vector<double>
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

        // Generate Gauss-Lobatto points for the current segment
        auto gauss_lobatto_tspan = gaussLobattoPoints(n_points, t_start, t_end);

        auto segment_start_time = std::chrono::high_resolution_clock::now();

        // Run ODE78 with Gauss-Lobatto points for the current segment
        auto segment_result = ode78(func, gauss_lobatto_tspan, y0, b, bh, c_vector, a, tol, rtol);

        auto segment_end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> segment_elapsed = segment_end_time - segment_start_time;

        // Store execution time for the segment
        execution_times.push_back(segment_elapsed.count());

        // Process results for each Gauss-Lobatto point
        for (const auto& result : segment_result) {
            double t_point = result[0];
            std::vector<double> current_position(result.begin() + 1, result.begin() + 4);  // Extract position

            time_points.push_back(t_point);            // Store time point
            final_positions.push_back(current_position); // Store final position
        }

        // Update `y0` for the next segment with the last state of this segment
        y0[0] = segment_result.back()[1];
        y0[1] = segment_result.back()[2];
        y0[2] = segment_result.back()[3];
        y0[3] = segment_result.back()[4];
        y0[4] = segment_result.back()[5];
        y0[5] = segment_result.back()[6];
    }

    // Record the total elapsed time
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    // Save execution times and positions to CSV files
    save_results_to_csv("Satellite", "ODE78", time_points, execution_times);
    save_final_positions_to_csv("Satellite", "ODE78", time_points, final_positions);

    // Print total execution time
    std::cout << "ODE78 Execution time: " << elapsed.count() << " seconds\n";

    // Print the final position
    std::cout << "ODE78 Final Position: X=" << y0[0] << " Y=" << y0[1] << " Z=" << y0[2] << std::endl;
}


void compare_ODE113_algorithm_with_segments(const std::string& algorithm_name,
                                            std::vector<std::vector<double>> (*algorithm_func)(
                                                    std::vector<double>(*)(double, const std::vector<double>&, double),
                                                    const std::vector<double>&,
                                                    const std::vector<double>&,
                                                    double, double, double, double),
                                            const std::vector<double>& r0, const std::vector<double>& v0,
                                            double total_time,
                                            double tol, double hmax, double hmin, double mu) {
    double segment_duration = total_time / NUM_SEGMENTS; // Duration of each segment
    std::vector<double> y0 = r0; // Initial position
    y0.insert(y0.end(), v0.begin(), v0.end()); // Append initial velocity to form state vector

    double total_execution_time = 0.0;

    // Vectors to store execution times, final positions, and time points
    std::vector<double> execution_times;
    std::vector<std::vector<double>> final_positions;
    std::vector<double> time_points;

    for (int segment = 0; segment < NUM_SEGMENTS; ++segment) {
        double t_start = segment * segment_duration;
        double t_end = (segment + 1) * segment_duration;

        // Generate Gauss-Lobatto points for the current segment
        std::vector<double> segment_time_points = gaussLobattoPoints(NUM_GAUSS_LOBATTO_POINTS, t_start, t_end);

        auto segment_start_time = std::chrono::high_resolution_clock::now();

        // Execute the ODE113 algorithm for all Gauss-Lobatto points in this segment
        std::vector<std::vector<double>> results = algorithm_func(satellite_motion, segment_time_points, y0, tol, hmax, hmin, mu);

        auto segment_end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> segment_elapsed = segment_end_time - segment_start_time;

        // Store execution time for the segment
        execution_times.push_back(segment_elapsed.count());

        // Process results for each Gauss-Lobatto point
        for (size_t i = 0; i < segment_time_points.size(); ++i) {
            time_points.push_back(segment_time_points[i]); // Store time point
            final_positions.push_back({results[i][0], results[i][1], results[i][2]}); // Store position
        }

        // Update y0 for the next segment with the last state of this segment
        y0 = results.back();
    }

    // Save execution times to a CSV file (segment-level granularity)
    save_results_to_csv("Satellite", algorithm_name, time_points, execution_times);

    // Save final positions to a CSV file (time point-level granularity)
    save_final_positions_to_csv("Satellite", algorithm_name, time_points, final_positions);

    // Print total execution time
    total_execution_time = std::accumulate(execution_times.begin(), execution_times.end(), 0.0);
    std::cout << algorithm_name << " total execution time: " << total_execution_time << " seconds\n";

    // Print the final position after the last segment
    std::cout << "Final position after " << total_time << " seconds: "
              << y0[0] << ", " << y0[1] << ", " << y0[2] << "\n";
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

    // File to save execution times
    std::ofstream exec_time_file("algorithm_execution_times.csv");
    if (!exec_time_file.is_open()) {
        std::cerr << "Error: Could not open file to save execution times.\n";
        return 1;
    }
    exec_time_file << "Algorithm,Execution Time (s)\n";

    // Timer for RK4 with segments
    auto start = std::chrono::high_resolution_clock::now();
    compare_RK4_algorithm_with_segments("RK4", rk4, r0, v0, total_time);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    exec_time_file << "RK4_with_segments," << elapsed.count() << "\n";
    std::cout << "RK4_with_segments execution time: " << elapsed.count() << " seconds\n";

    // Timer for RK8 with segments
    start = std::chrono::high_resolution_clock::now();
    double step_size = 0.1;  // Step size for RK8
    compare_RK8_algorithm_with_segments("RK8", RK8, r0, v0, total_time, step_size);
    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    exec_time_file << "RK8_with_segments," << elapsed.count() << "\n";
    std::cout << "RK8_with_segments execution time: " << elapsed.count() << " seconds\n";

    // Timer for ODE45 with segments
    start = std::chrono::high_resolution_clock::now();
    const double MU = 398600;    // Gravitational parameter
    const double RTOL = 1e-6;    // Relative tolerance
    const double ATOL = 1e-6;    // Absolute tolerance
    compare_ODE45_algorithm_with_segments("ODE45", ode45, r0, v0, total_time, MU, RTOL, ATOL);
    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    exec_time_file << "ODE45_with_segments," << elapsed.count() << "\n";
    std::cout << "ODE45_with_segments execution time: " << elapsed.count() << " seconds\n";

    // Timer for ODE78 with segments
    start = std::chrono::high_resolution_clock::now();
    double rtol = 1e-3;
    double atol = 1e-6;
    compare_ODE78_algorithm_with_segments(total_time);
    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    exec_time_file << "ODE78_with_segments," << elapsed.count() << "\n";
    std::cout << "ODE78_with_segments execution time: " << elapsed.count() << " seconds\n";

    // Timer for ODE113 with segments
    start = std::chrono::high_resolution_clock::now();
    compare_ODE113_algorithm_with_segments("ODE113", ODE113, r0, v0, total_time, tol, hmax, hmin, mu);
    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    exec_time_file << "ODE113_with_segments," << elapsed.count() << "\n";
    std::cout << "ODE113_with_segments execution time: " << elapsed.count() << " seconds\n";

    // Close the execution time file
    exec_time_file.close();

    return 0;
}
