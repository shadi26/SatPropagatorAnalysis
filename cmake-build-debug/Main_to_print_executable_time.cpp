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
const int NUM_GAUSS_LOBATTO_POINTS = 16;  // Number of points per segment, can be adjusted


void save_results_to_csv(const std::string& satellite_name, const std::string& algorithm_name,
                         const std::vector<double>& time_points, const std::vector<double>& position_differences) {
    std::string directory = "results";
    std::filesystem::create_directory(directory);

    std::string filename = directory + "/" + satellite_name + "_" + algorithm_name + "_results.csv";
    std::ofstream file(filename);

    if (file.is_open()) {
        file << "Time (s), Error (km)\n";
        for (size_t i = 0; i < time_points.size(); ++i) {
            file << std::fixed << std::setprecision(13)
                 << time_points[i] << "," << position_differences[i] << "\n";
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
            file << std::fixed << std::setprecision(13)
                 << time_points[i] << "," << positions[i][0] << ","
                 << positions[i][1] << "," << positions[i][2] << "\n";
        }
        std::cout << "Final positions saved to " << filename << std::endl;
    } else {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }

    file.close();
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

void compare_RK4_algorithm_with_segments(
        const std::string& algorithm_name,
        std::vector<std::vector<double>> (*algorithm_func)(
                std::vector<double>(*)(double, const std::vector<double>&),
                const std::vector<double>&,
                const std::vector<double>&),
        const std::vector<double>& r0,
        const std::vector<double>& v0,
        double total_time,
        double orbital_period) {

    // Calculate the number of orbits and segment duration
    int num_orbits = static_cast<int>(total_time / orbital_period); // Number of orbits
    int num_segments = 8; // Number of segments per orbit
    double segment_duration = orbital_period / num_segments; // Duration of each segment

    std::vector<double> y0 = r0; // Initial position
    y0.insert(y0.end(), v0.begin(), v0.end()); // Append initial velocity to form state vector

    double total_execution_time = 0.0;
    int total_points = 0; // Variable to count total Gauss-Lobatto points used

    // Vectors to store execution times, final positions, and time points
    std::vector<double> execution_times;
    std::vector<std::vector<double>> final_positions;
    std::vector<double> time_points;

    // Outer loop for orbits
    for (int orbit = 0; orbit < num_orbits; ++orbit) {
        for (int segment = 0; segment < num_segments; ++segment) {
            double t_start = orbit * orbital_period + segment * segment_duration;
            double t_end = t_start + segment_duration;

            // Generate Gauss-Lobatto points for the current segment
            std::vector<double> segment_time_points = gaussLobattoPoints(NUM_GAUSS_LOBATTO_POINTS, t_start, t_end);

            // Increment total_points by the number of Gauss-Lobatto points
            total_points += segment_time_points.size();

            // Start measuring execution time for this segment
            auto start_time = std::chrono::high_resolution_clock::now();

            // Execute the RK4 algorithm for all Gauss-Lobatto points in this segment
            std::vector<std::vector<double>> results = algorithm_func(satelliteMotion, segment_time_points, y0);

            // Stop measuring execution time for this segment
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

            // Store the execution time for the segment
            execution_times.push_back(execution_time);
        }
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

    // Print the total number of Gauss-Lobatto points used
    std::cout << "Total number of Gauss-Lobatto points used: " << total_points << "\n";
}






void compare_RK8_algorithm_with_segments(
        const std::string& algorithm_name,
        std::vector<std::vector<double>> (*algorithm_func)(
                std::vector<double>(*)(double, const std::vector<double>&, double),
                const std::vector<double>&,
                const std::vector<double>&,
                double),
        const std::vector<double>& r0,
        const std::vector<double>& v0,
        double total_time,
        double orbital_period,
        double h = 0.1) {


    // Calculate the number of orbits and segment duration
    int num_orbits = static_cast<int>(total_time / orbital_period); // Number of orbits
    int num_segments = 8; // Number of segments per orbit
    double segment_duration = orbital_period / num_segments; // Duration of each segment

    std::vector<double> y0 = r0; // Initial position
    y0.insert(y0.end(), v0.begin(), v0.end()); // Append initial velocity to form state vector

    double total_execution_time = 0.0;
    int total_points = 0; // Variable to count total Gauss-Lobatto points used

    // Vectors to store execution times, final positions, and time points
    std::vector<double> execution_times;
    std::vector<std::vector<double>> final_positions;
    std::vector<double> time_points;

    // Outer loop for orbits
    for (int orbit = 0; orbit < num_orbits; ++orbit) {
        for (int segment = 0; segment < num_segments; ++segment) {
            double t_start = orbit * orbital_period + segment * segment_duration;
            double t_end = t_start + segment_duration;

            // Generate Gauss-Lobatto points for the current segment
            std::vector<double> segment_time_points = gaussLobattoPoints(NUM_GAUSS_LOBATTO_POINTS, t_start, t_end);

            // Increment total_points by the number of Gauss-Lobatto points
            total_points += segment_time_points.size();

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

            // Store execution time for the segment
            execution_times.push_back(execution_time);
        }
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

    // Print the total number of Gauss-Lobatto points used
    std::cout << "Total number of Gauss-Lobatto points used: " << total_points << "\n";
}



void compare_ODE45_algorithm_with_segments(
        const std::string& algorithm_name,
        std::vector<std::vector<double>> (*algorithm_func)(
                std::vector<double>(*)(double, const std::vector<double>&, double),
                const std::vector<double>&,
                const std::vector<double>&,
                double, double, double),
        const std::vector<double>& r0,
        const std::vector<double>& v0,
        double total_time,
        double orbital_period,
        double mu,
        double rtol,
        double atol) {

    // Calculate the number of orbits and segment duration
    int num_orbits = static_cast<int>(total_time / orbital_period); // Number of orbits
    int num_segments = 8; // Number of segments per orbit
    double segment_duration = orbital_period / num_segments; // Duration of each segment

    std::vector<double> y0 = r0; // Initial position
    y0.insert(y0.end(), v0.begin(), v0.end()); // Append initial velocity to form state vector

    double total_execution_time = 0.0;
    int total_points = 0; // Variable to count total Gauss-Lobatto points used

    // Vectors to store execution times, final positions, and time points
    std::vector<double> execution_times;
    std::vector<std::vector<double>> final_positions;
    std::vector<double> time_points;

    // Outer loop for orbits
    for (int orbit = 0; orbit < num_orbits; ++orbit) {
        for (int segment = 0; segment < num_segments; ++segment) {
            double t_start = orbit * orbital_period + segment * segment_duration;
            double t_end = t_start + segment_duration;

            // Generate Gauss-Lobatto points for the current segment
            std::vector<double> segment_time_points = gaussLobattoPoints(NUM_GAUSS_LOBATTO_POINTS, t_start, t_end);

            // Increment total_points by the number of Gauss-Lobatto points
            total_points += segment_time_points.size();

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

            // Store execution time for the segment
            execution_times.push_back(execution_time);
        }
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

    // Print the total number of Gauss-Lobatto points used
    std::cout << "Total number of Gauss-Lobatto points used: " << total_points << "\n";
}



void compare_ODE78_algorithm_with_segments(
        std::vector<double> r0,
        std::vector<double> v0,
        double total_time1,
        double orbital_period) {

    if (std::isnan(total_time1)) {  // Check if total_time1 is NaN
        std::cerr << "Error: total_time1 is NaN." << std::endl;
        return;
    }

    double total_time = total_time1;
    int num_orbits = static_cast<int>(total_time / orbital_period); // Number of orbits
    int num_segments = 8;                // Number of segments per orbit
    double segment_time = orbital_period / num_segments; // Duration of each segment
    int n_points = 16;                    // Number of Gauss-Lobatto points per segment
    double mu = 398600;                   // Gravitational parameter
    double tol = 1e-16;                   // Tolerance
    double rtol = tol / 1000;
    double hmax = 0.01;
    double hmin = 1e-10;

    std::vector<double> y0 = r0;
    y0.insert(y0.end(), v0.begin(), v0.end());

    // Vectors to store results
    std::vector<double> time_points;
    std::vector<std::vector<double>> final_positions;
    std::vector<double> execution_times;

    // Convert c from map<int, double> to vector<double>
    std::vector<double> c_vector(c.size());
    for (const auto& pair : c) {
        c_vector[pair.first - 1] = pair.second;  // Convert map to vector, adjusting index
    }

    auto func = [&mu](double t, const std::vector<double>& y) { return satellite_motion(t, y, mu); };

    // Record the start time
    auto start = std::chrono::high_resolution_clock::now();

    // Outer loop for orbits
    for (int orbit = 0; orbit < num_orbits; ++orbit) {
        for (int segment = 0; segment < num_segments; ++segment) {
            double t_start = orbit * orbital_period + segment * segment_time;
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

            // Update y0 for the next segment with the last state of this segment
            y0[0] = segment_result.back()[1];
            y0[1] = segment_result.back()[2];
            y0[2] = segment_result.back()[3];
            y0[3] = segment_result.back()[4];
            y0[4] = segment_result.back()[5];
            y0[5] = segment_result.back()[6];
        }
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
    std::cout << "ODE78 Final Position: X=" << y0[0] << ", Y=" << y0[1] << ", Z=" << y0[2] << std::endl;

    // Print the total number of points used
    int total_points = num_orbits * num_segments * n_points;
    std::cout << "Total number of Gauss-Lobatto points used: " << total_points << "\n";
}





void compare_ODE113_algorithm_with_segments(
        const std::string& algorithm_name,
        const std::vector<double>& r0,
        const std::vector<double>& v0,
        double tol,
        double hmax,
        double hmin,
        double mu,
        double orbital_period,
        double total_time) {
    std::cout << "afs3 hyane hon" ;
    int num_orbits = static_cast<int>(total_time / orbital_period);
    int num_segments = NUM_SEGMENTS;
    double segment_duration = orbital_period / num_segments;

    std::vector<double> y0 = r0;
    y0.insert(y0.end(), v0.begin(), v0.end());

    double total_execution_time = 0.0;
    int total_points = 0;

    std::vector<double> execution_times;
    std::vector<std::vector<double>> final_positions;
    std::vector<double> time_points;

    for (int orbit = 0; orbit < num_orbits; ++orbit) {
        for (int segment = 0; segment < num_segments; ++segment) {
            double t_start = orbit * orbital_period + segment * segment_duration;
            double t_end = t_start + segment_duration;

            // Generate Gauss-Lobatto points
            std::vector<double> segment_time_points = gaussLobattoPoints(NUM_GAUSS_LOBATTO_POINTS, t_start, t_end);
            total_points += segment_time_points.size();

            auto segment_start_time = std::chrono::high_resolution_clock::now();

            // Execute ODE113 for the current segment
            ODE113Result results = ODE113(
                satellite_motion,
                segment_time_points,
                y0,
                mu,
                tol,
                tol, // Absolute tolerance
                hmax,
                hmin
            );

            auto segment_end_time = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> segment_elapsed = segment_end_time - segment_start_time;

            // Store execution time
            execution_times.push_back(segment_elapsed.count());
            total_execution_time += segment_elapsed.count();

            // Store results for each Gauss-Lobatto point
            for (size_t i = 0; i < results.time.size(); ++i) {
                time_points.push_back(results.time[i]);
                final_positions.push_back({results.values[i][0], results.values[i][1], results.values[i][2]});
            }

            // Update y0 for the next segment
            y0 = results.values.back();
        }
    }

    // Save results to CSV
    save_results_to_csv("Satellite", algorithm_name, time_points, execution_times);
    save_final_positions_to_csv("Satellite", algorithm_name, time_points, final_positions);

    // Print statistics
    std::cout << algorithm_name << " total execution time: " << total_execution_time << " seconds\n";
    std::cout << "Final position after " << total_time << " seconds: "
              << y0[0] << ", " << y0[1] << ", " << y0[2] << "\n";
    std::cout << "Total number of points used: " << total_points << "\n";
}





int main() {
    // Initial conditions
    std::vector<double> r0 = {  -19946.988367233487, 11684.423263957098, 43511.217134673760 };  // Initial position vector
    std::vector<double> v0 = { -0.803367495870, -1.762325494392, 0.200043952127  }; // Initial velocity vector

    // Define parameters for ODE113
    double tol = 1e-16;   // Tolerance
    double hmax = 0.01;   // Maximum step size
    double hmin = 1e-10;  // Minimum step size
    double mu = 398600;  // Gravitational parameter for Earth
    double total_time = 1000000;  // Total simulation time in seconds
    // Compute the orbital period
    double orbital_period =     57438.77;

    // File to save execution times
    std::ofstream exec_time_file("algorithm_execution_times.csv");
    if (!exec_time_file.is_open()) {
        std::cerr << "Error: Could not open file to save execution times.\n";
        return 1;
    }
    exec_time_file << "Algorithm,Execution Time (s)\n";

    // Timer for RK4 with segments
    auto start = std::chrono::high_resolution_clock::now();
    compare_RK4_algorithm_with_segments("RK4", rk4, r0, v0, total_time, orbital_period);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    exec_time_file << "RK4_with_segments," << elapsed.count() << "\n";
    std::cout << "RK4_with_segments execution time: " << elapsed.count() << " seconds\n";

    // Timer for RK8 with segments
    start = std::chrono::high_resolution_clock::now();
    double step_size = 0.1;  // Step size for RK8
    compare_RK8_algorithm_with_segments(
    "RK8",
    RK8,
    r0,
    v0,
    total_time,
    orbital_period, // Pass the orbital period here
    0.1             // Step size
);
    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    exec_time_file << "RK8_with_segments," << elapsed.count() << "\n";
    std::cout << "RK8_with_segments execution time: " << elapsed.count() << " seconds\n";

    // Timer for ODE45 with segments
    start = std::chrono::high_resolution_clock::now();
    const double MU = 398600;    // Gravitational parameter
    const double RTOL = 1e-16;    // Relative tolerance
    const double ATOL = 1e-16;    // Absolute tolerance
    compare_ODE45_algorithm_with_segments(
    "ODE45",
    ode45,
    r0,
    v0,
    total_time,
    orbital_period,
    398600,  // mu
    1e-6,    // rtol
    1e-6     // atol
);
    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    exec_time_file << "ODE45_with_segments," << elapsed.count() << "\n";
    std::cout << "ODE45_with_segments execution time: " << elapsed.count() << " seconds\n";

    // Timer for ODE78 with segments
    start = std::chrono::high_resolution_clock::now();
    double rtol = 1e-3;
    double atol = 1e-6;
    compare_ODE78_algorithm_with_segments(r0, v0, total_time, orbital_period);
    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    exec_time_file << "ODE78_with_segments," << elapsed.count() << "\n";
    std::cout << "ODE78_with_segments execution time: " << elapsed.count() << " seconds\n";

    // Timer for ODE113 with segments
    start = std::chrono::high_resolution_clock::now();
    compare_ODE113_algorithm_with_segments(
    "ODE113", r0, v0, 1e-6, 100, 1e-3, mu, orbital_period, 500000
);
    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    exec_time_file << "ODE113_with_segments," << elapsed.count() << "\n";
    std::cout << "ODE113_with_segments execution time: " << elapsed.count() << " seconds\n";

    // Close the execution time file
    exec_time_file.close();

    return 0;
}
