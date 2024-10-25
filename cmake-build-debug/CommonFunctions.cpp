// SatelliteUtils.cpp
#include "CommonFunctions.h"
#include <cmath>
#include <numeric>
#include <stdexcept>

double vectorNorm(const std::vector<double>& vec) {
    double sum = 0;
    for (double val : vec) {
        sum += val * val;
    }
    return std::sqrt(sum);
}

std::vector<double> vectorAdd(const std::vector<double>& a, const std::vector<double>& b) {
    std::vector<double> result(a.size());
    for (size_t i = 0; i < a.size(); ++i) {
        result[i] = a[i] + b[i];
    }
    return result;
}

std::vector<double> scalarMultiply(const std::vector<double>& vec, double scalar) {
    std::vector<double> result(vec.size());
    for (size_t i = 0; i < vec.size(); ++i) {
        result[i] = vec[i] * scalar;
    }
    return result;
}

std::vector<double> satelliteMotion(double t, const std::vector<double>& y) {
    std::vector<double> r = {y[0], y[1], y[2]};  // Position vector
    std::vector<double> v = {y[3], y[4], y[5]};  // Velocity vector

    double r_norm = vectorNorm(r);  // Magnitude of the position vector

    // Acceleration vector
    std::vector<double> a = {-mu / std::pow(r_norm, 3) * r[0],
                             -mu / std::pow(r_norm, 3) * r[1],
                             -mu / std::pow(r_norm, 3) * r[2]};

    std::vector<double> dydt(6);
    dydt[0] = v[0]; dydt[1] = v[1]; dydt[2] = v[2];  // Velocity components
    dydt[3] = a[0]; dydt[4] = a[1]; dydt[5] = a[2];  // Acceleration components

    return dydt;
}

std::vector<double> gaussLobattoPoints(int n, double a, double b) {
    std::vector<double> points(n);

    for (int i = 0; i < n; ++i) {
        points[i] = -std::cos(M_PI * i / (n - 1));
    }

    std::vector<double> scaled_points(n);
    for (int i = 0; i < n; ++i) {
        scaled_points[i] = 0.5 * (b - a) * (points[i] + 1) + a;
    }

    return scaled_points;
}

// Overload the + operator for vector addition (vector + vector)
std::vector<double> operator+(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) {
        throw std::invalid_argument("Vectors must be of the same length");
    }
    std::vector<double> result(a.size());
    for (size_t i = 0; i < a.size(); ++i) {
        result[i] = a[i] + b[i];
    }
    return result;
}

// Overload the + operator for scalar and vector addition (scalar + vector)
std::vector<double> operator+(double scalar, const std::vector<double>& v) {
    std::vector<double> result(v.size());
    for (size_t i = 0; i < v.size(); ++i) {
        result[i] = scalar + v[i];
    }
    return result;
}

// Overload the + operator for vector and scalar addition (vector + scalar)
std::vector<double> operator+(const std::vector<double>& v, double scalar) {
    return scalar + v;  // Use the scalar + vector overload
}

// Overload the - operator for vector subtraction (vector - vector)
std::vector<double> operator-(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) {
        throw std::invalid_argument("Vectors must be of the same length");
    }
    std::vector<double> result(a.size());
    for (size_t i = 0; i < a.size(); ++i) {
        result[i] = a[i] - b[i];
    }
    return result;
}

// Overload the - operator for scalar and vector subtraction (scalar - vector)
std::vector<double> operator-(double scalar, const std::vector<double>& v) {
    std::vector<double> result(v.size());
    for (size_t i = 0; i < v.size(); ++i) {
        result[i] = scalar - v[i];
    }
    return result;
}

// Overload the - operator for vector and scalar subtraction (vector - scalar)
std::vector<double> operator-(const std::vector<double>& v, double scalar) {
    std::vector<double> result(v.size());
    for (size_t i = 0; i < v.size(); ++i) {
        result[i] = v[i] - scalar;
    }
    return result;
}


// Overload the * operator for scalar multiplication (vector * scalar)
std::vector<double> operator*(const std::vector<double>& v, double scalar) {
    std::vector<double> result(v.size());
    for (size_t i = 0; i < v.size(); ++i) {
        result[i] = v[i] * scalar;
    }
    return result;
}

// Overload the * operator for scalar multiplication (scalar * vector)
std::vector<double> operator*(double scalar, const std::vector<double>& v) {
    return v * scalar;  // Use the other overloaded operator
}

// Overload the / operator for vector and scalar division (vector / scalar)
std::vector<double> operator/(const std::vector<double>& v, double scalar) {
    if (scalar == 0.0) {
        throw std::invalid_argument("Division by zero");
    }
    std::vector<double> result(v.size());
    for (size_t i = 0; i < v.size(); ++i) {
        result[i] = v[i] / scalar;
    }
    return result;
}

// Function to calculate satellite motion
std::vector<double> satelliteMotion(double t, const std::vector<double>& y, double mu) {
    std::vector<double> r(y.begin(), y.begin() + 3);  // Position vector
    std::vector<double> v(y.begin() + 3, y.end());    // Velocity vector

    double r_norm = vectorNorm(r);  // Magnitude of position vector

    // Acceleration vector
    std::vector<double> a(3);
    for (int i = 0; i < 3; ++i) {
        a[i] = -mu / std::pow(r_norm, 3) * r[i];
    }

    // Combine velocity and acceleration into a single vector
    std::vector<double> dydt(6);
    for (int i = 0; i < 3; ++i) {
        dydt[i] = v[i];     // Velocity
        dydt[i + 3] = a[i]; // Acceleration
    }

    return dydt;
}

// Satellite motion model differential equation
std::vector<double> satellite_motion(double t, const std::vector<double>& y, double mu) {
    std::vector<double> r(y.begin(), y.begin() + 3);  // Position vector
    std::vector<double> v(y.begin() + 3, y.end());    // Velocity vector
    double r_norm = std::sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);

    // Acceleration
    std::vector<double> a(3);
    for (int i = 0; i < 3; ++i) {
        a[i] = -mu / (r_norm * r_norm * r_norm) * r[i];
    }

    // Combine velocity and acceleration
    std::vector<double> dydt(6);
    for (int i = 0; i < 3; ++i) {
        dydt[i] = v[i];  // Velocity components
        dydt[i + 3] = a[i];  // Acceleration components
    }

    return dydt;
}

// Helper function for dot product (replacing the previous approach)
std::vector<double> dot_product(const std::vector<double>& weights, const std::vector<std::vector<double>>& k) {
    std::vector<double> result(k[0].size(), 0.0);
    for (size_t i = 0; i < weights.size(); ++i) {
        for (size_t j = 0; j < result.size(); ++j) {
            result[j] += weights[i] * k[i][j];
        }
    }
    return result;
}

// Helper function to scale a vector by a scalar
std::vector<double> scalar_multiply(double scalar, const std::vector<double>& vec) {
    std::vector<double> result(vec.size());
    for (size_t i = 0; i < vec.size(); ++i) {
        result[i] = scalar * vec[i];
    }
    return result;
}

// Norm function to calculate vector norm (like np.linalg.norm in Python)
double vector_norm(const std::vector<double>& vec) {
    return std::sqrt(std::inner_product(vec.begin(), vec.end(), vec.begin(), 0.0));
}
// Helper function to perform element-wise vector addition
std::vector<double> vector_add(const std::vector<double>& v1, const std::vector<double>& v2) {
    std::vector<double> result(v1.size());
    for (size_t i = 0; i < v1.size(); ++i) {
        result[i] = v1[i] + v2[i];
    }
    return result;
}