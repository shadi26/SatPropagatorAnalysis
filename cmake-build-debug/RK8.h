// RK8.h
#ifndef RK8_H
#define RK8_H

#include <vector>

// RK8 function declaration
std::vector<std::vector<double>> RK8(
    std::vector<double>(*f)(double, const std::vector<double>&, double),
    const std::vector<double>& t_gauss_lobatto, const std::vector<double>& Y0, double h = 0.1);

#endif // RK8_H
