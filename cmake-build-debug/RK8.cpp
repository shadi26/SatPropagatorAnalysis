#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include "CommonFunctions.h"
#include "RK8.h"

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

// Constants
const double MU_EARTH = 398600.0;  // Gravitational parameter for Earth in km^3/s^2



// Dormand-Prince 8th-order Butcher table (c, a, b coefficients)
const std::vector<double> c_DP8 = {0, 1.0 / 18, 1.0 / 12, 1.0 / 8, 5.0 / 16, 3.0 / 8, 59.0 / 400, 93.0 / 200,
                                   5490023248.0 / 9719169821.0, 13.0 / 20, 1201146811.0 / 1299019798.0, 1.0, 1.0};

const std::vector<std::vector<double>> a_DP8 = {
    {},
    {1.0 / 18},
    {1.0 / 48, 1.0 / 16},
    {1.0 / 32, 0, 3.0 / 32},
    {5.0 / 16, 0, -75.0 / 64, 75.0 / 64},
    {3.0 / 80, 0, 0, 3.0 / 16, 3.0 / 20},
    {29443841.0 / 614563906, 0, 0, 77736538.0 / 692538347, -28693883.0 / 1125000000, 23124283.0 / 1800000000},
    {16016141.0 / 946692911, 0, 0, 61564180.0 / 158732637, 22789713.0 / 633445777, 545815736.0 / 2771057229, -180193667.0 / 1043307555},
    {39632708.0 / 573591083, 0, 0, -433636366.0 / 683701615, -421739975.0 / 2616292301, 100302831.0 / 723423059,
     790204164.0 / 839813087, 800635310.0 / 3783071287},
    {246121993.0 / 1340847787, 0, 0, -37695042795.0 / 15268766246, -309121744.0 / 1061227803, -12992083.0 / 490766935,
     6005943493.0 / 2108947869, 393006217.0 / 1396673457, 123872331.0 / 1001029789},
    {-1028468189.0 / 846180014, 0, 0, 8478235783.0 / 508512852, 1311729495.0 / 1432422823, -10304129995.0 / 1701304382,
     -48777925059.0 / 3047939560, 15336726248.0 / 1032824649, -45442868181.0 / 3398467696, 3065993473.0 / 597172653},
    {185892177.0 / 718116043, 0, 0, -3185094517.0 / 667107341, -477755414.0 / 1098053517, -703635378.0 / 230739211,
     5731566787.0 / 1027545527, 5232866602.0 / 850066563, -4093664535.0 / 808688257, 3962137247.0 / 1805957418, 65686358.0 / 487910083},
    {403863854.0 / 491063109, 0, 0, -5068492393.0 / 434740067, -411421997.0 / 543043805, 652783627.0 / 914296604,
     11173962825.0 / 925320556, -13158990841.0 / 6184727034, 3936647629.0 / 1978049680, -160528059.0 / 685178525, 248638103.0 / 1413531060}
};

const std::vector<double> b_DP8 = {14005451.0 / 335480064, 0, 0, 0, 0, -59238493.0 / 1068277825, 181606767.0 / 758867731,
                                   561292985.0 / 797845732, -1041891430.0 / 1371343529, 760417239.0 / 1151165299,
                                   118820643.0 / 751138087, -528747749.0 / 2220607170, 1.0 / 4};

// Dormand-Prince 8th-order method using Gauss-Lobatto points
std::vector<std::vector<double>> RK8(
    std::vector<double>(*f)(double, const std::vector<double>&, double),
    const std::vector<double>& t_gauss_lobatto, const std::vector<double>& Y0, double h ) {

    std::vector<double> tout = t_gauss_lobatto;
    std::vector<std::vector<double>> yout(tout.size(), std::vector<double>(Y0.size()));
    yout[0] = Y0;

    std::vector<double> y = Y0;

    for (size_t i = 1; i < tout.size(); ++i) {
        double h_step = tout[i] - tout[i - 1];
        std::vector<std::vector<double>> k(c_DP8.size(), std::vector<double>(Y0.size()));

        // Calculate k values using Butcher table coefficients
        for (size_t j = 0; j < c_DP8.size(); ++j) {
            std::vector<double> y_temp;
            if (j == 0) {
                y_temp = y;
            } else {
                std::vector<double> sum(Y0.size(), 0.0);
                for (size_t l = 0; l < std::min(j, a_DP8[j].size()); ++l) {
                    sum = vectorAdd(sum, scalarMultiply(k[l], a_DP8[j][l]));
                }
                y_temp = vectorAdd(y, scalarMultiply(sum, h_step));
            }
            k[j] = f(tout[i - 1] + c_DP8[j] * h_step, y_temp, MU_EARTH);
        }

        // Update state vector y
        std::vector<double> sum(Y0.size(), 0.0);
        for (size_t j = 0; j < k.size(); ++j) {
            sum = vectorAdd(sum, scalarMultiply(k[j], b_DP8[j]));
        }
        y = vectorAdd(y, scalarMultiply(sum, h_step));
        yout[i] = y;
    }

    return yout;
}


