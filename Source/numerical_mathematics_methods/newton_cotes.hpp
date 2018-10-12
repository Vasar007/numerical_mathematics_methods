// Copyright (C) 2018 Vasily Vasilyev (vasar007@yandex.ru)

#pragma once

#include <iostream>
#include <cmath>
#include <array>

#include "utils.hpp"


namespace vv
{

/// ===== CONSTANT SECTION =====
double nc_weight_moment_0(const double x)
{
    return -7.0 / 3.0 * pow(2.9 - x, 3.0 / 7.0);
}

double nc_weight_moment_1(const double x)
{
    return -0.7 * (x + 203.0 / 30.0) * pow(2.9 - x, 3.0 / 7.0);
}

double nc_weight_moment_2(const double x)
{
    return -0.4117647058823529 * (x * x + 4.06 * x + 27.47266666666667) * pow(2.9 - x, 3.0 / 7.0);
}

template <class Function>
double nc_calculate_weight_moment(Function&& nc_weight_moment, const double a, const double b)
{
    return nc_weight_moment(b) - nc_weight_moment(a);
}

constexpr std::array<double, 3> nc_coefficients(const std::array<double, 3>& x,
                                                const std::array<double, 3>& u) noexcept
{
    std::array<double, 3> A
    {
        -(-u.at(2) + u.at(1) * x.at(1) + u.at(1) * x.at(2) - u.at(0) * x.at(1) * x.at(2))
          / (( x.at(0) - x.at(1)) * ( x.at(0) - x.at(2))),

        -(-u.at(2) + u.at(1) * x.at(0) + u.at(1) * x.at(2) - u.at(0) * x.at(0) * x.at(2))
          / ((-x.at(0) + x.at(1)) * ( x.at(1) - x.at(2))),

        -( u.at(2) - u.at(1) * x.at(0) - u.at(1) * x.at(1) + u.at(0) * x.at(0) * x.at(1))
          / (( x.at(1) - x.at(2)) * (-x.at(0) + x.at(2)))
    };
    return A;
}


/// ===== FUNCTION SECTION =====
template <class Function, std::size_t Size>
constexpr double newton_cotes(Function&& f, const double a, const double b,
                              const std::size_t segments, const std::size_t degree,
                              const std::size_t divisor, const std::array<double, Size>& Ci)
{
    const double seg_step = (b - a) / segments;
    
    double f_val = 0.0;
    for (std::size_t i = 0; i < segments; ++i)
    {
        double seg_a = a + seg_step * i;
        double seg_b = seg_a + seg_step;

        double f_subval = 0.0;
        for (std::size_t j = 0; j < degree + 1; ++j)
        {
            f_subval += Ci.at(j) * f(seg_a + j * (seg_b - seg_a) / degree);
        }
        f_val += (seg_b - seg_a) / divisor * f_subval;
    }
    return f_val;
}


template <class Function>
constexpr double newton_cotes(Function&& f, const double a, const double b,
                              const std::size_t segments, const std::size_t degree)
{
    const double seg_step = (b - a) / segments;
    
    double f_val = 0.0;
    for (std::size_t i = 0; i < segments; ++i)
    {
        double seg_a = a + seg_step * i;
        double seg_b = seg_a + seg_step;

        if (i == segments - 1) seg_b = b;
        const std::array<double, 3> Xi{ seg_a, (seg_a + seg_b) / 2.0, seg_b };
        const std::array<double, 3> u
        {
            nc_calculate_weight_moment(nc_weight_moment_0, seg_a, seg_b),
            nc_calculate_weight_moment(nc_weight_moment_1, seg_a, seg_b),
            nc_calculate_weight_moment(nc_weight_moment_2, seg_a, seg_b)
        };
        const std::array<double, 3> Ci = nc_coefficients(Xi, u);

        double f_subval = 0.0;
        for (std::size_t j = 0; j < degree + 1; ++j)
        {
            f_subval += Ci.at(j) * f(Xi.at(j));
        }
        f_val += f_subval;
    }
    return f_val;
}


// Simpsons Rule
template <class Function>
double simpsons(Function&& f, const double a, const double b, const std::size_t segments)
{ 
    constexpr std::size_t degree = 2;
    constexpr std::size_t divisor = 6;
    constexpr std::array<double, 3> Ci{ 1.0, 4.0, 1.0 };

    return newton_cotes(std::forward<Function>(f), a, b, segments, degree, divisor, Ci);
}


// Simpsons 3/8 Rule
template <class Function>
double simpsons_3_8(Function&& f, const double a, const double b, const std::size_t segments)
{
    constexpr std::size_t degree = 3;
    constexpr std::size_t divisor = 8;
    constexpr std::array<double, 4> Ci{ 1.0, 3.0, 3.0, 1.0 };

    return newton_cotes(std::forward<Function>(f), a, b, segments, degree, divisor, Ci);
}


// Booles Rule
template <class Function>
double booles(Function&& f, const double a, const double b, const std::size_t segments)
{ 
    constexpr std::size_t degree = 4;
    constexpr std::size_t divisor = 90;
    constexpr std::array<double, 5> Ci{ 7.0, 32.0, 12.0, 32.0, 7.0 };

    return newton_cotes(std::forward<Function>(f), a, b, segments, degree, divisor, Ci);
}


// Own weight.
template <class Function>
double nc_calculate_with_own_weight(Function&& f, const double a, const double b,
                                    const std::size_t segments)
{
    constexpr std::size_t degree = 2;

    return newton_cotes(std::forward<Function>(f), a, b, segments, degree);
}


constexpr double const_test(const double) noexcept
{
    return 42.0;
}

double pythagorian_test(const double t)
{
    return sin(t) * sin(t) + cos(t) * cos(t);
}


/// ===== TEST SECTION =====
#if ENABLE_TESTS_NEWTON_COTES_INTEGRAL

TEST_METHOD(newton_cotes_test_calculation_integral)
{
    std::cout << "Newton-Cotes method for integrals:\n\n";

    constexpr std::size_t L = 2;
    constexpr std::size_t M = 3;
    constexpr std::size_t N = 1;
    constexpr std::size_t N1 = N * L;
    constexpr std::size_t N2 = N1 * L;

    constexpr double a = 1.8;
    constexpr double b = 2.9;

    const auto function = [](const double x)
    {
        return 4.0 * cos(2.5 * x) * exp(4.0 * x / 7.0) + 2.5 * sin(5.5 * x) * exp(-3.0 * x / 5.0)
             + 4.3 * x;
    };

    double integral = nc_calculate_with_own_weight(function, a, b, N);
    const double S1 = integral;
    const double S2 = nc_calculate_with_own_weight(function, a, b, N1);
    const double S3 = nc_calculate_with_own_weight(function, a, b, N2);
    const double relative_error1 = (S2 - S1) / (1.0 - 1.0 / pow(L, M));
    const double relative_error2 = (S3 - S2) / (1.0 - 1.0 / pow(L, M));
    double convergence_rate = -log((S3 - S2) / (S2 - S1)) / log(L);

    constexpr double true_result = 57.48462064655287;
    constexpr double methodical_error = 4.8909; // Calculated formula value by Wolfram.

    constexpr double guarantee_multiplier = 0.95;
    double N_opt1 = std::ceil(N / pow(kEps / abs(relative_error1), 1.0 / M) / guarantee_multiplier);
    double N_opt2 = std::ceil(N1 / pow(kEps / abs(relative_error2), 1.0 / M) / guarantee_multiplier);
    const std::size_t N_opt = static_cast<std::size_t>(std::ceil(std::max(N_opt1, N_opt2) / 2.0) * 2.0);

    std::cout << "True result:       " << true_result << "\n\n";

    std::cout << "Result with n1:    " << S1 << '\n';
    std::cout << "Result with n2:    " << S2 << '\n';
    std::cout << "Result with n3:    " << S3 << '\n';
    std::cout << "Methodical error:  " << methodical_error << '\n';
    std::cout << "Relative error R1: " << relative_error1 << '\n';
    std::cout << "Relative error R2: " << relative_error2 << '\n';
    std::cout << "Convergence rate:  " << convergence_rate << '\n';
    std::cout << "Optimal step n1:   " << N_opt1 << '\n';
    std::cout << "Optimal step n2:   " << N_opt2;
    std::cout << "\n\n------------------------------------\n\n";
    
    const double integral_opt1 = nc_calculate_with_own_weight(function, a, b, N_opt);
    const double integral_opt2 = nc_calculate_with_own_weight(function, a, b, N_opt / 2);
    const double relative_error_opt = (integral_opt1 - integral_opt2) / (pow(L, M) - 1.0);

    std::cout << "Optimal step:      " << N_opt << '\n'; 
    std::cout << "Result with Nopt:  " <<  integral_opt1 << '\n';
    std::cout << "Relative error Ropt: " <<  relative_error_opt << "\n\n";

    double I1;
    double I2;
    double I3;
    double relative_error = 1.0;
    std::size_t Mi = 1;
    while (abs(relative_error) >= kEps)
    {
        I1 = nc_calculate_with_own_weight(function, a, b, Mi);
        I2 = nc_calculate_with_own_weight(function, a, b, Mi * L);
        I3 = nc_calculate_with_own_weight(function, a, b, Mi * L * L);
        convergence_rate = -log((I3 - I2) / (I2 - I1)) / log(L);
        std::cout << "Convergence rate:  " << convergence_rate << '\n';

        relative_error = (I3 - I2) / (pow(L, M) - 1.0);
        Mi *= 2;
    }

    std::cout << "Optimal step:      " << Mi * L << '\n'; 
    std::cout << "Result with Nopt:  " <<  I2 << '\n';
    std::cout << "Relative error Ropt: " <<  relative_error;
    std::cout << "\n\n------------------------------------\n\n";
}

#endif // ENABLE_TESTS_NEWTON_COTES_INTEGRAL

} // namespace vv

/*
 * Usage:
 * int main()
 * {
 *     constexpr std::size_t N = 1000;
 *     constexpr double a = 1.8;
 *     constexpr double b = 2.9;
 *     std::cout << "Simpsons:     " << simpsons(pythagorian_test, a, b, N) << '\n';
 *     std::cout << "Simpsons 3/8: " << simpsons_3_8(pythagorian_test, a, b, N) << '\n';
 *     std::cout << "Booles:       " << booles(pythagorian_test, a, b, N) << '\n';
 *     return 0;
 * }
 */
