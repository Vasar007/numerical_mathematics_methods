#pragma once

#include <iostream>
#include <cmath>

#include "poly34.h"
#include "utils.hpp"


namespace vv
{

// RESULT = 57.4846

/// ===== CONSTANT SECTION =====
double gauss_weight_moment_0(const double x)
{
    return -7.0 / 3.0 * pow(2.9 - x, 3.0 / 7.0);
}

double gauss_weight_moment_1(const double x)
{
    return -0.7 * (x + 203.0 / 30.0) * pow(2.9 - x, 3.0 / 7.0);
}

double gauss_weight_moment_2(const double x)
{
    return -0.4117647058823529 * (x * x + 4.06 * x + 27.47266666666667) * pow(2.9 - x, 3.0 / 7.0);
}

double gauss_weight_moment_3(const double x)
{
    return -0.2916666666666667 * (x * x * x + 3.582352941176471 * x * x + 14.54435294117648 * x
            + 98.4167882352941) * pow(2.9 - x, 3.0 / 7.0);
}

double gauss_weight_moment_4(const double x)
{
    return -0.2258064516129032 * (x * x * x * x + 3.383333333333334 * x * x * x
            + 12.12029411764707 * x * x + 49.20839411764709 * x + 332.9768001960785)
            * pow(2.9 - x, 3.0 / 7.0);
}

double gauss_weight_moment_5(const double x)
{
    return -(0.1842105263157895 * x * x * x * x * x + 0.6031409168081494 * x * x * x * x
             + 2.040626768534239 * x * x * x + 7.310245306102067 * x * x + 29.67959594277439 * x
             + 200.8319325461067) * pow(2.9 - x, 3.0 / 7.0);
}

template <class Function>
double gauss_calculate_weight_moment(Function&& gauss_weight_moment, const double a, const double b)
{
    return gauss_weight_moment(b) - gauss_weight_moment(a);
}

constexpr std::array<double, 3> gauss_coefficients_polynom(const std::array<double, 6>& u) noexcept
{
    std::array<double, 3> a
    {
        -((u.at(3) * u.at(3) * u.at(3) - 2.0 * u.at(2) * u.at(3) * u.at(4) + u.at(1) * u.at(4)
            * u.at(4) + u.at(2) * u.at(2) * u.at(5) - u.at(1) * u.at(3) * u.at(5)) / (u.at(2)
            * u.at(2) * u.at(2) - 2.0 * u.at(1) * u.at(2) * u.at(3) + u.at(0) * u.at(3) * u.at(3)
            + u.at(1) * u.at(1) * u.at(4) - u.at(0) * u.at(2) * u.at(4))),

        -((-u.at(2) * u.at(3) * u.at(3) + u.at(2) * u.at(2) * u.at(4) + u.at(1) * u.at(3) * u.at(4)
            - u.at(0) * u.at(4) * u.at(4) - u.at(1) * u.at(2) * u.at(5) + u.at(0) * u.at(3)
            * u.at(5)) / (u.at(2) * u.at(2) * u.at(2) - 2.0 * u.at(1) * u.at(2) * u.at(3) + u.at(0)
            * u.at(3) * u.at(3) + u.at(1) * u.at(1) * u.at(4) - u.at(0) * u.at(2) * u.at(4))),

        -((u.at(2) * u.at(2) * u.at(3) - u.at(1) * u.at(3) * u.at(3) - u.at(1) * u.at(2) * u.at(4)
            + u.at(0) * u.at(3) * u.at(4) + u.at(1) * u.at(1) * u.at(5) - u.at(0) * u.at(2)
            * u.at(5)) / (u.at(2) * u.at(2) * u.at(2) - 2.0 * u.at(1) * u.at(2) * u.at(3) + u.at(0)
            * u.at(3) * u.at(3) + u.at(1) * u.at(1) * u.at(4) - u.at(0) * u.at(2) * u.at(4)))
    };
    /*
    a0 -> -((u3^3 - 2 u2 u3 u4 + u1 u4^2 + u2^2 u5 - u1 u3 u5) / (u2^3 - 2 u1 u2 u3 + u0 u3^2 + u1^2 u4 - u0 u2 u4)),

    a1 -> -((-u2 u3^2 + u2^2 u4 + u1 u3 u4 - u0 u4^2 - u1 u2 u5 + u0 u3 u5)/(u2^3 - 2 u1 u2 u3 + u0 u3^2 + u1^2 u4 - u0 u2 u4)),

    a2 -> -((u2^2 u3 - u1 u3^2 - u1 u2 u4 + u0 u3 u4 + u1^2 u5 - u0 u2 u5)/(u2^3 - 2 u1 u2 u3 + u0 u3^2 + u1^2 u4 - u0 u2 u4))
    */
    return a;
};

constexpr std::array<double, 3> gauss_coefficients(const std::array<double, 3>& x,
                                                   const std::array<double, 6>& u) noexcept
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
template <class Function>
constexpr double gauss_integral(Function&& f, const double a, const double b,
                                const std::size_t segments)
{
    constexpr std::size_t degree = 3;
    constexpr std::array<double, 3> Xi{ -0.7745967, 0.0,       0.7745967 };
    constexpr std::array<double, 3> Ci{  0.5555556, 0.8888889, 0.5555556 };

    const double seg_step = (b - a) / segments;
    
    double f_val = 0.0;
    for (std::size_t i = 0; i < segments; ++i)
    {
        double seg_a = a + i * seg_step;
        double seg_b = seg_a + seg_step;
        double average = (seg_a + seg_b) / 2.0;

        double f_subval = 0.0;
        for (std::size_t j = 0; j < degree; ++j)
        {
            double Q = average + seg_step * Xi.at(j);
            f_subval += Ci.at(j) * f(Q);
        }
        f_val += seg_step * f_subval;
    }
    return f_val;
}


template <class Function>
constexpr double gauss_integral(Function&& f, const double a, const double b,
                                const std::size_t segments, const std::size_t degree)
{
    const double seg_step = (b - a) / segments;
    
    double f_val = 0.0;
    for (std::size_t i = 0; i < segments; ++i)
    {
        double seg_a = a + i * seg_step;
        double seg_b = seg_a + seg_step;

        if (i == segments - 1) seg_b = b;
        const std::array<double, 6> u
        {
            gauss_calculate_weight_moment(gauss_weight_moment_0, seg_a, seg_b),
            gauss_calculate_weight_moment(gauss_weight_moment_1, seg_a, seg_b),
            gauss_calculate_weight_moment(gauss_weight_moment_2, seg_a, seg_b),
            gauss_calculate_weight_moment(gauss_weight_moment_3, seg_a, seg_b),
            gauss_calculate_weight_moment(gauss_weight_moment_4, seg_a, seg_b),
            gauss_calculate_weight_moment(gauss_weight_moment_5, seg_a, seg_b)
        };
        const std::array<double, 3> coeff_a = gauss_coefficients_polynom(u);

        std::array<double, 3> Xi{};
        poly34::SolveP3(Xi.data(), coeff_a.at(2), coeff_a.at(1), coeff_a.at(0));
        
        utils::println(std::cout, Xi);

        const std::array<double, 3> Ci = gauss_coefficients(Xi, u);

        double f_subval = 0.0;
        for (std::size_t j = 0; j < degree; ++j)
        {
            f_subval += Ci.at(j) * f(Xi.at(j));
        }
        f_val += f_subval;
    }
    return f_val;
}


template <class Function>
double gauss_calculate_with_own_weight(Function&& f, const double a, const double b,
                                       const std::size_t segments)
{
    constexpr std::size_t degree = 3;

    return gauss_integral(std::forward<Function>(f), a, b, segments, degree);
}


/// ===== TEST SECTION =====
#if ENABLE_TESTS_GAUSS_INTEGRAL

TEST_METHOD(gauss_test_calculation_integral)
{
    std::cout << "Gauss method for integral:\n\n";

    constexpr std::size_t L = 2;
    constexpr std::size_t M = 6;
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

    const double integral = gauss_calculate_with_own_weight(function, a, b, N);
    const double S1 = integral;
    const double S2 = gauss_calculate_with_own_weight(function, a, b, N1);
    const double S3 = gauss_calculate_with_own_weight(function, a, b, N2);
    const double relative_error1 = (S2 - S1) / (1.0 - 1.0 / pow(L, M));
    const double relative_error2 = (S3 - S2) / (1.0 - 1.0 / pow(L, M));
    const double convergence_rate = -log((S3 - S2) / (S2 - S1)) / log(L);

    constexpr double true_result = 57.48462064655287;
    constexpr double methodical_error = 4.8909; // Calculated formula value by Wolfram.

    constexpr double guarantee_multiplier = 0.95;
    const double N_opt1 = std::ceil(N / pow(kEps / abs(relative_error1), 1.0 / M) / guarantee_multiplier);
    const double N_opt2 = std::ceil(N1 / pow(kEps / abs(relative_error2), 1.0 / M) / guarantee_multiplier);
    const std::size_t N_opt = static_cast<std::size_t>(std::ceil(std::max(N_opt1, N_opt2) / 2.0) * 2.0);

    std::cout << "True result:       " << true_result << "\n\n";

    std::cout << "Result with n1:    " << S1 << '\n';
    std::cout << "Result with n2:    " << S2 << '\n';
    std::cout << "Result with n3:    " << S3 << '\n';
    std::cout << "Result with 8:     " << gauss_calculate_with_own_weight(function, a, b, 8) << '\n';    
    std::cout << "Methodical error:  " << methodical_error << '\n';
    std::cout << "Relative error R1: " << relative_error1 << '\n';
    std::cout << "Relative error R2: " << relative_error2 << '\n';
    std::cout << "Convergence rate:  " << convergence_rate << '\n';
    std::cout << "Optimal step n1:   " << N_opt1 << '\n';
    std::cout << "Optimal step n2:   " << N_opt2;
    std::cout << "\n\n------------------------------------\n\n";
    
    std::cout << "Optimal step:      " << N_opt << '\n';
    const double integral_opt1 = gauss_calculate_with_own_weight(function, a, b, N_opt);
    const double integral_opt2 = gauss_calculate_with_own_weight(function, a, b, N_opt / 2);
    const double relative_error_opt = (integral_opt1 - integral_opt2) / (pow(L, M) - 1.0);
    
    std::cout << "Result with Nopt:  " <<  integral_opt1 << '\n';
    std::cout << "Result with Nopt / 2:  " <<  integral_opt2 << '\n';
    std::cout << "Relative error Ropt: " <<  relative_error_opt;
    std::cout << "\n\n------------------------------------\n\n";
}

#endif // ENABLE_TESTS_GAUSS_INTEGRAL

} // namespace vv

/*
 * Usage:
 * int main()
 * {
 *     constexpr std::size_t N = 100;
 *     constexpr double a = 1.0;
 *     constexpr double b = 2.0;
 *     constexpr double s = gauss_integral(function, a, b, N);
 *     std::cout << "I = " << s << '\n';
 * 
 *     return 0;
 * }
 */