// Copyright (C) 2018 Vasily Vasilyev (vasar007@yandex.ru)

#pragma once

#include <iostream>
#include <array>
#include <cassert>
#include <iterator>
#include <algorithm>

#include "cx_math.h"
#include "utils.hpp"


namespace vv
{
/// ===== CONSTANT SECTION =====
namespace
{

constexpr double epsilon   = 0.000001;
constexpr double epsilon_u = epsilon / 5.26;
constexpr double epsilon_v = epsilon / 4.04;
constexpr double epsilon_w = epsilon / 4.04;
constexpr double epsilon_z = epsilon / 4;

constexpr double PI_DIV_4  = static_cast<double>(cx::detail::pi()) / 4.0;
constexpr double PI_DIV_2  = static_cast<double>(cx::detail::pi()) / 2.0;
constexpr double PI3_DIV_4 = static_cast<double>(cx::detail::pi()) * 3.0 / 4.0;
constexpr double PI        = static_cast<double>(cx::detail::pi());
constexpr double PI5_DIV_4 = static_cast<double>(cx::detail::pi()) * 5.0 / 4.0;
constexpr double PI3_DIV_2 = static_cast<double>(cx::detail::pi()) * 3.0 / 2.0;
constexpr double PI7_DIV_4 = static_cast<double>(cx::detail::pi()) * 7.0 / 4.0;
constexpr double PI2       = static_cast<double>(cx::detail::pi()) * 2.0;

} // anonymous namespace

enum class function_type : uint8_t
{
    SIN,
    COS
};

struct always_false : std::false_type {};


/// ===== FUNCTION SECTION =====
constexpr double sin(const double arg, const double eps = epsilon)
{
    double result = arg;
    double n_member = arg;
    int n = 1;

    while (cx::abs(n_member) >= eps)
    {
        n_member *= -arg * arg / ((2.0 * n) * (2.0 * n + 1.0));
        result += n_member;
        ++n;
    }

    return result;
}


constexpr double cos(const double arg, const double eps = epsilon)
{
    double result = 1.0;
    double n_member = 1.0;
    int n = 1;

    while (cx::abs(n_member) >= eps)
    {
        n_member *= -arg * arg / ((2.0 * n) * (2.0 * n - 1.0));
        result += n_member;
        ++n;
    }

    return result;
}


constexpr double sqrt(const double arg, const double eps = epsilon)
{
    if (arg < 0.0)
    {
        return std::numeric_limits<double>::quiet_NaN();
    }

    double sum_n = 0.0;
    double sum_n1 = std::max(arg, 1.0);

    while (cx::abs(sum_n1 - sum_n) >= eps)
    {
        sum_n = sum_n1;
        sum_n1 = 0.5 * (sum_n + arg / sum_n);
    }

    return sum_n;
}


constexpr double exp(const double arg, const double eps = epsilon)
{
    double result = 1.0;
    double n_member = 1.0;
    int n = 1;

    while (cx::abs(n_member) >= eps)
    {
        n_member *= arg / n;
        result += n_member;
        ++n;
    }

    return result;
}


constexpr double reduce(const double arg, const double eps = epsilon)
{
    return (cx::abs(arg) - cx::floor(cx::abs(arg / eps)) * cx::abs(eps))
            * (arg >= 0.0 ? 1.0 : -1.0);
}


constexpr double deduce(const double arg, const function_type flag, const double eps = epsilon)
{
    if (arg > PI2) return deduce(reduce(arg, PI2), flag, eps);
    if (arg < -PI_DIV_4) return deduce(reduce(arg, PI_DIV_4), flag, eps);

    switch (flag)
    {
        case function_type::SIN:
            if (-PI_DIV_4 <= arg && arg <= PI_DIV_4)  return  vv::sin(arg,             eps);
            if ( PI_DIV_4 <= arg && arg <= PI_DIV_2)  return  vv::cos(PI_DIV_2 - arg,  eps);
            if ( PI_DIV_2 <= arg && arg <= PI3_DIV_4) return  vv::cos(arg - PI_DIV_2,  eps);
            if (PI3_DIV_4 <= arg && arg <= PI)        return  vv::sin(PI - arg,        eps);
            if (       PI <= arg && arg <= PI5_DIV_4) return -vv::sin(arg - PI,        eps);
            if (PI5_DIV_4 <= arg && arg <= PI3_DIV_2) return -vv::cos(PI3_DIV_2 - arg, eps);
            if (PI3_DIV_2 <= arg && arg <= PI7_DIV_4) return -vv::cos(arg - PI3_DIV_2, eps);
            if (PI7_DIV_4 <= arg && arg <= PI2)       return  vv::sin(arg - PI2,       eps);
            break;
        
        case function_type::COS:
            if (-PI_DIV_4 <= arg && arg <= PI_DIV_4)  return  vv::cos(arg,             eps);
            if ( PI_DIV_4 <= arg && arg <= PI_DIV_2)  return  vv::sin(PI_DIV_2 - arg,  eps);
            if ( PI_DIV_2 <= arg && arg <= PI3_DIV_4) return -vv::sin(arg - PI_DIV_2,  eps);
            if (PI3_DIV_4 <= arg && arg <= PI)        return -vv::cos(PI - arg,        eps);
            if (       PI <= arg && arg <= PI5_DIV_4) return -vv::cos(arg - PI,        eps);
            if (PI5_DIV_4 <= arg && arg <= PI3_DIV_2) return -vv::sin(PI3_DIV_2 - arg, eps);
            if (PI3_DIV_2 <= arg && arg <= PI7_DIV_4) return  vv::sin(arg - PI3_DIV_2, eps);
            if (PI7_DIV_4 <= arg && arg <= PI2)       return  vv::cos(arg - PI2,       eps);
            break;

        default:
            assert(false);
            break;
    }

    return 0.0;
}

constexpr size_t NUMBER_OF_ARRAYS = 8;

template <size_t N>
constexpr std::array<std::array<double, N>, NUMBER_OF_ARRAYS>
    do_work(const std::array<double, N>& data)
{
    std::array<double, N> my_results{};
    std::array<double, N> expected_results{};
    std::array<double, N> my_u_values{};
    std::array<double, N> my_v_values{};
    std::array<double, N> my_w_values{};
    std::array<double, N> u_values{};
    std::array<double, N> v_values{};
    std::array<double, N> w_values{};

    for (size_t i = 0; i < data.size(); ++i)
    {
        double u = vv::sqrt(1.0 + data.at(i) * data.at(i), epsilon_u);
        //double v = vv::sin(3.0 * data.at(i) + 0.1, epsilon_v);
        //double w = vv::cos(2.0 * data.at(i) + 0.3, epsilon_w);
        double v = vv::deduce(3.0 * data.at(i) + 0.1, vv::function_type::SIN, epsilon_v);
        double w = vv::deduce(2.0 * data.at(i) + 0.3, vv::function_type::COS, epsilon_w);
        my_results.at(i) = u * (v + w);
        my_u_values.at(i) = u; my_v_values.at(i) = v; my_w_values.at(i) = w;

        u = cx::sqrt(1.0 + data.at(i) * data.at(i));
        v = cx::sin(3.0 * data.at(i) + 0.1);
        w = cx::cos(2.0 * data.at(i) + 0.3);
        expected_results.at(i) = u * (v + w);
        u_values.at(i) = u; v_values.at(i) = v; w_values.at(i) = w;
    }

    return { my_results, expected_results, my_u_values, my_v_values, my_w_values, u_values,
             v_values, w_values };
}


/// ===== TEST SECTION =====
#if ENABLE_TESTS_FUNCTIONS

TEST_METHOD(functions_test)
{
    std::cout << "Constexpr math functions:\n\n";

    constexpr std::size_t num_of_values = 11;
    constexpr std::array<double, num_of_values> values = 
    { 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3 };

    constexpr auto results = vv::do_work(values);
    
    std::cout.precision(10);
    std::cout << std::fixed;
    std::cout << "NUMBER: |";
    std::cout << "RESULT:       | ";
    std::cout << "REAL RESULT: | ";
    std::cout << "DELTA: \n";
    for (std::size_t i = 0; i < num_of_values; ++i)
    {
        std::cout << i <<  (i != 10 ? "       | " : "      | ");
        std::cout << results.at(0).at(i) << " | ";
        std::cout << results.at(1).at(i) << " | ";
        std::cout << cx::abs(results.at(1).at(i) - results.at(0).at(i)) << '\n';
    }
    std::cout << "\n\n";
}

#endif // ENABLE_TESTS_FUNCTIONS

} // namesapce vv
