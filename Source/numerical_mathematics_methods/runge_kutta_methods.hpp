#pragma once

#include <iostream>
#include <cmath>
#include <vector>
#include <tuple>
#include <cassert>

#include "utils.hpp"
#include "matrix.hpp"


namespace vv
{
/// ===== CONSTANT SECTION =====
namespace
{

constexpr double A = -1.0;
constexpr double B = -2.0;
constexpr double C = -1.0;
constexpr std::size_t kNumber_F = 4;

} // anonymous namespace

detail::matrix::container<double> f(const double x)
{
    detail::matrix::container<double> result
    {
        std::exp(sin(x * x)),
        std::exp(B * std::sin(x * x)),
        C * std::sin(x * x) + A,
        std::cos(x * x)
    };
    return result;
}

std::tuple<matrix<double>, matrix<double>>
calc_true(const double a, const double b, const double h, const matrix<double>& x0,
          const matrix<double>& y0)
{
    const std::size_t n = static_cast<std::size_t>((b - a) / h);

    matrix<double> Y(n, kNumber_F);
    matrix<double> X(n, 1);

    X(0, 0) = x0(0, 0);
    Y(0) = y0(0);
    for (std::size_t i = 1; i < n; ++i)
    {
        X(i, 0) = a + i * h;
        Y(i) = f(X(i, 0));
    }
    return { X, Y };
}

detail::matrix::container<double>
F(const double x, const detail::matrix::container<double>& y)
{
    const detail::matrix::container<double> result
    {
        2.0 * x * std::pow(y.at(1), 1.0 / B) * y.at(3),
        2.0 * B * x * std::exp(B / C * (y.at(2) - A)) * y.at(3),
        2.0 * C * x * y.at(3),
        -2.0 * x * std::log(y.at(0))
    };
    return result;
}


/// ===== FUNCTION SECTION =====
namespace detail
{

void runga_kutta_2(const double c2, const double a21, const double b1,
                   const double b2, const double h,
                   double& X, matrix::container<double>& K1, matrix::container<double>& K2,
                   matrix::container<double>& Y1, const matrix::container<double>& Y)
{
    X += h;
    K1 = h * F(X, Y);
    K2 = h * F(X + c2 * h, Y + a21 * K1);
    Y1 = Y + (b1 * K1 + b2 * K2);
}


void runga_kutta_2(const double c2, const double a, const double a21, const double b1,
                   const double b2, const double h, const std::size_t i, ::vv::matrix<double>& X,
                   ::vv::matrix<double>& K1, ::vv::matrix<double>& K2, ::vv::matrix<double>& Y)
{
    X(i, 0) = a + i * h;
    K1(i) = h * F(X(i - 1, 0), Y(i - 1));
    K2(i) = h * F(X(i - 1, 0) + c2 * h, Y(i - 1) + a21 * K1(i));
    Y(i) = Y(i - 1) + (b1 * K1(i) + b2 * K2(i));
}


void runga_kutta_4(const double h, double& X, matrix::container<double>& K1,
                   matrix::container<double>& K2, matrix::container<double>& K3,
                   matrix::container<double>& K4, matrix::container<double>& Y1,
                   const matrix::container<double>& Y)
{
    X += h;
    K1 = h * F(X, Y);
    K2 = h * F(X + h / 3.0, Y + K1 / 3.0);
    K3 = h * F(X + 2.0 * h / 3.0, Y - K1 / 3.0 + K2);
    K4 = h * F(X + h, Y + K1 - K2 + K3);
    Y1 = Y + (K1 + 3.0 * K2 + 3.0 * K3 + K4) / 8.0;
}

void runga_kutta_4(const double a, const double h, const std::size_t i, ::vv::matrix<double>& X,
                   ::vv::matrix<double>& K1, ::vv::matrix<double>& K2, ::vv::matrix<double>& K3,
                   ::vv::matrix<double>& K4, ::vv::matrix<double>& Y)
{
    X(i, 0) = a + i * h;
    K1(i) = h * F(X(i - 1, 0), Y(i - 1));
    K2(i) = h * F(X(i - 1, 0) + h / 3.0, Y(i - 1) + K1(i) / 3.0);
    K3(i) = h * F(X(i - 1, 0) + 2.0 * h / 3.0, Y(i - 1) - K1(i) / 3.0 + K2(i));
    K4(i) = h * F(X(i - 1, 0) + h, Y(i - 1) + K1(i) - K2(i) + K3(i));
    Y(i) = Y(i - 1) + (K1(i) + 3.0 * K2(i) + 3.0 * K3(i) + K4(i)) / 8.0;
}

} // namespace detail


std::tuple<matrix<double>, matrix<double>>
runge_kutta_2(const double c2, const double a, const double b, const double h,
              const matrix<double>& x0, const matrix<double>& y0)
{
    assert(c2 != 0.0);
    const std::size_t n = static_cast<std::size_t>((b - a) / h);

    const double a21 = c2;
    const double b2 = 1.0 / (2.0 * c2);
    const double b1 = 1.0 - b2;
    //std::cout << "Info: " << a21 << ' ' << b1 << ' ' << b2 << '\n';

    matrix<double> K1(n + 1, kNumber_F);
    matrix<double> K2(n + 1, kNumber_F);
    matrix<double> Y(n + 1, kNumber_F);
    matrix<double> X(n + 1, 1);

    X(0, 0) = x0(0, 0);
    Y(0) = y0(0);
    for (std::size_t i = 1; i <= n; ++i)
    {
        detail::runga_kutta_2(c2, a, a21, b1, b2, h, i, X, K1, K2, Y);
    }
    return { X, Y };
}


std::tuple<matrix<double>, matrix<double>>
runge_kutta_4(const double a, const double b, const double h, const matrix<double>& x0,
              const matrix<double>& y0)
{
    const std::size_t n = static_cast<std::size_t>((b - a) / h);

    matrix<double> K1(n + 1, kNumber_F);
    matrix<double> K2(n + 1, kNumber_F);
    matrix<double> K3(n + 1, kNumber_F);
    matrix<double> K4(n + 1, kNumber_F);
    matrix<double> Y(n + 1, kNumber_F);
    matrix<double> X(n + 1, 1);

    X(0, 0) = x0(0, 0);
    Y(0) = y0(0);
    for (std::size_t i = 1; i <= n; ++i)
    {
        detail::runga_kutta_4(a, h, i, X, K1, K2, K3, K4, Y);
    }
    return { X, Y };
}


double norm(const detail::matrix::container<double>& x)
{
    double result = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i)
    {
        result += x.at(i) * x.at(i);
    }

    return sqrt(result);
}


matrix<double> extract(const matrix<double>& x, const std::size_t n)
{
    if (n >= x.get_columns_number())
    {
        return matrix<double>{};
    }

    matrix<double> result(x.get_rows_number(), 1);
    for (std::size_t i = 0; i < x.get_rows_number(); ++i)
    {
        result(i, 0) = x(i, n);
    }
    return result;
}


std::vector<std::pair<double, double>> process_data(const matrix<double>& x,
                                                    const matrix<double>& y,
                                                    const matrix<double>& y1)
{
    std::vector<std::pair<double, double>> result;
    if (x.get_rows_number() != y.get_rows_number() || y.get_rows_number() != y1.get_rows_number())
    {
        std::cout << "\nERROR: Matrices have different size! Cannot process data.\n";
        return result;
    }
    
    result.reserve(x.get_rows_number());
    for (std::size_t i = 0; i < x.get_rows_number(); ++i)
    {
        result.emplace_back(x(i, 0), norm(y(i) - y1(i)));
    }

    return result;
}

matrix<double> data_x_true;
matrix<double> data_y_true;

matrix<double> data_x_true_opp;
matrix<double> data_y_true_opp;

std::vector<std::pair<double, double>> graphic_data;
std::vector<std::pair<double, double>> graphic_data_opp;

std::vector<std::pair<double, double>> tol_graphic_data;
std::vector<std::pair<double, double>> tol_graphic_data_opp;

std::vector<std::pair<double, double>> tol_graphic_data_auto;
std::vector<std::pair<double, double>> tol_graphic_data_auto_opp;

matrix<double> graphic_data_x;
matrix<double> graphic_data_y;

matrix<double> graphic_data_x_opp;
matrix<double> graphic_data_y_opp;

std::vector<std::pair<double, double>> quality;
std::vector<std::pair<double, double>> quality_opp;

std::vector<std::pair<double, double>> h_acc;
std::vector<std::pair<double, double>> h_not_acc;

std::vector<std::pair<double, double>> h_acc_opp;
std::vector<std::pair<double, double>> h_not_acc_opp;

std::vector<std::pair<double, double>> count_F;
std::vector<std::pair<double, double>> count_F_opp;


double calculate(const double c2, const double a, const double b, const double h,
                 const matrix<double>& x0, const matrix<double>& y0)
{
    const std::size_t n = static_cast<std::size_t>((b - a) / h);
    const std::size_t n2 = n * 2;

    const auto [X1, Y1] = runge_kutta_2(c2, a, b, h, x0, y0);
    const auto [X1_2, Y1_2] = runge_kutta_2(c2, a, b, h / 2.0, x0, y0);
    
    const detail::matrix::container<double> yn = Y1(n);
    const detail::matrix::container<double> y2n = Y1_2(n2);

    constexpr std::size_t p = 2;
    const double R2n = norm(y2n - yn) / (std::pow(2, p) - 1.0);

    std::cout << "Step = " << h << "\nR2n = " << R2n << '\n';

    return R2n;
}


std::tuple<matrix<double>, matrix<double>>
calculate_last(const double c2, const double a, const double b, const double h,
               const matrix<double>& x0, const matrix<double>& y0)
{
    const std::size_t n = static_cast<std::size_t>((b - a) / h);
    const std::size_t n2 = n / 2;

    const auto[X1, Y1] = runge_kutta_2(c2, a, b, h, x0, y0);
    const auto[X1_2, Y1_2] = runge_kutta_2(c2, a, b, h * 2.0, x0, y0);

    const detail::matrix::container<double> y2n = Y1(n);
    const detail::matrix::container<double> yn = Y1_2(n2);

    constexpr std::size_t p = 2;
    const double R2n = norm(y2n - yn) / (std::pow(2.0, p) - 1.0);
    std::cout << "Step = " << h << "\nR2n = " << R2n << '\n';

    return { X1, Y1 };
}


double calculate_opp(const double a, const double b, const double h,
                     const matrix<double>& x0, const matrix<double>& y0)
{
    const std::size_t n = static_cast<std::size_t>((b - a) / h);
    const std::size_t n2 = n * 2;

    const auto [X2, Y2] = runge_kutta_4(a, b, h, x0, y0);
    const auto [X2_2, Y2_2] = runge_kutta_4(a, b, h / 2.0, x0, y0);

    const detail::matrix::container<double> yn_opp = Y2(n);
    const detail::matrix::container<double> y2n_opp = Y2_2(n2);

    constexpr std::size_t p = 4;
    const double R2n = norm(y2n_opp - yn_opp) / (std::pow(2, p) - 1.0);

    std::cout << "Opponent step = " << h << "\nR2n = " << R2n << '\n';

    return R2n;
}

std::tuple<matrix<double>, matrix<double>>
calculate_last_opp(const double a, const double b, const double h, const matrix<double>& x0,
                   const matrix<double>& y0)
{
    const std::size_t n = static_cast<std::size_t>((b - a) / h);
    const std::size_t n2 = n / 2;

    const auto [X1, Y1] = runge_kutta_4(a, b, h, x0, y0);
    const auto [X1_2, Y1_2] = runge_kutta_4(a, b, h * 2.0, x0, y0);

    const detail::matrix::container<double> y2n_opp = Y1(n);
    const detail::matrix::container<double> yn_opp = Y1_2(n2);

    constexpr std::size_t p = 4;
    const double R2n = norm(y2n_opp - yn_opp) / (std::pow(2.0, p) - 1.0);
    std::cout << "Opponent step = " << h << "\nR2n = " << R2n << '\n';

    return { X1, Y1 };
}


double calculate_h0(const double rtol, const double atol, const std::size_t p,
                    const matrix<double>& x0, const matrix<double>& y0)
{
    const auto abs_compare = [](const double a1, const double a2)
                             { return std::abs(a1) < std::abs(a2); };
    const double tol = norm(y0(0)) * rtol + atol;
    // Step 1:
    detail::matrix::container<double> f0 = F(x0(0, 0), y0(0));

    // Step 2:
    double max_x = *std::max_element(std::begin(x0(0)), std::end(x0(0)), abs_compare);
    double delta = std::pow(1.0 / (max_x + 1.0), p + 1) + std::pow(norm(f0), p + 1);

    // Step 3:
    double h1 = std::pow(tol / delta, 1.0 / (p + 1.0));

    // Step 4:
    double ux = x0(0, 0) + h1;
    detail::matrix::container<double> uy = y0(0) + h1 * F(x0(0, 0), y0(0));

    // Step 5:
    detail::matrix::container<double> f0_ = F(ux, uy);
    double delta_ = std::pow(1.0 / ux, p + 1) + std::pow(norm(f0_), p + 1);
    double h1_ = std::pow(tol / delta_, 1.0 / (p + 1.0));
    return std::min(h1, h1_);
}


std::tuple<matrix<double>, matrix<double>>
calculate_with_h_auto(const double c2, const double a, const double b, const matrix<double>& x0,
                      const matrix<double>& y0, const double rtol, const double atol)
{
    assert(c2 != 0.0);

    const double a21 = c2;
    const double b2 = 1.0 / (2.0 * c2);
    const double b1 = 1.0 - b2;

    constexpr std::size_t p = 2;

    const double h = calculate_h0(rtol, atol, p, x0, y0);
    const std::size_t n = static_cast<std::size_t>(b - a);

    matrix<double> Y(n, kNumber_F, mtx::mode::dynamically_expandable);
    matrix<double> X(n, 1, mtx::mode::dynamically_expandable);

    detail::matrix::container<double> K1_1;
    detail::matrix::container<double> K2_1;
    detail::matrix::container<double> Y_1 = y0(0);
    double X_1 = x0(0, 0);

    detail::matrix::container<double> Y_12;
    double X_12;

    detail::matrix::container<double> K1_2;
    detail::matrix::container<double> K2_2;
    detail::matrix::container<double> Y_2 = y0(0);
    double X_2 = x0(0, 0);

    X(0, 0) = x0(0, 0);
    Y(0) = y0(0);

    std::vector<double> h_vec{ h };
    h_vec.reserve(n);
    quality.reserve(n);
    h_acc.reserve(n);
    h_not_acc.reserve(n);

    std::size_t i = 1;
    while (X(i - 1, 0) + h_vec.at(i - 1) <= b)
    {
        std::cout << "Step = " << (i - 1) << " \th = " << h_vec.at(i - 1) << '\t';
        detail::runga_kutta_2(c2, a21, b1, b2, h_vec.at(i - 1), X_1, K1_1, K2_1, Y_1, Y(i - 1));

        detail::runga_kutta_2(c2, a21, b1, b2, h_vec.at(i - 1) / 2.0, X_2, K1_2, K2_2, Y_2, Y(i - 1));
        X_12 = X_2;
        Y_12 = Y_2;
        detail::runga_kutta_2(c2, a21, b1, b2, h_vec.at(i - 1) / 2.0, X_2, K1_2, K2_2, Y_2, Y_2);

        double max1 = 0.0;
        double max2 = 0.0;
        double max3 = 0.0;
        for (std::size_t j = 0; j < kNumber_F; ++j)
        {
            max1 += Y(i - 1, j) * Y(i - 1, j);
            max2 += Y_1.at(j) * Y_1.at(j);
            max3 += Y_2.at(j) * Y_2.at(j);
        }
        const double tol = rtol * std::sqrt(std::max({ max1, max2, max3 })) + atol;

        const double rn = norm(Y_2 - Y_1) / (1.0 -  1.0 / std::pow(2, p));

        if (rn > tol * std::pow(2.0, p))
        {
            h_not_acc.emplace_back(X_1, h_vec.at(i - 1));

            h_vec.at(i - 1) /= 2.0;
            Y_1 = Y_12;

            std::cout << "Chosen 1 way.\n";
        }
        else if (rn > tol)
        {
            h_not_acc.emplace_back(X_1, h_vec.at(i - 1));

            h_vec.emplace_back(h_vec.at(i - 1) / 2.0);
            Y(i) = Y_2;
            X(i, 0) = X_1;
            ++i;

            quality.emplace_back(X_1, norm(f(X_1) - Y_2) / rn);
            std::cout << "Chosen 2 way.\n";
        }
        else if (rn >= tol / std::pow(2.0, p + 1))
        {
            h_acc.emplace_back(X_1, h_vec.at(i - 1));

            h_vec.emplace_back(h_vec.at(i - 1));
            Y(i) = Y_1;
            X(i, 0) = X_1;
            ++i;

            quality.emplace_back(X_1, norm(f(X_1) - Y_1) / rn);
            std::cout << "Chosen 3 way.\n";
        }
        else
        {
            h_not_acc.emplace_back(X_1, h_vec.at(i - 1));

            const double h_max = *std::max_element(std::begin(h_vec), std::end(h_vec));
            if (h_max != h_vec.at(i - 1))
            {
                h_vec.push_back(std::min(h_vec.at(i - 1) * 2.0, h_max));
            }
            else
            {
                h_vec.push_back(h_vec.at(i - 1) * 2.0);
            }
            Y(i) = Y_1;
            X(i, 0) = X_1;
            ++i;

            quality.emplace_back(X_1, norm(f(X_1) - Y_1) / rn);
            std::cout << "Chosen 4 way.\n";
        }
    }
    constexpr double addressing = 3.0;
    constexpr double usage_F_rkm_2 = 2.0;
    count_F.emplace_back(-std::log10(rtol), std::log10((i - 1) * addressing * usage_F_rkm_2));

    return { X, Y };
}


std::tuple<matrix<double>, matrix<double>>
calculate_with_h_auto_opp(const double a, const double b, const matrix<double>& x0,
                          const matrix<double>& y0, const double rtol, const double atol)
{
    constexpr std::size_t p = 4;

    const double h = calculate_h0(rtol, atol, p, x0, y0);
    const std::size_t n = static_cast<std::size_t>(b - a);

    matrix<double> Y(n, kNumber_F, mtx::mode::dynamically_expandable);
    matrix<double> X(n, 1, mtx::mode::dynamically_expandable);

    detail::matrix::container<double> K1_1;
    detail::matrix::container<double> K2_1;
    detail::matrix::container<double> K3_1;
    detail::matrix::container<double> K4_1;
    detail::matrix::container<double> Y_1 = y0(0);
    double X_1 = x0(0, 0);

    detail::matrix::container<double> Y_12;
    double X_12;

    detail::matrix::container<double> K1_2;
    detail::matrix::container<double> K2_2;
    detail::matrix::container<double> K3_2;
    detail::matrix::container<double> K4_2;
    detail::matrix::container<double> Y_2 = y0(0);
    double X_2 = x0(0, 0);

    X(0, 0) = x0(0, 0);
    Y(0) = y0(0);

    std::vector<double> h_vec{ h };
    h_vec.reserve(n);
    quality_opp.reserve(n);
    h_acc_opp.reserve(n);
    h_not_acc_opp.reserve(n);

    std::size_t i = 1;
    while (X(i - 1, 0) + h_vec.at(i - 1) <= b)
    {
        std::cout << "Step = " << (i - 1) << " \th = " << h_vec.at(i - 1) << '\t';
        detail::runga_kutta_4(h_vec.at(i - 1), X_1, K1_1, K2_1, K3_1, K4_1, Y_1, Y(i - 1));

        detail::runga_kutta_4(h_vec.at(i - 1) / 2.0, X_2, K1_2, K2_2, K3_2, K4_2, Y_2, Y(i - 1));
        X_12 = X_2;
        Y_12 = Y_2;
        detail::runga_kutta_4(h_vec.at(i - 1) / 2.0, X_2, K1_2, K2_2, K3_2, K4_2, Y_2, Y_2);

        double max1 = 0.0;
        double max2 = 0.0;
        double max3 = 0.0;
        for (std::size_t j = 0; j < kNumber_F; ++j)
        {
            max1 += Y(i - 1, j) * Y(i - 1, j);
            max2 += Y_1.at(j) * Y_1.at(j);
            max3 += Y_2.at(j) * Y_2.at(j);
        }
        const double tol = rtol * std::sqrt(std::max({ max1, max2, max3 })) + atol;

        const double rn = norm(Y_2 - Y_1) / (1.0 -  1.0 / std::pow(2, p));

        if (rn > tol * std::pow(2.0, p))
        {
            h_not_acc_opp.emplace_back(X_1, h_vec.at(i - 1));

            h_vec.at(i - 1) /= 2.0;
            Y_1 = Y_12;

            std::cout << "Chosen 1 way.\n";
        }
        else if (rn > tol)
        {
            h_not_acc_opp.emplace_back(X_1, h_vec.at(i - 1));

            h_vec.emplace_back(h_vec.at(i - 1) / 2.0);
            Y(i) = Y_2;
            X(i, 0) = X_1;
            ++i;

            quality_opp.emplace_back(X_1, norm(f(X_1) - Y_2) / rn);
            std::cout << "Chosen 2 way.\n";
        }
        else if (rn >= tol / std::pow(2.0, p + 1))
        {
            h_acc_opp.emplace_back(X_1, h_vec.at(i - 1));

            h_vec.emplace_back(h_vec.at(i - 1));
            Y(i) = Y_1;
            X(i, 0) = X_1;
            ++i;

            quality_opp.emplace_back(X_1, norm(f(X_1) - Y_1) / rn);
            std::cout << "Chosen 3 way.\n";
        }
        else
        {
            h_not_acc_opp.emplace_back(X_1, h_vec.at(i - 1));

            const double h_max = *std::max_element(std::begin(h_vec), std::end(h_vec));
            if (h_max != h_vec.at(i - 1))
            {
                h_vec.push_back(std::min(h_vec.at(i - 1) * 2.0, h_max));
            }
            else
            {
                h_vec.push_back(h_vec.at(i - 1) * 2.0);
            }
            Y(i) = Y_1;
            X(i, 0) = X_1;
            ++i;

            quality_opp.emplace_back(X_1, norm(f(X_1) - Y_1) / rn);
            std::cout << "Chosen 4 way.\n";
        }
    }
    constexpr double addressing = 3.0;
    constexpr double usage_F_rkm_4 = 4.0;
    count_F_opp.emplace_back(-std::log10(rtol), std::log10((i - 1) * addressing * usage_F_rkm_4));

    return { X, Y };
}


/// ===== TEST SECTION =====
#if ENABLE_TESTS_RUNGE_KUTTA

TEST_METHOD(runge_kutta_calculation)
{
    std::cout << "Runge-Kutta method with p=2:\n\n";

    constexpr double a = 0.0;
    constexpr double b = 5.0;
    constexpr double c2 = 0.65;
    const matrix<double> x0{ 0.0, 0.0, 0.0, 0.0 };
    const matrix<double> y0{ 1.0, 1.0, A,   1.0 };

    constexpr std::size_t p = 2;
    constexpr double tol = 1e-5;
    
    double R2n;
    double R2n_prev;
    double h;
    for (std::size_t i = 5; i <= 13; ++i)
    {
        h = 1.0 / std::pow(2.0, i);
        R2n = calculate(c2, a, b, h, x0, y0);

        if (i == 5)
        {
            R2n_prev = R2n;
            continue;
        }
        graphic_data.emplace_back(std::log2(h), std::log2(R2n));
        R2n_prev = R2n;
    }
    //sconst double h_opt = h * std::pow(tol / std::abs(R2n), 1.0 / p) / 2.0;

    //std::cout << "\nh_opt = " << h_opt << "\n\n";
    
    //const auto [X1, Y1] = calculate_last(c2, a, b, h_opt, x0, y0);
    //const auto [X, Y] = calc_true(a, b, h_opt, x0, y0);
    
    //tol_graphic_data = process_data(X, Y, Y1);

    std::cout << "\nCalculations are finished.\n\n";
}


TEST_METHOD(runge_kutta_calculation_opp)
{
    std::cout << "Runge-Kutta method with p=4:\n\n";

    constexpr double a = 0.0;
    constexpr double b = 5.0;
    const matrix<double> x0{ 0.0, 0.0, 0.0, 0.0 };
    const matrix<double> y0{ 1.0, 1.0, A,   1.0 };

    constexpr std::size_t p = 4;
    constexpr double tol = 1e-5;

    double R2n;
    double R2n_prev;
    double h;
    for (std::size_t i = 5; i <= 13; ++i)
    {
        h = 1.0 / std::pow(2.0, i);
        R2n = calculate_opp(a, b, h, x0, y0);
        
        
        if (i == 5)
        {
            R2n_prev = R2n;
            continue;
        }
        graphic_data_opp.emplace_back(std::log2(h), std::log2(R2n));
        R2n_prev = R2n;
    }
    //const double h_opt = h * std::pow(tol / std::abs(R2n), 1.0 / p) / 2.0;

    //std::cout << "\nOpponent h_opt = " << h_opt << "\n\n";
    
    //const auto [X1, Y1] = calculate_last_opp(a, b, h_opt, x0, y0);
    //const auto [X, Y] = calc_true(a, b, h_opt, x0, y0);
    
    //tol_graphic_data_opp = process_data(X, Y, Y1);

    std::cout << "\nCalculations are finished.\n\n";
}


TEST_METHOD(runge_kutta_calculation_auto)
{
    std::cout << "Runge-Kutta method with p=2 and h_auto:\n\n";

    constexpr double a = 0.0;
    constexpr double b = 5.0;
    constexpr double c2 = 0.65;
    const matrix<double> x0{ 0.0, 0.0, 0.0, 0.0 };
    const matrix<double> y0{ 1.0, 1.0, A,   1.0 };

    constexpr std::size_t p = 2;
    constexpr double rtol = 1e-6;
    constexpr double atol = 1e-12;

    std::tie(graphic_data_x, graphic_data_y) = calculate_with_h_auto(c2, a, b, x0, y0, rtol, atol);

    const double h_mid = 1.0 * (b - a) / graphic_data_x.get_rows_number();
    std::cout << "\nCalculate true functions with h_mid = " << h_mid << '\n';

    const auto [X, Y] = calc_true(a, b, h_mid, x0, y0);
    data_x_true = X;
    data_y_true = Y;
    tol_graphic_data_auto = process_data(X, Y, graphic_data_y);

    std::cout << "\nCalculations are finished.\n\n";
}


TEST_METHOD(runge_kutta_calculation_auto_opp)
{
    std::cout << "Runge-Kutta method with p=4 and h_auto:\n\n";

    constexpr double a = 0.0;
    constexpr double b = 5.0;
    const matrix<double> x0{ 0.0, 0.0, 0.0, 0.0 };
    const matrix<double> y0{ 1.0, 1.0, A,   1.0 };

    constexpr std::size_t p = 4;
    constexpr double rtol = 1e-6;
    constexpr double atol = 1e-12;

    std::tie(graphic_data_x_opp, graphic_data_y_opp) = calculate_with_h_auto_opp(a, b, x0, y0, rtol, atol);

    const double h_mid = 1.0 * (b - a) / graphic_data_x_opp.get_rows_number();
    std::cout << "\nCalculate true functions with h_mid = " << h_mid << '\n';

    const auto [X, Y] = calc_true(a, b, h_mid, x0, y0);
    data_x_true_opp = X;
    data_y_true_opp = Y;
    tol_graphic_data_auto_opp = process_data(X, Y, graphic_data_y_opp);

    std::cout << "\nCalculations are finished.\n\n";
}


TEST_METHOD(runge_kutta_calculation_auto_rtol)
{
    std::cout << "Runge-Kutta method with p=2 and h_auto, count addressing:\n\n";

    constexpr double a = 0.0;
    constexpr double b = 5.0;
    constexpr double c2 = 0.65;
    const matrix<double> x0{ 0.0, 0.0, 0.0, 0.0 };
    const matrix<double> y0{ 1.0, 1.0, A,   1.0 };

    constexpr std::size_t p = 2;
    constexpr double atol = 1e-12;

    double rtol;
    for (int i = 4; i <= 8; ++i)
    {
        rtol = std::pow(10.0, -i);
        calculate_with_h_auto(c2, a, b, x0, y0, rtol, atol);
    }
    
    std::cout << "\nCalculations are finished.\n\n";
}


TEST_METHOD(runge_kutta_calculation_auto_rtol_opp)
{
    std::cout << "Runge-Kutta method with p=4 and h_auto, count addressing:\n\n";

    constexpr double a = 0.0;
    constexpr double b = 5.0;
    const matrix<double> x0{ 0.0, 0.0, 0.0, 0.0 };
    const matrix<double> y0{ 1.0, 1.0, A,   1.0 };

    constexpr std::size_t p = 4;
    constexpr double atol = 1e-12;

    double rtol;
    for (int i = 4; i <= 8; ++i)
    {
        rtol = std::pow(10.0, -i);
        calculate_with_h_auto_opp(a, b, x0, y0, rtol, atol);
    }

    std::cout << "\nCalculations are finished.\n\n";
}

#endif // ENABLE_TESTS_RUNGE_KUTTA


} // namespace vv
