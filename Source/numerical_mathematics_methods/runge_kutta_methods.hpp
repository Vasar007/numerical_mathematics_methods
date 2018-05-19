#pragma once

#include <iostream>
#include <cmath>
#include <vector>
#include <tuple>
#include <cassert>
#include <fstream>
#include <string>
#include <string_view>

#include "utils.hpp"
#include "matrix.hpp"



namespace utils
{

void out_data(const std::string& file_name, const std::string_view mode,
              const std::string_view param, const std::string_view title,
              const std::string_view x_label, const std::string_view y_label,
              const std::vector<std::pair<double, double>>& data)
{
    if (data.empty())
    {
        std::cout << "ERROR: empty data to process for file " << file_name << ".\n";
        return;
    }

    std::ofstream out_file(file_name);
    out_file << mode << '|' << param << '|' << title << '|' << x_label << '|' << y_label << '\n';
    for (const auto& [x, y] : data)
    {
        out_file << x << ' ' << y << '\n';
    }
}


void out_data(const std::string& file_name, const std::string_view mode,
              const std::string_view param, const std::string_view title,
              const std::string_view x_label, const std::string_view y_label,
              const std::vector<std::pair<double, double>>& data_1,
              const std::vector<std::pair<double, double>>& data_2)
{
    if (data_1.empty() || data_2.empty())
    {
        std::cout << "ERROR: empty data to process for file " << file_name << ".\n";
        return;
    }

    std::ofstream out_file(file_name);
    out_file << mode << '|' << param << '|' << title << '|' << x_label << '|' << y_label << '\n';
    for (const auto& [x, y] : data_1)
    {
        out_file << x << ' ' << y << '\n';
    }
    out_file << "#\n";
    for (const auto [x, y] : data_2)
    {
        out_file << x << ' ' << y << '\n';
    }
}


void out_data(const std::string& file_name, const std::string_view mode,
              const std::string_view param, const std::string_view title,
              const std::string_view x_label, const std::string_view y_label,
              const vv::matrix<double>& x, const vv::matrix<double>& y)
{
    if (x.empty() || y.empty())
    {
        std::cout << "ERROR: empty data to process for file " << file_name << ".\n";
        return;
    }
    if (x.get_rows_number() != y.get_rows_number())
    {
        std::cout << "\nERROR: X and Y have different size! Cannot process data.\n";
        return;
    }

    std::ofstream out_file(file_name);
    out_file << mode << '|' << param << '|' << title << '|' << x_label << '|' << y_label << '\n';
    for (std::size_t i = 0; i < x.get_rows_number(); ++i)
    {
        out_file << x.at(i, 0) << ' ' << y.at(i, 0) << ' ' << y.at(i, 1) << ' ' << y.at(i, 2)
            << ' ' << y.at(i, 3) << '\n';
    }
}


void out_data(const std::string& file_name, const std::string_view mode,
              const std::string_view param, const std::string_view title,
              const std::string_view x_label, const std::string_view y_label,
              const vv::matrix<double>& x1, const vv::matrix<double>& y1,
              const vv::matrix<double>& x2, const vv::matrix<double>& y2)
{
    if (x1.empty() || y1.empty() || x2.empty() || y2.empty())
    {
        std::cout << "ERROR: empty data to process for file " << file_name << ".\n";
        return;
    }
    if (x1.get_rows_number() != y1.get_rows_number()
     || x2.get_rows_number() != y2.get_rows_number())
    {
        std::cout << "\nERROR: X and Y have different size! Cannot process data.\n";
        return;
    }

    std::ofstream out_file(file_name);
    out_file << mode << '|' << param << '|' << title << '|' << x_label << '|' << y_label << '\n';
    for (std::size_t i = 0; i < x1.get_rows_number(); ++i)
    {
        out_file << x1.at(i, 0) << ' ' << y1.at(i, 0) << ' ' << y1.at(i, 1) << ' ' << y1.at(i, 2)
            << ' ' << y1.at(i, 3) << '\n';
    }
    out_file << "#\n";
    for (std::size_t i = 0; i < x2.get_rows_number(); ++i)
    {
        out_file << x2.at(i, 0) << ' ' << y2.at(i, 0) << ' ' << y2.at(i, 1) << ' ' << y2.at(i, 2)
            << ' ' << y2.at(i, 3) << '\n';
    }
}

} // namespace utils

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

    matrix<double> Y(n + 1, kNumber_F);
    matrix<double> X(n + 1, 1);

    X.at(0, 0) = x0.at(0, 0);
    Y.at(0) = y0.at(0);
    for (std::size_t i = 1; i <= n; ++i)
    {
        X.at(i, 0) = a + i * h;
        Y.at(i) = f(X.at(i, 0));
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

matrix::container<double> runga_kutta_2(const double c2, const double a21, const double b1,
                                        const double b2, const double h, double& X,
                                        matrix::container<double>& K1,
                                        matrix::container<double>& K2,
                                        matrix::container<double>& Y1,
                                        const matrix::container<double>& Y)
{
    const matrix::container<double> F_xy = F(X, Y);
    K1 = h * F_xy;
    K2 = h * F(X + c2 * h, Y + a21 * K1);
    Y1 = Y + (b1 * K1 + b2 * K2);
    X += h;
    return F_xy;
}

void runga_kutta_2(const double c2, const double a21, const double b1, const double b2,
                   const double h, double& X, matrix::container<double>& K1,
                   matrix::container<double>& K2, matrix::container<double>& Y1,
                   const matrix::container<double>& Y, const matrix::container<double>& F_xy)
{
    K1 = h * F_xy;
    K2 = h * F(X + c2 * h, Y + a21 * K1);
    Y1 = Y + (b1 * K1 + b2 * K2);
    X += h;
}


void runga_kutta_2(const double c2, const double a, const double a21, const double b1,
                   const double b2, const double h, const std::size_t i, ::vv::matrix<double>& X,
                   ::vv::matrix<double>& K1, ::vv::matrix<double>& K2, ::vv::matrix<double>& Y)
{
    X.at(i, 0) = a + i * h;
    K1.at(i) = h * F(X.at(i - 1, 0), Y.at(i - 1));
    K2.at(i) = h * F(X.at(i - 1, 0) + c2 * h, Y.at(i - 1) + a21 * K1.at(i));
    Y.at(i) = Y.at(i - 1) + (b1 * K1.at(i) + b2 * K2.at(i));
}


matrix::container<double> runga_kutta_4(const double h, double& X, matrix::container<double>& K1,
                                        matrix::container<double>& K2,
                                        matrix::container<double>& K3,
                                        matrix::container<double>& K4,
                                        matrix::container<double>& Y1,
                                        const matrix::container<double>& Y)
{
    const matrix::container<double> F_xy = F(X, Y);
    K1 = h * F_xy;
    K2 = h * F(X + h / 3.0, Y + K1 / 3.0);
    K3 = h * F(X + 2.0 * h / 3.0, Y - K1 / 3.0 + K2);
    K4 = h * F(X + h, Y + K1 - K2 + K3);
    Y1 = Y + (K1 + 3.0 * K2 + 3.0 * K3 + K4) / 8.0;
    X += h;
    return F_xy;
}


void runga_kutta_4(const double h, double& X, matrix::container<double>& K1,
                   matrix::container<double>& K2, matrix::container<double>& K3,
                   matrix::container<double>& K4, matrix::container<double>& Y1,
                   const matrix::container<double>& Y, const matrix::container<double>& F_xy)
{
    K1 = h * F_xy;
    K2 = h * F(X + h / 3.0, Y + K1 / 3.0);
    K3 = h * F(X + 2.0 * h / 3.0, Y - K1 / 3.0 + K2);
    K4 = h * F(X + h, Y + K1 - K2 + K3);
    Y1 = Y + (K1 + 3.0 * K2 + 3.0 * K3 + K4) / 8.0;
    X += h;
}


void runga_kutta_4(const double a, const double h, const std::size_t i, ::vv::matrix<double>& X,
                   ::vv::matrix<double>& K1, ::vv::matrix<double>& K2, ::vv::matrix<double>& K3,
                   ::vv::matrix<double>& K4, ::vv::matrix<double>& Y)
{
    X.at(i, 0) = a + i * h;
    K1.at(i) = h * F(X.at(i - 1, 0), Y.at(i - 1));
    K2.at(i) = h * F(X.at(i - 1, 0) + h / 3.0, Y.at(i - 1) + K1.at(i) / 3.0);
    K3.at(i) = h * F(X.at(i - 1, 0) + 2.0 * h / 3.0, Y.at(i - 1) - K1.at(i) / 3.0 + K2.at(i));
    K4.at(i) = h * F(X.at(i - 1, 0) + h, Y.at(i - 1) + K1.at(i) - K2.at(i) + K3.at(i));
    Y.at(i) = Y.at(i - 1) + (K1.at(i) + 3.0 * K2.at(i) + 3.0 * K3.at(i) + K4.at(i)) / 8.0;
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

    X.at(0, 0) = x0.at(0, 0);
    Y.at(0) = y0.at(0);
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

    X.at(0, 0) = x0.at(0, 0);
    Y.at(0) = y0.at(0);
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

    return std::sqrt(result);
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
        result.emplace_back(x.at(i, 0), norm(y.at(i) - y1.at(i)));
    }

    return result;
}


matrix<double> data_x_true;
matrix<double> data_y_true(5, kNumber_F, mtx::mode::dynamically_expandable);

matrix<double> data_x_true_opp;
matrix<double> data_y_true_opp(5, kNumber_F, mtx::mode::dynamically_expandable);

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
    assert(b > a);
    const std::size_t n = static_cast<std::size_t>((b - a) / h);
    const std::size_t n2 = n * 2;

    const auto [X1, Y1] = runge_kutta_2(c2, a, b, h, x0, y0);
    const auto [X1_2, Y1_2] = runge_kutta_2(c2, a, b, h / 2.0, x0, y0);
    
    const detail::matrix::container<double> yn = Y1.at(n);
    const detail::matrix::container<double> y2n = Y1_2.at(n2);

    constexpr std::size_t p = 2;
    const double R2n = norm(y2n - yn) / (std::pow(2, p) - 1.0);

    std::cout << "Step = " << h << "\nR2n = " << R2n << '\n';

    return R2n;
}


std::tuple<matrix<double>, matrix<double>>
calculate_last(const double c2, const double a, const double b, const double h,
               const matrix<double>& x0, const matrix<double>& y0)
{
    assert(b > a);
    const std::size_t n = static_cast<std::size_t>((b - a) / h);
    const std::size_t n2 = n * 2;

    const auto[X1, Y1] = runge_kutta_2(c2, a, b, h, x0, y0);
    const auto[X1_2, Y1_2] = runge_kutta_2(c2, a, b, h / 2.0, x0, y0);

    const detail::matrix::container<double> yn = Y1.at(n);
    const detail::matrix::container<double> y2n = Y1_2.at(n2);

    constexpr std::size_t p = 2;
    const double Rn = norm(y2n - yn) / (1.0 - 1.0 / std::pow(2.0, p));
    std::cout << "Step = " << h << "\nRn = " << Rn << '\n';

    return { X1, Y1 };
}


double calculate_opp(const double a, const double b, const double h,
                     const matrix<double>& x0, const matrix<double>& y0)
{
    assert(b > a);
    const std::size_t n = static_cast<std::size_t>((b - a) / h);
    const std::size_t n2 = n * 2;

    const auto [X2, Y2] = runge_kutta_4(a, b, h, x0, y0);
    const auto [X2_2, Y2_2] = runge_kutta_4(a, b, h / 2.0, x0, y0);

    const detail::matrix::container<double> yn_opp = Y2.at(n);
    const detail::matrix::container<double> y2n_opp = Y2_2.at(n2);

    constexpr std::size_t p = 4;
    const double R2n = norm(y2n_opp - yn_opp) / (std::pow(2, p) - 1.0);

    std::cout << "Opponent step = " << h << "\nR2n = " << R2n << '\n';

    return R2n;
}

std::tuple<matrix<double>, matrix<double>>
calculate_last_opp(const double a, const double b, const double h, const matrix<double>& x0,
                   const matrix<double>& y0)
{
    assert(b > a);
    const std::size_t n = static_cast<std::size_t>((b - a) / h);
    const std::size_t n2 = n * 2;

    const auto [X1, Y1] = runge_kutta_4(a, b, h, x0, y0);
    const auto [X1_2, Y1_2] = runge_kutta_4(a, b, h / 2.0, x0, y0);

    const detail::matrix::container<double> yn_opp = Y1.at(n);
    const detail::matrix::container<double> y2n_opp = Y1_2.at(n2);

    constexpr std::size_t p = 4;
    const double Rn = norm(y2n_opp - yn_opp) / (1.0 - 1.0 / std::pow(2.0, p));
    std::cout << "Opponent step = " << h << "\nRn = " << Rn << '\n';

    return { X1, Y1 };
}


double calculate_h0(const double rtol, const double atol, const std::size_t p,
                    const matrix<double>& x0, const matrix<double>& y0)
{
    const auto abs_compare = [](const double a1, const double a2)
                             { return std::abs(a1) < std::abs(a2); };
    const double tol = norm(y0.at(0)) * rtol + atol;
    // Step 1:
    detail::matrix::container<double> f0 = F(x0.at(0, 0), y0.at(0));

    // Step 2:
    double max_x = *std::max_element(std::begin(x0.at(0)), std::end(x0.at(0)), abs_compare);
    double delta = std::pow(1.0 / (max_x + 1.0), p + 1) + std::pow(norm(f0), p + 1);

    // Step 3:
    double h1 = std::pow(tol / delta, 1.0 / (p + 1.0));

    // Step 4:
    double ux = x0.at(0, 0) + h1;
    detail::matrix::container<double> uy = y0.at(0) + h1 * F(x0.at(0, 0), y0.at(0));

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
    assert(b > a);
    assert(c2 != 0.0);

    const double a21 = c2;
    const double b2 = 1.0 / (2.0 * c2);
    const double b1 = 1.0 - b2;

    constexpr std::size_t p = 2;

    matrix<double> Y(1, kNumber_F, mtx::mode::dynamically_expandable);
    matrix<double> X(1, 1, mtx::mode::dynamically_expandable);

    detail::matrix::container<double> K1_1;
    detail::matrix::container<double> K2_1;
    detail::matrix::container<double> Y_1 = y0.at(0);
    double X_1 = x0.at(0, 0);

    detail::matrix::container<double> Y_12;
    double X_12;

    detail::matrix::container<double> K1_2;
    detail::matrix::container<double> K2_2;
    detail::matrix::container<double> Y_2 = y0.at(0);
    double X_2 = x0.at(0, 0);

    X.at(0, 0) = x0.at(0, 0);
    Y.at(0) = y0.at(0);
    data_y_true.at(0) = y0.at(0);

    std::vector<double> h_vec{ calculate_h0(rtol, atol, p, x0, y0) };

    std::size_t i = 1;
    while (X.at(i - 1, 0) + h_vec.at(i - 1) <= b)
    {
        std::cout << "Step = " << (i - 1) << " \th = " << h_vec.at(i - 1) << '\t';
        const detail::matrix::container<double> F_xy = detail::runga_kutta_2(c2, a21, b1, b2,
                                                                h_vec.at(i - 1), X_1, K1_1, K2_1,
                                                                Y_1, Y.at(i - 1));

        detail::runga_kutta_2(c2, a21, b1, b2, h_vec.at(i - 1) / 2.0, X_2, K1_2, K2_2, Y_2,
                              Y.at(i - 1), F_xy);
        X_12 = X_2;
        Y_12 = Y_2;
        detail::runga_kutta_2(c2, a21, b1, b2, h_vec.at(i - 1) / 2.0, X_2, K1_2, K2_2, Y_2, Y_2);

        const double tol = rtol * std::max({ norm(Y.at(i - 1)), norm(Y_1), norm(Y_2) }) + atol;
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
            Y.at(i) = Y_2;
            X.at(i, 0) = X_1;
            data_y_true.at(i) = f(X_1);

            quality.emplace_back(X_1, norm(data_y_true.at(i) - Y_2) / rn);
            ++i;

            std::cout << "Chosen 2 way.\n";
        }
        else if (rn >= tol / std::pow(2.0, p + 1))
        {
            h_acc.emplace_back(X_1, h_vec.at(i - 1));

            h_vec.emplace_back(h_vec.at(i - 1));
            Y.at(i) = Y_1;
            X.at(i, 0) = X_1;
            data_y_true.at(i) = f(X_1);

            quality.emplace_back(X_1, norm(data_y_true.at(i) - Y_1) / rn);
            ++i;

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
            Y.at(i) = Y_1;
            X.at(i, 0) = X_1;
            data_y_true.at(i) = f(X_1);

            quality.emplace_back(X_1, norm(data_y_true.at(i) - Y_1) / rn);
            ++i;

            std::cout << "Chosen 4 way.\n";
        }
    }
    constexpr double addressing = 3.0;
    constexpr double usage_F_rkm_2 = 2.0;
    constexpr double s = addressing * usage_F_rkm_2 - 1.0;
    count_F.emplace_back(-std::log10(rtol), std::log10((i - 1) * s));

    return { X, Y };
}


std::tuple<matrix<double>, matrix<double>>
calculate_with_h_auto_opp(const double a, const double b, const matrix<double>& x0,
                          const matrix<double>& y0, const double rtol, const double atol)
{
    assert(b > a);
    constexpr std::size_t p = 4;

    matrix<double> Y(1, kNumber_F, mtx::mode::dynamically_expandable);
    matrix<double> X(1, 1, mtx::mode::dynamically_expandable);

    detail::matrix::container<double> K1_1;
    detail::matrix::container<double> K2_1;
    detail::matrix::container<double> K3_1;
    detail::matrix::container<double> K4_1;
    detail::matrix::container<double> Y_1 = y0.at(0);
    double X_1 = x0.at(0, 0);

    detail::matrix::container<double> Y_12;
    double X_12;

    detail::matrix::container<double> K1_2;
    detail::matrix::container<double> K2_2;
    detail::matrix::container<double> K3_2;
    detail::matrix::container<double> K4_2;
    detail::matrix::container<double> Y_2 = y0.at(0);
    double X_2 = x0.at(0, 0);

    X.at(0, 0) = x0.at(0, 0);
    Y.at(0) = y0.at(0);
    data_y_true_opp.at(0) = y0.at(0);

    std::vector<double> h_vec{ calculate_h0(rtol, atol, p, x0, y0) };

    std::size_t i = 1;
    while (X.at(i - 1, 0) + h_vec.at(i - 1) <= b)
    {
        std::cout << "Step = " << (i - 1) << " \th = " << h_vec.at(i - 1) << '\t';
        const detail::matrix::container<double> F_xy = detail::runga_kutta_4(h_vec.at(i - 1), X_1,
                                                            K1_1, K2_1, K3_1, K4_1, Y_1,
                                                            Y.at(i - 1));

        detail::runga_kutta_4(h_vec.at(i - 1) / 2.0, X_2, K1_2, K2_2, K3_2, K4_2, Y_2, Y.at(i - 1),
                              F_xy);
        X_12 = X_2;
        Y_12 = Y_2;
        detail::runga_kutta_4(h_vec.at(i - 1) / 2.0, X_2, K1_2, K2_2, K3_2, K4_2, Y_2, Y_2);

        const double tol = rtol * std::max({ norm(Y.at(i - 1)), norm(Y_1), norm(Y_2) }) + atol;
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
            Y.at(i) = Y_2;
            X.at(i, 0) = X_1;
            data_y_true_opp.at(i) = f(X_1);

            quality_opp.emplace_back(X_1, norm(data_y_true_opp.at(i) - Y_2) / rn);
            ++i;

            std::cout << "Chosen 2 way.\n";
        }
        else if (rn >= tol / std::pow(2.0, p + 1))
        {
            h_acc_opp.emplace_back(X_1, h_vec.at(i - 1));

            h_vec.emplace_back(h_vec.at(i - 1));
            Y.at(i) = Y_1;
            X.at(i, 0) = X_1;
            data_y_true_opp.at(i) = f(X_1);

            quality_opp.emplace_back(X_1, norm(data_y_true_opp.at(i) - Y_1) / rn);
            ++i;

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
            Y.at(i) = Y_1;
            X.at(i, 0) = X_1;
            data_y_true_opp.at(i) = f(X_1);

            quality_opp.emplace_back(X_1, norm(data_y_true_opp.at(i) - Y_1) / rn);
            ++i;

            std::cout << "Chosen 4 way.\n";
        }
    }
    constexpr double addressing = 3.0;
    constexpr double usage_F_rkm_4 = 4.0;
    constexpr double s = addressing * usage_F_rkm_4 - 1.0;
    count_F_opp.emplace_back(-std::log10(rtol), std::log10((i - 1) * s));

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
    double h;
    for (std::size_t i = 5; i <= 13; ++i) // If i < 5 may appear NAN because of log.
    {
        h = 1.0 / std::pow(2.0, i);
        R2n = calculate(c2, a, b, h, x0, y0);

        graphic_data.emplace_back(std::log2(h), std::log2(R2n));
    }
    constexpr double guarantee_multiplier = 0.95;
    const double h_opt = guarantee_multiplier * h * std::pow(tol / std::abs(R2n), 1.0 / p) / 2.0;

    std::cout << "\nh_opt = " << h_opt << "\n\n";
    
    const auto [X1, Y1] = calculate_last(c2, a, b, h_opt, x0, y0);
    const auto [X, Y] = calc_true(a, b, h_opt, x0, y0);
    
    tol_graphic_data = process_data(X, Y, Y1);

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
    double h;
    for (std::size_t i = 5; i <= 13; ++i) // If i < 5 may appear NAN because of log.
    {
        h = 1.0 / std::pow(2.0, i);
        R2n = calculate_opp(a, b, h, x0, y0);

        graphic_data_opp.emplace_back(std::log2(h), std::log2(R2n));
    }
    constexpr double guarantee_multiplier = 0.95;
    const double h_opt = guarantee_multiplier * h * std::pow(tol / std::abs(R2n), 1.0 / p) / 2.0;
    
    std::cout << "\nOpponent h_opt = " << h_opt << "\n\n";
    
    const auto [X1, Y1] = calculate_last_opp(a, b, h_opt, x0, y0);
    const auto [X, Y] = calc_true(a, b, h_opt, x0, y0);
    
    tol_graphic_data_opp = process_data(X, Y, Y1);

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

    data_x_true = graphic_data_x;
    tol_graphic_data_auto = process_data(data_x_true, data_y_true, graphic_data_y);

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

    data_x_true_opp = graphic_data_x_opp;
    tol_graphic_data_auto_opp = process_data(data_x_true_opp, data_y_true_opp, graphic_data_y_opp);

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


TEST_METHOD(do_runge_kutta_tests)
{
    // Calculate RKM with p = 2 and then calculate with h_opt.
    runge_kutta_calculation();
    
    // Calculate RKM with p = 4 and then calculate with h_opt.
    runge_kutta_calculation_opp();

    // Calculate RKM_auto with p = 2 and also save plot and quality with h_acc and h_not_acc.
    runge_kutta_calculation_auto();
    
    // Calculate RKM_auto with p = 4 and also save plot and quality with h_acc and h_not_acc.
    runge_kutta_calculation_auto_opp();


    utils::out_data("rkm_h.txt", "2p", "def",
                    "RKM with h = 1/2^k", "log2(h)", "log2(R2n)",
                    graphic_data, graphic_data_opp);
    utils::out_data("rkm_h_opt.txt", "2p", "def", "RKM with h_opt", "x", "y - y1",
                    tol_graphic_data, tol_graphic_data_opp);
    
    utils::out_data("rkm_2_auto.txt", "1p", "def",
                    "Tolerance for RKM_auto with p = 2", "x", "y - y1",
                    tol_graphic_data_auto);
    utils::out_data("rkm_2_auto_plot.txt", "8p", "r;r;r;r;b;b;b;b",
                    "Plot from RKM_auto with p = 2", "x", "y",
                    graphic_data_x, graphic_data_y,  data_x_true_opp, data_y_true_opp);
    utils::out_data("rkm_2_auto_quality.txt", "1p", "def",
                    "Quality for RKM_auto with p = 2", "x", "Quality",
                    quality);
    utils::out_data("rkm_2_auto_h.txt", "2p", "r.;b.",
                     "Choice 'h' for RKM_auto with p = 2", "x", "h",
                     h_acc, h_not_acc);
    
    utils::out_data("rkm_4_auto.txt", "1p", "def",
                    "Tolerance for RKM_auto with p = 4", "x", "y - y1",
                    tol_graphic_data_auto_opp);
    utils::out_data("rkm_4_auto_plot.txt", "8p", "r;r;r;r;b;b;b;b",
                    "Plot from RKM_auto with p = 4", "x", "y",
                    graphic_data_x_opp, graphic_data_y_opp, data_x_true_opp, data_y_true_opp);
    utils::out_data("rkm_4_auto_quality.txt", "1p", "def",
                    "Quality for RKM_auto with p = 4", "x", "Quality",
                    quality_opp);
    utils::out_data("rkm_4_auto_h.txt", "2p", "r.;b.",
                    "Choice 'h' for RKM_auto with p = 4", "x", "h",
                    h_acc_opp, h_not_acc_opp);


    // Calculate RKM_auto and count addressing to F.
    count_F.clear();
    runge_kutta_calculation_auto_rtol();
    // And for opponent method.
    count_F_opp.clear();
    runge_kutta_calculation_auto_rtol_opp();


    utils::out_data("rkm_auto_addressing.txt", "2p", "def",
                    "Counting addresses to F", "-log10(rtol)", "log10(Addressing_F)",
                    count_F, count_F_opp);
}

#endif // ENABLE_TESTS_RUNGE_KUTTA


} // namespace vv
