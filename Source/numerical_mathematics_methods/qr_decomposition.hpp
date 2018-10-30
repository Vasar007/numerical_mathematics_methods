// Copyright (C) 2018 Vasily Vasilyev (vasar007@yandex.ru)

#pragma once

#include <algorithm>
#include <array>
#include <stdexcept>
#include <type_traits>

#include "cx_matrix.hpp"
#include "constants.hpp"
#include "cx_math.h"


namespace vv
{

/// ===== FUNCTION SECTION =====
template <typename Type, std::size_t Rows, std::size_t Columns>
constexpr std::array<cx_matrix<Type, Rows, Columns>, 2>
    qr_decompose(const cx_matrix<Type, Rows, Columns>& mat, const Type eps = kDefault_eps<Type>)
{
    using size_type = typename cx_matrix<Type, Rows, Columns>::size_type;

    // Create identity cx_matrix.
    auto I = cx_matrix<Type, Rows, Rows>::create_identity();
    cx_matrix<Type, Rows, Rows> P = I;
    cx_matrix<Type, Rows, Rows> Q = I;
    
    cx_matrix<Type, Rows, Columns> R = mat;

    // Using Householder reflections.
    for (size_type i = 0; i < Columns; ++i)
    {
        cx_matrix<Type, Rows, 1> u{};
        cx_matrix<Type, Rows, 1> v{};
        
        Type norm{};
        for (size_type j = i; j < Rows; ++j)
        {
            u.at(j, 0) = R.at(j, i); // [j * Columns + i] 
            norm += u.at(j, 0) * u.at(j, 0);
        }
        norm = cx::sqrt(norm);
        
        Type alpha = u.at(i, 0) < Type{} ? -norm : norm;

        norm = Type{};
        for (size_type j = i; j < Rows; ++j)
        {
            v.at(j, 0) = (j == i ? u.at(j, 0) + alpha : u.at(j, 0));
            norm += v.at(j, 0) * v.at(j, 0);
        }
        norm = cx::sqrt(norm);

        if (norm < eps) continue;

        for (size_type j = i; j < Rows; ++j)
        {
            v.at(j, 0) /= norm;
        }

        P = I - (v * v.transpose()) * static_cast<Type>(2);

        R = P * R;
        Q = Q * P;

    } // for (size_type i = 0; i < Columns; ++i)
    
    return { Q, R };
}


/// ===== TEST SECTION =====
#if ENABLE_TESTS_QR

TEST_METHOD(qr_test_decomposition)
{
    std::cout << "QR-decompose:\n\n";
    
    constexpr auto qr_value = qr_decompose(mat_A);
    std::cout << "Matrix A:\n" << mat_A << "\n\n";
    std::cout << "Matrix Q:\n" << qr_value.at(0) << "\n\n";
    std::cout << "Matrix R:\n" << qr_value.at(1);
    std::cout << "\n\n------------------------------------\n\n";

    std::cout << "Checking QR-decompose:\n\n";
    constexpr auto qr_product = qr_value.at(0) * qr_value.at(1);
    std::cout << "Matrix (Q * R):\n" << qr_product;
    std::cout << "\n\n------------------------------------\n\n";
}

#endif // ENABLE_TESTS_QR

#if ENABLE_TESTS_QR_SOLVE

TEST_METHOD(qr_test_solve)
{
    std::cout << "QR-solve:\n\n";
    
    constexpr auto qr_value = qr_decompose(mat_A);
    constexpr auto rx = qr_value.at(1);
    constexpr auto qTb = qr_value.at(0).transpose() * vec_b;
    constexpr auto vecX_qr = vv::cx_matrix<double>::backward_substitution(rx, qTb, kEps);
    std::cout << "Matrix solved with Gaussian (include eps):\n\n";
    std::cout << "Matrix (QT * b):\n" << qTb << "\n\n";
    std::cout << "Matrix x:\n" << vecX_qr;
    std::cout << "\n\n------------------------------------\n\n";

    std::cout << "Checking QR-solve:\n\n";
    constexpr auto qr_solve_check = (rx * vecX_qr) - qTb;
    std::cout << "Result of (R * x - QT * b):\n" << qr_solve_check;
    std::cout << "\n\n------------------------------------\n\n";
}

#endif // ENABLE_TESTS_QR_SOLVE

} // namespace vv
