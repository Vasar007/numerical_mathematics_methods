// Copyright (C) 2018 Vasily Vasilyev (vasar007@yandex.ru)

#pragma once

#include <algorithm>
#include <array>
#include <stdexcept>
#include <type_traits>

#include "cx_matrix.hpp"
#include "constants.hpp"


namespace vv
{

/// ===== FUNCTION SECTION =====
template <typename Type, std::size_t Size>
constexpr std::array<cx_matrix<Type, Size, Size>, 2>
    lu_decompose(const cx_matrix<Type, Size, Size>& mat)
{
    using size_type = typename cx_matrix<Type, Size, Size>::size_type;

    cx_matrix<Type, Size, Size> lower{};
    cx_matrix<Type, Size, Size> upper{};

    // Decomposing cx_matrix into Upper and Lower triangular cx_matrix.
    for (size_type i = 0; i < Size; ++i)
    {
        // Upper Triangular.
        for (size_type k = i; k < Size; ++k)
        {
            // Summation of L(i, j) * U(j, k).
            Type sum{};
            for (size_type j = 0; j < i; ++j)
            {
                sum += lower.at(i, j) * upper.at(j, k);
            }
 
            // Evaluating U(i, k).
            upper.at(i, k) = mat.at(i, k) - sum;
        }
 
        // Lower Triangular.
        for (size_type k = i; k < Size; ++k)
        {
            if (i == k)
            {
                // Diagonal as 1.
                lower.at(i, i) = static_cast<Type>(1);
            }
            else
            {
                // Summation of L(k, j) * U(j, i).
                Type sum{};
                for (size_type j = 0; j < i; ++j)
                {
                    sum += lower.at(k, j) * upper.at(j, i);
                }
                // Evaluating L(k, i).
                lower.at(k, i) = (mat.at(k, i) - sum) / upper.at(i, i);
            }
        }
    } // for (size_type i = 0; i < Size; ++i)

    return { lower, upper };
}


template <typename Type, std::size_t Size>
constexpr long rank(const cx_matrix<Type, Size, Size>& mat)
{
    using size_type = typename cx_matrix<Type, Size, Size>::size_type;

    long rank = 0;
    const auto [L, U] = lu_decompose(mat);
    for (size_type i = 0; i < Size; ++i)
    {
        if (U.at(i, Size - 1) != Type{})
        {
            ++rank;
        }
    }
    return rank;
}


template <typename Type, std::size_t Size>
constexpr Type lu_determenant(const cx_matrix<Type, Size, Size>& mat)
{
    using size_type = typename cx_matrix<Type, Size, Size>::size_type;

    const auto [L, U] = lu_decompose(mat);

    Type detL = L.at(0, 0);
    Type detU = U.at(0, 0);
    for (size_type i = 1; i < Size; ++i)
    {
        detL *= L.at(i, i);
        detU *= U.at(i, i);
    }

    return detL * detU;
}


/// ===== TEST SECTION =====
#if ENABLE_TESTS_LU

TEST_METHOD(lu_test_decompose)
{
    std::cout << "LU-decompose:\n\n";
    
    constexpr auto lu_value = lu_decompose(data_eps);
    std::cout << "Matrix A:\n" << data_eps << "\n\n";
    std::cout << "Matrix L:\n" << lu_value.at(0) << "\n\n";
    std::cout << "Matrix U:\n" << lu_value.at(1);
    std::cout << "\n\n------------------------------------\n\n";

    constexpr auto lu_value_check = (lu_value.at(0) * lu_value.at(1)) - data_eps;
    std::cout << "Checking LU-decomposition:\n\n";
    std::cout << "Matrix (L * U) - A:\n" << lu_value_check;
    std::cout << "\n\n------------------------------------\n\n";
}


TEST_METHOD(lu_test_determenant)
{
    std::cout << "Determenant by LU:\n\n";
    
    constexpr auto lu_det = lu_determenant(data_eps);
    std::cout << "Determenant of A:\n" << lu_det;
    std::cout << "\n\n------------------------------------\n\n";
}


TEST_METHOD(lu_test_rank)
{
    std::cout << "LU-decomposition:\n\n";

    constexpr auto lu_value = vv::lu_decompose(data_eps);
    constexpr auto rank = vv::rank(data_eps);
    std::cout << "Matrix A:\n" << data_eps << "\n\n";
    std::cout << "Matrix L:\n" << lu_value.at(0) << "\n\n";
    std::cout << "Matrix U:\n" << lu_value.at(1) << "\n\n";
    std::cout << "Rank A:\n" << rank;
    std::cout << "\n\n------------------------------------\n\n";
}

#endif // ENABLE_TESTS_LU

} // namespace vv
