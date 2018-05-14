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
template <class Type, std::size_t Rows_A, std::size_t Columns_A>
constexpr std::tuple<cx_matrix<Type, Rows_A, Columns_A>, cx_matrix<Type, Rows_A, 1>, long>
    lup_decompose(cx_matrix<Type, Rows_A, Columns_A> A, const Type eps = kDefault_eps<Type>)
{
    // Alias the size_type.
    using size_type = typename cx_matrix<Type, Rows_A, Columns_A>::size_type;

    long permutations_counter = static_cast<long>(Rows_A);
    cx_matrix<Type, Rows_A, 1> P{};
    // Unit permutation cx_matrix, P[N] initialized with N.
    for (size_type i = 0; i < Rows_A; ++i)
    {
        P(i, 0) = static_cast<Type>(i);
    }
    
    for (size_type i = 0; i < Columns_A; ++i)
    {
        Type max_A{};
        size_type i_max = i;
        for (size_type k = i; k < Rows_A; ++k)
        {
            const auto absA = cx::abs(A(k, i));
            if (absA > max_A)
            { 
                max_A = absA;
                i_max = k;
            }
        }

        if (max_A < cx::abs(eps))
        { 
            // Failure, cx_matrix is degenerate.
            //throw std::domain_error("Matrix is singular!");
            return { cx_matrix<Type, Rows_A, Columns_A>{}, cx_matrix<Type, Rows_A, 1>{}, -1 };
        }

        if (i_max != i)
        {
            // Pivoting P.
            P.swap_rows(i, i_max);
            // Pivoting rows of A.
            A.swap_rows(i, i_max);
            // Counting pivots starting from N (for determinant).
            ++permutations_counter;
        }

        for (size_type j = i + 1; j < Rows_A; ++j)
        {
            A(j, i) /= A(i, i);
            for (size_type k = i + 1; k < Columns_A; ++k)
            {
                A(j, k) -= A(j, i) * A(i, k);
            }
        }
    } // for (size_type i = 0u; i < N; ++i)
    
    return { A, P, permutations_counter };
}


template <class Type, std::size_t Rows, std::size_t Columns_A, std::size_t Columns_b>
constexpr auto lup_solve(const cx_matrix<Type, Rows, Columns_A>& mat,
                         const cx_matrix<Type, Rows, Columns_b>& b,
                         const Type eps = kDefault_eps<Type>)
{
    static_assert(Columns_b == 1, "Matrix contains more than one columns in right "
                                  "hand side vector!");
                                   
    using size_type = typename cx_matrix<Type, Rows, Columns_A>::size_type;

    const auto [A, P, permutations_counter] = lup_decompose(mat, eps);

    cx_matrix<Type, Columns_A, 1> x{};
    for (size_type i = 0; i < Rows && i < Columns_A; ++i)
    {
        x(i, 0) = b(static_cast<size_type>(P(i, 0)), 0);

        for (size_type k = 0; k < i; ++k)
        {
            x(i, 0) -= A(i, k) * x(k, 0);
        }
    }

    for (size_type i = Columns_A - 1; i < Rows && i < Columns_A; --i)
    {
        for (size_type k = i + 1; k < Rows && k < Columns_A; ++k)
        {
            x(i, 0) -= A(i, k) * x(k, 0);
        }

        x(i, 0) /= A(i, i);
        if (i == 0u) break;
    }
    
    return x;
}


template <class Type, std::size_t Rows, std::size_t Columns>
constexpr Type lup_determenant(const cx_matrix<Type, Rows, Columns>& mat,
                           const Type eps = kDefault_eps<Type>)
{
    static_assert(Rows == Columns, "Matrix is not quadrant!");

    using size_type = typename cx_matrix<Type, Rows, Columns>::size_type;
    constexpr std::size_t N = Rows;

    const auto [A, P, permutations_counter] = lup_decompose(mat, eps); 
    Type det = A(0, 0);
    for (size_type i = 1; i < N; ++i)
    {
        det *= A(i, i);
    }

    return ((permutations_counter - N) % 2 == 0) ? det : -det;
}


template <class Type, std::size_t Rows, std::size_t Columns>
constexpr auto lup_invert(const cx_matrix<Type, Rows, Columns>& mat,
                          const Type eps = kDefault_eps<Type>)
{
    static_assert(Rows == Columns, "Matrix is not quadrant!");

    constexpr std::size_t N = Rows;
    using size_type = typename cx_matrix<Type, Rows, Columns>::size_type;

    const auto [A, P, permutations_counter] = lup_decompose(mat, eps);

    cx_matrix<Type, Rows, Columns> IA{};
    for (size_type j = 0; j < N; ++j)
    {
        for (size_type i = 0; i < N; ++i)
        {
            if (P(i, 0) == j)
            { 
                IA(i, j) = static_cast<Type>(1);
            }
            else
            {
                IA(i, j) = static_cast<Type>(0);
            }

            for (size_type k = 0; k < i; ++k)
            {
                IA(i, j) -= A(i, k) * IA(k, j);
            }
        }

        for (size_type i = N - 1; ; --i)
        {
            for (size_type k = i + 1; k < N; ++k)
            {
                IA(i, j) -= A(i, k) * IA(k, j);
            }

            IA(i, j) = IA(i, j) / A(i, i);

            if (i == 0) break;
        }
    } //  for (size_type j = 0; j < N; ++j)

    return IA;
}


/// ===== TEST SECTION =====
#if ENABLE_TESTS_LUP

TEST_METHOD(lup_test_decompose)
{
    std::cout << "LUP-decompose:\n\n";
    
    constexpr auto lup_tuple = lup_decompose(mat_A);
    constexpr auto mat_C = std::get<0>(lup_tuple);
    constexpr auto mat_P = std::get<1>(lup_tuple);
    constexpr auto pivots_number = std::get<2>(lup_tuple);
    std::cout << "Changed cx_matrix A:\n" << mat_C << "\n\n";
    std::cout << "Matrix P:\n" << mat_P << "\n\n";
    std::cout << "Number of pivots:\n" << pivots_number;
    std::cout << "\n\n------------------------------------\n\n";
}


TEST_METHOD(lup_test_solve)
{
    std::cout << "LUP-solve:\n\n";

    constexpr auto vec_x = lup_solve(mat_A, vec_b);
    std::cout << "Matrix A:\n" << mat_A << "\n\n";
    std::cout << "Matrix b:\n" << vec_b << "\n\n";
    std::cout << "Matrix X:\n" << vec_x;
    std::cout << "\n\n------------------------------------\n\n";

    std::cout << "Checking LUP-solve:\n\n";
    constexpr auto lup_solve_check = (mat_A * vec_x) - vec_b;
    std::cout << "Result of (A * x - b):\n" << lup_solve_check;
    std::cout << "\n\n------------------------------------\n\n";
}


#endif // ENABLE_TESTS_LUP

#if ENABLE_TESTS_LUP_ONLY_SQUARE

TEST_METHOD(lup_test_determenant)
{
    std::cout << "Determenant by LUP:\n\n";    

    constexpr auto lup_det = lup_determenant(mat_A);
    std::cout << "Determenant of A:\n" << lup_det;
    std::cout << "\n\n------------------------------------\n\n";
}


TEST_METHOD(lup_test_inverse)
{
    std::cout << "Invert by LUP:\n\n";

    constexpr auto lup_invert_mat = lup_invert(mat_A);
    std::cout << "Matrix A:\n" << mat_A << "\n\n";
    std::cout << "Inverted cx_matrix to A:\n" << lup_invert_mat;
    std::cout << "\n\n------------------------------------\n\n";

    std::cout << "Checking inverse:\n\n";
    constexpr auto inverse_check1 = mat_A * lup_invert_mat;
    std::cout << "Result of (A * IA):\n" << inverse_check1 << "\n\n";

    constexpr auto inverse_check2 = lup_invert_mat * mat_A;    
    std::cout << "Result of (IA * A):\n" << inverse_check2;
    std::cout << "\n\n------------------------------------\n\n";
}

#endif // ENABLE_TESTS_LUP_ONLY_SQUARE

} // namespace vv
