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

namespace detail::lupq
{

/// ===== ADDITION FUNCTIONAL SECTION =====
template <class Type, std::size_t Rows, std::size_t Columns>
constexpr ::vv::cx_matrix<Type, Rows, Columns>
    lupq_get_u_solution_impl(const ::vv::cx_matrix<Type, Rows, Columns>& C,
                             const std::size_t rank) noexcept
{
    using size_type = typename ::vv::cx_matrix<Type, Rows, Columns>::size_type;

    ::vv::cx_matrix<Type, Rows, Columns> result_U{};

    // Singular matrix.
    if (rank < Rows || rank < Columns)
    {
        for (size_type i = 0; i < rank; ++i)
        {
            for (size_type j = i; j < rank; j++)
            {
                result_U.at(i, j) = C.at(i, j);
            }
            for (size_type j = rank; j < Columns; j++)
            {
                result_U.at(i, j) = Type{};
            }
        }
        for (size_type i = rank; i < Rows; i++)
        {
            for (size_type j = i; j < Columns; j++)
            {
                result_U.at(i, j) = Type{};
            }
        }
    }
    // Regular matrix.
    else
    {
        for (size_type i = 0; i < Rows; ++i)
        {
            for (size_type j = i; j < Columns; ++j)
            {
                result_U.at(i, j) = C.at(i, j);
            }
        }
    }

    return result_U;
}


template <class Type, std::size_t Rows, std::size_t Columns>
constexpr ::vv::cx_matrix<Type, Rows, Rows>
    lupq_get_l_solution_impl(const ::vv::cx_matrix<Type, Rows, Columns>& C,
                             const std::size_t rank) noexcept
{
    using size_type = typename ::vv::cx_matrix<Type, Rows, Columns>::size_type;
    
    ::vv::cx_matrix<Type, Rows, Rows> result_L{};

    // Singular matrix.
    if (rank < Rows)
    {
        for (size_type i = 0; i < rank; ++i)
        {
            result_L.at(i, i) = static_cast<Type>(1);

            if (i == 0) continue;

            for (size_type j = i - 1; ; --j)
            {
                result_L.at(i, j) = C.at(i, j);

                if (j == 0) break;
            }
        }
        for (size_type i = rank; i < Rows; ++i)
        {
            if (i == 0) continue;

            for (size_type j = i; ; --j)
            {
                result_L.at(i, j) = Type{};

                if (j == 0) break;
            }
        }
    }
    // Regular matrix.
    else
    {
        for (size_type i = 0; i < Rows; ++i)
        {
            result_L.at(i, i) = static_cast<Type>(1);

            if (i == 0) continue;

            for (size_type j = i - 1; ; --j)
            {
                result_L.at(i, j) = C.at(i, j);

                if (j == 0) break;
            }
        }
    }

    return result_L;
}


template <class Type, std::size_t Rows, std::size_t Columns>
constexpr ::vv::cx_matrix<Type, Rows, Columns>
    lupq_get_u_impl(const ::vv::cx_matrix<Type, Rows, Columns>& C) noexcept
{
    using size_type = typename ::vv::cx_matrix<Type, Rows, Columns>::size_type;

    ::vv::cx_matrix<Type, Rows, Columns> result_U{};
    for (size_type i = 0; i < Rows; ++i)
    {
        for (size_type j = i; j < Columns; ++j)
        {
            result_U.at(i, j) = C.at(i, j);
        }
    }

    return result_U;
}

template <class Type, std::size_t Rows, std::size_t Columns>
constexpr ::vv::cx_matrix<Type, Rows, Rows>
    lupq_get_l_impl(const ::vv::cx_matrix<Type, Rows, Columns>& C) noexcept
{
    using size_type = typename ::vv::cx_matrix<Type, Rows, Columns>::size_type;

    ::vv::cx_matrix<Type, Rows, Rows> result_L{};

    result_L.at(0, 0) = static_cast<Type>(1);
    for (size_type i = 1; i < Rows; ++i)
    {
        size_type temp_index = i >= Columns ? Columns - 1 : i;
        for (size_type j = temp_index - 1; ; --j)
        {
            if (Columns < Rows)
            {
                result_L.at(i, j) = Type{};
            }
            else
            {
                result_L.at(i, j) = C.at(i, j);
            }

            if (j == 0) break;
        }
    }
    
    return result_L;
}


template <class Type, std::size_t Rows>
constexpr ::vv::cx_matrix<Type, Rows, Rows>
    lupq_get_p_impl(const ::vv::cx_matrix<Type, Rows, 1>& P) noexcept
{
    using size_type = typename ::vv::cx_matrix<Type, Rows, Rows>::size_type;

    ::vv::cx_matrix<Type, Rows, Rows> result_P{};
    for (size_type i = 0; i < Rows; ++i)
    {
        result_P.at(i, P.at(i, 0)) = static_cast<Type>(1);
    }
    
    return result_P;
}


template <class Type, std::size_t Columns>
constexpr ::vv::cx_matrix<Type, Columns, Columns>
    lupq_get_q_impl(const ::vv::cx_matrix<Type, Columns, 1>& Q) noexcept
{
    using size_type = typename ::vv::cx_matrix<Type, Columns, Columns>::size_type;

    ::vv::cx_matrix<Type, Columns, Columns> result_Q{};
    for (size_type i = 0; i < Columns; ++i)
    {
        result_Q.at(Q.at(i, 0), i) = static_cast<Type>(1);
    }
    return result_Q;
}

} // namespace detail::lupq


/// ===== FUNCTION SECTION =====
template <class Type, std::size_t Rows, std::size_t Columns>
constexpr std::tuple<cx_matrix<Type, Rows, Columns>, cx_matrix<Type, Rows, 1>,
                     cx_matrix<Type, Columns, 1>, std::size_t, long>
    lupq_decompose(cx_matrix<Type, Rows, Columns> mat)
{
    using size_type = typename cx_matrix<Type, Rows, Columns>::size_type;

    long permutations_counter = 0;

    cx_matrix<Type, Rows, 1> P{};
    cx_matrix<Type, Columns, 1> Q{};

    for (size_type i = 0; i < Rows; ++i)
    {
        P.at(i, 0) = static_cast<Type>(i);
    }

    for (size_type i = 0; i < Columns; ++i)
    {
        Q.at(i, 0) = static_cast<Type>(i);
    }

    size_type row_index{};
    size_type column_index{};

    // Calculate cx_matrix C = L + U - E, matrix P, matrix Q, rank and permutations_counter.
    while (!mat.is_empty(row_index, column_index) && row_index < Rows && column_index < Columns)
    {
        const auto indices_pair = mat.find_max_element(row_index, column_index);
        if (row_index != indices_pair.first)
        {
            mat.swap_rows(row_index, indices_pair.first);
            P.swap_rows(row_index, indices_pair.first);
            ++permutations_counter;
        }
        
        if (column_index != indices_pair.second)
        {
            mat.swap_columns(column_index, indices_pair.second);
            Q.swap_rows(column_index, indices_pair.second);
            ++permutations_counter;
        }

        for (size_type j = row_index + 1; j < Rows; ++j)
        {
            mat.at(j, column_index) /= mat.at(row_index, column_index);
            for (size_type k = column_index + 1; k < Columns; ++k)					
            {
                mat.at(j, k) -= mat.at(row_index, k) * mat.at(j, column_index);
            }
        }
        ++row_index;
        ++column_index;
    }
    size_type rank = std::max(row_index, column_index);

    return { mat, P, Q, rank, permutations_counter };
}


template <class Type, std::size_t Rows, std::size_t Columns>
constexpr std::tuple<cx_matrix<Type, Rows, Columns>, cx_matrix<Type, Rows, Rows>,
                     cx_matrix<Type, Rows, Rows>, cx_matrix<Type, Columns, Columns>, long>
    lupq_get_components(const cx_matrix<Type, Rows, Columns>& mat)
{
    const auto [C, P, Q, rank, permutations_counter] = lupq_decompose(mat);

    cx_matrix<Type, Rows, Columns> result_U = detail::lupq::lupq_get_u_impl(C);
    cx_matrix<Type, Rows, Rows> result_L = detail::lupq::lupq_get_l_impl(C);
    cx_matrix<Type, Rows, Rows> result_P = detail::lupq::lupq_get_p_impl(P);
    cx_matrix<Type, Columns, Columns> result_Q = detail::lupq::lupq_get_q_impl(Q);

    return { result_U, result_L, result_P, result_Q, rank };
}


template <class Type, std::size_t Rows, std::size_t Columns_A>
constexpr cx_matrix<Type, Columns_A, 1>
    lupq_solve(const cx_matrix<Type, Rows, Columns_A>& A,
               const cx_matrix<Type, Rows, 1>& b,
               const Type eps = kDefault_eps<Type>)
{
    const auto lupq_tuple = lupq_decompose(A);

    return lupq_solve(A, b, lupq_tuple, eps);
}


template <class Type, std::size_t Rows, std::size_t Columns_A>
constexpr cx_matrix<Type, Columns_A, 1>
    lupq_solve(const cx_matrix<Type, Rows, Columns_A>& A,
               const cx_matrix<Type, Rows, 1>& b,
               const std::tuple<cx_matrix<Type, Rows, Columns_A>, cx_matrix<Type, Rows, 1>,
                                cx_matrix<Type, Columns_A, 1>, std::size_t, long>& lupq_tuple,
               const Type eps = kDefault_eps<Type>)
{
    using size_type = typename cx_matrix<Type, Rows, Columns_A>::size_type;

    const auto [C, _P, _Q, rank, permutations_counter] = lupq_tuple;
    cx_matrix<Type, Rows, Rows> P = detail::lupq::lupq_get_p_impl(_P);
    cx_matrix<Type, Columns_A, Columns_A> Q = detail::lupq::lupq_get_q_impl(_Q);

    cx_matrix<Type, Columns_A, 1> solution{};
    solution = P * b;

    // Solve equation Lx = b.
    cx_matrix<Type, Rows, Rows> L = detail::lupq::lupq_get_l_solution_impl(C, rank);
    for (size_type i = 0; i < rank - 1; ++i)
    {
        for (size_type j = i + 1; j < rank; ++j)
        {
            solution.at(j, 0) -= solution.at(i, 0) * L.at(j, i);
        }
    }
    
    // Solve system Uy = b, xQ = y.
    cx_matrix<Type, Rows, Columns_A> U = detail::lupq::lupq_get_u_solution_impl(C, rank);
    for (size_type i = rank - 1; i > 0; --i)
    {
        solution.at(i, 0) /= U.at(i, i);

        for (size_type j = i - 1; ; --j)
        {
            solution.at(j, 0) -= solution.at(i, 0) * U.at(j, i);

            if (j == 0) break;
        }
    }
    solution.at(0, 0) /= U.at(0, 0); // To make U(0, 0) == Type{1}

    // Checking the compatibility of the system of equations.
    if (rank < Rows)
    {
        cx_matrix<Type, Rows, Columns_A> origin = P * A * Q;
        for (size_type i = rank; i < Rows && i < Columns_A; ++i)
        {
            Type x = solution.at(i, 0);
            for (size_type j = 0; j < rank && j < Columns_A; ++j)
            {
                x -= origin.at(i, j) * solution.at(j, 0);
            }

            if (cx::abs(x) > eps)
            {
                return cx_matrix<Type, Rows, 1>::get_error_matrix();
            }
        }
    }

    for (size_type i = rank; i < Columns_A; ++i)
    {
        solution.at(i, 0) = Type{};
    }

    return Q * solution;
}


/// ===== TEST SECTION =====
#if ENABLE_TESTS_LUPQ

TEST_METHOD(lupq_test_decompose)
{
    std::cout << "LUPQ-decompose:\n\n";
    
    constexpr auto lupq_tuple = lupq_get_components(mat_A);
    constexpr auto mat_U = std::get<0>(lupq_tuple);
    constexpr auto mat_L = std::get<1>(lupq_tuple);
    constexpr auto mat_P = std::get<2>(lupq_tuple);
    constexpr auto mat_Q = std::get<3>(lupq_tuple);
    constexpr auto rank = std::get<4>(lupq_tuple);
    std::cout << "Matrix A:\n" << mat_A << "\n\n";
    std::cout << "Matrix U:\n" << mat_U << "\n\n";
    std::cout << "Matrix L:\n" << mat_L << "\n\n";
    std::cout << "Matrix P:\n" << mat_P << "\n\n";
    std::cout << "Matrix Q:\n" << mat_Q << "\n\n";
    std::cout << "Rank A:\n" << rank;
    std::cout << "\n\n------------------------------------\n\n";

    std::cout << "Checking LUPQ-decompose:\n\n";
    constexpr auto lupq_value_check = (mat_L * mat_U) - (mat_P * mat_A * mat_Q);
    std::cout << "Matrix (L * U) - (P * A * Q):\n" << lupq_value_check;
    std::cout << "\n\n------------------------------------\n\n";
}


TEST_METHOD(lupq_test_solve)
{
    std::cout << "LUPQ-solve:\n\n";
    
    constexpr auto vec_x = lupq_solve(mat_A, vec_b);
    std::cout << "Matrix A:\n" << mat_A << "\n\n";
    std::cout << "Matrix b:\n" << vec_b << "\n\n";
    std::cout << "Matrix X:\n" << vec_x;
    std::cout << "\n\n------------------------------------\n\n";
    
    std::cout << "Checking LUPQ-solve:\n\n";
    constexpr auto lupq_solve_check = (mat_A * vec_x) - vec_b;
    std::cout << "Result of (A * x - b):\n" << lupq_solve_check;
    std::cout << "\n\n------------------------------------\n\n";
}

#endif // ENABLE_TESTS_LUPQ

} //namesace vv
