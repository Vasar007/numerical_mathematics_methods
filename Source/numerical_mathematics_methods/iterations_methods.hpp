#pragma once

#include <algorithm>
#include <array>
#include <stdexcept>
#include <type_traits>
#include <thread>
#include <chrono>

#include "utils.hpp"
#include "cx_loops.hpp"
#include "cx_random.hpp"
#include "cx_matrix.hpp"
#include "constants.hpp"
#include "cx_math.h"


namespace vv
{

/// ===== CONSTANT SECTION =====
namespace
{

constexpr long kMax_iterations_methods = 10'000;

} // anonymous namespace


/// ===== FUNCTION SECTION =====
template <class Type, std::size_t Rows, std::size_t Columns_A, std::size_t Columns_b>
constexpr std::pair<cx_matrix<Type, Rows, Columns_b>, long>
    jacobi_solve(const cx_matrix<Type, Rows, Columns_A>& A,
                 const cx_matrix<Type, Rows, Columns_b>& b,
                 const Type eps = kDefault_eps<Type>) noexcept
{
    static_assert(Rows == Columns_A, "Matrix is not quadrant!");
    static_assert(Columns_b == 1, "Matrix contains more than one columns in right "
                                  "hand side vector!");

    using size_type = typename cx_matrix<Type, Rows, Columns_A>::size_type;
    constexpr std::size_t N = Rows;

    long iterations_counter = 0;

    Type norm{};
    cx_matrix<Type, Rows, 1> x{};
    cx_matrix<Type, Rows, 1> previous_x{};
    do
    {
        for (size_type i = 0; i < N; ++i)
        {
            previous_x(i, 0u) = b(i, 0);
            for (size_type j = 0; j < N; ++j)
            {
                if (i != j)
                {
                    if (cx::abs(x(j, 0)) > eps)
                    {
                        previous_x(i, 0) -= A(i, j) * x(j, 0);
                    }
                }
            }
            previous_x(i, 0) /= A(i, i);
        }

        norm = cx::abs(x(0, 0) - previous_x(0, 0));
        for (size_type k = 0; k < N; ++k)
        {
            if (cx::abs(x(k, 0) - previous_x(k, 0)) > norm)
            {
                norm = cx::abs(x(k, 0) - previous_x(k, 0));
            }
            x(k, 0) = previous_x(k, 0);
        }

        ++iterations_counter;
    }
    while (norm > eps && iterations_counter < kMax_iterations_methods);

    return { x, iterations_counter };
}


template <class Type, std::size_t Rows, std::size_t Columns_A, std::size_t Columns_b>
constexpr std::pair<cx_matrix<Type, Rows, Columns_b>, long>
    seidel_solve(const cx_matrix<Type, Rows, Columns_A>& A,
                 const cx_matrix<Type, Rows, Columns_b>& b,
                 const Type eps = kDefault_eps<Type>) noexcept
{
    static_assert(Rows == Columns_A, "Matrix is not quadrant!");
    static_assert(Columns_b == 1, "Matrix contains more than one columns in right "
                                  "hand side vector!");

    using size_type = typename cx_matrix<Type, Rows, Columns_A>::size_type;
    constexpr std::size_t N = Rows;

    long iterations_counter = 0;

    constexpr auto converge = [](const auto& x_k, const auto& x_k_prev, const Type epsilon)
    {
        Type norm{};
        for (size_type i = 0; i < N; ++i)
        {
            norm += (x_k(i, 0) - x_k_prev(i, 0)) * (x_k(i, 0) - x_k_prev(i, 0));
        }

        return (cx::sqrt(norm) < epsilon);
    };

    cx_matrix<Type, Rows, 1> x{};
    cx_matrix<Type, Rows, 1> previous_x{};
    do
    {
        previous_x = x;

        for (size_type i = 0; i < N; ++i)
        {
            Type var{};
            for (size_type j = 0; j < i; ++j)
            {
                var += (A(i, j) * x(j, 0));
            }

            for (size_type j = i + 1; j < N; ++j)
            {
                var += (A(i, j) * previous_x(j, 0));
            }

            x(i, 0) = (b(i, 0) - var) / A(i, i);
        }

        ++iterations_counter;
    }
    while (!converge(x, previous_x, eps) && iterations_counter < kMax_iterations_methods);

    return { x, iterations_counter };
}


template <class Type, std::size_t Rows, std::size_t Columns>
constexpr cx_matrix<Type, Rows, Columns> generate_matrix() noexcept
{
    // Alias the size_type.
    using size_type = typename cx_matrix<Type, Rows, Columns>::size_type;

    cx_matrix<Type, Rows, Columns> result{};
    for (size_type i = 0; i < Rows; ++i)
    {
        for (size_type j = 0; j < Columns; ++j)
        {
            const Type number = static_cast<Type>(cx_random::get_random(i * Rows + j + 1));

            result(i, j) += number - static_cast<Type>(i * Rows + j + 1);
            if (i + 1u < Rows)
            {
                result(i + 1, j) = number * static_cast<Type>(0.001);
            }
        }
    }

    return result;
}


template <class Type, std::size_t Rows, std::size_t Columns>
constexpr cx_matrix<Type, Rows, Columns> generate_matrix_with_diagonal_predominance() noexcept
{
    static_assert(Rows == Columns, "Matrix is not quadrant!");

    // Alias the size_type.
    using size_type = typename cx_matrix<Type, Rows, Columns>::size_type;

    cx_matrix<Type, Rows, Columns> result{};
    for (size_type i = 0; i < Rows; ++i)
    {
        Type sum{};
        for (size_type j = 0; j < Columns; ++j)
        {
            const Type number = static_cast<Type>(cx_random::get_random(i * Rows + j));

            result(i, j) += number - static_cast<Type>(i * Rows + j);
            if (i + 1 < Rows)
            {
                result(i + 1, j) = number * static_cast<Type>(0.001);
            }
            sum += result(i, j);
        }

        result(i, i) += sum;
    }

    return result;
}


template <class Type, std::size_t Rows, std::size_t Columns>
constexpr cx_matrix<Type, Rows, Columns> generate_matrix_simmetrical() noexcept
{
    auto result = generate_matrix<Type, Rows, Columns>();
    result = result * result.transpose();

    return result;
}


/// ===== TEST SECTION =====
#if ENABLE_TESTS_ITERATIONS

TEST_METHOD(jacobi_test_solve)
{
    std::cout << "Matrix solved with Jacobi (include eps = " << kEps << "):\n\n";
    
    constexpr auto vec_x = jacobi_solve(mat_A, vec_b, kEps);
    std::cout << "Matrix x:\n" << vec_x.first << "\n\n";
    std::cout << "Number of iterations Jacobi:\n" << vec_x.second;
    std::cout << "\n\n------------------------------------\n\n";

    std::cout << "Checking Jacobi solve (include eps = " << kEps << "):\n\n";
    constexpr auto jacobi_solve_check = (mat_A * vec_x.first) - vec_b;
    std::cout << "Result of (A * x - b):\n" << jacobi_solve_check;
    std::cout << "\n\n------------------------------------\n\n";
}

TEST_METHOD(seidel_test_solve)
{
    std::cout << "Matrix solved with Seidel (include eps = " << kEps << "):\n\n";
    
    constexpr auto vec_x = seidel_solve(mat_A, vec_b, kEps);
    std::cout << "Matrix x:\n" << vec_x.first << "\n\n";
    std::cout << "Number of iterations Seidel:\n" << vec_x.second;
    std::cout << "\n\n------------------------------------\n\n";

    std::cout << "Checking Seidel solve (include eps = " << kEps << "):\n\n";
    constexpr auto seidel_solve_check = (mat_A * vec_x.first) - vec_b;
    std::cout << "Result of (A * x - b):\n" << seidel_solve_check;
    std::cout << "\n\n------------------------------------\n\n";
}

#endif // ENABLE_TESTS_ITERATIONS

#if ENABLE_TESTS_ITERATIONS_DIAGONAL

TEST_METHOD(gen_test_diagonal_predominance)
{
    std::cout << "Diagonal predominance cx_matrix solved with Seidel and Jacobi (include eps = " << kEps << "):\n\n";
    
    cx_loops::static_for<10>([] (const auto Index)
    {
        if constexpr (Index > 2)
        {
            constexpr auto A = generate_matrix_with_diagonal_predominance<double, Index, Index>();
            constexpr auto b = generate_matrix<double, Index, 1>();
            constexpr auto x1 = jacobi_solve(A, b, kEps);
            constexpr auto x2 = seidel_solve(A, b, kEps);
            std::cout << "Matrix A" << Index << ":\n" << A << "\n\n";
            std::cout << "Matrix b" << Index << ":\n" << b << "\n\n";
            std::cout << "Matrix x" << Index << ":\n" << x1.first << "\n\n";
            std::cout << "Number of iterations Jacobi" << Index << ":\n" << x1.second << "\n\n";
            std::cout << "Matrix x" << Index << ":\n" << x2.first << "\n\n";
            std::cout << "Number of iterations Seidel" << Index << ":\n" << x2.second;
            std::cout << "\n\n------------------------------------\n\n";

            if constexpr (Index % 5 == 0)
            {
                utils::pause_clear("\nPress Enter to continue ELIMINATION...");
            }            
        }
        else
        {
            using namespace std::chrono_literals;
            std::cout << "Initiate ELIMINATION " << (3 - Index) << "...\n";
            std::this_thread::sleep_for(1s);
        }
    });
}

#endif // ENABLE_TESTS_ITERATIONS_DIAGONAL

#if ENABLE_TESTS_ITERATONS_SIMMETRICAL

TEST_METHOD(gen_test_simmetrical)
{
    std::cout << "Simmetrical cx_matrix solved with Seidel and Jacobi (include eps = " << kEps << "):\n\n";
    
    cx_loops::static_for<10>([] (const auto Index)
    {
        if constexpr (Index > 2)
        {
            constexpr auto A = generate_matrix_simmetrical<double, Index, Index>();
            constexpr auto b = generate_matrix<double, Index, 1>();
            constexpr auto x1 = jacobi_solve(A, b, kEps);
            constexpr auto x2 = seidel_solve(A, b, kEps);
            std::cout << "Matrix A" << Index << ":\n" << A << "\n\n";
            std::cout << "Matrix b" << Index << ":\n" << b << "\n\n";
            std::cout << "Matrix x" << Index << ":\n" << x1.first << "\n\n";
            std::cout << "Number of iterations Jacobi" << Index << ":\n" << x1.second << "\n\n";
            std::cout << "Matrix x" << Index << ":\n" << x2.first << "\n\n";
            std::cout << "Number of iterations Seidel" << Index << ":\n" << x2.second;
            std::cout << "\n\n------------------------------------\n\n";

            if constexpr (Index % 5 == 0)
            {
                utils::pause_clear("\nPress Enter to continue ELIMINATION...");
            }            
        }
        else
        {
            using namespace std::chrono_literals;
            std::cout << "Initiate ELIMINATION " << (3 - Index) << "...\n";
            std::this_thread::sleep_for(1s);
        }
    });
}

#endif // ENABLE_TESTS_ITERATONS_SIMMETRICAL

} // namespace vv
