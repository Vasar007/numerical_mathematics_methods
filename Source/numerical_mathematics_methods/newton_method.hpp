#pragma once

#include <array>
#include <tuple>
#include <chrono>
#include <cassert>
#include <cmath>
#include <vector>

#include "cx_math.h"


namespace vv
{

/// ===== CONSTANT SECTION =====
namespace 
{

constexpr long kMax_iterations_newton = 1'000;

template <long N>
constexpr long kComplexity_lupq = 2 * N * N * N / 3;

template <long N>
constexpr long kComplexity_solve = N * N;

template <long N>
constexpr long kComplexity_calc_jacobi = N * N;

template <long N>
constexpr long kComplexity_calc_F = N;

constexpr long kMax = 10;

} // anonymous namespace

// ===== Data 1-1 (Nonlinear system) =====
// Not using constexpr because of C-functions <cmath> calculate better then <cx_math.h>
cx_matrix<double, 10, 1> nonlinear_matrix(const cx_matrix<double, 10, 1>& x)
{
    cx_matrix<double, 10, 1> result
    {
        cos(x(0, 0) * x(1, 0)) - exp(-3. * x(2, 0)) + x(3, 0) * x(4, 0) * x(4, 0) -
            x(5, 0) - sinh(2. * x(7, 0)) * x(8, 0) + 2. * x(9, 0) + 2.0004339741653854440,

        sin(x(0, 0) * x(1, 0)) + x(2, 0) * x(8, 0) * x(6, 0) - exp(-x(9, 0) + x(5, 0)) +
            3. * x(4, 0) * x(4, 0) - x(5, 0) * (x(7, 0) + 1.) + 10.886272036407019994,
        
        x(0, 0) - x(1, 0) + x(2, 0) - x(3, 0) + x(4, 0) - x(5, 0) + x(6, 0) - x(7, 0) + x(8, 0) -
            x(9, 0) - 3.1361904761904761904,
        
        2. * cos(-x(8, 0) + x(3, 0)) + x(4, 0) / (x(2, 0) + x(0, 0)) - sin(x(1, 0) *
            x(1, 0)) + cos(x(6, 0) * x(9, 0)) * cos(x(6, 0) * x(9, 0)) - x(7, 0) -
            0.170747270502230475,
        
        sin(x(4, 0)) + 2. * x(7, 0) * (x(2, 0) + x(0, 0)) - exp(-x(6, 0) *
            (-x(9, 0) + x(5, 0))) + 2. * cos(x(1, 0)) - 1. / (x(3, 0) - x(8, 0)) -
            0.368589627310127786,
        
        exp(x(0, 0) - x(3, 0) - x(8, 0)) + x(4, 0) * x(4, 0) / x(7, 0) + 0.5 *
            cos(3. * x(9, 0) * x(1, 0)) - x(5, 0) * x(2, 0) + 2.049108601677187511,
        
        x(1, 0) * x(1, 0) * x(1, 0) * x(6, 0) - sin(x(9, 0) / x(4, 0) + x(7, 0)) +
            (x(0, 0) - x(5, 0)) * cos(x(3, 0)) + x(2, 0) - 0.738043007620279801,
        
        x(4, 0) * (x(0, 0) - 2. * x(5, 0)) * (x(0, 0) - 2. * x(5, 0)) - 2. * sin(-x(8, 0) +
            x(2, 0)) + 1.5 * x(3, 0) - exp(x(1, 0) * x(6, 0) + x(9, 0)) + 3.566832198969380904,
        
        7. / x(5, 0) + exp(x(4, 0) + x(3, 0)) - 2 * x(1, 0) * x(7, 0) * x(9, 0) * x(6, 0) +
            3 * x(8, 0) - 3 * x(0, 0) - 8.439473450838325749,
        
        x(9, 0) * x(0, 0) + x(8, 0) * x(1, 0) - x(7, 0) * x(2, 0) + sin(x(3, 0) + x(4, 0) +
            x(5, 0)) * x(6, 0) - 0.7823809523809523809
    };
    return result;
}

constexpr cx_matrix<double, 10, 1> vector_b
{
     -2.0004339741653854440,
    -10.886272036407019994,
      3.1361904761904761904,
      0.170747270502230475,
      0.368589627310127786,
     -2.049108601677187511,
      0.738043007620279801,
     -3.566832198969380904,
      8.439473450838325749,
      0.7823809523809523809
};

// ===== Data 1-2-1 (vector x0) =====
constexpr cx_matrix<double, 10, 1> vector_x
{
     0.5,
     0.5,
     1.5,
    -1.0,
    -0.5,
     1.5,
     0.5,
    -0.5,
     1.5,
    -1.5
};

// ===== Data 1-2-2 (vector x0, extra task) =====
constexpr cx_matrix<double, 10, 1> vector_x_extra
{
     0.5,
     0.5,
     1.5,
    -1.0,
    -0.2,
     1.5,
     0.5,
    -0.5,
     1.5,
    -1.5
};

// ===== Data 1-3 (Jacobi cx_matrix) =====
cx_matrix<double, 10, 10> jacobi_matrix(const cx_matrix<double, 10, 1>& x)
{
    cx_matrix<double, 10, 10> result
    {
        -sin(x(0, 0) * x(1, 0)) * x(1, 0),
        -sin(x(0, 0) * x(1, 0)) * x(0, 0),
        3. * exp(-3. * x(2, 0)),
        x(4, 0) * x(4, 0),
        2. * x(3, 0) * x(4, 0),
        -1.,
        0.,
        -2. * cosh(2. * x(7, 0)) * x(8, 0),
        -sinh(2. * x(7, 0)),
        2.,
        cos(x(0, 0) * x(1, 0)) * x(1, 0),
        cos(x(0, 0) * x(1, 0)) * x(0, 0),
        x(8, 0) * x(6, 0),
        0.,
        6. * x(4, 0),
        -exp(-x(9, 0) + x(5, 0)) - x(7, 0) - 1.,
        x(2, 0) * x(8, 0),
        -x(5, 0),
        x(2, 0) * x(6, 0),
        exp(-x(9, 0) + x(5, 0)),
        1.,
        -1.,
        1.,
        -1.,
        1.,
        -1.,
        1.,
        -1.,
        1.,
        -1.,
        -x(4, 0) / ((x(2, 0) + x(0, 0)) * (x(2, 0) + x(0, 0))),
        -2. * cos(x(1, 0) * x(1, 0)) * x(1, 0),
        -x(4, 0) / ((x(2, 0) + x(0, 0)) * (x(2, 0) + x(0, 0))),
        -2. * sin(-x(8, 0) + x(3, 0)),
        1. / (x(2, 0) + x(0, 0)),
        0.,
        -2. * cos(x(6, 0) * x(9, 0)) * sin(x(6, 0) * x(9, 0)) * x(9, 0),
        -1,
        2. * sin(-x(8, 0) + x(3, 0)),
        -2. * cos(x(6, 0) * x(9, 0)) * sin(x(6, 0) * x(9, 0)) * x(6, 0),
        2. * x(7, 0),
        -2. * sin(x(1, 0)),
        2. * x(7, 0),
        1. / ((-x(8, 0) + x(3, 0)) * (-x(8, 0) + x(3, 0))),
        cos(x(4, 0)),
        x(6, 0) * exp(-x(6, 0) * (-x(9, 0) + x(5, 0))),
        -(x(9, 0) - x(5, 0)) * exp(-x(6, 0) * (-x(9, 0) + x(5, 0))),
        2. * x(2, 0) + 2. * x(0, 0),
        -1. / ((-x(8, 0) + x(3, 0)) * (-x(8, 0) + x(3, 0))),
        -x(6, 0) * exp(-x(6, 0) * (-x(9, 0) + x(5, 0))),
        exp(x(0, 0) - x(3, 0) - x(8, 0)),
        -1.5 * sin(3. * x(9, 0) * x(1, 0)) * x(9, 0),
        -x(5, 0),
        -exp(x(0, 0) - x(3, 0) - x(8, 0)),
        2. * x(4, 0) / x(7, 0),
        -x(2, 0),
        0.,
        -x(4, 0) * x(4, 0) / (x(7, 0) * x(7, 0)),
        -exp(x(0, 0) - x(3, 0) - x(8, 0)),
        -1.5 * sin(3. * x(9, 0) * x(1, 0)) * x(1, 0),
        cos(x(3, 0)),
        3. * x(1, 0) * x(1, 0) * x(6, 0),
        1.,
        -(x(0, 0) - x(5, 0)) * sin(x(3, 0)),
        cos(x(9, 0) / x(4, 0) + x(7, 0)) * x(9, 0) / (x(4, 0) * x(4, 0)),
        -cos(x(3, 0)),
        x(1, 0) * x(1, 0) * x(1, 0),
        -cos(x(9, 0) / x(4, 0) + x(7, 0)),
        0.,
        -cos(x(9, 0) / x(4, 0) + x(7, 0)) / x(4, 0),
        2. * x(4, 0) * (x(0, 0) - 2. * x(5, 0)),
        -x(6, 0) * exp(x(1, 0) * x(6, 0) + x(9, 0)),
        -2. * cos(-x(8, 0) + x(2, 0)),
        1.5,
        (x(0, 0) - 2. * x(5, 0)) * (x(0, 0) - 2. * x(5, 0)),
        -4. * x(4, 0) * (x(0, 0) - 2. * x(5, 0)),
        -x(1, 0) * exp(x(1, 0) * x(6, 0) + x(9, 0)),
        0.,
        2. * cos(-x(8, 0) + x(2, 0)),
        -exp(x(1, 0) * x(6, 0) + x(9, 0)),
        -3.,
        -2. * x(7, 0) * x(9, 0) * x(6, 0),
        0.,
        exp(x(4, 0) + x(3, 0)),
        exp(x(4, 0) + x(3, 0)),
        -7. / (x(5, 0) * x(5, 0)),
        -2. * x(1, 0) * x(7, 0) * x(9, 0),
        -2. * x(1, 0) * x(9, 0) * x(6, 0),
        3.,
        -2. * x(1, 0) * x(7, 0) * x(6, 0),
        x(9, 0),
        x(8, 0),
        -x(7, 0),
        cos(x(3, 0) + x(4, 0) + x(5, 0)) * x(6, 0),
        cos(x(3, 0) + x(4, 0) + x(5, 0)) * x(6, 0),
        cos(x(3, 0) + x(4, 0) + x(5, 0)) * x(6, 0),
        sin(x(3, 0) + x(4, 0) + x(5, 0)),
        -x(2, 0),
        x(1, 0),
        x(0, 0)
    };
    return result;
}

// ===== Data 2 (Test system) =====
cx_matrix<double, 2, 1> test_system(const cx_matrix<double, 2, 1>& x)
{
    cx_matrix<double, 2, 1> result
    {
        sin(2. * x(0, 0) - x(1, 0)) - 1.2 * x(0, 0) - 0.4,
        0.8 * x(0, 0) * x(0, 0) + 1.5 * x(1, 0) * x(1, 0) - 1.
    };
    return result;
}

cx_matrix<double, 2, 2> test_diff(const cx_matrix<double, 2, 1>& x)
{
    cx_matrix<double, 2, 2> result
    {
        2. * cos(2. * x(0, 0) - x(1, 0)) - 1.2,
        1.6 * x(0, 0),
        -cos(2. * x(0, 0) - x(1, 0)),
        3. * x(1, 0),
    };
    return result;
}

constexpr cx_matrix<double, 2, 1> test_x
{
     0.4,
    -0.75
};


/// ===== FUNCTION SECTION =====
template <class Type, class Function>
constexpr Type bisection_method(Type xn, Type xk, Function&& f,
                                const Type eps = kDefault_eps<Type>)
{
    if (f(xn) == static_cast<Type>(0))
    {
        return xn;
    }
    if (f(xk) == static_cast<Type>(0))
    {
        return xk;
    }

    Type xi{};
    while (xk - xn > eps)
    {
        Type dx = (xk - xn) / 2.0;
        xi = xn + dx;
        if (f(xn) * f(xk) < static_cast<Type>(0))
        {
            xk = xi;
        }
        else
        {
            xn = xi;
        }              
    }
    return xi;
}


template <class Type, class Function>
constexpr std::pair<Type, Type> bisection_localize(Type xn, Type xk, Function&& f,
                                                   const Type eps = kRough_eps<Type>)
{
    Type xi{};
    while (xk - xn > eps)
    {
        xi = (xn + xk) / 2.0;
        if (f(xk) * f(xi) < static_cast<Type>(0))
        {
            xn = xi;
            break;
        }
        
        xk = xi;
    }
    return { xn, xi };
}


template <class Type, class Function, class Derivative>
constexpr std::tuple<Type, Type, long> newton_method(Function&& f, Derivative&& df, Type x0,
                                                     const Type eps = kDefault_eps<Type>)
{
    Type x1 = x0 - f(x0) / df(x0);
    long iterations_counter = 1;

    while (cx::abs(x1 - x0) > eps && iterations_counter < kMax_iterations_newton)
    {
        x0 = x1;
        x1 -= f(x1) / df(x1);
        ++iterations_counter;
    }
    return { x1, x0, iterations_counter };
}


template <class Type, class Function, class Derivative>
constexpr std::tuple<Type, Type, long> newton_method(Function&& f, Derivative&& df, Type xn, Type xk,
                                               const Type eps = kDefault_eps<Type>)
{
    auto [x0, _x1] = bisection_localize(xn, xk, f);

    Type x1 = x0 - f(x0) / df(x0);
    long iterations_counter = 1;

    while (cx::abs(x1 - x0) > eps && iterations_counter < kMax_iterations_newton)
    {
        x0 = x1;
        x1 -= f(x1) / df(x1);
        ++iterations_counter;
    }
    return { x1, x0, iterations_counter };
}


template <class Type, std::size_t Row, std::size_t Columns = 1>
constexpr Type difference(const cx_matrix<Type, Row, Columns>& lhs,
                          const cx_matrix<Type, Row, Columns>& rhs) noexcept
{
    static_assert(Columns == 1, "Vectors have more than one columns!");

    Type result{};
    for (std::size_t i = 0; i < Row; ++i)
    {
        result += lhs(i, 0) - rhs(i, 0);
        //result = std::max(cx::abs(lhs(i, 0) - rhs(i, 0)), result);
    }
    return cx::abs(result);
}


template <class Type, std::size_t N, class Function, class Derivative>
constexpr std::tuple<cx_matrix<Type, N, 1>, long, long>
    newton_method_system(Function&& F, Derivative&& dF, cx_matrix<Type, N, 1> x0,
                         const Type eps = kDefault_eps<Type>)
{
    cx_matrix<Type, N, 1> x1 = x0 + lupq_solve(dF(x0), -1.0 * F(x0), eps);
    long iterations_counter = 1;
    long complexity = kComplexity_calc_jacobi<N> + kComplexity_calc_F<N> +
                      kComplexity_lupq<N> + kComplexity_solve<N>;

    while (difference(x1, x0) > eps && iterations_counter < kMax_iterations_newton)
    {
        x0 = x1;
        x1 += lupq_solve(dF(x0), -1.0 * F(x0), eps);
        ++iterations_counter;
        complexity += kComplexity_calc_jacobi<N> + kComplexity_calc_F<N> +
                      kComplexity_lupq<N> + kComplexity_solve<N>;
    }
    return { x1, iterations_counter, complexity };
}


template <class Type, std::size_t N, class Function, class Derivative>
constexpr std::tuple<cx_matrix<Type, N, 1>, long, long>
    mod_newton_method_system(Function&& F, Derivative&& dF, cx_matrix<Type, N, 1> x0, long k,
                             const Type eps = kDefault_eps<Type>)
{
    assert(k > 0);

    cx_matrix<Type, N, N> calc_jacobi = dF(x0);
    --k;
    long complexity = kComplexity_calc_jacobi<N>;

    auto lupq_tuple = lupq_decompose(calc_jacobi, eps);
    complexity += kComplexity_lupq<N>;

    cx_matrix<Type, N, 1> x1 = x0 + lupq_solve(calc_jacobi, -1.0 * F(x0), lupq_tuple, eps);
    long iterations_counter = 1;
    complexity += kComplexity_solve<N> + kComplexity_calc_F<N>;

    while (difference(x1, x0) > eps && iterations_counter < kMax_iterations_newton)
    {
        x0 = x1;

        if (k != 0)
        {
            --k;
            calc_jacobi = dF(x0);
            lupq_tuple = lupq_decompose(calc_jacobi, eps);
            complexity += kComplexity_calc_jacobi<N> + kComplexity_lupq<N>;
        }

        x1 += lupq_solve(calc_jacobi, -1.0 * F(x0), lupq_tuple, eps);
        ++iterations_counter;
        complexity += kComplexity_solve<N> + kComplexity_calc_F<N>;
    }
    return { x1, iterations_counter, complexity };
}


// TODO: fix cycling in automatic calculations.
template <class Type, std::size_t N, class Function, class Derivative>
std::tuple<cx_matrix<Type, N, 1>, long, long>
    mod_newton_method_system(Function&& F, Derivative&& dF, cx_matrix<Type, N, 1> x0, long k,
                             const bool, const Type eps = kDefault_eps<Type>)
{
    assert(k > 0);

    cx_matrix<Type, N, N> calc_jacobi = dF(x0);
    --k;
    long complexity = kComplexity_calc_jacobi<N>;

    auto lupq_tuple = lupq_decompose(calc_jacobi, eps);
    complexity += kComplexity_lupq<N>;

    cx_matrix<Type, N, 1> x1 = x0 + lupq_solve(calc_jacobi, -1.0 * F(x0), lupq_tuple, eps);
    long iterations_counter = 1;
    complexity += kComplexity_solve<N> + kComplexity_calc_F<N>;

    Type deviation = difference(x1, x0);
    Type prev_deviation = difference(x1, x0);
    std::vector<cx_matrix<Type, N, 1>> prev_x{ x0 };
    cx_matrix<Type, N, 1> prev_x0 = x0;
    prev_x.reserve(kMax_iterations_newton / 10);
    bool diverge = false;
    bool again_diverge = false;
    long cycling_heler = 0;
    while (deviation > eps && iterations_counter < kMax_iterations_newton)
    {
        x0 = x1;

        if (k != 0 || diverge)
        {
            k -= k == 0 ? 0 : 1;
            calc_jacobi = dF(x0);
            lupq_tuple = lupq_decompose(calc_jacobi, eps);
            complexity += kComplexity_calc_jacobi<N> + kComplexity_lupq<N>;
        }

        x1 += lupq_solve(calc_jacobi, -1.0 * F(x0), lupq_tuple, eps);
        ++iterations_counter;
        complexity += kComplexity_solve<N> + kComplexity_calc_F<N>;

        deviation = difference(x1, x0);
        if (deviation > prev_deviation)
        {
            std::cout << iterations_counter << " Diverge!\n";

            if (diverge)
            {
                again_diverge = true;
                if (!prev_x.empty())
                {
                    std::cout << "Diverge again!\n";
                    prev_x0 = std::move(prev_x.back());
                    prev_x.pop_back();
                }
            }
            diverge = true;

            // Not sure that it's work. Begin.
            if (cycling_heler == iterations_counter - 2)
            {
                if (!prev_x.empty())
                {
                    std::cout << "Cycling detected!\n";
                    prev_x.pop_back();
                    prev_x0 = std::move(prev_x.back());
                    prev_x.pop_back();
                }
            }
            cycling_heler = iterations_counter;
            // End.

            x1 = prev_x0;            
        }
        else
        {
            diverge = false;
            again_diverge = false;
            prev_x0 = x0;            
            prev_x.emplace_back(std::move(x0));
        }
        prev_deviation = deviation;
    }
    return { x1, iterations_counter, complexity };
}


template <class Type, std::size_t N, class Function, class Derivative>
constexpr std::tuple<cx_matrix<Type, N, 1>, long, long>
    hybrid_newton_method_system(Function&& F, Derivative&& dF, cx_matrix<Type, N, 1> x0,
                                const long k, const Type eps = kDefault_eps<Type>)
{
    assert(k > 0);

    cx_matrix<Type, N, N> calc_jacobi = dF(x0);
    long k_iter = 1;
    long complexity = kComplexity_calc_jacobi<N>;

    auto lupq_tuple = lupq_decompose(calc_jacobi, eps);
    complexity += kComplexity_lupq<N>;

    cx_matrix<Type, N, 1> x1 = x0 + lupq_solve(calc_jacobi, -1.0 * F(x0), lupq_tuple, eps);
    long iterations_counter = 1;
    complexity += kComplexity_solve<N> + kComplexity_calc_F<N>;

    while (difference(x1, x0) > eps && iterations_counter < kMax_iterations_newton)
    {
        x0 = x1;

        if (k_iter == k)
        {
            calc_jacobi = dF(x0);
            lupq_tuple = lupq_decompose(calc_jacobi, eps);
            k_iter = 1;
            complexity += kComplexity_calc_jacobi<N> + kComplexity_lupq<N>;
        }
        else
        {
            ++k_iter;
        }

        x1 += lupq_solve(calc_jacobi, -1.0 * F(x0), eps);
        ++iterations_counter;
        complexity += kComplexity_solve<N> + kComplexity_calc_F<N>;
    }
    return { x1, iterations_counter, complexity };
}


/// ===== TEST SECTION =====
#if ENABLE_TESTS_NEWTON

TEST_METHOD(newton_test_solve)
{
    std::cout << "Newton method for nonlinear equation:\n\n";

    // Function and derivative.
    constexpr auto f = [](const double x) { return cx::tan(x) - cx::cos(x) + 0.1; };
    constexpr auto df = [](const double x)
        { return (1.0 / (cx::cos(x) * cx::cos(x))) + cx::sin(x); };

    constexpr auto newton_solve = newton_method(f, df, -1.0, 1.0, kEps);
    constexpr auto answer = std::get<0>(newton_solve);
    constexpr auto x0 = std::get<1>(newton_solve);    
    constexpr auto iterations = std::get<2>(newton_solve);
    std::cout << "Epsilon = " << kEps << '\n';
    std::cout << "x0 = " << x0 << '\n';
    std::cout << "Answer = " << answer << '\n';
    std::cout << "Number of iterations = " << iterations;
    std::cout << "\n\n------------------------------------\n\n";

    std::cout << "Checking Newton solve (include eps = " << kEps << "):\n\n";
    constexpr auto newton_solve_check = f(answer);
    std::cout << "Result of f(Answer): " << newton_solve_check;
    std::cout << "\n\n------------------------------------\n\n";
}

#endif // ENABLE_TESTS_NEWTON

#if ENABLE_TESTS_NEWTON_SYSTEM

TEST_METHOD(newton_test_solve_system)
{
    std::cout << "Newton method for nonlinear system:\n\n";

    auto start = std::chrono::steady_clock::now();
    const auto newton_solve_system = newton_method_system(nonlinear_matrix, jacobi_matrix,
                                                          vector_x, kEps);
    std::chrono::duration<double> duration = std::chrono::steady_clock::now() - start;
    std::cout << "Calculation time: " << duration.count() << " ms\n";

    const auto answer = std::get<0>(newton_solve_system);
    const auto iterations = std::get<1>(newton_solve_system);
    const auto complexity = std::get<2>(newton_solve_system);
    std::cout << "Epsilon = " << kEps << '\n';
    std::cout << "x0:\n" << vector_x << "\n\n";
    std::cout << "Answer:\n" << answer << "\n\n";
    std::cout << "Number of iterations = " << iterations << '\n';
    std::cout << "Number of operations = " << complexity;
    std::cout << "\n\n------------------------------------\n\n";

    std::cout << "Checking Newton solve (include eps = " << kEps << "):\n\n";
    const auto newton_solve_check = nonlinear_matrix(answer);
    std::cout << "Result of F(Answer): " << newton_solve_check;
    std::cout << "\n\n------------------------------------\n\n";
}


#endif // ENABLE_TESTS_NEWTON_SYSTEM

#if ENABLE_TESTS_MOD_NEWTON_SYSTEM

TEST_METHOD(mod_newton_test_solve_system)
{
    std::cout << "Modified Newton method for nonlinear system:\n\n";

    for (long k = 1; k <= kMax; ++k)
    {
        auto start = std::chrono::steady_clock::now();
        const auto mod_newton_solve_system = mod_newton_method_system(nonlinear_matrix,
                                                                      jacobi_matrix, vector_x,
                                                                      k, true, kEps);
        std::chrono::duration<double> duration = std::chrono::steady_clock::now() - start;
        std::cout << k << " Calculation time: " << duration.count() << " ms\n";

        const auto answer = std::get<0>(mod_newton_solve_system);
        const auto iterations = std::get<1>(mod_newton_solve_system);
        const auto complexity = std::get<2>(mod_newton_solve_system);
        //std::cout << "Epsilon = " << kEps << '\n';
        //std::cout << "x0:\n" << vector_x << "\n\n";
        std::cout << "Answer:\n" << answer << "\n\n";
        std::cout << "Number of iterations = " << iterations << '\n';
        std::cout << "Number of operations = " << complexity;
        std::cout << "\n\n------------------------------------\n\n";

        std::cout << "Checking modified Newton solve (include eps = " << kEps << "):\n\n";
        const auto newton_solve_check = nonlinear_matrix(answer);
        std::cout << "Result of F(Answer): " << newton_solve_check;
        std::cout << "\n\n------------------------------------\n\n";
    }
}


#endif // ENABLE_TESTS_MOD_NEWTON_SYSTEM

#if ENABLE_TESTS_HYBRID_NEWTON_SYSTEM

TEST_METHOD(hybrid_newton_test_solve_system)
{
    std::cout << "Hybrid Newton method for nonlinear system:\n\n";

    for (long k = 1; k <= kMax; ++k)
    {
        if (k == 2) continue;
        auto start = std::chrono::steady_clock::now();
        const auto hybrid_newton_solve_system = hybrid_newton_method_system(nonlinear_matrix,
                                                                            jacobi_matrix,
                                                                            vector_x,
                                                                            k, kEps);
        std::chrono::duration<double> duration = std::chrono::steady_clock::now() - start;
        std::cout << k << " Calculation time: " << duration.count() << " ms\n";

        const auto answer = std::get<0>(hybrid_newton_solve_system);
        const auto iterations = std::get<1>(hybrid_newton_solve_system);
        const auto complexity = std::get<2>(hybrid_newton_solve_system);
        //std::cout << "Epsilon = " << kEps << '\n';
        //std::cout << "x0:\n" << vector_x << "\n\n";
        std::cout << "Answer:\n" << answer << "\n\n";
        std::cout << "Number of iterations = " << iterations << '\n';
        std::cout << "Number of operations = " << complexity;
        std::cout << "\n\n------------------------------------\n\n";

        //std::cout << "Checking hybrid Newton solve (include eps = " << kEps << "):\n\n";
        //const auto newton_solve_check = nonlinear_matrix(answer);
        //std::cout << "Result of F(Answer): " << newton_solve_check;
        //std::cout << "\n\n------------------------------------\n\n";
    }
}


#endif // ENABLE_TESTS_HYBRID_NEWTON_SYSTEM

} //namesace vv
