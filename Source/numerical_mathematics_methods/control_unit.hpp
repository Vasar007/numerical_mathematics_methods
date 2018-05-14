#pragma once

#include <iostream>
#include <stdexcept>


namespace vv
{

/// ===== MACROS SECTION =====
#define TRY_BLOCK(expr) try { expr } \
                        catch(std::exception& ex) { std::cout << ex.what() << "\n\n"; } \
                        catch(const char* ex) { std::cout << ex << "\n\n"; } \
                        catch(const std::string& ex) { std::cout << ex << "\n\n"; } \
                        catch(...) { std::cout << "Not handled exception!" << "\n\n"; }                        

#define TEST_METHOD_RUNTIME(name) template <class Type, std::size_t Rows, std::size_t Columns> \
                        void name(const matrix<Type, Rows, Columns>& mat_A, \
                                  const matrix<Type, Rows, 1>& vec_b = matrix<Type, Rows, 1>{}, \
                                  const Type eps = kDefault_eps<Type>)

#define TEST_METHOD(name) void name()


TEST_METHOD(general_test_condition_number);

TEST_METHOD(gaussian_test_for_matrix_N);
TEST_METHOD(gaussian_test_for_matrix_eps_N);
TEST_METHOD(gaussian_test_inverse);

TEST_METHOD(functions_test);

TEST_METHOD(lu_test_decompose);
TEST_METHOD(lu_test_determenant);
TEST_METHOD(lu_test_rank);

TEST_METHOD(lup_test_decompose);
TEST_METHOD(lup_test_solve);

TEST_METHOD(lup_test_determenant);
TEST_METHOD(lup_test_inverse);

TEST_METHOD(lupq_test_decompose);
TEST_METHOD(lupq_test_solve);

TEST_METHOD(qr_test_decomposition);
TEST_METHOD(qr_test_solve);

TEST_METHOD(jacobi_test_solve);
TEST_METHOD(seidel_test_solve);

TEST_METHOD(gen_test_diagonal_predominance);
TEST_METHOD(gen_test_simmetrical);

TEST_METHOD(newton_test_solve);
TEST_METHOD(newton_test_solve_system);
TEST_METHOD(mod_newton_test_solve_system);
TEST_METHOD(hybrid_newton_test_solve_system);

TEST_METHOD(newton_cotes_test_calculation_integral);
TEST_METHOD(gauss_test_calculation_integral);

TEST_METHOD(runge_kutta_calculation);


#define ENABLE_TESTS_GENERAL               0
#define ENABLE_TESTS_GAUSSIAN              0
#define ENABLE_TESTS_FUNCTIONS             1
#define ENABLE_TESTS_LU                    0
#define ENABLE_TESTS_LUP                   0
#define ENABLE_TESTS_LUP_ONLY_SQUARE       0
#define ENABLE_TESTS_LUPQ                  0
#define ENABLE_TESTS_QR                    0
#define ENABLE_TESTS_QR_SOLVE              0
#define ENABLE_TESTS_ITERATIONS            0
#define ENABLE_TESTS_ITERATIONS_DIAGONAL   0
#define ENABLE_TESTS_ITERATONS_SIMMETRICAL 0
#define ENABLE_TESTS_NEWTON                0
#define ENABLE_TESTS_NEWTON_SYSTEM         0
#define ENABLE_TESTS_MOD_NEWTON_SYSTEM     0
#define ENABLE_TESTS_HYBRID_NEWTON_SYSTEM  0
#define ENABLE_TESTS_NEWTON_COTES_INTEGRAL 0
#define ENABLE_TESTS_GAUSS_INTEGRAL        0
#define ENABLE_TESTS_RUNGE_KUTTA           0


void start_tests()
{
#if ENABLE_TESTS_GENERAL
    general_test_condition_number();
#endif // ENABLE_TESTS_GENERAL


#if ENABLE_TESTS_GAUSSIAN
    gaussian_test_for_matrix_N();
    gaussian_test_for_matrix_eps_N();
    gaussian_test_inverse();
#endif // ENABLE_TESTS_GAUSSIAN


#if ENABLE_TESTS_FUNCTIONS
    functions_test();
#endif // ENABLE_TESTS_FUNCTIONS


#if ENABLE_TESTS_LU
    lu_test_decompose();
    lu_test_determenant();
    lu_test_rank();
#endif // ENABLE_TESTS_LU


#if ENABLE_TESTS_LUP
    lup_test_decompose();
    lup_test_solve();
#endif // ENABLE_TESTS_LUP

#if ENABLE_TESTS_LUP_ONLY_SQUARE
    lup_test_determenant();
    lup_test_inverse();
#endif // ENABLE_TESTS_LUP_ONLY_SQUARE


#if ENABLE_TESTS_LUPQ
    lupq_test_decompose();
    lupq_test_solve();
#endif // ENABLE_TESTS_LUPQ


#if ENABLE_TESTS_QR
    qr_test_decomposition();
#endif // ENABLE_TESTS_QR


#if ENABLE_TESTS_QR_SOLVE
    qr_test_solve();
#endif // ENABLE_TESTS_QR_SOLVE


#if ENABLE_TESTS_ITERATIONS
    jacobi_test_solve();
    seidel_test_solve();
#endif // ENABLE_TESTS_ITERATIONS


#if ENABLE_TESTS_ITERATIONS_DIAGONAL
    gen_test_diagonal_predominance();
#endif // ENABLE_TESTS_ITERATIONS_DIAGONAL


#if ENABLE_TESTS_ITERATONS_SIMMETRICAL
    gen_test_simmetrical();
#endif // ENABLE_TESTS_ITERATONS_SIMMETRICAL


#if ENABLE_TESTS_NEWTON
    newton_test_solve();
#endif // ENABLE_TESTS_NEWTON


#if ENABLE_TESTS_NEWTON_SYSTEM
    newton_test_solve_system();
#endif // ENABLE_TESTS_NEWTON_SYSTEM


#if ENABLE_TESTS_MOD_NEWTON_SYSTEM
    mod_newton_test_solve_system();
#endif // ENABLE_TESTS_MOD_NEWTON_SYSTEM


#if ENABLE_TESTS_HYBRID_NEWTON_SYSTEM
    hybrid_newton_test_solve_system();
#endif // ENABLE_TESTS_HYBRID_NEWTON_SYSTEM

#if ENABLE_TESTS_NEWTON_COTES_INTEGRAL
    newton_cotes_test_calculation_integral();
#endif // ENABLE_TESTS_NEWTON_COTES_INTEGRAL

#if ENABLE_TESTS_GAUSS_INTEGRAL
    gauss_test_calculation_integral();
#endif // ENABLE_TESTS_GAUSS_INTEGRAL

#if ENABLE_TESTS_RUNGE_KUTTA
    runge_kutta_calculation();
#endif // ENABLE_TESTS_RUNGE_KUTTA
}

} // namespace vv

#include "general_tests.hpp"
#include "cx_functions.hpp"
#include "lu_decomposition.hpp"
#include "lup_decomposition.hpp"
#include "lupq_decomposition.hpp"
#include "qr_decomposition.hpp"
#include "iterations_methods.hpp"
#include "newton_method.hpp"
#include "newton_cotes.hpp"
#include "gauss_integral.hpp"
#include "runge_kutta_methods.hpp"
