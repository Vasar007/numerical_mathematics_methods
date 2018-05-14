#pragma once

#include <iostream>

#include "cx_matrix.hpp"
#include "constants.hpp"


namespace vv
{

/// ===== TEST SECTION =====
#if ENABLE_TESTS_GENERAL

TEST_METHOD(general_test_condition_number)
{
    std::cout << "Calculated condition number:\n\n";
    
    constexpr auto condition_number = mat_A.calculate_condition_number();
    std::cout << "Matrix A:\n" << mat_A << "\n\n";
    std::cout << "Condition number of matrix A:\n" << condition_number;
    std::cout << "\n\n------------------------------------\n\n";
}

#endif // ENABLE_TESTS_GENERAL

#if ENABLE_TESTS_GAUSSIAN

TEST_METHOD(gaussian_test_for_matrix_N)
{
    std::cout << "Matrix solved with Gaussian:\n\n";
    
    constexpr auto vec_x = cx_matrix<double>::solve(mat_A, vec_b);
    std::cout << "Matrix x:\n" << vec_x;
    std::cout << "\n\n------------------------------------\n\n";

    std::cout << "Checking Gaussian solve:\n\n";
    constexpr auto gaussian_check = (mat_A * vec_x) - vec_b;
    std::cout << "Result of (A * x - b):\n" << gaussian_check;
    std::cout << "\n\n------------------------------------\n\n";
}


TEST_METHOD(gaussian_test_for_matrix_eps_N)
{
    std::cout << "Matrix solved with Gaussian (include kEps = " << kEps <<  ")):\n\n";
    
    constexpr auto vec_x = cx_matrix<double>::solve(mat_A, vec_b, kEps);
    std::cout << "Matrix x:\n" << vec_x;
    std::cout << "\n\n------------------------------------\n\n";

    std::cout << "Checking Gaussian solve (include kEps = " << kEps <<  ")):\n\n";
    constexpr auto gaussian_solve_check = (mat_A * vec_x) - vec_b;
    std::cout << "Result of (Aeps * x - b):\n" << gaussian_solve_check;
    std::cout << "\n\n------------------------------------\n\n";
}


TEST_METHOD(gaussian_test_inverse)
{
    std::cout << "Invert by matrix method (Gaussian):\n\n";
    
    constexpr auto casual_invert = mat_A.inverse();
    std::cout << "Matrix A:\n" << mat_A << "\n\n";
    std::cout << "Inverted matrix to A:\n" << casual_invert;
    std::cout << "\n\n------------------------------------\n\n";

    std::cout << "Checking inverse:\n\n";
    constexpr auto inverse_check1 = mat_A * casual_invert;
    std::cout << "Result of (A * IA):\n" << inverse_check1 << "\n\n";

    constexpr auto inverse_check2 = casual_invert * mat_A;    
    std::cout << "Result of (IA * A):\n" << inverse_check2;
    std::cout << "\n\n------------------------------------\n\n";
}

#endif // ENABLE_TESTS_GAUSSIAN

} // namespace vv