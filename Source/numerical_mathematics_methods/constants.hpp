#pragma once

#include "cx_matrix.hpp"


namespace vv
{

/// ===== CONSTANT SECTION =====
// ===== DATA 0 (base constants) =====
template <class Type>
constexpr Type kDefault_eps = static_cast<Type>(1e-6);

template <class Type>
constexpr Type kRough_eps = static_cast<Type>(1e-3);

constexpr int N = 12;
constexpr double kEps = kDefault_eps<double>;

// ===== DATA 1 (base tests) =====
constexpr cx_matrix<double, 3, 3> data_N{ N + 2, 1,     1,
                                           1,     N + 4, 1,
                                           1,     1,     N + 6 };

constexpr cx_matrix<double, 3, 1> vec_b_N{ N + 4, N + 6, N + 8 };

// ===== DATA 2 (bad case) =====
constexpr cx_matrix<double, 3, 3> data{ 1, -1, -1,
                                        0,  1, -1,
                                        0,  0,  1  };

constexpr cx_matrix<double, 3, 3> eps_mat{ kEps * N, -kEps * N, -kEps * N,
                                           kEps * N,  kEps * N, -kEps * N,
                                           kEps * N,  kEps * N,  kEps * N };

constexpr cx_matrix<double, 3, 3> data_eps{ 1 + kEps * N, -1 - kEps * N,  -1 - kEps * N,
                                            kEps * N,      1 + kEps * N, -1 - kEps * N,
                                            kEps * N,      kEps * N,      1 + kEps * N };

constexpr cx_matrix<double, 3, 1> vec_b_eps{ -1, -1, -1 };

// ===== DATA 3 (singular cx_matrix, no solutions) =====
constexpr cx_matrix<double, 4, 4> data_singular{  2,  -1, 4, 7,
                                                   4,  -2, 8, 5,
                                                  -2,   1, 6, 9,
                                                  10,  -5, 1, 3  };

constexpr cx_matrix<double, 4, 1> vec_b_singular{ 4, 8, 6, 1 };

// ===== DATA 4 (for rank test) =====
constexpr cx_matrix<double, 4, 4> data_rank{  2.0,  1.0, 11.0,  2.0,
                                               1.0,  0.0,  4.0, -1.0,
                                              11.0,  4.0, 56.0,  5.0,
                                               2.0, -1.0,  5.0, -6.0  };

// ===== DATA 5 (not square cx_matrix, no solutions) =====
constexpr cx_matrix<double, 5, 3> data_not_square_singular{   2.0,  1.0, 11.0,
                                                               2.0,  1.0,  1.0,
                                                               0.0,  4.0, -1.0,
                                                               3.0, 11.0,  4.0,
                                                              56.0,  5.0,  5.0  };

constexpr cx_matrix<double, 5, 1> vec_b_not_square_singular{ 1.0, 0.0, 2.0, 3.0, 7.0 };

// ===== DATA 6 (not square cx_matrix) =====
constexpr cx_matrix<double, 5, 3> data_not_square_row{  3.0,  1.0,  1.0,
                                                         2.0,  2.0,  1.0,
                                                         7.0,  5.0,  2.0,
                                                         1.0, 11.0,  7.0,
                                                        19.0, 17.0, 16.0  };

constexpr cx_matrix<double, 5, 1> vec_b_not_square_row{ 5.0, 5.0, 14.0, 19.0, 52.0 };

// ===== DATA 7 (not square cx_matrix) =====
constexpr cx_matrix<double, 3, 5> data_not_square_col{  3.0, 1.0, 1.0, 4.0, 3.0,
                                                         2.0, 2.0, 1.0, 5.0, 4.0,
                                                         7.0, 5.0, 2.0, 1.0, 9.0  };

constexpr cx_matrix<double, 3, 1> vec_b_not_square_col{ 12.0, 14.0, 24.0 };

#define mat_A data_eps
#define vec_b vec_b_eps

void show_constants()
{
    std::cout << "Matrix A:\n" << vv::mat_A << "\n\n";
    std::cout << "Matrix b:\n" << vv::vec_b;
    std::cout << "\n\n------------------------------------\n\n";
}

} // namespace vv