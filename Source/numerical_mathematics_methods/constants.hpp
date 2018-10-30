// Copyright (C) 2018 Vasily Vasilyev (vasar007@yandex.ru)

#pragma once

#include "cx_matrix.hpp"


namespace vv
{

/// ===== CONSTANT SECTION =====
// ===== DATA 0 (base constants) =====
template <typename Type>
constexpr Type kDefault_eps = static_cast<Type>(1e-6);

template <typename Type>
constexpr Type kRough_eps = static_cast<Type>(1e-3);

constexpr int N = 12;
constexpr double kEps = kDefault_eps<double>;

// ===== DATA 1 (base tests) =====
constexpr cx_matrix<double, 3, 3> data_N{ N + 2.0,     1.0,     1.0,
                                              1.0, N + 4.0,     1.0,
                                              1.0,     1.0, N + 6.0  };

constexpr cx_matrix<double, 3, 1> vec_b_N{ N + 4.0, N + 6.0, N + 8.0 };

// ===== DATA 2 (bad case) =====
constexpr cx_matrix<double, 3, 3> data{ 1.0, -1.0, -1.0,
                                        0.0,  1.0, -1.0,
                                        0.0,  0.0,  1.0  };

constexpr cx_matrix<double, 3, 3> eps_mat{ kEps * N, -kEps * N, -kEps * N,
                                           kEps * N,  kEps * N, -kEps * N,
                                           kEps * N,  kEps * N,  kEps * N };

constexpr cx_matrix<double, 3, 3> data_eps{ 1.0 + kEps * N, -1.0 - kEps * N, -1.0 - kEps * N,
                                                  kEps * N,  1.0 + kEps * N, -1.0 - kEps * N,
                                                  kEps * N,        kEps * N,  1.0 + kEps * N };

constexpr cx_matrix<double, 3, 1> vec_b_eps{ -1, -1, -1 };

// ===== DATA 3 (singular cx_matrix, no solutions) =====
constexpr cx_matrix<double, 4, 4> data_singular{  2.0, -1.0, 4.0, 7.0,
                                                  4.0, -2.0, 8.0, 5.0,
                                                 -2.0,  1.0, 6.0, 9.0,
                                                 10.0, -5.0, 1.0, 3.0  };

constexpr cx_matrix<double, 4, 1> vec_b_singular{ 4.0, 8.0, 6.0, 1.0 };

// ===== DATA 4 (for rank test) =====
constexpr cx_matrix<double, 4, 4> data_rank{ 2.0,   1.0, 11.0,  2.0,
                                              1.0,  0.0,  4.0, -1.0,
                                             11.0,  4.0, 56.0,  5.0,
                                              2.0, -1.0,  5.0, -6.0  };

// ===== DATA 5 (not square cx_matrix, no solutions) =====
constexpr cx_matrix<double, 5, 3> data_not_square_singular{  2.0,  1.0,  11.0,
                                                             2.0,  1.0,   1.0,
                                                             0.0,  4.0,  -1.0,
                                                             3.0, 11.0,   4.0,
                                                            56.0,  5.0,   5.0  };

constexpr cx_matrix<double, 5, 1> vec_b_not_square_singular{ 1.0, 0.0, 2.0, 3.0, 7.0 };

// ===== DATA 6 (not square cx_matrix) =====
constexpr cx_matrix<double, 5, 3> data_not_square_row{  3.0,  1.0, 1.0,
                                                        2.0,  2.0, 1.0,
                                                        7.0,  5.0, 2.0,
                                                        1.0, 11.0, 7.0,
                                                       19.0, 17.0, 16.0 };

constexpr cx_matrix<double, 5, 1> vec_b_not_square_row{ 5.0, 5.0, 14.0, 19.0, 52.0 };

// ===== DATA 7 (not square cx_matrix) =====
constexpr cx_matrix<double, 3, 5> data_not_square_col{ 3.0, 1.0, 1.0, 4.0, 3.0,
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