#pragma once

#include <iterator>
#include <algorithm>
#include <array>
#include <stdexcept>
#include <initializer_list>
#include <string>
#include <cmath>
#include <cassert>
#include <iterator>

#include "cx_math.h"
#include "cx_algorthms.hpp"


/*
 *
 * As base for cx_matrix class using https://github.com/akalicki/matrix
 *
 */

namespace vv
{
/// ===== CONSTEXPR MATRIX IMLEMENTATION SECTION =====
template <class Type = double, std::size_t Rows = 1, std::size_t Columns = 1>
class cx_matrix
{
public:
    using value_type                    = Type;
    using size_type                     = std::size_t;
    using difference_type               = std::ptrdiff_t;

    template <class ConType = value_type, std::size_t N = Rows>
    using container                     = std::array<ConType, N>;
    using row_container                 = container<value_type, Columns>;
    using row_container_reference       = container<value_type, Columns>&;    
    using const_row_container_reference = const container<value_type, Columns>&;

    using pointer                       = value_type*;
    using const_pointer                 = const value_type*;

    using reference                     = value_type&;
    using const_reference               = const value_type&;


    static constexpr value_type EPS = static_cast<value_type>(1e-10);

    static_assert(std::is_arithmetic_v<value_type>, "Matrix elements type has to be arithmetic!");
    static_assert(Rows > 0 && Columns > 0, "Incorrect size parameters!");

    constexpr cx_matrix() noexcept = default;

    constexpr cx_matrix(const value_type value) noexcept
    : _data()
    {
        detail::fill(std::begin(_data), std::end(_data), value);
    }

    constexpr cx_matrix(const std::initializer_list<value_type> list) noexcept
    : _data()
    {
        size_type row_counter = 0;
        size_type col_counter = 0;
        for (const auto elem : list)
        {
            _data.at(row_counter).at(col_counter) = elem;
            ++col_counter;

            if (col_counter == Columns)
            {
                ++row_counter;
                if (row_counter == Rows && col_counter == Columns)
                {
                    break;
                }
                col_counter = 0;
            }
        }
    }
    

    ~cx_matrix() noexcept = default;
    constexpr cx_matrix(const cx_matrix& other) = default;
    constexpr cx_matrix& operator=(const cx_matrix& other) = default;
    constexpr cx_matrix(cx_matrix&& other) noexcept = default;
    constexpr cx_matrix& operator=(cx_matrix&& other) noexcept = default;



    std::string get_dimension() const
    {
        return std::to_string(get_rows_number()) + std::string("x")
               + std::to_string(get_columns_number());
    }


    constexpr size_type get_rows_number() const noexcept
    {
        return Rows;
    }


    constexpr size_type get_columns_number() const noexcept
    {
        return Columns;
    }

    constexpr bool empty() const noexcept
    {
        // Without any checks because class has static_assert and code won't be compiled
        // if number of Rows or Colums is equal to 0.
        return false;
    }


    constexpr const_reference operator()(const size_type i, const size_type j) const
    {
        return _data.at(i).at(j);
    }


    constexpr reference operator()(const size_type i, const size_type j)
    { 
        return _data.at(i).at(j);
    }


    constexpr const_row_container_reference& operator()(const size_type i) const
    {
        return _data.at(i);
    }


    constexpr row_container_reference& operator()(const size_type i)
    {
        return _data.at(i);
    }


    constexpr cx_matrix& operator+=(const cx_matrix& rhs) noexcept
    {
        for (size_type i = 0; i < Rows; ++i)
        {
            for (size_type j = 0; j < Columns; ++j)
            {
                _data.at(i).at(j) += rhs(i, j);
            }
        }
        return *this;
    }
    

    constexpr cx_matrix& operator-=(const cx_matrix& rhs) noexcept
    {
        for (size_type i = 0; i < Rows; ++i)
        {
            for (size_type j = 0; j < Columns; ++j)
            {
                _data.at(i).at(j) -= rhs(i, j);
            }
        }
        return *this;
    }
    

    constexpr cx_matrix& operator*=(const value_type value) noexcept
    {
        for (auto& row : _data)
        {
            for (auto& elem : row)
            {
                elem *= value;
            }
        }
        return *this;
    }
    

    constexpr cx_matrix& operator/=(const value_type value)
    {
        assert(value != value_type{});

        for (auto& row : _data)
        {
            for (auto& elem : row)
            {
                elem /= value;
            }
        }
        return *this;
    }


    constexpr auto operator^(const int value) const noexcept
    {
        cx_matrix temp(*this);
        return _exp_helper(temp, value);
    }


    constexpr void swap_rows(const size_type row_1, const size_type row_2) noexcept
    {
        detail::swap(_data.at(row_1), _data.at(row_2));
    }


    constexpr void swap_columns(const size_type column_1, const size_type column_2) noexcept
    {
        for (size_type i = 0; i < Rows; ++i)
        {
            detail::swap(_data.at(i).at(column_1), _data.at(i).at(column_2));
        }
    }
    

    constexpr cx_matrix<value_type, Columns, Rows> transpose() const noexcept
    {
        cx_matrix<value_type, Columns, Rows> transp{};
        for (size_type i = 0; i < Rows; ++i)
        {
            for (size_type j = 0; j < Columns; ++j)
            {
                transp(j, i) = _data.at(i).at(j);
            }
        }
        return transp;
    }


    constexpr bool is_empty(const size_type start_row_index = 0,
                            const size_type start_column_index = 0,
                            const value_type eps = EPS) const noexcept
    {
        if (start_row_index > Rows || start_column_index > Columns)
        {
            return false;
        }

        for (size_type i = start_row_index; i < Rows; ++i)
        {
            for (size_type j = start_column_index; j < Columns; ++j)
            {
                if (cx::abs(_data.at(i).at(j)) >= eps)
                {
                    return false;
                }
            }
        }

        return true;        
    }


    bool is_error_matrix() const noexcept
    {
        //const auto is_nan = [](const value_type x) { return std::isnan(x); };
        for (const auto& row : _data)
        {
            if (!detail::all_of(std::begin(row), std::end(row), std::isnan<value_type>))
            {
                return false;
            }
        }
        return true;
    }


    constexpr std::pair<size_type, size_type>
        find_max_element(const size_type start_row_index = 0,
                         const size_type start_column_index = 0) const noexcept
    {
        if (start_row_index > Rows || start_column_index > Columns)
        {
            return { start_row_index, start_column_index };
        }

        value_type abs_max{};
        size_type row = start_row_index;
        size_type column = start_column_index;

        for (size_type i = start_row_index; i < Rows; ++i)
        {
            for (size_type j = start_column_index; j < Columns; ++j)
            {
                if (cx::abs(_data.at(i).at(j)) > abs_max)
                {
                    abs_max = cx::abs(_data.at(i).at(j));
                    row = i;
                    column = j;
                }
            }
        }
        return { row, column };
    }


    constexpr value_type calculate_condition_number() const noexcept
    {
        constexpr auto abs_plus = [](const value_type a, const value_type b)
        {
            return cx::abs(a) + cx::abs(b);
        };

        value_type condition_number{};
        for (const auto& row: _data)
        {
            const auto summ_in_row = detail::accumulate(std::begin(row), std::end(row),
                                                        value_type{}, abs_plus);
            condition_number = std::max(condition_number, summ_in_row);
        }
        return condition_number;
    }


    template <class T = value_type, std::size_t Size = Rows>
    static constexpr cx_matrix<T, Size, Size> create_identity() noexcept
    {
        cx_matrix<T, Size, Size> temp{};
        for (size_type i = 0; i < Size; ++i)
        {
            for (size_type j = 0; j < Size; ++j)
            {
                if (i == j)
                {
                    temp(i, j) = static_cast<T>(1);
                }
                else
                {
                    temp(i, j) = value_type{};
                }
            }
        }
        return temp;
    }


    template <std::size_t Rows_, std::size_t Columns_A, std::size_t Columns_b>
    static constexpr std::pair<cx_matrix<value_type, Rows_, Columns_A>,
                               cx_matrix<value_type, Rows_, Columns_b>>
        forward_substitution(cx_matrix<value_type, Rows_, Columns_A> A,
                             cx_matrix<value_type, Rows_, Columns_b> b,
                             const value_type eps = EPS)
    {
        static_assert(Columns_b == 1, "Matrix contains more than one columns in right "
                                       "hand side vector!");

        // Gaussian elimination.
        for (size_type i = 0; i < Rows_; ++i)
        {
            if (cx::abs(A(i, i)) < eps)
            {
                // Pivot ~ 0 => throw error.
                throw std::domain_error("Error: the coefficient cx_matrix has 0 as a pivot.");
            }

            for (size_type j = i + 1; j < Rows_; ++j)
            {
                for (size_type k = i + 1; k < Columns_A; ++k)
                {
                    A(j, k) -= A(i, k) * (A(j, i) / A(i, i));
                    if (cx::abs(A(j, k)) < eps)
                    {
                        A(j, k) = value_type{};
                    }
                }

                b(j, 0) -= b(i, 0) * (A(j, i) / A(i, i));
                if (cx::abs(A(j, 0)) < eps)
                {
                    A(j, 0) = value_type{};
                }
                A(j, i) = value_type{};
            }
        } // for (size_type i = 0; i < Rows_; ++i)

        return { A, b };
    }


    template <std::size_t Rows_, std::size_t Columns_A, std::size_t Columns_b>
    static constexpr cx_matrix<value_type, Rows_, Columns_b>
        backward_substitution(const cx_matrix<value_type, Rows_, Columns_A>& A,
                              const cx_matrix<value_type, Rows_, Columns_b>& b,
                              const value_type eps = EPS)
    {
        static_assert(Columns_b == 1, "Matrix contains more than one columns in right "
                                       "hand side vector!");

        // Back substitution.
        cx_matrix<value_type, Rows_, 1> x{};
        for (size_type i = Rows_ - 1; ; --i)
        {
            value_type sum{};
            for (size_type j = i + 1; j < Rows_; ++j)
            {
                sum += A(i, j) * x(j, 0);
            }

            x(i, 0) = (b(i, 0) - sum) / A(i, i);
            if (cx::abs(x(i, 0)) < eps)
            {
                x(i, 0) = value_type{};
            }

            if (i == 0)
            {
                break;
            }
        } // for (size_type i = Rows_ - 1; ; --i)
        
        return x;
    }

    template <std::size_t Rows_, std::size_t Columns_A, std::size_t Columns_b>
    static constexpr cx_matrix<value_type, Rows_, Columns_b>
        solve(const cx_matrix<value_type, Rows_, Columns_A>& A,
              const cx_matrix<value_type, Rows_, Columns_b>& b, const value_type eps = EPS)
    {
        static_assert(Columns_b == 1, "Matrix contains more than one columns in right "
                                      "hand side vector!");

        const auto result_of_forw_subs = forward_substitution(A, b, eps);

        return backward_substitution(result_of_forw_subs.first, result_of_forw_subs.second);
    }


    template <std::size_t Rows_A, std::size_t Columns_A,
              std::size_t Rows_b, std::size_t Columns_b>
    static constexpr cx_matrix<value_type, Rows_b, 1>
        band_solve(cx_matrix<value_type, Rows_A, Columns_A> A,
                   cx_matrix<value_type, Rows_b, Columns_b> b,
                   const value_type coef)
    {
        // Optimized Gaussian elimination.
        value_type bandsBelow = static_cast<value_type>((coef - 1) / 2.0);
        for (size_type i = 0; i < Rows_A; ++i)
        {
            if (A(i, i) == value_type{})
            {
                // Pivot ~0 => throw exception.
                throw std::domain_error("Error: the coefficient cx_matrix has 0 as a "
                                        "pivot. Please fix the input and try again.");
            }
            for (size_type j = i + 1; j < Rows_A && j <= i + bandsBelow; ++j)
            {
                size_type k = i + 1;
                while (k < Columns_A && A(j, k))
                {
                    A(j, k) -= A(i, k) * (A(j, i) / A(i, i));
                    ++k;
                }
                b(j, 0u) -= b(i, 0u) * (A(j, i) / A(i, i));
                A(j, i) = value_type{};
            }
        }

        // Back substitution.
        cx_matrix<value_type, Rows_b, 1> x{};
        x(Rows_b - 1, 0) = b(Rows_b - 1, 0) / A(Rows_b - 1, Rows_b - 1);
        for (size_type i = Rows_b - 2; ; --i)
        {
            value_type sum{};
            for (size_type j = i + 1; j < Rows_b; ++j)
            {
                sum += A(i, j) * x(j, 0);
            }
            x(i, 0u) = (b(i, 0) - sum) / A(i, i);

            if (i == 0)
            {
                break;
            }
        }

        return x;
    }


    // Functions on vectors.
    template <std::size_t Rows_, std::size_t Columns_A, std::size_t Columns_B>
    static constexpr value_type dot_product(const cx_matrix<value_type, Rows_, Columns_A>& A,
                                            const cx_matrix<value_type, Rows_, Columns_B>& B,
                                            const size_type column) noexcept
    {
        value_type sum{};
        for (size_type i = 0; i < Rows_; ++i)
        {
            sum += (A(i, column) * B(i, column));
        }
        return sum;
    }


    // Functions on augmented matrices.
    template <std::size_t Rows_A, std::size_t Columns_A,
              std::size_t Rows_B, std::size_t Columns_B>
    static constexpr cx_matrix<value_type, Rows_A, Columns_A + Columns_B>
        augment(const cx_matrix<value_type, Rows_A, Columns_A>& A,
                const cx_matrix<value_type, Rows_B, Columns_B>& B) noexcept
    {
        cx_matrix<value_type, Rows_A, Columns_A + Columns_B> AB{};
        for (size_type i = 0; i < Rows_A; ++i)
        {
            for (size_type j = 0; j < Columns_A + Columns_B; ++j)
            {
                if (j < Columns_A)
                {
                    AB(i, j) = A(i, j);
                }
                else
                {
                    AB(i, j) = B(i, j - Columns_B);
                }
            }
        }
        return AB;
    }

    constexpr cx_matrix gaussian_eliminate(const value_type eps = EPS) const
    {
        cx_matrix Ab(*this);
        int rows = Rows;
        int cols = Columns;
        int Acols = cols - 1;

        int i = 0; // Row tracker.
        int j = 0; // Column tracker.

        // Iterate through the rows.
        while (i < rows)
        {
            // Find a pivot for the row.
            bool pivot_found = false;
            while (j < Acols && !pivot_found)
            {
                if (Ab(i, j) != value_type{})
                { // Pivot not equal to 0.
                    pivot_found = true;
                }
                else
                { // Check for a possible swap.
                    int max_row = i;
                    value_type max_val{};
                    for (int k = i + 1; k < rows; ++k)
                    {
                        value_type cur_abs = cx::abs(Ab(k, j));
                        if (cur_abs > max_val)
                        {
                            max_row = k;
                            max_val = cur_abs;
                        }
                    }
                    if (max_row != i)
                    {
                        Ab.swap_rows(max_row, i);
                        pivot_found = true;
                    }
                    else
                    {
                        ++j;
                    }
                }
            }

            // Perform elimination as normal if pivot was found.
            if (pivot_found)
            {
                for (int t = i + 1; t < rows; ++t)
                {
                    for (int s = j + 1; s < cols; ++s)
                    {
                        Ab(t, s) = Ab(t, s) - Ab(i, s) * (Ab(t, j) / Ab(i, j));
                        if (Ab(t, s) < eps && Ab(t, s) > -1 * eps)
                        {
                            Ab(t, s) = value_type{};
                        }
                    }
                    Ab(t, j) = value_type{};
                }
            }

            ++i;
            ++j;
        }

        return Ab;
    }

    constexpr cx_matrix row_reduce_from_gaussian(const value_type eps = EPS) const
    {
        cx_matrix result(*this);
        int rows = Rows;
        int cols = Columns;

        int i = rows - 1; // Row tracker.
        int j = cols - 1; // Column tracker.

        // Iterate through every row.
        while (i >= 0 && j >= 0)
        {
            // Find the pivot column.
            int k = j - 1;
            while (k >= 0)
            {
                if (result(i, k) != value_type{})
                {
                    j = k;
                }
                --k;
            }

            // Zero out elements above pivots if pivot not 0.
            if (result(i, j) != value_type{})
            {
            
                for (int t = i - 1; t >= 0; --t)
                {
                    for (int s = 0; s < cols; ++s)
                    {
                        if (s != j)
                        {
                            result(t, s) = result(t, s) - result(i, s) 
                                            * (result(t, j) / result(i, j));
                            if (result(t, s) < eps && result(t, s) > -1 * eps)
                            {
                                result(t, s) = value_type{};
                            }
                        }
                    }
                    result(t, j) = value_type{};
                }

                // Divide row by pivot.
                for (int l = j + 1; l < cols; ++l)
                {
                    result(i, l) = result(i, l) / result(i, j);
                    if (result(i, l) < eps && result(i, l) > -1 * eps)
                    {
                        result(i, l) = value_type{};
                    }
                }
                result(i, j) = static_cast<value_type>(1);

            }

            --i;
            --j;
        }

        return result;
    }


    // Not tested this method.
    void read_solutions_from_RREF(std::ostream& os)
    {
        cx_matrix result(*this);

        // Print number of solutions.
        bool has_solutions = true;
        bool dones_searching = false;
        size_type i = 0;
        while (!dones_searching && i < Rows)
        {
            bool allZeros = true;
            for (size_type j = 0; j < Columns - 1; ++j)
            {
                if (result(i, j) != value_type{})
                {
                    allZeros = false;
                }
            }
            if (allZeros && result(i, Columns - 1) != value_type{})
            {
                has_solutions = false;
                os << "NO SOLUTIONS\n\n";
                dones_searching = true;
            }
            else if (allZeros && result(i, Columns - 1) == value_type{})
            {
                os << "INFINITE SOLUTIONS\n\n";
                dones_searching = true;
            }
            else if (Rows < Columns - 1)
            {
                os << "INFINITE SOLUTIONS\n\n";
                dones_searching = true;
            }
            ++i;
        }
        if (!dones_searching)
        {
            os << "UNIQUE SOLUTION\n\n";
        }

        // Get solutions if they exist.
        if (has_solutions)
        {
            cx_matrix<value_type, Columns - 1, 1> particular{};
            cx_matrix<value_type, Columns - 1, 1> special{};

            for (size_type j = 0; j < Rows; ++j)
            {
                bool pivotFound = false;
                bool specialCreated = false;
                for (size_type k = 0; k < Columns - 1; ++k)
                {
                    if (result(k, k) != value_type{})
                    {
                        // If pivot variable, add b to particular.
                        if (!pivotFound)
                        {
                            pivotFound = true;
                            particular(k, 0) = result(k, Columns - 1);
                        }
                        else
                        { // Otherwise, add to special solution.
                            if (!specialCreated)
                            {
                                special = cx_matrix<value_type, Columns - 1, 1>{};
                                specialCreated = true;
                            }
                            special(k, 0) = -result(k, k);
                        }
                    }
                }
                os << "Special solution\n" << special << '\n';
            }
            os << "Particular solution:\n" << particular << '\n';
        }
    }


    constexpr cx_matrix inverse() const
    {
        auto I = create_identity();
        auto AI = augment(*this, I);
        auto U = AI.gaussian_eliminate();
        auto IAInverse = U.row_reduce_from_gaussian();

        cx_matrix AInverse{};
        for (size_type i = 0; i < Rows; ++i)
        {
            for (size_type j = 0; j < Columns; ++j)
            {
                AInverse(i, j) = IAInverse(i, j + Columns);
            }
        }
        return AInverse;
    }

    
    constexpr value_type calculate_elements_sum() const noexcept
    {
        value_type sum{};
        for (const auto& row: _data)
        {
            const auto sum_in_row = detail::accumulate(std::begin(row), std::end(row),
                                                       value_type{});
            sum += sum_in_row;
        }
        return sum;
    }


    template <class T = value_type, std::size_t Rows_E = Rows, std::size_t Columns_E = Columns>
    static constexpr cx_matrix<T, Rows_E, Columns_E> get_error_matrix() noexcept
    {
        cx_matrix<T, Rows_E, Columns_E> err_matrix{};
        for (auto& row: err_matrix._data)
        {
            detail::fill(std::begin(row), std::end(row), std::numeric_limits<T>::quiet_NaN());
        }
        return err_matrix;
    }


    friend std::ostream& operator<<(std::ostream& os, const cx_matrix& mat)
    {
        os << "[" << mat.get_dimension() << "]\n";
        for (const auto& row : mat._data)
        {
            std::copy(std::begin(row), std::end(row),
                      std::ostream_iterator<value_type>(os, " "));
            os << '\n';
        }
        return os;
    }


    friend std::istream& operator>>(std::istream& is, cx_matrix& mat)
    {
        for (auto& row : mat._data)
        {
            for (auto& elem : row)
            {
                is >> elem;
            }
        }
        return is;
    }


private:
    constexpr auto _exp_helper(const cx_matrix& mat, const int value) const noexcept
    {
        if (value == 0)
        { 
            return create_identity<value_type, Rows>();
        }
        if (value == 1)
        {
            return mat;
        }
        if (value % 2 == 0)
        {
            // value is even.
            return _exp_helper(mat * mat, value / 2);
        }
        // value is odd.
        return mat * _exp_helper(mat * mat, (value - 1) / 2);
    }

    container<container<value_type, Columns>, Rows> _data;
};


namespace detail::cx_matrix
{
    template <class Type, std::size_t N>
    using container = std::array<Type, N>;
} // namespace detail::matrix


template <class value_type, std::size_t Rows, std::size_t Columns>
constexpr cx_matrix<value_type, Rows, Columns> operator+(
    const cx_matrix<value_type, Rows, Columns>& lhs,
    const cx_matrix<value_type, Rows, Columns>& rhs) noexcept
{
    cx_matrix<value_type, Rows, Columns> temp(lhs);
    return (temp += rhs);
}


template <class value_type, std::size_t Rows, std::size_t Columns>
constexpr cx_matrix<value_type, Rows, Columns> operator-(
    const cx_matrix<value_type, Rows, Columns>& lhs,
    const cx_matrix<value_type, Rows, Columns>& rhs) noexcept
{
    cx_matrix<value_type, Rows, Columns> temp(lhs);
    return (temp -= rhs);
}


template <class value_type, std::size_t Rows, std::size_t Columns>
constexpr cx_matrix<value_type, Rows, Columns> operator*(
    const cx_matrix<value_type, Rows, Columns>& mat, const value_type value) noexcept
{
    cx_matrix<value_type, Rows, Columns> temp(mat);
    return (temp *= value);
}


template <class value_type, std::size_t Rows, std::size_t Columns>
constexpr cx_matrix<value_type, Rows, Columns> operator*(
    const value_type value, const cx_matrix<value_type, Rows, Columns>& mat) noexcept
{
    return (mat * value);
}


template <class value_type, std::size_t Rows_lhs, std::size_t Mid_dimension,
          std::size_t Columns_rhs>
constexpr cx_matrix<value_type, Rows_lhs, Columns_rhs> operator*(
    const cx_matrix<value_type, Rows_lhs, Mid_dimension>& lhs,
    const cx_matrix<value_type, Mid_dimension, Columns_rhs>& rhs) noexcept
{
    using size_type = typename cx_matrix<value_type, Rows_lhs, Columns_rhs>::size_type;

    cx_matrix<value_type, Rows_lhs, Columns_rhs> result{};
    detail::cx_matrix::container<value_type, Mid_dimension> thatColumn{};

    for (size_type j = 0; j < Columns_rhs; ++j)
    {
        for (size_type k = 0; k < Mid_dimension; ++k)
        {
            thatColumn.at(k) = rhs(k, j);
        }

        for (size_type i = 0; i < Rows_lhs; ++i)
        {
            const auto thisRow = lhs(i);
            value_type summand{};
            for (size_type k = 0; k < Mid_dimension; ++k)
            {
                summand += thisRow.at(k) * thatColumn.at(k);
            }
            result(i, j) = summand;
        }
    }
    return result;
}


template <class value_type, std::size_t Rows, std::size_t Columns>
constexpr cx_matrix<value_type, Rows, Columns> operator/(
    const cx_matrix<value_type, Rows, Columns>& mat, const value_type value)
{
    assert(value != value_type{});

    cx_matrix<value_type, Rows, Columns> temp(mat);
    return (temp /= value);
}


template <class value_type, std::size_t Rows, std::size_t Columns>
constexpr bool operator==(const cx_matrix<value_type, Rows, Columns>& lhs,
                          const cx_matrix<value_type, Rows, Columns>& rhs) noexcept
{
    using size_type = typename cx_matrix<value_type, Rows, Columns>::size_type;

    for (size_type row_index = 0; row_index < Rows; ++row_index)
    {
        if (!detail::equal(std::begin(lhs(row_index)), std::end(lhs(row_index)),
            std::begin(rhs(row_index))))
        {
            return false;
        }
    }
    return true;
}


template <class value_type, std::size_t Rows, std::size_t Columns>
constexpr bool operator!=(const cx_matrix<value_type, Rows, Columns>& lhs,
                          const cx_matrix<value_type, Rows, Columns>& rhs) noexcept
{
    return !(lhs == rhs);
}


/// Helpers operation
template <class value_type, std::size_t N>
detail::cx_matrix::container<value_type, N>
operator+(const detail::cx_matrix::container<value_type, N>& lhs,
          const detail::cx_matrix::container<value_type, N>& rhs) noexcept
{
    detail::cx_matrix::container<value_type, N> temp(lhs);
    for (std::size_t i = 0; i < temp.size(); ++i)
    {
        temp.at(i) += rhs.at(i);
    }
    return temp;
}

template <class value_type, std::size_t N>
detail::cx_matrix::container<value_type, N>
operator+(const detail::cx_matrix::container<value_type, N>& cont, const value_type& value) noexcept
{
    detail::cx_matrix::container<value_type, N> temp(cont);
    for (auto& elem : temp)
    {
        elem += value;
    }
    return temp;
}

template <class value_type, std::size_t N>
detail::cx_matrix::container<value_type, N>
operator+(const value_type& value, const detail::cx_matrix::container<value_type, N>& cont) noexcept
{
    return operator+(cont, value);
}

template <class value_type, std::size_t N>
detail::cx_matrix::container<value_type, N>
operator-(const detail::cx_matrix::container<value_type, N>& lhs,
          const detail::cx_matrix::container<value_type, N>& rhs) noexcept
{
    detail::cx_matrix::container<value_type, N> temp(lhs);
    for (std::size_t i = 0; i < temp.size(); ++i)
    {
        temp.at(i) -= rhs.at(i);
    }
    return temp;
}

template <class value_type, std::size_t N>
detail::cx_matrix::container<value_type, N>
operator-(const detail::cx_matrix::container<value_type, N>& cont, const value_type& value) noexcept
{
    detail::cx_matrix::container<value_type, N> temp(cont);
    for (auto& elem : temp)
    {
        elem -= value;
    }
    return temp;
}

template <class value_type, std::size_t N>
detail::cx_matrix::container<value_type, N>
operator-(const value_type& value, const detail::cx_matrix::container<value_type, N>& cont) noexcept
{
    return operator-(cont, value);
}

template <class value_type, std::size_t N>
detail::cx_matrix::container<value_type, N>
operator*(const detail::cx_matrix::container<value_type, N>& cont, const value_type& value) noexcept
{
    detail::cx_matrix::container<value_type, N> temp(cont);
    for (auto& elem : temp)
    {
        elem *= value;
    }
    return temp;
}

template <class value_type, std::size_t N>
detail::cx_matrix::container<value_type, N>
operator*(const value_type& value, const detail::cx_matrix::container<value_type, N>& cont) noexcept
{
    return operator*(cont, value);
}

template <class value_type, std::size_t N>
detail::cx_matrix::container<value_type, N>
operator/(const detail::cx_matrix::container<value_type, N>& cont, const value_type& value)
{
    assert(value != value_type{});

    detail::cx_matrix::container<value_type, N> temp(cont);
    for (auto& elem : temp)
    {
        elem /= value;
    }
    return temp;
}

template <class value_type, std::size_t N>
detail::cx_matrix::container<value_type, N>
operator/(const value_type& value, const detail::cx_matrix::container<value_type, N>& cont)
{
    return operator/(cont, value);
}

} // namespace vv
