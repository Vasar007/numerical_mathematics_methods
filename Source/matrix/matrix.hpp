#pragma once

#include <iterator>
#include <algorithm>
#include <vector>
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
 * As base for matrix class using https://github.com/akalicki/matrix
 *
 */

namespace vv
{
/// ===== MATRIX IMLEMENTATION SECTION =====

namespace mtx
{

enum class mode
{
    normal,
    dynamically_expandable
};

} // namespace mtx

    
template <class Type = double>
class matrix
{
public:
    using value_type                    = Type;
    using size_type                     = std::size_t;
    using difference_type               = std::ptrdiff_t;

    template <class ConType = value_type>
    using container                     = std::vector<ConType>;
    using row_container                 = container<value_type>;
    using row_container_reference       = container<value_type>&;    
    using const_row_container_reference = const container<value_type>&;

    using pointer                       = value_type*;
    using const_pointer                 = const value_type*;

    using reference                     = value_type&;
    using const_reference               = const value_type&;


    static constexpr value_type EPS = static_cast<value_type>(1e-10);

    static_assert(std::is_arithmetic_v<value_type>, "Matrix elements type has to be arithmetic!");


    matrix() = default;

    matrix(const size_type rows, const size_type columns,
           const mtx::mode work_mode = mtx::mode::normal)
    : _rows(rows)
    , _columns(columns)
    , _data(rows, container<value_type>(columns, value_type{}))
    , _mode(work_mode)
    { }

    matrix(const std::initializer_list<value_type> list)
    : _data{ list }
    , _mode(mtx::mode::normal)
    {
        _rows = _data.size();
        _columns = _data.at(0).size();
    }
    
    ~matrix() noexcept = default;
    matrix(const matrix& other) = default;
    matrix& operator=(const matrix& other) = default;
    matrix(matrix&& other) noexcept = default;
    matrix& operator=(matrix&& other) noexcept = default;


    std::string get_dimension() const
    {
        return std::to_string(get_rows_number()) + std::string("x")
               + std::to_string(get_columns_number());
    }


    size_type get_rows_number() const noexcept
    {
        return _rows;
    }


    size_type get_columns_number() const noexcept
    {
        return _columns;
    }


    mtx::mode mode() const noexcept
    {
        return _mode;
    }

    bool empty() const noexcept
    {
        return (_rows == 0) || (_columns == 0);
    }


    const_reference operator()(const size_type i, const size_type j) const
    {
        return _data.at(i).at(j);
    }


    reference operator()(const size_type i, const size_type j)
    { 
        if (_mode != mtx::mode::dynamically_expandable || i < _rows) return _data.at(i).at(j);

        if (i == _rows)
        {
            _data.push_back(container<value_type>(_columns, value_type{}));
            ++_rows;
        }
        return _data.at(i).at(j);
    }


    const_row_container_reference& operator()(const size_type i) const
    {
        return _data.at(i);
    }


    row_container_reference& operator()(const size_type i)
    {
        if (_mode != mtx::mode::dynamically_expandable || i < _rows) return _data.at(i);

        if (i == _rows)
        {
            _data.push_back(container<value_type>(_columns, value_type{}));
            ++_rows;
        }
        return _data.at(i);
    }


    matrix& operator+=(const matrix& rhs) noexcept
    {
        for (size_type i = 0; i < _rows; ++i)
        {
            for (size_type j = 0; j < _columns; ++j)
            {
                _data.at(i).at(j) += rhs(i, j);
            }
        }
        return *this;
    }
    

    matrix& operator-=(const matrix& rhs) noexcept
    {
        for (size_type i = 0; i < _rows; ++i)
        {
            for (size_type j = 0; j < _columns; ++j)
            {
                _data.at(i).at(j) -= rhs(i, j);
            }
        }
        return *this;
    }
    

    matrix& operator*=(const value_type value) noexcept
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
    

    matrix& operator/=(const value_type value)
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


    auto operator^(const int value) const
    {
        matrix temp(*this);
        return _exp_helper(temp, value);
    }


    void swap_rows(const size_type row_1, const size_type row_2) noexcept
    {
        detail::swap(_data.at(row_1), _data.at(row_2));
    }


    void swap_columns(const size_type column_1, const size_type column_2) noexcept
    {
        for (size_type i = 0; i < _rows; ++i)
        {
            detail::swap(_data.at(i).at(column_1), _data.at(i).at(column_2));
        }
    }
    

    matrix<value_type> transpose() const
    {
        matrix<value_type> transp(_columns, _rows);
        for (size_type i = 0; i < _rows; ++i)
        {
            for (size_type j = 0; j < _columns; ++j)
            {
                transp(j, i) = _data.at(i).at(j);
            }
        }
        return transp;
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


    value_type calculate_condition_number() const noexcept
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


    template <class T = value_type>
    static matrix<T> create_identity(const size_type size)
    {
        matrix<T> temp(size, size);
        for (size_type i = 0; i < size; ++i)
        {
            temp(i, i) = static_cast<T>(1);
        }
        return temp;
    }

    
    value_type calculate_elements_summ() const noexcept
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


    friend std::ostream& operator<<(std::ostream& os, const matrix& mat)
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


    friend std::istream& operator>>(std::istream& is, matrix& mat)
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
    auto _exp_helper(const matrix& mat, const int value) const
    {
        if (value == 0)
        { 
            return create_identity<value_type>();
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


    std::size_t _rows;
    std::size_t _columns;

    container<container<value_type>> _data;

    mtx::mode _mode;
};


namespace detail::matrix
{
    template <class Type>
    using container = std::vector<Type>;
} // namespace detail::matrix


template <class value_type>
matrix<value_type> operator+(const matrix<value_type>& lhs, const matrix<value_type>& rhs)
{
    matrix<value_type> temp(lhs);
    return (temp += rhs);
}


template <class value_type>
matrix<value_type> operator-(const matrix<value_type>& lhs, const matrix<value_type>& rhs)
{
    matrix<value_type> temp(lhs);
    return (temp -= rhs);
}


template <class value_type>
matrix<value_type> operator*(const matrix<value_type>& mat, const value_type value)
{
    matrix<value_type> temp(mat);
    return (temp *= value);
}


template <class value_type>
matrix<value_type> operator*(const value_type value, const matrix<value_type>& mat)
{
    return (mat * value);
}

template <class value_type>
matrix<value_type> operator*(const matrix<value_type>& lhs, const matrix<value_type>& rhs)
{
    using size_type = typename matrix<value_type>::size_type;

    if (lhs.get_columns_number() != rhs.get_rows_number())
    {
        return matrix<value_type>{};
    }

    const size_type mid_dimension = lhs.get_columns_number();
    matrix<value_type> result(lhs.get_rows_number(), rhs.get_columns_number());
    detail::matrix::container<value_type> thatColumn(mid_dimension);

    for (size_type j = 0; j < rhs.get_columns_number(); ++j)
    {
        for (size_type k = 0; k < mid_dimension; ++k)
        {
            thatColumn.at(k) = rhs(k, j);
        }

        for (size_type i = 0; i < lhs.get_rows_number(); ++i)
        {
            const auto thisRow = lhs(i);
            value_type summand{};
            for (size_type k = 0; k < mid_dimension; ++k)
            {
                summand += thisRow.at(k) * thatColumn.at(k);
            }
            result(i, j) = summand;
        }
    }
    return result;
}


template <class value_type>
matrix<value_type> operator/(const matrix<value_type>& mat, const value_type value)
{
    assert(value != value_type{});

    matrix<value_type> temp(mat);
    return (temp /= value);
}


template <class value_type>
bool operator==(const matrix<value_type>& lhs, const matrix<value_type>& rhs) noexcept
{
    using size_type = typename matrix<value_type>::size_type;

    if (lhs.get_rows_number() != rhs.get_rows_number())
    {
        return false;
    }

    for (size_type row_index = 0; row_index < lhs.get_rows_number(); ++row_index)
    {
        if (!detail::equal(std::begin(lhs(row_index)), std::end(lhs(row_index)),
            std::begin(rhs(row_index))))
        {
            return false;
        }
    }
    return true;
}


template <class value_type>
bool operator!=(const matrix<value_type>& lhs, const matrix<value_type>& rhs) noexcept
{
    if (lhs.get_rows_number() != rhs.get_rows_number())
    {
        return false;
    }
    return !(lhs == rhs);
}


/// Helpers operation
template <class value_type>
detail::matrix::container<value_type>
operator+(const detail::matrix::container<value_type>& lhs,
          const detail::matrix::container<value_type>& rhs)
{
    using size_type = typename matrix<value_type>::size_type;

    detail::matrix::container<value_type> temp(lhs);
    for (size_type i = 0; i < temp.size(); ++i)
    {
        temp.at(i) += rhs.at(i);
    }
    return temp;
}

template <class value_type>
detail::matrix::container<value_type>
operator+(const detail::matrix::container<value_type>& cont, const value_type& value)
{
    detail::matrix::container<value_type> temp(cont);
    for (auto& elem : temp)
    {
        elem += value;
    }
    return temp;
}

template <class value_type>
detail::matrix::container<value_type>
operator+(const value_type& value, const detail::matrix::container<value_type>& cont)
{
    return operator+(cont, value);
}

template <class value_type>
detail::matrix::container<value_type>
operator-(const detail::matrix::container<value_type>& lhs,
          const detail::matrix::container<value_type>& rhs)
{
    using size_type = typename matrix<value_type>::size_type;

    detail::matrix::container<value_type> temp(lhs);
    for (size_type i = 0; i < temp.size(); ++i)
    {
        temp.at(i) -= rhs.at(i);
    }
    return temp;
}

template <class value_type>
detail::matrix::container<value_type>
operator-(const detail::matrix::container<value_type>& cont, const value_type& value)
{
    detail::matrix::container<value_type> temp(cont);
    for (auto& elem : temp)
    {
        elem -= value;
    }
    return temp;
}

template <class value_type>
detail::matrix::container<value_type>
operator-(const value_type& value, const detail::matrix::container<value_type>& cont)
{
    return operator-(cont, value);
}

template <class value_type>
detail::matrix::container<value_type>
operator*(const detail::matrix::container<value_type>& cont, const value_type& value)
{
    detail::matrix::container<value_type> temp(cont);
    for (auto& elem : temp)
    {
        elem *= value;
    }
    return temp;
}

template <class value_type>
detail::matrix::container<value_type>
operator*(const value_type& value, const detail::matrix::container<value_type>& cont)
{
    return operator*(cont, value);
}

template <class value_type>
detail::matrix::container<value_type>
operator/(const detail::matrix::container<value_type>& cont, const value_type& value)
{
    assert(value != value_type{});

    detail::matrix::container<value_type> temp(cont);
    for (auto& elem : temp)
    {
        elem /= value;
    }
    return temp;
}

template <class value_type>
detail::matrix::container<value_type>
operator/(const value_type& value, const detail::matrix::container<value_type>& cont)
{
    return operator/(cont, value);
}

} // namespace vv
