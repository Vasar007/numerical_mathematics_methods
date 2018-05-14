#pragma once

#include <vector>

#include <SFML/Graphics.hpp>

#include "matrix.hpp"


namespace statistic
{

class transmitter final : public sf::Drawable
{
public:
    transmitter(const int max_x, const int max_y,
                const std::vector<std::pair<double, double>>& data,
                const sf::PrimitiveType type = sf::Lines, const sf::Color color = sf::Color::Black)
    : _limit_x(max_x)
    , _limit_y(max_y)
    , _color(color)
    , _mode(Mode::One)
    , _coord_system(sf::Lines, 4)
    , _shapes(sf::Lines)
    , _points(type)
    , _points_y1(type)
    , _points_y2(type)
    , _points_y3(type)
    , _points_y4(type)
    {
        initialize_coord_system();
        add_points(data);
    }

    transmitter(const int max_x, const int max_y, const vv::matrix<double>& x,
                const vv::matrix<double>& y, const sf::PrimitiveType type = sf::Lines,
                const sf::Color color = sf::Color::Black)
    : _limit_x(max_x)
    , _limit_y(max_y)
    , _color(color)
    , _mode(Mode::Four)
    , _coord_system(sf::Lines, 4)
    , _shapes(sf::Lines)
    , _points(type)
    , _points_y1(type)
    , _points_y2(type)
    , _points_y3(type)
    , _points_y4(type)
    {
        initialize_coord_system();
        add_points(x, y);
    }

    void reset_data(const std::vector<std::pair<double, double>>& data,
                    const sf::PrimitiveType type = sf::LineStrip,
                    const sf::Color color = sf::Color::Black)
    {
        _color = color;
        _points = sf::VertexArray(type);
        _points_y1.clear();
        _points_y2.clear();
        _points_y3.clear();
        _points_y4.clear();
        _mode = Mode::One;
        add_points(data);
    }

    void reset_data(const vv::matrix<double>& x, const vv::matrix<double>& y,
                    const sf::PrimitiveType type = sf::Lines,
                    const sf::Color color = sf::Color::Black)
    {
        _color = color;
        _points.clear();
        _points_y1 = sf::VertexArray(type);
        _points_y2 = sf::VertexArray(type);
        _points_y3 = sf::VertexArray(type);
        _points_y4 = sf::VertexArray(type);
        _mode = Mode::Four;
        add_points(x, y);
    }

private:
    void draw(sf::RenderTarget& target, sf::RenderStates) const override
    {
        target.draw(_coord_system);
        target.draw(_shapes);

        if (_mode == Mode::One)
        {
            target.draw(_points);
        }
        else if (_mode == Mode::Four)
        {
            target.draw(_points_y1);
            target.draw(_points_y2);
            target.draw(_points_y3);
            target.draw(_points_y4);
        }
    }

    void initialize_coord_system()
    {
        _coord_system[0] = { sf::Vector2f(1.0f, 1.0f), sf::Color::Red };
        _coord_system[1] = { sf::Vector2f(1.0f, static_cast<float>(_limit_y)), sf::Color::Red };
        _coord_system[2] = { sf::Vector2f(1.0f, _limit_y - 1.0f - _limit_y / 2.0f), sf::Color::Red };
        _coord_system[3] = { sf::Vector2f(_limit_x + 1.0f, _limit_y - 1.0f - _limit_y / 2.0f), sf::Color::Red };

        for (int i = 0; i < segments; ++i)
        {
            _shapes.append({ sf::Vector2f(1.0f * (i + 1) * _limit_x / segments, _limit_y - 6.0f - _limit_y / 2.0f), sf::Color::Red });
            _shapes.append({ sf::Vector2f(1.0f * (i + 1) * _limit_x / segments, _limit_y - 1.0f - _limit_y / 2.0f), sf::Color::Red });
        }
        for (int i = 0; i < segments; ++i)
        {
            _shapes.append({ sf::Vector2f(1.0f, 1.0f * (i + 1) * _limit_y / segments), sf::Color::Red });
            _shapes.append({ sf::Vector2f(6.0f, 1.0f * (i + 1) * _limit_y / segments), sf::Color::Red });
        }
    }

    void add_points(const std::vector<std::pair<double, double>>& data)
    {
        for (const auto point : data)
        {
            const auto coord_x = static_cast<float>(point.first * 50.0f);
            const auto coord_y = static_cast<float>(_limit_y - point.second * 10.0f) - 1.0f - _limit_y / 2.0f;

            _points.append({ sf::Vector2f(coord_x, coord_y), _color });
        }
    }

    void add_points(const vv::matrix<double>& x, const vv::matrix<double>& y)
    {
        if (x.get_rows_number() != y.get_rows_number())
        {
            std::cout << "\nERROR: X and Y have different size! Cannot process data.\n";
            return;
        }

        for (std::size_t i = 0; i < x.get_rows_number(); ++i)
        {
            const auto coord_x = static_cast<float>(x(i, 0) * 50.0f);
            const auto coord_y1 = static_cast<float>(_limit_y - y(i, 0) * 10.0f) - 1.0f - _limit_y / 2.0f;
            const auto coord_y2 = static_cast<float>(_limit_y - y(i, 1) * 10.0f) - 1.0f - _limit_y / 2.0f;
            const auto coord_y3 = static_cast<float>(_limit_y - y(i, 2) * 10.0f) - 1.0f - _limit_y / 2.0f;
            const auto coord_y4 = static_cast<float>(_limit_y - y(i, 3) * 10.0f) - 1.0f - _limit_y / 2.0f;

            _points_y1.append({ sf::Vector2f(coord_x, coord_y1), _color });
            _points_y2.append({ sf::Vector2f(coord_x, coord_y2), _color });
            _points_y3.append({ sf::Vector2f(coord_x, coord_y3), _color });
            _points_y4.append({ sf::Vector2f(coord_x, coord_y4), _color });
        }
    }


    static constexpr int segments = 50;

    enum class Mode
    {
        One,
        Four
    };

    const int _limit_x;

    const int _limit_y;

    Mode _mode;

    sf::Color _color;

    sf::VertexArray _coord_system;

    sf::VertexArray _shapes;

    sf::VertexArray _points;
    sf::VertexArray _points_y1;
    sf::VertexArray _points_y2;
    sf::VertexArray _points_y3;
    sf::VertexArray _points_y4;

};

} // namespace stat
