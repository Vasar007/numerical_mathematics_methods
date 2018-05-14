#include <iostream>
#include <fstream>
#include <string>
#include <string_view>
#include <vector>

#include <SFML/Graphics.hpp>

#include "utils.hpp"
#include "control_unit.hpp"
#include "transmitter.hpp"


void out_data(const std::string& file_name, const std::string_view mode,
              const std::string_view param, const std::string_view title,
              const std::string_view x_label, const std::string_view y_label,
              const std::vector<std::pair<double, double>>& data)
{
    if (data.empty())
    {
        std::cout << "ERROR: empty data to process for file " << file_name << ".\n";
        return;
    }

    std::ofstream out_file(file_name);
    out_file << mode << '|' << param << '|' << title << '|' << x_label << '|' << y_label << '\n';
    for (const auto [x, y] : data)
    {
        out_file << x << ' ' << y << '\n';
    }
}


void out_data(const std::string& file_name, const std::string_view mode,
              const std::string_view param, const std::string_view title,
              const std::string_view x_label, const std::string_view y_label,
              const std::vector<std::pair<double, double>>& data_1,
              const std::vector<std::pair<double, double>>& data_2)
{
    if (data_1.empty() || data_2.empty())
    {
        std::cout << "ERROR: empty data to process for file " << file_name << ".\n";
        return;
    }

    std::ofstream out_file(file_name);
    out_file << mode << '|' << param << '|' << title << '|' << x_label << '|' << y_label << '\n';
    for (const auto [x, y] : data_1)
    {
        out_file << x << ' ' << y << '\n';
    }
    out_file << "#\n";
    for (const auto [x, y] : data_2)
    {
        out_file << x << ' ' << y << '\n';
    }
}


void out_data(const std::string& file_name, const std::string_view mode,
              const std::string_view param, const std::string_view title,
              const std::string_view x_label, const std::string_view y_label,
              const vv::matrix<double>& x, const vv::matrix<double>& y)
{
    if (x.empty() || y.empty())
    {
        std::cout << "ERROR: empty data to process for file " << file_name << ".\n";
        return;
    }
    if (x.get_rows_number() != y.get_rows_number())
    {
        std::cout << "\nERROR: X and Y have different size! Cannot process data.\n";
        return;
    }

    std::ofstream out_file(file_name);
    out_file << mode << '|' << param << '|' << title << '|' << x_label << '|' << y_label << '\n';
    for (std::size_t i = 0; i < x.get_rows_number(); ++i)
    {
        out_file << x(i, 0) << ' ' << y(i, 0) << ' ' << y(i, 1) << ' ' << y(i, 2) << ' ' << y(i, 3)
            << '\n';
    }
}


void out_data(const std::string& file_name, const std::string_view mode,
              const std::string_view param, const std::string_view title,
              const std::string_view x_label, const std::string_view y_label,
              const vv::matrix<double>& x1, const vv::matrix<double>& y1,
              const vv::matrix<double>& x2, const vv::matrix<double>& y2)
{
    if (x1.empty() || y1.empty() || x2.empty() || y2.empty())
    {
        std::cout << "ERROR: empty data to process for file " << file_name << ".\n";
        return;
    }
    if (x1.get_rows_number() != y1.get_rows_number()
     || x2.get_rows_number() != y2.get_rows_number())
    {
        std::cout << "\nERROR: X and Y have different size! Cannot process data.\n";
        return;
    }

    std::ofstream out_file(file_name);
    out_file << mode << '|' << param << '|' << title << '|' << x_label << '|' << y_label << '\n';
    for (std::size_t i = 0; i < x1.get_rows_number(); ++i)
    {
        out_file << x1(i, 0) << ' ' << y1(i, 0) << ' ' << y1(i, 1) << ' ' << y1(i, 2) << ' '
            << y1(i, 3) << '\n';
    }
    out_file << "#\n";
    for (std::size_t i = 0; i < x2.get_rows_number(); ++i)
    {
        out_file << x2(i, 0) << ' ' << y2(i, 0) << ' ' << y2(i, 1) << ' ' << y2(i, 2) << ' '
            << y2(i, 3) << '\n';
    }
}


int main()
{
    utils::pause("Press Enter to start computations\n");

    std::cout.precision(15);
    std::cout << std::fixed;

    //TRY_BLOCK(vv::start_tests();)

    constexpr std::size_t width = 1600;
    constexpr std::size_t height = 900;
    sf::RenderWindow window(sf::VideoMode(width, height), "Numerical Mathematics Methods");
    window.setFramerateLimit(60);

    std::vector<statistic::transmitter> graph;
    while (window.isOpen())
    {
        sf::Event event{};
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
            {
                window.close();
            }
            else if (event.type == sf::Event::KeyPressed)
            {
                switch (event.key.code)
                {
                    /// Part 1.
                    case sf::Keyboard::Q: // x = log2(h), y = Rn with h = 1/2^k steps.
                        graph.clear();
                        graph.emplace_back(width, height, vv::graphic_data,
                                           sf::PrimitiveType::LineStrip, sf::Color::Black);
                        graph.emplace_back(width, height, vv::graphic_data_opp,
                                           sf::PrimitiveType::LineStrip, sf::Color::Magenta);
                        break;

                    case sf::Keyboard::W: // tolerance with h_opt.
                        graph.clear();
                        graph.emplace_back(width, height, vv::tol_graphic_data,
                                           sf::PrimitiveType::LineStrip, sf::Color::Black);
                        graph.emplace_back(width, height, vv::tol_graphic_data_opp,
                                           sf::PrimitiveType::LineStrip, sf::Color::Magenta);
                        break;

                    /// Part 2.
                    case sf::Keyboard::E: // tolerance with h_auto for p = 2.
                        graph.clear();
                        graph.emplace_back(width, height, vv::tol_graphic_data_auto,
                                           sf::PrimitiveType::LineStrip, sf::Color::Black);
                        break;

                    case sf::Keyboard::R: // plot for my solution for p = 2.
                        graph.clear();
                        graph.emplace_back(width, height, vv::graphic_data_x, vv::graphic_data_y,
                                           sf::PrimitiveType::Lines, sf::Color::Black);
                        break;

                    case sf::Keyboard::T: // quality for p = 2.
                        graph.clear();
                        graph.emplace_back(width, height, vv::quality,
                                           sf::PrimitiveType::LineStrip, sf::Color::Black);
                        break;

                    case sf::Keyboard::Y: // h_acc and h_not_acc for p = 2.
                        graph.clear();
                        graph.emplace_back(width, height, vv::h_acc,
                                           sf::PrimitiveType::Points, sf::Color::Black);
                        graph.emplace_back(width, height, vv::h_not_acc,
                                           sf::PrimitiveType::Points, sf::Color::Magenta);
                        break;

                    /// Part 3.
                    case sf::Keyboard::U: // tolerance with h_auto for p = 4.
                        graph.clear();
                        graph.emplace_back(width, height, vv::tol_graphic_data_auto_opp,
                                           sf::PrimitiveType::LineStrip, sf::Color::Black);
                        break;

                    case sf::Keyboard::I: // plot for my solution for p = 4.
                        graph.clear();
                        graph.emplace_back(width, height, vv::graphic_data_x_opp, vv::graphic_data_y_opp,
                                           sf::PrimitiveType::Lines, sf::Color::Black);
                        break;

                    case sf::Keyboard::O: // quality for p = 4.
                        graph.clear();
                        graph.emplace_back(width, height, vv::quality_opp,
                                           sf::PrimitiveType::LineStrip, sf::Color::Black);
                        break;

                    case sf::Keyboard::P: // h_acc and h_not_acc for p = 4.
                        graph.clear();
                        graph.emplace_back(width, height, vv::h_acc_opp,
                                           sf::PrimitiveType::Points, sf::Color::Black);
                        graph.emplace_back(width, height, vv::h_not_acc_opp,
                                           sf::PrimitiveType::Points, sf::Color::Magenta);
                        break;

                    case sf::Keyboard::L: // count addressing to F.
                        graph.clear();
                        graph.emplace_back(width, height, vv::count_F,
                                           sf::PrimitiveType::LineStrip, sf::Color::Black);
                        graph.emplace_back(width, height, vv::count_F_opp,
                                           sf::PrimitiveType::LineStrip, sf::Color::Magenta);
                        break;

                    /// Calculation settings.
                    case sf::Keyboard::Num1: // calculate RKM with p = 2 and then calculate with h_opt.
                        vv::graphic_data.clear();
                        vv::runge_kutta_calculation();
                        break;

                    case sf::Keyboard::Num2: // calculate RKM with p = 4 and then calculate with h_opt.
                        vv::graphic_data_opp.clear();
                        vv::runge_kutta_calculation_opp();
                        break;

                    case sf::Keyboard::Num3: // calculate RKM_auto with p = 2 and also save plot and quality with h_acc and h_not_acc.
                        vv::quality.clear();
                        vv::h_acc.clear();
                        vv::h_not_acc.clear();
                        vv::runge_kutta_calculation_auto();
                        break;

                    case sf::Keyboard::Num4: // calculate RKM_auto with p = 4 and also save plot and quality with h_acc and h_not_acc.
                        vv::quality_opp.clear();
                        vv::h_acc_opp.clear();
                        vv::h_not_acc_opp.clear();
                        vv::runge_kutta_calculation_auto_opp();
                        break;

                    case sf::Keyboard::Num5: // calculate RKM_auto and count addressing to F.
                        vv::count_F.clear();
                        vv::runge_kutta_calculation_auto_rtol();
                        vv::quality.clear();
                        vv::h_acc.clear();
                        vv::h_not_acc.clear();

                        vv::count_F_opp.clear();
                        vv::runge_kutta_calculation_auto_rtol_opp();
                        vv::quality_opp.clear();
                        vv::h_acc_opp.clear();
                        vv::h_not_acc_opp.clear();
                        break;

                    case sf::Keyboard::S:
                        out_data("rkm_h.txt", "2p", "def",
                                 "RKM with h = 1/2^k", "-log2(h)", "y",
                                 vv::graphic_data, vv::graphic_data_opp);
                        out_data("rkm_h_opt.txt", "2p", "def", "RKM with h_opt", "x", "y - y1",
                                 vv::tol_graphic_data, vv::tol_graphic_data_opp);

                        out_data("rkm_2_auto.txt", "1p", "def",
                                 "Tolerance for RKM_auto with p = 2", "x", "y - y1",
                                 vv::tol_graphic_data_auto);
                        out_data("rkm_2_auto_plot.txt", "8p", "r;r;r;r;b;b;b;b",
                                 "Plot from RKM_auto with p = 2", "x", "y",
                                 vv::graphic_data_x, vv::graphic_data_y,
                                 vv::data_x_true_opp, vv::data_y_true_opp);
                        out_data("rkm_2_auto_quality.txt", "1p", "def",
                                 "Quality for RKM_auto with p = 2", "x", "Quality",
                                 vv::quality);
                        out_data("rkm_2_auto_h.txt", "2p", "r,;b,",
                                 "Choice 'h' for RKM_auto with p = 2", "x", "h",
                                 vv::h_acc, vv::h_not_acc);

                        out_data("rkm_4_auto.txt", "1p", "def",
                                 "Tolerance for RKM_auto with p = 4", "x", "y - y1",
                                 vv::tol_graphic_data_auto_opp);
                        out_data("rkm_4_auto_plot.txt", "8p", "r;r;r;r;b;b;b;b",
                                 "Plot from RKM_auto with p = 4", "x", "y",
                                 vv::graphic_data_x_opp, vv::graphic_data_y_opp,
                                 vv::data_x_true_opp, vv::data_y_true_opp);
                        out_data("rkm_4_auto_quality.txt", "1p", "def",
                                 "Quality for RKM_auto with p = 4", "x", "Quality",
                                 vv::quality_opp);
                        out_data("rkm_4_auto_h.txt", "2p", "r,;b,",
                                 "Choice 'h' for RKM_auto with p = 4", "x", "h",
                                 vv::h_acc_opp, vv::h_not_acc_opp);
                        break;

                    case sf::Keyboard::D:
                        out_data("rkm_auto_addressing.txt", "2p", "def",
                                 "Counting addresses to F", "-log10(rtol)", "log10(addressing_F)",
                                 vv::count_F, vv::count_F_opp);
                        break;

                    default:
                        break;
                }
            }
        }

        window.clear(sf::Color::White);
        for (const auto& elem : graph)
        {
            window.draw(elem);
        }
        window.display();
    }
    return 0;
}
