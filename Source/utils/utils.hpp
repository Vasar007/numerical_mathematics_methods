// Copyright (C) 2018 Vasily Vasilyev (vasar007@yandex.ru)

#pragma once

#include <iostream>
#include <iterator>
#include <string_view>


namespace utils
{

template <class OutputStream, class Container>
void print(OutputStream& out, const Container& container)
{
    std::copy(std::begin(container), std::end(container),
              std::ostream_iterator<typename Container::value_type>(out, " "));
}


template <class OutputStream, class Container>
void println(OutputStream& out, const Container& container)
{
    std::copy(std::begin(container), std::end(container),
              std::ostream_iterator<typename Container::value_type>(out, " "));
    std::cout << '\n';
}


void pause(const std::string_view message = "\nPress the Enter key to continue...")
{
    do
    {
        std::cout << message;
    }
    while (std::cin.get() != '\n');
}


void pause_clear(const std::string_view message = "Press ENTER to continue...")
{
    std::cout << message << std::flush;
    std::cin.seekg(0u, std::ios::end);
    std::cin.clear();
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
}

} // namespace utils
