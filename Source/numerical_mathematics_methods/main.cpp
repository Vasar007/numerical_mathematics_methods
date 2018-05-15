#include <iostream>

#include "utils.hpp"
#include "control_unit.hpp"


int main()
{
    std::cout << "Numerical Mathematics Methods\n";
    utils::pause("Press Enter to start computations\n");

    std::cout.precision(15);
    std::cout << std::fixed;

    TRY_BLOCK(vv::start_tests();)

    std::cout << "\n\n";
    utils::pause();
    return 0;
}
