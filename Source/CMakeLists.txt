cmake_minimum_required(VERSION 3.8)

project(numerical_mathematics_methods)


if(WIN32)
    set(CMAKE_PREFIX_PATH "C:\\Users\\vasar\\Documents\\Visual Studio 2017\\Projects\\numerical_mathematics_methods\\SFML-2.4.2\\")
endif()

# Set compiler options.
if (MSVC)
    add_compile_options(/W4)
else ()
    add_compile_options(-Wall -Wextra -pedantic)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


include_directories(
    ${PROJECT_BINARY_DIR}
    cx_math
    utils
    cx_options
    poly34
    cx_matrix
    matrix
)

set(EXECUTABLE_NAME "numerical_mathematics_methods")

add_executable(${EXECUTABLE_NAME}
    cx_math/cx_math.h
    utils/utils.hpp
    cx_options/cx_loops.hpp
    cx_options/cx_random.hpp
    cx_options/cx_algorthms.hpp
    poly34/poly34.h
    cx_matrix/cx_matrix.hpp
    matrix/matrix.hpp


    numerical_mathematics_methods/constants.hpp
    numerical_mathematics_methods/control_unit.hpp

    numerical_mathematics_methods/general_tests.hpp

    numerical_mathematics_methods/lu_decomposition.hpp
    numerical_mathematics_methods/lup_decomposition.hpp
    numerical_mathematics_methods/lupq_decomposition.hpp

    numerical_mathematics_methods/qr_decomposition.hpp
    
    numerical_mathematics_methods/iterations_methods.hpp

    numerical_mathematics_methods/newton_method.hpp

    numerical_mathematics_methods/newton_cotes.hpp
    numerical_mathematics_methods/gauss_integral.hpp

    numerical_mathematics_methods/runge_kutta_methods.hpp

    numerical_mathematics_methods/transmitter.hpp

    numerical_mathematics_methods/main.cpp
)

# Detect and add SFML
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules" ${CMAKE_MODULE_PATH})
#Find any version 2.X of SFML
#See the FindSFML.cmake file for additional details and instructions
find_package(SFML 2.4 REQUIRED graphics window system)
if(SFML_FOUND)
    include_directories(${SFML_INCLUDE_DIR})
    target_link_libraries(${EXECUTABLE_NAME} ${SFML_LIBRARIES} ${SFML_DEPENDENCIES})
endif()


# Install target
install(TARGETS ${EXECUTABLE_NAME} DESTINATION bin)
