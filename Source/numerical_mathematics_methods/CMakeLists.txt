cmake_minimum_required(VERSION 3.8)

project(numerical_mathematics_methods_impl)

# Add source to this project's executable.
add_executable(${PROJECT_NAME}
    constants.hpp
    control_unit.hpp

    general_tests.hpp

    lu_decomposition.hpp
    lup_decomposition.hpp
    lupq_decomposition.hpp

    qr_decomposition.hpp
    
    iterations_methods.hpp

    newton_method.hpp

    newton_cotes.hpp
    gauss_integral.hpp

    main.cpp
)

target_link_libraries(${PROJECT_NAME}
    cx_math
    utils
    cx_options
    poly34
    cx_matrix
)


# Old script
cmake_minimum_required(VERSION 3.8)

project(numerical_mathematics_methods)


# Set compiler options.
if (MSVC)
    add_compile_options(/W4)
else ()
    add_compile_options(-Wall -Wextra -pedantic)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Add cx_math library.
add_library(cx_math STATIC
    cx_math/cx_math.h
)
target_include_directories(cx_math PUBLIC
    cx_math
)
set_target_properties(cx_math PROPERTIES LINKER_LANGUAGE CXX)

# Add utils library.
add_library(utils STATIC
    utils/utils.hpp
)
target_include_directories(utils PUBLIC
    utils
)
set_target_properties(utils PROPERTIES LINKER_LANGUAGE CXX)

# Add cx_options library.
add_library(cx_options STATIC
    cx_options/cx_loops.hpp
    cx_options/cx_random.hpp
    cx_options/cx_algorthms.hpp
)
target_include_directories(cx_options PUBLIC
    cx_options
)
set_target_properties(cx_options PROPERTIES LINKER_LANGUAGE CXX)

# Add poly34 library.
add_library(poly34 STATIC
    poly34/poly34.h
)
target_include_directories(poly34 PUBLIC
    poly34
)
set_target_properties(poly34 PROPERTIES LINKER_LANGUAGE CXX)

# Add cx_matrix library.
add_library(cx_matrix STATIC
    cx_matrix/cx_matrix.hpp
)
target_include_directories(cx_matrix PUBLIC
    cx_matrix
)
set_target_properties(cx_matrix PROPERTIES LINKER_LANGUAGE CXX)


add_executable(${PROJECT_NAME}
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

    numerical_mathematics_methods/main.cpp
)

target_link_libraries(${PROJECT_NAME}
    cx_math
    utils
    cx_options
    poly34
    cx_matrix
)


# Include sub-projects.
#add_subdirectory(numerical_mathematics_methods)
