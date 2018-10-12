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

    numerical_mathematics_methods/cx_functions.hpp

    numerical_mathematics_methods/lu_decomposition.hpp
    numerical_mathematics_methods/lup_decomposition.hpp
    numerical_mathematics_methods/lupq_decomposition.hpp

    numerical_mathematics_methods/qr_decomposition.hpp
    
    numerical_mathematics_methods/iterations_methods.hpp

    numerical_mathematics_methods/newton_method.hpp

    numerical_mathematics_methods/newton_cotes.hpp
    numerical_mathematics_methods/gauss_integral.hpp

    numerical_mathematics_methods/runge_kutta_methods.hpp

    numerical_mathematics_methods/main.cpp
)
