cmake_minimum_required(VERSION 3.31.5)
project(ucDE VERSION 0.1.0 LANGUAGES C)

add_executable(ucDE ucDE_endCom.c mt64.h mt19937-64.c erf.h erf.c compact_lib.h compact_lib.c cec17_test_func.c)

find_library(MATH_LIBRARY m)
if(MATH_LIBRARY)
    target_link_libraries(ucDE PUBLIC ${MATH_LIBRARY})
endif()