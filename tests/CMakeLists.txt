cmake_minimum_required(VERSION 3.16)

project(tests LANGUAGES CXX)

find_package(Boost REQUIRED COMPONENTS unit_test_framework)

add_executable(tests 
    slae/direct/test_tridiagonal.cpp
)

target_link_libraries(tests 
    PRIVATE 
        slae
        calc
        Boost::unit_test_framework
)

enable_testing()

add_test(NAME Tests COMMAND tests)