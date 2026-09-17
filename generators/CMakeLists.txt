cmake_minimum_required(VERSION 3.16)

project(generators LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(generators STATIC
    area/uniform_generator.cpp

    bathymetry/uniform_generator.cpp

    dz/uniform_generator.cpp
    dz/triple_point_generator.cpp

    viscosity/uniform_generator.cpp
    viscosity/linear_generator.cpp
    viscosity/wind_speed_generator.cpp

    wind/uniform_generator.cpp
)

target_include_directories(generators PUBLIC ${CMAKE_SOURCE_DIR})

target_link_libraries(generators PRIVATE
    calc
    utils
)