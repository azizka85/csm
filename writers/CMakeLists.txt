cmake_minimum_required(VERSION 3.16)

project(writers LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(writers STATIC
    utils.cpp

    statistics.cpp
    surface.cpp
    volume.cpp
)

target_include_directories(writers PUBLIC ${CMAKE_SOURCE_DIR})

target_link_libraries(writers PRIVATE
    calc
    utils         
)