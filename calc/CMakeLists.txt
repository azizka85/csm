cmake_minimum_required(VERSION 3.16)

project(calc LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(calc STATIC
    grid.cpp
    time.cpp
)

target_include_directories(calc PUBLIC ${CMAKE_SOURCE_DIR})

target_link_libraries(calc PRIVATE
    utils
)