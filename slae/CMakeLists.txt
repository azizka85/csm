cmake_minimum_required(VERSION 3.16)

project(slae LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(slae STATIC
    direct/tridiagonal.cpp
)

target_include_directories(slae PUBLIC ${CMAKE_SOURCE_DIR})

target_link_libraries(slae PRIVATE
    utils
)