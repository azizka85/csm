cmake_minimum_required(VERSION 3.16)

project(utils LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(utils STATIC 
	checks/bathymetry.cpp
	checks/vector.cpp
	checks/time.cpp
	checks/phys.cpp
)

target_include_directories(utils PUBLIC ${CMAKE_SOURCE_DIR})
