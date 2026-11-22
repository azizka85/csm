cmake_minimum_required(VERSION 3.28)

project(utils LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(MSVC)
    set(OpenMP_RUNTIME_MSVC "llvm")
endif()

find_package(OpenMP REQUIRED)

add_library(utils STATIC 
	checks/bathymetry.cpp
	checks/vector.cpp
	checks/time.cpp
	checks/phys.cpp
	checks/area.cpp

	checks/fs.cpp

	data.cpp
	bathymetry.cpp
	viscosity.cpp
)

target_include_directories(utils PUBLIC ${CMAKE_SOURCE_DIR})

target_link_libraries(utils PUBLIC OpenMP::OpenMP_CXX)
