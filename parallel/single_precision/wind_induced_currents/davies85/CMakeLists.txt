cmake_minimum_required(VERSION 3.28)

project(parallel_single_precision_wind_induced_currents_davies_85 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(OpenCL CONFIG REQUIRED)

add_executable(parallel_single_precision_wind_induced_currents_davies_85 
    main.cpp
    solver.cpp
    kernel.cpp
)

target_copy_opencl(parallel_single_precision_wind_induced_currents_davies_85)
target_copy_bathymetry_data(parallel_single_precision_wind_induced_currents_davies_85)
target_copy_wind_data(parallel_single_precision_wind_induced_currents_davies_85)

target_finalize_asset_copy(parallel_single_precision_wind_induced_currents_davies_85)

target_add_manifest(parallel_single_precision_wind_induced_currents_davies_85)

target_link_libraries(parallel_single_precision_wind_induced_currents_davies_85 PRIVATE
    utils
    slae
    calc
    generators
    writers
    OpenCL::OpenCL
)
