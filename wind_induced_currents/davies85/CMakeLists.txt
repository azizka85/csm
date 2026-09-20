cmake_minimum_required(VERSION 3.16)

project(wind_induced_currents_davies_85 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(wind_induced_currents_davies_85 
    main.cpp
    solver.cpp   
)

target_copy_bathymetry_data(wind_induced_currents_davies_85)
target_copy_wind_data(wind_induced_currents_davies_85)

target_finalize_asset_copy(wind_induced_currents_davies_85)

target_add_manifest(wind_induced_currents_davies_85)

target_link_libraries(wind_induced_currents_davies_85 PRIVATE
    utils
    slae
    calc
    generators
    writers
)
