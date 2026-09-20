cmake_minimum_required(VERSION 3.16)

project(tests LANGUAGES CXX)

find_package(doctest CONFIG REQUIRED)

add_executable(tests 
    test_main.cpp

    slae/direct/test_tridiagonal.cpp

    calc/test_grid.cpp
    calc/test_time.cpp
    calc/test_turbulence.cpp
    calc/test_wind.cpp

    generators/test_utils.cpp
    
    generators/area/test_uniform_generator.cpp
    generators/area/test_proj_generator.cpp

    generators/bathymetry/test_uniform_generator.cpp
    generators/bathymetry/test_gebco_generator.cpp

    generators/dz/test_uniform_generator.cpp
    generators/dz/test_triple_point_generator.cpp

    generators/viscosity/test_uniform_generator.cpp
    generators/viscosity/test_linear_generator.cpp
    generators/viscosity/test_wind_speed_generator.cpp

    generators/wind/test_uniform_generator.cpp
)

target_copy_bathymetry_data(tests)

target_finalize_asset_copy(tests)

target_add_manifest(tests)

target_link_libraries(tests 
    PRIVATE 
        slae
        calc
        generators
        doctest::doctest
)