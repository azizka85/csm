cmake_minimum_required(VERSION 3.16)

project(generators LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package (GeographicLib CONFIG REQUIRED)
find_package(netCDFCxx CONFIG REQUIRED)
find_package(boost_math CONFIG REQUIRED)

add_library(generators STATIC
    utils.cpp

    area/uniform_generator.cpp
    area/proj_generator.cpp

    bathymetry/uniform_generator.cpp
    bathymetry/cosine_generator.cpp
    bathymetry/gebco_generator

    dz/uniform_generator.cpp
    dz/triple_point_generator.cpp

    viscosity/uniform_generator.cpp
    viscosity/linear_generator.cpp
    viscosity/wind_speed_generator.cpp

    wind/uniform_generator.cpp
    wind/ecmwf_generator.cpp
)

target_include_directories(generators PUBLIC ${CMAKE_SOURCE_DIR})

target_link_libraries(generators PRIVATE
    calc
    utils    
    netCDF::netcdf netCDF::netcdf-cxx4
    Boost::math
    ${GeographicLib_LIBRARIES}    
)