cmake_minimum_required(VERSION 3.28)

project(parallel_wind_induced_currents_davies_85 LANGUAGES CXX CUDA)

set(CMAKE_CUDA_ARCHITECTURES 100)

add_compile_options($<$<COMPILE_LANGUAGE:CUDA>:-allow-unsupported-compiler>)

find_package(CUDAToolkit REQUIRED)
find_package(Thrust REQUIRED)

add_executable(parallel_wind_induced_currents_davies_85 main.cpp solver.cpp kernel.cu)

target_copy_bathymetry_data(parallel_wind_induced_currents_davies_85)
target_copy_wind_data(parallel_wind_induced_currents_davies_85)

target_finalize_asset_copy(parallel_wind_induced_currents_davies_85)

target_add_manifest(parallel_wind_induced_currents_davies_85)

target_include_directories(parallel_wind_induced_currents_davies_85 PRIVATE ${CUDAToolkit_INCLUDE_DIRS})

target_link_libraries(parallel_wind_induced_currents_davies_85 PRIVATE 
    utils
    calc
    generators
    writers
    parallel_slae
    CUDA::cudart
    CUDA::cublas
)

set_target_properties(parallel_wind_induced_currents_davies_85 PROPERTIES
    CXX_STANDARD 20
    CUDA_STANDARD 20
    CUDA_STANDARD_REQUIRED ON
    CMAKE_CXX_STANDARD_REQUIRED ON
    CUDA_SEPARABLE_COMPILATION ON
)