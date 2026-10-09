cmake_minimum_required(VERSION 3.28)

project(parallel_slae LANGUAGES CXX CUDA)

set(CMAKE_CUDA_ARCHITECTURES 100)

add_compile_options($<$<COMPILE_LANGUAGE:CUDA>:-allow-unsupported-compiler>)

find_package(CUDAToolkit REQUIRED)

add_library(parallel_slae STATIC 
	direct/tridiagonal.cu
)

target_include_directories(parallel_slae PRIVATE ${CUDAToolkit_INCLUDE_DIRS})

target_include_directories(parallel_slae PUBLIC ${CMAKE_SOURCE_DIR})

set_target_properties(parallel_slae PROPERTIES
	CXX_STANDARD 20
    CUDA_STANDARD 20
    CUDA_STANDARD_REQUIRED ON
    CMAKE_CXX_STANDARD_REQUIRED ON
    CUDA_SEPARABLE_COMPILATION ON
)
