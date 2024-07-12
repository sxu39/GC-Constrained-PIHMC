cmake_minimum_required(VERSION 3.21.4)
project(CPIHMC)

add_compile_definitions(HIGH_PREC)

add_definitions(-w -O3)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../bin)

include_directories(include)

include_directories(/home/KawaSaki/software/deepmd-kit/2.2.1_cuda_10.2/include)

link_directories(/usr/local/lib64)

link_directories(/home/KawaSaki/software/deepmd-kit/2.2.1_cuda_10.2/lib)

aux_source_directory(source CPP_LIST)

add_executable(CPIHMC ${CPP_LIST})

target_link_libraries(CPIHMC deepmd_cc tensorflow_cc)
