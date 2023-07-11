project(SiprosEnsemble)
cmake_minimum_required(VERSION 3.1)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

message("${CMAKE_BUILD_TYPE}")

# set(CMAKE_C_COMPILER "gcc-11")
# set(CMAKE_CXX_COMPILER "g++-11")
message("${CMAKE_CXX_COMPILER}")
set(CMAKE_CXX_STANDARD 17)
find_package(OpenMP REQUIRED)

aux_source_directory(${PROJECT_SOURCE_DIR}/src SRC_LIST)
include_directories(${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/MSToolkit/include)
include_directories(${PROJECT_SOURCE_DIR}/MSToolkit/src/zlib-1.2.11/include)
include_directories(${PROJECT_SOURCE_DIR}/MSToolkit/src/expat-2.2.9/include)

SET(EXECUTABLE_OUTPUT_PATH "${PROJECT_SOURCE_DIR}/bin")

add_subdirectory(MSToolkit)
add_subdirectory(openmp)
add_subdirectory(mpi)

