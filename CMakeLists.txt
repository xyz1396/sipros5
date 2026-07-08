cmake_minimum_required(VERSION 3.12)
project(Sipros)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

message("${CMAKE_BUILD_TYPE}")
message("${CMAKE_CXX_COMPILER}")

set(CMAKE_CXX_STANDARD 17)
find_package(OpenMP REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS C CXX)

if(TARGET hdf5::hdf5-static AND TARGET hdf5::hdf5_cpp-static)
  set(SIPROS_HDF5_LIBRARIES hdf5::hdf5-static hdf5::hdf5_cpp-static)
elseif(TARGET hdf5::hdf5-shared AND TARGET hdf5::hdf5_cpp-shared)
  set(SIPROS_HDF5_LIBRARIES hdf5::hdf5-shared hdf5::hdf5_cpp-shared)
elseif(TARGET hdf5-static AND TARGET hdf5_cpp-static)
  set(SIPROS_HDF5_LIBRARIES hdf5-static hdf5_cpp-static)
elseif(TARGET hdf5-shared AND TARGET hdf5_cpp-shared)
  set(SIPROS_HDF5_LIBRARIES hdf5-shared hdf5_cpp-shared)
elseif(TARGET HDF5::HDF5)
  set(SIPROS_HDF5_LIBRARIES HDF5::HDF5)
else()
  set(SIPROS_HDF5_LIBRARIES ${HDF5_LIBRARIES})
endif()

set(SIPROS_INCLUDE_DIRS
  ${PROJECT_SOURCE_DIR}/sipros/include
)
include_directories(${SIPROS_INCLUDE_DIRS})

SET(EXECUTABLE_OUTPUT_PATH "${PROJECT_SOURCE_DIR}/bin")

function(sipros_compile_definitions target_name)
  if(WIN32)
    target_compile_definitions(${target_name} PUBLIC _NO_THERMORAW _NOSQLITE)
  else()
    target_compile_definitions(${target_name} PUBLIC GCC _FILE_OFFSET_BITS=64 _NOSQLITE)
  endif()
endfunction()

function(sipros_optimization_options target_name)
  if(MSVC)
    target_compile_options(${target_name} PUBLIC "/fp:fast" "/openmp:experimental" "/openmp:llvm" "/arch:AVX2")
  elseif(BUILD_CONDA)
    target_compile_options(${target_name} PUBLIC -ffast-math -march=native)
  else()
    target_compile_options(${target_name} PUBLIC -ffast-math)
  endif()
endfunction()

add_subdirectory(sipros)

if(NOT MSVC AND NOT BUILD_CONDA)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  message(CMAKE_CXX_FLAGS: "${CMAKE_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
  message(CMAKE_EXE_LINKER_FLAGS: "${CMAKE_EXE_LINKER_FLAGS}")
endif()

if(MSVC)
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
  set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} /MT")
endif()

