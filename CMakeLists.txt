cmake_minimum_required(VERSION 3.12)
project(Sipros)
include(CTest)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

message("${CMAKE_BUILD_TYPE}")
message("${CMAKE_CXX_COMPILER}")

set(CMAKE_CXX_STANDARD 17)

# Release and Conda builds dynamically link their external runtime libraries.
option(BUILD_CONDA "Build with the Conda compiler suite" OFF)
option(SIPROS_BUILD_MPI "Build the siprosMPI executable" ON)

find_package(OpenMP REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS C CXX)

if(TARGET hdf5::hdf5-shared AND TARGET hdf5::hdf5_cpp-shared)
  set(SIPROS_HDF5_LIBRARIES hdf5::hdf5-shared hdf5::hdf5_cpp-shared)
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

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

function(sipros_optimization_options target_name)
  if(MSVC)
    target_compile_options(${target_name} PUBLIC
      "/fp:fast" "/arch:AVX2" "/utf-8")
  elseif(BUILD_CONDA)
    target_compile_options(${target_name} PUBLIC -ffast-math -march=native)
  else()
    target_compile_options(${target_name} PUBLIC -ffast-math)
  endif()
endfunction()

function(sipros_windows_conda_openmp target_name)
  if(NOT (MSVC AND BUILD_CONDA))
    return()
  endif()

  set(conda_openmp_import "$ENV{CONDA_PREFIX}/Library/lib/libomp.lib")
  if(NOT EXISTS "${conda_openmp_import}")
    message(FATAL_ERROR
      "The Conda OpenMP import library is required for Windows packaging: ${conda_openmp_import}")
  endif()
  # /openmp:llvm embeds /DEFAULTLIB:libomp.lib in each OpenMP object. Ignore
  # that implicit MSVC library and explicitly use Conda's redistributable
  # import library, whose runtime name is libomp.dll.
  target_link_options(${target_name} PRIVATE "/NODEFAULTLIB:libomp.lib")
  target_link_libraries(${target_name} PRIVATE "${conda_openmp_import}")
endfunction()

add_subdirectory(sipros)
add_subdirectory(aerith)
add_subdirectory(wf33)

if(MSVC)
  # Conda's HDF5 and PyTorch packages use the dynamic MSVC runtime. Mixing
  # /MT objects with those /MD libraries causes allocator and ABI failures.
  set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")
endif()
