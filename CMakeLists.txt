cmake_minimum_required(VERSION 3.12)
project(Sipros)
include(CTest)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

message("${CMAKE_BUILD_TYPE}")
message("${CMAKE_CXX_COMPILER}")

set(CMAKE_CXX_STANDARD 17)

# The system-toolchain build publishes a portable, fully static non-MPI
# executable.  The Conda build links all runtime libraries dynamically.
option(BUILD_CONDA "Build with the Conda compiler suite" OFF)
option(SIPROS_STATIC_HDF5 "Require static HDF5 C and C++ libraries" OFF)

# FindOpenMP otherwise prefers libgomp.so even when the final link uses
# -static, which makes the link fail with "attempted static link of dynamic
# object".  Select GCC's matching static OpenMP runtime up front.
set(SIPROS_FULLY_STATIC OFF)
if(NOT MSVC AND NOT BUILD_CONDA)
  set(SIPROS_FULLY_STATIC ON)
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    execute_process(
      COMMAND ${CMAKE_CXX_COMPILER} -print-file-name=libgomp.a
      OUTPUT_VARIABLE SIPROS_STATIC_LIBGOMP
      OUTPUT_STRIP_TRAILING_WHITESPACE
      RESULT_VARIABLE SIPROS_STATIC_LIBGOMP_RESULT
    )
    if(NOT SIPROS_STATIC_LIBGOMP_RESULT EQUAL 0 OR
       NOT IS_ABSOLUTE "${SIPROS_STATIC_LIBGOMP}" OR
       NOT EXISTS "${SIPROS_STATIC_LIBGOMP}")
      message(FATAL_ERROR "A static libgomp.a is required for the static Sipros build")
    endif()
    set(OpenMP_gomp_LIBRARY "${SIPROS_STATIC_LIBGOMP}" CACHE FILEPATH
      "Static GNU OpenMP runtime used by Sipros" FORCE)
  endif()
endif()

find_package(OpenMP REQUIRED)
if(SIPROS_STATIC_HDF5)
  # The standalone vcpkg HDF5 config exports these dependency names but does
  # not load them itself when it is consumed without the vcpkg toolchain.
  find_package(Threads REQUIRED)
  find_package(ZLIB CONFIG REQUIRED)
  find_package(libaec CONFIG REQUIRED)
  find_package(HDF5 CONFIG REQUIRED COMPONENTS C CXX static)
else()
  find_package(HDF5 REQUIRED COMPONENTS C CXX)
endif()

if(SIPROS_STATIC_HDF5 AND
   NOT ((TARGET hdf5::hdf5-static AND TARGET hdf5::hdf5_cpp-static) OR
        (TARGET hdf5-static AND TARGET hdf5_cpp-static)))
  message(FATAL_ERROR
    "SIPROS_STATIC_HDF5=ON, but the selected HDF5 package does not provide static C and C++ targets")
endif()

if(BUILD_CONDA AND TARGET hdf5::hdf5-shared AND TARGET hdf5::hdf5_cpp-shared)
  set(SIPROS_HDF5_LIBRARIES hdf5::hdf5-shared hdf5::hdf5_cpp-shared)
elseif(BUILD_CONDA AND TARGET hdf5-shared AND TARGET hdf5_cpp-shared)
  set(SIPROS_HDF5_LIBRARIES hdf5-shared hdf5_cpp-shared)
elseif(TARGET hdf5::hdf5-static AND TARGET hdf5::hdf5_cpp-static)
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

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

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

function(sipros_static_link_options target_name)
  if(SIPROS_FULLY_STATIC)
    # LINK_FLAGS is used instead of target_link_options to retain the declared
    # CMake 3.12 compatibility.
    set_property(TARGET ${target_name} APPEND_STRING PROPERTY LINK_FLAGS " -static")
  endif()
endfunction()

add_subdirectory(sipros)
add_subdirectory(aerith)

if(MSVC)
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
  set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} /MT")
endif()
