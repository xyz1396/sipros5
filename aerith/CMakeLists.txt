cmake_minimum_required(VERSION 3.12)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  project(aerith LANGUAGES CXX)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(aerith_filter
  src/filter.cpp
  src/output.cpp
  src/pin.cpp
  src/rt.cpp
  src/spectrum_entropy.cpp
  src/statistics.cpp
  src/svm.cpp
)
target_include_directories(aerith_filter PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
if(NOT TARGET OpenMP::OpenMP_CXX)
  find_package(OpenMP REQUIRED)
endif()
target_link_libraries(aerith_filter PRIVATE OpenMP::OpenMP_CXX)

if(BUILD_CONDA)
  set(AERITH_TORCH_DEFAULT ON)
else()
  set(AERITH_TORCH_DEFAULT OFF)
endif()
option(AERITH_ENABLE_TORCH "Enable native DIA-NN spectrum entropy prediction"
       ${AERITH_TORCH_DEFAULT})
if(AERITH_ENABLE_TORCH)
  find_package(Torch REQUIRED)
endif()
if(Torch_FOUND)
  if(NOT HDF5_FOUND)
    find_package(HDF5 REQUIRED COMPONENTS C CXX)
  endif()
  target_compile_definitions(aerith_filter PRIVATE AERITH_WITH_TORCH=1)
  target_include_directories(aerith_filter PRIVATE ${HDF5_INCLUDE_DIRS})
  target_link_libraries(aerith_filter PRIVATE ${TORCH_LIBRARIES})
  if(DEFINED SIPROS_HDF5_LIBRARIES)
    target_link_libraries(aerith_filter PRIVATE ${SIPROS_HDF5_LIBRARIES})
  elseif(TARGET HDF5::HDF5)
    target_link_libraries(aerith_filter PRIVATE HDF5::HDF5)
  else()
    target_link_libraries(aerith_filter PRIVATE ${HDF5_LIBRARIES})
  endif()
  if(TORCH_CXX_FLAGS)
    separate_arguments(AERITH_TORCH_CXX_FLAGS NATIVE_COMMAND "${TORCH_CXX_FLAGS}")
    target_compile_options(aerith_filter PRIVATE ${AERITH_TORCH_CXX_FLAGS})
  endif()
elseif(AERITH_ENABLE_TORCH)
  message(FATAL_ERROR "AERITH_ENABLE_TORCH requires a dynamic LibTorch installation")
else()
  message(STATUS "LibTorch not found: Aerith will build without --spectra support")
endif()
if(MSVC)
  target_compile_options(aerith_filter PRIVATE /W4 /O2)
else()
  target_compile_options(aerith_filter PRIVATE -Wall -Wextra -Wpedantic -O3 -ftree-vectorize)
endif()

add_executable(aerith src/main.cpp)
target_link_libraries(aerith PRIVATE aerith_filter)
if(Torch_FOUND)
  get_filename_component(AERITH_DIANN_MODEL
    "${CMAKE_CURRENT_SOURCE_DIR}/../tools/diann-2.6.1-fragmentation.pt" ABSOLUTE)
  if(NOT EXISTS "${AERITH_DIANN_MODEL}")
    message(FATAL_ERROR "Missing DIA-NN fragmentation model: ${AERITH_DIANN_MODEL}")
  endif()
  add_custom_command(TARGET aerith POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${AERITH_DIANN_MODEL}"
            "$<TARGET_FILE_DIR:aerith>/diann-2.6.1-fragmentation.pt")
endif()

include(CTest)
if(BUILD_TESTING)
  add_executable(aerith_tests tests/test_filter.cpp)
  target_link_libraries(aerith_tests PRIVATE aerith_filter)
  add_test(NAME aerith_unit_tests COMMAND aerith_tests)
endif()
