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
  src/statistics.cpp
  src/svm.cpp
)
target_include_directories(aerith_filter PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
if(NOT TARGET OpenMP::OpenMP_CXX)
  find_package(OpenMP REQUIRED)
endif()
target_link_libraries(aerith_filter PRIVATE OpenMP::OpenMP_CXX)
if(MSVC)
  target_compile_options(aerith_filter PRIVATE /W4 /O2)
else()
  target_compile_options(aerith_filter PRIVATE -Wall -Wextra -Wpedantic -O3 -ftree-vectorize)
endif()

add_executable(aerith src/main.cpp)
target_link_libraries(aerith PRIVATE aerith_filter)

include(CTest)
if(BUILD_TESTING)
  add_executable(aerith_tests tests/test_filter.cpp)
  target_link_libraries(aerith_tests PRIVATE aerith_filter)
  add_test(NAME aerith_unit_tests COMMAND aerith_tests)
endif()
