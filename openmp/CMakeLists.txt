# set(CMAKE_C_COMPILER "gcc-11")
# set(CMAKE_CXX_COMPILER "g++-11")

add_executable(Sipros ${SRC_LIST} ${PROJECT_SOURCE_DIR}/openmp/main.cpp)

IF(MSVC)
  find_package(OpenMP REQUIRED)
  target_compile_options(Sipros PUBLIC "/fp:fast" "/openmp:experimental" "/openmp:llvm" "/arch:AVX2")
  # make sure to use the static runtime library for release build
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
  set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} /MT")
  target_link_libraries(Sipros PUBLIC mstoolkit OpenMP::OpenMP_CXX)
ELSE()
  target_compile_options(Sipros PUBLIC -ffast-math -march=native)
# target_compile_options(Sipros PUBLIC -ffast-math)
  target_link_libraries(Sipros PUBLIC mstoolkit -fopenmp -static)
ENDIF()

# for conda only have libgomp.so
# target_link_libraries(Sipros PUBLIC OpenMP::OpenMP_CXX MStoolkit gomp dl)

IF(WIN32)
  target_compile_definitions(Sipros PUBLIC _NO_THERMORAW)
ElSE()
  target_compile_definitions(Sipros PUBLIC GCC _FILE_OFFSET_BITS=64)
ENDIF()
