# Add executable
file(GLOB SOURCES_Aerith "src/*.cpp")
add_executable(aerithFeatureExtractor ${SOURCES_Aerith} ${PROJECT_SOURCE_DIR}/src/averagine.cpp 
    ${PROJECT_SOURCE_DIR}/src/proNovoConfig.cpp
    ${PROJECT_SOURCE_DIR}/src/isotopologue.cpp)

# Include directories
include_directories(include)

IF(MSVC)
  target_link_libraries(aerithFeatureExtractor PUBLIC OpenMP::OpenMP_CXX)
  target_compile_options(aerithFeatureExtractor PUBLIC "/fp:fast" "/openmp:experimental" "/openmp:llvm" "/arch:AVX2")
ELSE()
  IF(BUILD_CONDA)
    target_link_libraries(aerithFeatureExtractor PUBLIC OpenMP::OpenMP_CXX)
    target_compile_options(aerithFeatureExtractor PUBLIC -ffast-math -march=native)
    # target_compile_options(aerithFeatureExtractor PUBLIC -ffast-math)
  ELSE()
    # this flag is not working for static link uper than gcc 11
    # target_link_libraries(aerithFeatureExtractor PUBLIC mstoolkit OpenMP::OpenMP_CXX -static -march=native)
    # target_compile_options(aerithFeatureExtractor PUBLIC -ffast-math -march=native)
    target_compile_options(aerithFeatureExtractor PUBLIC -ffast-math)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    message(CMAKE_CXX_FLAGS: "${CMAKE_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
    message(CMAKE_EXE_LINKER_FLAGS: "${CMAKE_EXE_LINKER_FLAGS}")
  ENDIF()
ENDIF()