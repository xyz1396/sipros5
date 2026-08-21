cmake_minimum_required(VERSION 3.12)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  project(aerith LANGUAGES CXX)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(MSVC)
  # CMake injects /W3 into its MSVC defaults. Aerith selects /W4 below, so
  # remove the inherited level in this directory instead of asking cl.exe to
  # override it for every translation unit.
  string(REGEX REPLACE "(^| )[/-]W[0-4]( |$)" " "
    CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
endif()

add_library(aerith_filter
  src/filter.cpp
  src/isotope.cpp
  src/output.cpp
  src/pin.cpp
  src/prediction_cache.cpp
  src/protein_assembly.cpp
  src/predicted_rt.cpp
  src/quantification.cpp
  src/rt.cpp
  src/spectrum_entropy.cpp
  src/statistics.cpp
  src/svm.cpp
)
target_include_directories(aerith_filter PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
if(TARGET sipros_core)
  sipros_optimization_options(aerith_filter)
  target_include_directories(aerith_filter PRIVATE
    ${PROJECT_SOURCE_DIR}/sipros/include)
  target_link_libraries(aerith_filter PRIVATE sipros_core)
else()
  message(FATAL_ERROR
    "Aerith requires the Sipros core isotope-chemistry target")
endif()
if(NOT TARGET OpenMP::OpenMP_CXX)
  find_package(OpenMP REQUIRED)
endif()
target_link_libraries(aerith_filter PRIVATE OpenMP::OpenMP_CXX)

if(BUILD_CONDA)
  set(AERITH_TORCH_DEFAULT ON)
else()
  set(AERITH_TORCH_DEFAULT OFF)
endif()
option(AERITH_ENABLE_TORCH "Enable native DIA-NN spectrum and RT prediction"
       ${AERITH_TORCH_DEFAULT})
option(AERITH_TORCH_CPU_ONLY
       "Link only the dynamic CPU LibTorch libraries" OFF)
if(AERITH_ENABLE_TORCH)
  if(AERITH_TORCH_CPU_ONLY)
    set(AERITH_TORCH_ROOT "" CACHE PATH
        "Root of the system CPU LibTorch installation")
    if(NOT AERITH_TORCH_ROOT)
      message(FATAL_ERROR
        "AERITH_TORCH_CPU_ONLY requires AERITH_TORCH_ROOT")
    endif()
    find_path(AERITH_TORCH_INCLUDE_DIR
      NAMES torch/script.h
      PATHS
        "${AERITH_TORCH_ROOT}/include"
        "${AERITH_TORCH_ROOT}/include/torch/csrc/api/include"
      NO_DEFAULT_PATH)
    find_path(AERITH_TORCH_API_INCLUDE_DIR
      NAMES torch/torch.h
      PATHS "${AERITH_TORCH_ROOT}/include/torch/csrc/api/include"
      NO_DEFAULT_PATH)
    find_library(AERITH_LIBTORCH_CPU_LIBRARY
      NAMES torch_cpu PATHS "${AERITH_TORCH_ROOT}/lib" NO_DEFAULT_PATH)
    find_library(AERITH_LIBC10_LIBRARY
      NAMES c10 PATHS "${AERITH_TORCH_ROOT}/lib" NO_DEFAULT_PATH)
    if(NOT AERITH_TORCH_INCLUDE_DIR OR
       NOT AERITH_TORCH_API_INCLUDE_DIR OR
       NOT AERITH_LIBTORCH_CPU_LIBRARY OR
       NOT AERITH_LIBC10_LIBRARY)
      message(FATAL_ERROR
        "Incomplete CPU LibTorch installation under ${AERITH_TORCH_ROOT}")
    endif()
    set(Torch_FOUND TRUE)
    set(AERITH_TORCH_LIBRARIES
      "-Wl,--no-as-needed"
      "${AERITH_LIBTORCH_CPU_LIBRARY}"
      "${AERITH_LIBC10_LIBRARY}"
      "-Wl,--as-needed")
    set(AERITH_TORCH_INCLUDE_DIRS
      "${AERITH_TORCH_INCLUDE_DIR}"
      "${AERITH_TORCH_API_INCLUDE_DIR}")
    message(STATUS
      "Aerith: using CPU-only release LibTorch from ${AERITH_TORCH_ROOT}")
  else()
    if(DEFINED AERITH_TORCH_PACKAGE_VERSION)
      message(STATUS
        "Aerith: configuring with PyTorch ${AERITH_TORCH_PACKAGE_VERSION}")
    endif()
    # Conda's compiler root-path policy can hide NVTX3 from a system CUDA
    # toolkit. Locate the header from either the active Conda environment or
    # the toolkit containing nvcc before loading TorchConfig.cmake.
    if(NOT nvtx3_dir)
      find_path(AERITH_NVTX3_INCLUDE_DIR
        NAMES nvtx3/nvToolsExt.h
        HINTS "$ENV{CONDA_PREFIX}/include")
      if(NOT AERITH_NVTX3_INCLUDE_DIR AND AERITH_ENABLE_CUDA)
        find_program(AERITH_NVCC_EXECUTABLE NAMES nvcc)
        if(AERITH_NVCC_EXECUTABLE)
          get_filename_component(AERITH_CUDA_BIN_DIR
            "${AERITH_NVCC_EXECUTABLE}" DIRECTORY)
          get_filename_component(AERITH_CUDA_ROOT
            "${AERITH_CUDA_BIN_DIR}" DIRECTORY)
          find_path(AERITH_NVTX3_INCLUDE_DIR
            NAMES nvtx3/nvToolsExt.h
            HINTS
              "${AERITH_CUDA_ROOT}/include"
              "${AERITH_CUDA_ROOT}/targets/x86_64-linux/include"
            NO_DEFAULT_PATH)
        endif()
      endif()
      if(AERITH_NVTX3_INCLUDE_DIR)
        set(nvtx3_dir "${AERITH_NVTX3_INCLUDE_DIR}" CACHE PATH
            "Directory containing the NVTX3 include directory" FORCE)
        message(STATUS "Aerith: found NVTX3 in ${nvtx3_dir}")
      endif()
    endif()
    # Some shared Conda libtorch builds enable Kineto in TorchConfig.cmake but
    # fold its required symbols into the shared Torch libraries instead of
    # packaging a standalone libkineto.  Upstream TorchConfig emits a warning
    # and appends nothing in that case.  Pre-seed its lookup with a removable
    # sentinel so this known optional packaging mismatch does not obscure real
    # CMake warnings.
    set(AERITH_KINETO_LINK_SENTINEL "")
    if(DEFINED Torch_DIR)
      get_filename_component(AERITH_TORCH_PREFIX
        "${Torch_DIR}/../../.." ABSOLUTE)
      find_library(AERITH_KINETO_LIBRARY
        NAMES kineto
        PATHS "${AERITH_TORCH_PREFIX}/lib")
      if(NOT AERITH_KINETO_LIBRARY)
        set(AERITH_KINETO_LINK_SENTINEL
          "AERITH_OPTIONAL_KINETO_NOT_PACKAGED")
        set(kineto_LIBRARY "${AERITH_KINETO_LINK_SENTINEL}")
      endif()
    endif()
    # Torch uses TORCH_CUDA_ARCH_LIST and emits a warning whenever CMake's
    # separate architecture variable is visible. Preserve the project value
    # around TorchConfig so both systems retain their intended settings.
    if(AERITH_ENABLE_CUDA AND NOT CMAKE_CUDA_COMPILER_LOADED)
      enable_language(CUDA)
    endif()
    set(AERITH_CUDA_ARCHITECTURES_DEFINED FALSE)
    if(DEFINED CMAKE_CUDA_ARCHITECTURES)
      set(AERITH_CUDA_ARCHITECTURES_DEFINED TRUE)
      set(AERITH_CUDA_ARCHITECTURES "${CMAKE_CUDA_ARCHITECTURES}")
    endif()
    get_property(AERITH_CUDA_ARCHITECTURES_CACHE_TYPE CACHE
      CMAKE_CUDA_ARCHITECTURES PROPERTY TYPE)
    if(AERITH_CUDA_ARCHITECTURES_CACHE_TYPE)
      get_property(AERITH_CUDA_ARCHITECTURES_CACHE_VALUE CACHE
        CMAKE_CUDA_ARCHITECTURES PROPERTY VALUE)
      unset(CMAKE_CUDA_ARCHITECTURES CACHE)
    endif()
    unset(CMAKE_CUDA_ARCHITECTURES)
    find_package(Torch 2.12.1 REQUIRED)
    unset(CMAKE_CUDA_ARCHITECTURES)
    if(AERITH_CUDA_ARCHITECTURES_CACHE_TYPE)
      set(CMAKE_CUDA_ARCHITECTURES
        "${AERITH_CUDA_ARCHITECTURES_CACHE_VALUE}"
        CACHE ${AERITH_CUDA_ARCHITECTURES_CACHE_TYPE}
        "CUDA architectures selected outside TorchConfig" FORCE)
    elseif(AERITH_CUDA_ARCHITECTURES_DEFINED)
      set(CMAKE_CUDA_ARCHITECTURES "${AERITH_CUDA_ARCHITECTURES}")
    endif()
    if(AERITH_KINETO_LINK_SENTINEL)
      list(REMOVE_ITEM TORCH_LIBRARIES "${AERITH_KINETO_LINK_SENTINEL}")
      unset(kineto_LIBRARY)
      message(STATUS
        "Aerith: optional standalone Kineto library is not packaged; omitting it")
    endif()
  endif()
endif()
if(Torch_FOUND)
  if(NOT HDF5_FOUND)
    find_package(HDF5 REQUIRED COMPONENTS C CXX)
  endif()
  target_compile_definitions(aerith_filter PRIVATE AERITH_WITH_TORCH=1)
  if(AERITH_TORCH_CPU_ONLY)
    target_compile_definitions(aerith_filter PRIVATE
      AERITH_TORCH_CPU_ONLY=1 _GLIBCXX_USE_CXX11_ABI=1)
    target_include_directories(aerith_filter SYSTEM PRIVATE
      ${AERITH_TORCH_INCLUDE_DIRS})
  endif()
  target_include_directories(aerith_filter PRIVATE ${HDF5_INCLUDE_DIRS})
  if(AERITH_TORCH_CPU_ONLY)
    target_link_libraries(aerith_filter PRIVATE ${AERITH_TORCH_LIBRARIES})
  else()
    target_link_libraries(aerith_filter PRIVATE ${TORCH_LIBRARIES})
  endif()
  if(DEFINED SIPROS_HDF5_LIBRARIES)
    target_link_libraries(aerith_filter PRIVATE ${SIPROS_HDF5_LIBRARIES})
  elseif(TARGET HDF5::HDF5)
    target_link_libraries(aerith_filter PRIVATE HDF5::HDF5)
  else()
    target_link_libraries(aerith_filter PRIVATE ${HDF5_LIBRARIES})
  endif()
  if(TORCH_CXX_FLAGS)
    separate_arguments(AERITH_TORCH_CXX_FLAGS NATIVE_COMMAND "${TORCH_CXX_FLAGS}")
    if(MSVC)
      # Aerith owns its warning level; do not let Torch's /W3 conflict with /W4.
      list(FILTER AERITH_TORCH_CXX_FLAGS EXCLUDE REGEX "^/W[0-4]$")
    endif()
    target_compile_options(aerith_filter PRIVATE ${AERITH_TORCH_CXX_FLAGS})
  endif()
elseif(AERITH_ENABLE_TORCH)
  message(FATAL_ERROR "AERITH_ENABLE_TORCH requires a dynamic LibTorch installation")
else()
  message(STATUS "LibTorch not found: Aerith will build without predicted-spectrum or predicted-RT features")
endif()
if(MSVC)
  target_compile_definitions(aerith_filter PRIVATE _CRT_SECURE_NO_WARNINGS)
  target_compile_options(aerith_filter PRIVATE
    /W4 /O2 /utf-8 /external:anglebrackets /external:W0)
  # These translation units instantiate PyTorch templates that emit narrowing
  # and unreachable-code diagnostics from MSVC/STL internals. Keep the same
  # warnings enabled for Aerith's own non-template code.
  set_source_files_properties(
    src/predicted_rt.cpp src/spectrum_entropy.cpp
    PROPERTIES COMPILE_OPTIONS "/wd4267;/wd4702")
else()
  target_compile_options(aerith_filter PRIVATE -Wall -Wextra -Wpedantic -O3 -ftree-vectorize)
endif()

add_executable(aerith src/main.cpp)
target_link_libraries(aerith PRIVATE aerith_filter)
sipros_windows_conda_openmp(aerith)
if(AERITH_TORCH_CPU_ONLY AND NOT MSVC)
  # The packaged CPU release resolves its complete non-glibc runtime from the
  # adjacent lib directory.
  set_target_properties(aerith PROPERTIES
    BUILD_WITH_INSTALL_RPATH TRUE
    INSTALL_RPATH "\$ORIGIN/lib")
endif()
if(Torch_FOUND)
  get_filename_component(AERITH_DIANN_MODEL
    "${CMAKE_CURRENT_SOURCE_DIR}/../tools/diann-2.6.1-fragmentation.pt" ABSOLUTE)
  if(NOT EXISTS "${AERITH_DIANN_MODEL}")
    message(FATAL_ERROR "Missing DIA-NN fragmentation model: ${AERITH_DIANN_MODEL}")
  endif()
  get_filename_component(AERITH_DIANN_RT_MODEL
    "${CMAKE_CURRENT_SOURCE_DIR}/../tools/diann-2.6.1-retention-time.pt" ABSOLUTE)
  if(NOT EXISTS "${AERITH_DIANN_RT_MODEL}")
    message(FATAL_ERROR "Missing DIA-NN RT model: ${AERITH_DIANN_RT_MODEL}")
  endif()
  add_custom_command(TARGET aerith POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${AERITH_DIANN_MODEL}"
            "$<TARGET_FILE_DIR:aerith>/diann-2.6.1-fragmentation.pt")
  add_custom_command(TARGET aerith POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${AERITH_DIANN_RT_MODEL}"
            "$<TARGET_FILE_DIR:aerith>/diann-2.6.1-retention-time.pt")
endif()

include(CTest)
if(BUILD_TESTING)
  add_executable(aerith_tests tests/test_filter.cpp)
  target_link_libraries(aerith_tests PRIVATE aerith_filter)
  if(Torch_FOUND)
    target_compile_definitions(aerith_tests PRIVATE AERITH_TEST_WITH_TORCH=1)
  endif()
  add_test(NAME aerith_unit_tests COMMAND aerith_tests)
endif()
