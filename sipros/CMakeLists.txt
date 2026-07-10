find_package(MPI REQUIRED COMPONENTS CXX)

set(SIPROS_COMMON_SOURCE ${PROJECT_SOURCE_DIR}/sipros/src/SiprosSearchRunner.cpp)
set(SIPROS_WORKFLOW_SOURCES
  ${PROJECT_SOURCE_DIR}/sipros/src/DatabaseSearchWorkflow.cpp
  ${PROJECT_SOURCE_DIR}/sipros/src/TheoreticalSpectraWorkflow.cpp
  ${PROJECT_SOURCE_DIR}/sipros/src/ExperimentalSpectraWorkflow.cpp
  ${PROJECT_SOURCE_DIR}/sipros/src/SearchSpectraWorkflow.cpp
)
set(SIPROS_ENTRYPOINT_SOURCES
  ${PROJECT_SOURCE_DIR}/sipros/src/omp_main.cpp
  ${PROJECT_SOURCE_DIR}/sipros/src/mpi_main.cpp
)

file(GLOB SIPROS_ALL_SOURCES CONFIGURE_DEPENDS ${PROJECT_SOURCE_DIR}/sipros/src/*.cpp)
set(SIPROS_CORE_SOURCES ${SIPROS_ALL_SOURCES})
list(REMOVE_ITEM SIPROS_CORE_SOURCES
  ${SIPROS_COMMON_SOURCE}
  ${SIPROS_WORKFLOW_SOURCES}
  ${SIPROS_ENTRYPOINT_SOURCES}
)

add_library(sipros_core STATIC ${SIPROS_CORE_SOURCES})
target_include_directories(sipros_core PUBLIC
  ${SIPROS_INCLUDE_DIRS}
  ${HDF5_INCLUDE_DIRS}
)
target_link_libraries(sipros_core PUBLIC OpenMP::OpenMP_CXX ${SIPROS_HDF5_LIBRARIES})
sipros_compile_definitions(sipros_core)
sipros_optimization_options(sipros_core)

add_library(sipros_common STATIC ${SIPROS_COMMON_SOURCE})
target_include_directories(sipros_common PUBLIC ${SIPROS_INCLUDE_DIRS})
target_link_libraries(sipros_common PUBLIC sipros_core)
sipros_compile_definitions(sipros_common)
sipros_optimization_options(sipros_common)

add_library(sipros_openmp_tools STATIC ${SIPROS_WORKFLOW_SOURCES})
target_include_directories(sipros_openmp_tools PUBLIC
  ${PROJECT_SOURCE_DIR}/sipros/include
  ${HDF5_INCLUDE_DIRS}
)
target_link_libraries(sipros_openmp_tools PUBLIC
  sipros_common
  sipros_core
  ${SIPROS_HDF5_LIBRARIES}
)
sipros_compile_definitions(sipros_openmp_tools)
sipros_optimization_options(sipros_openmp_tools)

add_executable(sipros ${PROJECT_SOURCE_DIR}/sipros/src/omp_main.cpp)
target_include_directories(sipros PUBLIC ${PROJECT_SOURCE_DIR}/sipros/include)
target_link_libraries(sipros PUBLIC sipros_openmp_tools)
sipros_compile_definitions(sipros)
sipros_optimization_options(sipros)

add_executable(siprosMPI ${PROJECT_SOURCE_DIR}/sipros/src/mpi_main.cpp)
target_include_directories(siprosMPI PUBLIC ${PROJECT_SOURCE_DIR}/sipros/include)
target_link_libraries(siprosMPI PUBLIC sipros_common sipros_core MPI::MPI_CXX OpenMP::OpenMP_CXX)
target_compile_options(siprosMPI PUBLIC -ffast-math -march=native -fPIC)
sipros_compile_definitions(siprosMPI)

if(BUILD_TESTING)
  add_executable(ms1_isotope_envelope_test
    ${PROJECT_SOURCE_DIR}/sipros/tests/test_ms1_isotope_envelope.cpp
  )
  target_include_directories(ms1_isotope_envelope_test PUBLIC ${SIPROS_INCLUDE_DIRS})
  target_link_libraries(ms1_isotope_envelope_test PUBLIC sipros_common sipros_core)
  sipros_compile_definitions(ms1_isotope_envelope_test)
  sipros_optimization_options(ms1_isotope_envelope_test)
  add_test(
    NAME ms1_isotope_envelope
    COMMAND ms1_isotope_envelope_test ${PROJECT_SOURCE_DIR}/configTemplates/SIP.cfg
  )
endif()
