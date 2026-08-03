if(SIPROS_BUILD_MPI)
  find_package(MPI REQUIRED COMPONENTS CXX)
endif()

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

if(SIPROS_BUILD_MPI)
  add_executable(siprosMPI ${PROJECT_SOURCE_DIR}/sipros/src/mpi_main.cpp)
  target_include_directories(siprosMPI PUBLIC ${PROJECT_SOURCE_DIR}/sipros/include)
  target_link_libraries(siprosMPI PUBLIC sipros_common sipros_core MPI::MPI_CXX OpenMP::OpenMP_CXX)
  target_compile_options(siprosMPI PUBLIC -ffast-math -march=native -fPIC)
  sipros_compile_definitions(siprosMPI)
endif()

if(BUILD_TESTING)
  add_executable(builtin_config_test
    ${PROJECT_SOURCE_DIR}/sipros/tests/test_builtin_config.cpp
  )
  target_include_directories(builtin_config_test PUBLIC ${SIPROS_INCLUDE_DIRS})
  target_link_libraries(builtin_config_test PUBLIC sipros_core)
  sipros_compile_definitions(builtin_config_test)
  sipros_optimization_options(builtin_config_test)
  add_test(
    NAME builtin_config
    COMMAND builtin_config_test
  )

  add_executable(ms1_isotope_envelope_test
    ${PROJECT_SOURCE_DIR}/sipros/tests/test_ms1_isotope_envelope.cpp
  )
  target_include_directories(ms1_isotope_envelope_test PUBLIC ${SIPROS_INCLUDE_DIRS})
  target_link_libraries(ms1_isotope_envelope_test PUBLIC sipros_common sipros_core)
  sipros_compile_definitions(ms1_isotope_envelope_test)
  sipros_optimization_options(ms1_isotope_envelope_test)
  add_test(
    NAME ms1_isotope_envelope
    COMMAND ms1_isotope_envelope_test
  )

  add_executable(pin_writer_schema_test
    ${PROJECT_SOURCE_DIR}/sipros/tests/test_pin_writer_schema.cpp
  )
  target_include_directories(pin_writer_schema_test PUBLIC ${SIPROS_INCLUDE_DIRS})
  target_link_libraries(pin_writer_schema_test PUBLIC sipros_core)
  sipros_compile_definitions(pin_writer_schema_test)
  sipros_optimization_options(pin_writer_schema_test)
  add_test(
    NAME pin_writer_schema
    COMMAND pin_writer_schema_test
  )

  add_executable(fragment_index_test
    ${PROJECT_SOURCE_DIR}/sipros/tests/test_fragment_index.cpp
  )
  target_include_directories(fragment_index_test PUBLIC ${SIPROS_INCLUDE_DIRS})
  target_link_libraries(fragment_index_test PUBLIC sipros_core)
  sipros_compile_definitions(fragment_index_test)
  sipros_optimization_options(fragment_index_test)
  add_test(
    NAME fragment_index
    COMMAND fragment_index_test
  )

  add_executable(spectra_index_test
    ${PROJECT_SOURCE_DIR}/sipros/tests/test_spectra_index.cpp
  )
  target_include_directories(spectra_index_test PUBLIC ${SIPROS_INCLUDE_DIRS})
  target_link_libraries(spectra_index_test PUBLIC sipros_core)
  sipros_compile_definitions(spectra_index_test)
  sipros_optimization_options(spectra_index_test)
  add_test(
    NAME spectra_index
    COMMAND spectra_index_test
  )

  add_executable(native_precursor_test
    ${PROJECT_SOURCE_DIR}/sipros/tests/test_native_precursors.cpp
  )
  target_include_directories(native_precursor_test PUBLIC ${SIPROS_INCLUDE_DIRS})
  target_link_libraries(native_precursor_test PUBLIC sipros_core)
  sipros_compile_definitions(native_precursor_test)
  sipros_optimization_options(native_precursor_test)
  add_test(
    NAME native_precursor
    COMMAND native_precursor_test
  )

  add_executable(protein_n_terminal_peptides_test
    ${PROJECT_SOURCE_DIR}/sipros/tests/test_protein_n_terminal_peptides.cpp
  )
  target_include_directories(protein_n_terminal_peptides_test PUBLIC ${SIPROS_INCLUDE_DIRS})
  target_link_libraries(protein_n_terminal_peptides_test PUBLIC sipros_core)
  sipros_compile_definitions(protein_n_terminal_peptides_test)
  sipros_optimization_options(protein_n_terminal_peptides_test)
  add_test(
    NAME protein_n_terminal_peptides
    COMMAND protein_n_terminal_peptides_test
  )

  add_executable(search_policy_test
    ${PROJECT_SOURCE_DIR}/sipros/tests/test_search_policy.cpp
  )
  target_include_directories(search_policy_test PUBLIC ${SIPROS_INCLUDE_DIRS})
  target_link_libraries(search_policy_test PUBLIC sipros_common sipros_core)
  sipros_compile_definitions(search_policy_test)
  sipros_optimization_options(search_policy_test)
  add_test(
    NAME search_policy
    COMMAND search_policy_test
  )
endif()
