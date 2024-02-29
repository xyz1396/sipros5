set(CMAKE_C_COMPILER "mpicc")
set(CMAKE_CXX_COMPILER "mpicxx")

add_executable(SiprosMPI ${SRC_LIST} ${PROJECT_SOURCE_DIR}/mpi/mpimain.cpp)
target_link_libraries(SiprosMPI PUBLIC mstoolkit -fopenmp)
target_compile_options(SiprosMPI PUBLIC -ffast-math -march=native -fPIC)
target_compile_definitions(SiprosMPI PUBLIC GCC _FILE_OFFSET_BITS=64)
