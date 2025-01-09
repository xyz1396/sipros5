set(CMAKE_C_COMPILER "mpicc")
set(CMAKE_CXX_COMPILER "mpicxx")

add_executable(siprosMPI ${SRC_LIST} ${PROJECT_SOURCE_DIR}/mpi/mpimain.cpp)
target_link_libraries(siprosMPI PUBLIC mstoolkit -fopenmp)
target_compile_options(siprosMPI PUBLIC -ffast-math -march=native -fPIC)
target_compile_definitions(siprosMPI PUBLIC GCC _FILE_OFFSET_BITS=64 _NOSQLITE)
