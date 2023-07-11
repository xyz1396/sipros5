# this one is important
set(CMAKE_SYSTEM_NAME Linux)
set(CMAKE_PLATFORM Linux)
#this one not so much
set(CMAKE_SYSTEM_VERSION 3)

# specify the cross compiler
# set(CMAKE_C_COMPILER ${CC})

# where is the target environment
# set(CMAKE_FIND_ROOT_PATH ${CONDA_PREFIX} ${CONDA_BUILD_SYSROOT}/usr/)
set(CMAKE_FIND_ROOT_PATH "/home/yixiong/miniconda3/envs/mpi" "/home/yixiong/miniconda3/envs/mpi/x86_64-conda-linux-gnu/sysroot/usr/")
message("${CMAKE_FIND_ROOT_PATH}")
# SET(CMAKE_SYSROOT ${CONDA_BUILD_SYSROOT})
# SET(CMAKE_SYSROOT "/home/yixiong/miniconda3/envs/mpi/x86_64-conda-linux-gnu/sysroot/usr/")

# search for programs in the build host directories
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM ONLY)
# for libraries and headers in the target directories
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)

# god-awful hack because it seems to not run correct tests to determine this:
set(__CHAR_UNSIGNED___EXITCODE 1)