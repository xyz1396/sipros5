aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/src/expat-2.2.9/src expatSRC)
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/src/MSToolkit MSToolkitSRC)
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/src/mzIMLTools mzIMLToolsSRC)
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/src/mzParser mzParserSRC)
## add sqlite
# aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/src/sqlite-3.32.1 sqliteSRC)
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/src/zlib-1.2.11/src zlibSRC)

## add sqlite
# add_library(MStoolkit ${expatSRC} ${MSToolkitSRC} ${mzIMLToolsSRC} ${mzParserSRC} ${sqliteSRC} ${zlibSRC})
add_library(MStoolkit ${expatSRC} ${MSToolkitSRC} ${mzIMLToolsSRC} ${mzParserSRC} ${zlibSRC})
## for .so lib
# set_target_properties(MStoolkit PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)
## for .a lib
# set_target_properties(MStoolkit PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)
## for Linux + glibc >=2.25 has getrandom in <sys/random.h>
# target_compile_definitions(MStoolkit PUBLIC GCC HAVE_GETRANDOM _FILE_OFFSET_BITS=64 _NO_THERMORAW _NOSQLITE)
## for kernel version < 3.17 without SYS_getrandom
# target_compile_definitions(MStoolkit PUBLIC GCC HAVE_SYSCALL_GETRANDOM _FILE_OFFSET_BITS=64 _NO_THERMORAW _NOSQLITE)
# target_compile_definitions(MStoolkit PUBLIC GCC XML_POOR_ENTROPY _FILE_OFFSET_BITS=64 _NO_THERMORAW _NOSQLITE)
target_compile_definitions(MStoolkit PUBLIC GCC XML_DEV_URANDOM _FILE_OFFSET_BITS=64 _NO_THERMORAW _NOSQLITE)
## add sqlite
# target_compile_definitions(MStoolkit PUBLIC GCC HAVE_GETRANDOM _FILE_OFFSET_BITS=64 _NO_THERMORAW)

