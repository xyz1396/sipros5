
# aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/src/expat-2.2.9/src expatSRC)
# aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/src/MSToolkit MSToolkitSRC)
# aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/src/mzIMLTools mzIMLToolsSRC)
# aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/src/mzParser mzParserSRC)

## add sqlite
# aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/src/sqlite-3.32.1 sqliteSRC)
# aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/src/zlib-1.2.11/src zlibSRC)

## add sqlite
# add_library(MStoolkit ${expatSRC} ${MSToolkitSRC} ${mzIMLToolsSRC} ${mzParserSRC} ${sqliteSRC} ${zlibSRC})
# add_library(MStoolkit ${expatSRC} ${MSToolkitSRC} ${mzIMLToolsSRC} ${mzParserSRC} ${zlibSRC})

## for .so lib
# set_target_properties(MStoolkit PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)

## for .a lib
# set_target_properties(MStoolkit PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)

## for Linux + glibc >=2.25 has getrandom in <sys/random.h>
# target_compile_definitions(MStoolkit PUBLIC GCC HAVE_GETRANDOM _FILE_OFFSET_BITS=64 _NO_THERMORAW _NOSQLITE)

## for kernel version < 3.17 without SYS_getrandom
# target_compile_definitions(MStoolkit PUBLIC GCC HAVE_SYSCALL_GETRANDOM _FILE_OFFSET_BITS=64 _NO_THERMORAW _NOSQLITE)
# target_compile_definitions(MStoolkit PUBLIC GCC XML_POOR_ENTROPY _FILE_OFFSET_BITS=64 _NO_THERMORAW _NOSQLITE)
# target_compile_definitions(MStoolkit PUBLIC GCC XML_DEV_URANDOM _FILE_OFFSET_BITS=64 _NO_THERMORAW _NOSQLITE)

## add sqlite
# target_compile_definitions(MStoolkit PUBLIC GCC HAVE_GETRANDOM _FILE_OFFSET_BITS=64 _NO_THERMORAW)


# Source directories
set(MSTOOLKIT_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/MSToolkit)
set(EXPAT_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/expat-2.2.9/lib)
set(MZID_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/mzIMLTools)
set(MZP_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/mzParser)
set(ZLIB_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/zlib-1.2.11)

# Library output directory
set(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/lib)

# for Windows and linux cross platform

IF (WIN32)
  # set stuff for windows
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /permissive- /GS /GL /W3 /Gy /Zc:wchar_t /Zi /Gm- /sdl /Zc:inline /fp:precise /errorReport:prompt /WX- /Zc:forScope /Gd /Oi /MD /FC /EHsc /nologo /diagnostics:column") 
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /permissive- /GS /GL /W3 /Gy /Zc:wchar_t /Zi /Gm- /sdl /Zc:inline /fp:precise /errorReport:prompt /WX- /Zc:forScope /Gd /Oi /MD /FC /EHsc /nologo /diagnostics:column")
  add_definitions(-DWIN32 -D_CRT_SECURE_NO_WARNINGS -D_CRT_NONSTDC_NO_DEPRECATE -DXML_STATIC -D_MBCS -D_NOSQLITE -D_NO_THERMORAW)
  # make sure to use the static runtime library for release build
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
  set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} /MT")
ELSE()
  # for linux
  add_definitions(-DGCC -DXML_DEV_URANDOM -D_FILE_OFFSET_BITS=64 -D_NO_THERMORAW -D_NOSQLITE)
ENDIF()

# Add library targets
file(GLOB MSTOOLKIT_SOURCES "${MSTOOLKIT_SRC_DIR}/*.cpp")
add_library(mstoolkit STATIC ${MSTOOLKIT_SOURCES})
file(GLOB EXPAT_SOURCES "${EXPAT_SRC_DIR}/*.c")
add_library(expat STATIC ${EXPAT_SOURCES})
file(GLOB MZID_SOURCES "${MZID_SRC_DIR}/*.cpp")
add_library(mzid STATIC ${MZID_SOURCES})
file(GLOB MZP_SOURCES "${MZP_SRC_DIR}/*.cpp")
add_library(mzparser STATIC ${MZP_SOURCES})
file(GLOB ZLIB_SOURCES "${ZLIB_SRC_DIR}/*.c")
add_library(zlib STATIC ${ZLIB_SOURCES})

# Link libraries
target_link_libraries(mzparser PRIVATE expat)
target_link_libraries(mstoolkit PRIVATE expat mzid mzparser zlib)
