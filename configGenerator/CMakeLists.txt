include_directories(src)
file(GLOB SOURCES_configGenerator "src/*.cpp")
add_executable(configGenerator ${SOURCES_configGenerator})
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
