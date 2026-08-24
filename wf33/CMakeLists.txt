find_package(Threads REQUIRED)
find_package(OpenGL REQUIRED)
find_package(glfw3 CONFIG REQUIRED)
find_package(imgui CONFIG REQUIRED)

if(MSVC)
  # CMake's Windows defaults contain /W3. Remove it before applying /W4 to
  # avoid command-line warning D9025 about one warning level overriding another.
  string(REGEX REPLACE "([/-])W[0-4]" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
endif()

add_library(siproswf_core STATIC
  src/workflow.cpp
  src/process.cpp
)
target_include_directories(siproswf_core PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/include
)
target_link_libraries(siproswf_core PUBLIC Threads::Threads)
target_compile_features(siproswf_core PUBLIC cxx_std_17)
target_compile_definitions(siproswf_core PUBLIC
  SIPROS_VERSION="${SIPROS_VERSION}"
)

if(MSVC)
  target_compile_options(siproswf_core PRIVATE /utf-8 /W4)
  target_compile_definitions(siproswf_core PRIVATE _CRT_SECURE_NO_WARNINGS)
else()
  target_compile_options(siproswf_core PRIVATE -Wall -Wextra -Wpedantic)
endif()

set(SIPROSWF_SOURCES
  src/main.cpp
  src/gui.cpp
  src/native_dialog.cpp
)
if(WIN32)
  list(APPEND SIPROSWF_SOURCES siproswf.rc siproswf.ico siproswf.manifest)
  add_executable(siproswf WIN32 ${SIPROSWF_SOURCES})
else()
  add_executable(siproswf ${SIPROSWF_SOURCES})
endif()
target_link_libraries(siproswf PRIVATE
  siproswf_core
  imgui::imgui
  glfw
  OpenGL::GL
)
target_compile_features(siproswf PRIVATE cxx_std_17)

if(MSVC)
  target_compile_options(siproswf PRIVATE /utf-8 /W4)
  target_compile_definitions(siproswf PRIVATE _CRT_SECURE_NO_WARNINGS)
elseif(NOT WIN32)
  target_compile_options(siproswf PRIVATE -Wall -Wextra -Wpedantic)
endif()

if(WIN32)
  target_link_libraries(siproswf PRIVATE ole32 shell32 uuid delayimp)
  target_link_options(siproswf PRIVATE
    /DELAYLOAD:imgui.dll
    /DELAYLOAD:glfw3.dll)
endif()

if(BUILD_TESTING)
  add_executable(siproswf_fake_tool tests/fake_tool.cpp)
  if(MSVC)
    target_compile_definitions(siproswf_fake_tool PRIVATE _CRT_SECURE_NO_WARNINGS)
  endif()
  add_executable(siproswf_test tests/test_workflow.cpp)
  add_dependencies(siproswf_test siproswf_fake_tool)
  target_link_libraries(siproswf_test PRIVATE siproswf_core)
  add_test(NAME siproswf_core COMMAND siproswf_test $<TARGET_FILE:siproswf_fake_tool>)
endif()
