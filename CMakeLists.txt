cmake_minimum_required(VERSION 3.15)
project(${SKBUILD_PROJECT_NAME} LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find Python
find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)

# Include directories
set(CPYPP_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/deps/cpypp/include")
set(FBITSET_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/deps/fbitset/include")
set(LIBPARENTH_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/deps/libparenth/include")

# Platform-specific compiler flags
if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
  message(STATUS "Configuring for macOS (Darwin)")
  set(CMAKE_OSX_DEPLOYMENT_TARGET "10.9")
  add_compile_options(-stdlib=libc++)
# elseif(MSVC)
#   # MSVC compiler flags
#   add_compile_options(/EHsc /bigobj /wd4996 /wd4267 /Zc:twoPhase-)
endif()

# parenth extension
Python_add_library(_parenth MODULE
  gristmill/_parenth.cpp
  WITH_SOABI
)
target_include_directories(_parenth PRIVATE
  ${CPYPP_INCLUDE_DIR}
  ${FBITSET_INCLUDE_DIR}
  ${LIBPARENTH_INCLUDE_DIR}
)
install(TARGETS _parenth LIBRARY DESTINATION gristmill)
