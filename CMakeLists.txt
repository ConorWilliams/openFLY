cmake_minimum_required(VERSION 3.16...3.22)

# Determine if openFLY is built as a subproject (using add_subdirectory) or if it is the master
# project.
if(NOT DEFINED FLY_MASTER_PROJECT)
  set(FLY_MASTER_PROJECT OFF)
  if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    set(FLY_MASTER_PROJECT ON)
    message(STATUS "CMake version: ${CMAKE_VERSION}")
  endif()
endif()

if(FLY_MASTER_PROJECT AND NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
  message(
    STATUS "Choose the type of build, options are: Debug Release RelWithDebInfo or MinSizeRel"
  )
endif()

project(openFLY CXX)

option(FLY_ALL "Build everything all at once" Off)

# Options that control generation of various targets.
option(FLY_APPS "Build the openFLY application." ${FLY_ALL})
option(FLY_EXAMPLES "Build the examples." ${FLY_ALL})
option(FLY_TESTS "Build the tests." ${FLY_ALL})
option(FLY_DOCS "Build the documentation." ${FLY_ALL})
option(FMT_PEDANTIC "Enable extra warnings." ${FLY_ALL})

# Get version from version.h
file(READ "include/libfly/version.hpp" version)

if(NOT version MATCHES "VERSION_MAJOR ([0-9]+)")
  message(FATAL_ERROR "Cannot get FLY_VERSION from version.hpp.")
else()
  math(EXPR CPACK_PACKAGE_VERSION_MAJOR ${CMAKE_MATCH_1})
endif()

if(NOT version MATCHES "VERSION_MINOR ([0-9]+)")
  message(FATAL_ERROR "Cannot get VERSION_MINOR from version.hpp.")
else()
  math(EXPR CPACK_PACKAGE_VERSION_MINOR ${CMAKE_MATCH_1})
endif()

if(NOT version MATCHES "VERSION_PATCH ([0-9]+)")
  message(FATAL_ERROR "Cannot get VERSION_PATCH from version.hpp.")
else()
  math(EXPR CPACK_PACKAGE_VERSION_PATCH ${CMAKE_MATCH_1})
endif()

set(FLY_VERSION
    "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}"
)

message(STATUS "openFLY version: ${FLY_VERSION}")

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

if(NOT CMAKE_RUNTIME_OUTPUT_DIRECTORY)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
  message(STATUS "Binaries will be written to: ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
endif()

# libFly is here
add_subdirectory(src)

if(FLY_APPS)
  add_subdirectory(apps)
endif()

if(FLY_EXAMPLES)
  add_subdirectory(examples)
endif()

if(FLY_TESTS)
  include(CTest)
  enable_testing()
  add_subdirectory(tests)
endif()

if(FLY_DOCS)
  add_subdirectory(docs)
endif()
