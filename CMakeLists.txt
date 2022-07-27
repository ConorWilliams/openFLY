cmake_minimum_required(VERSION 3.16...3.22)

# Determine if openFLY is built as a subproject (using add_subdirectory) or if it is the master project.
if (NOT DEFINED FLY_MASTER_PROJECT)
  set(FLY_MASTER_PROJECT OFF)
  if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    set(FLY_MASTER_PROJECT ON)
    message(STATUS "CMake version: ${CMAKE_VERSION}")
  endif ()
endif ()

if (FLY_MASTER_PROJECT)
  # Ensure -std=c++xx instead of -std=g++xx
  set(CMAKE_CXX_EXTENSIONS OFF)
endif()

project(openFLY CXX)

# Options that control generation of various targets.
option(FLY_APPS "Build the openFLY application." ${FLY_MASTER_PROJECT})
option(FLY_EXAMPLES "Build the examples." ${FLY_MASTER_PROJECT})
option(FLY_TESTS "Build the tests." ${FLY_MASTER_PROJECT})
option(FMT_PEDANTIC "Enable extra warnings." ${FLY_MASTER_PROJECT})


# Get version from version.h
file(READ "include/libfly/version.hpp" version)

if (NOT version MATCHES "VERSION_MAJOR ([0-9]+)")
  message(FATAL_ERROR "Cannot get FLY_VERSION from version.hpp.")
else()
  math(EXPR CPACK_PACKAGE_VERSION_MAJOR ${CMAKE_MATCH_1})
endif ()

if (NOT version MATCHES "VERSION_MINOR ([0-9]+)")
  message(FATAL_ERROR "Cannot get VERSION_MINOR from version.hpp.")
else ()
  math(EXPR CPACK_PACKAGE_VERSION_MINOR ${CMAKE_MATCH_1})
endif ()

if (NOT version MATCHES "VERSION_PATCH ([0-9]+)")
  message(FATAL_ERROR "Cannot get VERSION_PATCH from version.hpp.")
else()
  math(EXPR CPACK_PACKAGE_VERSION_PATCH ${CMAKE_MATCH_1})
endif ()

set(FLY_VERSION "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")

message(STATUS "openFLY version: ${FLY_VERSION}")

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

if (NOT CMAKE_RUNTIME_OUTPUT_DIRECTORY)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
endif ()

# libFly is here
add_subdirectory(src)

if (FLY_EXAMPLES)
  add_subdirectory(examples)
endif ()

if (FLY_TESTS)
  include(CTest)
  enable_testing()
  add_subdirectory(tests)
endif ()

# # Control fuzzing independent of the unit tests.
# if (FLY_FUZZ)
#   add_subdirectory(test/fuzzing)

#   # The FLY_FUZZ macro is used to prevent resource exhaustion in fuzzing
#   # mode and make fuzzing practically possible. It is similar to
#   # FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION but uses a different name to
#   # avoid interfering with fuzzing of projects that use {fly}.
#   # See also https://llvm.org/docs/LibFuzzer.html#fuzzer-friendly-build-mode.
#   target_compile_definitions(fly PUBLIC FLY_FUZZ)
# endif ()

# set(gitignore ${PROJECT_SOURCE_DIR}/.gitignore)
# if (FLY_MASTER_PROJECT AND EXISTS ${gitignore})
#   # Get the list of ignored files from .gitignore.
#   file (STRINGS ${gitignore} lines)
#   list(REMOVE_ITEM lines /doc/html)
#   foreach (line ${lines})
#     string(REPLACE "." "[.]" line "${line}")
#     string(REPLACE "*" ".*" line "${line}")
#     set(ignored_files ${ignored_files} "${line}$" "${line}/")
#   endforeach ()
#   set(ignored_files ${ignored_files}
#     /.git /breathe /format-benchmark sphinx/ .buildinfo .doctrees)

#   set(CPACK_SOURCE_GENERATOR ZIP)
#   set(CPACK_SOURCE_IGNORE_FILES ${ignored_files})
#   set(CPACK_SOURCE_PACKAGE_FILE_NAME fly-${FLY_VERSION})
#   set(CPACK_PACKAGE_NAME fly)
#   set(CPACK_RESOURCE_FILE_README ${PROJECT_SOURCE_DIR}/README.rst)
#   include(CPack)
# endif ()