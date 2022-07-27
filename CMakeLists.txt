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
option(FLY_DOC "Generate the doc target." ${FLY_MASTER_PROJECT})
option(FLY_INSTALL "Generate the install target." ${FLY_MASTER_PROJECT})
option(FLY_APPS "Generate the apps target." ${FLY_MASTER_PROJECT})
option(FLY_EXAMPLE "Generate the example target." ${FLY_MASTER_PROJECT})
option(FLY_TEST "Generate the test target." ${FLY_MASTER_PROJECT})
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

# For ${CMAKE_INSTALL_INCLUDEDIR}
include(GNUInstallDirs)

# Note: CMake is better at globbing now...
file(GLOB_RECURSE FLY_HEADERS CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp")
file(GLOB_RECURSE FLY_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")

# Define the openFLY library, libFLY, its includes and the needed defines.
add_library(fly ${FLY_SOURCES} ${FLY_HEADERS})
add_library(fly::fly ALIAS fly)

if (FLY_PEDANTIC)
  target_compile_options(fly PRIVATE
    $<$<CXX_COMPILER_ID:MSVC>:/W4>
    $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall -Wextra -Wpedantic>
  )
endif ()

target_include_directories(fly PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# Debug library postfix.
set(FLY_DEBUG_POSTFIX d)

set_target_properties(fly PROPERTIES
  VERSION ${FLY_VERSION} 
  SOVERSION ${CPACK_PACKAGE_VERSION_MAJOR}
  PUBLIC_HEADER "${FLY_HEADERS}"
  DEBUG_POSTFIX "${FLY_DEBUG_POSTFIX}"
)

# Set FLY_LIB_NAME for pkg-config fly.pc. W
set(FLY_LIB_NAME fly)

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(FLY_LIB_NAME ${FLY_LIB_NAME}${FLY_DEBUG_POSTFIX})
endif ()




# # Install targets.
# if (FLY_INSTALL)
#   include(CMakePackageConfigHelpers)
#   set_verbose(FLY_CMAKE_DIR ${CMAKE_INSTALL_LIBDIR}/cmake/fly CACHE STRING
#               "Installation directory for cmake files, a relative path that "
#               "will be joined with ${CMAKE_INSTALL_PREFIX} or an absolute "
#               "path.")
#   set(version_config ${PROJECT_BINARY_DIR}/fly-config-version.cmake)
#   set(project_config ${PROJECT_BINARY_DIR}/fly-config.cmake)
#   set(pkgconfig ${PROJECT_BINARY_DIR}/fly.pc)
#   set(targets_export_name fly-targets)

#   set_verbose(FLY_LIB_DIR ${CMAKE_INSTALL_LIBDIR} CACHE STRING
#               "Installation directory for libraries, a relative path that "
#               "will be joined to ${CMAKE_INSTALL_PREFIX} or an absolute path.")

#   set_verbose(FLY_PKGCONFIG_DIR ${CMAKE_INSTALL_LIBDIR}/pkgconfig CACHE PATH
#               "Installation directory for pkgconfig (.pc) files, a relative "
#               "path that will be joined with ${CMAKE_INSTALL_PREFIX} or an "
#               "absolute path.")

#   # Generate the version, config and target files into the build directory.
#   write_basic_package_version_file(
#     ${version_config}
#     VERSION ${FLY_VERSION}
#     COMPATIBILITY AnyNewerVersion)

#   join_paths(libdir_for_pc_file "\${exec_prefix}" "${FLY_LIB_DIR}")
#   join_paths(includedir_for_pc_file "\${prefix}" "${FLY_INC_DIR}")

#   configure_file(
#     "${PROJECT_SOURCE_DIR}/support/cmake/fly.pc.in"
#     "${pkgconfig}"
#     @ONLY)
#   configure_package_config_file(
#     ${PROJECT_SOURCE_DIR}/support/cmake/fly-config.cmake.in
#     ${project_config}
#     INSTALL_DESTINATION ${FLY_CMAKE_DIR})

#   set(INSTALL_TARGETS fly fly-header-only)

#   # Install the library and headers.
#   install(TARGETS ${INSTALL_TARGETS} EXPORT ${targets_export_name}
#           LIBRARY DESTINATION ${FLY_LIB_DIR}
#           ARCHIVE DESTINATION ${FLY_LIB_DIR}
#           PUBLIC_HEADER DESTINATION "${FLY_INC_DIR}/fly"
#           FRAMEWORK DESTINATION "."
#           RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

#   # Use a namespace because CMake provides better diagnostics for namespaced
#   # imported targets.
#   export(TARGETS ${INSTALL_TARGETS} NAMESPACE fly::
#          FILE ${PROJECT_BINARY_DIR}/${targets_export_name}.cmake)

#   # Install version, config and target files.
#   install(
#     FILES ${project_config} ${version_config}
#     DESTINATION ${FLY_CMAKE_DIR})
#   install(EXPORT ${targets_export_name} DESTINATION ${FLY_CMAKE_DIR}
#           NAMESPACE fly::)

#   install(FILES $<TARGET_PDB_FILE:${INSTALL_TARGETS}>
#           DESTINATION ${FLY_LIB_DIR} OPTIONAL)
#   install(FILES "${pkgconfig}" DESTINATION "${FLY_PKGCONFIG_DIR}")
# endif ()

# if (FLY_DOC)
#   add_subdirectory(doc)
# endif ()

# if (FLY_TEST)
#   enable_testing()
#   add_subdirectory(test)
# endif ()

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