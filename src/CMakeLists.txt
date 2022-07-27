
# For ${CMAKE_INSTALL_INCLUDEDIR}
include(GNUInstallDirs)

# Note: CMake is better at globbing now...
file(GLOB_RECURSE FLY_HEADERS CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/../include/*.hpp")
file(GLOB_RECURSE FLY_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

# Define the openFLY library, libFLY, its includes and the needed defines.
add_library(fly ${FLY_SOURCES} ${FLY_HEADERS})
add_library(fly::fly ALIAS fly)

if (FLY_PEDANTIC)
  target_compile_options(fly PRIVATE
    $<$<CXX_COMPILER_ID:MSVC>:/W4>
    $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall -Wextra -Wpedantic>
  )
endif ()

target_compile_features(fly PUBLIC cxx_std_17)

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