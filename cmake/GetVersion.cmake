# Get version from version.h
file(READ "include/libfly/utility/version.hpp" version)

if(NOT version MATCHES "VERSION_MAJOR ([0-9]+)")
    message(FATAL_ERROR "Cannot get FLY_VERSION from version.hpp.")
else()
    math(EXPR PROJECT_VERSION_MAJOR ${CMAKE_MATCH_1})
endif()

if(NOT version MATCHES "VERSION_MINOR ([0-9]+)")
    message(FATAL_ERROR "Cannot get VERSION_MINOR from version.hpp.")
else()
    math(EXPR PROJECT_VERSION_MINOR ${CMAKE_MATCH_1})
endif()

if(NOT version MATCHES "VERSION_PATCH ([0-9]+)")
    message(FATAL_ERROR "Cannot get VERSION_PATCH from version.hpp.")
else()
    math(EXPR PROJECT_VERSION_PATCH ${CMAKE_MATCH_1})
endif()

set(FLY_VERSION "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}")
