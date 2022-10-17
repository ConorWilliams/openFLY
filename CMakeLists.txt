cmake_minimum_required(VERSION 3.16...3.22)

# In-source guard
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
    message(
        FATAL_ERROR
        "In-source builds are not supported. "
        "Please read the BUILDING document before trying to build this project. "
        "You may need to delete 'CMakeCache.txt' and 'CMakeFiles/' first."
    )
endif()

project(
    openFLY
    DESCRIPTION "The open off-lattice kniteic Monte Carlo simulation framework"
    HOMEPAGE_URL "https://github.com/ConorWilliams/openFLY"
    LANGUAGES CXX C
)

# This variable is set by project() in CMake 3.21+
if(NOT PROJECT_IS_TOP_LEVEL)
    string(COMPARE EQUAL "${CMAKE_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}" PROJECT_IS_TOP_LEVEL)
endif()

# Developer mode enables targets and code paths in the CMake scripts that are only relevant for the
# developer(s) of openFLY Targets necessary to build the project must be provided unconditionally,
# so consumers can trivially build and package the project
option(FLY_DEVELOPER_MODE "Enable developer mode" PROJECT_IS_TOP_LEVEL)

include(cmake/GetVersion.cmake)

message(STATUS "openFLY version: ${FLY_VERSION}")

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# Get the sources
file(GLOB_RECURSE FLY_SOURCES_C CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/src/*.c")
file(GLOB_RECURSE FLY_SOURCES_CPP CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/src/*.cpp")

# Declare the lib, Manually include GSD here as its cmake is not for consumption
add_library(libFLY OBJECT ${FLY_SOURCES_C} ${FLY_SOURCES_CPP} "${PROJECT_SOURCE_DIR}/external/gsd/gsd/gsd.c")

target_include_directories(libFLY PUBLIC "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>")

target_include_directories(libFLY PRIVATE "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>")
target_include_directories(libFLY PRIVATE "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/external>")

target_compile_features(libFLY PUBLIC cxx_std_17)

# Link libraries
include(ExternalProject)

# System
find_package(OpenMP REQUIRED)

# All managed by vcpkg
find_package(fmt CONFIG REQUIRED)
find_package(Eigen3 CONFIG REQUIRED)
find_package(span-lite CONFIG REQUIRED)
find_package(xxHash CONFIG REQUIRED)

# Ancient voodoo
ExternalProject_Add(
    nauty
    URL "https://pallini.di.uniroma1.it/nauty27r4.tar.gz"
    URL_HASH "SHA256=b810c85a6fe299f3b4c9f24aaf929cac7f9546c2f35c20e1dd0adbc7408848a6"
    PREFIX ${CMAKE_BINARY_DIR}/nauty
    CONFIGURE_COMMAND ${CMAKE_BINARY_DIR}/nauty/src/nauty/configure --enable-tls
    BUILD_COMMAND make -j4
    TEST_COMMAND make checks
    INSTALL_COMMAND ""
    BUILD_IN_SOURCE 1
    LOG_CONFIGURE ON
    LOG_DOWNLOAD ON
    LOG_INSTALL ON
    LOG_BUILD ON
    LOG_TEST ON
    LOG_OUTPUT_ON_FAILURE ON
)

add_library(Nauty STATIC IMPORTED ${CMAKE_BINARY_DIR}/nauty/src/nauty/nauty.a)

set_target_properties(Nauty PROPERTIES
    IMPORTED_LOCATION ${CMAKE_BINARY_DIR}/nauty/src/nauty/nauty.a
)

# Key headers are: "nauty/nauty.h", "nauty/nautinv.h", "nauty/naututil.h", "nauty/nausparse.h"
target_include_directories(Nauty INTERFACE
    ${CMAKE_BINARY_DIR}/nauty/src
)

add_dependencies(Nauty nauty)

target_link_libraries(libFLY PRIVATE Nauty xxHash::xxhash)

target_link_libraries(libFLY PUBLIC Eigen3::Eigen fmt::fmt nonstd::span-lite OpenMP::OpenMP_CXX)

# Declare executable
add_executable(openFLY apps/main.cpp)
add_executable(openfly::exe ALIAS openFLY)
set_target_properties(openFLY PROPERTIES OUTPUT_NAME openfly)
target_compile_features(openFLY PRIVATE cxx_std_17)
target_link_libraries(openFLY PRIVATE libFLY)
target_link_libraries(openFLY PRIVATE fmt::fmt)

add_executable(geoFLY apps/geo.cpp)
target_compile_features(geoFLY PRIVATE cxx_std_17)
target_link_libraries(geoFLY PRIVATE libFLY Nauty)

# Dev mode
if(FLY_DEVELOPER_MODE)
    include(CTest)
    add_subdirectory(tests)
    add_subdirectory(examples)
endif()

option(FLY_DOCS " Build documentation using ... " OFF)

if(FLY_DOCS)
    add_subdirectory(docs)
endif()