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
    LANGUAGES CXX
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

# Declare library
file(GLOB_RECURSE FLY_SOURCES CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/src/*.cpp")

add_library(libFLY OBJECT ${FLY_SOURCES})

target_include_directories(libFLY PUBLIC "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>")

target_compile_features(libFLY PUBLIC cxx_std_17)

find_package(fmt REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(span-lite REQUIRED)

# target_link_libraries(libFLY PRIVATE )
target_link_libraries(libFLY PUBLIC Eigen3::Eigen fmt::fmt nonstd::span-lite)

# Declare executable
add_executable(openFLY apps/main.cpp)

add_executable(openfly::exe ALIAS openFLY)

set_target_properties(openFLY PROPERTIES OUTPUT_NAME openfly)

target_compile_features(openFLY PRIVATE cxx_std_17)

target_link_libraries(openFLY PRIVATE libFLY)

target_link_libraries(openFLY PRIVATE fmt::fmt)

# Dev mode
if(FLY_DEVELOPER_MODE)
    include(CTest)

    add_subdirectory(tests)

    add_subdirectory(examples)
endif()

option(FLY_DOCS "Build documentation using ..." OFF)

if(FLY_DOCS)
    add_subdirectory(docs)
endif()