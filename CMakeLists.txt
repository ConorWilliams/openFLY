cmake_minimum_required(VERSION 3.19)

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

# Add the cmake folder so custom
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# This variable is set by project() in CMake 3.21+
if(NOT PROJECT_IS_TOP_LEVEL)
    string(COMPARE EQUAL "${CMAKE_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}" PROJECT_IS_TOP_LEVEL)
endif()

# Avoid warning about DOWNLOAD_EXTRACT_TIMESTAMP in CMake 3.24:
if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.24.0")
    cmake_policy(SET CMP0135 NEW)
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

target_compile_features(libFLY PUBLIC cxx_std_20)

# Find depenedancies
include(ExternalProject)

# System
find_package(OpenMP REQUIRED)

find_package(KIM-API 2.3.0 CONFIG PATHS ${CMAKE_SOURCE_DIR}/build/kim-api/kim-api/install/share/cmake/kim-api)

# ${KIM_API_}
message(STATUS "Using KIM-API @\"${kim-api_BINARY_DIR}\"")

if(NOT KIM-API_FOUND)
    message(FATAL_ERROR "The KIM-API could not be found, install a local copy with:\n\"./scripts/kim-api/install.sh\"")
endif()

# Lib Backward is easiest as a submodule/subdirectory.
add_subdirectory(${CMAKE_SOURCE_DIR}/external/backward-cpp)

# All managed by vcpkg
find_package(fmt CONFIG REQUIRED)
find_package(Eigen3 CONFIG REQUIRED)
find_package(span-lite CONFIG REQUIRED)
find_package(xxHash CONFIG REQUIRED)
find_package(cereal CONFIG REQUIRED)

# Manually install lib Nauty. This is a bit of a pain, but it is the only way to get it to work.

# Ancient voodoo that is autoconf
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

# Link libraries
target_link_libraries(libFLY PRIVATE
    Nauty
    xxHash::xxhash
    KIM-API::kim-api
)

target_link_libraries(libFLY PUBLIC
    Eigen3::Eigen
    fmt::fmt
    nonstd::span-lite
    OpenMP::OpenMP_CXX
    cereal::cereal
)

# Declare executable(s)
# add_executable(skmcFLY apps/skmc.cpp ${BACKWARD_ENABLE})
# target_compile_features(skmcFLY PRIVATE cxx_std_17)
# target_link_libraries(skmcFLY PRIVATE libFLY)
# add_backward(skmcFLY)
add_executable(hydrogenFLY apps/hydrogen.cpp ${BACKWARD_ENABLE})
target_compile_features(hydrogenFLY PRIVATE cxx_std_20)
target_link_libraries(hydrogenFLY PRIVATE libFLY)
add_backward(hydrogenFLY)

# Dev mode
if(FLY_DEVELOPER_MODE)
    find_program(CCACHE "ccache")

    if(CCACHE)
        message(STATUS "Using ccache: ${CCACHE}")
        set(CMAKE_C_COMPILER_LAUNCHER ${CCACHE})
        set(CMAKE_CXX_COMPILER_LAUNCHER ${CCACHE})
    else()
        message(STATUS "Could not find ccache")
    endif()

    include(CTest)
    add_subdirectory(tests)

    option(FLY_SKIP_EXAMPLES "Disable examples" OFF)

    if(NOT FLY_SKIP_EXAMPLES)
        add_subdirectory(examples)
    endif()
endif()

option(FLY_DOCS "Build documentation using doxygen and sphinx " OFF)

if(FLY_DOCS)
    add_subdirectory(docs)
endif()
