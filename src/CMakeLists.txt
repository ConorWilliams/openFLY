# Note: CMake is better at globbing now...
file(GLOB_RECURSE FLY_HEADERS CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/../include/*.hpp")
file(GLOB_RECURSE FLY_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

# Define the openFLY library, libFLY, its includes and the needed defines.
add_library(libFLY ${FLY_SOURCES} ${FLY_HEADERS})
add_library(libFLY::libFLY ALIAS libFLY)

if(FLY_PEDANTIC)
  target_compile_options(
    libFLY PRIVATE $<$<CXX_COMPILER_ID:MSVC>:/W4> $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall -Wextra
                   -Wpedantic>
  )
endif()

target_compile_features(libFLY PUBLIC cxx_std_17)

target_include_directories(libFLY PUBLIC ${PROJECT_SOURCE_DIR}/include)

export(
  TARGETS libFLY
  NAMESPACE libFLY::
  FILE LibFlyTargets.cmake
)
