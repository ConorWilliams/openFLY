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

target_include_directories(fly PUBLIC ${PROJECT_SOURCE_DIR}/include)

