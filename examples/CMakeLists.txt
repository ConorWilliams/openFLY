# Note: CMake is better at globbing now...
file(GLOB_RECURSE FLY_TEST_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

add_executable(examplesFLY ${FLY_TEST_SOURCES})

target_compile_options(examplesFLY PRIVATE
  $<$<CXX_COMPILER_ID:MSVC>:/W4>
  $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall -Wextra -Wpedantic>
)

target_compile_features(examplesFLY PRIVATE cxx_std_17)

target_link_libraries(examplesFLY PRIVATE fly::fly)


