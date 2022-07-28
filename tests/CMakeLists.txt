find_package(Catch2 REQUIRED)

# Note: CMake is better at globbing now...
file(GLOB_RECURSE FLY_TEST_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

add_executable(testFLY ${FLY_TEST_SOURCES})

target_compile_options(testFLY PRIVATE
  $<$<CXX_COMPILER_ID:MSVC>:/W4>
  $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall -Wextra -Wpedantic>
)

target_compile_features(testFLY PRIVATE cxx_std_17)

target_link_libraries(testFLY PRIVATE libFLY::libFLY Catch2::Catch2WithMain)

add_test(NAME testFLY COMMAND testFLY) 

