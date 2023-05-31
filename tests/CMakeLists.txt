# ---- Dependencies ----

find_package(Catch2 REQUIRED)

include(Catch)

# ---- Tests ----
file(GLOB_RECURSE FLY_TEST_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

add_executable(testFLY ${FLY_TEST_SOURCES})

target_link_libraries(testFLY PRIVATE libFLY Catch2::Catch2WithMain Nauty)

target_compile_features(testFLY PRIVATE cxx_std_20)

catch_discover_tests(testFLY)
