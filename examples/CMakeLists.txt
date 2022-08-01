file(GLOB_RECURSE FLY_EXAMPLE_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

add_executable(examplesFLY ${FLY_EXAMPLE_SOURCES})

target_compile_features(examplesFLY PRIVATE cxx_std_17)

target_link_libraries(examplesFLY PRIVATE libFLY)
