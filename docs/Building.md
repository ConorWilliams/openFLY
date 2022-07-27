# Building

The requirements are:

- CMake 
- A C++17 compatible compiler
- Git

To configure:

```bash
cmake -S . -B build -DCMAKE_TOOLCHAIN_FILE=[path to vcpkg]/scripts/buildsystems/vcpkg.cmake 
```

Add `-GNinja` if you have Ninja.

To build:

```bash
cmake --build build
```

To run the tests (`--target` can be written as `-t` in CMake 3.15+):

```bash
ctest --test-dir build  
```
or run directly
```bash
./build/bin/testFLY
```

test $math = i$