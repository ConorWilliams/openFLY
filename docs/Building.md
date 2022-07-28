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

To run the tests with ctest:
```bash
ctest --test-dir build  
```
or to run them directly:
```bash
./build/bin/testFLY
```
