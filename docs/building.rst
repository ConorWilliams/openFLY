Compiling 
=========

Building libFLY
---------------

LibFLY is the library component of openFLY and designed to be reusable and composable.

Build requirements:

* cmake
* A C++17 compatible compiler that supports openMP
* git

First obtain a copy of the openFLY source code:

.. code:: console

    git clone git@github.com:ConorWilliams/openFLY.git --recursive

and enter the directory:

.. code:: console

    cd openFLY 

OpenFLY uses `vcpkg <https://github.com/microsoft/vcpkg>`_ to manage its dependencies; to bootstrap vcpkg:

.. code:: console

    ./vcpkg/bootstrap-vcpkg.sh   

.. note::
    If you already use vcpkg you can omit the ``--recursive`` from the ``git clone ...`` call and skip the bootstrapping step, in the configuration step you will need to modify ``-DCMAKE_TOOLCHAIN_FILE=...`` to point to your installation's ``vcpkg.cmake``.

Now we can configure our build, the simplest option is to run:

.. code:: console

    cmake -S . -B build -DCMAKE_TOOLCHAIN_FILE=vcpkg/scripts/buildsystems/vcpkg.cmake 

You can pass in many optional cmake parameters at this point to customize the compiler, linker, compiler-options, etc. If you have ninja installed try appending ``-GNinja``. 

Finally, to build libFLY:

.. code:: console

    cmake --build build

Building the openFLY application(s)
-----------------------------------

Build requirements:

* LibFLY (see above)
* MPI
* ...


Examples and tests
------------------

Documentation
------------------

Requirements
~~~~~~~~~~~~~~~~~~

The requirements are:

- CMake 
- 
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
