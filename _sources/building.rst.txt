Building 
========

Prerequisite
------------

OpenFLY requires:

* A C++17 compatible compiler that **supports openMP**
* cmake
* git
* MPI

Additionally openFLY has a number of dependencies (please refer to `vcpkg.json <https://github.com/ConorWilliams/openFLY/blob/master/vcpkg.json>`_ for a complete list). You can choose to download and install these manually but, the recommended method is to use `vcpkg <https://github.com/microsoft/vcpkg>`_ to manage them automatically.

.. tip:: 

   To get started with vcpkg in Unix, clone the repo:

   .. code:: console

      git clone  https://vcpkg.io/en/getting-started.html

   then bootstrap using the supplied script:

   .. code:: console

      ./vcpkg/bootstrap-vcpkg.sh

   finally set the ``VCPKG_ROOT`` environmental variable:

   .. code:: console

       export VCPKG_ROOT=<path to vcpkg>  

   For more detail see the official page: https://github.com/microsoft/vcpkg.


Compiling
---------

To build openFLY in release mode with a single-configuration generator, like the Unix Makefiles one:

.. code:: console

   cmake -S . -B build -DCMAKE_TOOLCHAIN_FILE=$VCPKG_ROOT/scripts/buildsystems/vcpkg.cmake -DCMAKE_BUILD_TYPE=Release
   cmake --build build
    

Similarly with a multi-configuration generator, like the Visual Studio ones:

.. code:: console

   cmake -S . -B build
   cmake --build build --config Release

.. note::
    MSVC by default is not standards compliant and you need to pass some flags to make it behave properly. See the ``flags-windows`` preset in the `CMakePresets.json <https://github.com/ConorWilliams/openFLY/blob/master/CMakePresets.json>`_  file for the flags and with what variable to provide them to CMake during configuration.


Hacking
=======

Here is some wisdom to help you build and test this project as a developer and potential contributor. If you'd like to learn more about openFLY's internals, refer to the :ref:`library api reference`.

Developer mode
--------------

Build system targets that are only useful for developers of this project are hidden if the ``FLY_DEVELOPER_MODE`` option is disabled. Enabling this option makes tests and other developer targets and options available. Not enabling this option means that you are a consumer of this project and thus you have no need for these targets and options.

Developer mode is always set to on in CI workflows.

Presets
~~~~~~~

This project makes use of `presets <https://cmake.org/cmake/help/latest/manual/cmake-presets.7.html>`_ to simplify the process of configuring the project. As a developer, you are recommended to always have the `latest CMake version <https://cmake.org/download/>`_ installed to make use of the latest Quality-of-Life additions.

You have a few options to pass ``DEV_DEVELOPER_MODE`` to the configure command, but this project prefers to use presets.

As a developer, you should create a ``CMakeUserPresets.json`` file at the root of the project:

.. code:: json

    {
    "version": 2,
    "configurePresets": [
        {
        "name": "dev",
        "binaryDir": "${sourceDir}/build/dev",
        "inherits": [
            "dev-mode",
            "vcpkg",
            "ci-<os>"
        ],
        "cacheVariables": {
            "FLY_DOCS": "OFF",
            "CMAKE_BUILD_TYPE": "Debug",
            "CMAKE_EXPORT_COMPILE_COMMANDS": "ON"
        }
        }
    ],
    "buildPresets": [
        {
        "name": "dev",
        "configurePreset": "dev",
        "configuration": "Debug"
        }
    ],
    "testPresets": [
        {
        "name": "dev",
        "configurePreset": "dev",
        "configuration": "Debug",
        "output": {
            "outputOnFailure": true
        }
        }
    ]
    }

You should replace ``<os>`` in your newly created presets file with the name of the operating system you have, which may be ``win64`` or ``unix``. You can see what these correspond to in the `CMakePresets.json <https://github.com/ConorWilliams/openFLY/blob/master/CMakePresets.json>`_ file.

``CMakeUserPresets.json`` is also the perfect place in which you can put all sorts of things that you would otherwise want to pass to the configure command in the terminal.

The above preset will make use of vcpkg. Make sure the ``VCPKG_ROOT`` environment variable is pointing at the directory where the vcpkg executable is. On Windows, you might also want to inherit from the `vcpkg-win64-static` preset, which will make vcpkg install the dependencies as static libraries. This is only necessary if you don't want to setup ``PATH`` to run tests.

Configure, build and test
~~~~~~~~~~~~~~~~~~~~~~~~~

If you followed the above instructions, then you can configure, build and test the project respectively with the following commands from the project root on any operating system with any build system:

.. code:: console

    cmake --preset=dev
    cmake --build --preset=dev
    ctest --preset=dev

.. note::

    Both the build and test commands accept a ``-j`` flag to specify the number of jobs to use, which should ideally be specified to the number of threads your CPU has. You may also want to add that to your preset using the ``jobs`` property, see the `presets documentation <https://cmake.org/cmake/help/latest/manual/cmake-presets.7.html>`_ for more details.

Developer mode targets
~~~~~~~~~~~~~~~~~~~~~~

In developer mode the test and examples will be build automatically. Additionally if you want to build a local version of the documentation you will require:

* doxygen
* python3 and packages:
    - sphinx
    - breathe
    - furo

Then you can build the documentation by supplying the ``-DFLY_DOCS=ON`` alongside the configure command or in your presets file.

