Hacking
===============

Here is some wisdom to help you build and test this project as a developer and/or potential contributor. 

The build system targets (examples and test) that are only useful for developers of this project are hidden if the ``FLY_DEVELOPER_MODE`` option is disabled. Enabling this option makes tests and other developer targets/options available. Not enabling this option means that you are a consumer of this project and thus you have no need for these targets/options.

Presets
--------

This project makes use of `presets <https://cmake.org/cmake/help/latest/manual/cmake-presets.7.html>`_ to simplify the process of configuring the project. As a developer, you are recommended to always have the `latest CMake version <https://cmake.org/download/>`_ installed to make use of the latest Quality-of-Life additions like presets.

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

You should replace ``<os>`` in your newly created ``CMakeUserPresets.json`` file with the name of the operating system you have, which may be ``win64`` or ``unix``. You can see what these correspond to in the `CMakePresets.json <https://github.com/ConorWilliams/openFLY/blob/master/CMakePresets.json>`_ file.

``CMakeUserPresets.json`` is also the perfect place in which you can put all sorts of things that you would otherwise want to pass to the configure command in the terminal.

The above preset will make use of vcpkg. Make sure the ``VCPKG_ROOT`` environment variable is pointing at the directory where the vcpkg executable is. On Windows, you might also want to inherit from the ``vcpkg-win64-static`` preset, which will make vcpkg install the dependencies as static libraries. This is only necessary if you don't want to setup ``PATH`` to run tests.

Configure, build and test
--------------------------

In developer mode the test and examples will be built automatically. If you followed the above instructions, then you can configure, build and test the project respectively with the following commands from the project root on any operating system with any build system:

.. code:: console

    cmake --preset=dev
    cmake --build --preset=dev
    ctest --preset=dev

.. note::

    Both the build and test commands accept a ``-j`` flag to specify the number of jobs to use, which should ideally be specified to the number of threads your CPU has. You may also want to add that to your preset using the ``jobs`` property, see the `presets documentation <https://cmake.org/cmake/help/latest/manual/cmake-presets.7.html>`_ for more details.


