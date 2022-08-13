Building 
========

Prerequisites
-------------

OpenFLY requires:

* A C++17/20 compatible compiler that supports `openMP 4.5+ <https://www.openmp.org/specifications/>`_
* cmake 3.14+
* git
* MPI installation
* doxygen 1.9+ [optional]
* python 3.6+ [optional]

Additionally openFLY has a number of dependencies (please refer to `vcpkg.json <https://github.com/ConorWilliams/openFLY/blob/master/vcpkg.json>`_ for a complete list). You can choose to download and install these manually but, the recommended method is to use `vcpkg <https://github.com/microsoft/vcpkg>`_ to manage them automatically.

.. tip:: 

   To get started with vcpkg in Unix, clone the repo:

   .. code:: console

      git clone https://github.com/microsoft/vcpkg

   then bootstrap using the supplied script:

   .. code:: console

      ./vcpkg/bootstrap-vcpkg.sh

   finally set the ``VCPKG_ROOT`` environmental variable:

   .. code:: console

       export VCPKG_ROOT=<path to vcpkg>  

   For more detail see the official page: https://github.com/microsoft/vcpkg.


Compiling openFLY
------------------

To build openFLY in release mode with a single-configuration generator, like the Unix Makefiles one:

.. code:: console

    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
    cmake --build build


Similarly with a multi-configuration generator, like the Visual Studio ones:

.. code:: console

   cmake -S . -B build
   cmake --build build --config Release

If you are using vcpkg to manage your dependancies then you will need to append:

.. code:: console

    -DCMAKE_TOOLCHAIN_FILE=$VCPKG_ROOT/scripts/buildsystems/vcpkg.cmake

to the configure step.

.. note::
    MSVC by default is not standards compliant and you need to pass some flags to make it behave properly. See the ``flags-windows`` preset in the `CMakePresets.json <https://github.com/ConorWilliams/openFLY/blob/master/CMakePresets.json>`_  file for the flags and with what variable to provide them to CMake during configuration.


Documentation
--------------------------------

If you want to build a local version of the documentation you will require the optional Prerequisites_

First you must install the required python packages:

.. code:: console

    pip3 install -r docs/requirements.txt

Then the documentation will build automatically when supplying the ``-DFLY_DOCS=ON`` to the configure step or in your ``CMakeUserPresets.json`` file if using :ref:`Presets`


