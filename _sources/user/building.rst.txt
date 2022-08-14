Building 
========

Prerequisites
-------------

OpenFLY requires:

* A C++17/20 compatible compiler that supports `openMP 4.5+ <https://www.openmp.org/specifications/>`_
* cmake 3.14+
* git
* MPI installation
* doxygen 1.9+ [optional, for documentation]
* python 3.6+ [optional, for documentation]

Additionally openFLY has a number of dependencies which are managed automatically through  `vcpkg <https://github.com/microsoft/vcpkg>`_ and git submodules.


Getting the source code
----------------------------------

First, download a copy of the source code using git:

.. tab:: https

    .. code:: console

        git clone --recursive https://github.com/ConorWilliams/openFLY.git

.. tab:: ssh

    .. code:: console

        git clone --recursive git@github.com:ConorWilliams/openFLY.git

then `cd` into the newly created folder:

.. code:: console

    cd openFLY

Compiling openFLY
------------------

To build openFLY in release mode with:

.. tab:: Single-configuration generator i.e. make, ninja, etc.

    First configure cmake (this will also fetch and build the vcpkg dependancies):

    .. code:: console

        cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE=./external/vcpkg/scripts/buildsystems/vcpkg.cmake

    Then compile the code:

    .. code:: console

        cmake --build build

.. tab:: Multi-configuration generator, i.e. Visual Studio

    First configure cmake (this will also fetch and build the vcpkg dependancies):
    
    .. code:: console

        cmake -S . -B build -DCMAKE_TOOLCHAIN_FILE=./external/vcpkg/scripts/buildsystems/vcpkg.cmake

    Then compile the code:

    .. code:: console

        cmake --build build --config Release 




.. note::
    MSVC by default is not standards compliant and you need to pass some flags to make it behave properly. See the ``flags-windows`` preset in the `CMakePresets.json <https://github.com/ConorWilliams/openFLY/blob/master/CMakePresets.json>`_  file for the flags and with what variable to provide them to CMake during configuration.


Documentation
--------------------------------

If you want to build a local version of the documentation you will require the optional Prerequisites_

First you must install the required python packages:

.. code:: console

    pip3 install -r docs/requirements.txt

Then the documentation will build automatically when supplying the ``-DFLY_DOCS=ON`` to the configure step or in your ``CMakeUserPresets.json`` file if using :ref:`Presets`


