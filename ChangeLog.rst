Changelog
============================


Unreleased
-------------------------------

Added
~~~~~~~~~

- ``AdjacentCells`` class
- ``Vector`` class
- ``Xoshiro`` class as the PRNG
- ``Box`` has a new ``get()`` method to fetch the underlying ``std::variant`` 
- ``visit`` utility function.

Changed
~~~~~~~~~~

- ``Property``'s ``array_ref_t`` is now an ``Eigen::Map`` to disallow resizing of ``SoA``'s individual arrays.
- ``SoA``'s ``destructive_resize`` now returns a boolean.
- Allow zero length ``SoA``s
- Crystal systems ``gen_image`` marked ``const``.
- ``BinaryFile`` internals reworked to support clang.

Fixes
~~~~~~~~~~~~~~
- ``SoA``'s ``operator()`` was broken.
- ``SoA``'s ``resize`` was broken.

Removed
~~~~~~~~~

- ``Box`` default constructor

Meta 
~~~~~~~~~~~~~~~~~~~~~~~~~

- vcpkg + gsd are now submodules so all dependencies are tracked by Dependabot!


Version 0.2.0
--------------------------------

The first released alpha version of openFLY! This is a minimal feature-set release with just the base classes that underpin openFLY as well as binary IO.

Added
~~~~~~~~~

- Binary IO using the GSD format through the ``BinaryFile`` class.

- ``Property`` base class template.
- ``TypeMap`` class template.
- ``Supercell`` class template.
- ``SoA`` class template.
- ``VoS`` class template.
- ``Atom`` class template.
- ``Box`` class and specialised crystal systems that it is built on.

- The utility.hpp file containing many utilities.

- New CI workflow now includes C++20 and Intel compilers, MSVC removed due to compiler bug.

Changed
~~~~~~~~~~

- GPL-2.0 -> GPL-3.0-or-later

Removed
~~~~~~~~~

Meta 
~~~~~~~~~~~~~~~~~~~~~~~~~

- Hdoc is no longer used to build the documentation.

Version 0.1.0-alpha
---------------------------

Continuous pre-release, not currently in a usable state.