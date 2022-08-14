Changelog
============================


Unreleased
-------------------------------

Added
~~~~~~~~~

- Adjacent cells class

Changed
~~~~~~~~~~

- vcpkg + gsd are now submodules so all dependencies are tracked by dependabot!
- ``Property``'s ``array_ref_t`` is now an ``Eigen::Map`` to disallow resizing of ``SoA``'s individual arrays.

Fixes
~~~~~~~~~~~~~~

- SoA resize was broken

Removed
~~~~~~~~~

Version 0.2.0
--------------------------------

The first released alpha version of openFLY! This is a minimal feature-set release with just the base classes that underpin openFLY as well as binary IO.

Added
~~~~~~~~~

- Binary IO using the GSD format through the ``FileGSD`` class.

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

- Hdoc is no longer used to build the documentation.

Version 0.1.0-alpha
---------------------------

Continuous pre-release, not currently in a usable state.