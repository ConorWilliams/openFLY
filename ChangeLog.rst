Changelog
============================

.. Unreleased
.. -------------------------------
.. Added
.. ~~~~~
.. Changed
.. ~~~~~~~
.. Removed
.. ~~~~~~~
.. Bugfixes
.. ~~~~~~~~
.. Meta 
.. ~~~~


Unreleased
----------
Added
~~~~~

- New ``Delta`` property.
- New ``StepLBFGS`` class.

Changed
~~~~~~~

- ``neigh::List``'s ``update()`` API changed

Removed
~~~~~~~
Bugfixes
~~~~~~~~
Meta 
~~~~

Version 0.3.0
------------------------

The second alpha version of openFLY, this release brings neighbour-list support to libFLY.

Added
~~~~~~~~~

- Neighbour-list support via the ``neigh::List`` class.
- Internal ``Vector`` class to replace ``std::vector``.
- ``Xoshiro`` class as the PRNG.
- ``Box`` has a new ``get()`` method to fetch the underlying ``std::variant``.
- ``visit`` utility function.
- ``neighbour::sort`` function to optimise ordering for neighbour operations.
- ``operator=`` for the ``Atom`` class
- New ``template_for`` utility function.
- Added ``min_width`` member to crystal specialisations.

Changed
~~~~~~~~~~

- ``Property``'s ``array_ref_t`` is now an ``Eigen::Map`` to disallow resizing of ``SoA``'s individual arrays.
- ``SoA``'s ``destructive_resize`` now returns a boolean.
- Allow zero length ``SoA`` s.
- ``BinaryFile`` internals reworked to support clang.
- ``VoS`` uses the ``Vector`` class.

Removed
~~~~~~~~~

- Removed the ``Orthorombic``'s deprecated member min-image. 
- Removed ``Box`` default constructor.

Bugfixes
~~~~~~~~~~~~~~

- Crystal systems ``gen_image`` marked ``const``.
- Fixes to test and examples that assumed 3D. 
- ``SoA``'s ``operator()`` was broken.
- ``SoA``'s ``resize`` was broken.


Meta 
~~~~~~~~~~~~~~~~~~~~~~~~~

- vcpkg + gsd are now submodules so all dependencies are tracked by Dependabot!
- Many documentation engancements.

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

- GPL-2.0 -> GPL-3.0-or-later.

Removed
~~~~~~~~~

Meta 
~~~~~~~~~~~~~~~~~~~~~~~~~

- Hdoc is no longer used to build the documentation.

Version 0.1.0 pre-release
---------------------------

Continuous pre-release, not currently in a usable state.