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

.. Version is specified here, in vcpkg.json, docs/index.rst and libfly/utility/version.hpp

Unreleased
-------------------------------
Added
~~~~~
Changed
~~~~~~~
Removed
~~~~~~~
Bugfixes
~~~~~~~~
Meta 
~~~~


Version 0.9.0
-------------------------------
Added
~~~~~

- New ``viewSoA`` alias.
- ``Supercell`` has a ``set_box`` method.
- New ``SKMC`` class for running OLKMC simulations.
- New ``SuperBasin`` class for managing a collection of basins.
- New ``SuperCache`` class for managing a collection of superbasins.
- New ``min_delta_max`` option for ``Catalogue``.
- ``Dimer`` has a ``use_history`` option.
- ``Dimer``'s ``find_sp`` accepts the number of frozen atoms.
- ``Hessian`` exposes a ``Matrix`` member type and methods to calculate eigen vectors.
- New ``centoid`` and ``centroid_align`` functions in the lattice namespace.
- New ``DetectVacacies`` class for explicitly detecting vacancies.
- ``Catalogue`` counts false positives and has an ``optimise`` method.

Changed
~~~~~~~


- ``Dimer`` fixes the COM of the cell during SP searches by projection translational components out.
- During SP searches reconstruction we deterministically chose the dimer to axis to improve reproducibility.
- During SP searches we utilise the fixed COM in the SP comparison step.
- Potential's hessian function is specified as non-mass weighted and a mass-weight function is included.
- Increased default ``MAX_GHOST_RATIO`` for ``neigh::List``.
- Changed catalogue to write in Cereal's portable_binary format.
- Split poisoned mechanisms into SP/Minima poisoned.
- Renamed ``xise`` to ``ssize``!
- New options for EAM parsing to support debugging and alternative tabulations.


Bugfixes
~~~~~~~~

- EAM Hessian bug fixed (was calculating slightly wrong Hessians!).
- EAM potential parser handles tabulations with ``numP != numR`` correctly.
- EAM potential does not require the same number of types in-use/available just that in-use is a subset of available. 

Meta 
~~~~

- Cmake's DOWNLOAD_EXTRACT_TIMESTAMP Warnings fixed.
- Removed CI macOS builds due to openMP bug on on GitHub actions runners.

Version 0.8.0
-------------------------------

First release to actually do any KMC! This release introduces the ``Basin`` class for performing N-fold way KMC as well as: serialisation support for the catalogue; an adaptive catalogue and a ``Master`` s API changed (again!).

Added
~~~~~

- New ``Basin`` class for performing KMC simulations
- New lattice functions: ``add_atoms()``, ``motif_to_lattice()``, ``remove_atoms()``, ``remove_sphere()``.
- New ``refine_tol()`` function in Catalogue.
- New ``reconstruct()`` method in Catalogue.
- New ``dprint`` utility function.
- Catalogue has serialisation support (hence all serialised members of the catalogue do to).
- Catalogue environments expose their delta_max.
- Catalogue has a method to calculate the (approximate) symmetries a geometry has.
- Dimer::Exit has a new return code.
- Many new options in Master::Options.
- New Xoshiro constructor to properly seed from ``std::random_device``.
- Box has a new min_image calculator method.

Changed
~~~~~~~

- ``template_for()`` enhanced.
- Catalogue now sets delta_max of environments as minimum of Options.delta_max and the appropriate function of Fingerprint.r_min() enabling a more general default.
- Catalogue environment's ``set_mech`` function moved to a Catalogue method.
- Mechanism's members tweaked to encode more reconstruction info.
- Hessians now contain zeros for frozen atoms instead of only being n_active x n_active to facilitate future conditioning and make them more like gradients.
- Master's default option mech_tol increased.
- Master's API has changed to facilitate new features + MPI.
- Changed the base type of Hash

Bugfixes
~~~~~~~~
- Typo in SoA::operator=
- SoA::rebind made const correct.
- Environments in the catalogue are allowed to have no mechanisms.
- Catalogue works correctly on the second pass (memoization bug). 

Meta 
~~~~
- New dependency Cereal
- Now uses ccache in developer mode.
- Fixed openMP compilation in MacOS cloud.

Version 0.7.0
-------------------------------

Overhauled saddle-point finding including: returning the new ``Mechanism`` class, automatic mechanisms symmetry identification and discovery and history dependant dimer searches.

Added
~~~~~

- New ``Mechanism`` class.
- New ``rebuild_geo_from_nl()`` geometry function.
- The ``Master`` class added which is an improved version of the ``MasterFinder`` class.

Changed
~~~~~~~

- ``Dimer``'s interface changed.
- Generalised some geometry functions: ``rmsd()``, ``grmsd()``.
- Potentials now compute mass weighted Hessians.
- EAM hessian computation made more cache efficient.

Removed
~~~~~~~

- The ``MasterFinder`` class was removed.

Bugfixes
~~~~~~~~

- Return wrong value in ``examples/env/geometry.cpp``.
- Fix EAM parser bug not reading masses correctly.
- Fix bug in finder which forgot to ``std::swap`` two states.



Version 0.6.0
-------------------------------

This release introduces local environments, the catalogue and implementation of our invariant and tolerant matching algorithm.


Added
~~~~~

- New ``Catalogue`` class.
- New ``Fingerprint`` class.
- New ``canon_hash()`` function.
- New geometry functions: ``centroid()``, ``rmsd()``, ``grmsd()``, ``for_equiv_perms()``.
- New ``ortho_onto()`` function.
- New ``Colour`` property.
- New ``Geometry`` class.
- New property ``Hash``.
- New meta programming utility ``is_narrowing_conversion_v``.
- ``Atom`` has a default constructor.
- ``VoS::atom_t`` exposes the underlying atom type. 
- New internal graph class.

Meta 
~~~~
- New dependencies Nauty and xxHash


Version 0.5.0
-------------------------------

This release introduce saddle-point finding and min->sp->min pathway finding. The concept of a generic potential was made more concrete to prevent a template explosion.


Added
~~~~~

- New ``Dimer`` saddle-point finder.
- New ``Rotor`` class.
- New ``perturb()`` function
- New ``MasterFinder`` class.
- ``SoA`` has a new rebind method. 

Changed
~~~~~~~

- Generalised ``StepLBFGS``'s ``.newton_step()``.
- ``Generic`` potential API + constructor changes
- Unified minimiser, saddle finder and dimer return codes to follow C conventions (truthy on failure);
- ``Spline`` methods clamp input.

Removed
~~~~~~~

- ``LBFGS`` no longer has special handling for dimer classes. 

Bugfixes
~~~~~~~~

- Const-corrected ``Generic::gradient``.
- Padded spline with terminator to fix-up floating point rounding errors.


Version 0.4.0
--------------

This release introduces generic potentials and the first concrete potential into openFLY, EAM. The EAM implementation includes support for analytic Hessians and is fully openMP parallelised. Additionally, an efficient parallel implementation of the LBFGS minimiser is included.

Added
~~~~~

- New ``Delta`` property.
- New ``StepLBFGS`` class.
- New ``Spline`` utility class.
- New ``DataEAM`` class with eam/fs parsing.
- New ``potential::Generic`` class.
- New ``EAM`` class.
- New ``xize`` utility function.
- New ``LBFGS`` class.
- New ``Hessian`` class.
- ``Frozen`` property has a tag to enable GSD IO.

Changed
~~~~~~~

- ``fly::near()`` now has customizable tolerances.
- ``neigh::List``'s ``update()`` API changed.
- ``SoA``'s converting constructors are now SFINE friendly.
- ``SoA``'s base classes are now public.
- ``TypeMap``'s converting constructor is now SFINE friendly.
- De-generalised ``SoA``'s converting constructors to allow implicit casts.
- Read methods on ``BinaryFile`` are ``const``.
- ``LBFGS`` force tolerance default tightened.
- ``Spline`` clamps interval.

Bugfixes
~~~~~~~~

- Box (Ortho and Triclinic, valid bounds now include zero).

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
- Many documentation enhancements.

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

Meta 
~~~~~~~~~~~~~~~~~~~~~~~~~

- Hdoc is no longer used to build the documentation.

Version 0.1.0 pre-release
---------------------------

Continuous pre-release, not currently in a usable state.