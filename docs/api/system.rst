System
====================================

Directory: ``libfly/system``

The ``system`` folder contains a collection of data structures for representing atoms, groups of atoms and the space the atoms exist in. Everything is contained within the namespace ``fly::system``.

Simulation space
--------------------------------

File: ``libfly/system/box.hpp``

.. doxygenfile:: libfly/system/box.hpp
    :sections: briefdescription detaileddescription

In theory libFLY could detect the most efficient crystal-system for any basis set and produce optimal code. Unfortunately, this is a lot of work hence, libFLY uses a smaller set - suitable for the majority of cases. If you would like to read more about these specialised boxes see:

.. toctree::
   
   crystal_systems.rst

Generalised box
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. doxygenclass:: fly::system::Box
    :members:
    :undoc-members:


Atoms
----------------------

File: ``libfly/system/atom.hpp``

.. doxygenfile:: libfly/system/atom.hpp
    :sections: briefdescription detaileddescription

Atom class
~~~~~~~~~~

.. doxygenstruct:: fly::system::Atom
    :members:
    :undoc-members:


MemTag
~~~~~~

.. doxygenstruct:: fly::system::MemTag
    :members:
    :undoc-members:


Built-in members
~~~~~~~~~~~~~~~~~~

.. doxygennamespace:: fly::builtin_m
  

Generic data structures for Atoms
-----------------------------------

The building block of many of the types in libFLY.

SoA
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

File: ``libfly/system/SoA.hpp``

.. doxygenfile:: libfly/system/SoA.hpp
    :sections: briefdescription detaileddescription

.. doxygenclass:: fly::system::SoA
    :members:
    :undoc-members:

VoS 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

File: ``libfly/system/VoS.hpp``

.. doxygenfile:: libfly/system/VoS.hpp
    :sections: briefdescription detaileddescription

.. doxygenclass:: fly::system::VoS
    :members:
    :undoc-members:


