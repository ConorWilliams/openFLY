System
====================================

Directory: ``libfly/system``

The ``system`` folder contains a collection of data structures for representing atoms, groups of atoms and the space the atoms exist in. Everything in the systems folder is contained within the namespace ``fly::system``.

Atoms
----------------------

File: ``libfly/system/boxes/atom.hpp``

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


Bultin members
~~~~~~~~~~~~~~~

.. doxygennamespace:: fly::builtin_m
  

Generic data structures for Atoms
-----------------------------------

The building block of many of the types in libFLY.

SoA
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

File: ``libfly/system/boxes/SoA.hpp``

.. doxygenfile:: libfly/system/SoA.hpp
    :sections: briefdescription detaileddescription

.. doxygenclass:: fly::system::SoA
    :members:
    :undoc-members:

VoS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

File: ``libfly/system/boxes/VoS.hpp``

.. doxygenfile:: libfly/system/VoS.hpp
    :sections: briefdescription detaileddescription

.. doxygenclass:: fly::system::VoS
    :members:
    :undoc-members:


Simulation space
--------------------------------

File: ``libfly/system/boxes/box.hpp``

.. doxygenfile:: libfly/system/box.hpp
    :sections: briefdescription detaileddescription

If you would like to read about the specialised crystal systems use the links below:

.. toctree::
   
   crystal_systems.rst

