System (``libfly/system``)
====================================

The ``system`` folder and namespace contains a collection of data structures for representing atoms, groups of atoms and the space the atoms exist in.

Atoms (``atom.hpp``)
----------------------

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
  

Generic data structures
-----------------------

Structure of Arrays (``SoA.hpp``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. doxygenclass:: fly::system::SoA
    :members:
    :undoc-members:

Vector of Structures (``VoS.hpp``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. doxygenclass:: fly::system::VoS
    :members:
    :undoc-members:


Simulation box (``box.hpp``)
--------------------------------

.. doxygenclass:: fly::system::AdjacentCells
    :members:
    :undoc-members:

Crystal Systems 
~~~~~~~~~~~~~~~~~~~~

Orthorhombic 
````````````````````````````````````````````

.. doxygenclass:: fly::system::Orthorhombic
    :members:
    :undoc-members:

Triclinic 
````````````````````````````````````````````