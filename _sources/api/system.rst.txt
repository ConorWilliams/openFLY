System
====================================

Directory: ``libfly/system``

The ``system`` folder contains a collection of data structures for representing atoms, groups of atoms and the space the atoms exist in. Everything is contained within the namespace ``fly::system``.

Simulation space
--------------------------------

File: ``libfly/system/box.hpp``

.. doxygenfile:: libfly/system/box.hpp
    :sections: briefdescription detaileddescription

Generalised box
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. doxygenclass:: fly::system::Box
    :members:
    :undoc-members:

Specialised boxes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In theory libFLY could detect the most efficient crystal-system for any basis-set and produce optimal code. Unfortunately, this is a lot of work hence, libFLY uses a smaller set - suitable for the majority of cases. If you would like to read more about these specialised boxes see:

.. toctree::
   
   crystal_systems.rst


Properties
------------------------

File: ``libfly/system/property.hpp``

.. doxygenfile:: libfly/system/property.hpp
    :sections: briefdescription detaileddescription


.. doxygenstruct:: fly::system::Property
    :members:
    :undoc-members:

Built-in 
~~~~~~~~~~~~~~~~~~

.. doxygennamespace:: fly::builtin_properties



Property generic data-structures
----------------------------------------

Atom class
~~~~~~~~~~

File: ``libfly/system/atom.hpp``

.. doxygenfile:: libfly/system/atom.hpp
    :sections: briefdescription detaileddescription

.. doxygenstruct:: fly::system::Atom
    :members:
    :undoc-members:

Type map
~~~~~~~~~~

File: ``libfly/system/typemap.hpp``

.. doxygenfile:: libfly/system/typemap.hpp
    :sections: briefdescription detaileddescription

.. doxygenclass:: fly::system::TypeMap
    :members:
    :undoc-members:


Structure of arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

File: ``libfly/system/SoA.hpp``

.. doxygenfile:: libfly/system/SoA.hpp
    :sections: briefdescription detaileddescription

.. doxygenclass:: fly::system::SoA
    :members:
    :undoc-members:

Array of structures
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

File: ``libfly/system/VoS.hpp``

.. doxygenfile:: libfly/system/VoS.hpp
    :sections: briefdescription detaileddescription

.. doxygenclass:: fly::system::VoS
    :members:
    :undoc-members:


Supercells
---------------------------

File: ``libfly/system/supercell.hpp``

.. doxygenfile:: libfly/system/supercell.hpp
    :sections: briefdescription detaileddescription

.. doxygenclass:: fly::system::Supercell
    :members:
    :undoc-members:

.. doxygenfunction:: fly::system::make_supercell



Hessians
-----------------------

File: ``libfly/system/hessian.hpp``

.. doxygenfile:: libfly/system/hessian.hpp
    :sections: briefdescription detaileddescription

.. doxygenclass:: fly::system::Hessian
    :members:
    :undoc-members:
