Local environments
======================


Directory: ``libfly/env``

The ``env`` folder contains routines working with local environments. Everything is contained within the namespace ``fly::env``. The local environment (LE) of an atom is the set of atoms within some neighbourhood.


Geometry
----------------

File: ``libfly/env/geometry.hpp``

.. doxygenfile:: libfly/env/geometry.hpp
    :sections: briefdescription detaileddescription


Functions
~~~~~~~~~~~

.. doxygenfunction:: fly::env::centroid

.. doxygenfunction:: fly::env::rmsd

.. doxygenfunction:: fly::env::grmsd

.. doxygenfunction:: fly::env::ortho_onto

.. doxygenfunction:: fly::env::for_equiv_perms


Geometry class
~~~~~~~~~~~~~~~~~~~~~~~

.. doxygenclass:: fly::env::Geometry
    :members:
    :undoc-members:


Local 
--------------------------

File: ``libfly/env/local.hpp``

.. doxygenfile:: libfly/env/local.hpp
    :sections: briefdescription detaileddescription

Fingerprints
~~~~~~~~~~~~~~~~~~

.. doxygenclass:: fly::env::Fingerprint
    :members:
    :undoc-members:


Local class
~~~~~~~~~~~~~~~~~~

.. doxygenclass:: fly::env::Local
    :members:
    :undoc-members:

LocalList class
~~~~~~~~~~~~~~~~~~

.. doxygenclass:: fly::env::LocalList
    :members:
    :undoc-members:

