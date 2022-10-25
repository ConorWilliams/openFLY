Environments
======================


Directory: ``libfly/env``

The ``env`` folder contains routines working with local environments. Everything is contained within the namespace ``fly::env``. The local environment (LE) of an atom is the set of atoms within some neighbourhood. It is assumed in OLKMC that the mechanisms accessible to an atom are completely contained-within and solely a-function-of its LE.


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


Heuristics
--------------------------

File: ``libfly/env/heuristics.hpp``

.. doxygenfile:: libfly/env/heuristics.hpp
    :sections: briefdescription detaileddescription

Fingerprints
~~~~~~~~~~~~~~~~~~

.. doxygenclass:: fly::env::Fingerprint
    :members:
    :undoc-members:


Graph hash
~~~~~~~~~~~~~~~~~~

.. doxygenfunction:: fly::env::canon_hash


Mechanisms
---------------------------

File: ``libfly/env/mechanisms.hpp``

.. doxygenfile:: libfly/env/mechanisms.hpp
    :sections: briefdescription detaileddescription

.. doxygenclass:: fly::env::Mechanism
    :members:
    :undoc-members:


The catalogue
--------------------------

File: ``libfly/env/catalogue.hpp``

.. doxygenfile:: libfly/env/catalogue.hpp
    :sections: briefdescription detaileddescription

.. doxygenclass:: fly::env::Catalogue
    :members:
    :undoc-members: