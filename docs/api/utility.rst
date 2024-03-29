Utility
======================================

Directory: ``libfly/utility``

The ``utility`` folder contains a collection of non-specific utilities, it is documented by file.

Versioning
---------------------------------------

File: ``libfly/version.hpp``

.. doxygenfile:: version.hpp

Core 
------------------------------------------

File: ``libfly/core.hpp``

.. doxygenfile:: core.hpp
    :sections: briefdescription detaileddescription

Configuration macros
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. doxygendefine:: FLY_SPATIAL_DIMS

Error handling
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. doxygenstruct:: fly::RuntimeError
    :members:
    :undoc-members:

-------------------

.. doxygenfunction:: error

-------------------

.. doxygenfunction:: verify

-------------------

.. doxygendefine:: ASSERT


Defines, variables, etc.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. doxygenvariable:: spatial_dims

-------------------

.. doxygenenum:: Sign

-------------------

.. doxygentypedef:: Vec

-------------------

.. doxygentypedef:: Mat

-------------------

.. doxygentypedef:: Arr

-------------------

Meta programming
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. doxygentypedef:: first_t

-------------------

.. doxygentypedef:: remove_cref_t

-------------------

.. doxygenvariable:: always_false

-------------------

.. doxygenvariable:: is_detected_v

-------------------

.. doxygentypedef:: detected_or_t

-------------------

.. doxygenvariable:: is_narrowing_conversion_v

Small functions
~~~~~~~~~~~~~~~~~~~~~~

.. doxygenfunction:: dprint

-------------------

.. doxygenfunction:: safe_cast

-------------------

.. doxygenfunction:: visit

-------------------

.. doxygenfunction:: template_for

-----------------------

.. doxygenfunction:: fly::ssize

Mathematical functions
~~~~~~~~~~~~~~~~~~~~~~~~

.. doxygenfunction:: near

-------------------

.. doxygenfunction:: product_scan

-------------------

.. doxygenfunction:: ipow

-------------------

.. doxygenfunction:: gdot

-------------------

.. doxygenfunction:: gnorm

-------------------

.. doxygenfunction:: gnorm_sq

-------------------

.. doxygenfunction:: hyperplane_normal


Classes
~~~~~~~~~~~~~~~~~~~~~~~~
 
.. doxygenclass:: fly::Defer
    :members:
    :undoc-members:

Timing 
~~~~~~~~~~~~~~~~~~~~~~~~

File: ``libfly/timeit.hpp``

.. doxygenfunction:: timeit


Random numbers 
----------------------

File ``libfly/random.hpp``

.. doxygenfile:: random.hpp
    :sections: briefdescription detaileddescription

.. doxygenclass:: fly::Xoshiro
    :members:
    :undoc-members:

Natural Splines
---------------------

File ``libfly/spline.hpp``

.. doxygenfile:: spline.hpp
    :sections: briefdescription detaileddescription

.. doxygenclass:: fly::Spline
    :members:
    :undoc-members:

Lattices
---------------------

File ``libfly/lattice.hpp``

.. doxygenfile:: lattice.hpp
    :sections: briefdescription detaileddescription

Building/manipulating
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. doxygenfunction:: motif_to_lattice

.. doxygenfunction:: add_atoms
    
.. doxygenfunction:: remove_atoms
    
.. doxygenfunction:: remove_sphere

.. doxygenfunction:: fly::centroid

.. doxygenfunction:: centroid_align

Vacancies
~~~~~~~~~~

.. doxygenclass:: fly::DetectVacancies
    :members:
    :undoc-members:

.. doxygenfunction:: kruskal_max