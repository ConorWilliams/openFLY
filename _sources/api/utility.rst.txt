Utility
======================================

The ``utility`` folder contains a collection of non-specific utilities, it is documented by file.

Versioning
---------------------------------------

File: ``version.hpp``

.. doxygenfile:: version.hpp

Core 
------------------------------------------

File: ``core.hpp``

.. doxygenfile:: core.hpp
    :sections: briefdescription detaileddescription

Macros
~~~~~~~~~

.. doxygendefine:: ASSERT

.. doxygendefine:: VERIFY    

.. doxygendefine:: FLY_SPATIAL_DIMS


Defines, variables, etc.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. doxygenvariable:: spatial_dims

.. doxygenvariable:: max_atomic_num

.. doxygenenum:: Sign

.. doxygentypedef:: Vec

.. doxygentypedef:: Arr

.. doxygentypedef:: Mat

.. doxygentypedef:: first_t

.. doxygentypedef:: remove_cref_t

.. doxygenvariable:: always_false


Mathematical functions
~~~~~~~~~~~~~~~~~~~~~~~~

.. doxygenfunction:: near

.. doxygenfunction:: product_scan

.. doxygenfunction:: ipow

.. doxygenfunction:: gdot

.. doxygenfunction:: gnorm

.. doxygenfunction:: gnorm_sq

.. doxygenfunction:: hyperplane_normal


Classes
~~~~~~~~~~~~~~~~~~~~~~~~


.. doxygenclass:: fly::Defer
    :members:
    :undoc-members:

Timing 
------------------------------------------

File ``timeit.hpp``

.. doxygenfile:: timeit.hpp
    :sections: briefdescription detaileddescription

.. doxygenfunction:: timeit

