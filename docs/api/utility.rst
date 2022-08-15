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

Configuration macros
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. doxygendefine:: FLY_SPATIAL_DIMS

Error handling
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. doxygenstruct:: fly::RuntimeError
    :members:
    :undoc-members:

.. doxygenfunction:: error

.. doxygenfunction:: verify

.. doxygendefine:: ASSERT


Defines, variables, etc.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. doxygenvariable:: spatial_dims

.. doxygenenum:: Sign

.. doxygentypedef:: Vec

.. doxygentypedef:: Mat

.. doxygentypedef:: Arr

.. doxygentypedef:: first_t

.. doxygentypedef:: remove_cref_t

.. doxygenvariable:: always_false


Small functions
~~~~~~~~~~~~~~~~~~~~~~

.. doxygenfunction:: safe_cast

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

.. doxygenclass:: fly::Vector
    :members:
    :undoc-members:
    
.. doxygenclass:: fly::Defer
    :members:
    :undoc-members:

Timing 
~~~~~~~~~~~~~~~~~~~~~~~~

File ``timeit.hpp``

.. doxygenfunction:: timeit


Random numbers 
----------------------

File ``random.hpp``

.. doxygenfile:: random.hpp
    :sections: briefdescription detaileddescription


.. doxygenclass:: fly::Xoshiro
    :members:
    :undoc-members:


