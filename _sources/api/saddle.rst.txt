Saddle-point searching
======================================

Directory: ``libfly/saddle``

The ``saddle`` folder contains a collection of utilities for saddle-point finding. Everything is contained within the namespace ``fly::saddle``.

A (strictly first-order) saddle-point is a stationary point where the Hessian matrix has exactly one negative eigenvalue.


Dimer
---------------------

File: ``libfly/saddle/dimer.hpp``

.. doxygenfile:: libfly/saddle/dimer.hpp
    :sections: briefdescription detaileddescription

.. doxygenclass:: fly::saddle::Dimer
    :members:
    :undoc-members:



Rotor
---------------------

File: ``libfly/saddle/rotor.hpp``

.. doxygenfile:: libfly/saddle/rotor.hpp
    :sections: briefdescription detaileddescription

.. doxygenclass:: fly::saddle::Rotor
    :members:
    :undoc-members:


Perturb
------------------

File: ``libfly/saddle/perturb.hpp``

.. doxygenfile:: libfly/saddle/perturb.hpp
    :sections: briefdescription detaileddescription

.. doxygenfunction:: fly::saddle::perturb



