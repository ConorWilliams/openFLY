Minimisers
======================================

Directory: ``libfly/minimise``

The ``minimise`` folder contains a collection of potential energy minimisers. Everything is contained within the namespace ``fly::minimise``.


LBFGS
----------------

Newton stepper
~~~~~~~~~~~~~~~~~~~~

File: ``libfly/minimise/LBFGS/lbfgs_core.hpp``

.. doxygenfile:: libfly/minimise/LBFGS/lbfgs_core.hpp
    :sections: briefdescription detaileddescription


.. doxygenclass:: fly::minimise::StepLBFGS
    :members:
    :undoc-members:


Minimiser
~~~~~~~~~~~~~

File: ``libfly/minimise/LBFGS/lbfgs.hpp``

.. doxygenfile:: libfly/minimise/LBFGS/lbfgs.hpp
    :sections: briefdescription detaileddescription


.. doxygenclass:: fly::minimise::LBFGS
    :members:
    :undoc-members:
