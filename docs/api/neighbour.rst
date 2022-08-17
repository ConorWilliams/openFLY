Neighbour
======================================

Directory: ``libfly/neighbour``

The ``neigh`` folder contains a collection of data structures and algorithms for generating and manipulating *neighbours*. Neighbours of an atom are atoms withing some cut-off radius, ``r_cut``. Everything is contained within the namespace ``fly::neigh``.

Lists
----------------

File: ``libfly/neigh/box.hpp``

.. doxygenfile:: libfly/neigh/list.hpp
    :sections: briefdescription detaileddescription

.. doxygenclass:: fly::neigh::List
    :members:
    :undoc-members: