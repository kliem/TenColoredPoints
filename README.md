# TenColoredPoints
Scripts proving that the optimal colored Tverberg problem holds for 10 points in the plane.

Given 10 points in general position in the plane.
They induce an acyclic chirotope of rank 3 on those 10 points that is
nowhere zero.
We construct a `k`-partite graph according to
- Kliem, Jonathan: A new k-partite graph k-clique iterator and the optimal colored Tverberg problem for ten colored points. Preprint (2021).
A counter example of the optimal colored Tverberg problem would imply
the existence of k-clique in such a graph.

We use the `k`-clique iterator
- https://github.com/kliem/KPartiteKClique
to verify whether it has a `k`-clique.

Iterating over all such chirotopes with the python package
- `pseudo_order_types`
this package verifies the optimal colored Tverberg problem for
10 points in the plane.

Alternatively, we can iterate with the list
- http://www.ist.tugraz.at/staff/aichholzer/research/rp/triangulations/ordertypes/
over all such realizable chirotopes.

## Dependencies

- `Cython`
- `memory_allocator`
