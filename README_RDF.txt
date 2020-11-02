RDF branch
11/2/2020   Rob Riggleman

This branch calculates g(r) for the specified particle types. 

It requires more command-line arguments to function:

./postproc-rdf [input.lammpstrj] [frame 0] [final frame] [bin size] [type A]
[type B]

The first three arguments are commonplace.

bin size = the resolution with which to calculate g(r)

type A: one of the types involved in the g(r)
type B: other type involved. 

Type A can be the same as type B, and both can be set to "all"

