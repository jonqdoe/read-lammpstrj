c++ implementation of read_lammpstrj
10/28/2020          R. Riggleman

This code uses classes (strings/vectors) to read in a lammpstrj file. It is
more general than previous implementations in that it can deal with arbitrary
ordering of columns and data other than "x", "y", "z". 

It assumes that id, type, mol, x, y, z, if present, will all come as the first
six columns, and that any computes will follow the z coordinate. 

To do:
- include the neighborlist functionality

- include particle-to-mesh functionality
