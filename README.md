depletion
=========

Date of Creation: 18 April 2013.
Authors: Kevin Manalo (GT), Michael Chin (GT)

Purpose: Creation of modules used in development of openmc-depletion.

Compilers: nagfor, ifort

Subroutines needed: 

Compressed Row Storage Insert subroutine - A subroutine that can insert data into arbitrary positions of a CSR matrix.
Symbolic Gaussian Elimination Vertex Fill algorithm subroutine - A subroutine that adds fill verticies to a CSR matrix.
      This subroutine has dependencies: Compressed Row Storage Insert Subroutine
