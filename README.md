# Introduction to Scientific Computing and Programming in Chemistry
Scripts developed in the course of Introduction to Scientific Computing and Programming in Chemistry at the Institute of Chemistry at the University of São Paulo (IQ-USP)

## z_matrix.py
Python script to convert xyz cartesian coordinates in z-matrix coordinates

How to use:

./z_matrix.py xyz_file.xyz

This script will read the xyz coordinates from the xyz_file and write a new file with the same name but with .zmat extension and z-matrix coordinates.

## Spin_coupling.sh
Bash script for extracting spin-spin coupling constants from a Gaussian program output file.

How to use:

It takes at least 4 arguments, the first will be the ouput file, the second will be the constant type (KFC, JFC, KSD, JSD, KDSO, JDSO, KPSO, JDSO, KTOTAL, JTOTAL), the third will be the number of first atom of the coupling, from the fourth onwards will be the numbers of the other atoms that the first one will couple.
Example:

./Spin_Coupling.sh output_file.log JDSO 14 19 22 1 7
This will give the JDSO spin-spin coupling constant to the pairs 14 19, 14 22, 14 1, and 14 7.

## boltzman.pl
Perl script to calculate the energy differences between different conformations of the same molecule and the boltzman populations for each conformation from the outputs generated by the geometry and frequency optimization calculation generated by the Gausssian 16 program.

How to use:

./boltzman.pl conformer1.log conformer2.log conformer3.log ...

The energies diferences and boltzman populations will be calculated and printed in screen 

## geom.f90
Source code writed in fortran90 to calculate the atoms distances from a xyz file and printing just the distances between atoms that are bonded

How to use:

Copile the source code with some compiler. Example: gfrotran -o geom.x geom.f90
Run the executable geom.x 
Type the name of the xyz file that the distances will be calculated

This script will read the xyz coordinates from the xyz_file and print the distances in terminal.

## geom.pl
Perl script to calculate the atoms distances from a xyz file.

How to use:

./geom.pl xyz_file.xyz

This script will read the xyz coordinates from the xyz_file and print the distances in terminal.
