#!/bin/python3
#
# matriz_z.py
#
# ************************************************************************************************
# *                        Conversion xyz coordinates to z-matrix coordinates                    *
# ************************************************************************************************
#
# https://github.com/danielarley/Introduction-to-scientific-computing-and-programming-in-chemistry
#
# Script to convert xyz cartesian coordinates in z-matrix coordinates
#
# How to use:
# ./matriz_z.py xyz_file.xyz
#
# This script will read the xyz coordinates from the xyz_file and write a new file with the same name
# but with .zmat extension and z-matrix coordinates.
#

import numpy as np
import os
import argparse
import math

### Functions ###

def open_xyzfile(file):
    """ Open .xyz file and stores symbols and coordinates in arrays.
    
    Parameters:
    file - xyz file with cartesian coordinates of the system
    
    Return:
    symbols - Element symbols (1N array)
    coordinates - the coordinates of each atom (ND array, where N= number of atoms)
    num_atoms - number of atoms in the system (int)
    """
    xyz_file = np.genfromtxt(fname=file,skip_header=2,dtype='unicode') # Skips the first two lines and opens the file like an N-dimensional array
    symbols = xyz_file[:,0] # Store the atoms symbols
    coordinates = (xyz_file[:,1:]) # each array dimension represents the coordinates of an atom 
    coordinates = coordinates.astype(np.float64)
    num_atoms=len(symbols)
    return symbols, coordinates, num_atoms

def calculate_distance(atom1_coord, atom2_coord):
    """ Calculate the distance between two atoms.
    
    Parameters:
    atom1_coord - xyz coordinates of atom1 (1N array)
    atom2_coord - xyz coordinates of atom2 (1N array)
    
    Return:
    bond_lenght - Distance between atom1 and atom2 (float)
    """
    delta_x = atom1_coord[0] - atom2_coord[0]
    delta_y = atom1_coord[1] - atom2_coord[1]
    delta_z = atom1_coord[2] - atom2_coord[2]
    bond_length = math.sqrt(delta_x**2+delta_y**2+delta_z**2)
    return bond_length

def calculate_parameter(atom1_symb,atom2_symb):
    """ Calculates the parameter that is used to define wheter atom1 and atom2 are bonded.
    
    Parameters:
    atom1_symb - atom1 symbol (str)
    atom2_symb - atom2 symbol (str)
    
    Return:
    parameter - parameter that define wheter atom1 and atom2 is bonded (float)
    """
    elem=["H","He","Li","Be","B","C","N","O","F","Ne","Na",\
          "Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti",\
          "V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As",\
          "Se","Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru",\
          "Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs",\
          "Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy",\
          "Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt",\
          "Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th",\
          "Pa","U","Np","Pu","Am","Cm"]
    cov_radii=[0.310,0.280,1.280,0.960,0.840,0.760,0.710,0.660,0.570,\
               0.580,1.660,1.410,1.210,1.110,1.070,1.050,1.020,1.060,\
               2.030,1.760,1.700,1.600,1.530,1.390,1.500,1.420,1.380,\
               1.240,1.320,1.220,1.220,1.200,1.190,1.200,1.200,1.160,\
               2.200,1.950,1.900,1.750,1.640,1.540,1.470,1.460,1.420,\
               1.390,1.450,1.440,1.420,1.390,1.390,1.380,1.390,1.400,\
               2.440,2.150,2.070,2.040,2.030,2.010,1.990,1.980,1.980,\
               1.960,1.940,1.920,1.920,1.890,1.900,1.870,1.870,1.750,\
               1.700,1.620,1.510,1.440,1.410,1.360,1.360,1.320,1.450,\
               1.460,1.480,1.400,1.500,1.500,2.600,2.210,2.150,2.060,\
               2.000,1.960,1.900,1.870,1.800,1.690]
    for k in range(0,95):
        if elem[k] == str(atom1_symb):
            radii1=cov_radii[k]
        if elem[k] == str(atom2_symb):
            radii2=cov_radii[k]
    vdw_sum=radii1+radii2
    parameter=vdw_sum*1.1
    return parameter

def bond_check(atom_distance,parameter):
    """ Check if a distance is a bond based on sum of vdw radii.
    
    Parameters:
    atom_distance - distance between two atoms (float)
    parameters - Parameters that defines a bond: vdw_sum*1.3 (float)
    
    Return: 
    True - If less or equal than parameter (logical)
    False - If greater than parameter or zero (logical)
    """
    if atom_distance > 0 and atom_distance <= parameter:
        return True
    else:
        return False

def calculate_angle_v1_v2(v1,v2):
    """ Calculate the angle between two three-dimensional vectors
    
    Parameters:
    v1- vector 1 (1N array)
    v2- vector 2 (1N array)
    
    Return: 
    angle_v1_v2 - angle between v1 and v2 (float)
    """
    scalar_v1_v2=np.dot(v1,v2) #Scalar product
    mod_v1=np.sqrt(np.dot(v1,v1)) # v1 module
    mod_v2=np.sqrt(np.dot(v2,v2)) # v2 module
    cos_ang=scalar_v1_v2/(mod_v1*mod_v2) # cosine of the angle
    angle_v1_v2=math.acos(cos_ang)*57.2958 # Angle calculate
    return angle_v1_v2

def calculate_angle(atom1_coord,atom2_coord,atom3_coord):
    """"Calculate the angle between three atoms.
    
    Parameters:
    atom1_coord - xyz coordinates of atom1 (1N array)
    atom2_coord - xyz coordinates of atom2 (1N array)
    atom3_coord - xyz coordinates of atom3 (1N array)
    
    Return:
    angle - angle between atom1, atom2, and atom3 (float)
    """
    v1=[atom2_coord[0] - atom1_coord[0],atom2_coord[1] - atom1_coord[1],atom2_coord[2] - atom1_coord[2]] # atom2 - atom1
    v2=[atom3_coord[0] - atom1_coord[0],atom3_coord[1] - atom1_coord[1],atom3_coord[2] - atom1_coord[2]] #atom3 - atom1
    angle=calculate_angle_v1_v2(v1,v2)
    return angle

def calculate_dihedral(atom1_coord,atom2_coord,atom3_coord,atom4_coord):
    """"Calculate the dihedral between four atoms.
    
    Parameters:
    atom1_coord - xyz coordinates of atom1 (1N array)
    atom2_coord - xyz coordinates of atom2 (1N array)
    atom3_coord - xyz coordinates of atom3 (1N array)
    atom4_coord - xyz coordinates of atom4 (1N array)
    
    Return:
    dihedral - dihedral between atom1, atom2, atom3, and atom4 (float)
    """
    # Calculate coordinates for vectors v21, v32 and v43
    v21=[atom2_coord[0] - atom1_coord[0],atom2_coord[1] - atom1_coord[1],atom2_coord[2] - atom1_coord[2]] # a2 - a1
    v32=[atom3_coord[0] - atom2_coord[0],atom3_coord[1] - atom2_coord[1],atom3_coord[2] - atom2_coord[2]] # a3 - a2
    v43=[atom4_coord[0] - atom3_coord[0],atom4_coord[1] - atom3_coord[1],atom4_coord[2] - atom3_coord[2]] # a4 - a3    
    # Calculate cross vectors
    v21_x_v32=np.cross(v21,v32)
    v32_x_v43=np.cross(v32,v43)   
    # Calculate normal vectros
    n1=v21_x_v32/(np.sqrt(np.dot(v21_x_v32,v21_x_v32)))
    n2=v32_x_v43/(np.sqrt(np.dot(v32_x_v43,v32_x_v43)))    
    # Unit vectors and y and x
    m1 = n2
    m2 = v32/(np.sqrt(np.dot(v32,v32)))
    m3 = np.cross(m2,m1)
    y=np.dot(n1,m3)
    x=np.dot(n1,m1)    
    # Calculate Dihedral
    dihedral = -math.atan2(y,x)
    dihedral = np.degrees(dihedral)
    return dihedral


## Get the arguments.
parser = argparse.ArgumentParser(description="This script analyzes a user given xyz file and write a new z-matrix file")
parser.add_argument("xyz_file", help="The filepath for the xyz file.")
args = parser.parse_args()

## Open file
file = os.path.join(args.xyz_file)
symbols, coords, num_atoms = open_xyzfile(file)
file_name=os.path.basename(file)
split_file=file.split('.')
name=str(split_file[0])
name=name+'.zmat'

## Distances calculations
number1=[] # Second colum in z-matrix
distances=[] # Stores distances values
for i in range(3, num_atoms):
    for j in range(0, i):
        distance=calculate_distance(coords[i],coords[j])
        parameter=calculate_parameter(symbols[i],symbols[j])
# Check if a distance is a bond based on sum of vdw radii
        if bond_check(distance,parameter) is True:
            number1.append(j+1)
            distances.append(distance)
            break
# Checks if the atom is not connected to any other according to the vdw criteria, if true, its distance from the previous atom is calculated
    if len(number1) != i-2:
        distance=calculate_distance(coords[i],coords[i-1])
        number1.append(i)
        distances.append(distance)

## Angles calculations
number2=[] # Fourth colum in z-matrix
angles=[] # Stores angles values
for i in range(3, num_atoms):
    for j in range(0, i):
        if number1[i-3] - 1 != j: 
            distance=calculate_distance(coords[number1[i-3]-1],coords[j])
            parameter=calculate_parameter(symbols[number1[i-3]-1],symbols[j])
# Check if a distance is a bond based on sum of vdw radii
            if bond_check(distance,parameter) is True:      
                angle=calculate_angle(coords[number1[i-3]-1],coords[j],coords[i])
                number2.append(j+1)
                angles.append(angle)
                break                
# Checks if the atom is not connected to any other according to the vdw criteria, if true, its angle with the previous atom is calculated
    if len(number2) != i-2:
        for j in range(i-1,-1,-1):
            if number1[i-3] - 1 != j:
                number2.append(j+1)
                angle=calculate_angle(coords[number1[i-3]-1],coords[j],coords[i])
                angles.append(angle)
                break

## Dihedrals calculations
number3=[] # Sixth colum in z-matrix
dihedrals=[] # Stores adihedrals values
for i in range(3, num_atoms):
    for j in range(0, i):
        if number1[i-3] - 1 != j and number2[i-3] - 1 != j : 
            distance=calculate_distance(coords[number2[i-3]-1],coords[j])
            parameter=calculate_parameter(symbols[number2[i-3]-1],symbols[j])
# Check if a distance is a bond based on sum of vdw radii            
            if bond_check(distance,parameter) is True:      
                dihedral=calculate_dihedral(coords[i],coords[number1[i-3]-1],coords[number2[i-3]-1],coords[j])
                number3.append(j+1)
                dihedrals.append(dihedral)
                break
# Checks if the atom is not connected to any other according to the vdw criteria, if true, its dihedral with the previous atom is calculated
    if len(number3) != i-2:
        for j in range(i-1,-1,-1):
            if number1[i-3] - 1 != j and number2[i-3] - 1 != j :
                number3.append(j+1)
                dihedral=calculate_dihedral(coords[i],coords[number1[i-3]-1],coords[number2[i-3]-1],coords[j])
                dihedrals.append(dihedral)
                break

## Writing results 
zmat_file=open(name,'w+')
# Writing first 2 lines
zmat_file.write(F'{symbols[0]}\n')
d1=calculate_distance(coords[0],coords[1])
zmat_file.write(F'{symbols[1]}\t1\t{d1:.3f}\n')
# Writing the third line
d2=calculate_distance(coords[0],coords[2])
p1=calculate_parameter(symbols[0],symbols[2])
if bond_check(d2,p1) is True:    
    a1=calculate_angle(coords[0],coords[1],coords[2])
    zmat_file.write(F'{symbols[2]}\t1\t{d2:.3f}\t2\t{a1:.2f}\n')
else:
    a2=calculate_angle(coords[1],coords[0],coords[2])
    d3=calculate_distance(coords[1],coords[2])
    zmat_file.write(F'{symbols[2]}\t2\t{d3:.3f}\t1\t{a2:.2f}\n')
# Writing all other lines
for i in range(0,num_atoms-3):
    zmat_file.write(F'{symbols[i+3]}\t{number1[i]}\t{float(distances[i]):.3f}\t{number2[i]}\t{float(angles[i]):.2f}\t{number3[i]}\t{float(dihedrals[i]):.2f}\n') 
zmat_file.close()
