#!/bin/python3
#
# dynamics_analysis.py
#
# ************************************************************************************************
# *                          Developed by Daniel Arley Santos Oliveira                           *
# *                             Last version finished in 31/12/2022                              *
# ************************************************************************************************
#
# https://github.com/danielarley/Introduction-to-scientific-computing-and-programming-in-chemistry
#
# This script reads and performs different analyzes on an xyz file containing the cartesian 
# coordinates of a solute molecule and several solvent molecules obtained through molecular dynamics. 
#
# How to use:
# ./dynamics_analysis.py xyz_file.xyz num_solvent_molecules
#
# Parameters:
# xyz_file.xyz - file containing solute and solvent coordinates
# num_solvent_molecules - number of solvent molecules to write in trajectory file
# 

import numpy as np
import os
import argparse
import math
import sys

## Functions ##

def open_xyzfile(file):
    """ Open .xyz file and stores symbols and coordinates in arrays.
    
    Parameters:
    file - xyz file with cartesian coordinates of the system
    
    Return:
    symbols - Element symbols (1N array)
    coordinates - Coordinates of each atom (ND array, where N= number of atoms)
    num_atoms - Number of atoms in the system (int)
    cell_lenght- Box size (float)
    """
    # Skips the first two lines and opens the file like an N-dimensional array
    xyz_file = np.genfromtxt(fname=file,skip_header=2,dtype='unicode') 
    symbols = xyz_file[:,0] # Store the atoms symbols
    coordinates = (xyz_file[:,1:]) # Each array dimension represents the coordinates of an atom 
    coordinates = coordinates.astype(np.float64)
    num_atoms=len(symbols)
    with open(file,'r') as f:
        f.readline()
        cell_lenght=float(f.readline()) # Getting the box size
    return symbols, coordinates, num_atoms, cell_lenght

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
    atom1_symb - Atom1 symbol (str)
    atom2_symb - Atom2 symbol (str)
    
    Return:
    parameter - Parameter that define wheter atom1 and atom2 is bonded (float)
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
    atom_distance - Distance between two atoms (float)
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
    v1- Vector 1 (1N array)
    v2- Vector 2 (1N array)
    
    Return: 
    angle_v1_v2 - Angle between v1 and v2 (float)
    """
    scalar_v1_v2=np.dot(v1,v2) #Scalar product
    mod_v1=np.sqrt(np.dot(v1,v1)) # v1 module
    mod_v2=np.sqrt(np.dot(v2,v2)) # v2 module
    cos_ang=scalar_v1_v2/(mod_v1*mod_v2) # cosine of the angle
    angle_v1_v2=math.acos(cos_ang)*57.2958 # Angle calculate
    return angle_v1_v2

def calculate_angle(atom1_coord,atom2_coord,atom3_coord):
    """" Calculate the angle between three atoms.
    
    Parameters:
    atom1_coord - xyz coordinates of atom1 (1N array)
    atom2_coord - xyz coordinates of atom2 (1N array)
    atom3_coord - xyz coordinates of atom3 (1N array)
    
    Return:
    angle - Angle between atom1, atom2, and atom3 (float)
    """
    v1=[atom2_coord[0] - atom1_coord[0],atom2_coord[1] - atom1_coord[1],atom2_coord[2] - atom1_coord[2]] # atom2 - atom1
    v2=[atom3_coord[0] - atom1_coord[0],atom3_coord[1] - atom1_coord[1],atom3_coord[2] - atom1_coord[2]] #atom3 - atom1
    angle=calculate_angle_v1_v2(v1,v2)
    return angle

def calculate_dihedral(atom1_coord,atom2_coord,atom3_coord,atom4_coord):
    """" Calculate the dihedral between four atoms.
    
    Parameters:
    atom1_coord - xyz coordinates of atom1 (1N array)
    atom2_coord - xyz coordinates of atom2 (1N array)
    atom3_coord - xyz coordinates of atom3 (1N array)
    atom4_coord - xyz coordinates of atom4 (1N array)
    
    Return:
    dihedral - Dihedral between atom1, atom2, atom3, and atom4 (float)
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

def molecule_sep(symbols,coords,num_atoms):
    """ Sorts atoms in molecules.
    
    Parameters:
    symbols - Element symbols (1N array)
    coords - Coordinates of each atom (ND array, where N= number of atoms)
    num_atoms - Number of atoms in the system (int)
    
    Return:
    molecule_number- Indices of the molecules of each atom (list with N elements, where N=number of atoms)
    num_molecules - Number of molecules in the system (int)
    """
    molecule_number=[0]
    for i in range(0,num_atoms):
            for j in range(0, i):
                distance=calculate_distance(coords[i],coords[j])
                parameter=calculate_parameter(symbols[i],symbols[j])
    # Check if a distance is a bond based on sum of vdw radii
                if bond_check(distance,parameter) is True:
    #if atom i is bonded with atom j, then molecule indice j is given to atom i
                    molecule_number.append(molecule_number[j]) 
                    break
    # if atom i and atom j are not bonded, then then molecule indice j + 1 is given to atom i       
            if len(molecule_number) != i+1:   
                molecule_number.append(molecule_number[j]+1)
    #Getting the number of molecules in the system
    num_molecules=max(molecule_number)
    return molecule_number, num_molecules

def get_solute_solvents(n_mol,molecule_number):
    """ Identifies the solute molecule and the closest solvent molecules
    
    Parameters:
    n_mol - number of closest molecules desired (int)
    molecule_number - Indices of the molecules of each atom (list with N elements, where N=number of atoms)
    
    Return:
    solute_atoms - List with the indices of the solute atoms (list)
    distances - List of distances between the closest solvent molecules and solute molecule (list)
    nearby_atoms - List with the indices of the atoms closest to the solute (list)
    nearby_molecules - List with the indices of the solvent molecules closest to the solute (list)
    """
    # Getting solute atoms
    solute_atoms=[]
    count=0
    for i in molecule_number:
        if i == 0:
            solute_atoms.append(count)
        count+=1
    # Stores the shortest distance of each solvent molecule to the solute molecule in a dictionary    
    distances_dic={}
    for i in range(0,num_atoms):
        a=10000
        if i not in solute_atoms:
            for j in solute_atoms:    
                distance=calculate_distance(coords[i],coords[j])
                if distance < a:
                    a=distance
            distances_dic[a]=i
    # Getting the n_mol closest solvent molecules
    a=10000 
    nearby_atoms=[] # Atoms closest to the solute
    nearby_molecules=[] # Solvent molecules closest to the solute
    distances=[] # Distances between the closest solvent molecules and solute molecule
    for i in range(0,n_mol):
        for j in distances_dic.keys():        
            b=molecule_number[distances_dic[j]]
    # Checks if the atom is not already part of one of the nearest molecules
            if j < a and b not in nearby_molecules:
                a=j
        distances.append(a)        
        nearby_atoms.append(distances_dic[a])
        nearby_molecules.append(molecule_number[distances_dic[a]])
        # Deletes the closest atom to the solute from the dictionary (prevents storing the same molecule twice)
        del distances_dic[a] 
        a=10000     
    return solute_atoms, distances, nearby_atoms, nearby_molecules

def solute_and_solvent(name,molecule_number,symbols,coords,nearby_molecules):
    """ Writes an xyz file with the solute and the n_mol closest solvent molecules
    
    Parameters:
    name - File name (str)
    molecule_number - Indices of the molecules of each atom (list with N elements, where N=number of atoms)
    symbols - Element symbols (1N array)
    coords - Coordinates of each atom (ND array, where N= number of atoms)
    nearby_molecules - List with the indices of the solvent molecules closest to the solute (list)
    """
    file_name=name+'_solute_and_closest_solvent_molecules.xyz'
    solute_and_solvent=open(file_name,'w+')
    # Getting solute coordinates
    count=0
    data=[]
    for i,j in enumerate(molecule_number):
        if j == 0:
            a=(F'{symbols[i]}\t{coords[i,0]}\t{coords[i,1]}\t{coords[i,2]}\n')
            data.append(a)
            count+=1
    # Getting the n_mol closest solvent molecules coordinates
    for i in nearby_molecules:
        for j,k in enumerate(molecule_number):
            if i==k:
                b=(F'{symbols[j]}\t{coords[j,0]}\t{coords[j,1]}\t{coords[j,2]}\n')
                data.append(b)
                count+=1
    solute_and_solvent.write(F'{count}\n\n')
    solute_and_solvent.writelines(data)  
    solute_and_solvent.close()

def trajectory(name,molecule_number,symbols,coords,nearby_molecules,distances):
    """ Writes a xyz trajectory file with the solute coordinates and n_mol closest solvent molecules
    
    Parameters:
    name - File name (str)
    molecule_number - Indices of the molecules of each atom (list with N elements, where N=number of atoms)
    symbols - Element symbols (1N array)
    coords - Coordinates of each atom (ND array, where N= number of atoms)
    nearby_molecules - List with the indices of the solvent molecules closest to the solute (list)
    distances - List of distances between the closest solvent molecules and solute molecule (list)
    """
    name_traj=name+'_trajectory.xyz'
    trajectory=open(name_traj,'w+')
    #Writes solute coordinates
    count=0
    data=[]
    for i,j in enumerate(molecule_number):
        if j == 0:
            a=(F'{symbols[i]}\t{coords[i,0]}\t{coords[i,1]}\t{coords[i,2]}\n')
            data.append(a)
            count+=1
    trajectory.write(F'{count}\nAdding Molecule 1 with {count} atom(s) R= 0.00000\n')
    trajectory.writelines(data)
    #Writes the solvent coordinates in trajectory format
    count2=0
    for i in nearby_molecules:
        count3=0
        for j,k in enumerate(molecule_number):
            if i==k:
                b=(F'{symbols[j]}\t{coords[j,0]}\t{coords[j,1]}\t{coords[j,2]}\n')
                data.append(b)
                count+=1
                count3+=1
        trajectory.write(F'{count}\nAdding Molecule {count2+2} with {count3} atom(s) R={distances[count2]:.5f} \n')
        trajectory.writelines(data)
        count2+=1
    trajectory.close()

def putting_in_box(num_atoms,coords,cell_lenght): 
    """ Put atoms outside the box into the box and writes a xyz file with the new coordinates
    
    Parameters:
    num_atoms - Number of atoms in the system (int)
    coords - Coordinates of each atom (ND array, where N= number of atoms)
    cell_lenght - Box size (float)
    
    Return:
    new_coords - Coordinates of each atom with solute in the origin of cartesian coordinates (ND array, where N= number of atoms)
    coords_inside - Coordinates of each atom already inside the box (ND array, where N= number of atoms)
    """
    # Putting the solute in the origin of cartesian coordinates
    new_coords=coords.copy()
    x0=coords[0,0]
    y0=coords[0,1]
    z0=coords[0,2]
    for i in range(0,num_atoms):
        new_coords[i]=[coords[i,0]-x0,coords[i,1]-y0,coords[i,2]-z0] 
    #Putting atoms outside the box into the box
    coords_inside=new_coords.copy()
    for i in range(0,num_atoms):
        for j in range (0,3):
            if coords_inside[i,j] > cell_lenght/2 :
                coords_inside[i,j]-=cell_lenght
            elif coords_inside[i,j] < -cell_lenght/2 :
                coords_inside[i,j]+=cell_lenght
    # Writing a xyz file with all atoms inside the box
    name_inside=name+'_inside.xyz'
    inside=open(name_inside,'w+')
    inside.write(F'{num_atoms}\n\n')
    data=[]
    for i in range(0,num_atoms):
        a=(F'{symbols[i]}\t{coords_inside[i,0]}\t{coords_inside[i,1]}\t{coords_inside[i,2]}\n')
        data.append(a)
    inside.writelines(data)
    inside.close()            
    return new_coords, coords_inside

def create_super_cell(new_coords,cell_lenght,symbols,name):
    """ Creates a super cell three times the size of the first and return coordinates and symbols of the super cell
    
    Parameters:
    new_coords - Coordinates of each atom with solute in the origin of cartesian coordinates (ND array, where N= number of atoms)
    cell_lenght - Box size (float) 
    symbols - Element symbols (1N array)
    name - File name (str)
    
    Return:
    coords_super_cell - Element symbols of each atom of the super cell
    symbols_super_cell - Coordinates of each atom of the super cell
    """
    # Copying the central coordinates to the other positions of the super cube
    coords_super_cell=new_coords.copy()
    for i in range(-1,2):
        for j in range(-1,2):
            for k in range(-1,2):
                position=[i,j,k]
                add_coords=new_coords.copy()
                ini_position=[0,0,0]
                if position != ini_position:
                    for l in range(0,3):
                        if position[l] == 1:
                            for m in range(0,len(add_coords)):
                                add_coords[m,l]+=cell_lenght
                        if position[l] == -1:
                            for m in range(0,len(add_coords)):
                                add_coords[m,l]-=cell_lenght
                    coords_super_cell=np.vstack([coords_super_cell,add_coords])
    # Creating a list with the symbols of the super cube atoms
    symbols_s_c=list(symbols)
    symbols_super_cell=symbols_s_c.copy()
    for i in range(0,26):
        symbols_super_cell=symbols_super_cell+symbols_s_c
    # Writing a xyz file with super cell coordinates
    name_super_cell=name+'_super_cell.xyz'
    super_cell=open(name_super_cell,'w+')
    super_cell.write(F'{str(len(coords_super_cell))}\n\n')
    data=[]
    for i in range(len(coords_super_cell)):
        a=(F'{symbols_super_cell[i]}\t{coords_super_cell[i,0]}\t{coords_super_cell[i,1]}\t{coords_super_cell[i,2]}\n')
        data.append(a)
    super_cell.writelines(data)
    super_cell.close()
    return coords_super_cell, symbols_super_cell

def print_bonds(num_atoms,coords,symbols):
    """ Print bond distances
    
    Parameters:
    symbols - Element symbols (1N array)
    coords - Coordinates of each atom (ND array, where N= number of atoms)
    num_atoms - Number of atoms in the system (int)
    """
    for i in range(0, num_atoms):
        for j in range(0, num_atoms):
            if i < j:
                distance=calculate_distance(coords[i],coords[j])
                parameter=calculate_parameter(symbols[i],symbols[j])
    # Check if a distance is a bond based on sum of vdw radii
                if bond_check(distance,parameter) is True:
                    print(F'R({i+1}-{symbols[i]},{j+1}-{symbols[j]}) = {distance:.4f}')
                
def print_angles_dihedrals(num_atoms,coords,symbols):
    """ Print angles and dihedrals
    
    Parameters:
    symbols - Element symbols (1N array)
    coords - Coordinates of each atom (ND array, where N= number of atoms)
    num_atoms - Number of atoms in the system (int)
    """
    angles_data=[] # list that will stores angles values 
    dihedrals_data=[] # list that will stores dihedrals values 
    for i in range(0, num_atoms):
        for j in range(0, num_atoms):
            if i < j:
                distance=calculate_distance(coords[i],coords[j])
                parameter=calculate_parameter(symbols[i],symbols[j])
    # Check if a distance is a bond based on sum of vdw radii
                if bond_check(distance,parameter) is True:
                    for k in range(0,num_atoms):
                        if k != i and k != j:
                            distance2=calculate_distance(coords[k],coords[j])
                            parameter2=calculate_parameter(symbols[k],symbols[j])
                            distance3=calculate_distance(coords[k],coords[i])
                            parameter3=calculate_parameter(symbols[k],symbols[i])
    # Check if a distance is a bond based on sum of vdw radii
                            if bond_check(distance2,parameter2) is True or bond_check(distance3,parameter3) is True :
                                angle=calculate_angle(coords[i],coords[k],coords[j])
                                A=(F'A({i+1}-{symbols[i]},{k+1}-{symbols[k]},{j+1}-{symbols[j]}) = {angle:.4f}')
                                angles_data.append(A)
                                for l in range(0,num_atoms):
                                    if l != i and l != j and l !=k:
                                        distance4=calculate_distance(coords[l],coords[i])
                                        parameter4=calculate_parameter(symbols[l],symbols[j])
                                        distance5=calculate_distance(coords[l],coords[i])
                                        parameter5=calculate_parameter(symbols[l],symbols[j])
                                        distance6=calculate_distance(coords[l],coords[k])
                                        parameter6=calculate_parameter(symbols[l],symbols[k])
    # Check if a distance is a bond based on sum of vdw radii
                                        if bond_check(distance4,parameter4) is True or bond_check(distance5,parameter5) is True or bond_check(distance6,parameter6) is True :
                                                dihedral=calculate_dihedral(coords[i],coords[j],coords[k],coords[l])
                                                D=(F'D({i+1}-{symbols[i]},{j+1}-{symbols[j]},{k+1}-{symbols[k]},{l+1}-{symbols[l]}) = {dihedral:.4f}')
                                                dihedrals_data.append(D)
    # Printing angles
    print('Angles Data:')
    for i in angles_data:
        print(i)
    # Printing Dihedrals
    print('\nDihedrals Data:')
    for i in dihedrals_data:
        print(i)

## End of functions        

## Open file
file = os.path.join(sys.argv[1])
symbols, coords, num_atoms, cell_lenght = open_xyzfile(file)
split_file=file.split('.')
name=str(split_file[0]) # Get file name

## Printing bonds Distances
print('Printing distances values... \n')
print_bonds(num_atoms,coords,symbols)

## Printing angles and dihedrals
print('Printing angles and dihedrals... \n')
print_angles_dihedrals(num_atoms,coords,symbols)

## Sorting atoms in molecules
print('Sorting atoms in molecules... \n')
molecule_number, num_molecules=molecule_sep(symbols,coords,num_atoms)
print(F'There are {num_molecules} molecules\n')

## Identifying the solute molecule and the closest solvent molecules
print('Identifying the solute molecule and the closest solvent molecules... \n')
solute_atoms, distances, nearby_atoms, nearby_molecules=get_solute_solvents(int(sys.argv[2]),molecule_number)

## Writing an xyz file with the solute and the closest solvent molecules
print('Writing an xyz file with the solute and the closest solvent molecules... \n')
solute_and_solvent(name,molecule_number,symbols,coords,nearby_molecules)

##  Writing a xyz trajectory file with the solute coordinates and the closest solvent molecules 
print('Writing a xyz trajectory file with the solute coordinates and the closest solvent molecules... \n')
trajectory(name,molecule_number,symbols,coords,nearby_molecules,distances)

## Putting atoms outside the box into the box
print('Putting atoms outside the box into the box.. \n')
new_coords, coords_inside= putting_in_box(num_atoms,coords,cell_lenght)

## Printing bonds distances for atoms already inside the box
print('Printing distances values for atoms already inside the box... \n')
print_bonds(num_atoms,coords_inside,symbols)

## Printing angles and dihedrals for the atoms already inside the box
print('Printing angles and dihedrals for the atoms already inside the box... \n')
print_angles_dihedrals(num_atoms,coords_inside,symbols)

## Sorting atoms into molecules for the atoms already inside the box
print('Sorting atoms in molecules for the atoms already inside the box... \n')
molecule_number_inside, num_molecules_inside=molecule_sep(symbols,coords_inside,num_atoms)
print(F'There are {num_molecules_inside} molecules inside the box\n')

## Identifying the solute molecule and the closest solvent molecules for the atoms already inside the box
print('Identifying the solute molecule and the closest solvent molecules for the atoms already inside the box... \n')
solute_atoms_inside, distances_inside, nearby_atoms_inside, nearby_molecules_inside=get_solute_solvents(int(sys.argv[2]),molecule_number_inside)

## Writing an xyz file with the solute and the closest solvent molecules for the atoms already inside the box
print('Writing an xyz file with the solute and the closest solvent molecules for the atoms already inside the box... \n')
name_inside=name+'_inside'
solute_and_solvent(name_inside,molecule_number_inside,symbols,coords_inside,nearby_molecules_inside)

##  Writing a xyz trajectory file with the solute coordinates and the closest solvent molecules for the atoms already inside the box
print('''Writing a xyz trajectory file with the solute coordinates and the closest solvent molecules for the atoms 
already inside the box... \n''')
trajectory(name_inside,molecule_number_inside,symbols,coords_inside,nearby_molecules_inside,distances_inside)

## Writing a xyz file with super cell coordinates
print('Writing a xyz file with super cell coordinates using initial coordinates... \n')
coords_super_cell, symbols_super_cell=create_super_cell(new_coords,cell_lenght,symbols,name)

## Printing bonds distances for atoms of the super cell
print('Printing distances values for atoms of the super cell... \n')
print_bonds(len(coords_super_cell),coords_super_cell,symbols_super_cell)

## Sorting atoms into molecules for the super cell atoms
print('Sorting atoms in molecules for the super cell atoms... \n')
molecule_number_sp, num_molecules_sp=molecule_sep(symbols_super_cell,coords_super_cell,len(coords_super_cell))
print(F'There are {num_molecules_sp} molecules inside the super cell\n')

## Identifying the solute molecule and the closest solvent molecules for the super cell atoms
print('Identifying the solute molecule and the closest solvent molecules for the super cell atoms... \n')
solute_atoms_sp, distances_sp, nearby_atoms_sp, nearby_molecules_sp=get_solute_solvents(int(sys.argv[2]),molecule_number_sp)

## Writing an xyz file with the solute and the closest solvent molecules for the super cell atoms
print('Writing an xyz file with the solute and the closest solvent molecules for the super cell atoms... \n')
name_sp=name+'_super_cel'
solute_and_solvent(name_sp,molecule_number_sp,symbols_super_cell,coords_super_cell,nearby_molecules_sp)

##  Writing a xyz trajectory file with the solute coordinates and the closest solvent molecules for the super cell atoms
print('Writing a xyz trajectory file with the solute coordinates and the closest solvent molecules for the super cell atoms... \n')
trajectory(name_sp,molecule_number_sp,symbols_super_cell,coords_super_cell,nearby_molecules_sp,distances_sp)