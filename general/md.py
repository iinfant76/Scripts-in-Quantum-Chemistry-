__all__ = ['make_bond_matrix', 'get_numberofatoms', 'read_atomlist',
           'read_xyz_coordinates', 'rdf']

import numpy as np 
import subprocess
import pandas as pd 
from .md_f import make_bond_matrix_f, make_angle_matrix_f

def make_bond_matrix(n_atoms, coords):
    #Build a tensor made (n_atoms, axes, n_atoms), where axes = x-x0, y-y0, z-z0
    dist = np.stack(
            np.stack(
                    (coords[i_atm, i_ax] - coords[:, i_ax]) ** 2 for i_ax in range(coords.shape[1])
                    )
            for i_atm in range(n_atoms)
            )
    # Builds the bond distance matrix between atoms 
    r = np.stack(
            np.sqrt(np.sum(dist[i_atm, :, :], axis=0)) for i_atm in range(n_atoms)
            )
    return r

def get_numberofatoms(fn): 
    # Retrieve number of lines in trajectory file
    cmd = "wc -l {}".format(fn)
    l = subprocess.check_output(cmd.split()).decode()
    n_lines = int(l.split()[0]) # Number of lines in traj file

    # Read number of atoms in the molecule. It is usually the first line in a xyz file  
    with open(fn) as f:
       l = f.readline()
       n_atoms = int(l.split()[0])

    # Get the number of frames in the trajectory file
    n_frames = int(int(n_lines)/(n_atoms+2))
    
    return n_atoms, n_frames 

def read_atomlist(fn, n_atoms):
    # Read atomic list from xyz file
    atoms = pd.read_csv(fn, nrows = n_atoms, delim_whitespace=True, header=None, 
                        skiprows = 2, usecols=[0]).astype(str).values
    return atoms     


def read_xyz_coordinates(fn, n_atoms, iframe): 
    # Read xyz coordinate from a (trajectory) xyz file. 
    coords = pd.read_csv(fn, nrows = n_atoms, delim_whitespace=True, header=None, 
             skiprows = (2 + (n_atoms + 2) * (iframe - 1)), usecols=(1,2,3)).astype(float).values
    
    return coords  

def rdf(fn, atoms_i, atoms_j, dr, rmax, start, stride): 
    # Get the number of atoms and number of frames in trajectory file
    n_atoms, n_frames = get_numberofatoms(fn) 

    # Create grid of r values to compute the pair correlation function g_ij 
    r_grid = np.arange(0, rmax, dr)
    
    # Compute the volume in a concentric sphere of size dr at a distance r from the origin  
    volume = 4 * np.pi * np.power(r_grid, 2) * dr # elemental volume dV = 4*pi*r^2*dr
    vol_tot = np.sum(volume) 
    
    # Read the order of atoms from the first frame of the trajectory. This order is unchanged. 
    atoms = read_atomlist(fn, n_atoms) 

    # Find the indexes in the bond_matrix of the atomic types involved in the calculation of g_ij
    index_i = np.where( atoms == atoms_i )
    index_j = np.where( atoms == atoms_j )

    for iframe in range(start, n_frames, stride): 
        # Read coordinates from iframe
        coords = read_xyz_coordinates(fn, n_atoms, iframe) 
        # Compute bond distance matrix for iframe
        bond_mtx = make_bond_matrix_f(coords, n_atoms)
        # Slice the bond_matrix with only the atom types
        sliced_mtx = np.triu(bond_mtx[np.ix_(index_i[0], index_j[0])]) 
        # Count the number of atoms within r and r+dr  
        rho_ab = np.stack( 
                    np.where( (sliced_mtx > r_grid[idr]) & (sliced_mtx < r_grid[idr+1]) )[0].size 
                    for idr in range(r_grid.size-1))
        # Compute the total density of pair of atoms
        rho_tot = np.sum(rho_ab) 
        rho_bulk = rho_tot / vol_tot  
#        density = counts / volume[1:] # skip the first element dV which is 0 
        g_ab = rho_ab / rho_bulk 
        # Store density in g_ij
        if (iframe == start):
            g_ij = g_ab
        else:
            g_ij = np.column_stack((g_ij,g_ab))

    # Average over all frames
    g_ij_av = np.average(g_ij, axis=1)

    return r_grid[1:], g_ij_av 

def adf(fn, atoms_i, atoms_j, atoms_k, start, stride): 
    # Get the number of atoms and number of frames in trajectory file
    n_atoms, n_frames = get_numberofatoms(fn) 

    # Create grid of r values to compute the pair correlation function g_ij 
    ang_grid = np.arange(0, 360, 1)
    
    # Read the order of atoms from the first frame of the trajectory. This order is unchanged. 
    atoms = read_atomlist(fn, n_atoms) 

    # Find the indexes in the bond_matrix of the atomic types involved in the calculation of g_ij
    index_i = np.where( atoms == atoms_i )
    index_j = np.where( atoms == atoms_j )
    index_k = np.where( atoms == atoms_k )
    
    for iframe in range(start, n_frames, stride): 
        # Read coordinates from iframe
        coords = read_xyz_coordinates(fn, n_atoms, iframe) 
        # Compute bond distance matrix for iframe
        bond_mtx = make_bond_matrix_f(coords, n_atoms)
        # Compute the angle matrix and convert it in degrees 
        teta_mtx = np.degrees(make_angle_matrix_f(bond_mtx, coords, n_atoms))
        # Slice the bond_matrix with only the atom types
        sliced_mtx = np.triu(teta_mtx[np.ix_(index_i[0], index_j[0], index_k[0])]) 
        # Count the number of atoms within r and r+dr  
        rho_ab = np.stack( 
                    np.where( (sliced_mtx > ang_grid[idr]) & (sliced_mtx < ang_grid[idr+1]) )[0].size 
                    for idr in range(ang_grid.size-1))
        # Compute the total density of pair of atoms
        g_abc = rho_ab 
        # Store density in g_ijk
        if (iframe == start):
            g_ijk = g_abc
        else:
            g_ijk = np.column_stack((g_ijk,g_abc))

    # Average over all frames
    g_ijk_av = np.average(g_ijk, axis=1)

    return ang_grid[1:], g_ijk_av 



