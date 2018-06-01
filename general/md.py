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
    n_bins = r_grid.size # Used for the histogram 

    # Read the order of atoms from the first frame of the trajectory. This order is unchanged. 
    atoms = read_atomlist(fn, n_atoms) 

    # Find the indexes in the bond_matrix of the atomic types involved in the calculation of g_ij
    index_i = np.where( atoms == atoms_i )
    index_j = np.where( atoms == atoms_j )

    rdf_tot = np.zeros(n_bins)
    for iframe in range(start, n_frames, stride): 
        # Read coordinates from iframe
        coords = read_xyz_coordinates(fn, n_atoms, iframe) 
        # Compute bond distance matrix for iframe
        bond_mtx = make_bond_matrix_f(coords, n_atoms)
        # Slice the bond_matrix with only the atom types
        sliced_mtx = np.triu(bond_mtx[np.ix_(index_i[0], index_j[0])]) 
        # Count the number of atoms within r and r+dr  
        rdf, r_grid = np.histogram(sliced_mtx, bins = n_bins, range=(0, rmax)) 
        # Sum over all rdf 
        rdf_tot += rdf

    # Center radii 
    r = 0.5 * ( r_grid[1:] + r_grid[:-1] ) 
    # Compute the volume in a concentric sphere of size dr at a distance r from the origin  
    volume = (4/3) * np.pi * ( np.power(r_grid[1:], 3) - np.power(r_grid[:-1], 3) )  # elemental volume dV = 4*pi*r^2*dr
    vol_tot = np.sum(volume)
    # Compute the rho_ab(r) : 
    n_steps = int(n_frames/stride) # Number of frames used to compute rdf 
    rho_ab = rdf_tot / ( volume * n_steps ) # This is the pair density distribution for each frame   
    # Compute the bulk density 
    n_tot = np.sum(rdf_tot[1:]) / n_steps # skip first element otherwise it counts an atom with itself  
    rho_bulk = n_tot / vol_tot  
    # Density = counts / volume[1:] # skip the first element dV which is 0 
    g_ab = rho_ab / rho_bulk 
    # Compute also the potential of mean force 
    w_ab = -np.log(g_ab) 
   
    return r[1:], g_ab[1:], w_ab[1:] # we skip the first element at r=0, which in the histogram is counted many times: counts an atom with itself    

def adf(fn, atoms_i, atoms_j, atoms_k, bond_max_ij, bond_max_jk, start, stride): 
    # Get the number of atoms and number of frames in trajectory file
    n_atoms, n_frames = get_numberofatoms(fn) 

    # Create grid of r values to compute the pair correlation function g_ij 
    ang_grid = np.arange(1, 181, 0.5)
    n_bins = ang_grid.size  
    
    # Read the order of atoms from the first frame of the trajectory. This order is unchanged. 
    atoms = read_atomlist(fn, n_atoms) 

    # Find the indexes in the bond_matrix of the atomic types involved in the calculation of g_ij
    index_i = np.where( atoms == atoms_i )
    index_j = np.where( atoms == atoms_j )
    index_k = np.where( atoms == atoms_k )
    
    adf_tot = np.zeros(n_bins - 1)     
    for iframe in range(start, n_frames, stride): 
        # Read coordinates from iframe
        coords = read_xyz_coordinates(fn, n_atoms, iframe) 
        # Compute bond distance matrix for iframe
        bond_mtx = make_bond_matrix_f(coords, n_atoms)
        # Compute the angle matrix and convert it in degrees 
        teta_mtx = np.degrees(make_angle_matrix_f(bond_mtx, coords, bond_max_ij, bond_max_jk, n_atoms))
        # Slice the bond_matrix with only the atom types
        sliced_mtx = teta_mtx[np.ix_(index_i[0], index_j[0], index_k[0])] 
        i, j, k = np.indices(sliced_mtx.shape)
        condition = (i != j) & (j != k) & (i < k)
        arr = np.extract(condition, sliced_mtx)
        # Count the number of atoms within r and r+dr  
        adf, ang_grid = np.histogram(arr, bins = n_bins - 1, range=(1, 180))
        # Compute the total density of pair of atoms
        adf_tot += adf 

    # Average over all frames
    ang = 0.5 * (ang_grid[1:] + ang_grid[:-1])  
    n_steps = int(n_frames/stride) 
    adf_av = adf_tot / n_steps 

    return ang, adf_av  



