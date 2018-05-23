import numpy as np 
from general.md import rdf 
from noodles import gather, lift 
from qmflows import Settings, cp2k, run, templates 

# Fixed Parameters for computing the RDF  
start = 1 
stride = 50
rmax = 15.0
dr = 0.05

# Reference qm trajectory
fn_qm = '/Users/iinfante/Documents/University/Scripts/test/CsPbBr3_2.4nm_PBE_MD_NVE_2.5fs-pos-1.xyz'

# List of pairs used to compute rdf. More pairs better match with QM data.  
pairs = ['Pb Br', 'Cs Br', 'Cs Pb', 'Cs Cs', 'Pb Pb', 'Br Br']

# Starting parameters for the MM trajectory. Here we use LJ only and electrostatics, but can be anything. It has to match the list of pairs order. 
# We will make a dictionary. Best way of using this. 
# Atomic charges
atomic_charge = {
      'Cs' : 1.00, 
      'Pb' : 2.00,
      'Br' : -1.00 } 

# Lennard-Jones parameters in Ke (eps) and angstrom (sigma) 
eps : {
      'Cs Pb' : 4.92033,
      'Cs Br' : 9.481468,
      'Pb Br' : 185.4128,  
      'Pb Pb' : 96.2185, 
      'Cs Cs' : 0.2516114, 
      'Br Br' : 357.29}

sigma : {
      'Cs Pb' : 4.31277,
      'Cs Br' : 9.481468,
      'Pb Br' : 3.58329,
      'Pb Pb' : 3.0, 
      'Cs Cs' : 6.20, 
      'Br Br' : 4.28}

# Compute the RDF from the reference QM trajectory for each pair of atoms  
atoms_i = np.stack(pairs[i].split()[0] for i in range(len(pairs))) # Rearrange the pairs in two lists. It is useful for using the rdf function
atoms_j = np.stack(pairs[j].split()[1] for j in range(len(pairs)))
rdf_qm = np.stack(rdf(fn_qm, i, j, dr, rmax, start, stride)[1] for i, j in zip(atoms_i, atoms_j))

# Compute the MM trajectory 
# ======================
# Settings for CP2k
# ======================

# Set path for basis set 

# Settings specifics  
s = Settings()
#s.basis = "DZVP-MOLOPT-SR-GTH"
#s.potential = "GTH-PBE"
#s.cell_parameters = 42.0

# =======================
job = cp2k(templates.mm.overlay(s), fn_qm) 

energies = gather(*job.energy) 
results = run(energies, folder='.') 

#fn_mm = '/Users/iinfante/Documents/University/Scripts/test/CsPbBr3_2.4nm_MM-pos-1.xyz'
#rdf_mm = np.stack(rdf(fn_mm, i, j, dr, rmax, start, stride)[1] for i, j in zip(atoms_i, atoms_j))  

#delta_gij = rdf_qm - rdf_mm 

#e_qmmm = np.sum(np.sqrt(np.sum(np.power(delta_gij, 2), axis = 1) ) )

 



