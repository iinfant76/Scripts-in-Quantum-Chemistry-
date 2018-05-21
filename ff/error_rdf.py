import numpy as np 
from general.md_f import rdf 

# Fixed Parameters for computing the RDF  
start = 1 
stride = 50
rmax = 15.0
dr = 0.05

# List of pairs used to compute rdf. More pairs better match with QM data.  
pairs = ['Pb Br', 'Cs Br', 'Cs Pb', 'Cs Cs', 'Pb Pb', 'Br Br']

# Starting parameters for the MM trajectory. Here we use LJ only, but can be anything. It has to match the list of pairs order. 
eps = 

# Compute the RDF from the reference QM trajectory 
fn_qm = '/Users/iinfante/Documents/University/Scripts/test/CsPbBr3_2.4nm_PBE_MD_NVE_2.5fs-pos-1.xyz'
atoms_i = np.stack(pairs[i].split()[0] for i in range(len(pairs))) # Rearrange the pairs in two lists. It is useful for using the rdf function
atoms_j = np.stack(pairs[j].split()[1] for j in range(len(pairs)))
rdf_qm = np.stack(rdf(fn_qm, i, j, dr, rmax, start, stride)[1] for i, j in zip(atoms_i, atoms_j))

# Compute the MM trajectory 

fn_mm = '/Users/iinfante/Documents/University/Scripts/test/CsPbBr3_2.4nm_MM-pos-1.xyz'
rdf_mm = np.stack(rdf(fn_mm, i, j, dr, rmax, start, stride)[1] for i, j in zip(atoms_i, atoms_j))  

delta_gij = rdf_qm - rdf_mm 

e_qmmm = np.sum(np.sqrt(np.sum(delta_gij, axis = 1)))
 

