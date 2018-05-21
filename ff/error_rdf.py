import numpy as np 
from general.md import rdf 

fn_qm = '/Users/iinfante/Documents/University/Scripts/test/CsPbBr3_2.4nm_PBE_MD_NVE_2.5fs-pos-1.xyz'
start = 1 
stride = 50
rmax = 15.0
dr = 0.05

pairs = ['Pb Br', 'Cs Br', 'Cs Pb', 'Cs Cs', 'Pb Pb', 'Br Br']

atoms_i = np.stack(pairs[i].split()[0] for i in range(len(pairs)))
atoms_j = np.stack(pairs[j].split()[1] for j in range(len(pairs)))

x, rdf_qm = np.stack(rdf(fn_qm, i, j, dr, rmax, start, stride) for i, j in zip(atoms_i, atoms_j))
  


