**************************
Tutorial for MD in CP2k
**************************

This tutorial will help you in setting up input files for CP2k to run classical molecular dynamics simulations. 
There are some tedious steps that need to be done in order to run such simulations, especially if you do not have a .pdb file available and need to make your own. This tutorial will make several assumptions (e.g. you already have a force-field available for your system), however it will certainly help you in preparing your input files quickly.
Because my field of research is in colloidal nanocrystals, then the set-up proposed here is mostly for those type of systems. However, the scripts are generally applicable to other type of systems as well. 

Pre-requirements 
================
* Force-field parameters

It is required that you already have force-field (FF) parameters for your system already available. 

In the case of a nanocrystal passivated with organic ligands immersed in a solvent, usually you can obtain the parameters for the organic ligands from the `CGenFF <https://cgenff.paramchem.org/>` online library . The capping ligands are indeed made by organic molecules, e.g. oleate, alkyl-ammonium, etc., which have been widely employed in biosystems, and it will be rather straightforward to obtain good starting FF parameters. The same applies to solvents, for which force-fields have been widely tested. 

In case of semiconductor II-VI, IV-VI, perovksite nanocrystals, the availability of FF parameters is less certain and it is very likely that you need to build your own. For example, in case of CdSe nanocrystals we have made a parameterization ourselves that you can find in this publication: `https://pubs.acs.org/doi/pdfplus/10.1021/acs.jctc.6b01089'. For other systems, you need to rely to what is already available in the literature and fine-tune to your goals. 

* XYZ coordinates for the whole system

We assume that you have built your system (NC+ligands+solvents) in such a way that each fragment is stacked one over another. This means that your xyz file should look like this: 

CsPbBr3 with methylammonium passivated at the surface in toluene::

    # Coordinates of the nanocrystal core
    # All atoms of one type are grouped together
    # Cs atoms
    Cs x.xxxx y.yyyy z.zzzz
    Cs x.xxxx y.yyyy z.zzzz
    ...
    Cs x.xxxx y.yyyy z.zzzz
    # Pb atoms 
    Pb x.xxxx y.yyyy z.zzzz
    Pb x.xxxx y.yyyy z.zzzz
    ...
    Pb x.xxxx y.yyyy z.zzzz
    # Br atoms 
    Br x.xxxx y.yyyy z.zzzz
    Br x.xxxx y.yyyy z.zzzz
    ...
    Br x.xxxx y.yyyy z.zzzz
    # Then coordinates of the ligands 
    # Ligand 1: methylammonium 
    N  x.xxxx y.yyyy z.zzzz
    C  x.xxxx y.yyyy z.zzzz
    H  x.xxxx y.yyyy z.zzzz
    H  x.xxxx y.yyyy z.zzzz
    H  x.xxxx y.yyyy z.zzzz
    # Ligand 2: methylammonium (same atomic order as Ligand 1) 
    N  x.xxxx y.yyyy z.zzzz
    C  x.xxxx y.yyyy z.zzzz
    H  x.xxxx y.yyyy z.zzzz
    H  x.xxxx y.yyyy z.zzzz
    H  x.xxxx y.yyyy z.zzzz   
    ...
    # Ligand N: methylammonium (same atomic order as Ligand N)
    N  x.xxxx y.yyyy z.zzzz
    C  x.xxxx y.yyyy z.zzzz
    H  x.xxxx y.yyyy z.zzzz
    H  x.xxxx y.yyyy z.zzzz
    H  x.xxxx y.yyyy z.zzzz   
    # Then coordinates of the solvent 
    # Solvent 1: toluene     
    ...
    # Solvent 2: toluene (same atomic order as Solvent 1)
    ...
    # Solvent M: toluene (same atomic order as Solvent 2)
    ...



    




