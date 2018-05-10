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

In the case of a nanocrystal passivated with organic ligands immersed in a solvent, usually you can obtain the parameters for the organic ligands from the `CGenFF <https://cgenff.paramchem.org/>` online library . The capping ligands are indeed made by organic molecules, e.g. oleate, alkyl-ammonium, etc., which have been widely employed in biosystems, and it will be rather straightforward to obtain good FF parameters.  
The same applies to solvents, for which force-fields have been widely tested. 

In case of semiconductor II-VI, IV-VI, perovksite nanocrystals, the 


