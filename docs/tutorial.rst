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

We assume that you have built your system (NC+ligands+solvents) in such a way that each fragment is stacked one on top of each other. This means that your xyz file should look like this: 

CsPbBr\ :sub:`3`\ with methylammonium passivated at the surface in toluene (CsPbBr\ :sub:`3`\.xyz)::

    XXXX # Total number of atoms for the whole system 
    
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
    H  x.xxxx y.yyyy z.zzzz
    H  x.xxxx y.yyyy z.zzzz
    C  x.xxxx y.yyyy z.zzzz
    H  x.xxxx y.yyyy z.zzzz
    H  x.xxxx y.yyyy z.zzzz
    H  x.xxxx y.yyyy z.zzzz
    # Ligand 2: methylammonium (same atomic order as Ligand 1) 
    N  x.xxxx y.yyyy z.zzzz
    H  x.xxxx y.yyyy z.zzzz
    H  x.xxxx y.yyyy z.zzzz
    C  x.xxxx y.yyyy z.zzzz
    H  x.xxxx y.yyyy z.zzzz
    H  x.xxxx y.yyyy z.zzzz
    H  x.xxxx y.yyyy z.zzzz   
    ...
    # Ligand N: methylammonium (same atomic order as Ligand N)
    N  x.xxxx y.yyyy z.zzzz
    H  x.xxxx y.yyyy z.zzzz
    H  x.xxxx y.yyyy z.zzzz
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

It is paramount that your system follows this order, otherwise the scripts needed to write psf and pdb files will not work. 

Building the PSF file for each molecular/ionic fragment  
======================
A PSF file of the whole system (NC+ligands+solvents) is mandatory in your MD simulations. This file contains the information of the atomic connectivity of your system. In order to build such file, we will make use of the script ``xyz2psf.py`` that is available in this repository. The best way to tackle this is to construct a .psf file for each of the fragment types present in your system.

* Building psf files for the nanocrystal core

To make it more practical, here, we will focus on the CsPbBr\ :sub:`3`\ system. The core of the nanocrystal is usually made by 2, 3 or more types of atoms. Usually, in a classical force-field, these atoms are treated only with non-bonded interactions (either Lennard-Jones or Buckingham type). For this reason, we need to generate a simple psf file, with no connection between atoms, for each atomic species. This will be done by the ``xyz2psf.py`` script. This script will read an xyz file for the fragment of interest. The xyz file, however, must be augumented with two extra information that are required in the molecular mechanics calculations: the atomic label of the atoms of interest as used in the force-field and the charge of each atom. For example, in the case of Cs ions, you need the following file (which you can call Cs.xyz) ::

    XXX # Number of atoms 
    
    Cs x.xxxx y.yyyy z.zzzz  1.00  Cs 
    Cs x.xxxx y.yyyy z.zzzz  1.00  Cs
    ...
    Cs x.xxxx y.yyyy z.zzzz  1.00  Cs

In the first column you have the atomic label of Cs as in the force-field (in this case it coincides with the actual atomic name), then in columns 2-4 you have the coordinates (in angstrom) of each atomic type. In column 5, you have the charges of the Cs ion used in the MD calculation and in the last column you have the atomic element name. Once you have generated this file, then you can use the ``xyz2psf.py`` script in the following way::

    xyz2psf.py -file Cs.xyz -id 1 -isolated 

-file Cs.xyz will read the xyz file containing the info explained above;

-id is a label to identify the residue number (you can choose anything, but it's always best that you use the same order of fragments used to build the xyz for the whole system), in this case 1;

-isolated means that there is no connection between atoms. 

The next step is to build the .psf files also for the other two atomic species of the nanocrystal core, Pb and Br. First, you have to build the Pb.xyz and Br.xyz as done for Cs. Pb.xyz will look like this::

    XXX # Number of Pb ions 
    
    Pb x.xxxx y.yyyy z.zzzz  2.00  Pb 
    Pb x.xxxx y.yyyy z.zzzz  2.00  Pb
    ...
    Pn x.xxxx y.yyyy z.zzzz  2.00  Pb

and Br.xyz ::

    XXX # Number of Br ions 
    
    Br x.xxxx y.yyyy z.zzzz  -1.00  Br 
    Br x.xxxx y.yyyy z.zzzz  -1.00  Br
    ...
    Br x.xxxx y.yyyy z.zzzz  -1.00  Br

Then, you can execute the ``xyz2psf.py`` script on both files::

    xyz2psf.py -file Pb.xyz -id 2 -isolated 
    xyz2psf.py -file Br.xyz -id 3 -isolated 

As you can see the only difference from before is the id number that is changed for each atomic species. 

* Building psf files for the ligands

Next step is to build the psf file for the each ligand type. If you only have one ligand type, e.g. methylammonium (MA), but several of them passivating the nanocrystal, it is enough to build one .psf file. The coordinates are only important to determine the connectivity, but they play no other role. Therefore, what you have to do is to take the structure of one of the MA molecules/ions of your system and build a new MA.xyz file as you have done for the nanocrystal core. It should look like this::

    7  # Number of atoms of one MA unit 
   
    NG3P  x.xxxx y.yyyy z.zzzz  -0.130 N
    HGA5  x.xxxx y.yyyy z.zzzz   0.250 H
    HGA5  x.xxxx y.yyyy z.zzzz   0.250 H
    CG21  x.xxxx y.yyyy z.zzzz   0.180 C
    HGA2  x.xxxx y.yyyy z.zzzz   0.150 H
    HGA2  x.xxxx y.yyyy z.zzzz   0.150 H
    HGA2  x.xxxx y.yyyy z.zzzz   0.150 H

This time in the first column you see clearly that the atomic name corresponds to the atomic label used in the force-field. In the fift column you find the atomic charge for each elemet: the total charge of this system is +1 as it should be for an ammonium species. Finally, in the last column you have the real atomic name. Note that the label names and charge values given in this example are not realistic and just used for the purpose of this tutorial. If you want to use correct MA parameters, find them in available databases. 

Once you have prepared your MA.xyz, you can then use ``xyz2psf.py`` script in the following way::

    xyz2psf.py -file MA.xyz -id 4 -bond_thresh 1.6 

The main differences with the previous utilization of this script are: (1) the lack of the -isolated flag, as in this case the atoms are connected; and (2) the presence of the -bond_thresh flag, which is a threshold, in Angstrom, that is used to check if atoms are connected. In this case this value corresponds to the default value and you can obtain the same result by executing the script like this ::

    xyz2psf.py -file MA.xyz -id 4

Usually the default value is fine for atoms like C, N, O, H, etc. However, for atoms like S and P you may need to use larger thresholds. You must be aware that with large thresholds, it is possible to find bonds between atoms that are not actually connected. Unfortunately, at the moment, the script is not intelligent enough to discern these situations. Please be aware of using this with caution. 

If you have anothe ligand type together with MA, you will need to create anore XX.xyz file in the same way as done for MA. 

* Building psf files for the solvent molecules

In this case you need to build the psf file for only one solvent type, e.g. toluene or water. The procedure is the same as the ligand type. 

Building the PDB file for the whole system
==========================================
With the same xyz files used to build the psf files, i.e. those that contain the atomic labels for the FF and the atomic charges, you can now construct the pdb file for the whole system. In this case what you need to do is to provide also the xyz file for the whole system, the one that you have prepared at the beginning (in this case you do not need to add charges and atomic labels). Then, you need to executed the ``xyz2pdb.py`` script in the following way::

    xyz2pdb.py -whole CsPbBr3.xyz -nc Cs.xyz Pb.xyz Br.xyz -ligands MA.xyz -n_ligands 24 -solvents toluene.xyz -n_solvents 521

-whole will read the xyz file of the whole system;
-nc is a list of xyz files (with charges and atomic labels) for each atomic species placed in the same order as they appear in the CsPbBr3.xyz file. 
-ligands is a list of xyz files (with charges and atomic labels) for the ligand types. In this case, you need to provide the xyz file for only one ligand type, coordinates are not important here as they are taken from the CsPbBr3.xyz file.
-n_ligands is the number of ligands that appear in your system, which are made of only one ligand type. This number should correspond to the number of ligands stacked in the CsPbBr3.xyz file of the whole system
-solvents is a list of xyz files for the solvent types. In this case, you need to provide the xyz file for only one ligand type, coordinates are not important here as they are taken from the CsPbBr3.xyz file.
-n_solvents is the number of solvents that appear in your system, which are made of only one solvent type. This number should correspond to the number of solvents stacked in the CsPbBr3.xyz file of the whole system

It is mandatory that in executing this script you follow the same order as in the main xyz file that you provide. The outcome is the generation of a pdb file for your entire system 

Running CP2k simulations for your system 
==========================================
Now that you have prepared your pdb and psf files, you can now perform your first simulation of this system using cp2k. In order to do so, you have to modify the &topology section of your cp2k input script in the following way::

    &TOPOLOGY
      CHARGE_BETA
      COORD_FILE_FORMAT PDB
      COORD_FILE_NAME ./CsPbBr3.xyz.pdb
      CONNECTIVITY MOL_SET
      &MOL_SET
        &MOLECULE
          NMOL 1
          CONN_FILE_NAME ./Cs.xyz.psf
        &END
        &MOLECULE
          NMOL 1
          CONN_FILE_NAME ./Pb.xyz.psf
        &END
        &MOLECULE
          NMOL 1
          CONN_FILE_NAME ./Br.xyz.psf
        &END
        &MOLECULE
          NMOL 24
          CONN_FILE_NAME ./MA.xyz.psf
        &END
        &MOLECULE
          NMOL 521
          CONN_FILE_NAME ./toluene.xyz.psf
        &END  
      &END
      &DUMP_PDB
      &END
      &DUMP_PSF
      &END
      &CENTER_COORDINATES
      &END
    &END TOPOLOGY

Note that the keyword charge_beta means that the charges are taken from the pdb file. However, it is preferrable that also the psf files contain the same charges, as it is done here. 
The dump_pdb and dump_psf will write a new psf and pdb script for the entire system and are updated at each steps of the MD simulations. These can then be used for subsequent calculations using CP2k. 











    




