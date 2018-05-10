************
Installation
************

qm_scripts does not require complicated installations and a recent version of Python 3, combined with numpy and matplotlib libraries is sufficient.
It is adviced that you create a conda environment and install qm_scripts directly from there using a python package. 
The adviced installation procedure is as follows: 

Conda Installation 
==================
- Download miniconda for python3: miniconda (also you can install the complete anaconda version).
- Install according to: installConda.
- Reopen terminal (or type source ``~/.bashrc``).
- Create a new virtual environment using the following commands

.. code:: bash

     conda create -n qm_scripts

- Activate the new virtual environment

.. code:: bash 

     source activate qm_scripts

To exit the virtual environment type ``source deactivate``.

Via Python Package
==================

Install the package (or add it to your ``requirements.txt`` file). In this case, the dependencies will be installed automatically. 

.. code:: bash

    pip install qm_scripts  

Via Git or Download
===================

You can also install this package by ``git clone`` from a git repository:

.. code:: bash

    git clone https://github.com/iinfant76/qm_scripts.git

In this case you have to ``pip install`` numpy and matplotlib separately or you need to have them installed already.  

If you are interested in the latest developments of this repository, then you can access the ``develop`` branch using:

.. code:: bash 

    git checkout develop 
    git pull 



