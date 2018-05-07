from setuptools import setup

setup(
    name='qm_scripts',
    version='0.1.1',
    description='A set of scripts useful for analyzing outputs and setup inputs in quantum chemistry',
    license='GPL',
    url='https://github.com/iinfant76/qm_scripts',
    author=['Ivan Infante'],
    author_email='iinfant76@gmail.com',
    keywords='MD cp2k',
    packages=[
        "general"],
    classifiers=[
        'License :: OSI Approved :: GPL License'
        'Intended Audience :: Science/Research',
        'programming language :: python :: 3.6',
        'development status :: 1- Alpha',
        'intended audience :: science/research',
        'topic :: scientific/engineering :: chemistry'
    ],
    install_requires=[
        'numpy', 'matplotlib', 'scipy', 'qmflows'], 
#    dependency_links=[
#            "https://github.com/iinfant76/Scripts-in-Quantum-Chemistry-/tarball/master#egg=qm_scripts"],
    scripts=[
        'cp2k_md/xyz2pdb.py',
        'cp2k_md/xyz2psf.py']
)

