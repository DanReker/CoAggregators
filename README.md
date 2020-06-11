# Code and data to study and predict drug-excipient co-aggregation for nanoparticle formation

For more information, please refer to the associated publications.
https://www.biorxiv.org/content/10.1101/786251v1


If you use this data or code, please kindly cite ` Reker et al. BioRxiv 2019 10.1101/786251v1 `

# Dependencies

The machine learning and the molecular simulations have been enabled through a broad range of available software packages. 
We have generated two distinct conda environments to enable the implementation of these two pipelines.

## Machine Learning
runs in Python 2.7 using the [scikit-learn library](https://scikit-learn.org/) as well as the [RDKit](http://rdkit.org/). 
Optional libraries from
[OpenPyXL](https://openpyxl.readthedocs.io/en/stable/) and [Joblib](https://joblib.readthedocs.io/en/latest/) enable
reading of self-aggregation data from the Shoichet lab as well as parallelization of code for more efficient processing 
(not implemented in this version of the repository). 
A fresh conda environment can be set up to run these experiments using


```
conda create -n python2 python=2.7 
conda activate python2 
conda install scikit-learn
conda install -c rdkit rdkit
conda install -c anaconda openpyxl
conda install -c anaconda joblib
```

## Simulations

Dependent on Python libraries [OpenMM](http://openmm.org/), [openmoltools](https://github.com/choderalab/openmoltools), [Open Babel](http://openbabel.org/), 
[ambermini](https://github.com/choderalab/ambermini) as a stripped down version of [Amber](http://ambermd.org/),
[MDTraj](http://mdtraj.org/), [ParmEd](http://parmed.github.io/ParmEd), and [RDKit](http://rdkit.org/). 
Also requires an installation of [PACKMOL](http://m3g.iqm.unicamp.br/packmol/home.shtml) 
A fresh installation of Ubuntu 18 and anaconda can be set up to run these experiments using

```
conda create --name simulations
conda activate simulations
conda install -c omnia -c conda-forge openmm
python -m simtk.testInstallation
conda install -c openbabel openbabel
conda install -c kyleabeauchamp ambermini
conda install -c conda-forge mdtraj
conda install -c omnia parmed
conda install -c omnia openmoltools
conda install -c rdkit rdkit

sudo apt-get install gfortran
wget https://github.com/m3g/packmol/archive/20.010.tar.gz
tar -zxvf 20.010.tar.gz 
cd packmol-20.010 
./configure
make
export PATH="/home/user/packmol-20.010:$PATH"
```


# Descriptions of folders

###  data 
contains structural information of the used molecules and the screening results

### machine learning 
contains Python code to run predictions of self-aggregation and co-aggregation. All scripts also contain evaluations and bechmarks as well as adversarial control calculations.

### simulations 
contains code to run and automatically analyze short MD simulations to assess interaction potential of drug molecules and excipients
