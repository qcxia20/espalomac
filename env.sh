#!/bin/bash
# To create a conda environment that contains all the following dependencies. Here is to install cpu-only
# Dependencies:
# openmm openmmtools openmmforcefields rdkit openff-toolkit dgl torch qcportal
conda create -n espaloma
conda activate espaloma
conda install -c conda-forge python==3.8.4
conda install -y rdkit -c conda-forge # rdkit2022
conda install -y pytorch==1.10.1 torchvision==0.11.2 torchaudio==0.10.1 cpuonly -c pytorch
conda install -y -c dglteam dgl
conda install -c conda-forge openmmforcefields #0.11 to match 3.8 version, openmm and openff-toolkit will be installed at the same time
conda install -y -c conda-forge openmmtools
pip install qcportal #qcportal is used to submit QCarchive jobs and parse

git clone https://github.com/choderalab/espaloma.git
cd espaloma
python setup.py install
