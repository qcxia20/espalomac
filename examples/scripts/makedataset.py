#!/pubhome/qcxia02/miniconda3/envs/espaloma/bin/python3.8
import os
import sys
from pathlib import Path, PosixPath
import glob
import json
import shutil
import argparse

from rdkit import Chem
import torch
import dgl
from ase.units import Bohr

import espaloma as esp
from espaloma.graphs.graph import Graph
from openff.toolkit.topology import Molecule

import utils

UNICON_EXE = "/pubhome/soft/ZBH2021/unicon/unicon"
CGENFF_EXE = "/pubhome/qcxia02/soft/cgenff/cgenff"
CHARMM = "/pubhome/qcxia02/soft/charmm/bin/charmm"
utilspath = "/pubhome/qcxia02/work/espalomac/examples/scripts/utils"

def reorgxyz(xyzs, atommaplist):
    newxyzs = []
    for i in atommaplist:
        newxyzs.append(xyzs[i])
    return newxyzs


if __name__ == "__main__":
    parser = argparse.ArgumentParser("make a espaloma dealable dataset (from a dataset of sdf molecules)")
    parser.add_argument("--prefix", type=str, help="prefix for all confs of a single molecule")
    parser.add_argument("--count", type=int, help="count as the folder name to create different working dirs for different molecules")
    parser.add_argument("--rawdataset", type=PosixPath, help="path of raw dataset, which contains molecules and related infos")
    args = parser.parse_args()

    prefix = args.prefix
    globprefix = glob.escape(prefix) # for Path.glob, [ will be transformed into [[]
    count = args.count
    # dataset_name = "TorsionNet500_mm_opt_geometries"
    # dataset_path = Path("/pubhome/qcxia02/dataset/TorsionNet500") / dataset_name
    dataset_path = args.rawdataset
    dataset_name = dataset_path.name

    mol2outputpath = Path("/pubhome/qcxia02/work/espalomac/examples/datasets/mol2s")
    crdoutputpath = Path("/pubhome/qcxia02/work/espalomac/examples/datasets/crds")
    dataoutputpath = Path("/pubhome/qcxia02/work/espalomac/examples/datasets") / dataset_name
    logpath = Path("/pubhome/qcxia02/work/espalomac/examples/datasets/logs")
    tmppath = Path(f"/tmp/qcxia02/{count}")
    
    if not dataoutputpath.exists():
        dataoutputpath.mkdir()

    mol2outputpath_prefix = mol2outputpath / prefix
    crdoutputpath_prefix = crdoutputpath / prefix
    logpath_prefix = logpath / prefix
    dataoutputpath_prefix = dataoutputpath / prefix

    for path in [mol2outputpath_prefix, crdoutputpath_prefix, dataoutputpath_prefix, logpath_prefix]:
        if path.exists():
            shutil.rmtree(path)
        path.mkdir()

    if not Path(f"/tmp/qcxia02").exists():
        Path(f"/tmp/qcxia02").mkdir()
    if tmppath.exists():
        shutil.rmtree(tmppath)
    tmppath.mkdir()

    shutil.copytree(f"{utilspath}/toppar", f"{tmppath}/toppar")
    shutil.copy(f"{utilspath}/energy.inp", tmppath)
    shutil.copy(f"{utilspath}/toppar.str", tmppath)
    shutil.copy(f"{utilspath}/toppar_omm.str", tmppath)

    allxyzs = []
    relmmElist, qmElist, relqmElist = [], [], [] # MM and QM energy in dataset
    totElist, totE_nodihe_list, Edihe_list = [], [], [] # MM energy calculated by CHARMM
    tmplist = []
    # mapdiffcount = 0

    for i, sdffile in enumerate(dataset_path.glob(f"{globprefix}*")):
        ## NOTE: 1. sdf2xyz, sdf2mol2, get energy from file
        rdmol = Chem.SDMolSupplier(str(sdffile), removeHs = False)[0]
        smi = Chem.MolToSmiles(rdmol)
        relmmElist.append(float(rdmol.GetProp("Rel_energy_MM")))
        qmElist.append(float(rdmol.GetProp("PSI4_ENERGY")))
        relqmElist.append(float(rdmol.GetProp("rel_energy_QM")))

        atommaplist = list(map(int, rdmol.GetProp("_smilesAtomOutputOrder")[1:-2].split(",")))
        if i == 0:
            tmplist = atommaplist
        elif atommaplist != tmplist:
            sys.exit(f">>> Atom order between different conformations of {prefix} are not the same. Please check")

        xyzs = utils.sdf2xyz(sdffile)
        newxyzs = reorgxyz(xyzs, atommaplist)
        allxyzs.append(newxyzs)
        mol2outfile = mol2outputpath_prefix/f"tmp{i}.mol2"
        mol2filename = mol2outfile.name
        mol2filenameprefix = mol2filename.split(".mol2")[0]
        os.system(f'{UNICON_EXE} -i "{sdffile}" -o "{mol2outfile}"')
        os.system(f'sed -i "2c {mol2filenameprefix}" {mol2outfile}')# change name in case of unexpected resname 

        if i == 0:
            ## NOTE: 2.0 get graph
            mol = Molecule.from_smiles(smi,allow_undefined_stereo=True,hydrogens_are_explicit=True)
            graph = Graph(mol)
            heterograph = graph.heterograph
            homograph = graph.homograph

            ## NOTE: 2.1 sdf2mol2, assign CGenFF parameters (get .str file)
            outstrfile = mol2outputpath_prefix / (mol2filenameprefix + ".str")
            outstrallfile = mol2outputpath_prefix / (mol2filenameprefix + "_all.str")
            os.system(f'{CGENFF_EXE} "{mol2outfile}" -i mol2 -f "{outstrfile}"')
            os.system(f'{CGENFF_EXE} "{mol2outfile}" -i mol2 -a -f "{outstrallfile}"')

        ## NOTE: 2.2 mol2tocrd
        mol2lines = mol2outfile.read_text().split("\n")
        crdlines, atomsum = utils.mol2tocrdlines(mol2lines)
        outname = mol2filenameprefix + ".crd"
        outcrdfile = crdoutputpath_prefix / outname
        utils.writecrdlines(outcrdfile, crdlines, atomsum)


        ## NOTE: 3. calculate CHARMM energy
        os.chdir(tmppath) 
        cwd = os.getcwd()
        shutil.copy(outstrfile, tmppath/"tmp.str")
        tmpname = outstrfile.name.split(".str")[0]
        shutil.copy(outcrdfile, (tmppath/"tmp.crd"))
        outfile = logpath_prefix / (tmpname+'.out')
        os.system(f'{CHARMM} -i energy.inp workdir={cwd} resi=tmp0 str={cwd}/tmp.str crd={cwd}/tmp.crd cgenff:true modcutoff:true -o energy.out')
        shutil.copy("energy.out", outfile)
        # os.remove(tmppath/"tmp.crd")
        # os.remove(tmppath/"tmp.str")
    
    
        ## NOTE: 4. parse charmm energy
        totE, totE_nodihe = utils.Cout2E(outfile)
        Edihe = totE - totE_nodihe
                
        totElist.append(totE)
        totE_nodihe_list.append(totE_nodihe)
        Edihe_list.append(Edihe)


    # 5. sum and write
    print(torch.tensor(allxyzs).transpose(0,1).shape)
    print(torch.tensor(allxyzs).transpose(0,1))

    print(torch.tensor(allxyzs).shape)

    heterograph.nodes["n1"].data["xyz"] = torch.tensor(allxyzs).transpose(0,1)
    heterograph.nodes["g"].data["u_mm"] = torch.tensor([relmmElist])
    heterograph.nodes["g"].data["u_qm"] = torch.tensor([qmElist])
    heterograph.nodes["g"].data["u_ref"] = torch.tensor([relqmElist])
    heterograph.nodes["g"].data["u_charmm_totE"] = torch.tensor([totElist])
    heterograph.nodes["g"].data["u_charmm_totE_nodihe"] = torch.tensor([totE_nodihe_list])
    heterograph.nodes["g"].data["u_charmm_Edihe"] = torch.tensor([Edihe_list])
    

    dgl.save_graphs(str(dataoutputpath_prefix/"heterograph.bin"), [heterograph])
    dgl.save_graphs(str(dataoutputpath_prefix/"homograph.bin"), [homograph])
    with open(str(dataoutputpath_prefix/"mol.json"), 'w') as f:
        json.dump(mol.to_json(), f)

    shutil.rmtree(mol2outputpath_prefix)
    shutil.rmtree(crdoutputpath_prefix)
    shutil.rmtree(logpath_prefix)
    shutil.rmtree(tmppath)
