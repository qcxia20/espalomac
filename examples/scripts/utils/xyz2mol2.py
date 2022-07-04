#!/pubhome/qcxia02/miniconda3/envs/espaloma/bin/python3.8
import os
import sys
from pathlib import Path
from ase.units import Bohr

# import espaloma as esp
# import torch
import dgl
import json
from openff.toolkit.topology import Molecule

UNICON_EXE = "/pubhome/soft/ZBH2021/unicon/unicon"

def xyz2mol2(onedatapath,outputpath_prefix,tmppath):
    path = Path(onedatapath)

    clustername = path.name
    path = str(path)
    clusnum = int(path.split("-")[-1])
    
    # homograph = dgl.load_graphs(path + "/homograph.bin")[0][0]
    heterograph = dgl.load_graphs(path + "/heterograph.bin")[0][0]
            
    with open(path + "/mol.json", "r") as f_handle:
        mol = json.load(f_handle)
            
    try:
        mol = Molecule.from_json(mol)
    except:
        mol = Molecule.from_dict(mol)
            
    rdmol = mol.to_rdkit()
    xyzs_groups = heterograph.nodes["n1"].data["xyz"].transpose(0,1)
            
    from rdkit import Chem
    from rdkit.Geometry import Point3D
    from rdkit.Chem import AllChem
            
    AllChem.EmbedMolecule(rdmol,useRandomCoords=True)
            
    for xyzs in xyzs_groups:
        conf = rdmol.GetConformer()
            
        for i in range(rdmol.GetNumAtoms()):
            x,y,z = xyzs[i]
            conf.SetAtomPosition(i,Point3D(float(x)*Bohr,float(y)*Bohr,float(z)*Bohr))
        rdmol.AddConformer(conf, assignId=True)
            
    # rdmol.SetProp("_Name", clustername)
    rdmol.SetProp("_Name", f"tmp{clusnum}")

            
    writer = Chem.SDWriter(f"{tmppath}/{clustername}.sdf")
    for i,conf in enumerate(rdmol.GetConformers()):
        if i>0:
            writer.write(rdmol, confId = conf.GetId())
    writer.close()
            
    os.system(f'{UNICON_EXE} -i "{tmppath}/{clustername}.sdf" -o "{outputpath_prefix/clustername}.mol2"')
    os.system(f'rm "{tmppath}/{clustername}.sdf"')
    return 

if __name__ == "__main__":
    prefix = sys.argv[1]

    dataset_name = "gen2"
    dataset_path = Path("/pubhome/qcxia02/git-repo/espaloma/example") / dataset_name
    outputpath = Path("/pubhome/qcxia02/work/espalomac/examples/datasets/mol2s")

    outputpath_prefix = outputpath / prefix
    if not outputpath_prefix.exists():
        outputpath_prefix.mkdir()

    for path in dataset_path.glob(f"{prefix}*"):
        xyz2mol2(path,outputpath_prefix)

    
        

