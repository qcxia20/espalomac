#!/pubhome/qcxia02/miniconda3/envs/espaloma/bin/python3.8
import os
import sys
from pathlib import Path
import glob
import json
import shutil

import torch
import dgl
from ase.units import Bohr

import espaloma as esp
from openff.toolkit.topology import Molecule

# from utils import *
import utils

UNICON_EXE = "/pubhome/soft/ZBH2021/unicon/unicon"
CGENFF_EXE = "/pubhome/qcxia02/soft/cgenff/cgenff"
CHARMM = "/pubhome/qcxia02/soft/charmm/bin/charmm"
utilspath = "/pubhome/qcxia02/work/espalomac/examples/scripts/utils"


if __name__ == "__main__":
    prefix = sys.argv[1]
    globprefix = glob.escape(prefix) # for Path.glob, [ will be transformed into [[]
    count = sys.argv[2]
    print(prefix)
    dataset_name = "gen2"
    dataset_path = Path("/pubhome/qcxia02/git-repo/espaloma/example") / dataset_name
    mol2outputpath = Path("/pubhome/qcxia02/work/espalomac/examples/datasets/mol2s")
    crdoutputpath = Path("/pubhome/qcxia02/work/espalomac/examples/datasets/crds")
    dataoutputpath = Path("/pubhome/qcxia02/work/espalomac/examples/datasets/gen2")
    logpath = Path("/pubhome/qcxia02/work/espalomac/examples/datasets/logs")
    tmppath = Path(f"/tmp/qcxia02/{count}")

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

    ## NOTE: 1. xyz2mol2
    for path in dataset_path.glob(f"{globprefix}-*"):
        utils.xyz2mol2(path, mol2outputpath_prefix, tmppath)
    
    ## NOTE: 2.1 derive first mol2 and assign CGenFF parameters (get .str file)
    # derive_first_mol2
    for mol2file in mol2outputpath_prefix.glob("OpenFF*.mol2"):
        mol2filename = mol2file.name
        clusnum = mol2filename.split(".mol2")[0].split("-")[-1]
        outmol2file = mol2outputpath_prefix / ("tmp" + clusnum + ".mol2")
        outstrfile = mol2outputpath_prefix / ("tmp" + clusnum + ".str")
        outstrallfile = mol2outputpath_prefix / ("tmp" + clusnum + "_all.str")
        os.system(f'{utilspath}/derive_first_mol2.py "{mol2file}" "{outmol2file}"')
        os.system(f'{CGENFF_EXE} "{outmol2file}" -i mol2 -f "{outstrfile}"')
        os.system(f'{CGENFF_EXE} "{outmol2file}" -i mol2 -a -f "{outstrallfile}"')

    ## NOTE: 2.2 mol2tocrd
    for mol2file in mol2outputpath_prefix.glob("OpenFF*.mol2"):
        allmol2lines = utils.next_mol2_lines(mol2file)
        allmol2lineslist = list(allmol2lines)
    
        if len(allmol2lineslist) == 1:
            mol2lines = allmol2lineslist[0]
            crdlines, atomsum = utils.mol2tocrdlines(mol2lines)
            
            outname = mol2file.name.split(".mol2")[0] + ".crd"
            utils.writecrdlines(crdoutputpath_prefix / outname, crdlines, atomsum)
        
        else:
            for i, mol2lines in enumerate(allmol2lineslist):
                crdlines, atomsum = utils.mol2tocrdlines(mol2lines)
                
                outname = mol2file.name.split(".mol2")[0] + f"_{i}" + ".crd"
                utils.writecrdlines(crdoutputpath_prefix / outname, crdlines, atomsum)

    ## NOTE: 3. calculate CHARMM energy
    countdict = {}
    os.chdir(tmppath)
    cwd = os.getcwd()
    for strfile in mol2outputpath_prefix.glob("tmp*.str"):
        if "all" not in strfile.name:
            shutil.copy(strfile, tmppath/"tmp.str")
            tmpname = strfile.name.split(".str")[0]
            clusnum = strfile.name.split(".str")[0].split("tmp")[1]
            for count, crdfile in enumerate(crdoutputpath_prefix.glob(f"OpenFF*-{clusnum}_*.crd")):
                j = crdfile.name.split("_")[-1].split(".")[0]
                shutil.copy(crdfile, (tmppath/"tmp.crd"))
                outfile = logpath_prefix / (tmpname+'_'+j+'.out')
                # # path length should not  be too long
                # os.system(f'{CHARMM} -i energy.inp workdir={cwd} resi={tmpname} str={cwd}/tmp.str crd={cwd}/tmp.crd cgenff:true modcutoff:true -o "{outfile}"') 

                os.system(f'{CHARMM} -i energy.inp workdir={cwd} resi={tmpname} str={cwd}/tmp.str crd={cwd}/tmp.crd cgenff:true modcutoff:true -o energy.out')
                shutil.copy("energy.out", outfile)
                
                os.remove(tmppath/"tmp.crd")
            os.remove(tmppath/"tmp.str")
            countdict[clusnum] = count+1


    ## NOTE: 4. parse charmm energy
    for path in dataset_path.glob(f"{globprefix}*"):
        shutil.copytree(path, (dataoutputpath_prefix / path.name))
    
    for clusnum in countdict.keys():
        totElist, totE_nodihe_list, Edihe_list =  [], [], []
        for clusidx in range(countdict[clusnum]):
            outfile = logpath_prefix / f"tmp{clusnum}_{clusidx}.out"
            totE, totE_nodihe = utils.Cout2E(outfile)
            Edihe = totE - totE_nodihe
            
            totElist.append(totE)
            totE_nodihe_list.append(totE_nodihe)
            Edihe_list.append(Edihe)

        hetgraphfile = dataoutputpath_prefix / (prefix + f"-{clusnum}") / "heterograph.bin"
        heterograph = dgl.load_graphs(str(hetgraphfile))[0][0]
        heterograph.nodes["g"].data["u_charmm_totE"] = torch.tensor([totElist])
        heterograph.nodes["g"].data["u_charmm_totE_nodihe"] = torch.tensor([totE_nodihe_list])
        heterograph.nodes["g"].data["u_charmm_Edihe"] = torch.tensor([Edihe_list])

        dgl.save_graphs(str(hetgraphfile), [heterograph])


    shutil.rmtree(mol2outputpath_prefix)
    shutil.rmtree(crdoutputpath_prefix)
    shutil.rmtree(logpath_prefix)

    shutil.rmtree(tmppath)
