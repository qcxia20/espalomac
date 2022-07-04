#!/pubhome/qcxia02/miniconda3/envs/espaloma/bin/python3.8
import argparse
from pathlib import Path, PosixPath
import sys

def next_mol2_lines(infile):
    """Method to return one mol2 block once."""
    lines = list()

    for line in open(infile):
        if "@<TRIPOS>MOLECULE" in line:
            if len(lines) == 0:
                lines.append(line)
            else:
                yield lines
                lines = list()
                lines.append(line)
        else:
            lines.append(line)

    yield lines

def mol2tocrdlines(mol2lines):
    atomsum = int(mol2lines[2].split()[0])
    resname = mol2lines[1].split("\n")[0]
    if len(resname) > 8:
        resname = "lig"

    crdlines = []
    for i,line in enumerate(mol2lines[6:6+atomsum]):
        atomname = line.split()[1]
        x,y,z = float(line.split()[2]), float(line.split()[3]), float(line.split()[4])
        crdline = f"{i+1:10d}{1:10d}  {resname:8s}  {atomname:8s}{x:20.10f}{y:20.10f}{z:20.10f}  {resname:8s}  {'1':8s}{0:20.10f}"
        crdlines.append(crdline) 
    
    return crdlines, atomsum

def writecrdlines(outputfilepathname, crdlines, atomsum):
    with open(outputfilepathname, 'w') as f:
        f.write("*Created by mol2tocrd\n*\n") # title
        f.write("%10d" % atomsum + "  EXT\n") # NATOM
        f.write("\n".join(crdlines)) # coordinate content
        f.write("\n")

"""
See https://charmm-gui.org/charmmdoc/io.html for more info about .crd file format
title
NATOM (I10)
ATOMNO RESNO   RES  TYPE  X     Y     Z   SEGID RESID Weighting
  I10   I10 2X A8 2X A8       3F20.10     2X A8 2X A8 F20.10
"""


if __name__ == "__main__":
    parser = argparse.ArgumentParser("convert mol2 to crd, note that resname <= 8")
    parser.add_argument("--input", type=PosixPath, help="input mol2 file path, absolute path is better")
    parser.add_argument("--outputpath", type=PosixPath, help="crd output path without file name")
    args = parser.parse_args()

    mol2file = args.input
    molname = str(mol2file).split("/")[-2]
    outpath = args.outputpath

    allmol2lines = next_mol2_lines(mol2file)
    allmol2lineslist = list(allmol2lines)
    # print(len(list(allmol2lines)))
    if not (outpath / molname).exists():
        (outpath / molname).mkdir()

    if len(allmol2lineslist) == 1:
        mol2lines = allmol2lineslist[0]
        crdlines, atomsum = mol2tocrdlines(mol2lines)
        
        outname = mol2file.name.split(".mol2")[0] + ".crd"
        writecrdlines(outpath / molname / outname, crdlines, atomsum)
    
    else:
        for i, mol2lines in enumerate(allmol2lineslist):
            crdlines, atomsum = mol2tocrdlines(mol2lines)
            
            outname = mol2file.name.split(".mol2")[0] + f"_{i}" + ".crd"
            writecrdlines(outpath / molname / outname, crdlines, atomsum)
