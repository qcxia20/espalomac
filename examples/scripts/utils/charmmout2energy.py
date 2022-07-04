#!/pubhome/qcxia02/miniconda3/envs/espaloma/bin/python3.8

def Cout2E(outfile):
    lines = outfile.read_text().split("\n")
    Elineidx = 0
    for i, line in enumerate(lines):
        if line.startswith("ENER>"):
            Elineidx = i
    
    Elines = lines[Elineidx:Elineidx+3]
    totE = float(Elines[0].split()[2])
    diheE = float(Elines[1].split()[5])
    totE_nodihe = totE-diheE
    
    return totE, totE_nodihe



if __name__ == "__main__":
    import sys
    from pathlib import Path
    outfile = Path(sys.argv[1])
    totE, totE_nodihe = Cout2E(outfile)



