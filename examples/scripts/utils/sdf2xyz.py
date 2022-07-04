from pathlib import Path

def sdf2xyz(sdffile):
    xyzs = []
    lines = Path(sdffile).read_text().split("\n")
    for i, line in enumerate(lines):
        linesplit = line.split()
        if len(linesplit) == 16: # atom xyz lines
            xyzs.append([float(linesplit[0]), float(linesplit[1]), float(linesplit[2])])
    
    return xyzs

