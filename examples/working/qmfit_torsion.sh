#$ -S /bin/bash
#$ -cwd
#$ -o log
#$ -e log
#$ -j y
#$ -q benz
#$ -V
#$ -pe benz 32


conda activate espaloma

# python qmfit_torsion.py \
# --Nepochs 10 \
# --datasetpath "/pubhome/qcxia02/work/espalomac/examples/datasets/gen2/OpenFFGen2OptSet4eMoleculesDiscrepancyc[c@@h]1cn(cc[c@@]1(c)o)c(=o)cn2nccn2" \
# --learningrate "1e-3"

# python qmfit_torsion.py \
# --Nepochs 10000 \
# --datasetpath "/pubhome/qcxia02/work/espalomac/examples/datasets/modgen2" \
# --learningrate "1e-3"

# python qmfit_torsion.py \
# --Nepochs 100 \
# --datasetpath "/pubhome/qcxia02/work/espalomac/examples/datasets/TorNet500_mm_subset" \
# --learningrate "1e-3"

python qmfit_torsion.py \
--Nepochs 1000 \
--datasetpath "/pubhome/qcxia02/work/espalomac/examples/datasets/TorsionNet500_mm_opt_geometries" \
--learningrate "1e-3"
