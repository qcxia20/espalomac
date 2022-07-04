strfile="/pubhome/qcxia02/work/espalomac/examples/datasets/mol2s/OpenFFGen2OptSet4eMoleculesDiscrepancyc#cc1cccc(c1)nc(=o)n2ccc3ccccc3c2/tmp0.str" \
crdfile="/pubhome/qcxia02/work/espalomac/examples/datasets/crds/OpenFFGen2OptSet4eMoleculesDiscrepancyc#cc1cccc(c1)nc(=o)n2ccc3ccccc3c2/OpenFFGen2OptSet4eMoleculesDiscrepancyc#cc1cccc(c1)nc(=o)n2ccc3ccccc3c2-0_0.crd"

cp $strfile ./test.str
cp $crdfile ./test.crd
/pubhome/qcxia02/soft/charmm/bin/charmm \
-i energy.inp \
workdir=`readlink -f .` \
resi=OpenFFGe \
cgenff:true \
modcutoff:true \
str=`readlink -f test.str` \
crd=`readlink -f test.crd`