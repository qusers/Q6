#!/bin/bash

sed 's/ \(\S\S\S\S\) \(\S\S\S\) /  \1\2 /' ../all_amino_acids.pdb > mol.pdb

qprep5 < prep.inp > prep.log 
qdyn5_r8 relax.inp > relax.log

if [[ $(grep "PDB file successfully written" prep.log) && $(grep "terminated normally" relax.log) ]]
then
    grep -m1 -A3 "Energy summary at step" relax.log
else
    echo "FAILED"
fi



