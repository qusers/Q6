#!/bin/bash

teLeap -f tleap.in \
       -I ff-amber14/cmd -I ff-amber14/lib \
       -I ff-amber14/parm -I ff-amber14/prep > tleap.out

sander -O -i relax.in -o relax.out -p mol.prmtop -c mol.inpcrd

echo "Energies:"
grep -m1 ANGLE -A2 relax.out

