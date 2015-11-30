#!/bin/python
load prep/lig_w.pdb, test
python
for i in range(1,6):
    cmd.load("eq%1d.dcd" % i, "test", discrete=0)
python end
python
for i in range(1,6):
    cmd.load("dc%1d.dcd" % i, "test", discrete=0)
python end

select waters, solvent
select ligand, not waters
hide lines, waters
show sticks, ligand
deselect
