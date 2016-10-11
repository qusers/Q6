### Amber ff14SB modification

This folder contains the ff14sb force field as found in AmberTools16 (the force
field is according to their license, in the public domain).

Additionally, several files have been created in the root directory:

**amber12_mod.lib** is the modified amber12.lib, with added ACE and NME groups
from amber12ct.lib and amber12nt.lib

**arn_94.prepi** (Neutral arginine in amber94) was obtained from:
http://signe.teokem.lu.se/~ulf/Methods/arn.html
http://pubs.acs.org/doi/suppl/10.1021/jm0608210

arn.prepi is arn_94.prepi with modified atom_types to match those found in
amino12.in (ff14SB) for ARG. The charges of library entries in force fields 
94 and 14 are the same, so they do not need modification.

**arn.lib** was created with tleap (AmberTools16):

> source leaprc.protein.ff14SB
> loadAmberPrep arn.prepi
> saveOff ARN arn.lib


