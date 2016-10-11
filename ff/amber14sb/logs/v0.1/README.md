### QAmber14SB v0.1

The whole *amino12.lib* Amber library, plus `ACE` and `NME` (from
*amino12nt.lib* and *amino12ct.lib*), as well as all the parameters 
found in *parm10.dat* and *ff14sb.frcmod* in *AmberTools16* were converted to Q
format with the use of QTools. 

The following residues are included:
```
ACE ALA ARG ASH ASN ASP 
CYM CYS CYX 
GLH GLN GLU GLY 
HID HIE HIP HYP 
ILE LEU LYN LYS 
MET 
NHE NME 
PHE PRO 
SER 
THR TRP TYR 
VAL
```


To test the parameters, a chain composed of 28 different amino acids, capped
with ACE and NME was constructed. Single point gas phase energies of Amber, Q
and QTools (both, in-situ converted and pre-made Q force fields are used) are 
then calculated and different energy contributions can be compared.

Additionally, neutral Arginine `ARN` and TIP3P water `HOH` were added manually.  
(see *conv_test/Amber/ff-amber14/README.md*)

