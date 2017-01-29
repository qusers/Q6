#!/usr/bin/env python2

import sys
import logging

try:
    from Qpyl.core.qparameter import QPrm
    from Qpyl.core.qlibrary import QLib
    from Qpyl.core.qstructure import QStruct
    from Qpyl.core.qtopology import QTopology
    from Qpyl.common import SpecialFormatter
except ImportError:
    print "QTools not installed"
    sys.exit(1)

logger = logging.getLogger('Qpyl')
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler(sys.stdout)
formatter = SpecialFormatter()
handler.setFormatter(formatter)
logger.addHandler(handler)

def log_dups(dups):
    for dup in dups:
        logger.info("Overwritten: {}".format(dup))

#
# Amber14FF to Qamber14
#
# Convert Amber14 lib (+prepin for impropers) and parm+frcmod to Q lib/prm
# Load all_amino_acids.pdb and all_amino_acids.mol2 and build topology
# separately.
#
# ignore_errors=True is needed because frcmod.ff14SB overwrites some parameters
#
qal = QLib("amber")
qap = QPrm("amber", ignore_errors=True)
qal.read_amber_lib("../Amber/ff-amber14/amber12_mod.lib")
qal.read_amber_lib("../Amber/ff-amber14/arn.lib")
qal.read_prepin_impropers("../Amber/ff-amber14/prep/amino12.in")
qal.read_prepin_impropers("../Amber/ff-amber14/arn.prepi")
log_dups(qap.read_amber_parm("../Amber/ff-amber14/parm/parm10.dat"))
log_dups(qap.read_amber_parm("../Amber/ff-amber14/parm/parm10.dat"))
log_dups(qap.read_amber_frcmod("../Amber/ff-amber14/parm/frcmod.ff14SB"))



# add options to parameters
for line in """name		Q-Amber14SB
type		AMBER
vdw_rule    arithmetic !vdW combination rule (geometric or arithmetic)
scale_14	0.8333 ! electrostatic 1-4 scaling factor
switch_atoms	off
improper_potential	periodic
improper_definition explicit""".splitlines():
    lf = line.split()
    qap.options[lf[0]] = " ".join(lf[1:])

# remove head from ACE and tail from NME
cons = qal.residue_dict["ACE"].connections
cons = [con for con in cons if "head" not in con]
qal.residue_dict["ACE"].connections = cons

cons = qal.residue_dict["NME"].connections
cons = [con for con in cons if "tail" not in con]
qal.residue_dict["NME"].connections = cons


open("qamber14_gen.lib", "w").write(qal.get_string())
open("qamber14_gen.prm", "w").write(qap.get_string())


print "# Topology with converted parameters (and mol2):"
qas1 = QStruct("../all_amino_acids.mol2", "mol2")
qat = QTopology(qal, qap, qas1)
q_tors = sum([len(list(tor.prm.get_prms())) for tor in qat.torsions])
print "Bonds: ", len(qat.bonds)
print "Angles: ", len(qat.angles)
print "Torsions: ", len(qat.torsions)
print "Q Torsions (diff parms): ", q_tors
print "Impropers: ", len(qat.impropers)
print "Bond energy: ", sum([bond.calc()[0] for bond in qat.bonds])
print "Angle energy: ", sum([ang.calc()[0] for ang in qat.angles])
print "Torsion energy: ", sum([tor.calc()[0] for tor in qat.torsions])
print "Improper energy: ", sum([imp.calc()[0] for imp in qat.impropers])


print
print "# Topology with pre-made parameters (and pdb):"
qal2 = QLib("amber")
qap2 = QPrm("amber")
qas2 = QStruct("../all_amino_acids.pdb", "pdb", ignore_errors=True)
qal2.read_lib("../Q/ff-qamber14/qamber14.lib")
log_dups(qap2.read_prm("../Q/ff-qamber14/qamber14.prm"))
qat = QTopology(qal2, qap2, qas2)
q_tors = sum([len(list(tor.prm.get_prms())) for tor in qat.torsions])
print "Bonds: ", len(qat.bonds)
print "Angles: ", len(qat.angles)
print "Torsions: ", len(qat.torsions)
print "Q Torsions (diff parms): ", q_tors
print "Impropers: ", len(qat.impropers)
print "Bond energy: ", sum([bond.calc()[0] for bond in qat.bonds])
print "Angle energy: ", sum([ang.calc()[0] for ang in qat.angles])
print "Torsion energy: ", sum([tor.calc()[0] for tor in qat.torsions])
print "Improper energy: ", sum([imp.calc()[0] for imp in qat.impropers])


