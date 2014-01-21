# makefile for the Q package:
# Qprep5/Qdyn5/Qdyn5p/Qfep5/Qdum5/Qcalc5
# version 5.01, 2003-08-26

SUFFIX=.f90 .F90

FC = ifort
MPIFC = mpiifort
FCFLAGS = -O0 -g -check uninit -traceback
LD = $(FC)
MPILD = $(MPIFC)
LDFLAGS = -traceback

QfepSource = mpiglob.f90 nrgy.f90 misc.f90 parse.f90
QfepFSource = qfep.F90

QprepSource = q_prep.f90 topo.f90 misc.f90 mpiglob.f90 parse.f90 prefs.f90 prep.f90 prmfile.f90 index.f90 mask.f90 trj.f90 sizes.f90 avetr.f90 maskmanip.f90 nrgy.f90

QcalcSource = calc_base.f90 calc_chemscore.f90 calc_fit.f90 calc_geom.f90 calc_pmfscore.f90 calc_com_ke.f90 calc_com.f90 calc_rdf.f90 calc_rms.f90 calc_rmsf.f90 calc_entropy.f90 calc_nb.f90 calc_xscore.f90 eigen.f90 index.f90 mask.f90 maskmanip.f90 misc.f90 mpiglob.f90 nrgy.f90 parse.f90 prmfile.f90 qatom.f90 qcalc.f90 sizes.f90 topo.f90 trj.f90

QdynFSource = md.F90 qdyn.F90
QdynSource = mask.f90 misc.f90 mpiglob.f90 nrgy.f90 prmfile.f90 qatom.f90 sizes.f90 topo.f90 trj.f90 index.f90

QdumFSource = qdyn.F90
QdumSource = mask.f90 misc.f90 mpiglob.f90 nrgy.f90 prmfile.f90 qatom.f90 sizes.f90 topo.f90 trj.f90 index.f90

QdynpSource = mask.f90 misc.f90 mpiglob.f90 nrgy.f90 prmfile.f90 qatom.f90 sizes.f90 topo.f90 trj.f90 index.f90

QfepObjects = $(QfepSource:.f90=.o) $(QfepFSource:.F90=.o)
QprepObjects = $(QprepSource:.f90=.o)
QcalcObjects = $(QcalcSource:.f90=.o)
QdynObjects = $(QdynSource:.f90=.o) $(QdynFSource:.F90=.o)
QdumObjects = $(QdumSource:.f90=.o) $(QdumFSource:.F90=.o) md_dum.o
QdynpObjects = $(QdynpSource:.f90=.o) md_mpi.o qdyn_mpi.o

PROGRAMS = Qfep5 Qdyn5p Qdyn5 Qprep5 Qcalc5 Qdum5

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

%.o: %.F90
	$(FC) $(FCFLAGS) -c $<

###########################################################
#
# pseudo-targets
#
###########################################################

all:	$(PROGRAMS)

clean:
	rm -f *.o *.mod *.M *.kmo *.il $(PROGRAMS)

###########################################################
#
# real build targets: programs
#
###########################################################

Qfep5: $(QfepObjects)
	$(LD) $(LDFLAGS) -o $@ $(QfepObjects) $(FLIBS)

Qprep5: $(QprepObjects)
	$(LD) $(LDFLAGS) -o $@ $(QprepObjects) $(FLIBS)

Qcalc5: $(QcalcObjects)
	$(LD) $(LDFLAGS) -o $@ $(QcalcObjects) $(FLIBS)

Qdyn5:	$(QdynObjects)
	$(LD) $(LDFLAGS) -o $@ $(QdynObjects) $(FLIBS)

Qdum5: $(QdumObjects)
	$(LD) $(LDFLAGS) -o $@ $(QdumObjects) $(FLIBS)

Qdyn5p: $(QdynpObjects)
	$(MPILD) $(LDFLAGS) -o $@ $(QdynpObjects) $(FLIBS)

# Some objects need special handling

md_dum.o: md.F90
	$(FC) $(FCFLAGS) -DDUM -c $< -o $@

md_mpi.o: md.F90
	$(MPIFC) $(FCFLAGS) -DUSE_MPI -c $< -o $@

qdyn_mpi.o: qdyn.F90
	$(MPIFC) $(FCFLAGS) -DUSE_MPI -c $< -o $@

# Include the fortran module deps
include mod-deps

# Add deps for the obejcts with special handliung
md_mpi.o: mpiglob.o qatom.o sizes.o trj.o
md_dum.o: mpiglob.o qatom.o sizes.o trj.o
qdyn_mpi.o: md_mpi.o mpiglob.o

