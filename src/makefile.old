# makefile for the Q package:
# Qprep5/Qdyn5/Qdyn5p/Qfep5/Qdum5/Qcalc5
# version 5.01, 2003-08-26

###########################################################
#
# use make without parameters for build instructions
#
###########################################################

default: what


###########################################################
#
# common declarations
#
###########################################################

OBJ_FLG=-o
OUTDIR="./"

###########################################################
#
# supported platforms
#
###########################################################

#
# Cray T3E
#

#F90FLG=-O2,unroll2,pipeline2 -dp -lmfastv

cray-t3e:
	@make F90="f90" \
	F90FLG="-O1" \
	F90LIBS= \
	PP_SUFFIX=".F90" \
	CPP_FLG= \
	OBJ_FLG="-b"\
	all

cray-t3e-debug:
	@make F90="f90" \
	F90FLG="-O1" \
	F90LIBS= \
	PP_SUFFIX=".F90" \
	CPP_FLG= \
	OBJ_FLG="-b"\
	all

cray-t3e-mpi:
	@make F90="f90" \
	F90FLG="-O1" \
	F90LIBS= \
	USE_MPI= \
	PP_SUFFIX=".F90" \
	CPP_FLG= \
	OBJ_FLG="-b"\
	all


#
# Digital alpha/OSF1
#
alpha-osf1-mpi:
	@make F90="f90" \
	F90FLG="-O1" \
	F90LIBS="-L/usr/opt/MPI104/lib -lmpi -lfmpi -limc -lrt -threads" \
	USE_MPI= \
	CPP_FLG="-cpp" \
	all

alpha-osf1-ev5:
	@make F90="f90" \
	F90FLG="-arch ev5 -tune ev5 -O4 -fast -pipeline -speculate by_routine" \
	F90LIBS= \
	CPP_FLG="-cpp" \
	all

alpha-osf1-ev6:
	@make F90="f90" \
	F90FLG="-arch ev6 -tune ev6 -O4 -fast -pipeline -speculate by_routine" \
	F90LIBS= \
	CPP_FLG="-cpp" \
	all

alpha-osf1-ev5-i:
	@make F90="f90" \
	F90FLG="-arch ev5 -tune ev5 -O4 -fast -pipeline -speculate by_routine" \
	F90LIBS= \
	CPP_FLG="-cpp" \
	CCFLG="-arch ev5 -fast" \
	Qdyn5i

#
# Sun Solaris 5.5
#
solaris:
	@make F90="f90" \
	F90FLG="-fast -xlibmopt" \
	F90LIBS= \
	PP_SUFFIX=".F90" \
	CPP_FLG= \
	all

#
# IRIX 
#
irix6.4:
	@make F90="f90" \
	F90FLG="-O2 -OPT:Olimit=0,fast_sqrt=ON,rsqrt=ON -n32 -mips4 -TARG:platform=ip27" \
	F90LIBS= \
	CPP_FLG="-cpp" \
	all

irix6.3:
	@make F90="f90" \
	F90FLG="-O2 -OPT:fast_sqrt=ON:rsqrt=ON -n32 -TARG:processor=r4000" \
	F90LIBS= \
	CPP_FLG= \
	all

irix6.2:
	@make F90="f90" \
	F90FLG="-O2 -OPT:fast_sqrt=ON -n32 -mips3 -TARG:platform=ip22" \
	F90LIBS= \
	CPP_FLG="-cpp" \
	all

#
# Linux Portland pgi ia-32 compiler
# (by Peter Vagedes, FU Berlin)

linux_pgi:
	@cp md.f90 md.F90
	@cp qdyn.f90 qdyn.F90
	@cp qcalc.f90 qcalc.F90
	@cp qfep.f90 qfep.F90
	@cp prep.f90 prep.F90
	@make F90="pgf77" \
	F90FLG="-fast" \
	CPP_FLG=  \
	F90LIBS= \
	CCFLAG="-O" \
	CC="cc" \
	PP_SUFFIX=".F90" \
	all

linux_abs:
	@cp md.f90 md.F90
	@cp qdyn.f90 qdyn.F90
	@cp qcalc.f90 qcalc.F90
	@cp qfep.f90 qfep.F90
	@cp prep.f90 prep.F90
	@make F90="f95 -O2 -cpu:athlon " \
	F90FLG="-DNO_FLUSH" \
	CPP_FLG=  \
	F90LIBS="-lU77" \
	CCFLAG="-O" \
	CC="cc" \
	PP_SUFFIX=".F90" \
	Qdyn5

#
#IBM AIX (UNIX) 
#

#For parallel execution on IBM SP2 at PDC with 'xlf90' compiler
power2-rs/6000:
	@make \
	CPP = -qsuffix=cpp=f90 \
	CPPFLAG= -WF,-D \
	F90FLG=-O2 -qsuffix=f=f90 -qintlog\
	F90=mpxlf90 \
	all

#
#Linux Red hat 7.3 on IA-32, NSC in Linköping
#


ifc:
	@cp md.f90 md.F90
	@cp qdyn.f90 qdyn.F90
	@cp qcalc.f90 qcalc.F90
	@cp qfep.f90 qfep.F90
	@cp prep.f90 prep.F90
	@make F90="ifc" \
	F90FLG="-w -DNO_FLUSH" \
	F90LIBS="-Vaxlib" \
	CPP_FLG="-cpp" \
	CCFLAG="-O" \
	CC="cc" \
	PP_SUFFIX="F90" \
	all

ifort_unopt:
	@cp md.f90 md.F90
	@cp qdyn.f90 qdyn.F90
	@cp qcalc.f90 qcalc.F90
	@cp qfep.f90 qfep.F90
	@cp prep.f90 prep.F90
	@make F90="ifort" \
	F90FLG="-w -DNO_FLUSH -O0 -g -debug extended -static -DPROFILING" \
	F90LIBS="-Vaxlib" \
	CPP_FLG="-cpp" \
	CCFLAG="-O" \
	CC="cc" \
	PP_SUFFIX="F90" \
	Qcalc5


pgf90:
	@cp md.f90 md.F90
	@cp qdyn.f90 qdyn.F90
	@cp qcalc.f90 qcalc.F90
	@cp qfep.f90 qfep.F90
	@cp prep.f90 prep.F90
	@make F90="pgf90" \
	F90FLG="-w -DNO_FLUSH -fast" \
	F90LIBS= \
	CPP_FLG= \
	CCFLAG="-O" \
	CC="cc" \
	PP_SUFFIX="F90" \
	all




ifc_scampi:
	@cp md.f90 md.F90
	@cp qdyn.f90 qdyn.F90
	@cp qcalc.f90 qcalc.F90
	@cp qfep.f90 qfep.F90
	@cp prep.f90 prep.F90
	@make F90="ifc" \
	F90FLG="-w -DNO_FLUSH -xW -ip -ipo -unroll -tpp7 -O3 -Nscampi" \
	F90LIBS="-Vaxlib" \
	CPP_FLG="-cpp" \
	CCFLAG="-O" \
	CC="cc" \
	PP_SUFFIX="F90" \
	Qdyn5p



ifc_scampidbg:
	@cp md.f90 md.F90
	@cp qdyn.f90 qdyn.F90
	@cp qcalc.f90 qcalc.F90
	@cp qfep.f90 qfep.F90
	@cp prep.f90 prep.F90
	@make F90="ifc" \
	F90FLG="-w -DNO_FLUSH -g -Nscampi" \
	F90LIBS="-Vaxlib" \
	CPP_FLG="-cpp" \
	CCFLAG="-O" \
	CC="cc" \
	PP_SUFFIX="F90" \
	Qdyn5p



ifc_mpich:
	@cp md.f90 md.F90
	@cp qdyn.f90 qdyn.F90
	@cp qcalc.f90 qcalc.F90
	@cp qfep.f90 qfep.F90
	@cp prep.f90 prep.F90
	@make F90="ifc -w -DMPICH -xW -ip -ipo -unroll -tpp7 -O3" \
	F90FLG="-w -DNO_FLUSH" \
	F90LIBS="-Vaxlib" \
	CPP_FLG="-cpp" \
	CCFLAG="-O" \
	CC="cc" \
	PP_SUFFIX="F90" \
	Qdyn5p

#
# Other...
#


generic:
	@make F90="f90" \
	F90FLG="-O" \
	F90LIBS= \
	CPP_FLG="-cpp" \
	all

#
# debug
#

debug:
	@make F90="f90" \
	F90FLG="-g" \
	F90LIBS= \
	CPP_FLG="-cpp" \
	all

###########################################################
#
# pseudo-targets
#
###########################################################

all:	Qfep5 Qprep5 Qdyn5 Qdum5 Qcalc5 

clean:
	rm -f *.o *.F90 *.mod *.M *.kmo *.il

nuke:
	rm -f *.o *.F90 *.mod *.M *.kmo *.il Qfep5 Qdyn5p Qdyn5 Qprep5 Qcalc5 Qdum5


${OUTDIR}:
	mkdir ${OUTDIR}

###########################################################
#
# real build targets: programs
#
###########################################################

Qfep5: qfep.o mpiglob.o nrgy.o misc.o parse.o
	${F90} ${F90FLG} -o Qfep5 qfep.o mpiglob.o nrgy.o misc.o parse.o ${F90LIBS}

Qprep5: q_prep.o topo.o misc.o mpiglob.o parse.o prefs.o prep.o prmfile.o index.o mask.o trj.o sizes.o avetr.o
	${F90} ${F90FLG} -o Qprep5 q_prep.o topo.o maskmanip.o misc.o mpiglob.o parse.o prefs.o prep.o prmfile.o index.o  mask.o trj.o sizes.o avetr.o ${F90LIBS}

Qcalc5: calc_base.o calc_chemscore.o calc_fit.o calc_geom.o calc_pmfscore.o calc_com_ke.o calc_com.o calc_rdf.o calc_rms.o calc_rmsf.o calc_entropy.o calc_nb.o calc_xscore.o eigen.o index.o mask.o maskmanip.o misc.o mpiglob.o nrgy.o parse.o prmfile.o qatom.o qcalc.o sizes.o topo.o trj.o
	${F90} ${F90FLG} -o Qcalc5 calc_base.o calc_chemscore.o calc_fit.o calc_geom.o calc_pmfscore.o calc_com_ke.o calc_com.o calc_rdf.o calc_rms.o calc_rmsf.o calc_entropy.o calc_nb.o calc_xscore.o eigen.o index.o mask.o maskmanip.o misc.o mpiglob.o nrgy.o parse.o prmfile.o qatom.o qcalc.o sizes.o topo.o trj.o ${F90LIBS}

Qdyn5:	md.o mask.o misc.o mpiglob.o nrgy.o prmfile.o qatom.o qdyn.o sizes.o topo.o trj.o index.o
	${F90} ${F90FLG} -o Qdyn5 md.o mask.o misc.o mpiglob.o nrgy.o prmfile.o qatom.o qdyn.o sizes.o topo.o trj.o index.o ${F90LIBS}

Qdum5: md_dum.o mask.o misc.o mpiglob.o nrgy.o prmfile.o qatom.o qdyn.o sizes.o topo.o trj.o index.o
	${F90} ${F90FLG} -o Qdum5 md_dum.o mask.o misc.o mpiglob.o nrgy.o prmfile.o qatom.o qdyn.o sizes.o topo.o trj.o index.o ${F90LIBS}

Qdyn5p: md_mpi.o mask.o misc.o mpiglob.o nrgy.o prmfile.o qatom.o qdyn_mpi.o sizes.o topo.o trj.o index.o
	${F90} ${F90FLG} mask.o misc.o mpiglob.o nrgy.o prmfile.o qatom.o qdyn_mpi.o sizes.o topo.o trj.o index.o md_mpi.o -o Qdyn5p ${F90LIBS}


###########################################################
#
# object modules
#
###########################################################

avetr.o: avetr.f90 prep.o
	${F90} ${F90FLG} -c avetr.f90

calc_base.o:calc_base.f90 topo.o
	${F90} ${F90FLG} -c calc_base.f90

calc_chemscore.o: calc_chemscore.f90 maskmanip.o trj.o prmfile.o index.o qatom.o
	${F90} ${F90FLG} -c calc_chemscore.f90

calc_entropy.o:calc_entropy.f90 calc_base.o maskmanip.o trj.o calc_fit.o
	${F90} ${F90FLG} -c calc_entropy.f90

calc_fit.o:calc_fit.f90 calc_base.o maskmanip.o
	${F90} ${F90FLG} -c calc_fit.f90

calc_geom.o:calc_geom.f90 calc_base.o
	${F90} ${F90FLG} -c calc_geom.f90

calc_nb.o:calc_nb.f90 calc_base.o maskmanip.o parse.o
	${F90} ${F90FLG} -c calc_nb.f90

calc_pmfscore.o: calc_pmfscore.f90 calc_base.o maskmanip.o trj.o topo.o prmfile.o index.o qatom.o misc.o
	${F90} ${F90FLG} -c calc_pmfscore.f90 

calc_rdf.o:calc_rdf.f90 calc_base.o parse.o maskmanip.o
	${F90} ${F90FLG} -c calc_rdf.f90

calc_rms.o:calc_rms.f90 calc_base.o maskmanip.o
	${F90} ${F90FLG} -c calc_rms.f90

calc_com_ke.o:calc_com_ke.f90 calc_base.o maskmanip.o
	${F90} ${F90FLG} -c calc_com_ke.f90

calc_com.o:calc_com.f90 calc_base.o maskmanip.o
	${F90} ${F90FLG} -c calc_com.f90

calc_rmsf.o:calc_rmsf.f90 calc_base.o maskmanip.o
	${F90} ${F90FLG} -c calc_rmsf.f90

calc_xscore.o: calc_xscore.f90 calc_base.o maskmanip.o trj.o topo.o prmfile.o index.o qatom.o misc.o
	${F90} ${F90FLG} -c calc_xscore.f90

eigen.o:eigen.f90
	${F90} ${F90FLG} -c eigen.f90

invsqrt_q.o:invsqrt_q.c
	$(CC) $(CCFLG) -c invsqrt_q.c

index.o:index.f90
	${F90} ${F90FLG} -c index.f90

mask.o:	mask.f90 topo.o
	${F90} ${F90FLG} -c mask.f90

maskmanip.o:maskmanip.f90 mask.o misc.o parse.o
	${F90} ${F90FLG} -c maskmanip.f90

md.o:	md.f90 mpiglob.o qatom.o sizes.o trj.o topo.o 
	${F90} ${F90FLG} ${CPP_FLG} -c md.F90

md_dum.o:md.f90 mpiglob.o qatom.o sizes.o topo.o
	${F90} ${F90FLG} ${CPP_FLG} -DDUM -c ${OBJ_FLG} md_dum.o md.F90

md_mpi.o: md.f90 mpiglob.o qatom.o sizes.o topo.o trj.o
	${F90} ${F90FLG} ${CPP_FLG} -DUSE_MPI -c ${OBJ_FLG} md_mpi.o md.F90

misc.o: misc.f90 sizes.o
	${F90} ${F90FLG} -c misc.f90

mpiglob.o: mpiglob.f90 sizes.o nrgy.o
	${F90} ${F90FLG} -c mpiglob.f90

nrgy.o: nrgy.f90 sizes.o
	${F90} ${F90FLG} -c nrgy.f90

parse.o: parse.f90 misc.o
	${F90} ${F90FLG} -c parse.f90

prefs.o: prefs.f90
	${F90} ${F90FLG} -c prefs.f90

prep.o: maskmanip.o sizes.o parse.o prmfile.o trj.o index.o prefs.o prep.f90
	${F90} ${F90FLG} -c prep.f90

prmfile.o: prmfile.f90 misc.o mpiglob.o
	${F90} ${F90FLG} -c prmfile.f90

q_prep.o: q_prep.f90 prep.o avetr.o
	${F90} ${F90FLG} -c q_prep.f90

qatom.o: qatom.f90 misc.o nrgy.o prmfile.o sizes.o index.o topo.o
	${F90} ${F90FLG} -c qatom.f90

qcalc.o: qcalc.f90 calc_chemscore.o calc_pmfscore.o calc_xscore.o trj.o calc_base.o calc_rms.o calc_fit.o calc_geom.o 
	${F90} ${F90FLG} ${CPP_FLG} -c qcalc.F90

qdyn.o: qdyn.f90 md.o mpiglob.o
	${F90} ${F90FLG} ${CPP_FLG} -c qdyn.F90

qdyn_dum.o: qdyn.f90 md.o mpiglob.o
	${F90} ${F90FLG} ${CPP_FLG} -DDUM -c ${OBJ_FLG} qdyn_dum.o qdyn.F90

qdyn_mpi.o: qdyn.f90 md_mpi.o mpiglob.o
	${F90} ${F90FLG} ${CPP_FLG} -DUSE_MPI -c ${OBJ_FLG} qdyn_mpi.o qdyn.F90

qfep.o: qfep.f90 nrgy.o parse.o
	${F90} ${F90FLG} ${CPP_FLG} -c qfep.F90

sizes.o: sizes.f90
	${F90} ${F90FLG} -c sizes.f90

topo.o:  topo.f90 misc.o mpiglob.o sizes.o
	${F90} ${F90FLG} -c topo.f90

trj.o:  trj.f90 mask.o misc.o 
	${F90} ${F90FLG} -c trj.f90


###########################################################
#
# default target: build instructions
#
###########################################################

what:
	@echo "usage: make <target>"
	@echo
	@echo "<target> is one of:"
	@echo "cray-t3e         Cray UNICOS (not tested recently)"
	@echo "alpha-osf1-ev5   Tru64 UNIX (Digital UNIX) v3.2 and up on EV5(21164) CPU"
	@echo "alpha-osf1-ev6   Tru64 UNIX (Digital UNIX) v4.0 and up on EV6(21264) CPU"
	@echo "alpha-osf1-mpi   MPI parallell on Tru64 UNIXv4.0 and up on EV5 (not tested)"
	@echo "alpha-osf1-ev5-i Tru64 UNIX v3.2+ on EV5 + use C routine for sqrt(1/x)"
	@echo "irix6.4"
	@echo "irix6.2"
	@echo "solaris"
	@echo "linux_pgi        Linux with Portland pgi ia-32 compiler"
	@echo "linux_absoft     Linux with Absoft compiler athlon opt "
	@echo "power2-rs/6000   AIX4.3 (UNIX) on IBM SP2"
	@echo "ifc              Linux with Intel compiler"
	@echo "ifc_scampi       Parallel version of Qdyn5 on Linux with Intel compiler and Scampi MPI"
	@echo "generic"
	@echo "debug"


