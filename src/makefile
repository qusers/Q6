################################################################################
#  Q V.6.0 Makefile                                                            #
#  code authors: Johan Aqvist, Martin Almlof, Martin Ander, Jens Carlson,      #
#  Isabella Feierberg, Peter Hanspers, Anders Kaplan, Karin Kolmodin,          #
#  Kajsa Ljunjberg, John Marelius, Martin Nervall                              #
#  keeper of the code:                                                         #
#  maintainers: Beat Amrein, Alexandre Barrozo                                 #
#  latest update: May 31, 2017                                                 #
#  make [option] [COMP=compiler]                                               #
################################################################################


################################################################################
#  Q Build Auto Increment                                                      #
################################################################################

BUILD_NUMBER_FILE=build-number.txt

################################################################################
# Intel Fortran ifort ( tested in linux and osx)
################################################################################
ifeq ($(COMP),ifort)
	PREC=		-DQDOUBLE
	FC=			ifort 
	FC_OPT=		-O3 -ip -ipo -xHost -inline-forceinline \
				-unroll-aggressive -finline-limit=1000\
				-no-prec-div -no-inline-max-per-compile \
				-no-inline-max-per-routine \
				-no-inline-max-total-size \
				-auto
	FFLAGS=		${PREC}
	FFLAGS_R8=	
	FC_DEBUG=	-unroll 
	FPP_FLG=	-fpp -DINTEL

	MPIFC=		mpif90
	MPIFLAGS=	-f90=${FC} -Nmpi
	OMPI_FC=	${FC}
	MPI_FC=		${FC}
	MPI_FCFLAGS=	${FFLAGS} ${FC_OPT}
	MPI=		${MPI_FCFLAGS}
	MPFLAG=		-openmp ${MPI}

	DEBUG=          -debug -traceback -g -O0  ${FFLAGS} -DDEBUG
	DEBUGMPI=       -debug -traceback -g -O0  ${FFLAGS} -DDEBUG
	DEBUGMP=        -openmp -debug -traceback -g -O0 -I${VT_ROOT}/include -L${VT_LIB_DIR} ${VT_ADD_LIBS} ${FFLAGS} -DDEBUG


################################################################################
# OSX Mavericks using gfortran
################################################################################
else ifeq ($(COMP),osx)
	PREC=		-DQDOUBLE
	FC= gfortran
	FC_OPT=         -O3
	FFLAGS=         -fstack-protector -DG95=1 \
                        -ffree-line-length-none -fcray-pointer -fall-intrinsics\
                        -std=legacy -Wall -Wtabs ${PREC}
	FFLAGS_R8=	-fdefault-real-8
	FPP_FLG=        -cpp -DOSX

	MPIFC=          mpif90
	OMPI_FC=        ${FC}
	OMPI_FCFLAGS=   ${FFLAGS} ${FC_OPT}
	MPI=		${OMPI_FCFLAGS}

	DEBUG=          -debug  -g ${FFLAGS}
	DEBUGMPI=       -debug  -g ${FFLAGS}

################################################################################
# PGI Compiler -- compiles at csb
################################################################################
else ifeq ($(COMP),pgi)
	PREC=		-DQDOUBLE
	FC=		pgf90
	FC_OPT=		-fastsse -Mipa=fast,inline -Minline=lrf_update,size:10000,levels:2 -O4 -Munroll -mcmodel=medium -Minfo=inline
	FFLAGS=		${PREC}
	FFLAGS_R8=	-r8
	FPP_FLG=	-Mpreprocess -DPGI

	MPIFC=		mpif90
	MPIFLAGS=	-lmpi
	MPI=		${MPIFLAGS} ${FFLAGS} ${FC_OPT}
	MPI_FC=		${FC}
	MPFLAG=		-mp ${MPI}

################################################################################
# IBM Compiler -- untested
################################################################################
else ifeq ($(COMP),ibm)
	PREC=		-DQDOUBLE
	FC=             xlf90
	FC_OPT=         -qsmp=omp
	FFLAGS=		-qrealsize=8 -DIBM ${PREC}

################################################################################
# GNU gcc-gfortran
################################################################################
else ifeq ($(COMP),gcc)
	PREC=		-DQDOUBLE
	FC=		gfortran
	FC_OPT=		-Ofast -unroll -finline-limit=10000000 -flto 
	FPP_FLG=	-cpp -DGCC

	MPIFC=		mpif90
	OMPI_FC=	${FC}
	OMPI_FCFLAGS=	-DG95=1 \
                        -ffree-line-length-none -fcray-pointer -fall-intrinsics -std=legacy 
	MPI=		${OMPI_FCFLAGS} ${FC_OPT} ${PREC}
	MPFLAG=		-fopenmp ${MPI}

	FFLAGS=		${OMPI_FCFLAGS} ${PREC}
	FFLAGS_R8=	
	DEBUG_FC=	${PREC}



	DEBUG=          -Og -fbacktrace -g ${DEBUG_FC}	-ffree-line-length-none -fcray-pointer -fall-intrinsics -std=legacy -DDEBUG -Wall -Wno-tabs -Wno-character-truncation
	DEBUGMPI=       -Og -fbacktrace -g ${DEBUG_FC}      -ffree-line-length-none -fcray-pointer -fall-intrinsics -std=legacy  -fstack-check -fstack-protector-all -DDEBUG -Wall -Wno-tabs -Wno-character-truncation
	DEBUGMP=        -Og -fopenmp -fbacktrace -g ${DEBUG_FC}      -ffree-line-length-none -fcray-pointer -fall-intrinsics -std=legacy  -fstack-check -fstack-protector-all -DDEBUG

	PROF=		-pg ${DEBUG}

else
	PREC=		-DQDOUBLE
	FC=		gfortran
	FC_OPT=		-O3 -unroll 
	FPP_FLG=	-cpp -DUNKNOWN

	MPI=		-DG95=1 \
                        -ffree-line-length-none -fcray-pointer -fall-intrinsics -std=legacy
	MPFLAG=		-fopenmp ${MPI}
	FFLAGS=		--static ${MPI} ${PREC}
	FFLAGS_R8=	
endif

################################################################################
# Display options
################################################################################
default: what

################################################################################
# Targets
################################################################################
all:    tools serial 
	

tools:
	@make FFLAGS="${FFLAGS} ${FC_OPT}" \
	qpi qfep qprep qdum qcalc qmin move1 move2

serial:	
	@make FFLAGS="${FFLAGS} ${FFLAGS_R8} ${FC_OPT}" \
	qdyn move1 moveqdyn

debug:
	@make FFLAGS="${DEBUG}" \
	qpi qdyn qfep qprep qdum qcalc qmin move1 move2 moveqdyn

prof:
	@make FFLAGS="${PROF}" \
	qfep qprep qdyn qdum qcalc move1 move2

mpi:
	@make FFLAGS="${MPI}" \
	qdynp qpip qminp move1 move3

mpidebug:
	@make FFLAGS="${DEBUGMPI}" \
	qdynp qpip qminp move1 move3 
hybrid:
	@make FFLAGS="${MPFLAG}" \
	qdynh move1 move4

hybriddebug:
	@make FFLAGS="${DEBUGMP}" \
	qdynh move1 move4

clean:
	-rm -f *.o  *.mod *.M *.kmo *.il *.oo

nuke:
	-rm -rf *.o  *.mod *.M *.kmo *.il *.oo qfep qdynh qdynp qdyn qprep qcalc qdum qmin qminp qpi qpip ../bin ../obj

qcalc qprep qfep: mpiglob.o

qcalc qdyn qdum qdynp qdynh qprep qfep qpi qpip qmin qminp: version.o misc.o math.o

qcalc qdyn qdum qdynp qdynh qprep qpi qpip qmin qminp: mask.o prmfile.o sizes.o topo.o trj.o index.o

qcalc qdyn qdum qdynp qdynh qfep qpi qpip qmin qminp: nrgy.o

qcalc: qatom.o

qcalc qprep qfep: parse.o

qcalc qprep: maskmanip.o

qdyn qdum: qdyn.o

qmin: qmin.o

qdyn qdum qdynp qdynh qpi qpip qmin qminp: qatom.o exc.o globals.o

qprep: q_prep.o prefs.o prep.o avetr.o nrgy.o bonded.o
	${FC} ${FFLAGS} ${FLIBS} $+ -o $@

qfep: qfep.o 
	${FC} ${FFLAGS} ${FLIBS} $+ -o $@

qcalc: calc_base.o calc_chemscore.o calc_fit.o calc_geom.o calc_pmfscore.o \
	calc_com_ke.o calc_com.o calc_rdf.o calc_rms.o calc_rmsf.o \
	calc_entropy.o calc_nb.o calc_xscore.o qcalc.o eigen.o globals.o
	${FC} ${FFLAGS} ${FLIBS} $+ -o $@

qdyn qdum qpi qmin: simprep.o qcp.o potene.o nonbondene.o nonbonded.o qalloc.o \
	bonded.o bondene.o mpiglob.o 

qdynp qdynh qpip qminp: simprep_mpi.o qcp_mpi.o potene_mpi.o nonbondene_mpi.o qalloc_mpi.o \
	nonbonded.o bonded.o bondene_mpi.o mpiglob_mpi.o

qdyn : md.o 
	${FC} ${FFLAGS} ${FLIBS} $+ -o $@

qmin : minim.o
	${FC} ${FFLAGS} ${FLIBS} $+ -o $@

qminp :qmin_mpi.o minim_mpi.o
	${MPIFC} ${MPIFLAGS} ${FFLAGS} ${FLIBS} $+ -o $@

qdum : md_dum.o
	${FC} ${FFLAGS} ${FLIBS} $+ -o $@

qdynp : md_mpi.o qdyn_mpi.o 
	${MPIFC} ${MPIFLAGS} ${FFLAGS} ${FLIBS} $+ -o $@

qdynh : md_mpi.o qdyn_mpi.o
	${MPIFC} ${MPIFLAGS} ${FFLAGS} ${FLIBS} $+ -o $@

qpi   : qpi.o
	${FC} ${FFLAGS} ${FLIBS} $+ -o $@

qpip  : qpi_mpi.o
	${MPIFC} ${MPIFLAGS} ${FFLAGS} ${FLIBS} $+ -o $@

################################################################################
# Object modules
################################################################################
avetr.o: avetr.f90 prep.o
	${FC} ${FFLAGS}  -c $<

bonded.o: bonded.f90 math.o sizes.o
	${FC} ${FFLAGS} -c $<

bondene.o: bondene.f90 math.o sizes.o bonded.o globals.o topo.o qatom.o qalloc.o
	${FC} ${FFLAGS} ${FPP_FLG} -c $<

bondene_mpi.o: bondene.f90 math.o sizes.o bonded.o globals.o topo.o qatom.o qalloc_mpi.o
	${MPIFC} ${MPIFLAGS} ${FFLAGS} ${FPP_FLG}  -c bondene.f90 -o bondene_mpi.o

calc_base.o:calc_base.f90 topo.o prmfile.o qatom.o
	${FC} ${FFLAGS} -c $<

calc_chemscore.o: calc_chemscore.f90 maskmanip.o trj.o prmfile.o index.o qatom.o
	${FC} ${FFLAGS} -c $<

calc_entropy.o:calc_entropy.f90 calc_base.o maskmanip.o trj.o calc_fit.o
	${FC} ${FFLAGS} -c $<

calc_fit.o:calc_fit.f90 calc_base.o maskmanip.o eigen.o
	${FC} ${FFLAGS} -c $<

calc_geom.o:calc_geom.f90 calc_base.o
	${FC} ${FFLAGS} -c $<

calc_nb.o:calc_nb.f90 calc_base.o maskmanip.o parse.o qatom.o prmfile.o
	${FC} ${FFLAGS} -c $<

calc_pmfscore.o: calc_pmfscore.f90 calc_base.o maskmanip.o trj.o topo.o prmfile.o index.o qatom.o misc.o
	${FC} ${FFLAGS} -c $<

calc_xscore.o: calc_xscore.f90 calc_base.o maskmanip.o trj.o topo.o prmfile.o index.o qatom.o misc.o
	${FC} ${FFLAGS} -c $<

calc_rdf.o:calc_rdf.f90 calc_base.o parse.o maskmanip.o
	${FC} ${FFLAGS} -c $<

calc_rms.o:calc_rms.f90 calc_base.o maskmanip.o
	${FC} ${FFLAGS} -c $<

calc_com_ke.o:calc_com_ke.f90 calc_base.o maskmanip.o
	${FC} ${FFLAGS} -c $<

calc_com.o:calc_com.f90 calc_base.o maskmanip.o
	${FC} ${FFLAGS} -c $<

calc_rmsf.o:calc_rmsf.f90 calc_base.o maskmanip.o
	${FC} ${FFLAGS} -c $<

eigen.o:eigen.f90
	${FC} ${FFLAGS} -c $<

exc.o:exc.f90 misc.o qatom.o maskmanip.o
	${FC} ${FFLAGS} -c $<

globals.o: globals.f90 mpiglob.o maskmanip.o math.o sizes.o
	${FC} ${FFLAGS} ${FPP_FLG} -c $<

index.o:index.f90
	${FC} ${FFLAGS} -c $<

mask.o:	mask.f90 topo.o
	${FC} ${FFLAGS} -c $<

maskmanip.o:maskmanip.f90 mask.o misc.o parse.o
	${FC} ${FFLAGS} -c $<

math.o: math.f90 sizes.o
	${FC} ${FFLAGS} ${FPP_FLG} -c $<

#arm grid by including this flag: -DUSE_GRID
md.o:	md.f90 mpiglob.o qatom.o sizes.o trj.o topo.o exc.o version.o \
	math.o globals.o simprep.o potene.o qcp.o
	${FC} ${FFLAGS} ${FPP_FLG} -c md.f90

md_dum.o:md.f90 mpiglob.o qatom.o sizes.o topo.o exc.o version.o \
	math.o globals.o simprep.o potene.o qcp.o
	${FC} ${FFLAGS} ${FPP_FLG} -DDUM -c md.f90 -o md_dum.o

md_mpi.o: md.f90 mpiglob_mpi.o qatom.o sizes.o topo.o trj.o exc.o version.o \
	math.o globals.o simprep_mpi.o potene_mpi.o qcp_mpi.o
	${MPIFC} ${MPIFLAGS} ${FFLAGS} ${FPP_FLG} -DUSE_MPI \
	-c md.f90 -o md_mpi.o

minim.o:  minim.f90 mpiglob.o qatom.o sizes.o trj.o topo.o exc.o version.o \
          math.o globals.o simprep.o potene.o qcp.o
	${FC} ${FFLAGS} ${FPP_FLG} -c minim.f90

minim_mpi.o: minim.f90 mpiglob_mpi.o qatom.o sizes.o topo.o trj.o exc.o version.o \
	math.o globals.o simprep_mpi.o potene_mpi.o qcp_mpi.o
	${MPIFC} ${MPIFLAGS} ${FFLAGS} ${FPP_FLG} -DUSE_MPI \
	-c minim.f90 -o minim_mpi.o

misc.o: misc.f90 sizes.o
	${FC} ${FFLAGS} -c $<

mpiglob.o: mpiglob.f90 sizes.o nrgy.o
	${FC} ${FFLAGS} ${FPP_FLG} -c $<

mpiglob_mpi.o: mpiglob.f90 sizes.o nrgy.o
	${MPIFC} ${MPIFLAGS} ${FPP_FLG} ${FFLAGS} -DUSE_MPI \
	-c mpiglob.f90 -o mpiglob_mpi.o

nonbonded.o : nonbonded.f90 sizes.o globals.o math.o qatom.o
	${FC} ${FFLAGS} -c $<

nonbondene.o: nonbondene.f90 sizes.o qatom.o qalloc.o globals.o math.o nonbonded.o exc.o mpiglob.o
	${FC} ${FFLAGS} ${FPP_FLG} -c nonbondene.f90

nonbondene_mpi.o: nonbondene.f90 sizes.o qatom.o qalloc_mpi.o globals.o math.o nonbonded.o exc.o mpiglob_mpi.o
	${MPIFC} ${MPIFLAGS} ${FFLAGS} ${FPP_FLG} -DUSE_MPI \
	-c nonbondene.f90 -o nonbondene_mpi.o

nrgy.o: nrgy.f90 sizes.o
	${FC} ${FFLAGS} -c $<

parse.o: parse.f90 misc.o
	${FC} ${FFLAGS} -c $<

potene.o: potene.f90 sizes.o globals.o nonbondene.o bondene.o mpiglob.o
	${FC} ${FFLAGS} ${FPP_FLG} -c potene.f90

potene_mpi.o: potene.f90 sizes.o globals.o bondene_mpi.o nonbondene_mpi.o mpiglob_mpi.o
	${MPIFC} ${MPIFLAGS} ${FFLAGS} ${FPP_FLG} -DUSE_MPI \
	-c potene.f90 -o potene_mpi.o

prefs.o: prefs.f90
	${FC} ${FFLAGS} -c $<

prep.o: prep.f90 maskmanip.o sizes.o parse.o prmfile.o trj.o index.o prefs.o
	${FC} ${FFLAGS} ${FPP_FLG} -c $<

prmfile.o: prmfile.f90 misc.o mpiglob.o
	${FC} ${FFLAGS} -c $<

qalloc.o: qalloc.f90 sizes.o mpiglob.o globals.o qatom.o math.o trj.o
	${FC} ${FFLAGS} ${FPP_FLG} -c qalloc.f90

qalloc_mpi.o: qalloc.f90 sizes.o mpiglob_mpi.o globals.o qatom.o math.o trj.o
	${MPIFC} ${MPIFLAGS} ${FFLAGS} ${FPP_FLG} -DUSE_MPI \
	-c qalloc.f90 -o qalloc_mpi.o

q_prep.o: q_prep.f90 prep.o avetr.o version.o
	${FC} ${FFLAGS} -c $<

qatom.o: qatom.f90 misc.o nrgy.o prmfile.o sizes.o index.o topo.o globals.o
	${FC} ${FFLAGS} -c $<

qcalc.o: qcalc.f90 calc_chemscore.o calc_pmfscore.o calc_xscore.o trj.o calc_base.o calc_rms.o calc_fit.o calc_geom.o version.o
	${FC} ${FFLAGS} ${FPP_FLG} -c qcalc.f90

qcp.o: 	qcp.f90 topo.o qatom.o sizes.o nrgy.o potene.o globals.o math.o mpiglob.o exc.o
	${FC} ${FFLAGS} ${FPP_FLG} -c qcp.f90

qcp_mpi.o: qcp.f90 topo.o qatom.o sizes.o nrgy.o potene_mpi.o globals.o math.o mpiglob_mpi.o exc.o
	${MPIFC} ${MPIFLAGS} ${FFLAGS} ${FPP_FLG} -DUSE_MPI \
	-c qcp.f90 -o qcp_mpi.o

qdyn.o: qdyn.f90 md.o mpiglob.o version.o simprep.o 
	${FC} ${FFLAGS} ${FPP_FLG} ${BUILD_NUMBER_LDFLAGS} -c qdyn.f90

qdyn_dum.o: qdyn.f90 md.o mpiglob.o version.o simprep.o
	${FC} ${FFLAGS} ${FPP_FLG}  ${BUILD_NUMBER_LDFLAGS} -DDUM -c qdyn.f90 -o qdyn_dum.o

qdyn_mpi.o: qdyn.f90 md_mpi.o mpiglob_mpi.o version.o simprep_mpi.o
	${MPIFC} ${MPIFLAGS} ${FFLAGS} ${FPP_FLG}  ${BUILD_NUMBER_LDFLAGS} -DUSE_MPI \
        -c qdyn.f90 -o qdyn_mpi.o

qfep.o: qfep.f90 nrgy.o parse.o version.o
	${FC} ${FFLAGS} ${FPP_FLG} -c qfep.f90

qmin.o: qmin.f90 minim.o mpiglob.o version.o simprep.o
	${FC} ${FFLAGS} ${FPP_FLG} ${BUILD_NUMBER_LDFLAGS} -c qmin.f90

qmin_mpi.o: qmin.f90 minim_mpi.o simprep_mpi.o mpiglob_mpi.o
	${MPIFC} ${MPIFLAGS} ${FFLAGS} ${FPP_FLG}  ${BUILD_NUMBER_LDFLAGS} -DUSE_MPI \
	-c qmin.f90 -o qmin_mpi.o

qpi.o: qpi.f90 qcp.o version.o simprep.o mpiglob.o
	${FC} ${FFLAGS} ${FPP_FLG} ${BUILD_NUMBER_LDFLAGS} -c qpi.f90

qpi_mpi.o: qpi.f90 qcp_mpi.o version.o simprep_mpi.o mpiglob_mpi.o
	${MPIFC} ${MPIFLAGS} ${FFLAGS} ${FPP_FLG}  ${BUILD_NUMBER_LDFLAGS} -DUSE_MPI \
	-c qpi.f90 -o qpi_mpi.o

simprep.o: simprep.f90 sizes.o trj.o mpiglob.o qatom.o version.o globals.o exc.o math.o qalloc.o qcp.o \
	nonbondene.o mpiglob.o
	${FC} ${FFLAGS} ${FPP_FLG} -c simprep.f90

simprep_mpi.o: simprep.f90 sizes.o trj.o mpiglob.o qatom.o version.o globals.o exc.o math.o qalloc_mpi.o qcp_mpi.o \
	nonbondene_mpi.o mpiglob_mpi.o
	${MPIFC} ${MPIFLAGS} ${FFLAGS} ${FPP_FLG} -DUSE_MPI \
	-c simprep.f90 -o simprep_mpi.o

sizes.o: sizes.f90
	${FC} ${FFLAGS} ${FPP_FLG} -c $<

topo.o:  topo.f90 misc.o mpiglob.o sizes.o version.o math.o
	${FC} ${FFLAGS} -c $<

trj.o:  trj.f90 mask.o misc.o
	${FC} ${FFLAGS} -c $<

version.o: version.f90
	${FC} ${FFLAGS} ${FPP_FLG} ${BUILD_NUMBER_LDFLAGS} -c $<

moveqdyn:
	mkdir -p ../bin ; mv qdyn ../bin/qdyn
move1:
	mkdir -p ../obj ; mv *.o *.mod ../obj/
	@if ! test -f $(BUILD_NUMBER_FILE); then echo 0 > $(BUILD_NUMBER_FILE); fi
	@echo $$(($$(cat $(BUILD_NUMBER_FILE)) + 1)) > $(BUILD_NUMBER_FILE)

move2:
	mkdir -p ../bin ; mv qfep qprep qdum qcalc qpi qmin ../bin/

move3:
	mkdir -p ../bin ; mv qdynp qpip qminp ../bin/

move4:
	mkdir -p ../bin ; mv qdynh ../bin/




include buildnumber.mak

################################################################################


################################################################################
# Build instructions
################################################################################
what:
	@echo "Use: make [target] [COMP=compiler]"
	@echo
	@echo "[target] is one of:"
	@echo "all              Everything, including optimized and optimized qdyn/qdynp"
	@echo "tools            Speed optimized versions of qdyn, qdum, qfep, qcalc and qprep"
	@echo "debug            Tools with debug information (stacktraces)."
	@echo "mpi              qdynp using a currently loaded MPI library"
	@echo "mpidebug         qdynp with debug information (stacktraces)."
	@echo "hybrid           qdynh using a currently loaded MPI library and OpenMP library"
	@echo "hybriddebug      qdynh with debug information (stacktraces)."
	@echo
	@echo "For compiler/os pick one of:"
	@echo "ifort            Intel Fortran compiler"
	@echo "osx              Gnu Compiler Collection (GCC) in mac"
	@echo "gcc              GCC in linux"
	@echo "pgi              Portland Group compiler"
################################################################################

