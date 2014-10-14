################################################################################
#  Q V.5.06 Makefile                                                           #
#  code authors: Johan Aqvist, Martin Almlof, Martin Ander, Jens Carlson,      #
#  Isabella Feierberg, Peter Hanspers, Anders Kaplan, Karin Kolmodin,          #
#  Kajsa Ljunjberg, John Marelius, Martin Nervall                              #
#  keeper of the code: Paul Bauer                                              #
#  maintainers: Beat Amrein, Alexandre Barrozo, Paul Bauer, Mauricio Esguerra  #
#  latest update: May 9, 2014                                                  #
#  make [option] [COMP=compiler]                                               #
################################################################################

################################################################################
# Intel Fortran ifort ( tested in linux and osx)
################################################################################
ifeq ($(COMP),ifort)
	FC=		ifort
	FC_OPT=         -O3 -ip -ipo -unroll -i-static
	FC_DEBUG=	-unroll 
	FFLAGS=         -fstack-protector ${FC_OPT}
	FPP_FLG=	-fpp

	MPIFC=		mpif90
	MPIFLAGS=	-f90=${FC} -Nmpi
	OMPI_FC=	${FC}
	MPI_FC=		${FC}
	MPI_FCFLAGS=	${FFLAGS}
	MPI=		${MPI_FCFLAGS}

	DEBUG=          -debug -traceback -g  -I${VT_ROOT}/include -L${VT_LIB_DIR} ${VT_ADD_LIBS} ${FC_OPT}
	DEBUGMPI=       -debug -traceback -g  -I${VT_ROOT}/include -L${VT_LIB_DIR} ${VT_ADD_LIBS} ${FC_OPT}

################################################################################
# OSX Mavericks using gfortran
################################################################################
else ifeq ($(COMP),osx)
	FC= gfortran
	FC_OPT=         -O3 -ffree-line-length-none -fcray-pointer -fall-intrinsics \
	                -std=legacy -Wall -Wtabs
	FFLAGS=         -fstack-protector ${FC_OPT} -DG95=1
	FPP_FLG=        -cpp

	MPIFC=          mpif90
	OMPI_FC=        ${FC}
	OMPI_FCFLAGS=   ${FFLAGS}
	MPI=		${OMPI_FCFLAGS}

	DEBUG=          -debug  -g ${FC_OPT}
	DEBUGMPI=       -debug  -g ${FC_OPT}

################################################################################
# PGI Compiler -- compiles at csb
################################################################################
else ifeq ($(COMP),pgi)
	FC=		pgf90
	FC_OPT=		-fastsse -Mipa=fast,inline -O3 -Munroll -mcmodel=medium
	FFLAGS=		${FC_OPT}
	FPP_FLG=	-Mpreprocess

	MPIFC=		mpif90
	MPIFLAGS=	-lmpi#-f90=${FC}
	MPI=		${MPIFLAGS} ${FFLAGS}
	MPI_FC=		${FC}

################################################################################
# IBM Compiler -- untested
################################################################################
else ifeq ($(COMP),ibm)
	FC=             xlf90
	FC_OPT=         -qsmp=omp
	FFLAGS=

################################################################################
# GNU gcc-gfortran
################################################################################
else ifeq ($(COMP),gcc)
	FC=		gfortran
	FC_OPT=		-O3 -unroll -ffree-line-length-none -fcray-pointer -fall-intrinsics -std=legacy
	FFLAGS=		-fstack-protector ${FC_OPT} -DG95=1
	FPP_FLG=	-cpp

	MPIFC=		mpif90
	OMPI_FC=	${FC}
	OMPI_FCFLAGS=	${FFLAGS}
	MPI=		${OMPI_FCFLAGS}

	DEBUG=          -fbacktrace -g ${DEBUG_FC}	
	DEBUGMPI=       -fbacktrace -g ${DEBUG_FC}

else
	FC=		gfortran
	FC_OPT=		-O3 -unroll -ffree-line-length-none -fcray-pointer -fall-intrinsics -std=legacy
	FFLAGS=		-fstack-protector ${FC_OPT} -DG95=1
	FPP_FLG=	-cpp

	MPI=		${FFLAGS}
endif

################################################################################
# Display options
################################################################################
default: what

################################################################################
# Targets
################################################################################
all:	qfep5 qprep5 qdyn5 qdum5 qcalc5 move1 move2

debug:
	@make FFLAGS="${DEBUG}" \
	qfep5 qprep5 qdyn5 qdum5 qcalc5 move1 move2

mpi:
	@make FFLAGS="${MPI}" \
	qdyn5p move1 move3

mpidebug:
	@make FFLAGS="${DEBUGMPI}" \
	qdyn5p move1 move3

clean:
	-rm -f *.o  *.mod *.M *.kmo *.il *.oo

nuke:
	-rm -rf *.o  *.mod *.M *.kmo *.il *.oo qfep5 qdyn5p qdyn5 qprep5 qcalc5 qdum5 ../bin ../obj

qcalc5 qdyn5 qdum5 qdyn5p qprep5 qfep5: misc.o mpiglob.o

qcalc5 qdyn5 qdum5 qdyn5p qprep5: mask.o prmfile.o sizes.o topo.o trj.o index.o

qcalc5 qdyn5 qdum5 qdyn5p qfep5: nrgy.o

qcalc5 qdyn5 qdum5 qdyn5p: qatom.o exc.o

qcalc5 qprep5 qfep5: parse.o

qcalc5 qprep5: maskmanip.o

qdyn5 qdum5: qdyn.o

qprep5: q_prep.o prefs.o prep.o avetr.o nrgy.o
	${FC} ${FFLAGS} ${FLIBS} $+ -o $@

qfep5: qfep.o
	${FC} ${FFLAGS} ${FLIBS} $+ -o $@

qcalc5: calc_base.o calc_chemscore.o calc_fit.o calc_geom.o calc_pmfscore.o \
	calc_com_ke.o calc_com.o calc_rdf.o calc_rms.o calc_rmsf.o \
	calc_entropy.o calc_nb.o calc_xscore.o eigen.o qcalc.o
	${FC} ${FFLAGS} ${FLIBS} $+ -o $@

qdyn5 : md.o
	${FC} ${FFLAGS} ${FLIBS} $+ -o $@

qdum5 : md_dum.o
	${FC} ${FFLAGS} ${FLIBS} $+ -o $@

qdyn5p : md_mpi.o qdyn_mpi.o
	${MPIFC} ${MPIFLAGS} ${FFLAGS} ${FLIBS} $+ -o $@

################################################################################
# Object modules
################################################################################
avetr.o: avetr.f90 prep.o
	${FC} ${FFLAGS} -c $<

calc_base.o:calc_base.f90 topo.o prmfile.o qatom.o
	${FC} ${FFLAGS} -c $<

calc_chemscore.o: calc_chemscore.f90 maskmanip.o trj.o prmfile.o index.o qatom.o
	${FC} ${FFLAGS} -c $<

calc_entropy.o:calc_entropy.f90 calc_base.o maskmanip.o trj.o calc_fit.o
	${FC} ${FFLAGS} -c $<

calc_fit.o:calc_fit.f90 calc_base.o maskmanip.o
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

index.o:index.f90
	${FC} ${FFLAGS} -c $<

mask.o:	mask.f90 topo.o
	${FC} ${FFLAGS} -c $<

maskmanip.o:maskmanip.f90 mask.o misc.o parse.o
	${FC} ${FFLAGS} -c $<

md.o:	md.f90 mpiglob.o qatom.o sizes.o trj.o topo.o exc.o
	${FC} ${FFLAGS} ${FPP_FLG} -c md.f90

md_dum.o:md.f90 mpiglob.o qatom.o sizes.o topo.o
	${FC} ${FFLAGS} ${FPP_FLG} -DDUM -c md.f90 exc.o -o md_dum.o

md_mpi.o: md.f90 mpiglob.o qatom.o sizes.o topo.o trj.o exc.o
	${MPIFC} ${MPIFLAGS} ${FFLAGS} ${FPP_FLG} -DUSE_MPI \
	-c md.f90 -o md_mpi.o

misc.o: misc.f90 sizes.o
	${FC} ${FFLAGS} -c $<

mpiglob.o: mpiglob.f90 sizes.o nrgy.o
	${FC} ${FFLAGS} -c $<

nrgy.o: nrgy.f90 sizes.o
	${FC} ${FFLAGS} -c $<

parse.o: parse.f90 misc.o
	${FC} ${FFLAGS} -c $<

prefs.o: prefs.f90
	${FC} ${FFLAGS} -c $<

prep.o: prep.f90 maskmanip.o sizes.o parse.o prmfile.o trj.o index.o prefs.o
	${FC} ${FFLAGS} -c $<

prmfile.o: prmfile.f90 misc.o mpiglob.o
	${FC} ${FFLAGS} -c $<

q_prep.o: q_prep.f90 prep.o avetr.o
	${FC} ${FFLAGS} -c $<

qatom.o: qatom.f90 misc.o nrgy.o prmfile.o sizes.o index.o topo.o
	${FC} ${FFLAGS} -c $<

qcalc.o: qcalc.f90 calc_chemscore.o calc_pmfscore.o calc_xscore.o trj.o calc_base.o calc_rms.o calc_fit.o calc_geom.o
	${FC} ${FFLAGS} ${FPP_FLG} -c qcalc.f90

qdyn.o: qdyn.f90 md.o mpiglob.o
	${FC} ${FFLAGS} ${FPP_FLG} -c qdyn.f90

qdyn_dum.o: qdyn.f90 md.o mpiglob.o
	${FC} ${FFLAGS} ${FPP_FLG} -DDUM -c qdyn.f90 -o qdyn_dum.o

qdyn_mpi.o: qdyn.f90 md_mpi.o mpiglob.o
	${MPIFC} ${MPIFLAGS} ${FFLAGS} ${FPP_FLG} -DUSE_MPI \
        -c qdyn.f90 -o qdyn_mpi.o

qfep.o: qfep.f90 nrgy.o parse.o
	${FC} ${FFLAGS} ${FPP_FLG} -c qfep.f90

sizes.o: sizes.f90
	${FC} ${FFLAGS} -c $<

topo.o:  topo.f90 misc.o mpiglob.o sizes.o
	${FC} ${FFLAGS} -c $<

trj.o:  trj.f90 mask.o misc.o
	${FC} ${FFLAGS} -c $<

move1:
	mkdir -p ../obj ; mv *.o *.mod ../obj/

move2:
	mkdir -p ../bin ; mv qfep5 qprep5 qdyn5 qdum5 qcalc5 ../bin/

move3:
	mkdir -p ../bin ; mv qdyn5p ../bin/

################################################################################


################################################################################
# Build instructions
################################################################################
what:
	@echo "Use: make [target] [COMP=compiler]"
	@echo
	@echo "[target] is one of:"
	@echo "all              Everything except the MPI parallel version of Q (qdyn5p)"
	@echo "debug            All with debug information (stacktraces)."
	@echo "mpi              qdyn5p using a currently loaded MPI library"
	@echo "mpidebug         qdyn5p with debug information (stacktraces)."
	@echo
	@echo "For compiler/os pick one of:"
	@echo "ifort            Intel Fortran compiler"
	@echo "osx              Gnu Compiler Collection (GCC) in mac"
	@echo "gcc              GCC in linux"
	@echo "pgi              Portland Group compiler"
################################################################################

