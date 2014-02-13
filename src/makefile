################################################################################
#  Q V.5 Makefile for linux                                                    #
#  main code author: Johan Aqvist                                              #
#  maintainers: Beat Amrein, Alexandre Barrozo, Paul Bauer, Mauricio Esguerra  #
#  latest update: January 29, 2014                                             #
################################################################################

# Specify compilers and their options in the following variables.
FC=		gfortran
CC=		gcc
MPIFC=		mpif90
FC_OPT=		-O3 -unroll --static -ffree-line-length-none -fcray-pointer -fall-intrinsics -std=legacy -Wall
CC_OPT=		-O3
#MPIFLAGS=	-f90=${FC}
OMPI_FC=	${FC}
OMPI_CC=	${CC}	
FFLAGS=		-fstack-protector ${FC_OPT} -DG95=1
CCFLAGS=	${CC_OPT}
OMPI_FCFLAGS=	${FFLAGS}
OMPI_CXXFLAGS=	${CCFLAGS}
#FLIBS=		-Vaxlib
FPP_FLG=	-cpp  
# 
###

###
# Display options
default: what
###

###
# Targets
all:	qfep5 qprep5 qdyn5 qdum5 qcalc5 move1 move2

debug:	
	@make FC_OPT='-O0 -unroll --static -ffree-line-length-none -fcray-pointer -fall-intrinsics -std=legacy' FFLAGS='-fstack-protector ${FC_OPT} -DG95=1'         \
	qfep5 qprep5 qdyn5 qdum5 qcalc5

mpi:	
	@make FC_OPT=-'O3 -unroll -ffree-line-length-none -fcray-pointer -fall-intrinsics -std=legacy' FFLAGS='-fstack-protector ${FC_OPT} -DG95=1' \
	qdyn5p

mpi_debug:	
	@make FC_OPT=-'O0 -unroll -ffree-line-length-none -fcray-pointer -fall-intrinsics -std=legacy -g'          \
	qdyn5p

clean:
	-rm -f *.o  *.mod *.M *.kmo *.il

nuke:
	-rm -f *.o  *.mod *.M *.kmo *.il qfep5 qdyn5p qdyn5 qprep5 qcalc5 qdum5

qcalc5 qdyn5 qdum5 qdyn5p qprep5 qfep5: misc.o mpiglob.o 

qcalc5 qdyn5 qdum5 qdyn5p qprep5: mask.o prmfile.o sizes.o topo.o trj.o index.o 

qcalc5 qdyn5 qdum5 qdyn5p qfep5: nrgy.o

qcalc5 qdyn5 qdum5 qdyn5p: qatom.o  

qcalc5 qprep5 qfep5: parse.o 

qcalc5 qprep5: maskmanip.o 

qdyn5 qdum5: qdyn.o

qprep5: q_prep.o prefs.o prep.o avetr.o
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

###
# Object modules
avetr.o: avetr.f90 prep.o
	${FC} ${FFLAGS} -c $<

calc_base.o:calc_base.f90 topo.o
	${FC} ${FFLAGS} -c $<

calc_chemscore.o: calc_chemscore.f90 maskmanip.o trj.o prmfile.o index.o qatom.o
	${FC} ${FFLAGS} -c $<

calc_entropy.o:calc_entropy.f90 calc_base.o maskmanip.o trj.o calc_fit.o
	${FC} ${FFLAGS} -c $<

calc_fit.o:calc_fit.f90 calc_base.o maskmanip.o
	${FC} ${FFLAGS} -c $<

calc_geom.o:calc_geom.f90 calc_base.o
	${FC} ${FFLAGS} -c $<

calc_nb.o:calc_nb.f90 calc_base.o maskmanip.o parse.o
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

invsqrt_q.o:invsqrt_q.c
	$(CC) $(CCFLAGS) -c $<

index.o:index.f90
	${FC} ${FFLAGS} -c $<

mask.o:	mask.f90 topo.o
	${FC} ${FFLAGS} -c $<

maskmanip.o:maskmanip.f90 mask.o misc.o parse.o
	${FC} ${FFLAGS} -c $<

md.o:	md.f90 mpiglob.o qatom.o sizes.o trj.o topo.o 
	${FC} ${FFLAGS} ${FPP_FLG} -c md.f90

md_dum.o:md.f90 mpiglob.o qatom.o sizes.o topo.o
	${FC} ${FFLAGS} ${FPP_FLG} -DDUM -c md.f90 -o md_dum.o 

md_mpi.o: md.f90 mpiglob.o qatom.o sizes.o topo.o trj.o
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
	mkdir ../obj ; mv *.o *.mod ../obj/

move2:
	mkdir ../bin ; mv qfep5 qprep5 qdyn5 qdyn5p qdum5 qcalc5 ../bin/
#
###

###
# Build instructions
what:
	@echo "usage: make <target>"
	@echo
	@echo "<target> is one of:"
	@echo "all              Everything except the MPI parallel version of Q (qdyn5p)"
	@echo "debug            'all' but with debug information (stacktraces etc)."
	@echo "mpi              qdyn5p using a loaded mpi library"
	@echo "mpi_debug        'mpi' but with debug information (stacktraces etc)."
	@echo
#
###
