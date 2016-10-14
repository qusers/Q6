! (C) 2014 Uppsala Molekylmekaniska HB, Uppsala, Sweden
! simprep.f90
! based on
! md.f90
! by Johan Åqvist, John Marelius, Anders Kaplan, Isabella Feierberg, Martin Nervall & Martin Almlöf
! general simulation setup and memory management

module MD

! used modules
!use PROFILING
use SIZES
use TRJ
use MPIGLOB
use QATOM
use EXC
use QCP
use VERSIONS
#if defined (_DF_VERSION_)
use DFPORT
use DFLIB
#endif
!$ use omp_lib
implicit none
#if defined (USE_MPI)
include "mpif.h"
#endif



!-----------------------------------------------------------------------
!	shared variables
!-----------------------------------------------------------------------
!	Constants
real(kind=prec)			::	pi, deg2rad	!set in sub startup
! no hard coded values for something like this
! rho will now be read from the topology, as it should be
real(kind=prec) 		::  rho_wat = -1.0_prec  ! molecules / A**3, set from topo
real(kind=prec), parameter		::  Boltz = 0.001986_prec !TODO: it is not place for that

!Read status
integer                     :: stat

!print temperature if it changed more than 2% in one time step
real(kind=prec), parameter				::	TEMP_PRINT_THRESHOLD=0.02_prec !TODO: it is not place for that


!	Memory management
integer						::	alloc_status


!	Topology information
! --- Atoms
integer						::	natom
integer                     :: nat3
! --- q-atom number or 0 for non-q
integer(TINY), allocatable		::	iqatom(:)
! --- water topology
! atoms of water, total water molecules, excluded water molecules
integer						::	nwat
real(kind=prec)					:: mu_w

!-------------------------------------------------------------------
!	Periodic box information
!-------------------------------------------------------------------

!******PWadded variable 2001-10-10
logical						:: box, rigid_box_centre
logical						:: put_solute_back_in_box, put_solvent_back_in_box

!variabels used for constant pressure algorithm
logical						:: constant_pressure = .false.
integer						:: volume_try, volume_acc
logical                                         :: atom_based_scaling = .false.
real(kind=prec)						:: pressure, max_vol_displ
integer						:: pressure_seed
real(kind=prec), allocatable		:: mol_mass(:)
real(kind=prec), allocatable		:: mass(:)



!variabels used when user change boxsize in inputfile
logical						:: control_box
real(kind=prec)						:: new_boxl(3) 

!-----------------------------------------------------------------------
!	Dynamics control information     
!-----------------------------------------------------------------------


! --- MD parameters
! friction, gkT, randf and randv are variables to be used for the Langevin thermostat - A. Barrozo
! numchain and nhchain are variables to be used for the Nosé-Hoover thermostat
integer						::	nsteps, istep
integer						::	iseed, numchain
logical						::	restart
real(kind=prec)						::	dt, dt2, Temp0, Tmaxw, tau_T, friction, gkT, randf, randv, kbT, nhq
logical						::	shake_solvent, shake_solute
logical						::	shake_hydrogens, shake_heavy
logical						::	separate_scaling
!new logical to select for gaussian or uniform random numbers in LANGEVIN
!themorstat
logical                                         ::      langevin_gauss_rand = .false.
!Paul does not like overflow errors, increased character array length
!Added ENUM for names to allow integer comparisson later on
character(len=80)				::	name_thermostat
ENUM, bind(c) 
	ENUMERATOR	:: BERENDSEN, LANGEVIN, NOSEHOOVER
END ENUM
integer						::	thermostat = -1
!Added support for different integration schemes
!For now leap-frog and velocity-verlet
!Alexandre Barrozo and Paul Bauer, October 2014
character(len=80)			:: name_integrator
integer						:: integrator = -1
ENUM, bind(c)
	ENUMERATOR	:: LEAPFROG,VELVERLET
END ENUM
real(kind=prec), allocatable				::	xnh(:), vnh(:), qnh(:), Gnh(:) 

! --- Non-bonded strategy
logical						::	use_LRF
integer						::	NBcycle
real(kind=prec)						::	Rcpp, Rcww, Rcpw, Rcq, RcLRF
! we do this now without limits
!integer, parameter			                ::	max_atyp = 255
integer                                                 ::      num_atyp = -1
integer,allocatable					::	ljcod(:,:)!max_atyp,max_atyp)


! --- Output control parameters
integer						::	itrj_cycle, iene_cycle
integer						::	itemp_cycle, iout_cycle
logical						::	force_rms
!******PWadded 2001-10-23
integer						::	ivolume_cycle
!PB added for box control			
integer						:: nsolvent = 0,nsolute = 0

integer,allocatable				:: mvd_mol(:)
! --- Protein boundary
logical						::	exclude_bonded
real(kind=prec)						::	fk_pshell
real(kind=prec)      :: fk_fix    = 200.0_prec

! --- Water sphere
real(kind=prec)						::	Dwmz, awmz, rwat_in
real(kind=prec)						::	fkwpol
logical						::	wpol_restr, wpol_Born
logical      :: freeze
real(kind=prec)						::	fk_wsphere, crgtot, crgQtot
integer(AI), allocatable	::	list_sh(:,:), nsort(:,:)
real(kind=prec),  allocatable		::	theta(:), theta0(:), tdum(:)
integer						::	nwpolr_shell, n_max_insh


type SHELL_TYPE
        real(kind=prec)					::	rout, dr, cstb
        real(kind=prec)					::	avtheta, avn_insh, theta_corr
        integer					::	n_insh
end type SHELL_TYPE

type OLD_SHELL_TYPE
        real                                 ::      rout, dr, cstb
        real                                 ::      avtheta, avn_insh, theta_corr
        integer                                 ::      n_insh
end type OLD_SHELL_TYPE

type SHELL_TYPE_SINGLE
        real(kind=singleprecision)                                 ::      rout, dr, cstb
        real(kind=singleprecision)                                 ::      avtheta, avn_insh, theta_corr
        integer                                 ::      n_insh
end type SHELL_TYPE_SINGLE

type SHELL_TYPE_DOUBLE
        real(kind=doubleprecision)                                 ::      rout, dr, cstb
        real(kind=doubleprecision)                                 ::      avtheta, avn_insh, theta_corr
        integer                                 ::      n_insh
end type SHELL_TYPE_DOUBLE
#ifndef PGI
type SHELL_TYPE_QUAD
        real(kind=quadprecision)                                 ::      rout, dr, cstb
        real(kind=quadprecision)                                 ::      avtheta, avn_insh, theta_corr
        integer                                 ::      n_insh
end type SHELL_TYPE_QUAD
#endif

type(SHELL_TYPE), allocatable::	wshell(:)
type(OLD_SHELL_TYPE), allocatable::old_wshell(:)
type(SHELL_TYPE_SINGLE), allocatable:: wshell_single(:)
type(SHELL_TYPE_DOUBLE), allocatable:: wshell_double(:)
#ifndef PGI
type(SHELL_TYPE_QUAD), allocatable:: wshell_quad(:)
#endif
! constants & default values
integer, parameter			::	itdis_update		= 100
real(kind=prec),  parameter			::	wpolr_layer			= 3.0001_prec !TODO: it is not place for that
real(kind=prec), parameter				::	drout				= 0.5_prec !TODO: it is not place for that
real(kind=prec), parameter			::	tau_T_default		= 100.0_prec
real(kind=prec), parameter			::	rcpp_default		= 10.0_prec
real(kind=prec), parameter			::	rcww_default		= 10.0_prec
real(kind=prec), parameter			::	rcpw_default		= 10.0_prec
real(kind=prec), parameter                      ::      rcq_default_pbc         = -1.0_prec
real(kind=prec), parameter                      ::      rcLRF_default_pbc       = -1.0_prec
real(kind=prec), parameter			::	rcq_default_sph		= 99.0_prec
real(kind=prec), parameter			::	rcLRF_default		= 99.0_prec
!recxl_i is set to rexcl_o * shell_default
real(kind=prec), parameter   :: fk_fix_default = 200.0_prec
real(kind=prec), parameter			::	shell_default		= 0.85_prec
real(kind=prec), parameter			::	fk_pshell_default	= 10.0_prec
integer, parameter			::	itrj_cycle_default	= 100
integer, parameter			::	iene_cycle_default	= 10
integer, parameter			::	iout_cycle_default	= 10
integer, parameter			::	nb_cycle_default	= 10
real(kind=prec), parameter			::	fkwpol_default		= 20.0_prec
real(kind=prec), parameter			::	fk_wsphere_default	= 60.0_prec
logical, parameter			::	wpol_restr_default	= .true.

integer, parameter			::	ivolume_cycle_default = 10

integer, parameter			::	ipressure_cycle_default = 100
 ! Yes, and use Born corr.
!constants in the sigmoid function giving default Dwmz as function of radius


! --- File names
character(len=200)			::	top_file
character(len=200)			::	restart_file
character(len=200)			::	xfin_file
character(len=200)			::	trj_file
character(len=200)			::	fep_file
character(len=200)			::	ene_file
character(len=200)			::	exrstr_file
character(len=200)			::	xwat_file

! --- Excluded groups, added 2014 Paul Bauer
integer					:: ngroups_gc
logical					:: use_excluded_groups
type(MASK_TYPE)                         :: allmask
type(Q_ENE_HEAD)			:: ene_header

! --- Restraints
integer						::	implicit_rstr_from_file
integer      :: nrstr_seq, nrstr_pos, nrstr_dist, nrstr_angl, nrstr_wall


type RSTRSEQ_TYPE
        integer(AI)				::	i,j
        real(kind=prec)					::	fk
        integer(TINY)			::	ih
        integer					::	to_centre !flag for restraining to geom. or mass centre
end type RSTRSEQ_TYPE


type RSTRPOS_TYPE
        integer(AI)				::	i
        integer(TINY)			::	ipsi
        real(kind=prec)					::	fk(3)
        real(kind=prec)					::	x(3)
end type RSTRPOS_TYPE


type RSTRDIS_TYPE
        integer(AI)				::	i,j
        integer(TINY)			::	ipsi
        real(kind=prec)					::	fk
        real(kind=prec)					::	d1, d2
        character(len=20)       ::  itext,jtext
end type RSTRDIS_TYPE

type RSTRANG_TYPE 
        integer(AI)    :: i,j,k
        integer(TINY)  :: ipsi
        real(kind=prec)        :: fk
        real(kind=prec)        :: ang
!       character(len=20)       ::  itext,jtext,ktext
end type RSTRANG_TYPE                                                                                                

type RSTRWAL_TYPE

        integer(AI)				::	i,j
        real(kind=prec)					::	d, fk, aMorse, dMorse
        integer(TINY)			::	ih
end type RSTRWAL_TYPE


type(RSTRSEQ_TYPE), allocatable::	rstseq(:)
type(RSTRPOS_TYPE), allocatable::	rstpos(:)
type(RSTRDIS_TYPE), allocatable::	rstdis(:)
type(RSTRANG_TYPE), allocatable:: rstang(:)
type(RSTRWAL_TYPE), allocatable::	rstwal(:)

! New ! Type for internal interactions of solvent atoms
! will be precomputed in prep sim
type SOLV_INT_TYPE
	integer(AI)         			:: i,j
	real(kind=prec)				:: vdWA,vdWB
	real(kind=prec)				:: elec
end type SOLV_INT_TYPE

type(SOLV_INT_TYPE),allocatable			:: nonbnd_solv_int(:)
integer						:: num_solv_int = -1

! New ! Type for precomputing interactions between arbitrary atoms pairs
type PRECOMPUTE_INT
        real(kind=prec)                         :: vdWA,vdWB,elec,score
        logical                                 :: set
end type PRECOMPUTE_INT
! special for qq, because soft pairs
type PRECOMPUTE_INT_QQ
        real(kind=prec)                         :: vdWA,vdWB,elec,score
        logical                                 :: set,soft
end type PRECOMPUTE_INT_QQ


type(PRECOMPUTE_INT),allocatable                :: pp_precomp(:,:)
type(PRECOMPUTE_INT),allocatable                :: pw_precomp(:,:)
type(PRECOMPUTE_INT),allocatable                :: ww_precomp(:,:)
type(PRECOMPUTE_INT),allocatable                :: qp_precomp(:,:,:)
type(PRECOMPUTE_INT),allocatable                :: qw_precomp(:,:,:)
type(PRECOMPUTE_INT_QQ),allocatable             :: qq_precomp(:,:,:)
!-----------------------------------------------------------------------
!	Coordinates, velocities, forces
!-----------------------------------------------------------------------
real(kind=prec), allocatable		::	d(:)
real(kind=prec), allocatable		::	x(:)
real(kind=prec), allocatable		::	xx(:) !for shake
real(kind=prec), allocatable		::	v(:)
real(kind=prec), allocatable			::	winv(:)
real(kind=prec)						::	grms !RMS force


!-----------------------------------------------------------------------
!	Energies , EQ is defined in qatom.f90
!-----------------------------------------------------------------------
type(ENERGIES)				::	E,old_E,previous_E
!added the variables for the MC_volume function as globals
real(kind=prec)						::	Tfree, Tfree_solvent, Tfree_solute, Temp
real(kind=prec)						::	Temp_solvent, Temp_solute, Texcl_solute, Texcl_solvent


!-----------------------------------------------------------------------
!	Nonbonded pair information
!-----------------------------------------------------------------------
type NB_TYPE
        integer(AI)				:: i, j
        integer					:: cgp_pair ! cgp_pair only used with periodic conditions
        real(kind=prec)                         :: vdWA,vdWB,elec
end type NB_TYPE

type CGP_PAIR_TYPE
        integer(AI)				:: i, j !switching atoms (or equal in case of no switching atoms) of the chargegroups
        real(kind=prec)				:: x, y, z !periodical shifts
end type CGP_PAIR_TYPE

type NBQP_TYPE
        integer(AI)				:: i, j
        integer					:: cgp_pair !this variable only used with periodic conditions
        real(kind=prec)                         :: vdWA,vdWB,elec,score
end type NBQP_TYPE


type NBQ_TYPE
        integer(AI)				:: iq, jq  !q-atom numbers
        logical                                 :: soft
        real(kind=prec)                         :: vdWA,vdWB,elec,score
end type NBQ_TYPE

#ifdef USE_GRID
! added support for grid partitioning of md system
! each grid is the size of NB-Cutoff for the interaction
! with the last grid having no atoms (PBC) or the sphere atoms (SPH)
! grids are build like this
! 1 n+1 2n+1 . n*m+1
! 2 n+2 2n+2 . n*m+2
! 3 n+3 2n+3 . n*m+3
! . .   .    . .
! n 2n  3n   . n*(m+1)
! Paul Bauer 2015

type NB_GRID_POS_TYPE
        real(kind=prec)			:: x,y,z           ! grid starting coordinates
        real(kind=prec)			:: xend,yend,zend  ! grid end coordinates
end type NB_GRID_POS_TYPE
#endif
integer						::	nbpp_pair !current no solute-solute pairs
type(NB_TYPE), allocatable, target::nbpp(:)

type(NB_TYPE), allocatable, target::old_nbpp(:)
integer	::	nbww_pair !current no solvent-solvent pairs, implicit and explicit
type(NB_TYPE), allocatable, target::nbww(:)
integer						::	old_nbww_pair,old_nbpw_pair,old_nbpw_cgp_pair
type(NB_TYPE), allocatable, target::old_nbww(:)
integer						::	nbpw_pair !current no solute-solvent pairs
type(NB_TYPE), allocatable, target::nbpw(:)

type(NB_TYPE), allocatable, target::old_nbpw(:)
integer						::	nbqq_max !max number of q-q pairs in any state
integer(TINY), allocatable	::	qconn(:,:,:) !Q-atom connectivity list
integer						:: nbqq_pair(max_states),nbqqp_pair(max_states)
type(NBQ_TYPE), allocatable ::nbqq(:,:)
type(NBQ_TYPE), allocatable ::nbqqp(:,:)

integer						::	nbqp_pair !current no of qatom-solute pairs
type(NBQP_TYPE), allocatable, target::nbqp(:,:)

type(NBQP_TYPE), allocatable, target::old_nbqp(:,:)
integer						::	nbqw_pair !current no of q-atom-water mol. pairs
type(NBQP_TYPE), allocatable	::	nbqw(:,:)

#ifdef USE_GRID
! new structures for grid 
! Paul Bauer 2015
type(NB_GRID_POS_TYPE),allocatable, target::grid_pp(:)
type(NB_GRID_POS_TYPE),allocatable, target::grid_pw(:)
type(NB_GRID_POS_TYPE),allocatable, target::grid_ww(:)

! new variables for number of grids
integer					::pp_gridnum,pw_gridnum,ww_gridnum
integer					::pp_ndim,pw_ndim,ww_ndim
integer,allocatable			:: ww_igrid(:)
integer,allocatable                     :: pw_igrid(:)
integer,allocatable                     :: pp_igrid(:)
integer,allocatable                     ::grid_pp_ngrp(:)
integer,allocatable                     ::grid_pw_ngrp(:)
integer,allocatable                     ::grid_ww_ngrp(:)
integer,allocatable                     ::grid_pp_grp(:,:)
integer,allocatable                     ::grid_pw_grp(:,:)
integer,allocatable                     ::grid_ww_grp(:,:)

logical,allocatable                     ::grid_pp_int(:,:,:,:)
logical,allocatable                     ::grid_pw_int(:,:,:,:)
logical,allocatable                     ::grid_ww_int(:,:,:,:)

integer                                 :: gridstor_pp,gridstor_pw,gridstor_ww
#endif
!these three used only under periodic conditions
integer						:: nbpp_cgp_pair !number of solute-solute chargegroups interacting
type(CGP_PAIR_TYPE), allocatable:: nbpp_cgp(:)
type(CGP_PAIR_TYPE), allocatable:: old_nbpp_cgp(:)
integer						:: nbpw_cgp_pair
type(CGP_PAIR_TYPE), allocatable :: nbpw_cgp(:)
type(CGP_PAIR_TYPE), allocatable :: old_nbpw_cgp(:)
integer						::	nbqp_cgp_pair
type(CGP_PAIR_TYPE), allocatable :: nbqp_cgp(:)
type(CGP_PAIR_TYPE), allocatable :: old_nbqp_cgp(:)

!Note: all the variables for the MC_volume function are now added here

!special monitoring of pairs
integer (TINY),allocatable  :: special_LJcod(:,:,:,:)

! LRF related variables
integer(AI), allocatable	::	iwhich_cgp(:)


type LRF_TYPE
        real(kind=prec)					::	cgp_cent(3)
        real(kind=prec)					::	phi0
        real(kind=prec)					::	phi1(3)
        real(kind=prec)					::	phi2(9)
        real(kind=prec)					::	phi3(27)
end type LRF_TYPE

type(LRF_TYPE), allocatable	::	lrf(:)
type(LRF_TYPE), allocatable	::	old_lrf(:)  !for constant pressure: MC_volume routine

type(node_assignment_type)	:: calculation_assignment


!shake types & variables
!convergence criterion (fraction of distance)
real(kind=prec), parameter			::	SHAKE_TOL = 0.0001_prec
integer, parameter			::	SHAKE_MAX_ITER = 1000


type SHAKE_BOND_TYPE
        integer(AI)				::	i,j
        real(kind=prec)					::	dist2
        logical					::	ready
end type SHAKE_BOND_TYPE


type SHAKE_MOL_TYPE
        integer					::	nconstraints
        type(SHAKE_BOND_TYPE), allocatable	:: bond(:)
end type SHAKE_MOL_TYPE


integer						::	shake_constraints, shake_molecules
type(SHAKE_MOL_TYPE), allocatable :: shake_mol(:)

!----END OF SHARED VARIABLES

!----START OF PUBLIC SUBROUTINES
contains


!----------------------------------------------------------------------


subroutine simprep_startup
! initialise used modules
call qatom_startup
call trj_startup


! initialise constants
pi = 4.0_prec*atan(one)
deg2rad = pi/180.0_prec

end subroutine simprep_startup


!----------------------------------------------------------------------

subroutine simprep_shutdown
! call used modules' shutdown subroutines
if (use_excluded_groups) then
call excluded_shutdown(ngroups_gc)
end if
call qatom_shutdown
call index_shutdown
call trj_shutdown
end subroutine simprep_shutdown

!----------------------------------------------------------------------

subroutine die_prep(cause)
! args
character(*), optional		:: cause


!
! exit with an error message
!


! local vars
integer						:: i
! flush stuff
integer(4), parameter			:: stdout_unit = 6
! external flush disabled for gfortran
! external flush

if (nodeid .eq. 0) then
        write(*,*)
        call centered_heading('ABNORMAL TERMINATION', '!')
        call close_output_files

        ! apologise
        write(*,'(a)') 'ABNORMAL TERMINATION of program'
        if (present(cause)) then
                write(*,'(79a)') 'Terminating due to ', cause
                endif
        write (*,'(79a)') ('!',i=1,79)
#if defined(CRAY)
        !Cray can't flush stdout...
#elif defined(NO_FLUSH)
		!When you can't flush
#else
        ! flush stdout
        call flush(stdout_unit)
#endif
end if	
! clean up
call md_deallocate


#if defined (USE_MPI)
! abort all processes with exit code 255
call MPI_Abort(MPI_COMM_WORLD, 255, ierr)
#else
! stop with a message to stderr
stop 'Program terminated abnormally'
#endif

end subroutine die

!----------------------------------------------------------------------

! --- Memory management routines
#ifdef USE_GRID
subroutine allocate_grid_pp
allocate(grid_pp(pp_gridnum),stat=alloc_status)
call check_alloc('md atom grid pp')
allocate(grid_pp_ngrp(pp_gridnum),stat=alloc_status)
call check_alloc('md atom grid pp ngroups')
allocate(grid_pp_grp(pp_gridnum,gridstor_pp),stat=alloc_status)
call check_alloc('md atom grid pp groups')
allocate(grid_pp_int(pp_gridnum,pp_ndim,pp_ndim,pp_ndim),stat=alloc_status)
call check_alloc('md atom grid pp interaction matrix')
end subroutine allocate_grid_pp


subroutine allocate_grid_pw
allocate(grid_pw(pw_gridnum),stat=alloc_status)
call check_alloc('md atom grid pw')
allocate(grid_pw_ngrp(pw_gridnum),stat=alloc_status)
call check_alloc('md atom grid pw ngroups')
allocate(grid_pw_grp(pw_gridnum,gridstor_pw),stat=alloc_status)
call check_alloc('md atom grid pw groups')
allocate(grid_pw_int(pw_gridnum,pw_ndim,pw_ndim,pw_ndim),stat=alloc_status)
call check_alloc('md atom grid pw interaction matrix')
end subroutine allocate_grid_pw

subroutine allocate_grid_ww
allocate(grid_ww(ww_gridnum),stat=alloc_status)
call check_alloc('md atom grid ww')
allocate(grid_ww_ngrp(ww_gridnum),stat=alloc_status)
call check_alloc('md atom grid ww ngroups')
allocate(grid_ww_grp(ww_gridnum,gridstor_ww),stat=alloc_status)
call check_alloc('md atom grid ww groups')
allocate(grid_ww_int(ww_gridnum,ww_ndim,ww_ndim,ww_ndim),stat=alloc_status)
call check_alloc('md atom grid ww interaction matrix')
end subroutine allocate_grid_ww
#endif

subroutine allocate_natom_arrays

allocate(x(natom*3), &
                xx(natom*3), &
                v(natom*3), &
                d(natom*3), &
                winv(natom), &
                iqatom(natom), &
                stat=alloc_status)
call check_alloc('atom data arrays')
!set arrays to zero after allocation
!we can not trust the compiler to do this
	x(:)=0.0
	xx(:)=0.0
	v(:)=0.0
	d(:)=0.0
	winv(:)=0.0
	iqatom(:)=0
end subroutine allocate_natom_arrays

!----------------------------------------------------------------------
!Allocate arrays to be used for the Nosé-Hoover thermostat
subroutine allocate_nhchain_arrays

allocate(xnh(numchain), &
                vnh(numchain), &
                qnh(numchain), &
		Gnh(numchain), &
                stat=alloc_status)
call check_alloc('Nose-Hoover data arrays')
end subroutine allocate_nhchain_arrays

!----------------------------------------------------------------------
!Allocate arrays that hold no. pairs per chargegroup.
subroutine allocate_nbxx_per_cgp

 allocate(nbpp_per_cgp(ncgp_solute), & 
          nbww_per_cgp(nwat), &
          nbqp_per_cgp(ncgp_solute), &
          nbqw_per_cgp(nwat), &
          nbpw_per_cgp(ncgp_solute), &
          stat=alloc_status)
 call check_alloc('MPI data arrays')

end subroutine allocate_nbxx_per_cgp

!----------------------------------------------------------------------

subroutine allocate_lrf_arrays

if (use_PBC .and. constant_pressure) then
	allocate(iwhich_cgp(natom), lrf(ncgp), old_lrf(ncgp), stat=alloc_status)
else
	allocate(iwhich_cgp(natom), lrf(ncgp), stat=alloc_status)
end if
call check_alloc('LRF arrays')
end subroutine allocate_lrf_arrays

!----------------------------------------------------------------------
#if defined(USE_MPI)
subroutine allocate_mpi

if(nodeid .eq. 0) then
 allocate(mpi_status(MPI_STATUS_SIZE,numnodes-1), & 
         request_recv(numnodes-1,3), &
         d_recv(natom*3,numnodes-1), &
         E_recv(numnodes-1), &
         EQ_recv(nstates,ene_header%arrays,numnodes-1), &
         stat=alloc_status)
 call check_alloc('MPI data arrays head node')
else
 allocate(E_send(1), &
          EQ_send(nstates,ene_header%arrays), &
          stat=alloc_status)
 call check_alloc('MPI data arrays slave nodes')
end if
#ifdef _OPENMP
allocate(qomp_elec(nstates,ene_header%arrays),qomp_vdw(nstates,ene_header%arrays),stat=alloc_status)
call check_alloc('OMP temporary energies')
#endif
end subroutine allocate_mpi
#endif
!----------------------------------------------------------------------

subroutine allocate_watpol_arrays

allocate(list_sh(n_max_insh, nwpolr_shell), &
nsort(n_max_insh,nwpolr_shell), &
theta(nwat), &
theta0(nwat), &
tdum(nwat), &
stat=alloc_status)
call check_alloc('water polarisation shell arrays')

end subroutine allocate_watpol_arrays
!----------------------------------------------------------------------

subroutine reallocate_watpol_arrays
!local variables
integer(AI), allocatable        ::      old_list_sh(:,:), old_nsort(:,:)
integer				::	old_n_max_insh

old_n_max_insh = n_max_insh
allocate(old_list_sh(old_n_max_insh,nwpolr_shell),old_nsort(n_max_insh, nwpolr_shell),stat=alloc_status)
call check_alloc('water polarisation shell arrays reallocation')

old_list_sh(1:old_n_max_insh,1:nwpolr_shell)=list_sh(1:old_n_max_insh,1:nwpolr_shell)
old_nsort(1:old_n_max_insh,1:nwpolr_shell)=nsort(1:old_n_max_insh,1:nwpolr_shell)

deallocate(list_sh,nsort)

n_max_insh = int(n_max_insh * 1.5_prec) ! make array a big larger, again ...
allocate(list_sh(n_max_insh,nwpolr_shell),nsort(n_max_insh, nwpolr_shell),stat=alloc_status)
call check_alloc('water polarisation shell arrays new allocation')

list_sh(1:old_n_max_insh,1:nwpolr_shell)=old_list_sh(1:old_n_max_insh,1:nwpolr_shell)
nsort(1:old_n_max_insh,1:nwpolr_shell)=old_nsort(1:old_n_max_insh,1:nwpolr_shell)

deallocate(old_list_sh,old_nsort)
call check_alloc('water polarisation shell arrays')

end subroutine reallocate_watpol_arrays

!----------------------------------------------------------------------

subroutine check_alloc(message)

! argument
character(*) message

! local var
character(120) allocmsg

if (alloc_status .ne. 0) then
allocmsg = '>>> Out of memory trying to allocate '//message
call die(allocmsg)
end if
end subroutine check_alloc

!----------------------------------------------------------------------

subroutine md_deallocate
! deallocatde this module's own arrays. Called by shutdown

! local variables for deallocation loop
	integer			:: ii, jj
! use status to avoid errors if not allocated

! atom arrays
deallocate (x, stat=alloc_status)
deallocate (xx, stat=alloc_status)
deallocate (v, stat=alloc_status)
deallocate (d, stat=alloc_status)
deallocate (winv, stat=alloc_status)
deallocate (iqatom, stat=alloc_status)

if(allocated(ljcod)) deallocate(ljcod, stat=alloc_status)

! shake stuff
if(allocated(shake_mol)) then
	do ii=1,nmol
		if (allocated(shake_mol(ii)%bond)) deallocate(shake_mol(ii)%bond)
	end do
	deallocate(shake_mol)
end if
! solvent stuff
if(allocated(aLJ_solv)) deallocate(aLJ_solv)
if(allocated(bLJ_solv)) deallocate(bLJ_solv)
if(allocated(chg_solv)) deallocate(chg_solv)
if(allocated(nonbnd_solv_int)) deallocate(nonbnd_solv_int)
! Nosé-Hoover array
if ( thermostat == NOSEHOOVER) then
deallocate (xnh, stat=alloc_status)
deallocate (vnh, stat=alloc_status)
deallocate (qnh, stat=alloc_status)
deallocate (Gnh, stat=alloc_status)
end if
! nonbond lists
deallocate(nbpp, nbpw, nbww, nbqq, nbqp, nbqw, qconn, stat=alloc_status)
! precompute arrays
if(allocated(pp_precomp)) deallocate(pp_precomp)
if(allocated(pw_precomp)) deallocate(pw_precomp)
if(allocated(ww_precomp)) deallocate(ww_precomp)
if(allocated(qp_precomp)) deallocate(qp_precomp)
if(allocated(qw_precomp)) deallocate(qw_precomp)
if(allocated(qq_precomp)) deallocate(qq_precomp)
! nonbond charge groups
if (allocated(nbpp_per_cgp)) deallocate(nbpp_per_cgp)
if (allocated(nbww_per_cgp)) deallocate(nbww_per_cgp)
if (allocated(nbqp_per_cgp)) deallocate(nbqp_per_cgp)
if (allocated(nbqw_per_cgp)) deallocate(nbqw_per_cgp)
if (allocated(nbpw_per_cgp)) deallocate(nbpw_per_cgp)
! nonbond monitor
if (allocated(monitor_group_int)) deallocate(monitor_group_int)

! stuff for PBC
if (allocated(nbqp_cgp)) deallocate(nbqp_cgp)
if (allocated(nbpp_cgp)) deallocate(nbpp_cgp)
if (allocated(nbpw_cgp)) deallocate(nbpw_cgp)
if (allocated(mvd_mol)) deallocate(mvd_mol)
if (allocated(old_lrf)) deallocate(old_lrf)
if (allocated(mol_mass)) deallocate(mol_mass)
if (allocated(mass)) deallocate(mass)

! watpol arrays
deallocate(wshell, stat=alloc_status)
deallocate(list_sh, stat=alloc_status)
deallocate(nsort, stat=alloc_status)
deallocate(theta, stat=alloc_status)
deallocate(theta0, stat=alloc_status)
deallocate(tdum, stat=alloc_status)

! LRF arrays
deallocate(iwhich_cgp, lrf, stat=alloc_status)

! restraints
deallocate(rstseq, stat=alloc_status)
deallocate(rstpos, stat=alloc_status)
deallocate(rstdis, stat=alloc_status)
deallocate(rstang, stat=alloc_status)
deallocate(rstwal, stat=alloc_status)

! excluded groups
deallocate(ST_gc, stat=alloc_status)
#ifdef _OPENMP
if (allocated(qomp_elec)) deallocate(qomp_elec)
if (allocated(qomp_vdw)) deallocate(qomp_vdw)
#endif
! new file header for energy file
if (allocated(ene_header%types)) then
deallocate(ene_header%types,ene_header%numres,ene_header%gcnum,ene_header%resid)
endif
! header data for energy file
if (allocated(ene_header%types)) deallocate(ene_header%types)
if (allocated(ene_header%numres)) deallocate(ene_header%numres)
if (allocated(ene_header%resid)) deallocate(ene_header%resid)
if (allocated(ene_header%gcnum)) deallocate(ene_header%gcnum)
!variables for MC_volume function
if (allocated(old_nbww)) deallocate(old_nbww)
if (allocated(old_nbpw)) deallocate(old_nbpw)
if (allocated(old_nbpp)) deallocate(old_nbpp)
if (allocated(old_nbqp)) deallocate(old_nbqp)

!more variables for MC_Volume ...
if (allocated(old_nbpp_cgp)) deallocate(old_nbpp_cgp)
if (allocated(old_nbpw_cgp)) deallocate(old_nbpw_cgp)
if (allocated(old_nbqp_cgp)) deallocate(old_nbqp_cgp)
#if defined (USE_MPI)
!MPI arrays
deallocate(nbpp_per_cgp ,stat=alloc_status)
deallocate(nbww_per_cgp ,stat=alloc_status)
deallocate(nbqp_per_cgp ,stat=alloc_status)
deallocate(nbqw_per_cgp ,stat=alloc_status)
deallocate(nbpw_per_cgp ,stat=alloc_status)


if(nodeid .eq. 0) then
   deallocate(d_recv, stat=alloc_status)
   deallocate(E_recv, stat=alloc_status)
   deallocate(EQ_recv, stat=alloc_status)
else
   deallocate(E_send, stat=alloc_status)
   deallocate(EQ_send, stat=alloc_status)
end if
#endif

end subroutine md_deallocate

!-----------------------------------------------------------------------

subroutine reallocate_nonbondlist_pp

! variables
type(NB_TYPE), allocatable	:: old_nbxx(:)
integer						:: old_max

! copy
old_max = calculation_assignment%pp%max
allocate(old_nbxx(old_max), stat=alloc_status)
call check_alloc('reallocating non-bonded interaction list')
old_nbxx(1:old_max) = nbpp(1:old_max)

! deallocate and copy back
deallocate(nbpp)
calculation_assignment%pp%max = int(calculation_assignment%pp%max * 1.05) + 200
allocate(nbpp(calculation_assignment%pp%max), stat = alloc_status)
call check_alloc('reallocating non-bonded interaction list')
nbpp(1:old_max) = old_nbxx(1:old_max)

! deallocate copy
deallocate(old_nbxx)

! tell the world
write (*,100) calculation_assignment%pp%max
100 format('>>> reallocating p-p pair list, new max is ', i8)

end subroutine reallocate_nonbondlist_pp
!----------------------------------------------------------------------
subroutine reallocate_nbpp_cgp

! variables
type(CGP_PAIR_TYPE), allocatable	:: old_nbxx(:)
integer						:: old_max, new_max

! copy
old_max = size(nbpp_cgp, 1)
allocate(old_nbxx(old_max), stat=alloc_status)
call check_alloc('reallocating non-bonded charge group list')
old_nbxx(1:old_max) = nbpp_cgp(1:old_max)

! deallocate and copy back
deallocate(nbpp_cgp)
new_max = int( old_max*1.05) + 200
allocate(nbpp_cgp(new_max), stat = alloc_status)
call check_alloc('reallocating non-bonded charge group list')
nbpp_cgp(1:old_max) = old_nbxx(1:old_max)

! deallocate copy
deallocate(old_nbxx)

! tell the world
write (*,100) new_max
100 format('>>> reallocating p-p charge group pair list, new max is ', i8)

end subroutine reallocate_nbpp_cgp
!----------------------------------------------------------------------
subroutine reallocate_nbpw_cgp

! variables
type(CGP_PAIR_TYPE), allocatable	:: old_nbxx(:)
integer						:: old_max, new_max

! copy
old_max = size(nbpw_cgp, 1)
allocate(old_nbxx(old_max), stat=alloc_status)
call check_alloc('reallocating non-bonded charge group list')
old_nbxx(1:old_max) = nbpw_cgp(1:old_max)

! deallocate and copy back
deallocate(nbpw_cgp)
new_max = int( old_max*1.05) + 200
allocate(nbpw_cgp(new_max), stat = alloc_status)
call check_alloc('reallocating non-bonded charge group list')
nbpw_cgp(1:old_max) = old_nbxx(1:old_max)

! deallocate copy
deallocate(old_nbxx)

! tell the world
write (*,100) new_max
100 format('>>> reallocating p-w charge group pair list, new max is ', i8)

end subroutine reallocate_nbpw_cgp
!----------------------------------------------------------------------
subroutine reallocate_nbqp_cgp

! variables
type(CGP_PAIR_TYPE), allocatable	:: old_nbxx(:)
integer						:: old_max, new_max

! copy
old_max = size(nbqp_cgp, 1)
allocate(old_nbxx(old_max), stat=alloc_status)
call check_alloc('reallocating non-bonded charge group list')
old_nbxx(1:old_max) = nbqp_cgp(1:old_max)

! deallocate and copy back
deallocate(nbqp_cgp)
new_max = int( old_max*1.05) + 200
allocate(nbqp_cgp(new_max), stat = alloc_status)
call check_alloc('reallocating non-bonded charge group list')
nbqp_cgp(1:old_max) = old_nbxx(1:old_max)

! deallocate copy
deallocate(old_nbxx)

! tell the world
write (*,100) new_max
100 format('>>> reallocating q-p charge group pair list, new max is ', i8)

end subroutine reallocate_nbqp_cgp

!----------------------------------------------------------------------

subroutine reallocate_nonbondlist_pw
! variables
type(NB_TYPE), allocatable	:: old_nbxx(:)
integer						:: old_max

! copy
old_max = calculation_assignment%pw%max
allocate(old_nbxx(old_max), stat=alloc_status)
call check_alloc('reallocating non-bonded interaction list')
old_nbxx(1:old_max) = nbpw(1:old_max)

! deallocate and copy back
deallocate(nbpw)
calculation_assignment%pw%max = int(calculation_assignment%pw%max * 1.05) + 200
allocate(nbpw(calculation_assignment%pw%max), stat = alloc_status)
call check_alloc('reallocating non-bonded interaction list')
nbpw(1:old_max) = old_nbxx(1:old_max)

! deallocate copy
deallocate(old_nbxx)

! tell the world
write (*,100) calculation_assignment%pw%max
100 format('>>> reallocating p-w pair list, new max is ', i8)

end subroutine reallocate_nonbondlist_pw

!----------------------------------------------------------------------

subroutine reallocate_nonbondlist_qp
! variables
type(NBQP_TYPE), allocatable	:: old_nbxx(:,:)
integer						:: old_max, is

old_max = calculation_assignment%qp%max
allocate(old_nbxx(old_max,nstates), stat=alloc_status)
call check_alloc('reallocating non-bonded interaction list')

do is=1,nstates
! copy
old_nbxx(1:old_max,is) = nbqp(1:old_max,is)
end do
! deallocate and copy back
deallocate(nbqp)
calculation_assignment%qp%max = int(calculation_assignment%qp%max * 1.05) + 200
allocate(nbqp(calculation_assignment%qp%max,nstates), stat = alloc_status)
call check_alloc('reallocating non-bonded interaction list')
do is=1,nstates
! copy
nbqp(1:old_max,is) = old_nbxx(1:old_max,is)
end do

! deallocate copy
deallocate(old_nbxx)

! tell the world
write (*,100) calculation_assignment%qp%max
100 format('>>> reallocating q-s pair list, new max is ', i8)

end subroutine reallocate_nonbondlist_qp

!----------------------------------------------------------------------
subroutine reallocate_nonbondlist_qw
! variables
type(NBQP_TYPE), allocatable    :: old_nbxx(:,:)
integer                                         :: old_max ,is

! copy
old_max = calculation_assignment%qw%max
allocate(old_nbxx(old_max,nstates), stat=alloc_status)
call check_alloc('reallocating non-bonded interaction list')
do is=1,nstates
! copy
old_nbxx(1:old_max,is) = nbqw(1:old_max,is)
end do

! deallocate and copy back
deallocate(nbqw)
calculation_assignment%qw%max = int(calculation_assignment%qw%max * 1.05) + 200
allocate(nbqw(calculation_assignment%qw%max,nstates), stat = alloc_status)
call check_alloc('reallocating non-bonded interaction list')
do is=1,nstates
! copy
nbqw(1:old_max,is) = old_nbxx(1:old_max,is)
end do

! deallocate copy
deallocate(old_nbxx)

! tell the world
write (*,100) calculation_assignment%qw%max
100 format('>>> reallocating q-w pair list, new max is ', i8)

end subroutine reallocate_nonbondlist_qw

subroutine reallocate_nonbondlist_ww
! variables
type(NB_TYPE), allocatable		:: old_nbxx(:)
integer						:: old_max

! copy
allocate(old_nbxx(calculation_assignment%ww%max), stat=alloc_status)
call check_alloc('reallocating non-bonded interaction list')
old_nbxx(1:calculation_assignment%ww%max) = nbww(1:calculation_assignment%ww%max)
old_max = calculation_assignment%ww%max


! deallocate and copy back
deallocate(nbww)
calculation_assignment%ww%max = int(calculation_assignment%ww%max * 1.05) + 200 + nwat
allocate(nbww(calculation_assignment%ww%max), stat = alloc_status)
call check_alloc('reallocating non-bonded interaction list')
nbww(1:old_max) = old_nbxx(1:old_max)

! deallocate copy
deallocate(old_nbxx)

! tell the world
write (*,100) calculation_assignment%ww%max
100 format('>>> reallocating w-w pair list, new max is ', i8)

end subroutine reallocate_nonbondlist_ww

!----------------------------------------------------------------------

#ifdef USE_GRID
subroutine create_grid_pp
! *** local variables
integer						:: i,j,k,ndim,num
integer						:: i3
real(kind=prec)					:: xtmp,ytmp,ztmp
real(kind=prec)					:: xmax,ymax,zmax
real(kind=prec)					:: xmin,ymin,zmin
!have 1.75 times the cutoff length to account for large charge groups in solute residues
real(kind=prec)					:: gridRc
integer						:: li,lj,lk,ui,uj,uk

gridRc = 1.75_prec * Rcpp

if(use_PBC) then
! when using period boundary conditions, the grids are generated accoring to the box centre and the desired  system size
! this depends on the volume, so we need to be careful with the coordinates to not miss any atom
! so the coordinate information needs to be checked/reset in populate grids
! Paul Bauer 2015

! one problem could be the large movement of molecules if they are not put back in the box
! this needs to be checked for, too
! also, the coordinates of each grid need to be corrected after a box movement

! number of one dimensional spacings cubed
pp_gridnum=int(ceiling(((boxlength(1)+boxlength(2)+boxlength(3))/3)/gridRc)**3)
ndim=ceiling(((boxlength(1)+boxlength(2)+boxlength(3))/3)/gridRc)
pp_ndim=ndim
if (ncgp_solute .gt. 500) then
gridstor_pp = (ncgp_solute*8/pp_gridnum)+100
else
gridstor_pp = ncgp_solute
endif
call allocate_grid_pp

! set starting coordinates for first grid
xtmp=boxcentre(1)-boxlength(1)/2
ytmp=boxcentre(2)-boxlength(2)/2
ztmp=boxcentre(3)-boxlength(3)/2

num=0
do i=1,ndim
 ytmp=boxcentre(2)-boxlength(2)/2
 do j=1,ndim
  ztmp=boxcentre(3)-boxlength(3)/2
  do k=1,ndim
   num=num+1

   grid_pp_int(num,:,:,:)=.false.
   grid_pp_int(num,i,j,k)=.true.

! stuff to account for PBC wraparound for the interaction matrix
   if(i+1.le.ndim) then
    ui=i+1
   else
    ui=1
   end if
   if(i-1.ge.1) then
    li=i-1
   else
    li=ndim
   end if
   if(j+1.le.ndim) then
    uj=j+1
   else
    uj=1
   end if
   if(j-1.ge.1) then
    lj=j-1
   else
    lj=ndim
   end if
   if(k+1.le.ndim) then
    uk=k+1
   else
    uk=1
   end if
   if(k-1.ge.1) then
    lk=k-1
   else
    lk=ndim
   end if
   
! generating interaction matrix
   grid_pp_int(num,ui,uj,uk)=.true.
   grid_pp_int(num,ui,uj,lk)=.true.
   grid_pp_int(num,ui,lj,uk)=.true.
   grid_pp_int(num,ui,lj,lk)=.true.
   grid_pp_int(num,li,uj,uk)=.true.
   grid_pp_int(num,li,uj,lk)=.true.
   grid_pp_int(num,li,lj,uk)=.true.
   grid_pp_int(num,li,lj,lk)=.true.
   grid_pp_int(num,ui,uj,k)=.true.
   grid_pp_int(num,ui,lj,k)=.true.
   grid_pp_int(num,li,uj,k)=.true.
   grid_pp_int(num,li,lj,k)=.true.
   grid_pp_int(num,ui,j,uk)=.true.
   grid_pp_int(num,ui,j,lk)=.true.
   grid_pp_int(num,li,j,uk)=.true.
   grid_pp_int(num,li,j,lk)=.true.
   grid_pp_int(num,i,uj,uk)=.true.
   grid_pp_int(num,i,uj,lk)=.true.
   grid_pp_int(num,i,lj,uk)=.true.
   grid_pp_int(num,i,lj,lk)=.true.
   grid_pp_int(num,ui,j,k)=.true.
   grid_pp_int(num,li,j,k)=.true.
   grid_pp_int(num,i,j,uk)=.true.
   grid_pp_int(num,i,j,lk)=.true.
   grid_pp_int(num,i,uj,k)=.true.
   grid_pp_int(num,i,lj,k)=.true.

   grid_pp(num)%x=xtmp
   grid_pp(num)%y=ytmp
   grid_pp(num)%z=ztmp

   grid_pp(num)%xend=xtmp+gridRc
   grid_pp(num)%yend=ytmp+gridRc
   grid_pp(num)%zend=ztmp+gridRc

   if (k .eq. ndim) then
    if (grid_pp(num)%zend .lt. (boxcentre(3)+(boxlength(3)/2))) then
     grid_pp(num)%zend=boxcentre(3)+(boxlength(3)/2)
    end if
   end if
   ztmp=ztmp+gridRc
  end do
  if (j .eq. ndim) then
   if(grid_pp(num)%yend .lt. (boxcentre(2)+(boxlength(2)/2))) then
    grid_pp(num)%yend=boxcentre(2)+(boxlength(2)/2)
   end if
  end if
  ytmp=ytmp+gridRc
 end do
 if (i .eq. ndim) then
  if(grid_pp(num)%xend .lt. (boxcentre(1)+(boxlength(1)/2))) then
   grid_pp(num)%xend=boxcentre(1)+(boxlength(1)/2)
  end if
 end if
 xtmp=xtmp+gridRc
end do

else
! we are using spherical boundaries, meaning that we now have some problems
! we need to know the largest and smallest values for all coordiantes to generate the grid
i = 1
i3 = i*3-3
xmax=x(i3+1)
ymax=x(i3+2)
zmax=x(i3+3)
xmin=xmax
ymin=ymax
zmin=zmax
! now iterate over all cooridnates to find min and max
! atom one is already read in, don't need to read it again
do i=2,natom
 i3=i*3-3
 if (xmin .gt. x(i3+1)) xmin=x(i3+1)
 if (ymin .gt. x(i3+2)) ymin=x(i3+2)
 if (zmin .gt. x(i3+3)) zmin=x(i3+3)
 if (xmax .lt. x(i3+1)) xmax=x(i3+1)
 if (ymax .lt. x(i3+2)) ymax=x(i3+2)
 if (zmax .lt. x(i3+3)) zmax=x(i3+3)
end do
! now add some buffer to each of them -> half of Rcpp
xmin=xmin-(gridRc/3)
ymin=ymin-(gridRc/3)
zmin=zmin-(gridRc/3)
xmax=xmax+(gridRc/3)
ymax=ymax+(gridRc/3)
zmax=zmax+(gridRc/3)
! the system of boxes will be cubic in any case, so now find the largest number of 
! grid spaces one any side
ndim=ceiling((xmax-xmin)/gridRc)
if (ceiling((ymax-ymin)/gridRc) .gt. ndim ) then
 ndim=ceiling((ymax-ymin)/gridRc)
else if (ceiling((zmax-zmin)/gridRc) .gt. ndim ) then
 ndim=(ceiling((zmax-zmin)/gridRc))
end if
pp_gridnum=ndim**3
if (ncgp_solute .gt. 500) then
gridstor_pp = (ncgp_solute*8/pp_gridnum)+100
else
gridstor_pp = ncgp_solute
endif
pp_ndim=ndim
call allocate_grid_pp

xtmp=xmin
ytmp=ymin
ztmp=zmin

num=1
do i=1,ndim
 ytmp=ymin
 do j=1,ndim
  ztmp=zmin
  do k=1,ndim
! now fill the interaction matrix
! no wraparound in this case, so be careful
! when there would be wraparound, just point to my own cell again

grid_pp_int(num,:,:,:)=.false.
grid_pp_int(num,i,j,k)=.true.


   if((i+1).le.ndim) then
    ui=i+1
   else
    ui=i
   end if
   if((i-1).ge.1) then
    li=i-1
   else
    li=i
   end if
   if((j+1).le.ndim) then
    uj=j+1
   else
    uj=j
   end if
   if((j-1).ge.1) then
    lj=j-1
   else
    lj=j
   end if
   if((k+1).le.ndim) then
    uk=k+1
   else
    uk=k
   end if
   if((k-1).ge.1) then
    lk=k-1
   else
    lk=k
   end if

   grid_pp_int(num,ui,uj,uk)=.true.
   grid_pp_int(num,ui,uj,lk)=.true.
   grid_pp_int(num,ui,lj,uk)=.true.
   grid_pp_int(num,ui,lj,lk)=.true.
   grid_pp_int(num,li,uj,uk)=.true.
   grid_pp_int(num,li,uj,lk)=.true.
   grid_pp_int(num,li,lj,uk)=.true.
   grid_pp_int(num,li,lj,lk)=.true.
   grid_pp_int(num,ui,uj,k)=.true.
   grid_pp_int(num,ui,lj,k)=.true.
   grid_pp_int(num,li,uj,k)=.true.
   grid_pp_int(num,li,lj,k)=.true.
   grid_pp_int(num,ui,j,uk)=.true.
   grid_pp_int(num,ui,j,lk)=.true.
   grid_pp_int(num,li,j,uk)=.true.
   grid_pp_int(num,li,j,lk)=.true.
   grid_pp_int(num,i,uj,uk)=.true.
   grid_pp_int(num,i,uj,lk)=.true.
   grid_pp_int(num,i,lj,uk)=.true.
   grid_pp_int(num,i,lj,lk)=.true.
   grid_pp_int(num,ui,j,k)=.true.
   grid_pp_int(num,li,j,k)=.true.
   grid_pp_int(num,i,j,uk)=.true.
   grid_pp_int(num,i,j,lk)=.true.
   grid_pp_int(num,i,uj,k)=.true.
   grid_pp_int(num,i,lj,k)=.true.

   grid_pp(num)%x=xtmp
   grid_pp(num)%y=ytmp
   grid_pp(num)%z=ztmp

   grid_pp(num)%xend=xtmp+gridRc
   grid_pp(num)%yend=ytmp+gridRc
   grid_pp(num)%zend=ztmp+gridRc

   num=num+1
   ztmp=ztmp+gridRc
  end do
  ytmp=ytmp+gridRc
 end do
 xtmp=xtmp+gridRc
end do

end if

end subroutine create_grid_pp 

subroutine create_grid_pw
! *** local variables
integer                                         :: i,j,k,ndim,num
integer                                         :: i3
real(kind=prec)                                 :: xtmp,ytmp,ztmp
real(kind=prec)                                 :: xmax,ymax,zmax
real(kind=prec)                                 :: xmin,ymin,zmin
! again make grid larger than cutoff to account for charge groups
! can be smaller this time because we only have water as one of them
real(kind=prec)					:: gridRc
integer                                         :: li,lj,lk,ui,uj,uk

gridRc = 1.5_prec * Rcpw

if(use_PBC) then
! when using period boundary conditions, the grids are generated accoring to the box centre and the desired  system size
! this depends on the volume, so we need to be careful with the coordinates to not miss any atom
! so the coordinate information needs to be checked/reset in populate grids
! Paul Bauer 2015

! one problem could be the large movement of molecules if they are not put back in the box
! this needs to be checked for, too
! also, the coordinates of each grid need to be corrected after a box movement

! number of one dimensional spacings cubed
pw_gridnum=int(ceiling(((boxlength(1)+boxlength(2)+boxlength(3))/3)/gridRc)**3)
ndim=ceiling(((boxlength(1)+boxlength(2)+boxlength(3))/3)/gridRc)
pw_ndim=ndim
if (nwat .gt. 500) then
gridstor_pw = (nwat*8/pw_gridnum)+100
else
gridstor_pw = nwat
endif
call allocate_grid_pw

! set starting coordinates for first grid
xtmp=boxcentre(1)-(boxlength(1)/2)
ytmp=boxcentre(2)-(boxlength(2)/2)
ztmp=boxcentre(3)-(boxlength(3)/2)


num=0
do i=1,ndim
 ytmp=boxcentre(2)-(boxlength(2)/2)
 do j=1,ndim
  ztmp=boxcentre(3)-(boxlength(3)/2)
  do k=1,ndim
   num=num+1

   grid_pw_int(num,:,:,:)=.false.
   grid_pw_int(num,i,j,k)=.true.

! stuff to account for PBC wraparound for the interaction matrix
   if(i+1.le.ndim) then
    ui=i+1
   else
    ui=1
   end if
   if(i-1.ge.1) then
    li=i-1
   else
    li=ndim
   end if
   if(j+1.le.ndim) then
    uj=j+1
   else
    uj=1
   end if
   if(j-1.ge.1) then
    lj=j-1
   else
    lj=ndim
   end if
   if(k+1.le.ndim) then
    uk=k+1
   else
    uk=1
   end if
   if(k-1.ge.1) then
    lk=k-1
   else
    lk=ndim
   end if

! generating interaction matrix
   grid_pw_int(num,ui,uj,uk)=.true.
   grid_pw_int(num,ui,uj,lk)=.true.
   grid_pw_int(num,ui,lj,uk)=.true.
   grid_pw_int(num,ui,lj,lk)=.true.
   grid_pw_int(num,li,uj,uk)=.true.
   grid_pw_int(num,li,uj,lk)=.true.
   grid_pw_int(num,li,lj,uk)=.true.
   grid_pw_int(num,li,lj,lk)=.true.
   grid_pw_int(num,ui,uj,k)=.true.
   grid_pw_int(num,ui,lj,k)=.true.
   grid_pw_int(num,li,uj,k)=.true.
   grid_pw_int(num,li,lj,k)=.true.
   grid_pw_int(num,ui,j,uk)=.true.
   grid_pw_int(num,ui,j,lk)=.true.
   grid_pw_int(num,li,j,uk)=.true.
   grid_pw_int(num,li,j,lk)=.true.
   grid_pw_int(num,i,uj,uk)=.true.
   grid_pw_int(num,i,uj,lk)=.true.
   grid_pw_int(num,i,lj,uk)=.true.
   grid_pw_int(num,i,lj,lk)=.true.
   grid_pw_int(num,ui,j,k)=.true.
   grid_pw_int(num,li,j,k)=.true.
   grid_pw_int(num,i,j,uk)=.true.
   grid_pw_int(num,i,j,lk)=.true.
   grid_pw_int(num,i,uj,k)=.true.
   grid_pw_int(num,i,lj,k)=.true.


   grid_pw(num)%x=xtmp
   grid_pw(num)%y=ytmp
   grid_pw(num)%z=ztmp

   grid_pw(num)%xend=xtmp+gridRc
   grid_pw(num)%yend=ytmp+gridRc
   grid_pw(num)%zend=ztmp+gridRc

   if (k .eq. ndim) then
    if (grid_pw(num)%zend .lt. (boxcentre(3)+(boxlength(3)/2))) then
     grid_pw(num)%zend=boxcentre(3)+(boxlength(3)/2)
    end if
   end if
   ztmp=ztmp+gridRc
  end do
  if (j .eq. ndim) then
   if(grid_pw(num)%yend .lt. (boxcentre(2)+(boxlength(2)/2))) then
    grid_pw(num)%yend=boxcentre(2)+(boxlength(2)/2)
   end if
  end if
  ytmp=ytmp+gridRc
 end do
 if (i .eq. ndim) then
  if(grid_pw(num)%xend .lt. (boxcentre(1)+(boxlength(1)/2))) then
   grid_pw(num)%xend=boxcentre(1)+(boxlength(1)/2)
  end if
 end if
 xtmp=xtmp+gridRc
end do

else
! we are using spherical boundaries, meaning that we now have some problems
! we need to know the largest and smallest values for all coordiantes to generate the grid
i = 1
i3 = i*3-3
xmax=x(i3+1)
ymax=x(i3+2)
zmax=x(i3+3)
xmin=xmax
ymin=ymax
zmin=zmax
! now iterate over all cooridnates to find min and max
! atom one is already read in, don't need to read it again
do i=2,natom
 i3=i*3-3
 if (xmin .gt. x(i3+1)) xmin=x(i3+1)
 if (ymin .gt. x(i3+2)) ymin=x(i3+2)
 if (zmin .gt. x(i3+3)) zmin=x(i3+3)
 if (xmax .lt. x(i3+1)) xmax=x(i3+1)
 if (ymax .lt. x(i3+2)) ymax=x(i3+2)
 if (zmax .lt. x(i3+3)) zmax=x(i3+3)
end do
! now add some buffer to each of them -> half of Rcpp
xmin=xmin-(gridRc/3)
ymin=ymin-(gridRc/3)
zmin=zmin-(gridRc/3)
xmax=xmax+(gridRc/3)
ymax=ymax+(gridRc/3)
zmax=zmax+(gridRc/3)
! the system of boxes will be cubic in any case, so now find the largest number of 
! grid spaces one any side
ndim=ceiling((xmax-xmin)/gridRc)
if (ceiling((ymax-ymin)/gridRc) .gt. ndim ) then
 ndim=ceiling((ymax-ymin)/gridRc)
else if (ceiling((zmax-zmin)/gridRc) .gt. ndim ) then
 ndim=(ceiling((zmax-zmin)/gridRc))
end if
pw_gridnum=ndim**3
pw_ndim=ndim
if (nwat .gt. 500) then
gridstor_pw = (nwat*8/pw_gridnum)+100
else
gridstor_pw = nwat
endif
call allocate_grid_pw

xtmp=xmin
ytmp=ymin
ztmp=zmin

num=1
do i=1,ndim
 ytmp = ymin
 do j=1,ndim
  ztmp = zmin
  do k=1,ndim
! now fill the interaction matrix
! no wraparound in this case, so be careful
! when there would be wraparound, just point to my own cell again

grid_pw_int(num,:,:,:)=.false.
grid_pw_int(num,i,j,k)=.true.


   if(i+1.le.ndim) then
    ui=i+1
   else
    ui=i
   end if
   if(i-1.ge.1) then
    li=i-1
   else
    li=i
   end if
   if(j+1.le.ndim) then
    uj=j+1
   else
    uj=j
   end if
   if(j-1.ge.1) then
    lj=j-1
   else
    lj=j
   end if
   if(k+1.le.ndim) then
    uk=k+1
   else
    uk=k
   end if
   if(k-1.ge.1) then
    lk=k-1
   else
    lk=k
   end if

   grid_pw_int(num,ui,uj,uk)=.true.
   grid_pw_int(num,ui,uj,lk)=.true.
   grid_pw_int(num,ui,lj,uk)=.true.
   grid_pw_int(num,ui,lj,lk)=.true.
   grid_pw_int(num,li,uj,uk)=.true.
   grid_pw_int(num,li,uj,lk)=.true.
   grid_pw_int(num,li,lj,uk)=.true.
   grid_pw_int(num,li,lj,lk)=.true.
   grid_pw_int(num,ui,uj,k)=.true.
   grid_pw_int(num,ui,lj,k)=.true.
   grid_pw_int(num,li,uj,k)=.true.
   grid_pw_int(num,li,lj,k)=.true.
   grid_pw_int(num,ui,j,uk)=.true.
   grid_pw_int(num,ui,j,lk)=.true.
   grid_pw_int(num,li,j,uk)=.true.
   grid_pw_int(num,li,j,lk)=.true.
   grid_pw_int(num,i,uj,uk)=.true.
   grid_pw_int(num,i,uj,lk)=.true.
   grid_pw_int(num,i,lj,uk)=.true.
   grid_pw_int(num,i,lj,lk)=.true.
   grid_pw_int(num,ui,j,k)=.true.
   grid_pw_int(num,li,j,k)=.true.
   grid_pw_int(num,i,j,uk)=.true.
   grid_pw_int(num,i,j,lk)=.true.
   grid_pw_int(num,i,uj,k)=.true.
   grid_pw_int(num,i,lj,k)=.true.

   grid_pw(num)%x=xtmp
   grid_pw(num)%y=ytmp
   grid_pw(num)%z=ztmp

   grid_pw(num)%xend=xtmp+gridRc
   grid_pw(num)%yend=ytmp+gridRc
   grid_pw(num)%zend=ztmp+gridRc

   num=num+1
   ztmp=ztmp+gridRc
  end do
  ytmp=ytmp+gridRc
 end do
 xtmp=xtmp+gridRc
end do

end if
end subroutine create_grid_pw

subroutine create_grid_ww
! *** local variables
integer                                         :: i,j,k,ndim,num
integer                                         :: i3
real(kind=prec)                                 :: xtmp,ytmp,ztmp
real(kind=prec)                                 :: xmax,ymax,zmax
real(kind=prec)                                 :: xmin,ymin,zmin
! again grids are slightly larger than cutoff, can be onle 1.25 because we only have water molecules
real(kind=prec)					:: gridRc 
integer                                         :: li,lj,lk,ui,uj,uk

gridRc = 1.25_prec * Rcww

if(use_PBC) then
! when using period boundary conditions, the grids are generated accoring to the box centre and the desired  system size
! this depends on the volume, so we need to be careful with the coordinates to not miss any atom
! so the coordinate information needs to be checked/reset in populate grids
! Paul Bauer 2015

! one problem could be the large movement of molecules if they are not put back in the box
! this needs to be checked for, too
! also, the coordinates of each grid need to be corrected after a box movement

! number of one dimensional spacings cubed
ww_gridnum=int(ceiling(((boxlength(1)+boxlength(2)+boxlength(3))/3)/gridRc)**3)
ndim=ceiling(((boxlength(1)+boxlength(2)+boxlength(3))/3)/gridRc)
ww_ndim=ndim
if (nwat .gt. 500) then
gridstor_ww = (nwat*8/ww_gridnum)+100
else
gridstor_ww = nwat
endif
call allocate_grid_ww
! set starting coordinates for first grid
xtmp=boxcentre(1)-(boxlength(1)/2)
ytmp=boxcentre(2)-(boxlength(2)/2)
ztmp=boxcentre(3)-(boxlength(3)/2)


num=0
do i=1,ndim
 ytmp=boxcentre(2)-(boxlength(2)/2)
 do j=1,ndim
  ztmp=boxcentre(3)-(boxlength(3)/2)
  do k=1,ndim
   num=num+1

   grid_ww_int(num,:,:,:)=.false.
   grid_ww_int(num,i,j,k)=.true.

! stuff to account for PBC wraparound for the interaction matrix
   if(i+1.le.ndim) then
    ui=i+1
   else
    ui=1
   end if
   if(i-1.ge.1) then
    li=i-1
   else
    li=ndim
   end if
   if(j+1.le.ndim) then
    uj=j+1
   else
    uj=1
   end if
   if(j-1.ge.1) then
    lj=j-1
   else
    lj=ndim
   end if
   if(k+1.le.ndim) then
    uk=k+1
   else
    uk=1
   end if
   if(k-1.ge.1) then
    lk=k-1
   else
    lk=ndim
   end if

! generating interaction matrix
   grid_ww_int(num,ui,uj,uk)=.true.
   grid_ww_int(num,ui,uj,lk)=.true.
   grid_ww_int(num,ui,lj,uk)=.true.
   grid_ww_int(num,ui,lj,lk)=.true.
   grid_ww_int(num,li,uj,uk)=.true.
   grid_ww_int(num,li,uj,lk)=.true.
   grid_ww_int(num,li,lj,uk)=.true.
   grid_ww_int(num,li,lj,lk)=.true.
   grid_ww_int(num,ui,uj,k)=.true.
   grid_ww_int(num,ui,lj,k)=.true.
   grid_ww_int(num,li,uj,k)=.true.
   grid_ww_int(num,li,lj,k)=.true.
   grid_ww_int(num,ui,j,uk)=.true.
   grid_ww_int(num,ui,j,lk)=.true.
   grid_ww_int(num,li,j,uk)=.true.
   grid_ww_int(num,li,j,lk)=.true.
   grid_ww_int(num,i,uj,uk)=.true.
   grid_ww_int(num,i,uj,lk)=.true.
   grid_ww_int(num,i,lj,uk)=.true.
   grid_ww_int(num,i,lj,lk)=.true.
   grid_ww_int(num,ui,j,k)=.true.
   grid_ww_int(num,li,j,k)=.true.
   grid_ww_int(num,i,j,uk)=.true.
   grid_ww_int(num,i,j,lk)=.true.
   grid_ww_int(num,i,uj,k)=.true.
   grid_ww_int(num,i,lj,k)=.true.

   grid_ww(num)%x=xtmp
   grid_ww(num)%y=ytmp
   grid_ww(num)%z=ztmp
   grid_ww(num)%xend=xtmp+gridRc
   grid_ww(num)%yend=ytmp+gridRc
   grid_ww(num)%zend=ztmp+gridRc

   if (k .eq. ndim) then
    if (grid_ww(num)%zend .lt. (boxcentre(3)+(boxlength(3)/2))) then
     grid_ww(num)%zend=boxcentre(3)+(boxlength(3)/2)
    end if
   end if
   ztmp=ztmp+gridRc
  end do
  if (j .eq. ndim) then
   if(grid_ww(num)%yend .lt. (boxcentre(2)+(boxlength(2)/2))) then
    grid_ww(num)%yend=boxcentre(2)+(boxlength(2)/2)
   end if
  end if
  ytmp=ytmp+gridRc
 end do
 if (i .eq. ndim) then
  if(grid_ww(num)%xend .lt. (boxcentre(1)+(boxlength(1)/2))) then
   grid_ww(num)%xend=boxcentre(1)+(boxlength(1)/2)
  end if
 end if
 xtmp=xtmp+gridRc
end do

else
! we are using spherical boundaries, meaning that we now have some problems
! we need to know the largest and smallest values for all coordiantes to generate the grid
i = 1
i3 = i*3-3
xmax=x(i3+1)
ymax=x(i3+2)
zmax=x(i3+3)
xmin=xmax
ymin=ymax
zmin=zmax
! now iterate over all cooridnates to find min and max
! atom one is already read in, don't need to read it again
do i=2,natom
 i3=i*3-3
 if (xmin .gt. x(i3+1)) xmin=x(i3+1)
 if (ymin .gt. x(i3+2)) ymin=x(i3+2)
 if (zmin .gt. x(i3+3)) zmin=x(i3+3)
 if (xmax .lt. x(i3+1)) xmax=x(i3+1)
 if (ymax .lt. x(i3+2)) ymax=x(i3+2)
 if (zmax .lt. x(i3+3)) zmax=x(i3+3)
end do
! now add some buffer to each of them -> half of Rcpp
xmin=xmin-(gridRc/3)
ymin=ymin-(gridRc/3)
zmin=zmin-(gridRc/3)
xmax=xmax+(gridRc/3)
ymax=ymax+(gridRc/3)
zmax=zmax+(gridRc/3)
! the system of boxes will be cubic in any case, so now find the largest number of 
! grid spaces one any side
ndim=ceiling((xmax-xmin)/gridRc)
if (ceiling((ymax-ymin)/gridRc) .gt. ndim ) then
 ndim=ceiling((ymax-ymin)/gridRc)
else if (ceiling((zmax-zmin)/gridRc) .gt. ndim ) then
 ndim=(ceiling((zmax-zmin)/gridRc))
end if
ww_gridnum=ndim**3
ww_ndim=ndim
if (nwat .gt. 500) then
gridstor_ww = (nwat*8/ww_gridnum)+100
else
gridstor_ww = nwat
endif
call allocate_grid_ww
xtmp=xmin
ytmp=ymin
ztmp=zmin

num=1
do i=1,ndim
 ytmp = ymin
 do j=1,ndim
  ztmp = zmin
  do k=1,ndim
! now fill the interaction matrix
! no wraparound in this case, so be careful
! when there would be wraparound, just point to my own cell again

grid_ww_int(num,:,:,:)=.false.
grid_ww_int(num,i,j,k)=.true.

   if(i+1.le.ndim) then
    ui=i+1
   else
    ui=i
   end if
   if(i-1.ge.1) then
    li=i-1
   else
    li=i
   end if
   if(j+1.le.ndim) then
    uj=j+1
   else
    uj=j
   end if
   if(j-1.ge.1) then
    lj=j-1
   else
    lj=j
   end if
   if(k+1.le.ndim) then
    uk=k+1
   else
    uk=k
   end if
   if(k-1.ge.1) then
    lk=k-1
   else
    lk=k
   end if

   grid_ww_int(num,ui,uj,uk)=.true.
   grid_ww_int(num,ui,uj,lk)=.true.
   grid_ww_int(num,ui,lj,uk)=.true.
   grid_ww_int(num,ui,lj,lk)=.true.
   grid_ww_int(num,li,uj,uk)=.true.
   grid_ww_int(num,li,uj,lk)=.true.
   grid_ww_int(num,li,lj,uk)=.true.
   grid_ww_int(num,li,lj,lk)=.true.
   grid_ww_int(num,ui,uj,k)=.true.
   grid_ww_int(num,ui,lj,k)=.true.
   grid_ww_int(num,li,uj,k)=.true.
   grid_ww_int(num,li,lj,k)=.true.
   grid_ww_int(num,ui,j,uk)=.true.
   grid_ww_int(num,ui,j,lk)=.true.
   grid_ww_int(num,li,j,uk)=.true.
   grid_ww_int(num,li,j,lk)=.true.
   grid_ww_int(num,i,uj,uk)=.true.
   grid_ww_int(num,i,uj,lk)=.true.
   grid_ww_int(num,i,lj,uk)=.true.
   grid_ww_int(num,i,lj,lk)=.true.
   grid_ww_int(num,ui,j,k)=.true.
   grid_ww_int(num,li,j,k)=.true.
   grid_ww_int(num,i,j,uk)=.true.
   grid_ww_int(num,i,j,lk)=.true.
   grid_ww_int(num,i,uj,k)=.true.
   grid_ww_int(num,i,lj,k)=.true.

   grid_ww(num)%x=xtmp
   grid_ww(num)%y=ytmp
   grid_ww(num)%z=ztmp

   grid_ww(num)%xend=xtmp+gridRc
   grid_ww(num)%yend=ytmp+gridRc
   grid_ww(num)%zend=ztmp+gridRc

   num=num+1
   ztmp=ztmp+gridRc
  end do
  ytmp=ytmp+gridRc
 end do
 xtmp=xtmp+gridRc
end do

end if
end subroutine create_grid_ww
#endif

!-----------------------------------------------------------------------


subroutine close_input_files
close (1)
if(restart) close (2)
if ( implicit_rstr_from_file .eq. 1 ) close (12)
close (13)

end subroutine close_input_files

!-----------------------------------------------------------------------

subroutine close_output_files
close (3)
if ( itrj_cycle .gt. 0 ) close (10)
if ( iene_cycle .gt. 0 ) close (11)

end subroutine close_output_files

!-----------------------------------------------------------------------

subroutine open_files
! --> restart file (2)
if(restart) then
open (unit=2, file=restart_file, status='old', form='unformatted', action='read', err=2)
end if

! --> final coords (3)
open (unit=3, file=xfin_file, status='unknown', form='unformatted', action='write', err=3)

! --> energy output file (11)
if ( iene_cycle .gt. 0 ) then
open (unit=11, file=ene_file, status='unknown', form='unformatted', action='write', err=11)
end if

! --> external file for implicit position restraints (12)
if ( implicit_rstr_from_file .eq. 1 ) then
open (unit=12, file=exrstr_file, status='old', form='unformatted', action='read', err=12)
end if


return

! crude error handling
2 call die('error opening restart file.')
3 call die('error opening final coordinates file.')
11 call die('error opening energy output file.')
12 call die('error opening position restraints file.')

end subroutine open_files

!-----------------------------------------------------------------------

subroutine gauss (am,sd,v,ig)
! arguments
real(kind=prec)					::	am,sd,v
integer					::	ig

! local variables
integer					::	i
real(kind=prec)					::	a,y

a=zero
do i=1,12
y=randm(ig)
a=a+y
end do
v=(a-6.0_prec)*sd+am
end subroutine gauss

!-----------------------------------------------------------------------
subroutine get_fep
! local variables
character(len=200)			::	libtext
character(len=2)			::	qaname
integer					::	i,j,k,iat
!temp. array for reallocating long-range exclusion list
integer(AI), pointer	::	tempexlong(:,:)

! --- # states, # q-atoms
if(.not. qatom_load_atoms(fep_file)) then
        call die('failure to load Q-atoms from FEP file.')
end if

! set flags
do i=1,nqat
        if(iqseq(i) > 0 .and. iqseq(i) <= nat_solute)  then
                iqatom(iqseq(i)) = i
        else if(iqseq(i) == 0) then
                write(*,10) i
        else
                write(*,20) i, iqseq(i)
                call die('invalid q-atom data')
        end if
end do
10	format('>>> WARNING: Q-atom no. ',i2,' is not associated with a topology atom.')
20	format('>>>>> ERROR: Q-atom no. ',i2,' has invalid topology number ',i5)
!allocate memory for qatom charges
allocate(qcrg(nqat,nstates), stat=alloc_status)
call check_alloc('Qatom charges')

! --- copy topology charges

do i=1,nqat
        do j=1,nstates
                qcrg(i,j)=crg(iqseq(i))
        end do
end do

!initialize softcore lookup array
allocate (sc_lookup(nqat,natyps+nqat,nstates))
sc_lookup(:,:,:)=zero

!load rest of fep file
if(.not. qatom_load_fep(fep_file)) then
        call die('failure to load FEP file.')
end if

!Adapt LJ parameters to topology
!If arithmetic combination rule take sqrt(epsilon) now
if (qvdw_flag .and. ivdw_rule .eq. 2 ) then
        qbvdw(:,1) = sqrt( qbvdw(:,1) )
        qbvdw(:,3) = sqrt( qbvdw(:,3) )
end if

!remove redefined bonded interactions from topology
if(nqbond > 0 .or. nqangle > 0 .or. nqtor > 0 .or. nqimp > 0 ) then
        write(*,*)
        call centered_heading('Removing redefined interactions from topology','-')
230		format('type',t10,' atom1 atom2 atom3 atom4')
        write(*,230)
231		format(a,t10,4i6)
        !remove bonds that were redefined
        do i=1,nbonds
                do j=1,nqbond
                        if ( (bnd(i)%i==qbnd(j)%i .and. bnd(i)%j==qbnd(j)%j) .or. &
                                (bnd(i)%i==qbnd(j)%j .and. bnd(i)%j==qbnd(j)%i) ) then
                                bnd(i)%cod = 0
                                write (*,231) 'bond',bnd(i)%i,bnd(i)%j
                        end if
                end do
        end do

        !remove angles that were redefined
        do i=1,nangles
                do j=1,nqangle
                        if((ang(i)%i.eq.qang(j)%i .and. ang(i)%j.eq.qang(j)%j .and. &
                                ang(i)%k.eq.qang(j)%k)                          .or. &
                                (ang(i)%i.eq.qang(j)%k .and. ang(i)%j.eq.qang(j)%j .and. &
                                ang(i)%k.eq.qang(j)%i) )                         then

                                ang(i)%cod = 0
                                write (*,231) 'angle',ang(i)%i,ang(i)%j,ang(i)%k
                        end if
                end do
        end do

        !remove torsions that were redefined
        do i=1,ntors
                do j=1,nqtor
			if(( (tor(i)%i.eq.qtor(j)%i .and. tor(i)%j.eq.qtor(j)%j .and. &
				tor(i)%k.eq.qtor(j)%k .and. tor(i)%l.eq.qtor(j)%l) .or. &
				(tor(i)%i.eq.qtor(j)%l .and. tor(i)%j.eq.qtor(j)%k .and. &
				tor(i)%k.eq.qtor(j)%j .and. tor(i)%l.eq.qtor(j)%i )) .and. &
				tor(i)%cod .ne. 0) then
                                tor(i)%cod = 0
                                write (*,231) 'torsion', tor(i)%i,tor(i)%j,tor(i)%k,tor(i)%l
						end if
                end do
        end do


        !remove impropers that were redefined
        select case(ff_type)
        case(FF_CHARMM) !special code for CHARMM
                do i=1,nimps
                        do j=1,nqimp
				if ( ((	(imp(i)%i.eq.qimp(j)%i) .or. &
					(imp(i)%i.eq.qimp(j)%l) .or. &
					(imp(i)%l.eq.qimp(j)%i)) .and. &
				     (	(imp(i)%j.eq.qimp(j)%l) .or. &
					(imp(i)%j.eq.qimp(j)%i) .or. &
					(imp(i)%j.eq.qimp(j)%j) .or. &
					(imp(i)%j.eq.qimp(j)%k)) .and. &
				     (	(imp(i)%k.eq.qimp(j)%i) .or. &
					(imp(i)%k.eq.qimp(j)%j) .or. &
					(imp(i)%k.eq.qimp(j)%k) .or. &
					(imp(i)%k.eq.qimp(j)%l))) .and. &
				     (	imp(i)%cod .ne. 0)) then
                                        imp(i)%cod = 0
                                        write (*,231) &
                                        'improper',imp(i)%i,imp(i)%j,imp(i)%k,imp(i)%l
                                end if
                        end do
                end do

        case default
        do i=1,nimps
                do j=1,nqimp
                        if(( (imp(i)%i.eq.qimp(j)%i) .and. (imp(i)%j.eq.qimp(j)%j) .and. &
                                (imp(i)%k.eq.qimp(j)%k) .and. (imp(i)%l.eq.qimp(j)%l)) .or. &
                                ((imp(i)%i.eq.qimp(j)%k) .and. (imp(i)%j.eq.qimp(j)%j) .and. &
                                (imp(i)%k.eq.qimp(j)%i) .and. (imp(i)%l.eq.qimp(j)%l)) .or. &
                                ((imp(i)%i.eq.qimp(j)%k) .and. (imp(i)%j.eq.qimp(j)%j) .and. &
                                (imp(i)%k.eq.qimp(j)%l) .and. (imp(i)%l.eq.qimp(j)%i)) .or. &
                                ((imp(i)%i.eq.qimp(j)%l) .and. (imp(i)%j.eq.qimp(j)%j) .and. &
                                (imp(i)%k.eq.qimp(j)%i) .and. (imp(i)%l.eq.qimp(j)%k)) .or. &
                                ((imp(i)%i.eq.qimp(j)%l) .and. (imp(i)%j.eq.qimp(j)%j) .and. &
                                (imp(i)%k.eq.qimp(j)%k) .and. (imp(i)%l.eq.qimp(j)%i))) then
                                imp(i)%cod = 0
                                write(*,231)'improper',imp(i)%i,imp(i)%j,imp(i)%k,imp(i)%l
                        end if
                end do
        end do
        end select
end if

!check special exclusions
!modify exclusion lists to inclue special exclusions between Q and non-Q
if(nexspec > 0) then
        allocate(tempexlong(2,nexlong+nexspec))
        tempexlong(:, 1:nexlong) = listexlong(:, 1:nexlong)
        deallocate(listexlong)
        listexlong => tempexlong
end if

do k = 1, nexspec
        i = exspec(k)%i
        j = exspec(k)%j
        if(i < 1 .or. i > nat_pro .or. j < 1 .or. j > nat_pro) then
                write(*, 592) k, i, j
                call die('invalid special exclusion data')
        end if
        !if one or more non-Q-atoms modify exclusion lists
        if(iqatom(i)==0 .or. iqatom(j)==0) then
                !With non-Q-atoms involved only accept all or no states
                if(any(exspec(k)%flag(1:nstates))) then
                        if(.not. all(exspec(k)%flag(1:nstates))) then
                                write(*,594) k
                                call die('invalid special exclusion data')
                        else !exlcude in all states
                                if(abs(j-i) <= max_nbr_range) then
                                        if(i < j) then
                                                listex(j-i,i) = .true.
                                        else
                                                listex(i-j,j) = .true.
                                        end if
                                else
                                        nexlong = nexlong + 1
                                        listexlong(1, nexlong) = i
                                        listexlong(2, nexlong) = j
                                end if
                        end if
                end if
        end if
end do
592	format('>>>>> ERROR: Special exclusion pair ',i2,' (',i5,1x,i5,') is invalid')
594	format('>>>>> ERROR: Non-Q-atom special excl. pair ',i2,' must be on in all or no states')
end subroutine get_fep

!-----------------------------------------------------------------------

subroutine get_fname (text,length,filnam)
! arguments
character(len=200)			::	text,filnam
integer					::	length
! local variables


integer					::	i

length=200
do i=1,200
if ( text(i:i) .eq. ' ' ) then
length=i-1
goto 10
end if
end do
10 filnam(1:length)=text(1:length)

end subroutine get_fname

!-----------------------------------------------------------------------

#if defined (USE_MPI)
!Defines and allocates variables needed in the md-calculations
!The node initiation is written for AI = 4. If changes are made to any size in
! sizes.f90 the MPI-code must be changed accordingly. It is not dynamically
! implemented yet.
subroutine init_nodes
!
! initialise slave nodes, sending to slaves:
!
! variables:
!  natom,nwat,nsteps,use_LRF,NBcycle,crg_ow,crg_hw,Rcpp,Rcww,Rcpw,Rcq,xpcent
!  nat_solute,ncgp,ncgp_solute,ivdw_rule,iuse_switch_atom,el14_scale,n14long
!  nexlong,natyps,nljtyp,rexcl_o,nstates,nqat,qvdw_flag,nqlib,RcLRF,
!  use_PBC, qswitch, nmol, nat_pro
!
! arrays:
!  x,v,iqatom,ljcod,qconn,iwhich_cgp,lrf,excl,iac,crg,cgpatom,cgp,iaclib
!  list14,listex,list14long,listexlong,iqseq,qiac,qcrg,qavdw,qbvdw,EQ(:)%lambda,
!  boxlength, inv_boxl, boxcentre, sc_lookup 
!

integer, parameter			:: vars = 40    !increment this var when adding data to broadcast in batch 1
integer				   	:: blockcnt(vars), ftype(vars)
integer(kind=MPI_ADDRESS_KIND)			   	:: fdisp(vars)
integer					:: mpitype_batch,mpitype_batch2
integer					:: nat3,j,jj,ii
!temp for shake
integer,allocatable                     :: shakeconsttmp(:)
integer,allocatable                     :: shakebondi(:)
integer,allocatable                     :: shakebondj(:)
integer                                 :: shakenumconst
real(kind=prec),allocatable             :: shakedist(:)
real(kind=wp8), allocatable		:: temp_lambda(:)
integer, parameter                      ::maxint=2147483647
real(kind=wp8), parameter                        ::maxreal=1E35_prec
integer  :: MPI_AI_INTEGER, MPI_TINY_INTEGER, i_loop

!external MPI_Address
!external MPI_Bcast

!**********
!2002-11-28 
!MN-> This will work with new implementations of MPI standard >= 2
!The MPI library at PDC does not support these definitions when I tried to use them.
!Using these routines will allow a change made to the sizes in sizes.f90 to
! affect the mpi. Without them the variables below marked (AI) and (TINY) will have to 
! be changed manually.
!When using this part make sure the vars marked with comments (AI) and (TINY) are 
! changed to MPI_AI_INTEGER and MPI_TINY_INTEGER.

!external MPI_Type_Create_F90_Integer
!external MPI_SizeOf

!Define data types
! This is wrong, the 1:st param is "Precision, in decimal digits", not bits
!call MPI_Type_Create_F90_Integer((8*AI-1),MPI_AI_INTEGER,ierr)
!call MPI_Type_Create_F90_Integer((8*TINY-1),MPI_TINY_INTEGER,ierr)
!To check the size in bytes of the new types use
!call MPI_SizeOf(MPI_AI_INTEGER,size,ierr)
!call MPI_SizeOf(MPI_TINY_INTEGER,size,ierr)



if (nodeid .eq. 0) call centered_heading('Distributing data to slave nodes', '-')
!***************************

! --- mandatory data, first batch ---

if (nodeid .eq. 0) write (*,'(80a)') 'MD data, first batch'

! initialise custom MPI data type; default is an integer scalar
! The variables below are sorted in modular order.
! Make sure to add next variable to the right module and 
! add +1 to 'vars'.
blockcnt(:) = 1
ftype(:) = MPI_INTEGER

! run control constants: natom, nwat, nsteps, NBmethod, NBcycle
call MPI_Bcast(natom, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast natom')
call MPI_Bcast(nwat, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast nwat')
call MPI_Bcast(nsteps, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast nsteps')
call MPI_Bcast(use_LRF, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast use_LRF')
call MPI_Bcast(NBcycle, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast NBcycle')

!some more new stuff for the thermostats/integrators
!now done in parallel
call MPI_Bcast(friction,1,MPI_REAL8,0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast friction')
call MPI_Bcast(thermostat,1,MPI_INTEGER,0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast thermostat')
call MPI_Bcast(dt,1,MPI_REAL8,0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast dt')
call MPI_Bcast(dt2,1,MPI_REAL8,0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast dt2')
call MPI_Bcast(Temp0,1,MPI_REAL8,0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast Temp0')
call MPI_Bcast(iseed,1,MPI_INTEGER,0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast dt2')
call MPI_Bcast(Temp0,1,MPI_REAL8,0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast iseed')

! water parameters: chg_solv array and solv_atom (used by nonbond_?w)
! also data type and internal interaction array
ftype(1) = MPI_INTEGER4
ftype(2) = MPI_INTEGER4
ftype(3) = MPI_REAL8
ftype(4) = MPI_REAL8
ftype(5) = MPI_REAL8
blockcnt(1) = 1
blockcnt(2) = 1
blockcnt(3) = 1
blockcnt(4) = 1
blockcnt(5) = 1
fdisp(1) = 0
fdisp(2) = 4
fdisp(3) = 4+4
fdisp(4) = 4+4+8
fdisp(5) = 4+4+8+8
call MPI_Type_create_struct(5, blockcnt, fdisp, ftype, mpitype_batch_solv_int, ierr)
if (ierr .ne. 0) call die('init_nodes MPI_Type_create_struct solv_int')
call MPI_Type_commit(mpitype_batch_solv_int, ierr)
if (ierr .ne. 0) call die('init_nodes MPI_Type_commit solv_int')

if (nwat .gt. 0 ) then
call MPI_Bcast(solvent_type,1,MPI_INTEGER,0,MPI_COMM_WORLD, ierr)
if (ierr .ne. 0 ) call die('init_nodes/MPI_Bcast solvent_type')
call MPI_Bcast(solv_atom,1,MPI_INTEGER,0,MPI_COMM_WORLD, ierr)
if (ierr .ne. 0 ) call die('init_nodes/MPI_Bcast solv_atom')
call MPI_Bcast(num_solv_int,1,MPI_INTEGER,0,MPI_COMM_WORLD, ierr)
if (ierr .ne. 0 ) call die('init_nodes/MPI_Bcast num_solv_int')

if (nodeid .ne. 0 ) then
allocate(chg_solv(solv_atom),stat=alloc_status)
call check_alloc('chg_solv')
allocate(aLJ_solv(solv_atom,3),stat=alloc_status)
call check_alloc('aLJ_solv')
allocate(bLJ_solv(solv_atom,3),stat=alloc_status)
call check_alloc('bLJ_solv')
allocate(nonbnd_solv_int(num_solv_int),stat=alloc_status)
call check_alloc('nonbnd_solv_int')
end if
call MPI_Bcast(chg_solv,solv_atom,MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast chg_solv array')
call MPI_Bcast(aLJ_solv,3*solv_atom,MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast chg_solv array')
call MPI_Bcast(bLJ_solv,3*solv_atom,MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast chg_solv array')
call MPI_Bcast(nonbnd_solv_int,num_solv_int,mpitype_batch_solv_int, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast nonbnd_solv_int struct')
end if

! cutoffs: Rcpp, Rcww, Rcpw, Rcq, RcLRF (used by pair list generating functions)
call MPI_Bcast(Rcpp, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast Rcpp')
call MPI_Bcast(Rcww, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast Rcww')
call MPI_Bcast(Rcpw, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast Rcpw')
call MPI_Bcast(Rcq, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast Rcq')
call MPI_Bcast(RcLRF, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast RcLRF')




!Periodic Boudary Condition
call MPI_Bcast(use_PBC, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast use_PBC')
call MPI_Bcast(boxcentre, 3, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast boxcentre')
call MPI_Bcast(boxlength, 3, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast boxlength')
call MPI_Bcast(inv_boxl, 3, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast inv_boxl')
call MPI_Bcast(qswitch, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast qswitch')
call MPI_Bcast(constant_pressure, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast constant_pressure')
call MPI_Bcast(ivolume_cycle, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast ivolume_cycle')
call MPI_Bcast(rigid_box_centre, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast rigid_box_centre')
call MPI_Bcast(put_solvent_back_in_box, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast put_solvent_back_in_box')
call MPI_Bcast(put_solute_back_in_box, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast put_solute_back_in_box')

! xpcent            from TOPO, needed for listgeneration
call MPI_Bcast(xpcent, 3, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast xpcent')


! a bunch of vars from the TOPO module
call MPI_Bcast(nat_solute, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast nat_solute')
call MPI_Bcast(ncgp, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast ncgp')
call MPI_Bcast(ncgp_solute, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast ncgp_solute')
call MPI_Bcast(ivdw_rule, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast ivdw_rule')
call MPI_Bcast(iuse_switch_atom, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast iuse_switch_atom')
call MPI_Bcast(el14_scale, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast el14_scale')
call MPI_Bcast(n14long, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast n14long')
call MPI_Bcast(nexlong, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast nexlong')
call MPI_Bcast(natyps, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast natyps')
call MPI_Bcast(rexcl_o, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast rexcl')
call MPI_Bcast(nmol, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast nmol')
call MPI_Bcast(nres_solute, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast nres_solute')
call MPI_Bcast(nat_pro, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast nat_pro')

! water pol stuff
call MPI_Bcast(wpol_restr, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast wpol_restr')
call MPI_Bcast(n_max_insh, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast n_max_insh')
call MPI_Bcast(nwpolr_shell, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast nwpolr_shell')
call MPI_Bcast(fk_pshell, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast fk_pshell')
call MPI_Bcast(fk_wsphere, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast fk_wsphere')
call MPI_Bcast(Dwmz, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast Dwmz')
call MPI_Bcast(awmz, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast awmz')
call MPI_Bcast(fkwpol, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast fkwpol')
call MPI_Bcast(rwat, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast rwat')



#ifdef USE_GRID
!first few grid variables
call MPI_Bcast(gridstor_pp, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast gridstor_pp')
call MPI_Bcast(gridstor_pw, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast gridstor_pw')
call MPI_Bcast(gridstor_ww, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast gridstor_ww')

! water grid lookup
if (nwat .gt. 0 ) then
if (nodeid .ne. 0 ) then
allocate(ww_igrid(nwat),stat=alloc_status)
call check_alloc('ww_igrid slave nodes')
allocate(pw_igrid(ncgp_solute),stat=alloc_status)
call check_alloc('pw_igrid slave nodes')
end if
call MPI_Bcast(ww_igrid,nwat,MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast ww_igrid')
call MPI_Bcast(pw_igrid,ncgp_solute,MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast pw_igrid')
end if
! protein grid lookup
if (ncgp_solute .gt. 0) then
if (nodeid .ne. 0 ) then
allocate(pp_igrid(ncgp_solute),stat=alloc_status)
call check_alloc('pp_igrid slave nodes')
end if
call MPI_Bcast(pp_igrid,ncgp_solute,MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast pp_igrid')
end if ! ncgp_solute .gt. 0
#endif

!vars from QATOM
call MPI_Bcast(nstates, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast nstates')
call MPI_Bcast(use_excluded_groups,1,MPI_LOGICAL, 0,MPI_COMM_WORLD,ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast use_excluded_groups')

if (use_excluded_groups) then
	call MPI_Bcast(ngroups_gc, 1, MPI_INTEGER, 0,MPI_COMM_WORLD,ierr)
	if (ierr .ne. 0) call die('init_nodes/MPI_Bcast ngroups_gc')
end if
call MPI_Bcast(nqat, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast nqat')
call MPI_Bcast(qvdw_flag, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast qvdw_flag')
call MPI_Bcast(nqlib, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast nqlib')

if (use_excluded_groups) then
        allocate(tempmask(nat_pro,ngroups_gc),stat=alloc_status)
        call check_alloc('temp mask arrays')
        if (nodeid .ne. 0) then
                allocate(ST_gc(ngroups_gc),stat=alloc_status)
                call check_alloc('node gc array')
        end if
        if (nodeid.eq.0) then
                do jj=1,ngroups_gc
                tempmask(1:nat_pro,jj)=ST_gc(jj)%gcmask%mask(1:nat_pro)
                end do
        end if
! Wait first for Master here
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call MPI_Bcast(tempmask,nat_pro*ngroups_gc,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
        if (ierr .ne. 0) call die('init_nodes/MPI_Bcast gc atom masks')
        if (nodeid.ne.0) then
                do jj=1,ngroups_gc
                call mask_initialize(ST_gc(jj)%gcmask)
                ST_gc(jj)%gcmask%mask(1:nat_pro)=tempmask(1:nat_pro,jj)
                end do
        end if
        deallocate(tempmask,stat=alloc_status)
	call MPI_Bcast(ST_gc(1:ngroups_gc)%caltype,ngroups_gc,MPI_INTEGER, 0,MPI_COMM_WORLD,ierr)
	if (ierr .ne. 0) call die('init_nodes/MPI_Bcast group contrib calc types')
end if

call MPI_Bcast(xwcent,3,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast wxcent')

! --- MD data, second batch ---

if (nodeid .eq. 0) write (*,'(80a)') 'MD data, second batch'

! allocate arrays
if (nodeid .ne. 0) then
 call allocate_natom_arrays
 if (thermostat == NOSEHOOVER) then
  call allocate_nhchain_arrays
 end if
end if

! broadcast x, v and winv
nat3 = 3*natom

call MPI_Bcast(x, nat3, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast x')
call MPI_Bcast(v, nat3, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast v')


call MPI_Bcast(winv,natom,MPI_REAL8,0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast winv')

!Broadcast iqatom
call MPI_Bcast(iqatom, natom, MPI_INTEGER2, 0, MPI_COMM_WORLD, ierr) !(TINY)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast iqatom')

call MPI_Bcast(num_atyp,1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast num_atyp')

!Broadcast ljcod
if (nodeid .ne. 0) then
allocate(ljcod(num_atyp,num_atyp))
end if
call MPI_Bcast(ljcod, size(ljcod), MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast ljcod')

!Broadcast qconn(nstates,nat_solute, nqat)
if (nodeid .ne. 0) then
allocate(qconn(nstates,nat_solute, nqat),stat=alloc_status)
call check_alloc('qconn')
end if
call MPI_Bcast(qconn, size(qconn), MPI_INTEGER2, 0, MPI_COMM_WORLD, ierr) ! (TINY)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast qconn')


! --- Periodic boundary condition data ---


if (nodeid.ne.0) then
allocate(mvd_mol(nmol))
end if


! --- lrf data ---

if (use_LRF) then
! lrf stuff

if (nodeid .eq. 0) write (*,'(80a)') 'lrf data'

! allocate arrays
if (nodeid .ne. 0) call allocate_lrf_arrays

!MPI_INTEGER4 is used instead of MPI_AI_INTEGER
!Change to mpi_type_create, see note above or note2 in sizes.f90
! iwhich_cgp
call MPI_Bcast(iwhich_cgp, natom, MPI_INTEGER4, 0, MPI_COMM_WORLD, ierr) !(AI)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast lrf parameters')


! lrf
! we now keep the data type of the lrf for later usage
ftype(:) = MPI_REAL8
blockcnt(1) = 3					! real(kind=prec) cgp_cent(3)
fdisp(1) = 0
blockcnt(2) = 1					! real(kind=prec) phi0
fdisp(2) = 3*8
blockcnt(3) = 3					! real(kind=prec) phi1(3)
fdisp(3) = 3*8 + 8
blockcnt(4) = 9					! real(kind=prec) phi2(9)
fdisp(4) = 3*8 + 8 + 3*8
blockcnt(5) = 27				! real(kind=prec) phi3(27)
fdisp(5) = 3*8 + 8 + 3*8 + 9*8
call MPI_Type_create_struct(5, blockcnt, fdisp, ftype, mpitype_batch_lrf, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Type_create_struct')
call MPI_Type_commit(mpitype_batch_lrf, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Type_commit')
call MPI_Bcast(lrf, ncgp, mpitype_batch_lrf, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast lrf parameters')
call MPI_Op_create(lrf_add,.false.,mpi_lrf_add,ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Op_create mpi_lrf_add')
call MPI_Op_create(lrf_cgp_rep,.false.,mpi_lrf_cgp_rep,ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Op_create mpi_lrf_cgp_rep')

end if !(use_LRF)

! --- data from the TOPO module ---

if (nodeid .eq. 0) write (*,'(80a)') 'TOPO data'

! allocate topology arrays
if (nodeid .ne. 0) then
! don't allocate memory for stuff we don't need
! these array size variables are actually used
max_cgp=ncgp
max_atyps = natyps
max_14long = n14long
max_exlong = nexlong
max_atom = natom

call topo_allocate_atom(alloc_status)
call check_alloc('topology arrays')
call topo_allocate_potential(alloc_status)
call check_alloc('topology arrays')
allocate(istart_mol(nmol+1), &
stat=alloc_status)
call check_alloc('topology arrays')
end if

! broadcast excl
call MPI_Bcast(excl, natom, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast excl')
! broadcast istart_mol
call MPI_Bcast(istart_mol, nmol+1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast istart_mol')

! Bcast iac, crg and cgpatom 
call MPI_Bcast(iac, natom, MPI_INTEGER2, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast iac')
call MPI_Bcast(crg, natom, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast crg')
call MPI_Bcast(cgpatom, natom, MPI_INTEGER4, 0, MPI_COMM_WORLD, ierr) !(AI)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast cgpatom')

call MPI_Bcast(xtop,nat3,MPI_REAL8,0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast xtop')
call MPI_Bcast(shell,natom, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast shell')


! cgp
!Use MPI_Type_create_struct here too
ftype(:) = MPI_INTEGER4 !(AI)
blockcnt(:) = 1
fdisp(1) = 0				! integer(AI) iswitch
fdisp(2) = AI				! integer(AI) first
fdisp(3) = AI + AI			! integer(AI) last
call MPI_Type_create_struct(3, blockcnt, fdisp, ftype, mpitype_batch, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Type_create_struct')
call MPI_Type_commit(mpitype_batch, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Type_commit')
call MPI_Bcast(cgp, ncgp, mpitype_batch, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast cgp')
call MPI_Type_free(mpitype_batch, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Type_free')

! iaclib
ftype(:) = MPI_REAL8
blockcnt(1) = 1					! real(kind=prec) mass
fdisp(1) = 0
blockcnt(2) = nljtyp				! real(kind=prec) avdw(nljtyp)
fdisp(2) = 8
blockcnt(3) = nljtyp				! real(kind=prec) bvdw(nljtyp)
fdisp(3) = 8 + 8*nljtyp
call MPI_Type_create_struct(3, blockcnt, fdisp, ftype, mpitype_batch, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Type_create_struct')
call MPI_Type_commit(mpitype_batch, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Type_commit')
call MPI_Bcast(iaclib, max_atyps, mpitype_batch, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast iaclib')
call MPI_Type_free(mpitype_batch, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Type_free')

! list14 and listex share the same format: logical listxx(max_nbr_range,max_atom)
call MPI_Bcast(list14, size(list14), MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast list14')
call MPI_Bcast(listex, size(listex), MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast listex')

! list14long and listexlong share the same format: integer(AI) listxxlong(2,max_nxxlong)
call MPI_Bcast(list14long, 2*n14long, MPI_INTEGER4, 0, MPI_COMM_WORLD, ierr) !(AI)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast list14long')
call MPI_Bcast(listexlong, 2*nexlong, MPI_INTEGER4, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast listexlong')

#ifdef USE_GRID

! data for the different grids
call MPI_Bcast(pp_gridnum,1,MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast pp_gridnum')
call MPI_Bcast(pp_ndim,1,MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast pp_ndim')
call MPI_Bcast(pw_gridnum,1,MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast pw_gridnum')
call MPI_Bcast(pw_ndim,1,MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast pw_ndim')
call MPI_Bcast(ww_gridnum,1,MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast ww_gridnum')
call MPI_Bcast(ww_ndim,1,MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast ww_ndim')

if (nodeid.ne.0) then
!prepare grids on the different nodes using the data from master
if (ncgp_solute .gt. 0) then
allocate(grid_pp(pp_gridnum),stat=alloc_status)
call check_alloc('MPI pp_grid')
allocate(grid_pp_grp(pp_gridnum,gridstor_pp),stat=alloc_status)
call check_alloc('MPI grid_pp_grp')
allocate(grid_pp_ngrp(pp_gridnum),stat=alloc_status)
call check_alloc('MPI grid_pp_ngrp')
allocate(grid_pp_int(pp_gridnum,pp_ndim,pp_ndim,pp_ndim),stat=alloc_status)
call check_alloc('MPI grid_pp_int')
end if ! ncgp_solute .gt. 0
if (nwat .gt. 0) then
allocate(grid_pw(pw_gridnum),stat=alloc_status)
call check_alloc('MPI pw_grid')
allocate(grid_ww(ww_gridnum),stat=alloc_status)
call check_alloc('MPI ww_grid')
allocate(grid_pw_grp(pw_gridnum,gridstor_pw),stat=alloc_status)
call check_alloc('MPI grid_pw_grp')
allocate(grid_pw_ngrp(pw_gridnum),stat=alloc_status)
call check_alloc('MPI grid_pw_ngrp')
allocate(grid_pw_int(pw_gridnum,pw_ndim,pw_ndim,pw_ndim),stat=alloc_status)
call check_alloc('MPI grid_pw_int')
allocate(grid_ww_grp(ww_gridnum,gridstor_ww),stat=alloc_status)
call check_alloc('MPI grid_ww_grp')
allocate(grid_ww_ngrp(ww_gridnum),stat=alloc_status)
call check_alloc('MPI grid_ww_ngrp')
allocate(grid_ww_int(ww_gridnum,ww_ndim,ww_ndim,ww_ndim),stat=alloc_status)
call check_alloc('MPI grid_ww_int')
end if !nwat > 0
end if ! nodeid .ne. 0

!make first process wait until nodes have allocated grid arrays
call MPI_BARRIER(MPI_COMM_WORLD, ierr)

! also use make struct for grids
! loop based assignment not working
ftype(1) = MPI_REAL8            ! real(kind=prec)       x
ftype(2) = MPI_REAL8            ! real(kind=prec)       y
ftype(3) = MPI_REAL8            ! real(kind=prec)       z
ftype(4) = MPI_REAL8            ! real(kind=prec)       xend
ftype(5) = MPI_REAL8            ! real(kind=prec)       yend
ftype(6) = MPI_REAL8            ! real(kind=prec)       zend
blockcnt(1) = 1
blockcnt(2) = 1
blockcnt(3) = 1
blockcnt(4) = 1
blockcnt(5) = 1
blockcnt(6) = 1
fdisp(1) = 0     
fdisp(2) = 1*8
fdisp(3) = 2*8
fdisp(4) = 3*8
fdisp(5) = 4*8
fdisp(6) = 5*8
call MPI_Type_create_struct(6, blockcnt, fdisp, ftype, mpitype_batch, ierr)
if (ierr .ne. 0) call die('init_nodes MPI_Type_create_struct grid')
call MPI_Type_commit(mpitype_batch, ierr)
if (ierr .ne. 0) call die('init_nodes MPI_Type_commit grid')
if (ncgp_solute .gt. 0) then
call MPI_Bcast(grid_pp, pp_gridnum, mpitype_batch, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes MPI_Bcast grid_pp')
end if
if (nwat .gt. 0 ) then
call MPI_Bcast(grid_pw, pw_gridnum, mpitype_batch, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes MPI_Bcast grid_pw')
call MPI_Bcast(grid_pw, pw_gridnum, mpitype_batch, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes MPI_Bcast grid_ww')
end if ! nwat > 0
call MPI_Type_free(mpitype_batch, ierr)
if (ierr .ne. 0) call die('init_nodes MPI_Type_free grid')

! now the rest of the individual arrays
! only need to send interaction matrix, rest is not constant and updated in
! populate grids
if (ncgp_solute .gt. 0 ) then
call MPI_Bcast(grid_pp_int,pp_gridnum*pp_ndim**3,MPI_LOGICAL,0, MPI_COMM_WORLD,ierr)
if (ierr .ne. 0) call die('init_nodes MPI_Bcast grid_pp_int')
! make grid_pp_grp mpi data structure here
call MPI_Type_contiguous(gridstor_pp, MPI_INTEGER, mpitype_batch_ppgrid, ierr)
if (ierr .ne. 0) call die('init_nodes MPI_Type_contiguous gridstor_pp')

endif
if (nwat .gt. 0) then
call MPI_Bcast(grid_pw_int,pw_gridnum*pw_ndim**3,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
if (ierr .ne. 0) call die('init_nodes MPI_Bcast grid_pw_int')
call MPI_Bcast(grid_ww_int,ww_gridnum*ww_ndim**3,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
if (ierr .ne. 0) call die('init_nodes MPI_Bcast grid_ww_int')

call MPI_Type_contiguous(gridstor_pw, MPI_INTEGER, mpitype_batch_pwgrid, ierr)
if (ierr .ne. 0) call die('init_nodes MPI_Type_contiguous gridstor_pw')
call MPI_Type_contiguous(gridstor_ww, MPI_INTEGER, mpitype_batch_wwgrid, ierr)
if (ierr .ne. 0) call die('init_nodes MPI_Type_contiguous gridstor_ww')

endif

call MPI_Op_create(grid_add,.false.,mpi_grid_add,ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Op_create mpi_grid_add')


#endif

! --- data from the QATOM module ---

if (nodeid .eq. 0) write (*,'(80a)') 'QATOM data'

! allocate memory
if (nodeid .ne. 0) then
allocate(iqseq(nqat), &
qiac(nqat,nstates), &
qcrg(nqat,nstates), &
qavdw(nqlib,nljtyp), &
qbvdw(nqlib,nljtyp), &
EQ(nstates), &
sc_lookup(nqat,natyps+nqat,nstates), &
stat=alloc_status)
!if we use PBC, we allocate the storage for the old energies here
if( use_PBC ) then
allocate(old_EQ(nstates))
end if
call check_alloc('Q-atom arrays')
end if
!Broadcast sc_lookup(nqat,natyps+nqat,nstates)
if (nstates.ne.0) then
call MPI_Bcast(sc_lookup, size(sc_lookup), MPI_REAL8, 0, MPI_COMM_WORLD,ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast sc_lookup')
else
call MPI_Bcast(sc_lookup,  nstates, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast sc_lookup')
end if

! integer(AI) ::  iqseq(nqat)
!Change to mpi_type_create  (AI)
call MPI_Bcast(iqseq, nqat, MPI_INTEGER4, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast iqseq')

!  integer ::  qiac(nqat,nstates)
if (nstates.ne.0) then
call MPI_Bcast(qiac, size(qiac), MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast qiac')
else
call MPI_Bcast(qiac, nstates, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast qiac')
end if
! real(4) ::  qcrg(nqat,nstates)
if (nstates.ne.0) then
call MPI_Bcast(qcrg, size(qcrg), MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast qcrg')
else
call MPI_Bcast(qcrg, nstates, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast qcrg')
end if

if(qvdw_flag) then
!MN20030409-> Havn't tried with qvdw_flag == .true.
! qavdw and qbvdw share the same format: real(kind=prec) qxvdw(nqlib,nljtyp)
if (nstates.ne.0) then
call MPI_Bcast(qavdw, size(qavdw), MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast qavdw')
call MPI_Bcast(qbvdw, size(qbvdw), MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast qbvdw')
else
call MPI_Bcast(qavdw, nljtyp, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast qavdw')
call MPI_Bcast(qbvdw, nljtyp, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast qbvdw')
end if
end if
if (nstates .gt. 0) then
! Broadcast EQ(:)%lambda
allocate(temp_lambda(1:nstates), stat=alloc_status)
call check_alloc('Q-atom energy array')
if (nodeid .eq. 0) temp_lambda(1:nstates) = EQ(1:nstates)%lambda
call MPI_Bcast(temp_lambda, nstates, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast EQ%lambda')
if (nodeid .ne. 0) EQ(1:nstates)%lambda = temp_lambda(1:nstates)
deallocate(temp_lambda)
end if
!Final Bcast for the group contribution header data type
call MPI_Bcast(ene_header%arrays,1,MPI_INTEGER,0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast group contrib arrays')
call MPI_Bcast(ene_header%totresid,1,MPI_INTEGER,0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast group contrib totresid')

if (nodeid .ne. 0 ) then
allocate(ene_header%types(ene_header%arrays),ene_header%numres(ene_header%arrays),&
	ene_header%resid(ene_header%totresid),ene_header%gcnum(ene_header%arrays),&
	stat=alloc_status)
call check_alloc('Group contrib type arrays')
end if
call MPI_Bcast(ene_header%types,ene_header%arrays,MPI_INTEGER,0,MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast group contrib arrays')
call MPI_Bcast(ene_header%numres,ene_header%arrays,MPI_INTEGER,0,MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast group contrib arrays')
call MPI_Bcast(ene_header%gcnum,ene_header%arrays,MPI_INTEGER,0,MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast group contrib arrays')
call MPI_Bcast(ene_header%resid,ene_header%totresid,MPI_INTEGER,0,MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast group contrib arrays')

!Now make sure the nodes have all the EQ arrays allocated, too
if (nodeid .ne. 0) then
do jj=1,nstates
allocate(EQ(jj)%qx(ene_header%arrays),EQ(jj)%qq(ene_header%arrays),&
	EQ(jj)%qp(ene_header%arrays),EQ(jj)%qw(ene_header%arrays),&
	EQ(jj)%total(ene_header%arrays),stat=alloc_status)
!same as above, storage for old arrays is allocated here
if( use_PBC ) then
allocate(old_EQ(jj)%qx(ene_header%arrays),old_EQ(jj)%qq(ene_header%arrays),&
        old_EQ(jj)%qp(ene_header%arrays),old_EQ(jj)%qw(ene_header%arrays),&
        old_EQ(jj)%total(ene_header%arrays),stat=alloc_status)
end if
call check_alloc('Group contrib EQ arrays slave nodes')
end do
end if


!send the shake data here now, deleted old reference to this above
!shake now done on each node for own set of atoms
!PB 2015
if (nodeid .eq. 0) write (*,'(80a)') 'Shake control data'
call MPI_Bcast(shake_hydrogens,1,MPI_LOGICAL,0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast shake hydrogens')
call MPI_Bcast(shake_heavy,1,MPI_LOGICAL,0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast shake heavy')
call MPI_Bcast(shake_solvent,1,MPI_LOGICAL,0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast shake solvent')
call MPI_Bcast(shake_solute,1,MPI_LOGICAL,0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast shake solute')

call MPI_Bcast(Ndegf, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)    !bara i div init_
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast Ndegf')
call MPI_Bcast(Ndegfree, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)    !bara i div init_
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast Ndegfree')

call MPI_Bcast(shake_molecules,1,MPI_INTEGER,0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast shake molecules')
call MPI_Bcast(shake_constraints,1,MPI_INTEGER,0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast shake constraints')

!now we have to do some nasty stuff to get the shake data from master to the
!slaves, we can not just send the shake arrays, as they contain pointers and are
!thus not MPI transferable
!so we first make temp arrays for shake on master and send those to the nodes
if (shake_molecules .gt. 0) then
allocate(shakeconsttmp(shake_molecules), stat=alloc_status)
call check_alloc('Temporary shake arrays')
if (nodeid .eq. 0) then
shakenumconst = 0
do jj=1,shake_molecules
shakeconsttmp(jj)=shake_mol(jj)%nconstraints
shakenumconst = shakenumconst + shake_mol(jj)%nconstraints
end do
end if
if (nodeid .ne.0 ) then
allocate(shake_mol(nmol),stat = alloc_status)
call check_alloc('Slave node shake main array')
end if
call MPI_Bcast(shakenumconst,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast shakenumconst')
call MPI_Bcast(shakeconsttmp,shake_molecules,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast shakeconsttmp')

allocate(shakebondi(shakenumconst),shakebondj(shakenumconst),shakedist(shakenumconst),&
                stat=alloc_status)
call check_alloc('Temporary shake constraint arrays')
if (nodeid .eq. 0 ) then
shakenumconst = 0
do jj = 1,shake_molecules
        ii = 1
        do while (ii .le. shake_mol(jj)%nconstraints)
        shakenumconst = shakenumconst + 1
        shakebondi(shakenumconst) = shake_mol(jj)%bond(ii)%i
        shakebondj(shakenumconst) = shake_mol(jj)%bond(ii)%j
        shakedist(shakenumconst) = shake_mol(jj)%bond(ii)%dist2
        ii = ii +1
        end do
end do
end if
call MPI_Bcast(shakebondi,shakenumconst,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast shakebondi')
call MPI_Bcast(shakebondj,shakenumconst,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast shakebondj')
call MPI_Bcast(shakedist,shakenumconst,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast shakedist')

!yeah, now the slaves have all the shake information, need to put it now back
!into the data structure
if (nodeid .ne.0 ) then
do jj = 1, shake_molecules
shakenumconst = 0
shake_mol(jj)%nconstraints = shakeconsttmp(jj)
allocate(shake_mol(jj)%bond(shakeconsttmp(jj)),stat=alloc_status)
call check_alloc('slave node shkae bond array')
        ii = 1
        do while (ii .le. shake_mol(jj)%nconstraints)
        shakenumconst = shakenumconst + 1
        shake_mol(jj)%bond(ii)%i = shakebondi(shakenumconst)
        shake_mol(jj)%bond(ii)%j = shakebondj(shakenumconst)
        shake_mol(jj)%bond(ii)%dist2 = shakedist(shakenumconst)
        ii = ii + 1
        end do
end do
end if

deallocate(shakeconsttmp,shakebondi,shakebondj,shakedist)

end if ! shake_mol .gt. 0


!And the last -> length of bytes to recieve from the EQ arrays
call MPI_Bcast(reclength,1,MPI_INTEGER,0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast EQ array byte length')

if (nodeid .eq. 0) then 
call centered_heading('End of initiation', '-')
print *
end if

!and we sync before the end :)
call MPI_BARRIER(MPI_COMM_WORLD, ierr)
!Finally allocate for  slaves:E_send, EQ_send
!For master :E_recv,d_recv
call allocate_mpi  

end subroutine init_nodes
#endif


!-----------------------------------------------------------------------

subroutine init_shake
!
! initialize shake constraints
!
!locals
integer						::	mol, b, ia, ja, constr, angle, i
real(kind=prec)						:: exclshk
integer						::	src, trg
integer						::	solute_shake_constraints
logical                                         ::      shaken

!allocate molecule list
allocate(shake_mol(nmol), stat=alloc_status)
call check_alloc('shake molecule array')

shake_mol(:)%nconstraints = 0
mol = 0
exclshk = zero

!count bonds to be constrained in each molecule
!also count shake constraints involving excluded atoms
do b=1,nbonds
        ia = bnd(b)%i
        ja = bnd(b)%j
        do while(ia >= istart_mol(mol+1))
                !new molecule
                mol = mol +1
        end do
        !skip redefined bonds
        if(bnd(b)%cod == 0) cycle
!new shake flag to spec if people want all atoms, all solvent hydrogens, solute
!hydrogens, solute atoms ...
        shaken = .false.
        if(shake_solute .and. (ia.le.nat_solute)) then
                if(.not. heavy(ia) .or. .not. heavy(ja)) then
                        if(shake_hydrogens) then
                                shake_mol(mol)%nconstraints = shake_mol(mol)%nconstraints + 1
                                shaken = .true.
                        end if
                else
                        if(shake_heavy) then
                                shake_mol(mol)%nconstraints = shake_mol(mol)%nconstraints + 1
                                shaken = .true.
                        end if
                end if
        end if
        if(shake_solvent .and. (ia.gt.nat_solute)) then
                if(.not. heavy(ia) .or. .not. heavy(ja)) then
                        if(shake_hydrogens) then
                                shake_mol(mol)%nconstraints = shake_mol(mol)%nconstraints + 1
                                shaken = .true.
                        end if
                else
                        if(shake_heavy) then
                                shake_mol(mol)%nconstraints = shake_mol(mol)%nconstraints + 1
                                shaken = .true.
                        end if
                end if
        end if

           if(shaken .and. ( .not. use_PBC) ) then
                if(excl(ia)) exclshk = exclshk + 0.5_prec
                if(excl(ja)) exclshk = exclshk + 0.5_prec
           end if

end do
!count extra shake constraints from fep file in appropriate molecule
do b = 1, nqshake
    ia=iqshake(b)
        mol = 1
        do while(ia >= istart_mol(mol+1))
                mol = mol + 1
        end do
        shake_mol(mol)%nconstraints = shake_mol(mol)%nconstraints + 1
end do

!allocate bond lists for each molecule
do mol = 1, nmol
        !allocate(sbtemp(nconstr(mol), status = alloc_status)
        allocate(shake_mol(mol)%bond(shake_mol(mol)%nconstraints), stat = alloc_status)
        call check_alloc('shake bond array')
        !shake_mol(mol)%bonds => sbtemp
end do

mol = 0
!add the constraint
do b=1,nbonds
        shaken = .false.
        ia = bnd(b)%i
        ja = bnd(b)%j
        do while(ia >= istart_mol(mol+1)) 
                !new molecule
                mol = mol +1
                shake_mol(mol)%nconstraints = 0
        end do
        !skip redefined bonds
        if(bnd(b)%cod == 0) cycle
        if(shake_solute .and. (ia.le.nat_solute)) then
                if(.not. heavy(ia) .or. .not. heavy(ja)) then
                        if(shake_hydrogens) then
                                shaken = .true.
                        end if
                else
                        if(shake_heavy) then
                                shaken = .true.
                        end if
                end if
        end if
        if(shake_solvent .and. (ia.gt.nat_solute)) then
                if(.not. heavy(ia) .or. .not. heavy(ja)) then
                        if(shake_hydrogens) then
                                shaken = .true.
                        end if
                else
                        if(shake_heavy) then
                                shaken = .true.
                        end if
                end if
        end if


        if(shaken) then
                shake_mol(mol)%nconstraints = shake_mol(mol)%nconstraints + 1
                shake_mol(mol)%bond(shake_mol(mol)%nconstraints)%i = ia
                shake_mol(mol)%bond(shake_mol(mol)%nconstraints)%j = ja
                shake_mol(mol)%bond(shake_mol(mol)%nconstraints)%dist2 = &
                        bondlib(bnd(b)%cod)%bnd0**2
                !set the bond code to -1 for shaken bonds
                !bnd(b) will be deleted by shrink_topology
                bnd(b)%cod = -1
        end if
end do

!add extra shake constraints from fep file to appropriate molecule
do b = 1, nqshake
    ia=iqshake(b)
        ja=jqshake(b)
        mol = 1
        do while(ia >= istart_mol(mol+1))
                mol = mol + 1
        end do
        !see if already shaken
        do constr = 1, shake_mol(mol)%nconstraints
                if((ia == shake_mol(mol)%bond(constr)%i .and. &
                        ja == shake_mol(mol)%bond(constr)%j) .or. &
                   (ja == shake_mol(mol)%bond(constr)%i .and. &
                        ia == shake_mol(mol)%bond(constr)%j)) then
                        !found it: will overwrite 
                        !also decrement number of constraints
                        shake_mol(mol)%nconstraints = shake_mol(mol)%nconstraints - 1
                        exit
                end if
        end do
        !constr now contains the right index
        shake_mol(mol)%bond(constr)%i = ia
        shake_mol(mol)%bond(constr)%j = ja
        shake_mol(mol)%bond(constr)%dist2 = &
                dot_product(EQ(1:nstates)%lambda,qshake_dist(b,1:nstates))**2
        shake_mol(mol)%nconstraints = shake_mol(mol)%nconstraints + 1
end do

!get total number of shake constraints in solute (used for separate scaling of temperatures)
solute_shake_constraints = sum(shake_mol(1:nmol-nwat)%nconstraints)


!remove molecules with zero constraints from list
trg = 1
src = 2
do while(src <= nmol)
        if(shake_mol(trg)%nconstraints == 0) then
                shake_mol(trg) = shake_mol(src)
                !clear source
                shake_mol(src)%nconstraints = 0
                src = src + 1
        else
                trg = trg + 1
                if(trg == src) src = src + 1
        end if
end do
shake_molecules = trg

!total number of constraints
shake_constraints = sum(shake_mol(1:shake_molecules)%nconstraints)
write(*,100) shake_constraints
write(*,101) shake_molecules
100	format(/,'Number of shake constraints             = ',i10)
101	format('No. molecules with shake constraints    = ',i10)
! calculate #degrees of freedom
Ndegf=3*natom-shake_constraints    !changed from Ndegf=3*natom-3-shake_constraints, center of mass position is NOT CONstrained in the simulation, but IS constrained for initial temperatures....
Ndegfree=Ndegf-3*nexats+exclshk

Ndegf_solvent = Ndegf - 3*nat_solute + solute_shake_constraints
Ndegf_solute = Ndegf - Ndegf_solvent

Ndegfree_solvent = 3*(natom - nat_solute) - (shake_constraints - solute_shake_constraints)
Ndegfree_solute = Ndegfree - Ndegfree_solvent

if (Ndegfree_solvent*Ndegfree_solute .eq. 0) then    ! if either solvent or solute have 0 degrees of freedom, turn off separate scaling (in case it's on) and do not print detailed temperatures
	detail_temps = .false.
	separate_scaling = .false.
else
	detail_temps = .true.
end if


!clear angles which are shaken (i and k atoms shaken)
do mol=1, shake_molecules
        do constr = 1, shake_mol(mol)%nconstraints
        ia = shake_mol(mol)%bond(constr)%i
        ja = shake_mol(mol)%bond(constr)%j
                do angle = 1, nangles
                        if((ang(angle)%i == ia .and. ang(angle)%k == ja) .or. &
                           (ang(angle)%i == ja .and. ang(angle)%k == ia)) then
                                ang(angle)%cod = 0
                                exit
                        end if
                end do
        end do
end do
! we can not initilaize the nose-hoover variables before
! as they need the shake information
! so we set it after shake has been started
if ( thermostat == NOSEHOOVER ) then
        do i=1,numchain
                qnh(i)=nhq
                vnh(i)=0
                xnh(i)=0
        end do
        qnh(1)=Ndegf*nhq
end if

end subroutine init_shake

!-----------------------------------------------------------------------

subroutine maxwell
! *** local variables
integer						:: i,j,k
real(kind=prec)						:: sd,vg,kT

!	Generate Maxwellian velocities 
kT = Boltz*Tmaxw

do i=1,natom
        sd = sqrt (kT/iaclib(iac(i))%mass)
        do j=1,3
                call gauss (zero,sd,vg,iseed)
                k=(i-1)*3+j
                v(k)=vg
        end do
end do

end subroutine maxwell

!-----------------------------------------------------------------------
! Subroutine for the Nosé-Hoover chain propagation
subroutine nh_prop

real(kind=prec)				:: dt4, dt8, expf, s
integer				:: i, i3

!initializing all the constants to be used
dt4=0.5_prec*dt2
dt8=0.5_prec*dt4

Gnh(1) = Ndegf*Boltz*(Temp-Temp0)/qnh(1)

do i=2,numchain
	Gnh(i)=(qnh(i-1)*vnh(i-1)**2 - kbT)/qnh(i)
end do

vnh(numchain)=vnh(numchain)+Gnh(numchain)*dt4
xnh(numchain)=xnh(numchain)+vnh(numchain)*dt2

do i=numchain-1,1,-1
	expf=exp(-dt8*vnh(i+1))
	vnh(i)=(vnh(i)*expf+Gnh(i)*dt4)*expf
	xnh(i)=xnh(i)+vnh(i)*dt2
end do

s=exp(-vnh(1)*dt2)
do i=1,natom*3
	v(i)=v(i)*s
end do

Gnh(1) = Ndegf*Boltz*(s*s*Temp - Temp0)/qnh(1)
vnh(1) = (vnh(1)*expf+Gnh(1)*dt4)*expf

do i=2,numchain-1
        expf=exp(-dt8*vnh(i+1))
        Gnh(i)=(qnh(i-1)*vnh(i-1)**2 - kbT)/qnh(i)
        vnh(i)=(vnh(i)*expf+Gnh(i)*dt4)*expf
end do

Gnh(numchain)=(qnh(numchain-1)*vnh(numchain-1)**2 - kbT)/qnh(numchain)
vnh(numchain)=vnh(numchain)+Gnh(numchain)*dt4

end subroutine nh_prop

!-----------------------------------------------------------------------
subroutine prep_coord


! local variables
integer(4)          :: i,nat3
logical             :: old_restart = .false.
TYPE(qr_vec),allocatable :: x_old(:),v_old(:)
TYPE(qr_vec)        :: old_boxlength, old_boxcentre
integer             :: headercheck,myprec
!new variables for differen tprecision restart files
TYPE(qr_vecs),allocatable :: x_single(:),v_single(:)
TYPE(qr_vecs)             :: boxl_single,boxc_single
TYPE(qr_vecd),allocatable :: x_double(:),v_double(:)
TYPE(qr_vecd)             :: boxl_double,boxc_double
#ifndef PGI
TYPE(qr_vecq),allocatable   :: x_quad(:),v_quad(:)
TYPE(qr_vecq)               :: boxl_quad,boxc_quad
#endif
if (prec .eq. singleprecision) then
myprec = -137
elseif (prec .eq. doubleprecision) then
myprec = -1337
#ifndef PGI
elseif (prec .eq. quadprecision) then
myprec = -13337
#endif
else
call die('No such precision')
end if



! --- Refresh topology coords. if needed (external restraints file)
if ( implicit_rstr_from_file .eq. 1 ) then
write (*,'(/,a,/)') 'Refreshing topology coords for restraining...'
read(12) headercheck
  if ((headercheck .ne. -137).and.(headercheck.ne.-1337).and.(headercheck.ne.-13337)) then
!old restart file without header canary
      rewind(12)
      old_restart = .true.
  end if
      read (12) nat3
      rewind(12)
        if(nat3 /= 3*natom) then
          write(*,100) nat3/3, natom
          call die('wrong number of atoms in restart file')
       end if
  if (.not.old_restart) then
     read (12) headercheck
     if (myprec .ne. headercheck) then
       write(*,*) '>>> WARNING: Using mismatched precision in restart file'
       if (headercheck .eq. -137) then
         allocate(x_single(natom))
         read (12,err=112,end=112) nat3, (x_single(i),i=1,natom)
         xtop(1:natom) = x_single(1:natom)
         deallocate(x_single)
       else if (headercheck .eq. -1337) then
         allocate(x_double(natom))
         read (12,err=112,end=112) nat3, (x_double(i),i=1,natom)
         xtop(1:natom) = x_double(1:natom)
         deallocate(x_double)
       else if (headercheck .eq. -13337) then
#ifndef PGI
         allocate(x_quad(natom))
         read (12,err=112,end=112) nat3, (x_quad(i),i=1,natom)
         xtop(1:natom) = x_quad(1:natom)
         deallocate(x_quad)
#else
         call die('Quadruple precision not supported in PGI')
#endif
       end if
     else
       read (12,err=112,end=112) nat3, (xtop(i),i=1,natom)
     end if
  else
     allocate(x_old(natom))
     read (12,err=112,end=112) nat3, (x_old(i),i=1,natom)
     xtop(1:natom) = x_old(1:natom)
     deallocate(x_old)
  end if

end if
old_restart =.false.
headercheck=0
!Assign restraints of kind res:atom their numerical atom numbers
do i=1,nrstr_dist
  if(rstdis(i)%itext .ne. 'nil') then
    if (scan(rstdis(i)%itext,':') .ne. 0) then
      rstdis(i)%i=get_atom_from_resnum_atnum(rstdis(i)%itext)
	else
      read(rstdis(i)%itext,*) rstdis(i)%i
    end if
    if (scan(rstdis(i)%jtext,':') .ne. 0) then
      rstdis(i)%j=get_atom_from_resnum_atnum(rstdis(i)%jtext)
    else
      read(rstdis(i)%jtext,*) rstdis(i)%j
    end if
  end if
end do

! --- Make spherical restraining shell lists based on
!     the xtop coords.
if (.not. use_PBC) then

	if(rexcl_i > rexcl_o) then
	  call die('inner radius of restrained shell must be < exclusion radius')
	end if
	!first find atoms in shell 
	if (rexcl_i >= 0) then      !if rexcl_i is defined...
	  if (rexcl_i <= 1.00) then   !if rexcl_i is defined as fraction of rexcl_o
	    rexcl_i = rexcl_i * rexcl_o  !translate to Angstrom
	  end if 
	  if(iuse_switch_atom == 1) then
	     call make_shell
	  else
	     call make_shell2
	  end if
	else
	  write (*,'(/,a,/)') 'Restrained shell not defined!'
	end if
else
	shell(:) = .false.
end if ! .not. use_PBC
! --- read restart file

call allocate_natom_arrays
 if (thermostat == NOSEHOOVER) then
  call allocate_nhchain_arrays
 end if
if(restart) then
        ! topology routine has determined nwat, natom and allocated storage
        call centered_heading('Reading restart file','-')
        read(2) headercheck
        if ((headercheck .ne. -137).and.(headercheck.ne.-1337).and.(headercheck.ne.-13337)) then
!old restart file without header canary
          rewind(2)
          old_restart = .true.
        end if
          read (2) nat3
          rewind(2)
        if(nat3 /= 3*natom) then
                write(*,100) nat3/3, natom
100			format('>>>>> ERROR:',i5,' atoms in restart file not equal to',i5,&
                        ' in topology.')
                call die('wrong number of atoms in restart file')
        end if
        if (.not.old_restart) then
          read(2) headercheck
          if (myprec .ne. headercheck) then
            write(*,*) '>>> WARNING: Using mismatched precision in restart file'
            if (headercheck .eq. -137) then
              allocate(x_single(natom),v_single(natom))
              read (2,err=112,end=112) nat3, (x_single(i),i=1,natom)
              read (2,err=112,end=112) nat3, (v_single(i),i=1,natom)
              x(1:natom) = x_single(1:natom)
              v(1:natom) = v_single(1:natom)
              deallocate(x_single,v_single)
              if( use_PBC) then
                 read(2,err=112,end=112) boxl_single
                 read(2,err=112,end=112) boxc_single
                 boxlength = boxl_single
                 boxcentre = boxc_single
              end if
            else if (headercheck .eq. -1337) then
              allocate(x_double(natom),v_double(natom))
              read (2,err=112,end=112) nat3, (x_double(i),i=1,natom)
              read (2,err=112,end=112) nat3, (v_double(i),i=1,natom)
              x(1:natom) = x_double(1:natom)
              v(1:natom) = v_double(1:natom)
              deallocate(x_double,v_double)
              if( use_PBC) then
                 read(2,err=112,end=112) boxl_double
                 read(2,err=112,end=112) boxc_double
                 boxlength = boxl_double
                 boxcentre = boxc_double
              end if
            else if (headercheck .eq. -13337) then
#ifndef PGI
              allocate(x_quad(natom),v_quad(natom))
              read (2,err=112,end=112) nat3, (x_quad(i),i=1,natom)
              read (2,err=112,end=112) nat3, (v_quad(i),i=1,natom)
              x(1:natom) = x_quad(1:natom)
              v(1:natom) = v_quad(1:natom)
              deallocate(x_quad,v_quad)
              if( use_PBC) then
                 read(2,err=112,end=112) boxl_quad
                 read(2,err=112,end=112) boxc_quad
                 boxlength = boxl_quad
                 boxcentre = boxc_quad
              end if
#else
              call die('Quadruple precision not supported in PGI')
#endif
            end if
          else
              read (2,err=112,end=112) nat3, (x(i),i=1,natom)
              read (2,err=112,end=112) nat3, (v(i),i=1,natom)
              if( use_PBC) then
                read(2,err=112,end=112) boxlength
                read(2,err=112,end=112) boxcentre
              end if
          end if
        else
          allocate(x_old(natom),v_old(natom))
          read (2,err=112,end=112) nat3, (x_old(i),i=1,natom)
          read (2,err=112,end=112) nat3, (v_old(i),i=1,natom)
          write(*,*) 'Read coordinates and velocities from previous version of qdyn'
          x(1:natom) = x_old(1:natom)
          v(1:natom) = v_old(1:natom)
          deallocate(x_old,v_old)
          if( use_PBC) then
             read(2,err=112,end=112) old_boxlength
             read(2,err=112,end=112) old_boxcentre
             write(*,*) 'Read boxlength and center from previous version of qdyn'
             boxlength = old_boxlength
             boxcentre = old_boxcentre
          end if
        end if
        write (*,'(a30,i8)')   'Total number of atoms        =',natom
        write (*,'(a30,i8,/)') 'Number of waters encountered =',nwat

        if( use_PBC) then
                write(*,*)
                write(*,'(a16,3f8.3)') 'Boxlength     =', boxlength(:)
                write(*,'(a16,3f8.3)') 'Centre of box =', boxcentre(:)
        end if
        !water polarisation data will be read from restart file in wat_shells
else
        x(1:natom) = xtop(1:natom)
end if

! clear iqatom atom array
iqatom(:) = 0

return
#if defined(USE_MPI)
112 call MPI_Abort(MPI_COMM_WORLD,1,ierr)
#else
112 stop 'Aborting due to errors reading restart file.'
#endif

end subroutine prep_coord

!-----------------------------------------------------------------------

!Sort out heavy atoms in restrained shell. Use protein center to calculate distance.
!Uses coordinates from topology unless 'implicit_rstr_from_file' is specified.
subroutine make_shell
! *** Local variables
	integer						::	i,ig,i3
	real(kind=prec)						::	rin2,r2

	nshellats = 0
	rin2  = rexcl_i**2

	shell(:) = .false.

	do ig=1,ncgp_solute
       if (.not. excl(cgp(ig)%iswitch) .and. heavy(cgp(ig)%iswitch)) then 
		i3 = 3*cgp(ig)%iswitch-3
		r2 = ( xtop(i3+1) - xpcent(1) )**2 &
			+( xtop(i3+2) - xpcent(2) )**2 &
			+( xtop(i3+3) - xpcent(3) )**2

		if(r2 > rin2) then
			do i=cgp(ig)%first, cgp(ig)%last
        nshellats = nshellats + 1 
				shell(cgpatom(i)) = .true.
			end do
		end if
   end if
	end do
	write(*,105) nshellats, rexcl_i, rexcl_o
105	format('Found   ',i6,' solute atoms in the restrained shell region (',f6.2,' to ',f6.2,')')
end subroutine make_shell

!------------------------------------------------------------------------

!Sort out heavy atoms in restrained shell. Use protein center to calculate distance.
!Use coordinates from topology unless 'implicit_rstr_from_file' is specified
subroutine make_shell2
! *** Local variables
	integer						::	i,ig,i3,k
	real(kind=prec)						::	rout2,rin2,r2
	real(kind=prec), allocatable		::	cgp_cent(:,:)
	nshellats = 0
	rin2  = rexcl_i**2

	shell(:) = .false.

	allocate(cgp_cent(3,ncgp+nwat))

	cgp_cent(:,:) = zero

	do ig=1,ncgp_solute
    if (.not. excl(cgp(ig)%iswitch) .and. heavy(cgp(ig)%iswitch)) then
                do i = cgp(ig)%first,cgp(ig)%last
			i3 = cgpatom(i)*3
			cgp_cent(:,ig) = cgp_cent(:,ig) + xtop(i3-2:i3)
		end do
        cgp_cent(:,ig) = cgp_cent(:,ig)/real(cgp(ig)%last - cgp(ig)%first +1)
		r2 = dot_product(cgp_cent(:,ig)-xpcent(:),cgp_cent(:,ig)-xpcent(:))

		if ( r2 .gt. rin2 ) then
			do i=cgp(ig)%first, cgp(ig)%last
				nshellats = nshellats + 1
				shell(cgpatom(i)) = .true.
			end do
		end if
   end if
	end do

	deallocate(cgp_cent) 
	write(*,105) nshellats, rexcl_i, rexcl_o
105	format('Found   ',i6,' solute atoms in the restrained shell region (',f6.2,' to ',f6.2,'Å)')
end subroutine make_shell2

!-----------------------------------------------------------------------

subroutine init_trj
!locals
integer						::	trj_atoms

!initialise trajectory atom mask
if(itrj_cycle > 0) then
        call trj_initialize(frames=nsteps/itrj_cycle, steps_before=itrj_cycle,&
                interval=itrj_cycle, steps=nsteps,	degf=ndegfree, &
                topfile=top_file)

        trj_atoms = trj_commit_mask()
        write(*,100) trj_atoms
        if(.not. trj_create(trj_file)) then
                call die('failure to open trajectory file')
        end if
end if

100	format('Coordinates for',i6,' atoms will be written to the trajectory.')
end subroutine init_trj

!-----------------------------------------------------------------------
subroutine precompute_interactions
! nice and small routine to call all the other stuff below
if (nat_solute .ne. 0) call pp_int_comp
if ((nwat .gt. 0) .and. (nat_solute .ne.0)) call pw_int_comp
if ((nat_solute .ne. 0) .and. (nqat .ne.0)) call qp_int_comp
if (nqat .ne. 0) call qq_int_comp
if ((nwat .gt. 0) .and. (nqat.ne.0)) call qw_int_comp
if (nwat .gt. 0) call ww_int_comp


end subroutine precompute_interactions

!-----------------------------------------------------------------------

subroutine pp_int_comp
! here comes the real stuff
! this one takes info from the list update routines to 
! make the full list of all possible interactions and have them ready for
! immidiate lookup

! locals
integer                         :: ig,jg,nl
integer                         :: ia,ja,i,j
allocate(pp_precomp(nat_solute,nat_solute),stat=alloc_status)
call check_alloc('Protein-Protein precomputation array')

! 
pp_precomp(:,:)%set   = .false.
pp_precomp(:,:)%score = zero
pp_precomp(:,:)%elec  = zero
pp_precomp(:,:)%vdWA  = zero
pp_precomp(:,:)%vdWB  = zero

igloop: do ig = calculation_assignment%pp%start, calculation_assignment%pp%end
        ia = cgp(ig)%iswitch
        if ( excl(ia) ) cycle igloop

jgloop: do jg = 1, ncgp_solute
                ja = cgp(jg)%iswitch
                if ( excl(ja) ) cycle jgloop
! count each charge group pair once only
                if ( ((ig .gt. jg) .and. (mod(ig+jg,2) .eq. 0)) .or. &
                        ((ig .lt. jg) .and. (mod(ig+jg,2) .eq. 1)) ) &
                        cycle jgloop
ialoop:         do ia = cgp(ig)%first, cgp(ig)%last
                        i = cgpatom(ia)
!             --- q-atom ? ---
                        if ( iqatom(i) .ne. 0 ) cycle ialoop

jaloop:                 do ja = cgp(jg)%first, cgp(jg)%last
                                j = cgpatom(ja)
!             --- q-atom ? ---
                                if ( iqatom(j).ne.0 ) cycle jaloop
! count once
                                if ( ig .eq. jg .and. i .ge. j ) cycle jaloop

                                if ( abs(j-i) .le. max_nbr_range ) then
                                        if ( i .lt. j ) then
                                                if ( listex(j-i,i) ) then
                                                        cycle jaloop
                                                else if (list14(j-i,i)) then
                                                        call precompute_set_values_pp(i,j,3)
                                                        cycle jaloop
                                                end if
                                        else
                                                if ( listex(i-j,j) ) then
                                                        cycle jaloop
                                                else if ( list14(i-j,j)) then
                                                        call precompute_set_values_pp(i,j,3)
                                                        cycle jaloop
                                                end if
                                        end if
                                else
                                        do nl = 1, nexlong
                                                if ( (listexlong(1,nl) .eq. i .and. &
                                                        listexlong(2,nl) .eq. j      ) .or. &
                                                        (listexlong(1,nl) .eq. j .and. &
                                                        listexlong(2,nl) .eq. i      ) ) then
                                                        cycle jaloop
                                                end if
                                        end do
                                        do nl = 1, n14long
                                                if ( (list14long(1,nl) .eq. i .and. &
                                                        list14long(2,nl) .eq. j      ) .or. &
                                                        (list14long(1,nl) .eq. j .and. &
                                                        list14long(2,nl) .eq. i      ) ) then
                                                        call precompute_set_values_pp(i,j,3)
                                                        cycle jaloop
                                                end if
                                        end do
                                end if
                                call precompute_set_values_pp(i,j,ljcod(iac(i),iac(j)))
                        end do jaloop
                end do ialoop
        end do jgloop
end do igloop
end subroutine pp_int_comp

!-----------------------------------------------------------------------

subroutine pw_int_comp
! locals
integer                         :: ig,jg,ia,a_ind,b_ind

allocate(pw_precomp(nat_solute,solv_atom),stat=alloc_status)
call check_alloc('Protein-Solvent precomputation array')

pw_precomp(:,:)%set = .false.
pw_precomp(:,:)%score = zero
pw_precomp(:,:)%elec  = zero
pw_precomp(:,:)%vdWA  = zero
pw_precomp(:,:)%vdWB  = zero

igloop:do ig = calculation_assignment%pw%start, calculation_assignment%pw%end
jgloop: do jg = 1, solv_atom
ialoop:         do ia = cgp(ig)%first, cgp(ig)%last
                        a_ind = cgpatom(ia)
                        b_ind = jg
                        if ( iqatom(a_ind) .ne. 0 ) cycle ialoop
                        call precompute_set_values_pw(a_ind,b_ind,ljcod(iac(a_ind),iac(nat_solute+b_ind)))
                end do ialoop
        end do jgloop
end do igloop
end subroutine pw_int_comp

!-----------------------------------------------------------------------

subroutine qp_int_comp
! locals
integer                         :: ig,jq,ia,a_ind,is,vdw

allocate(qp_precomp(nat_solute,nqat,nstates),stat=alloc_status)
call check_alloc('Protein-QAtom precomputation array')

qp_precomp(:,:,:)%set   = .false.
qp_precomp(:,:,:)%score = zero
qp_precomp(:,:,:)%elec  = zero
qp_precomp(:,:,:)%vdWA  = zero
qp_precomp(:,:,:)%vdWB  = zero

igloop:do ig = 1, nat_solute
        if(any(qconn(:,ig,:) .le. 3)) cycle igloop
        do jq = 1, nqat
                vdw = ljcod(iac(ig),iac(iqseq(jq)))
                do is = 1 ,nstates
                        if(qconn(is, ig, jq) .eq. 4) vdw = 3
                        call precompute_set_values_qp(jq,ig,is,vdw)
                end do
        end do
end do igloop
end subroutine qp_int_comp

!-----------------------------------------------------------------------

subroutine qq_int_comp
! locals
integer                         :: iq,jq,ia,ja,k,l,vdw,i,is
real(kind=prec)                 :: tmp_elscale
logical                         :: found
allocate(qq_precomp(nqat,nqat,nstates),stat=alloc_status)
call check_alloc('QAtom-QAtom precomputation array')

qq_precomp(:,:,:)%set   = .false.
qq_precomp(:,:,:)%soft  = .false.
qq_precomp(:,:,:)%score = zero
qq_precomp(:,:,:)%elec  = zero
qq_precomp(:,:,:)%vdWA  = zero
qq_precomp(:,:,:)%vdWB  = zero

do iq = 1, nqat - 1
        ia = iqseq(iq)
        do jq = iq + 1, nqat
                ja = iqseq(jq)
                do is = 1, nstates
                        if(qconn(is, ja, iq) .ge. 4) then
                                tmp_elscale = one
                                if (nel_scale .ne. 0) then
                                        i = 1
                                        found = .false.
                                        do while(( i .le. nel_scale) .and. &
                                                (.not.found))
                                                k=qq_el_scale(i)%iqat
                                                l=qq_el_scale(i)%jqat
                                                if ((iq .eq. k .and. jq .eq. l) .or. &
                                                        (iq .eq. l .and. jq .eq. k)) then
                                                        tmp_elscale = qq_el_scale(i)%el_scale(is)
                                                        found = .true.
                                                        ! from Masoud
                                                end if
                                                i = i + 1
                                        end do
                                end if
                                if(qconn(is, ja, iq) .eq. 4) then
                                        vdw = 3
                                elseif(.not. qvdw_flag) then
                                        vdw = ljcod(iac(ia),iac(ja))
                                else
                                        vdw = 1
                                        i = 1
                                        found = .false.
                                        do while(( i .le. nqexpnb) .and. &
                                                (.not.found))
                                                if ((iq .eq. iqexpnb(i) .and.  &
                                                        jq .eq. jqexpnb(i))  .or. &
                                                        ( jq .eq. iqexpnb(i) .and. &
                                                        iq .eq. jqexpnb(i)))  then
                                                        vdw  = 2
                                                        found = .true.
                                                end if
                                                i = i + 1
                                        end do
                                end if ! (qconn = 4)
                                call precompute_set_values_qq(iq,jq,is,vdw,tmp_elscale)
                        end if
                end do
        end do
end do
! now we put those that were before also on q-q, but are actually q-p, on the
! q-p list where they belong!
do ja = 1, nat_solute
        if(iqatom(ja) .ne. 0) cycle
        if(any(qconn(:,ja,:) .le. 3)) then
                !bonded or angled to at least one Q-atom
                do iq = 1, nqat
                        do is = 1, nstates
                                if(qconn(is, ja, iq) .ge. 4) then
                                        if(qconn(is, ja, iq) .eq. 4) then
                                                vdw = 3
                                        elseif(qvdw_flag) then
                                                vdw = 1
                                        else
                                                vdw = ljcod(iac(ia),iac(ja))
                                        end if
                                        call precompute_set_values_qp(iq,ja,is,vdw)
                                end if
                        end do
                end do
        end if
end do
! prepare q-atom nonbond lists that do not need updating
#ifdef USE_MPI
if(nodeid .eq.0) then
#endif
call nbqqlist
#ifdef USE_MPI
end if
#endif
end subroutine qq_int_comp

!-----------------------------------------------------------------------

subroutine qw_int_comp
! locals
integer                         :: iq,jg,ia,a_ind,b_ind

allocate(qw_precomp(nqat,solv_atom,nstates),stat=alloc_status)
call check_alloc('QAtom-Solvent precomputation array')

qw_precomp(:,:,:)%set = .false.
qw_precomp(:,:,:)%score = zero
qw_precomp(:,:,:)%elec  = zero
qw_precomp(:,:,:)%vdWA  = zero
qw_precomp(:,:,:)%vdWB  = zero

igloop:do iq = 1, nqat
jgloop: do jg = 1, solv_atom
                a_ind = iq
                b_ind = jg
                call precompute_set_values_qw(a_ind,b_ind,ljcod(iac(nat_solute+b_ind),iac(iqseq(iq))))
        end do jgloop
end do igloop
end subroutine qw_int_comp

!-----------------------------------------------------------------------

subroutine ww_int_comp
! locals
integer                         :: ig,jg,ia,a_ind,b_ind

allocate(ww_precomp(solv_atom,solv_atom),stat=alloc_status)
call check_alloc('Solvent-Solvent precomputation array')

ww_precomp(:,:)%set = .false.
ww_precomp(:,:)%score = zero
ww_precomp(:,:)%elec  = zero
ww_precomp(:,:)%vdWA  = zero
ww_precomp(:,:)%vdWB  = zero

igloop:do ig = 1, solv_atom
jgloop: do jg = 1, solv_atom
                a_ind = ig
                b_ind = jg
                call precompute_set_values_ww(a_ind,b_ind,ljcod(iac(nat_solute+a_ind),iac(nat_solute+b_ind)))
        end do jgloop
end do igloop

end subroutine ww_int_comp

!-----------------------------------------------------------------------

subroutine precompute_set_values_pp(i,j,vdw)
! arguments
integer                         :: i,j,vdw
! locals
real(kind=prec)                 :: tempA,tempB
select case(ivdw_rule)
case(VDW_GEOMETRIC)
tempA = iaclib(iac(i))%avdw(vdw) * iaclib(iac(j))%avdw(vdw)
tempB = iaclib(iac(i))%bvdw(vdw) * iaclib(iac(j))%bvdw(vdw)
pp_precomp(i,j)%vdWA = tempA
pp_precomp(i,j)%vdWB = tempB

case(VDW_ARITHMETIC)
tempA = iaclib(iac(i))%avdw(vdw) + iaclib(iac(j))%avdw(vdw)
tempA = tempA**2
tempA = tempA * tempA * tempA
tempB = iaclib(iac(i))%bvdw(vdw) * iaclib(iac(j))%bvdw(vdw)

pp_precomp(i,j)%vdWA = (tempA**2) * tempB
pp_precomp(i,j)%vdWB = 2.0_prec * tempA * tempB

end select
pp_precomp(i,j)%elec = crg(i) * crg(j)

if ( (pp_precomp(i,j)%vdWA.ne.zero) .or. (pp_precomp(i,j)%vdWB.ne.zero) &
        .or. (pp_precomp(i,j)%elec.ne.zero)) pp_precomp(i,j)%set  = .true.

if (vdw .eq. 3 ) pp_precomp(i,j)%elec = pp_precomp(i,j)%elec * el14_scale

! make sure both possible ways to get there are prepared
pp_precomp(j,i) = pp_precomp(i,j)

end subroutine precompute_set_values_pp

!-----------------------------------------------------------------------

subroutine precompute_set_values_pw(i,j,vdw)
! arguments
integer                         :: i,j,vdw
! locals
real(kind=prec)                 :: tempA,tempB
select case(ivdw_rule)
case(VDW_GEOMETRIC)
tempA = iaclib(iac(i))%avdw(vdw) * aLJ_solv(j,vdw)
tempB = iaclib(iac(i))%bvdw(vdw) * bLJ_solv(j,vdw)
pw_precomp(i,j)%vdWA = tempA
pw_precomp(i,j)%vdWB = tempB

case(VDW_ARITHMETIC)
tempA = iaclib(iac(i))%avdw(vdw) + aLJ_solv(j,vdw)
tempA = tempA**2
tempA = tempA * tempA * tempA
tempB = iaclib(iac(i))%bvdw(vdw) * bLJ_solv(j,vdw)

pw_precomp(i,j)%vdWA = (tempA**2) * tempB
pw_precomp(i,j)%vdWB = 2.0_prec * tempA * tempB

end select
pw_precomp(i,j)%elec = crg(i) * chg_solv(j)

if ( (pw_precomp(i,j)%vdWA.ne.zero) .or. (pw_precomp(i,j)%vdWB.ne.zero) &
        .or. (pw_precomp(i,j)%elec.ne.zero)) pw_precomp(i,j)%set  = .true.

end subroutine precompute_set_values_pw

!-----------------------------------------------------------------------

subroutine precompute_set_values_qp(iq,j,istate,vdw)
! arguments
integer                         :: iq,i,j,vdw,qvdw,istate
! locals
real(kind=prec)                 :: tempA,tempB

i = iqseq(iq)
if (vdw .eq. 2 ) then
qvdw = 1
else
qvdw = vdw
end if
! for reference, we put all (!!!!!!!!!einself) qp interactions here
! even those listed for some reason under nbqqlist
! now nbqq only (!!!einself) has the q-q interactions

select case(ivdw_rule)
case(VDW_GEOMETRIC)
if (qvdw_flag) then
tempA = qavdw(qiac(iq,istate),qvdw) * iaclib(iac(j))%avdw(vdw)
tempB = qbvdw(qiac(iq,istate),qvdw) * iaclib(iac(j))%bvdw(vdw)
else
tempA = iaclib(iac(i))%avdw(vdw) * iaclib(iac(j))%avdw(vdw)
tempB = iaclib(iac(i))%bvdw(vdw) * iaclib(iac(j))%bvdw(vdw)
end if
qp_precomp(j,iq,istate)%vdWA = tempA
qp_precomp(j,iq,istate)%vdWB = tempB
case(VDW_ARITHMETIC)
if (qvdw_flag) then
tempA = qavdw(qiac(iq,istate),qvdw) + iaclib(iac(j))%avdw(vdw)
tempB = qbvdw(qiac(iq,istate),qvdw) * iaclib(iac(j))%bvdw(vdw)
else
tempA = iaclib(iac(i))%avdw(vdw) + iaclib(iac(j))%avdw(vdw)
tempB = iaclib(iac(i))%bvdw(vdw) * iaclib(iac(j))%bvdw(vdw)
end if
tempA = tempA**2
tempA = tempA * tempA * tempA
qp_precomp(j,iq,istate)%vdWA = (tempA**2) * tempB
qp_precomp(j,iq,istate)%vdWB = 2.0_prec * tempA * tempB
end select
if(.not. qq_use_library_charges) then
qp_precomp(j,iq,istate)%elec = qcrg(iq,istate) * crg(j)
else
qp_precomp(j,iq,istate)%elec = crg(i) * crg(j)
end if

if ( (qp_precomp(j,iq,istate)%vdWA.ne.zero) .or. (qp_precomp(j,iq,istate)%vdWB.ne.zero) &
        .or. (qp_precomp(j,iq,istate)%elec.ne.zero)) qp_precomp(j,iq,istate)%set  = .true.

qp_precomp(j,iq,istate)%score = sc_lookup(iq,iac(j),istate)
if (vdw .eq. 3 ) qp_precomp(j,iq,istate)%elec = qp_precomp(j,iq,istate)%elec * el14_scale

end subroutine precompute_set_values_qp

!-----------------------------------------------------------------------

subroutine precompute_set_values_qq(iq,jq,istate,vdw,q_elscale)
! arguments
integer                         :: iq,jq,i,j,vdw,istate
real(kind=prec)                 :: q_elscale ! from Masoud
! locals
real(kind=prec)                 :: tempA,tempB

i = iqseq(iq)
j = iqseq(jq)

! this one needs special treatment, as the nbqqlist routine makes a list of both
! q-q and q-p interactions for some reason (because of bonded/angled
! interactions to at least one q atom)
! so we just split it because else this will lead to branch points later on

select case(ivdw_rule)
case(VDW_GEOMETRIC)
if (qvdw_flag) then
tempA = qavdw(qiac(iq,istate),vdw) * qavdw(qiac(jq,istate),vdw)
tempB = qbvdw(qiac(iq,istate),vdw) * qbvdw(qiac(jq,istate),vdw)
else
tempA = iaclib(iac(i))%avdw(vdw) * iaclib(iac(j))%avdw(vdw)
tempB = iaclib(iac(i))%bvdw(vdw) * iaclib(iac(j))%bvdw(vdw)
end if
qq_precomp(iq,jq,istate)%vdWA = tempA
qq_precomp(iq,jq,istate)%vdWB = tempB
case(VDW_ARITHMETIC)
if (qvdw_flag) then
! soft pair case needs special treatment !
if (vdw .eq. 2 ) then
tempA = qavdw(qiac(iq,istate),vdw) * qavdw(qiac(jq,istate),vdw)
else
tempA = qavdw(qiac(iq,istate),vdw) + qavdw(qiac(jq,istate),vdw)
end if
tempB = qbvdw(qiac(iq,istate),vdw) * qbvdw(qiac(jq,istate),vdw)
else
tempA = iaclib(iac(i))%avdw(vdw) + iaclib(iac(j))%avdw(vdw)
tempB = iaclib(iac(i))%bvdw(vdw) * iaclib(iac(j))%bvdw(vdw)
end if
! soft pair case needs special treatment !
if ((vdw .eq. 2 ) .and. (qvdw_flag)) then
qq_precomp(iq,jq,istate)%vdWA = tempA
qq_precomp(iq,jq,istate)%vdWB = tempB
else
tempA = tempA**2
tempA = tempA * tempA * tempA
qq_precomp(iq,jq,istate)%vdWA = (tempA**2) * tempB
qq_precomp(iq,jq,istate)%vdWB = 2.0_prec * tempA * tempB
end if
end select
if(.not. qq_use_library_charges) then
qq_precomp(iq,jq,istate)%elec = qcrg(iq,istate) * qcrg(jq,istate)
else
qq_precomp(iq,jq,istate)%elec = crg(i) * crg(j)
end if
qq_precomp(iq,jq,istate)%elec  = qq_precomp(iq,jq,istate)%elec * q_elscale
qq_precomp(iq,jq,istate)%score = sc_lookup(iq,natyps+jq,istate)

if ( (qq_precomp(iq,jq,istate)%vdWA.ne.zero) .or. (qq_precomp(iq,jq,istate)%vdWB.ne.zero) &
        .or. (qq_precomp(iq,jq,istate)%elec.ne.zero)) qq_precomp(iq,jq,istate)%set = .true.

if (vdw .eq. 3 ) qq_precomp(iq,jq,istate)%elec = qq_precomp(iq,jq,istate)%elec * el14_scale
if ((vdw .eq. 2 ).and.(qvdw_flag)) qq_precomp(iq,jq,istate)%soft = .true.
end subroutine precompute_set_values_qq

!-----------------------------------------------------------------------

subroutine precompute_set_values_qw(iq,j,vdw)
! arguments
integer                         :: iq,i,j,vdw,istate,iacj,qvdw
! locals
real(kind=prec)                 :: tempA,tempB
! need to check first if vdw is right, qvdw is 1,1,3 instead of normal 1,2,3
if (vdw .eq. 2 ) then
qvdw = 1 
else
qvdw = vdw
end if

i = iqseq(iq)
iacj = iac(nat_solute+j)
do istate = 1, nstates
select case(ivdw_rule)
case(VDW_GEOMETRIC)
if (qvdw_flag) then
tempA = qavdw(qiac(iq,istate),qvdw) * aLJ_solv(j,vdw)
tempB = qbvdw(qiac(iq,istate),qvdw) * bLJ_solv(j,vdw)
else
tempA = iaclib(iac(i))%avdw(vdw) * aLJ_solv(j,vdw)
tempB = iaclib(iac(i))%bvdw(vdw) * bLJ_solv(j,vdw)
end if
qw_precomp(iq,j,istate)%vdWA = tempA
qw_precomp(iq,j,istate)%vdWB = tempB
case(VDW_ARITHMETIC)
if (qvdw_flag) then
tempA = qavdw(qiac(iq,istate),qvdw) + aLJ_solv(j,vdw)
tempB = qbvdw(qiac(iq,istate),qvdw) * bLJ_solv(j,vdw)
else
tempA = iaclib(iac(i))%avdw(vdw) + aLJ_solv(j,vdw)
tempB = iaclib(iac(i))%bvdw(vdw) * bLJ_solv(j,vdw)
end if
tempA = tempA**2
tempA = tempA * tempA * tempA
qw_precomp(iq,j,istate)%vdWA = (tempA**2) * tempB
qw_precomp(iq,j,istate)%vdWB = 2.0_prec * tempA * tempB
end select
if(.not. qq_use_library_charges) then
qw_precomp(iq,j,istate)%elec = qcrg(iq,istate) * chg_solv(j)
else
qw_precomp(iq,j,istate)%elec = crg(i) * chg_solv(j)
end if

qw_precomp(iq,j,istate)%score = sc_lookup(iq,iacj,istate)
end do

end subroutine precompute_set_values_qw

!-----------------------------------------------------------------------

subroutine precompute_set_values_ww(i,j,vdw)
! arguments
integer                         :: i,j,vdw
! locals
real(kind=prec)                 :: tempA,tempB
select case(ivdw_rule)
case(VDW_GEOMETRIC)
tempA = aLJ_solv(i,vdw) * aLJ_solv(j,vdw)
tempB = bLJ_solv(i,vdw) * bLJ_solv(j,vdw)
ww_precomp(i,j)%vdWA = tempA
ww_precomp(i,j)%vdWB = tempB

case(VDW_ARITHMETIC)
tempA = aLJ_solv(i,vdw) + aLJ_solv(j,vdw)
tempA = tempA**2
tempA = tempA * tempA * tempA
tempB = bLJ_solv(i,vdw) * bLJ_solv(j,vdw)

ww_precomp(i,j)%vdWA = (tempA**2) * tempB
ww_precomp(i,j)%vdWB = 2.0_prec * tempA * tempB

end select

ww_precomp(i,j)%set  = .true.
ww_precomp(i,j)%elec = chg_solv(i) * chg_solv(j)

end subroutine precompute_set_values_ww

!-----------------------------------------------------------------------
subroutine prep_sim
! local variables
integer						:: i, j, ig, istate,runvar,iw,irc_solvent
type(Q_ENE_HEAD)				:: tempheader

if (nodeid .eq. 0) then	
        write(*,*)
        call centered_heading('Initialising dynamics', '-')
end if

! Set parameters (bonds, angles, charges,...) & restraints for water   
! also get the number of atoms in each solvent molecule
! just take all non solute atoms and divide by nwat
! allocate chg_solv array here after getting number of solvent atoms/molecule
if(nwat > 0) then
! solv atom is read in during topo read now :)
	allocate(chg_solv(solv_atom))
	allocate(aLJ_solv(solv_atom,3))
	allocate(bLJ_solv(solv_atom,3))
        select case (solvent_type)
	 case (SOLVENT_SPC,SOLVENT_ALLATOM)
		do iw=1,solv_atom
			chg_solv(iw)=crg(nat_solute+iw)
			aLJ_solv(iw,1:3)=iaclib(iac(nat_solute+iw))%avdw(1:3)
			bLJ_solv(iw,1:3)=iaclib(iac(nat_solute+iw))%bvdw(1:3)
		end do
        case(SOLVENT_GENERAL)
                !add appropriate code here
                call die('Topology contains mixed solvent. This feature is not implemented yet.')
        end select

        if( .not. use_PBC ) then
! sovent density also set in wat_sphere to value from topology
! default will be  0.0335, water density
                call wat_sphere
                if (wpol_restr) call wat_shells

        else !compute charges of the system for box case 
        !(done in subroutine wat_sphere for sphere case)
                !calc. total charge of non-Q-atoms
                crgtot = zero
                do i = 1, nat_solute
                        if ( iqatom(i)==0 ) crgtot = crgtot + crg(i)
                end do
                write (*,60) crgtot
60 format ('Total charge of non-Q atoms             = ',f10.2)


        !calc effective charge of whole system at this lambda
                crgQtot = zero
                do i = 1, nqat
                        do istate = 1, nstates
                        crgtot = crgtot + qcrg(i,istate)*EQ(istate)%lambda
                        crgQtot = crgQtot + qcrg(i,istate)*EQ(istate)%lambda
                        end do
                end do

                write (*,70) crgtot

70 format ('Total charge of system                  = ',f10.2)

        end if
end if
! set the charge group membership for every topology atom only if using LRF or PBC
if(use_LRF .or. use_PBC) then
        call allocate_lrf_arrays

        do ig = 1, ncgp
                do i = cgp(ig)%first, cgp(ig)%last
                        iwhich_cgp(cgpatom(i)) = ig
                end do
        end do
!set nsolvent and nsolute to values from topology
!for box control later
	nsolvent = nres - nres_solute
	nsolute = nmol - nsolvent
end if

!	Prepare an array of inverse masses
winv(:) = one/iaclib(iac(:))%mass

if(use_PBC) allocate(mvd_mol(nmol))

if( use_PBC ) then 
        !compute masses of all molecules
        allocate(mol_mass(nmol))
        mol_mass(:) = zero

        do i = 1,nmol-1 !all molecules but the last
                do j = istart_mol(i), istart_mol(i+1)-1 !all atoms of molecule
                        mol_mass(i) = mol_mass(i) + iaclib(iac(j))%mass
                end do
        end do

        do j = istart_mol(nmol), natom !last molecule
                mol_mass(nmol) = mol_mass(nmol) + iaclib(iac(j))%mass
        end do

        mol_mass(:) = one/mol_mass(:)

        !prepare array of masses
        allocate(mass(natom))
!make this a duplicate of the iaclib array to avoid truncation errors
!        mass(:) = one/winv(:)
        mass(:)=iaclib(iac(:))%mass
! moved here or it will segfault :(
        if (control_box) then
        boxlength(:) = new_boxl(:)
        if ( put_solute_back_in_box .or. put_solvent_back_in_box ) then !only put back in box if either solute or solvent should be put back (qdyn input option)
                call put_back_in_box
        end if
        write(*,'(a)') 'Boxsize changed. Equilibration may be needed'
        end if
end if

!initialization for the Nosé-Hoover chain thermostat
!is done after shake has started, here shake_constraints are still zero

!scale charges by sqrt(coulomb_constant) 
crg(:) = crg(:) * sqrt(coulomb_constant)
if(nwat > 0) then
chg_solv(:) = chg_solv(:) * sqrt(coulomb_constant)
end if
if(nqat > 0) then
        qcrg(:,:) = qcrg(:,:) * sqrt(coulomb_constant)
end if
if (use_excluded_groups) then
!start writing the energy file header with the needed information
!one array for each group + 1 for the default environment
	write(*,*)
	write(*,*) 'Preparing residue groups for group contribution calculation'
	ene_header%arrays = ngroups_gc + 1
	allocate(ene_header%types(ngroups_gc+1),ene_header%numres(ngroups_gc+1),&
			ene_header%gcnum(ngroups_gc+1))
	runvar = 2
        do i=1,ngroups_gc
                call mask_initialize(ST_gc(i)%gcmask)
                call gc_make_array(ST_gc(i))
		ene_header%numres(runvar)=ST_gc(i)%count
		select case (trim(ST_gc(i)%caltypen))
			case ('full')
			ene_header%types(runvar)=FULL
			ene_header%gcnum(runvar)=i
			case ('electro')
			ene_header%types(runvar)=ELECTRO
			ene_header%gcnum(runvar)=i
			case ('vdw')
			ene_header%types(runvar)=VDW
			ene_header%gcnum(runvar)=i
			case ('all')
!special case, do all the calculations together
!needs reallocation of header arrays
			allocate(tempheader%types(runvar),tempheader%numres(runvar),tempheader%gcnum(runvar))
			tempheader%types(1:runvar)=ene_header%types(1:runvar)
			tempheader%numres(1:runvar)=ene_header%numres(1:runvar)
			tempheader%gcnum(1:runvar)=ene_header%gcnum(1:runvar)
			deallocate(ene_header%types,ene_header%numres,ene_header%gcnum)
!add two more calculations
!			runvar = ene_header%arrays + 2
			allocate(ene_header%types(ene_header%arrays+2),&
			ene_header%numres(ene_header%arrays+2),&
			ene_header%gcnum(ene_header%arrays+2))
			ene_header%types(1:runvar)=tempheader%types(1:runvar)
			ene_header%types(runvar)=FULL
			ene_header%types(runvar+1)=ELECTRO
			ene_header%types(runvar+2)=VDW
			ene_header%numres(1:runvar)=tempheader%numres(1:runvar)
			ene_header%numres(runvar+1)=ene_header%numres(runvar)
			ene_header%numres(runvar+2)=ene_header%numres(runvar)
			ene_header%gcnum(1:runvar)=tempheader%gcnum(1:runvar)
			ene_header%gcnum(runvar)=i
			ene_header%gcnum(runvar+1)=i
			ene_header%gcnum(runvar+2)=i
			ene_header%arrays=ene_header%arrays+2
			runvar = runvar + 2
			deallocate(tempheader%types,tempheader%numres,tempheader%gcnum)


			case default
!should never ever reach this, end execution now!
			write(*,*) 'Unrecognized statement in group contribution case list'
			call die('group contribution case list')
		end select
                do j=1,ST_gc(i)%count
                	stat = mask_add(ST_gc(i)%gcmask,ST_gc(i)%sendtomask(j))
                end do
		runvar = runvar + 1
	write (*,*) 'Finished group ',i
        end do
	do i=2,ene_header%arrays
		ene_header%totresid = ene_header%totresid + ene_header%numres(i)
	end do
else
!no extra groups, write minimal info to header
	ene_header%arrays = 1
	ene_header%totresid = 0 ! to prevent allocation errors
	allocate(ene_header%types(1),ene_header%numres(1),ene_header%gcnum(1))
end if
	ene_header%totresid = ene_header%totresid + 1
	allocate(ene_header%resid(ene_header%totresid))
	ene_header%types(1)=NOGC
	ene_header%numres(1)=1
	ene_header%resid(1)=-1
	ene_header%gcnum(1)=-1
	runvar=2
	do i=2,ene_header%arrays
!give each residue its topology number
		do j=1,ene_header%numres(i)
			ene_header%resid(runvar)=get_from_mask(ST_gc(ene_header%gcnum(i)),j)
			runvar = runvar + 1
		end do
	end do
!now that we know how many arrays we need -> allocate them in EQ
	do i=1,nstates
	allocate(EQ(i)%qx(ene_header%arrays),EQ(i)%qq(ene_header%arrays),&
		EQ(i)%qp(ene_header%arrays),EQ(i)%qw(ene_header%arrays),&
		EQ(i)%total(ene_header%arrays))
	end do

!if using PBC, alredy allocate the old_EQ arrays so we don't have to do it
!every time we update the volume
if( use_PBC ) then
	allocate(old_EQ(nstates))
        do i=1,nstates
        allocate(old_EQ(i)%qx(ene_header%arrays),old_EQ(i)%qq(ene_header%arrays),&
                old_EQ(i)%qp(ene_header%arrays),old_EQ(i)%qw(ene_header%arrays),&
                old_EQ(i)%total(ene_header%arrays))
        end do
end if
!and we can write the header to the energy file :)
!write (11) canary,ene_header%arrays,ene_header%totresid,ene_header%types(1:ene_header%arrays),&
!	ene_header%numres(1:ene_header%arrays),ene_header%resid(1:ene_header%totresid),&
!	ene_header%gcnum(1:ene_header%arrays)

! the last thing to prepare are the MD grids
! so we call them here
! they will be filled later when the pair lists are created the first time
#ifdef USE_GRID
call create_grid_pp
if (nwat .gt. 0) then
call create_grid_pw
call create_grid_ww
end if

! now allocate the array needed to store the grid lookup information for
! the charge groups
allocate(pp_igrid(ncgp_solute))
allocate(pw_igrid(ncgp_solute))
allocate(ww_igrid(nwat))
#endif

! prepare internal nonbonded interactions for solvent
call prep_sim_precompute_solvent

#if defined (USE_MPI)
	reclength=nstates*((2*ene_header%arrays)+(2*ene_header%arrays))
#endif
end subroutine prep_sim

!-----------------------------------------------------------------------

subroutine prep_sim_precompute_solvent
! locals
integer,allocatable		:: interaction(:,:)
integer				:: i,j,na,nb,counter,last
real(kind=prec)			:: tempA,tempB
type(SOLV_INT_TYPE),allocatable	:: tmp_solv_int(:)

if (ntors .gt. ntors_solute) then
! we need a way to troll the solvent interactions without having to build a list for the 1-4
! interactions and save it in the topology, this could kill topology reading otherwise
! so we start trolling them here before shake and other stuff is assigned
! This means we know all solvent internal stuff
! the maximum possible number is n^2, but we exclude self interaction and nearest neighbors

! fill array with all possible interaction codes
allocate(interaction(solv_atom,solv_atom),stat=alloc_status)
call check_alloc('solvent self interaction array')
do i = 1 , solv_atom
do j = 1 , solv_atom
interaction(i,j) = ljcod(iac(nat_solute+i),iac(nat_solute+j))
end do
end do

! first, exclude all self interactions
do i = 1 , solv_atom
interaction(i,i) = 0
end do

! need to know kast atom of first solvent molecule, because we can stop searching there :)
last = nat_solute + solv_atom

! now, troll the bond list for interactions to exclude
do i = nbonds_solute+1 , nbonds
na = bnd(i)%i
nb = bnd(i)%j
if (na.gt.last .or. nb.gt.last) cycle
! get the actual index of the solvent atom
! and set this index to zero
na = na - nat_solute
nb = nb - nat_solute
interaction(na,nb) = 0
interaction(nb,na) = 0
end do
! same for angle list
do i = nangles_solute+1, nangles
na = ang(i)%i
nb = ang(i)%k
if (na.gt.last .or. nb.gt.last) cycle
na = na - nat_solute
nb = nb - nat_solute
interaction(na,nb) = 0
interaction(nb,na) = 0
end do
! different now for torsion
! set to type 3 (1-4 interaction) if not previously set to zero
! and set inverted interaction to zero to prevent double counting
do i = ntors_solute+1, ntors
na = tor(i)%i
nb = tor(i)%l
if (na.gt.last .or. nb.gt.last) cycle
na = na - nat_solute
nb = nb - nat_solute
if ((interaction(na,nb) .ne. 0) .or. (interaction(nb,na) .ne. 0) ) then
interaction(na,nb) = 3
interaction(nb,na) = 0
end if
end do
counter = 0
! now we know the remaining interactions that we need to precalculate
! allocate full size array and later shrink it to the number needed
allocate(tmp_solv_int(solv_atom**2),stat=alloc_status)
call check_alloc('tmp solvent internal nonbonds')
do i = 1, solv_atom
do j = 1, solv_atom
if (interaction(i,j) .ne. 0 ) then
counter = counter + 1
tmp_solv_int(counter)%i=i
tmp_solv_int(counter)%j=j
select case(ivdw_rule)
case(VDW_GEOMETRIC)
tmp_solv_int(counter)%vdWA = aLJ_solv(i,interaction(i,j)) * aLJ_solv(j,interaction(i,j))
tmp_solv_int(counter)%vdWB = bLJ_solv(i,interaction(i,j)) * bLJ_solv(j,interaction(i,j))
case(VDW_ARITHMETIC)
tempA = aLJ_solv(i,interaction(i,j)) + aLJ_solv(j,interaction(i,j))
tempA = tempA**2
tempA = tempA * tempA * tempA
tempB = bLJ_solv(i,interaction(i,j)) * bLJ_solv(j,interaction(i,j))

tmp_solv_int(counter)%vdWA = (tempA**2) * tempB
tmp_solv_int(counter)%vdWB = 2.0_prec * tempA * tempB

end select

tmp_solv_int(counter)%elec = chg_solv(i) * chg_solv(j)
if (interaction(i,j) .eq. 3) &
tmp_solv_int(counter)%elec = tmp_solv_int(counter)%elec * el14_scale
end if
end do
end do
if (counter .le. 0) then
num_solv_int = 1
tmp_solv_int(:)%i=-1
tmp_solv_int(:)%j=-1
tmp_solv_int(:)%vdWA=-1E35_prec
tmp_solv_int(:)%vdWB=-1E35_prec
tmp_solv_int(:)%elec=-1E35_prec
else
num_solv_int = counter
end if

! now make real array to shrink done to the actual number of interactions
allocate(nonbnd_solv_int(num_solv_int),stat=alloc_status)
call check_alloc('solvent internal nonbonds')
nonbnd_solv_int(1:num_solv_int) = tmp_solv_int(1:num_solv_int)

! done, clean up
deallocate(tmp_solv_int)
deallocate(interaction)

else
! needed to prevent MPI from crashing during broadcast
num_solv_int = 1
allocate(nonbnd_solv_int(num_solv_int),stat=alloc_status)
call check_alloc('solvent internal nonbonds')
nonbnd_solv_int(:)%i=-1
nonbnd_solv_int(:)%j=-1
nonbnd_solv_int(:)%vdWA=-1E35_prec
nonbnd_solv_int(:)%vdWB=-1E35_prec
nonbnd_solv_int(:)%elec=-1E35_prec
end if


end subroutine prep_sim_precompute_solvent

!-----------------------------------------------------------------------

subroutine prep_sim_version(version)
! arguments
character(*)	:: version

! local values
integer :: canary = 1337

if (prec .eq. singleprecision) then
canary = 137
else if (prec .eq. doubleprecision) then
canary = 1337
#ifndef PGI
else if (prec .eq. quadprecision) then
canary = 13337
#endif
else
call die('No such precision')
end if
! hard coded value
!now write version number so Miha stops complaining
ene_header%version = trim(version)

if ( iene_cycle .gt. 0 ) then
write (11) canary,ene_header%arrays,ene_header%totresid,ene_header%types(1:ene_header%arrays),&
ene_header%numres(1:ene_header%arrays),ene_header%resid(1:ene_header%totresid),&
ene_header%gcnum(1:ene_header%arrays),ene_header%version
end if

end subroutine prep_sim_version


!-----------------------------------------------------------------------

subroutine qangle (istate)
! arguments
integer						:: istate

! local variables
integer						:: ia,i,j,k,ic,i3,j3,k3,im,icoupl,ib
real(kind=prec)						:: bji,bjk,scp,ang,da,ae,dv,gamma
real(kind=prec)						:: rji(3),rjk(3),f1,di(3),dk(3)


do ia = 1, nqangle


ic = qang(ia)%cod(istate)
 !skip if angle not present (code 0)
if ( ic > 0 ) then

gamma = one
icoupl = 0

do im = 1, nang_coupl
   if ( iang_coupl(1,im) .eq. ia ) then
      icoupl = im
      ib     = iang_coupl(2,im)
      gamma = EMorseD(ib)
	  !couple improper to bond breaking not making
	  if ( iang_coupl(3,im) .eq. 1) gamma = one - gamma
	  exit
   end if
end do

i  = qang(ia)%i
j  = qang(ia)%j
k  = qang(ia)%k
i3 = 3*i-3
j3 = 3*j-3
k3 = 3*k-3

rji(1) = x(i3+1) - x(j3+1)
rji(2) = x(i3+2) - x(j3+2)
rji(3) = x(i3+3) - x(j3+3)
rjk(1) = x(k3+1) - x(j3+1)
rjk(2) = x(k3+2) - x(j3+2)
rjk(3) = x(k3+3) - x(j3+3)
bji = sqrt ( rji(1)**2 + rji(2)**2 + rji(3)**2 )
bjk = sqrt ( rjk(1)**2 + rjk(2)**2 + rjk(3)**2 )
scp = ( rji(1)*rjk(1) + rji(2)*rjk(2) + rji(3)*rjk(3) )
scp = scp / (bji*bjk)
if ( scp .gt.  one ) scp =  one
if ( scp .lt. -one ) scp = -one
ang = acos(scp)
da = ang - qanglib(ic)%ang0
ae = 0.5_prec*qanglib(ic)%fk*da**2
EQ(istate)%q%angle = EQ(istate)%q%angle + ae*gamma

dv = gamma*qanglib(ic)%fk*da*EQ(istate)%lambda
f1 = sin ( ang )
if ( abs(f1) .lt. 1.e-12_prec ) f1 = 1.e-12_prec
f1 =  -one / f1
di(1) = f1 * ( rjk(1)/(bji*bjk) - scp*rji(1)/bji**2 )
di(2) = f1 * ( rjk(2)/(bji*bjk) - scp*rji(2)/bji**2 )
di(3) = f1 * ( rjk(3)/(bji*bjk) - scp*rji(3)/bji**2 )
dk(1) = f1 * ( rji(1)/(bji*bjk) - scp*rjk(1)/bjk**2 )
dk(2) = f1 * ( rji(2)/(bji*bjk) - scp*rjk(2)/bjk**2 )
dk(3) = f1 * ( rji(3)/(bji*bjk) - scp*rjk(3)/bjk**2 )
d(i3+1) = d(i3+1) + dv*di(1)
d(i3+2) = d(i3+2) + dv*di(2)
d(i3+3) = d(i3+3) + dv*di(3)
d(k3+1) = d(k3+1) + dv*dk(1)
d(k3+2) = d(k3+2) + dv*dk(2)
d(k3+3) = d(k3+3) + dv*dk(3)
d(j3+1) = d(j3+1) - dv*( di(1) + dk(1) )
d(j3+2) = d(j3+2) - dv*( di(2) + dk(2) )
d(j3+3) = d(j3+3) - dv*( di(3) + dk(3) )

if ( icoupl .ne. 0 ) then

   i  = qbnd(ib)%i
   j  = qbnd(ib)%j
   i3 = 3*i-3
   j3 = 3*j-3

   d(i3+1) = d(i3+1) + dMorse_i(1,ib)*ae
   d(i3+2) = d(i3+2) + dMorse_i(2,ib)*ae
   d(i3+3) = d(i3+3) + dMorse_i(3,ib)*ae
   d(j3+1) = d(j3+1) + dMorse_j(1,ib)*ae
   d(j3+2) = d(j3+2) + dMorse_j(2,ib)*ae
   d(j3+3) = d(j3+3) + dMorse_j(3,ib)*ae

end if
end if
end do
end subroutine qangle

!-----------------------------------------------------------------------

subroutine qurey_bradley (istate)
! arguments
integer						:: istate



! local variables
integer						::	ia,i,j,k,ic,i3,j3,k3,im,icoupl,ib
real(kind=prec)						::	gamma
real(kind=prec)						::	rik(3), dik, du, ru, Eurey

do ia = 1, nqangle
        ic = qang(ia)%cod(istate)
        !skip if angle not present (code 0)
        if ( ic == 0  .or. qanglib(ic)%ureyfk == zero) cycle

gamma = one
icoupl = 0

do im = 1, nang_coupl
   if ( iang_coupl(1,im) .eq. ia ) then
      icoupl = im
      ib     = iang_coupl(2,im)
      gamma = EMorseD(ib)
	  !couple improper to bond breaking not making
	  if ( iang_coupl(3,im) .eq. 1) gamma = 1 - gamma
   end if
end do

i  = qang(ia)%i
j  = qang(ia)%j
k  = qang(ia)%k
i3 = 3*i-3
j3 = 3*j-3
k3 = 3*k-3
        rik(1) = x(k3+1) - x(i3+1)
        rik(2) = x(k3+2) - x(i3+2)
        rik(3) = x(k3+3) - x(i3+3)
        dik = sqrt(rik(1)*rik(1) + rik(2)*rik(2) + rik(3)*rik(3))
        ru = dik - qanglib(ic)%ureyr0
        Eurey = qanglib(ic)%ureyfk*ru**2
        EQ(istate)%q%angle = EQ(istate)%q%angle + Eurey*gamma
        du = gamma*2*(qanglib(ic)%ureyfk*ru/dik)*EQ(istate)%lambda

        d(k3+1) = d(k3+1) + du*rik(1)
        d(k3+2) = d(k3+2) + du*rik(2)
        d(k3+3) = d(k3+3) + du*rik(3)
        d(i3+1) = d(i3+1) - du*rik(1)
        d(i3+2) = d(i3+2) - du*rik(2)
        d(i3+3) = d(i3+3) - du*rik(3)

if ( icoupl .ne. 0 ) then

   i  = qbnd(ib)%i
   j  = qbnd(ib)%j
   i3 = 3*i-3
   j3 = 3*j-3

   d(i3+1) = d(i3+1) + dMorse_i(1,ib)*Eurey
   d(i3+2) = d(i3+2) + dMorse_i(2,ib)*Eurey
   d(i3+3) = d(i3+3) + dMorse_i(3,ib)*Eurey


   d(j3+1) = d(j3+1) + dMorse_j(1,ib)*Eurey
   d(j3+2) = d(j3+2) + dMorse_j(2,ib)*Eurey
   d(j3+3) = d(j3+3) + dMorse_j(3,ib)*Eurey


end if
end do
end subroutine qurey_bradley

!-----------------------------------------------------------------------

subroutine qbond (istate)
! arguments
integer						::	istate

! local variables
integer						::	ib,i,j,ic,i3,j3
real(kind=prec)						::	b,db,be,dv,fexp
real(kind=prec)						::	rij(3)

do ib = 1, nqbond

        ic = qbnd(ib)%cod(istate)
        !code 0 means bond not present
if ( ic > 0 ) then

i  = qbnd(ib)%i
j  = qbnd(ib)%j
i3 = 3*i-3
j3 = 3*j-3

rij(1) = x(j3+1) - x(i3+1)
rij(2) = x(j3+2) - x(i3+2)
rij(3) = x(j3+3) - x(i3+3)

b = sqrt ( rij(1)**2 + rij(2)**2 + rij(3)**2 )
db = b - qbondlib(ic)%r0

        fexp = exp ( -qbondlib(ic)%amz*db ) 
        be = qbondlib(ic)%Dmz*(fexp*fexp-2.0_prec*fexp) + 0.5_prec*qbondlib(ic)%fk*db**2
        EMorseD(ib) = -(fexp*fexp-2.0_prec*fexp)
EQ(istate)%q%bond = EQ(istate)%q%bond + be
dv = (2.0_prec*qbondlib(ic)%Dmz*qbondlib(ic)%amz*(fexp-fexp*fexp) + qbondlib(ic)%fk*db)*EQ(istate)%lambda/b

d(i3+1) = d(i3+1) - dv*rij(1)
d(i3+2) = d(i3+2) - dv*rij(2)
d(i3+3) = d(i3+3) - dv*rij(3)
d(j3+1) = d(j3+1) + dv*rij(1)
d(j3+2) = d(j3+2) + dv*rij(2)
d(j3+3) = d(j3+3) + dv*rij(3)

!Force scaling factor to be 1 when distance is smaller than r0
 if ( db > 0 ) then
 	EMorseD(ib) = -(fexp*fexp-2.0_prec*fexp)
 	dMorse_i(:,ib) = +2.0_prec*qbondlib(ic)%amz*(fexp-fexp*fexp)*EQ(istate)%lambda/b*rij(:)
 	dMorse_j(:,ib) = -2.0_prec*qbondlib(ic)%amz*(fexp-fexp*fexp)*EQ(istate)%lambda/b*rij(:)
 else
 	EMorseD(ib) = one
 	dMorse_i(:,ib) = zero
 	dMorse_j(:,ib) = zero
 end if

end if
end do
end subroutine qbond

!-----------------------------------------------------------------------

subroutine qimproper (istate)
! arguments
integer						::	istate

! local variables
integer						::	i,j,k,l,ip,ic,i3,j3,k3,l3
integer						::	icoupl,im,ib
real(kind=prec)						::	bj,bk,scp,phi,sgn,pe,dv,arg,f1,gamma
real(kind=prec)						::	rji(3),rjk(3),rkl(3),rnj(3),rnk(3)
real(kind=prec)						::	rki(3),rlj(3),dp(12),di(3),dl(3)

do ip = 1,nqimp

ic = qimp(ip)%cod(istate)

if ( ic > 0 ) then

gamma = one
icoupl = 0

do im = 1, nimp_coupl
   if ( iimp_coupl(1,im) .eq. ip ) then
      icoupl = im
      ib     = iimp_coupl(2,im)
      gamma = EMorseD(ib)
	  !couple improper to bond breaking not making
	  if ( iimp_coupl(3,im) .eq. 1) gamma = 1 - gamma
   end if 
end do



i = qimp(ip)%i
j = qimp(ip)%j
k = qimp(ip)%k
l = qimp(ip)%l

i3=i*3-3
j3=j*3-3
k3=k*3-3
l3=l*3-3
rji(1) = x(i3+1) - x(j3+1)
rji(2) = x(i3+2) - x(j3+2)
rji(3) = x(i3+3) - x(j3+3)
rjk(1) = x(k3+1) - x(j3+1)
rjk(2) = x(k3+2) - x(j3+2)
rjk(3) = x(k3+3) - x(j3+3)
rkl(1) = x(l3+1) - x(k3+1)
rkl(2) = x(l3+2) - x(k3+2)
rkl(3) = x(l3+3) - x(k3+3)
rnj(1) =  rji(2)*rjk(3) - rji(3)*rjk(2)
rnj(2) =  rji(3)*rjk(1) - rji(1)*rjk(3)
rnj(3) =  rji(1)*rjk(2) - rji(2)*rjk(1)


rnk(1) = -rjk(2)*rkl(3) + rjk(3)*rkl(2)
rnk(2) = -rjk(3)*rkl(1) + rjk(1)*rkl(3)
rnk(3) = -rjk(1)*rkl(2) + rjk(2)*rkl(1)
bj = sqrt ( rnj(1)**2 + rnj(2)**2 + rnj(3)**2 )
bk = sqrt ( rnk(1)**2 + rnk(2)**2 + rnk(3)**2 )
scp = (rnj(1)*rnk(1)+rnj(2)*rnk(2)+rnj(3)*rnk(3))/(bj*bk)
if ( scp .gt.  one ) scp =  one
if ( scp .lt. -one ) scp = -one
phi = acos ( scp )
sgn =  rjk(1)*(rnj(2)*rnk(3)-rnj(3)*rnk(2)) &
     +rjk(2)*(rnj(3)*rnk(1)-rnj(1)*rnk(3)) &
     +rjk(3)*(rnj(1)*rnk(2)-rnj(2)*rnk(1))
if ( sgn .lt. 0 ) phi = -phi

! ---       energy

arg = phi - qimplib(ic)%imp0
arg = arg - 2.0_prec*pi*nint(arg/(2.0_prec*pi))
dv = qimplib(ic)%fk*arg
pe  = 0.5_prec*dv*arg
EQ(istate)%q%improper = EQ(istate)%q%improper + pe*gamma
dv = dv*gamma*EQ(istate)%lambda

! ---       forces

f1 = sin ( phi ) 
if ( abs(f1) .lt. 1.e-12_prec ) f1 = 1.e-12_prec
f1 =  -one / f1
di(1) = f1 * ( rnk(1)/(bj*bk) - scp*rnj(1)/bj**2 )
di(2) = f1 * ( rnk(2)/(bj*bk) - scp*rnj(2)/bj**2 )
di(3) = f1 * ( rnk(3)/(bj*bk) - scp*rnj(3)/bj**2 )
dl(1) = f1 * ( rnj(1)/(bj*bk) - scp*rnk(1)/bk**2 )
dl(2) = f1 * ( rnj(2)/(bj*bk) - scp*rnk(2)/bk**2 )
dl(3) = f1 * ( rnj(3)/(bj*bk) - scp*rnk(3)/bk**2 )

rki(1) =  rji(1) - rjk(1)
rki(2) =  rji(2) - rjk(2)
rki(3) =  rji(3) - rjk(3)
rlj(1) = -rjk(1) - rkl(1)
rlj(2) = -rjk(2) - rkl(2)
rlj(3) = -rjk(3) - rkl(3)

dp(1)  = rjk(2)*di(3) - rjk(3)*di(2)
dp(2)  = rjk(3)*di(1) - rjk(1)*di(3)
dp(3)  = rjk(1)*di(2) - rjk(2)*di(1)
dp(4)  = rki(2)*di(3)-rki(3)*di(2)+rkl(2)*dl(3)-rkl(3)*dl(2)
dp(5)  = rki(3)*di(1)-rki(1)*di(3)+rkl(3)*dl(1)-rkl(1)*dl(3)
dp(6)  = rki(1)*di(2)-rki(2)*di(1)+rkl(1)*dl(2)-rkl(2)*dl(1)
dp(7)  = rlj(2)*dl(3)-rlj(3)*dl(2)-rji(2)*di(3)+rji(3)*di(2)
dp(8)  = rlj(3)*dl(1)-rlj(1)*dl(3)-rji(3)*di(1)+rji(1)*di(3)
dp(9)  = rlj(1)*dl(2)-rlj(2)*dl(1)-rji(1)*di(2)+rji(2)*di(1)
dp(10) = rjk(2)*dl(3) - rjk(3)*dl(2)
dp(11) = rjk(3)*dl(1) - rjk(1)*dl(3)


dp(12) = rjk(1)*dl(2) - rjk(2)*dl(1)

d(i3+1) = d(i3+1) + dv*dp(1)
d(i3+2) = d(i3+2) + dv*dp(2)
d(i3+3) = d(i3+3) + dv*dp(3)
d(j3+1) = d(j3+1) + dv*dp(4)
d(j3+2) = d(j3+2) + dv*dp(5)
d(j3+3) = d(j3+3) + dv*dp(6)
d(k3+1) = d(k3+1) + dv*dp(7)
d(k3+2) = d(k3+2) + dv*dp(8)
d(k3+3) = d(k3+3) + dv*dp(9)
d(l3+1) = d(l3+1) + dv*dp(10)
d(l3+2) = d(l3+2) + dv*dp(11)
d(l3+3) = d(l3+3) + dv*dp(12)

if ( icoupl .ne. 0 ) then

   i  = qbnd(ib)%i
   j  = qbnd(ib)%j
   i3 = 3*i-3
   j3 = 3*j-3

   d(i3+1) = d(i3+1) + dMorse_i(1,ib)*pe
   d(i3+2) = d(i3+2) + dMorse_i(2,ib)*pe
   d(i3+3) = d(i3+3) + dMorse_i(3,ib)*pe
   d(j3+1) = d(j3+1) + dMorse_j(1,ib)*pe
   d(j3+2) = d(j3+2) + dMorse_j(2,ib)*pe
   d(j3+3) = d(j3+3) + dMorse_j(3,ib)*pe

end if
end if
end do
end subroutine qimproper

!-----------------------------------------------------------------------

subroutine qtorsion (istate)
! arguments
integer						::	istate

! local variables
integer						::	i,j,k,l,ip,ic,i3,j3,k3,l3
integer						::	icoupl,im,ib
real(kind=prec)						::	bj,bk,scp,phi,sgn,pe,dv,arg,f1,gamma
real(kind=prec)						::	rji(3),rjk(3),rkl(3),rnj(3),rnk(3)
real(kind=prec)						::	rki(3),rlj(3),dp(12),di(3),dl(3)

do ip = 1,nqtor

ic = qtor(ip)%cod(istate)

if ( ic > 0 ) then

gamma = one
icoupl = 0

do im = 1, ntor_coupl
   if ( itor_coupl(1,im) .eq. ip ) then
      icoupl = im
      ib     = itor_coupl(2,im)
      gamma = EMorseD(ib)
	  !couple improper to bond breaking not making
	  if ( itor_coupl(3,im) .eq. 1) gamma = 1 - gamma
   end if
end do

i = qtor(ip)%i
j = qtor(ip)%j
k = qtor(ip)%k
l = qtor(ip)%l

i3=i*3-3
j3=j*3-3
k3=k*3-3
l3=l*3-3
rji(1) = x(i3+1) - x(j3+1)
rji(2) = x(i3+2) - x(j3+2)
rji(3) = x(i3+3) - x(j3+3)
rjk(1) = x(k3+1) - x(j3+1)
rjk(2) = x(k3+2) - x(j3+2)
rjk(3) = x(k3+3) - x(j3+3)
rkl(1) = x(l3+1) - x(k3+1)
rkl(2) = x(l3+2) - x(k3+2)
rkl(3) = x(l3+3) - x(k3+3)
rnj(1) =  rji(2)*rjk(3) - rji(3)*rjk(2)
rnj(2) =  rji(3)*rjk(1) - rji(1)*rjk(3)
rnj(3) =  rji(1)*rjk(2) - rji(2)*rjk(1)
rnk(1) = -rjk(2)*rkl(3) + rjk(3)*rkl(2)
rnk(2) = -rjk(3)*rkl(1) + rjk(1)*rkl(3)
rnk(3) = -rjk(1)*rkl(2) + rjk(2)*rkl(1)
bj = sqrt ( rnj(1)**2 + rnj(2)**2 + rnj(3)**2 )
bk = sqrt ( rnk(1)**2 + rnk(2)**2 + rnk(3)**2 )
scp = (rnj(1)*rnk(1)+rnj(2)*rnk(2)+rnj(3)*rnk(3))/(bj*bk)
if ( scp .gt.  one ) scp =  one
if ( scp .lt. -one ) scp = -one
phi = acos ( scp )
sgn =  rjk(1)*(rnj(2)*rnk(3)-rnj(3)*rnk(2)) &
     +rjk(2)*(rnj(3)*rnk(1)-rnj(1)*rnk(3)) &


    +rjk(3)*(rnj(1)*rnk(2)-rnj(2)*rnk(1))
if ( sgn .lt. 0 ) phi = -phi




! ---       energy

arg = qtorlib(ic)%rmult*phi-qtorlib(ic)%deltor
pe = qtorlib(ic)%fk*(one+cos(arg))
EQ(istate)%q%torsion = EQ(istate)%q%torsion + pe*gamma
dv = -qtorlib(ic)%rmult*qtorlib(ic)%fk*sin(arg)*gamma*EQ(istate)%lambda

! ---       forces

f1 = sin ( phi ) 
if ( abs(f1) .lt. 1.e-12_prec ) f1 = 1.e-12_prec
f1 =  -one / f1
di(1) = f1 * ( rnk(1)/(bj*bk) - scp*rnj(1)/bj**2 )
di(2) = f1 * ( rnk(2)/(bj*bk) - scp*rnj(2)/bj**2 )
di(3) = f1 * ( rnk(3)/(bj*bk) - scp*rnj(3)/bj**2 )
dl(1) = f1 * ( rnj(1)/(bj*bk) - scp*rnk(1)/bk**2 )
dl(2) = f1 * ( rnj(2)/(bj*bk) - scp*rnk(2)/bk**2 )
dl(3) = f1 * ( rnj(3)/(bj*bk) - scp*rnk(3)/bk**2 )

rki(1) =  rji(1) - rjk(1)
rki(2) =  rji(2) - rjk(2)


rki(3) =  rji(3) - rjk(3)
rlj(1) = -rjk(1) - rkl(1)
rlj(2) = -rjk(2) - rkl(2)
rlj(3) = -rjk(3) - rkl(3)

dp(1)  = rjk(2)*di(3) - rjk(3)*di(2)
dp(2)  = rjk(3)*di(1) - rjk(1)*di(3)
dp(3)  = rjk(1)*di(2) - rjk(2)*di(1)
dp(4)  = rki(2)*di(3)-rki(3)*di(2)+rkl(2)*dl(3)-rkl(3)*dl(2)
dp(5)  = rki(3)*di(1)-rki(1)*di(3)+rkl(3)*dl(1)-rkl(1)*dl(3)
dp(6)  = rki(1)*di(2)-rki(2)*di(1)+rkl(1)*dl(2)-rkl(2)*dl(1)
dp(7)  = rlj(2)*dl(3)-rlj(3)*dl(2)-rji(2)*di(3)+rji(3)*di(2)
dp(8)  = rlj(3)*dl(1)-rlj(1)*dl(3)-rji(3)*di(1)+rji(1)*di(3)
dp(9)  = rlj(1)*dl(2)-rlj(2)*dl(1)-rji(1)*di(2)+rji(2)*di(1)
dp(10) = rjk(2)*dl(3) - rjk(3)*dl(2)
dp(11) = rjk(3)*dl(1) - rjk(1)*dl(3)
dp(12) = rjk(1)*dl(2) - rjk(2)*dl(1)

d(i3+1) = d(i3+1) + dv*dp(1)
d(i3+2) = d(i3+2) + dv*dp(2)
d(i3+3) = d(i3+3) + dv*dp(3)
d(j3+1) = d(j3+1) + dv*dp(4)
d(j3+2) = d(j3+2) + dv*dp(5)
d(j3+3) = d(j3+3) + dv*dp(6)
d(k3+1) = d(k3+1) + dv*dp(7)
d(k3+2) = d(k3+2) + dv*dp(8)
d(k3+3) = d(k3+3) + dv*dp(9)
d(l3+1) = d(l3+1) + dv*dp(10)
d(l3+2) = d(l3+2) + dv*dp(11)
d(l3+3) = d(l3+3) + dv*dp(12)

if ( icoupl .ne. 0 ) then

   i  = qbnd(ib)%i
   j  = qbnd(ib)%j
   i3 = 3*i-3
   j3 = 3*j-3



  d(i3+1) = d(i3+1) + dMorse_i(1,ib)*pe
   d(i3+2) = d(i3+2) + dMorse_i(2,ib)*pe
   d(i3+3) = d(i3+3) + dMorse_i(3,ib)*pe
   d(j3+1) = d(j3+1) + dMorse_j(1,ib)*pe
   d(j3+2) = d(j3+2) + dMorse_j(2,ib)*pe
   d(j3+3) = d(j3+3) + dMorse_j(3,ib)*pe

end if
end if
end do
end subroutine qtorsion

!-----------------------------------------------------------------------


real(kind=prec) function randm (ig)
! arguments
integer					::	ig

! local variables
integer, parameter		::	m = 100000000
integer, parameter		::  m1 = 10000
integer, parameter		::	mult=31415821
integer					::	irandh,irandl,multh,multl
real(kind=prec)					::	r
integer, save				::	irand = 0
integer, save				::  new = 0

if (new .eq. 0) then
new = 1
irand = mod (iabs(ig),m)
end if

! --- multiply irand by mult, but take into account that overflow must
! --- be discarded, and do not generate an error.
irandh = irand / m1
irandl = mod(irand, m1)
multh = mult / m1
multl = mod(mult, m1)

irand = mod(irandh*multl + irandl*multh, m1) * m1 + irandl*multl
irand = mod(irand + 1, m)

! --- convert irand to a real random number between 0 and 1.
r = real(irand / 10, kind=prec) * 10 / real(m, kind=prec)
if ((r .le. 0.e0_prec) .or. (r .gt. 1.e0_prec)) r = 0.e0_prec
randm = r
ig = irand

end function randm

!-----------------------------------------------------------------------

integer function shake(xx, x)
!arguments
real(kind=prec)						::	xx(:), x(:)
!	returns no. of iterations

! *** local variables
integer						::	i,j,i3,j3,mol,ic,nits
real(kind=prec)						::	xij2,diff,corr,scp,xxij2
real(kind=prec)						::	xij(3),xxij(3)
logical                                                 ::      dead_shake =.false.

#if defined (PROFILING)
real(kind=prec)                                         :: start_loop_time
start_loop_time = rtime()
#endif

! reset niter
shake = 0

do mol=1,shake_molecules
        ! for every molecule:
        ! reset nits (iterations per molecule)
        nits = 0
        ! reset iready for every constraint
        shake_mol(mol)%bond(:)%ready = .false.
        do while (.not.dead_shake) !iteration loop
                do ic=1,shake_mol(mol)%nconstraints
                        ! for every constraint:

                        if (.not. shake_mol(mol)%bond(ic)%ready) then
                                ! repeat until done:

                                i = shake_mol(mol)%bond(ic)%i
                                j = shake_mol(mol)%bond(ic)%j
                                i3 = i*3-3
                                j3 = j*3-3
                                xij(1)  = x(i3+1) - x(j3+1)
                                xij(2)  = x(i3+2) - x(j3+2)
                                xij(3)  = x(i3+3) - x(j3+3)
                                xij2    = xij(1)**2+xij(2)**2+xij(3)**2
                                diff    = shake_mol(mol)%bond(ic)%dist2 - xij2
                                if(abs(diff) < SHAKE_TOL*shake_mol(mol)%bond(ic)%dist2) then
                                        shake_mol(mol)%bond(ic)%ready = .true. ! in range
                                end if
                                xxij(1) = xx(i3+1) - xx(j3+1)
                                xxij(2) = xx(i3+2) - xx(j3+2)
                                xxij(3) = xx(i3+3) - xx(j3+3)
                                scp = xij(1)*xxij(1)+xij(2)*xxij(2)+xij(3)*xxij(3)
                                corr = diff/(2.0_prec*scp*(winv(i)+winv(j)))

                                x(i3+1) = x(i3+1)+xxij(1)*corr*winv(i)
                                x(i3+2) = x(i3+2)+xxij(2)*corr*winv(i)
                                x(i3+3) = x(i3+3)+xxij(3)*corr*winv(i)
                                x(j3+1) = x(j3+1)-xxij(1)*corr*winv(j)
                                x(j3+2) = x(j3+2)-xxij(2)*corr*winv(j)
                                x(j3+3) = x(j3+3)-xxij(3)*corr*winv(j)
                        end if
                end do

            nits = nits+1

                ! see if every constraint is met
                if(all(shake_mol(mol)%bond(1:shake_mol(mol)%nconstraints)%ready)) then
                        exit !from iteration loop
                elseif(nits >= SHAKE_MAX_ITER) then
                        ! fail on too many iterations
                        do ic=1,shake_mol(mol)%nconstraints
                                if (.not. shake_mol(mol)%bond(ic)%ready) then
                                        ! for every failed constraint

                                        i = shake_mol(mol)%bond(ic)%i
                                        j = shake_mol(mol)%bond(ic)%j
                                        i3 = i*3-3
                                        j3 = j*3-3
                                        xxij(1) = xx(i3+1) - xx(j3+1)
                                        xxij(2) = xx(i3+2) - xx(j3+2)
                                        xxij(3) = xx(i3+3) - xx(j3+3)
                                        xxij2   = xxij(1)**2+xxij(2)**2+xxij(3)**2
                                        write (*,100) i,j,sqrt(xxij2),&
                                                sqrt(shake_mol(mol)%bond(ic)%dist2)
                                end if
                        end do
                        dead_shake = .true.
                end if
100         format ('>>> Shake failed, i,j,d,d0 = ',2i6,2f10.5)
        end do

        ! update niter
        shake = shake+nits
end do

! set niter to the average number of iterations per molecule
if (dead_shake) call die('shake failure')
shake=shake/nmol
#if defined (PROFILING)
profile(7)%time = profile(7)%time + rtime() - start_loop_time
#endif

end function shake

!-----------------------------------------------------------------------

subroutine shrink_topology
!get rid of bonds and angles where all atoms are excluded
!or where the code has been set to 0 due to q-[bonds|angles|...]

!locals
integer						::	i, removed

if(exclude_bonded) then
        call centered_heading &
                ('Eliminating torsions & impropers for excluded atoms', '-')
end if

10	format('Reduced number of ',a,t31,'from ',i8,' to ')
12	format(i8)

i = 1
removed = 0
do while(i <= nbonds)
        !if all atoms excluded
        if(bnd(i)%cod <= 0) then
                !bond code either 0 (bond redefined in FEP file) 
                !or -1 (bond removed by shake)
                if(i <= nbonds_solute) then
                        bnd(i) = bnd(nbonds_solute)
                        bnd(nbonds_solute) = bnd(nbonds)
                        nbonds_solute = nbonds_solute - 1
                else
                        bnd(i) = bnd(nbonds)
                endif
                nbonds = nbonds - 1
                cycle !don't change i now, 
        end if
        i = i + 1 
end do

i = 1
do while(i <= nangles)
        !if all atoms excluded
        if(ang(i)%cod == 0) then
                !move last angle to current position
                if(i <= nangles_solute) then
                        ang(i) = ang(nangles_solute)
                        ang(nangles_solute) = ang(nangles)
                        nangles_solute = nangles_solute - 1
                else
                        ang(i) = ang(nangles)
                endif
                nangles = nangles - 1
                cycle !don't change i now, 
        end if
        i = i + 1 
end do

if(exclude_bonded) write(*,10, advance='no') 'torsions', ntors
i = 1
do while(i <= ntors)
        !if all atoms excluded
        if((exclude_bonded .and. excl(tor(i)%i) .and. excl(tor(i)%j) &
                .and. excl(tor(i)%k) .and. excl(tor(i)%l)) &
                .or. tor(i)%cod == 0) then
                !move last bond to current position
                if(i <= ntors_solute) then
                        tor(i) = tor(ntors_solute)
                        tor(ntors_solute) = tor(ntors)
                        ntors_solute = ntors_solute - 1
                else
                        tor(i) = tor(ntors)
                endif
                ntors = ntors - 1
                cycle !don't change i now, 
        end if
        i = i + 1 
end do
if(exclude_bonded) write(*, 12) ntors

if(exclude_bonded) write(*,10, advance='no') 'impropers', nimps
i = 1
do while(i <= nimps)
        !if all atoms excluded
        if(exclude_bonded .and. excl(imp(i)%i) .and. excl(imp(i)%j) &
                .and. excl(imp(i)%k) .and. excl(imp(i)%l) &
                .or. imp(i)%cod == 0) then
                if(i <= nimps_solute) then
                        imp(i) = imp(nimps_solute)
                        imp(nimps_solute) = imp(nimps)
                        nimps_solute = nimps_solute - 1
                else
                        imp(i) = imp(nimps)
                endif
                nimps = nimps - 1
                cycle !don't change i now, 
        end if
        i = i + 1 
end do
if(exclude_bonded) write(*, 12) nimps

end subroutine shrink_topology

!--------------------------------------------------------------------


subroutine stop_cm_translation
! local variables
integer						::	i,j,k
real(kind=prec)						::	rmass,totmass
real(kind=prec)						::	vcm(3)

! calculate totmass and vcm
totmass = zero
do i=1,3
vcm(i) = zero
end do
do i=1,natom
rmass = iaclib(iac(i))%mass
totmass=totmass+rmass
do j=1,3
k=(i-1)*3+j
vcm(j)=vcm(j)+rmass*v(k)
end do
end do

! scale vcm
do j=1,3
vcm(j)=vcm(j)/totmass
end do

! update v
do i=1,natom
do j=1,3
k=(i-1)*3+j
v(k)=v(k)-vcm(j)
end do
end do
end subroutine stop_cm_translation

!-----------------------------------------------------------------------
subroutine topology
! local variables
integer					::	nat3
integer					::	i
real(kind=prec)					::	box_min, vtemp, vtemp1

!
! read topology
!
! will init:
!  natom
!  lots of stuff from topo_load
!  nwat
!  anglib, torlib, implib (conversion)
!  ljcod
!  [iaclib%bvdw] (conversion)

if(.not. topo_load(top_file, require_version=4.15)) then
        call die('Failed to load topology.')
end if
natom = nat_pro

nwat = (natom - nat_solute) / solv_atom
!add extra element to molecule start atom array to keep track of last atom
istart_mol(nmol+1) = nat_pro + 1

! abort if no atoms
if (natom .eq. 0) call die('zero particles to simulate')

! convert libraries from degrees to radians
anglib(1:nangcod)%ang0 = deg2rad*anglib(1:nangcod)%ang0
torlib(1:ntorcod)%paths = one/torlib(1:ntorcod)%paths
torlib(1:ntorcod)%deltor = deg2rad*torlib(1:ntorcod)%deltor
implib(1:nimpcod)%imp0 = deg2rad*implib(1:nimpcod)%imp0

num_atyp = maxval(iac)

allocate(ljcod(num_atyp,num_atyp))

ljcod(:,:) = 1
do i=1,nlj2
ljcod(lj2(i)%i, lj2(i)%j) = 2
ljcod(lj2(i)%j, lj2(i)%i) = 2
end do


!	If arithmetic combination rule (ivdw_rule=2) take sqrt(epsilon) now
if ( ivdw_rule .eq. 2 ) then
do i=1,natyps
  iaclib(i)%bvdw(:) = sqrt(abs(iaclib(i)%bvdw(:)))
end do
end if

!check if same boundary in topology and Qdyn input file
if(  ((.not. use_PBC) .and.  box) .or. ( use_PBC .and. (.not. box)) ) then
call die('Must have same boundary (sphere or box) in topology and input file')
end if

if( use_PBC ) then
if( any(boxlength(:) == zero ) ) then
        inv_boxl(:) = zero
else
        inv_boxl(:) = one/boxlength(:)
end if

!check cut-offs if periodic box used
box_min = min( boxlength(1), boxlength(2), boxlength(3) )

!Solute-solute cut-off radii
if( .not. (box_min .gt. Rcpp*2) ) then
        call die('Solute-solute cut-off radii too large')

!Solvent-solvent
else if( .not. (box_min .gt. Rcww*2) ) then
        call die('Solvent-solvent cut-off radii too large')

!Solute-solvent
else if( .not. (box_min .gt. Rcpw*2) ) then
        call die('Solute-solvent cut-off radii too large')

!Q-atom
else if( .not. (box_min .gt. rcq*2) ) then
        call die('Q-atom cut-off radii too large')
!LRF
else if( .not. (box_min .gt. RcLRF*2) ) then
        call die('LRF cut-off radii too large')
end if
end if
end subroutine topology

!-----------------------------------------------------------------------

real(kind=prec) function torsion(istart, iend)
!arguments
integer						::	istart, iend

! local variables
integer						::	ip
real(kind=prec)						::	scp,phi,dv,arg,f1
real(kind=prec)						::	bjinv, bkinv, bj2inv, bk2inv
real(kind=prec)					::	rji(3),rjk(3),rkl(3),rnj(3),rnk(3)
real(kind=prec)					::	rki(3),rlj(3),dp(12),di(3),dl(3)
type(TOR_TYPE), pointer		::  t
type(TORLIB_TYPE), pointer	::	lib
#ifdef _OPENMP
integer :: quotient, remainder
#endif

! global variables used:
!  tor, torlib, x, d

! calculate the total energy of all torsion angles
! updates d

torsion = zero

!$omp parallel default(none) shared(threads_num, istart, iend, istep, tor, torlib, x, d, torsion) private(quotient, remainder,ip,scp,phi,dv,arg,f1,bjinv,bkinv,bj2inv,bk2inv,rji,rjk,rkl,rnj,rnk,rki,rlj,dp,di,dl,t,lib)
#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (iend - istart + 1)/threads_num
remainder = MOD(iend - istart + 1, threads_num)
mp_start = thread_id * quotient + istart + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif

mp_real_tmp = zero

do ip = mp_end, mp_start,-1

#else
do ip = iend, istart,-1
#endif
        t => tor(ip)
        lib => torlib(t%cod)
rji(1) = x(t%i*3-2) - x(t%j*3-2)
rji(2) = x(t%i*3-1) - x(t%j*3-1)
rji(3) = x(t%i*3-0) - x(t%j*3-0)
rjk(1) = x(t%k*3-2) - x(t%j*3-2)
rjk(2) = x(t%k*3-1) - x(t%j*3-1)
rjk(3) = x(t%k*3-0) - x(t%j*3-0)
rkl(1) = x(t%l*3-2) - x(t%k*3-2)
rkl(2) = x(t%l*3-1) - x(t%k*3-1)
rkl(3) = x(t%l*3-0) - x(t%k*3-0)

rnj(1) =  rji(2)*rjk(3) - rji(3)*rjk(2)
rnj(2) =  rji(3)*rjk(1) - rji(1)*rjk(3)
rnj(3) =  rji(1)*rjk(2) - rji(2)*rjk(1)
rnk(1) = -rjk(2)*rkl(3) + rjk(3)*rkl(2)
rnk(2) = -rjk(3)*rkl(1) + rjk(1)*rkl(3)
rnk(3) = -rjk(1)*rkl(2) + rjk(2)*rkl(1)

        bj2inv = one/(rnj(1)**2 + rnj(2)**2 + rnj(3)**2 )
bk2inv = one/(rnk(1)**2 + rnk(2)**2 + rnk(3)**2 )
        bjinv = sqrt(bj2inv)
        bkinv = sqrt(bk2inv)

        ! calculate scp and phi
scp = (rnj(1)*rnk(1)+rnj(2)*rnk(2)+rnj(3)*rnk(3))*(bjinv*bkinv)
if ( scp .gt.  one ) then
                scp =  one
                phi = acos (one) ! const
        else if ( scp .lt. -one ) then
                scp = -one
                phi = acos (-one) ! const
        else
    phi = acos ( scp )
        end if
if(rjk(1)*(rnj(2)*rnk(3)-rnj(3)*rnk(2)) &
                +rjk(2)*(rnj(3)*rnk(1)-rnj(1)*rnk(3)) &
    +rjk(3)*(rnj(1)*rnk(2)-rnj(2)*rnk(1)) .lt. 0) then
                phi = -phi
        end if

! ---       energy

arg = lib%rmult*phi-lib%deltor
#ifdef _OPENMP
    mp_real_tmp = mp_real_tmp + lib%fk*(one+cos(arg))*lib%paths   !lib%paths is previously inverted 
#else    
    torsion = torsion + lib%fk*(one+cos(arg))*lib%paths   !lib%paths is previously inverted 
#endif
dv = -lib%rmult*lib%fk*sin(arg)*lib%paths

! ---       forces

f1 = sin ( phi ) 
if ( abs(f1) .lt. 1.e-12_prec ) f1 = 1.e-12_prec
f1 =  -one / f1
di(1) = f1 * ( rnk(1)*(bjinv*bkinv) - scp*rnj(1)*bj2inv )
di(2) = f1 * ( rnk(2)*(bjinv*bkinv) - scp*rnj(2)*bj2inv )
di(3) = f1 * ( rnk(3)*(bjinv*bkinv) - scp*rnj(3)*bj2inv )
dl(1) = f1 * ( rnj(1)*(bjinv*bkinv) - scp*rnk(1)*bk2inv )
dl(2) = f1 * ( rnj(2)*(bjinv*bkinv) - scp*rnk(2)*bk2inv )
dl(3) = f1 * ( rnj(3)*(bjinv*bkinv) - scp*rnk(3)*bk2inv )

rki(1) =  rji(1) - rjk(1)
rki(2) =  rji(2) - rjk(2)
rki(3) =  rji(3) - rjk(3)
rlj(1) = -rjk(1) - rkl(1)
rlj(2) = -rjk(2) - rkl(2)
rlj(3) = -rjk(3) - rkl(3)

dp(1)  = rjk(2)*di(3) - rjk(3)*di(2)
dp(2)  = rjk(3)*di(1) - rjk(1)*di(3)
dp(3)  = rjk(1)*di(2) - rjk(2)*di(1)
dp(4)  = rki(2)*di(3)-rki(3)*di(2)+rkl(2)*dl(3)-rkl(3)*dl(2)
dp(5)  = rki(3)*di(1)-rki(1)*di(3)+rkl(3)*dl(1)-rkl(1)*dl(3)
dp(6)  = rki(1)*di(2)-rki(2)*di(1)+rkl(1)*dl(2)-rkl(2)*dl(1)
dp(7)  = rlj(2)*dl(3)-rlj(3)*dl(2)-rji(2)*di(3)+rji(3)*di(2)
dp(8)  = rlj(3)*dl(1)-rlj(1)*dl(3)-rji(3)*di(1)+rji(1)*di(3)
dp(9)  = rlj(1)*dl(2)-rlj(2)*dl(1)-rji(1)*di(2)+rji(2)*di(1)
dp(10) = rjk(2)*dl(3) - rjk(3)*dl(2)
dp(11) = rjk(3)*dl(1) - rjk(1)*dl(3)
dp(12) = rjk(1)*dl(2) - rjk(2)*dl(1)

d(t%i*3-2) = d(t%i*3-2) + dv*dp(1)
d(t%i*3-1) = d(t%i*3-1) + dv*dp(2)
d(t%i*3-0) = d(t%i*3-0) + dv*dp(3)
d(t%j*3-2) = d(t%j*3-2) + dv*dp(4)
d(t%j*3-1) = d(t%j*3-1) + dv*dp(5)
d(t%j*3-0) = d(t%j*3-0) + dv*dp(6)
d(t%k*3-2) = d(t%k*3-2) + dv*dp(7)
d(t%k*3-1) = d(t%k*3-1) + dv*dp(8)
d(t%k*3-0) = d(t%k*3-0) + dv*dp(9)
d(t%l*3-2) = d(t%l*3-2) + dv*dp(10)
d(t%l*3-1) = d(t%l*3-1) + dv*dp(11)
d(t%l*3-0) = d(t%l*3-0) + dv*dp(12)
end do
#ifdef _OPENMP
!$omp atomic update
torsion = torsion + mp_real_tmp
#endif
!$omp end parallel

end function torsion

!-----------------------------------------------------------------------
subroutine restrain_solvent 
! local variables
integer						::	iw,i,i3,isolv,jsolv
real(kind=prec)						::	b,db,erst,dv,fexp
real(kind=prec) 					::	dr(3),dcent(3)
real(kind=prec)						::	shift

! global variables used:
!  E, Boltz, Tfree, fk_wsphere, nwat, nat_pro, x, xwcent, rwat, Dwmz, awmz, d

if(fk_wsphere .ne. zero) then
        shift = sqrt (Boltz*Tfree/fk_wsphere)
else
        shift = zero
end if
do iw = ncgp_solute + 1, ncgp
        i  = cgp(iw)%iswitch
        if (excl(i)) cycle ! skip excluded topology waters
        if(solv_atom.eq.3) then
!use old code
        i3 = 3*i-3
        dr(1) = x(i3+1) - xwcent(1)
        dr(2) = x(i3+2) - xwcent(2)
        dr(3) = x(i3+3) - xwcent(3)
        else
        dcent(:)=zero
        jsolv = iw - ncgp_solute
        i = nat_solute + solv_atom*jsolv - (solv_atom-1)
        do isolv = 0 , solv_atom - 1
                i3 = 3*(i+isolv)-3
                dcent(1) = dcent(1) + x(i3+1)
                dcent(2) = dcent(2) + x(i3+2)
                dcent(3) = dcent(3) + x(i3+3)
        end do
        dcent(:) = dcent(:)/solv_atom
        dr(:) = dcent(:) - xwcent(:)
        end if
        b = sqrt ( dr(1)**2 + dr(2)**2 + dr(3)**2 )
        db = b - (rwat - shift)

        ! calculate erst and dv
        if ( db > 0 ) then
                erst = 0.5_prec * fk_wsphere * db**2 - Dwmz
                dv = fk_wsphere*db/b
        else
                if (b > zero) then
                  fexp = exp ( awmz*db )
                  erst = Dwmz*(fexp*fexp-2.0_prec*fexp)
                  dv = -2.0_prec*Dwmz*awmz*(fexp-fexp*fexp)/b
                else
                  dv = 0
                  erst = 0
                end if
        end if
        if (solvent_type .gt. 0) then
                dv = dv / solv_atom
        end if

        ! update energy and forces
        E%restraint%solvent_radial = E%restraint%solvent_radial + erst
        if (solvent_type .gt. 0) then
                jsolv = iw - ncgp_solute
                i = nat_solute + solv_atom*jsolv - (solv_atom-1)
                do isolv = 0 , solv_atom - 1
                i3 = 3*(i+isolv)-3
                d(i3+1) = d(i3+1) + dv*dr(1)
                d(i3+2) = d(i3+2) + dv*dr(2)
                d(i3+3) = d(i3+3) + dv*dr(3)
                end do
        else
                d(i3+1) = d(i3+1) + dv*dr(1)
                d(i3+2) = d(i3+2) + dv*dr(2)
                d(i3+3) = d(i3+3) + dv*dr(3)
        end if
end do
end subroutine restrain_solvent

!-----------------------------------------------------------------------

subroutine wat_sphere
! local variables
integer					:: i,i3,kr,isort,int_wat,istate
real(kind=prec)					:: rc,rnwat
real(kind=prec), save				:: dr(3)
real(kind=prec)					:: crgexcl

!possibly override target sphere radius from topology
if(rwat_in > zero) rwat = rwat_in


!calc. total charge of non-excluded non-Q-atoms and excluded atoms
crgtot = zero
crgexcl = zero
rho_wat = topo_rho_wat
do i = 1, nat_solute
	if ( .not. excl(i) ) then
		if ( iqatom(i)==0 ) then
			crgtot = crgtot + crg(i)
		end if
	else
		crgexcl = crgexcl + crg(i)
	end if
end do
write (*,60) 'non-Q atoms', crgtot
write (*,60) 'excluded atoms', crgexcl
60 format ('Total charge of ',a,t41,'= ',f10.2)

!calc effective charge of simulation sphere at this lambda
crgQtot = zero
do i = 1, nqat
	do istate = 1, nstates
		crgtot = crgtot + qcrg(i,istate)*EQ(istate)%lambda
		crgQtot = crgQtot + qcrg(i,istate)*EQ(istate)%lambda
	end do
end do
write (*,70) crgtot
70 format ('Total charge of system                  = ',f10.2)

if (.not. wpol_born) crgQtot = zero !ignore total Q if Born corr. is off
if ( nwat .eq. 0 ) return


!	Set default values for unspecified optional parameters
if(fk_wsphere == -1) then
        !
        ! To be replaced by function of rc giving appropriate default for any sphere
        !
        fk_wsphere = fk_wsphere_default
end if
if(fkwpol == -1) then
        !
        ! To be replaced by function of rc giving appropriate default for any sphere
        !
        fkwpol = fkwpol_default
end if
if(Dwmz == -1) then !Use magic function to get suitable Dwmz
        Dwmz = 0.26_prec*exp(-0.19_prec*(rwat-15.0_prec))+0.74_prec
end if
if(awmz == -1) then !use magic for the reach of the Morse potential
        awmz = 0.2_prec/(one+exp(0.4_prec*(rwat-25.0_prec)))+0.3_prec
end if

write (*,90) rwat, fk_wsphere, Dwmz, awmz
90	format ('Target water sphere radius              = ',f10.2,/,&
                'Surface inward harmonic force constant  = ',f10.2,/,&
                'Surface attraction well depth           = ',f10.2,/,&
                'Surface attraction well width           = ',f10.2)
92	format ('Water polarisation restraints           : ',a)
if(.not. wpol_restr) then
        write(*,92) 'OFF'
else if(wpol_born) then
        write(*,92) 'ON, Born correction enabled'
        write(*, 100) fkwpol
else 
        write(*,92) 'ON, Born correction disabled'
        write(*, 100) fkwpol
end if
100	format('Radial polarisation force constant      = ',f10.2)

end subroutine wat_sphere

!-----------------------------------------------------------------------

subroutine wat_shells
! set up the shells for polarisation restraining

! local variables
real(kind=prec)						::	rout, dr, ri, Vshell, rshell, drs, eps_diel
integer						::	is, n_insh


integer						::	nwpolr_shell_restart, filestat
integer						::	bndcodw, angcodw

logical :: old_restart = .false.
real(8),allocatable :: x_old(:), v_old(:)
real(kind=singleprecision),allocatable :: x_1(:), v_1(:)
real(kind=doubleprecision),allocatable :: x_2(:), v_2(:)
#ifndef PGI
real(kind=quadprecision),allocatable :: x_3(:), v_3(:)
#endif
real(8) :: old_boxlength(3),old_boxcentre(3)
real(kind=singleprecision) :: boxlength_1(3),boxcentre_1(3)
real(kind=doubleprecision) :: boxlength_2(3),boxcentre_2(3)
#ifndef PGI
real(kind=quadprecision) :: boxlength_3(3),boxcentre_3(3)
#endif
integer :: headercheck,i,myprec,dummy

if (prec .eq. singleprecision) then
myprec = -137
elseif (prec .eq. doubleprecision) then
myprec = -1337
#ifndef PGI
elseif (prec .eq. quadprecision) then
myprec = -13337
#endif
else
call die('No such precision')
end if

!calc mu_w
!look up bond code for water
bndcodw = bnd(nbonds)%cod
angcodw = ang(nangles)%cod
!find charge of water O = charge of 1st solvent atom

mu_w = -chg_solv(1)*bondlib(bndcodw)%bnd0*cos(anglib(angcodw)%ang0/2)

! shell widths are drout, 2drout, 3drout
drs = wpolr_layer / drout !number of drouts 

! calc number of shells based on arithmetic series sum formula
nwpolr_shell = int(-0.5_prec + sqrt(2*drs + 0.25_prec)) 
allocate(wshell(nwpolr_shell), stat=alloc_status)
call check_alloc('water polarisation shell array')

write(*, 100) nwpolr_shell
100	format(/,'Setting up ', i1, ' water shells for polarisation restraints.')

if(restart) then !try to load theta_corr from restart file
!first, rewind file to find out if we have an old or new restart
   rewind(2)
   read(2) headercheck
   if ((headercheck .ne. -137).and.(headercheck .ne. -1337).and.(headercheck .ne. -13337)) then
     old_restart = .true.
     rewind(2)
     allocate(x_old(3*natom),v_old(3*natom))
     read(2) dummy,(x_old(i),i=1,nat3)
     read(2) dummy,(v_old(i),i=1,nat3)
     if( use_PBC) then
       read(2) old_boxlength(:)
       read(2) old_boxcentre(:)
     end if
     deallocate(x_old,v_old)
   else if (headercheck .eq. -137) then
     allocate(x_1(3*natom),v_1(3*natom))
     read(2) dummy,(x_1(i),i=1,nat3)
     read(2) dummy,(v_1(i),i=1,nat3)
     if( use_PBC) then
       read(2) boxlength_1(:)
       read(2) boxcentre_1(:)
     end if
     deallocate(x_1,v_1)
   else if (headercheck .eq. -1337) then
     allocate(x_2(3*natom),v_2(3*natom))
     read(2) dummy,(x_2(i),i=1,nat3)
     read(2) dummy,(v_2(i),i=1,nat3)
     if( use_PBC) then
       read(2) boxlength_2(:)
       read(2) boxcentre_2(:)
     end if
     deallocate(x_2,v_2)
   else if (headercheck .eq. -13337) then
#ifndef PGI
     allocate(x_3(3*natom),v_3(3*natom))
     read(2) dummy,(x_3(i),i=1,nat3)
     read(2) dummy,(v_3(i),i=1,nat3)
     if( use_PBC) then
       read(2) boxlength_3(:)
       read(2) boxcentre_3(:)
     end if
     deallocate(x_3,v_3)
#else
     call die('Quadruple precision not supported in PGI')
#endif
   end if
        read(2, iostat=filestat) nwpolr_shell_restart
        if(filestat .ne. 0 .or. nwpolr_shell_restart /= nwpolr_shell) then
                write(*,102) 
                wshell(:)%theta_corr = zero
        else
                backspace(2)
                if (old_restart) then
                   allocate(old_wshell(nwpolr_shell), stat=alloc_status)
                   read(2) nwpolr_shell_restart,old_wshell(:)%theta_corr
                   wshell(:)%theta_corr = old_wshell(:)%theta_corr
                   write(*,*) 'Loaded water polarization data from old qdyn restart file'
                   deallocate(old_wshell)
                else
                   if (headercheck .ne. myprec) then
                   write(*,*) '>>> WARNING: Using water polarisation from restart with mismatched precision'
                     if (headercheck .eq. -137) then
                       allocate(wshell_single(nwpolr_shell), stat=alloc_status)
                       read(2) nwpolr_shell_restart,wshell_single(:)%theta_corr
                       wshell(:)%theta_corr = wshell_single(:)%theta_corr
                       deallocate(wshell_single)
                     else if (headercheck .eq. -1337) then
                       allocate(wshell_double(nwpolr_shell), stat=alloc_status)
                       read(2) nwpolr_shell_restart,wshell_double(:)%theta_corr
                       wshell(:)%theta_corr = wshell_double(:)%theta_corr
                       deallocate(wshell_double)
                     else if (headercheck .eq. -13337) then
#ifndef PGI
                       allocate(wshell_quad(nwpolr_shell), stat=alloc_status)
                       read(2) nwpolr_shell_restart,wshell_quad(:)%theta_corr
                       wshell(:)%theta_corr = wshell_quad(:)%theta_corr
                       deallocate(wshell_quad)
#else
                       call die('Quadruple precision not supported in PGI')
#endif
                     end if
                   else
                     read(2) nwpolr_shell_restart,wshell(:)%theta_corr
                   end if
                end if
                write(*,103)
        end if
else
        wshell(:)%theta_corr = zero
end if

102	format('>>> WARNING: Failed to read polarisation restraint data from restart file.')
103	format('Loaded polarisation restraint data from restart file.')

write(*,'(a)') 'Shell #    outer radius    inner radius'
110	format(i7, 2f16.2)

rout = rwat
n_max_insh = 0
do is = 1, nwpolr_shell 
        wshell(is)%avtheta = 0
        wshell(is)%avn_insh = 0
        wshell(is)%rout = rout
    dr = drout*is
        ri = rout - dr
        wshell(is)%dr = dr
        Vshell = rout**3 - ri**3
        n_insh = int(4 * pi/3 * Vshell * rho_wat)
        if (n_insh > n_max_insh) n_max_insh = n_insh
rshell = (0.5_prec*(rout**3+ri**3))**(one/3.0_prec)


        ! --- Note below: 0.98750 = (1-1/epsilon) for water
        ! now deprecated because we actually read the dielectric from the
        ! solvent files as dielectric = epsilon * 1000
        eps_diel = real(dielectric/1000,kind=prec)
        eps_diel = one-(one/eps_diel)
!wshell(is)%cstb = crgQtot*0.98750_prec/(rho_wat*mu_w*4.0_prec*pi*rshell**2)
        wshell(is)%cstb = crgQtot*eps_diel/(rho_wat*mu_w*4.0_prec*pi*rshell**2)
        write(*, 110) is, rout, ri
        rout = rout - dr
end do

n_max_insh = n_max_insh * 1.5_prec !take largest and add some extra
call allocate_watpol_arrays

end subroutine wat_shells

!-----------------------------------------------------------------------

subroutine watpol
! local variables
integer						:: iw,is,i,i3,il,jl,jw,imin,jmin,isolv,jsolv,j3
real(kind=prec)						:: dr,rw,rshell,rm,rc,scp
real(kind=prec)						:: tmin,arg,avtdum,dv,f0
real(kind=prec), save					:: f1(3),f2(3),f3(3)
real(kind=prec) 					:: rmu(3),rcu(3), rmc(3)

! global variables used:
!  E, wshell, bndw0, deg2rad, angw0, nwat, theta, theta0, nat_pro, x, xwcent,
!  tdum, nwpolr_shell, list_sh, pi, nsort, istep, itdis_update, fkwpol, d

! reset wshell%n_insh
wshell(:)%n_insh = 0

! calculate theta(:), tdum(:), wshell%n_insh
do iw = 1, nwat
theta(iw)  = zero
theta0(iw) = zero

i  = nat_solute + solv_atom*iw - (solv_atom-1) 
if(excl(i)) cycle ! skip excluded topology solvent
i3 = i*3-3

! function needs to be rewritten as we no longer just have 3-point water
! molecules, meaning that the simple geometry won't work any longer
! we will only keep the old code for SPC/TIP3P like solvent

rmu(:)=zero ! solvent vector
rmc(:)=zero
if (solv_atom .eq. 3) then
!SPC/TIP3P solvent, keep old code
do isolv = 1, solv_atom - 1
jsolv = i + isolv
j3 = jsolv*3-3
rmu(1) = rmu(1) + x(j3+1)
rmu(2) = rmu(2) + x(j3+2)
rmu(3) = rmu(3) + x(j3+3)
end do
rmu(1) = rmu(1) - 2.0_prec*x(i3+1)   
rmu(2) = rmu(2) - 2.0_prec*x(i3+2)
rmu(3) = rmu(3) - 2.0_prec*x(i3+3)
else
! no longer 3-point solvent, create solvent vector using individual atom charges
! of the solvent molecule
do isolv = 0, solv_atom - 1
jsolv = i + isolv
j3 = jsolv*3-3
rmc(1) = rmc(1) + x(j3+1)
rmc(2) = rmc(2) + x(j3+2)
rmc(3) = rmc(3) + x(j3+3)
rmu(1) = rmu(1) + (chg_solv(isolv+1)*x(j3+1))
rmu(2) = rmu(2) + (chg_solv(isolv+1)*x(j3+2))
rmu(3) = rmu(3) + (chg_solv(isolv+1)*x(j3+3))
end do
rmc(:) = rmc(:)/solv_atom
end if

rm = sqrt ( rmu(1)**2 + rmu(2)**2 + rmu(3)**2 )
rmu(1) = rmu(1)/rm
rmu(2) = rmu(2)/rm
rmu(3) = rmu(3)/rm

if (solv_atom .eq. 3) then
!SPC/TIP3P solvent, keep old code
rcu(1) = x(i3+1) - xwcent(1)    !Radial vector to center solvent atom
rcu(2) = x(i3+2) - xwcent(2) 
rcu(3) = x(i3+3) - xwcent(3) 
else
! different solvent, use new code for geometric center
rcu(:) = rmc(:) - xwcent(:)
end if
rc = sqrt ( rcu(1)**2 + rcu(2)**2 + rcu(3)**2 )
rcu(1) = rcu(1)/rc
rcu(2) = rcu(2)/rc
rcu(3) = rcu(3)/rc

scp = rmu(1)*rcu(1)+rmu(2)*rcu(2)+rmu(3)*rcu(3)   !Calculate angle between solvent vector and radial vector
if ( scp .gt.  one ) scp =  one
if ( scp .lt. -one ) scp = -one
theta(iw) = acos( scp )
tdum(iw) = theta(iw)

if ( rc > wshell(nwpolr_shell)%rout-wshell(nwpolr_shell)%dr ) then
  do is = nwpolr_shell, 2, -1
        if(rc <= wshell(is)%rout) exit
  end do
  wshell(is)%n_insh = wshell(is)%n_insh + 1
  if (wshell(is)%n_insh .ge. n_max_insh) then
	call reallocate_watpol_arrays
  end if
  list_sh(wshell(is)%n_insh,is) = iw
end if
end do

! sort the solvent molecules according to theta
do is = 1, nwpolr_shell
imin = 0
do il = 1, wshell(is)%n_insh
tmin = 2.0_prec*pi
do jl = 1, wshell(is)%n_insh
jw = list_sh(jl,is)
if ( tdum(jw) .lt. tmin ) then
  jmin = jw
  tmin = theta(jw)
end if
end do
imin = imin+1
nsort(imin,is) = jmin
  tdum(jmin) = 99999.0_prec

end do

end do

! calculate energy and force
if ( istep .ne. 0 .and. mod(istep,itdis_update) .eq. 0) then
call centered_heading('Solvent polarisation restraint data', '-')
write(*,'(a)') 'shell    <n>    <theta>    theta_0 theta_corr'
do is = 1, nwpolr_shell
wshell(is)%avtheta = wshell(is)%avtheta / real (itdis_update, kind=prec)
wshell(is)%avn_insh = wshell(is)%avn_insh / real (itdis_update, kind=prec)
wshell(is)%theta_corr = wshell(is)%theta_corr + wshell(is)%avtheta-acos(wshell(is)%cstb)
write (*,10) is,wshell(is)%avn_insh,wshell(is)%avtheta/deg2rad, &
     acos(wshell(is)%cstb)/deg2rad,wshell(is)%theta_corr/deg2rad
10	  format(i5,1x,f6.1,3x,f8.3,3x,f8.3,3x,f8.3)
wshell(is)%avtheta = zero
wshell(is)%avn_insh = zero
end do
end if

do is = 1, nwpolr_shell
if(wshell(is)%n_insh == 0) cycle !skip empty shell
avtdum = zero
do il = 1, wshell(is)%n_insh
iw = nsort(il,is)
arg = one + (one - 2.0_prec*real(il, kind=prec))/real(wshell(is)%n_insh, kind=prec)
theta0(il) = acos ( arg )
theta0(il) = theta0(il)-3.0_prec*sin(theta0(il))*wshell(is)%cstb/2.0_prec
if ( theta0(il) .lt. zero ) theta0(il) = zero
if ( theta0(il) .gt. pi)   theta0(il) = pi

avtdum = avtdum + theta(iw)

E%restraint%water_pol = E%restraint%water_pol + 0.5_prec*fkwpol* &
     (theta(iw)-theta0(il)+wshell(is)%theta_corr)**2

dv = fkwpol*(theta(iw)-theta0(il)+wshell(is)%theta_corr)

i  = nat_solute + solv_atom*iw - (solv_atom-1)
i3 = i*3-3
rmu(:)=zero ! solvent vector
rmc(:)=zero
if (solv_atom .eq. 3) then
!SPC/TIP3P solvent, keep old code
do isolv = 1, solv_atom - 1
jsolv = i + isolv
j3 = jsolv*3-3
rmu(1) = rmu(1) + x(j3+1)
rmu(2) = rmu(2) + x(j3+2)
rmu(3) = rmu(3) + x(j3+3)
end do
rmu(1) = rmu(1) - 2.0_prec*x(i3+1)
rmu(2) = rmu(2) - 2.0_prec*x(i3+2)
rmu(3) = rmu(3) - 2.0_prec*x(i3+3)
else
! no longer 3-point solvent, create solvent vector based on individual charges of
! the solvent molecule
do isolv = 0, solv_atom - 1
jsolv = i + isolv
j3 = jsolv*3-3
rmc(1) = rmc(1) + x(j3+1)
rmc(2) = rmc(2) + x(j3+2)
rmc(3) = rmc(3) + x(j3+3)
rmu(1) = rmu(1) + (chg_solv(isolv+1)*x(j3+1))
rmu(2) = rmu(2) + (chg_solv(isolv+1)*x(j3+2))
rmu(3) = rmu(3) + (chg_solv(isolv+1)*x(j3+3))
end do
rmc(:) = rmc(:)/solv_atom
end if
rm = sqrt ( rmu(1)**2 + rmu(2)**2 + rmu(3)**2 )
rmu(1) = rmu(1)/rm
rmu(2) = rmu(2)/rm
rmu(3) = rmu(3)/rm
if (solv_atom .eq. 3) then
!SPC/TIP3P solvent, keep old code
rcu(1) = x(i3+1) - xwcent(1)
rcu(2) = x(i3+2) - xwcent(2)
rcu(3) = x(i3+3) - xwcent(3)
else
! different solvent, use new code for geometric center
rcu(:) = rmc(:) - xwcent(:)
end if
rc = sqrt ( rcu(1)**2 + rcu(2)**2 + rcu(3)**2 )
rcu(1) = rcu(1)/rc
rcu(2) = rcu(2)/rc
rcu(3) = rcu(3)/rc



scp = rmu(1)*rcu(1)+rmu(2)*rcu(2)+rmu(3)*rcu(3)
if ( scp .gt.  one ) scp =  one
if ( scp .lt. -one ) scp = -one
f0 = sin ( acos(scp) )
if ( abs(f0) .lt. 1.e-12_prec ) f0 = 1.e-12_prec
f0 = -one / f0
f0 = dv*f0

f1(1) = -2.0_prec*(rcu(1)-rmu(1)*scp)/rm
f1(2) = -2.0_prec*(rcu(2)-rmu(2)*scp)/rm
f1(3) = -2.0_prec*(rcu(3)-rmu(3)*scp)/rm
f3(1) = (rcu(1)-rmu(1)*scp)/rm
f3(2) = (rcu(2)-rmu(2)*scp)/rm
f3(3) = (rcu(3)-rmu(3)*scp)/rm

f2(1) = ( rmu(1)-rcu(1)*scp)/rc
f2(2) = ( rmu(2)-rcu(2)*scp)/rc
f2(3) = ( rmu(3)-rcu(3)*scp)/rc

d(i3+1) = d(i3+1) + f0 * ( f1(1) + f2(1) )
d(i3+2) = d(i3+2) + f0 * ( f1(2) + f2(2) )
d(i3+3) = d(i3+3) + f0 * ( f1(3) + f2(3) )
do isolv = 1, solv_atom - 1
jsolv = i + isolv
j3 = jsolv*3-3
d(j3+1) = d(j3+1) + f0 * ( f3(1) )
d(j3+2) = d(j3+2) + f0 * ( f3(2) )
d(j3+3) = d(j3+3) + f0 * ( f3(3) )
end do
end do

wshell(is)%avtheta = wshell(is)%avtheta + avtdum/real(wshell(is)%n_insh, kind=prec)
wshell(is)%avn_insh = wshell(is)%avn_insh + wshell(is)%n_insh
end do
end subroutine watpol

!----------------------------------------------------------------------------

subroutine write_out
! local variables
integer					::	i,istate

! header line
if(istep >= nsteps) then
write(*,3) 'Energy summary'
else
write(*,2) 'Energy summary', istep
end if
2 format('======================= ',A15,' at step ',i6,' ========================')
3 format('=========================== FINAL ',A15,' =============================')

! legend line
write(*,4) 'el', 'vdW' ,'bond', 'angle', 'torsion', 'improper'
4 format(16X, 6A10)

! row by row: solute, solvent, solute-solvent, LRF, q-atom
write(*,6) 'solute', E%pp%el, E%pp%vdw, E%p%bond, E%p%angle, E%p%torsion, E%p%improper
6 format(A,T17, 6F12.2)

if(nwat > 0) then
write(*,6) 'solvent', E%ww%el, E%ww%vdw, E%w%bond, E%w%angle, E%w%torsion, E%w%improper
end if

write(*,6) 'solute-solvent', E%pw%el, E%pw%vdw

if(use_LRF) then
write(*,6) 'LRF', E%LRF
end if

if(nqat .gt. 0) then
	write(*,6) 'Q-atom', E%qx%el, E%qx%vdw, E%q%bond, E%q%angle, E%q%torsion, E%q%improper
end if

! restraints
write(*,*)
write(*,4) 'total', 'fix', 'slvnt_rad', 'slvnt_pol', 'shell', 'solute'
write(*,6) 'restraints', E%restraint%total, E%restraint%fix, &
        E%restraint%solvent_radial, E%restraint%water_pol, E%restraint%shell, &
        E%restraint%protein
write(*,*)

! totals
if(force_rms) then
        grms = sqrt(dot_product(d(:), d(:))/(3*natom))
        write(*,4) 'total', 'potential', 'kinetic', '', 'RMS force'
        write(*,14) 'SUM', E%potential+E%kinetic, E%potential, E%kinetic, grms
else
        write(*,4) 'total', 'potential', 'kinetic'
        write(*,6) 'SUM', E%potential+E%kinetic, E%potential, E%kinetic
end if
14 format(A,T17, 3F10.2, 10X, F10.2)

! q-atom energies
if(nstates > 0) then
if(istep >= nsteps) then
  write(*,3) 'Q-atom energies'
else
  write(*,2) 'Q-atom energies', istep
end if

write(*,26) 'el', 'vdW' ,'bond', 'angle', 'torsion', 'improper'

do istate =1, nstates
  write (*,32) 'Q-Q', istate, EQ(istate)%lambda, EQ(istate)%qq(1)%el, EQ(istate)%qq(1)%vdw
end do
write(*,*)
if(nat_solute > nqat) then !only if there is something else than Q-atoms in topology
  do istate =1, nstates
        write (*,32) 'Q-prot', istate,EQ(istate)%lambda, EQ(istate)%qp(1)%el, EQ(istate)%qp(1)%vdw
  end do
  write(*,*)
end if

if(nwat > 0) then
  do istate =1, nstates
        write (*,32) 'Q-wat', istate, EQ(istate)%lambda, EQ(istate)%qw(1)%el, EQ(istate)%qw(1)%vdw
  end do
  write(*,*)
end if

do istate =1, nstates
  write (*,32) 'Q-surr.',istate, EQ(istate)%lambda, &
                EQ(istate)%qp(1)%el + EQ(istate)%qw(1)%el, EQ(istate)%qp(1)%vdw &
                + EQ(istate)%qw(1)%vdw
end do
write(*,*)

do istate = 1, nstates
  write (*,36) 'Q-any', istate, EQ(istate)%lambda, EQ(istate)%qx(1)%el,&
                EQ(istate)%qx(1)%vdw, EQ(istate)%q%bond, EQ(istate)%q%angle,&
                EQ(istate)%q%torsion, EQ(istate)%q%improper
end do
write(*,*)

write(*,22) 'total', 'restraint'
do istate = 1, nstates
  write (*,32) 'Q-SUM', istate, EQ(istate)%lambda,&
                EQ(istate)%total(1), EQ(istate)%restraint
end do
do i=1,noffd
  write (*,360) offd(i)%i, offd(i)%j, Hij(offd(i)%i, offd(i)%j), &
                offd2(i)%k, offd2(i)%l, offd(i)%rkl
360	  format ('H(',i2,',',i2,') =',f8.2,' dist. between Q-atoms',2i4, ' =',f8.2)
end do
end if

if(monitor_group_pairs > 0) then
        call centered_heading('Monitoring selected groups of nonbonded interactions','=')
        write (*,37,advance='no')
        write (*,38) (istate,istate, istate=1,nstates)
        do i=1,monitor_group_pairs
                write (*,39,advance='no') i,monitor_group_pair(i)%Vwsum, &
                        monitor_group_pair(i)%Vwel,monitor_group_pair(i)%Vwlj
                write (*,40) (monitor_group_pair(i)%Vel(istate), &
                        monitor_group_pair(i)%Vlj(istate), istate=1,nstates)
        end do
end if

write(*,'(80a)') '==============================================================================='


22	format('type   st lambda',2A10)
26	format('type   st lambda',6a10)
32	format (a,T8,i2,f7.4,2f10.2)
36	format (a,T8,i2,f7.4,6f10.2)
37  format ('pair   Vwsum    Vwel    Vwvdw')
38  format (3(i4,':Vel',i3,':Vvdw'))
39  format (i2,f10.2,f8.2,f9.2)
40  format (3(2f8.2))


if(use_PBC .and. constant_pressure .and. istep>=nsteps ) then
        write(*,*)
        write(*,'(a)') '=========================== VOLUME CHANGE SUMMARY ==========================='
        write(*,45) boxlength(1)*boxlength(2)*boxlength(3)
        write(*,*)
        write(*,46) 'total', 'accepted', 'ratio'
        write(*,47) 'Attempts', volume_try, volume_acc, real(volume_acc, kind=prec)/volume_try
write(*,'(80a)') '==============================================================================='
end if
45 format('Final volume: ', f10.3)
46 format(16X, 3A10)
47 format(A,T17, 2i10, f10.3)

end subroutine write_out

!-----------------------------------------------------------------------

subroutine write_trj

if(.not. trj_write(x)) then
        call die('failure to write to trajectory file')
end if

end subroutine write_trj

!-----------------------------------------------------------------------

subroutine write_xfin
! local variables
integer						::	i,nat3
integer        :: canary = -1337
nat3 = natom*3

if (prec .eq. singleprecision) then
canary = -137
elseif (prec .eq. doubleprecision) then
canary = -1337
#ifndef PGI
elseif (prec .eq. quadprecision) then
canary = -13337
#endif
else
call die('No such precision')
end if
rewind (3)
!new canary on top of file
write (3) canary
write (3) nat3, (x(i),i=1,nat3)
write (3) nat3, (v(i),i=1,nat3)
!save dynamic polarisation restraint data
	if(wpol_restr .and. allocated(wshell)) then
        write (3) nwpolr_shell, wshell(:)%theta_corr
end if

if( use_PBC )then
        write(3) boxlength(:)
        write(3) boxcentre(:)
end if
end subroutine write_xfin

!-----------------------------------------------------------------------
!Put molecules back in box for nice visualisation.
!Change boxcentre if rigid_box_centre is off.
!Update cgp_centers for LRF.
!-----------------------------------------------------------------------
subroutine put_back_in_box

real(kind=prec)				::	boxc(1:3),old_boxc(1:3)
integer				::	i, j, starten, slutet
!the borders of the periodic box
real(kind=prec)				::	x_max, x_min, y_max, y_min, z_max, z_min
real(kind=prec)				:: cm(1:3)
integer				::  k, ig
integer				::  pbib_start, pbib_stop

!the function can be run exclusivly on node 0
!no gain from doing it on each node
if (nodeid.eq.0) then
old_boxc(:)=boxcentre(:)
if( .not. rigid_box_centre ) then !if the box is allowed to float around, center on solute if present, otherwise center on solvent
        if( nat_solute > 0) then
                slutet = nat_solute
                !starten = ncgp_solute + 1
                
        else !if no solute present, centre box around solvent
                slutet = natom
                !starten = 1
        end if
        !find center
        boxc(:) = zero
        do i = 1,slutet
                boxc(1) = boxc(1) + x( 3*i-2 )
                boxc(2) = boxc(2) + x( 3*i-1 )
                boxc(3) = boxc(3) + x( 3*i   )
        end do
        boxc(:) = boxc(:)/slutet
        boxcentre(:) = boxc(:) !store new boxcentre
        ! starten = ncgp_solute + 1 !solute can not move around the corner
else
        !use boxcenter given in topology, ie. 'moving' solute
        boxc(:) = boxcentre(:)
        !starten = 1
end if

!calculate the borders of the periodic box
x_max = boxc(1) + boxlength(1)/2
x_min = boxc(1) - boxlength(1)/2
y_max = boxc(2) + boxlength(2)/2
y_min = boxc(2) - boxlength(2)/2
z_max = boxc(3) + boxlength(3)/2
z_min = boxc(3) - boxlength(3)/2

mvd_mol(:) = 0 

!pbib_start and pbib_stop are the starting and stopping molecule indexes of which molecules to Put Back In Box
pbib_start = 1
pbib_stop  = nmol
if ( .not. put_solute_back_in_box ) then !we're not putting solute back in box
!we now have nsolute to keep track of this stuff, don't use this old way
	pbib_start = nsolute + 1
end if
if ( .not. put_solvent_back_in_box ) then !we're not putting solvent back in box
!same here, we know the molecule numbers
	pbib_stop = nsolute
end if          


do i=pbib_start,pbib_stop
        cm(:) =zero
        do j = istart_mol(i),istart_mol(i+1)-1 !loop over all atoms in molecule i
                cm(1) = cm(1) + x(j*3-2)*mass(j)
                cm(2) = cm(2) + x(j*3-1)*mass(j)
                cm(3) = cm(3) + x(j*3  )*mass(j)
        end do
        cm(:) = cm(:) * mol_mass(i) !centre of mass of molecule i
!mol_mass is actually 1/mol_mass
        !x-direction
        if( cm(1) .gt. x_max) then !position of centre of mass

                        do j= istart_mol(i),istart_mol(i+1)-1 !move the molecule
						x(j*3-2) = x(j*3-2) - boxlength(1)
				        end do
						mvd_mol(i) = 1
		else if ( cm(1) .lt. x_min ) then
                        do j= istart_mol(i),istart_mol(i+1)-1 !move the molecule
						x(j*3-2) = x(j*3-2) + boxlength(1)
				        end do
						mvd_mol(i) = 1
        end if

       ! y-direction
        if( cm(2) .gt. y_max) then !position of centre of mass
                        do j= istart_mol(i),istart_mol(i+1)-1 !move the molecule
						x(j*3-1) = x(j*3-1) - boxlength(2)
				        end do
						mvd_mol(i) = 1
        else if ( cm(2) .lt. y_min ) then
                        do j= istart_mol(i),istart_mol(i+1)-1 !move the molecule
						x(j*3-1) = x(j*3-1) + boxlength(2)
				        end do
						mvd_mol(i) = 1
        end if

        !z-direction
        if( cm(3) .gt. z_max) then !position of centre of mass
                        do j= istart_mol(i),istart_mol(i+1)-1 !move the molecule
						x(j*3  ) = x(j*3  ) - boxlength(3)
				        end do
						mvd_mol(i) = 1
        else if ( cm(3) .lt. z_min ) then
                        do j= istart_mol(i),istart_mol(i+1)-1 !move the molecule
						x(j*3  ) = x(j*3  ) + boxlength(3)
				        end do
						mvd_mol(i) = 1
  		end if
end do !over molecules
end if !if(nodeid .eq. 0)
!now that node 0 has done the job, broadcast the changes to the slaves
#if defined(USE_MPI)
!only during MD, skip for first setup, mpitypes are 
!set after this, so save to use as proxy
if(mpitype_batch_lrf.ne.0) then
	call MPI_Bcast(x, nat3, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
	if (ierr .ne. 0) call die('init_nodes/MPI_Bcast x')
	call MPI_Bcast(boxcentre, 3, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
	if (ierr .ne. 0) call die('init_nodes/MPI_Bcast x')
	call MPI_Bcast(old_boxc, 3, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
	if (ierr .ne. 0) call die('init_nodes/MPI_Bcast x')
end if
#endif

#ifdef USE_GRID
! we now need to update the coordinates of all grid boxes
! to shift by the amount the new box has shifted
grid_pp(:)%x = grid_pp(:)%x - (old_boxc(1)-boxcentre(1))
grid_pp(:)%y = grid_pp(:)%y - (old_boxc(2)-boxcentre(2))
grid_pp(:)%z = grid_pp(:)%z - (old_boxc(3)-boxcentre(3))

grid_pp(:)%xend = grid_pp(:)%xend - (old_boxc(1)-boxcentre(1))
grid_pp(:)%yend = grid_pp(:)%yend - (old_boxc(2)-boxcentre(2))
grid_pp(:)%zend = grid_pp(:)%zend - (old_boxc(3)-boxcentre(3))

grid_pw(:)%x = grid_pw(:)%x - (old_boxc(1)-boxcentre(1))
grid_pw(:)%y = grid_pw(:)%y - (old_boxc(2)-boxcentre(2))
grid_pw(:)%z = grid_pw(:)%z - (old_boxc(3)-boxcentre(3))

grid_pw(:)%xend = grid_pw(:)%xend - (old_boxc(1)-boxcentre(1))
grid_pw(:)%yend = grid_pw(:)%yend - (old_boxc(2)-boxcentre(2))
grid_pw(:)%zend = grid_pw(:)%zend - (old_boxc(3)-boxcentre(3))

grid_ww(:)%x = grid_ww(:)%x - (old_boxc(1)-boxcentre(1))
grid_ww(:)%y = grid_ww(:)%y - (old_boxc(2)-boxcentre(2))
grid_ww(:)%z = grid_ww(:)%z - (old_boxc(3)-boxcentre(3))

grid_ww(:)%xend = grid_ww(:)%xend - (old_boxc(1)-boxcentre(1))
grid_ww(:)%yend = grid_ww(:)%yend - (old_boxc(2)-boxcentre(2))
grid_ww(:)%zend = grid_ww(:)%zend - (old_boxc(3)-boxcentre(3))

#endif
!LRF: if molecule moved update all cgp_centers from first charge group to the last one
if (use_LRF) then

!Broadcast mvd_mol(:) & x(:)
!#if defined(USE_MPI)
!call MPI_Bcast(mvd_mol, nmol, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
!if (ierr .ne. 0) call die('put_back_in_box MPI_BCast mvd_mol')
!!broadcast start and stop molecule so that each node can do its job
!call MPI_Bcast(pbib_start, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
!if (ierr .ne. 0) call die('put_back_in_box MPI_Bcast pbib_start')
!call MPI_Bcast(pbib_stop, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
!if (ierr .ne. 0) call die('put_back_in_box MPI_Bcast pbib_stop')
!#endif

if (nodeid .eq.0) then
do k=pbib_start,pbib_stop
		if (mvd_mol(k) == 1) then  
						do ig=iwhich_cgp(istart_mol(k)),iwhich_cgp(istart_mol(k+1)-1)
								lrf(ig)%cgp_cent(:) = 0
								do i  = cgp(ig)%first, cgp(ig)%last
										lrf(ig)%cgp_cent(:) = lrf(ig)%cgp_cent(:) + x(cgpatom(i)*3-2:cgpatom(i)*3)
								end do
						lrf(ig)%cgp_cent(:) = lrf(ig)%cgp_cent(:)/real(cgp(ig)%last - cgp(ig)%first +1, kind=prec)
						end do
		end if 
end do
end if ! nodeid .eq. 0
#if defined(USE_MPI)
!for now just bcast the whole thing
!only during MD, skip for first setup, mpitypes are 
!set after this, so save to use as proxy
if(mpitype_batch_lrf.ne.0) then
	call MPI_Bcast(lrf,ncgp,mpitype_batch_lrf,0,MPI_COMM_WORLD, ierr)
	if (ierr .ne. 0) call die('put_back_in_box MPI_Bcast lrf')
end if
! for later use of PBC update on each node
!call lrf_gather
#endif
end if

end subroutine put_back_in_box
!----------------------------------------------------------------------------
subroutine MC_volume()

real(kind=prec)									:: old_x(1:3*nat_pro), old_xx(1:3*nat_pro), x_move(1:3*nat_pro)
real(kind=prec)									:: old_boxl(1:3), old_inv(1:3)
real(kind=prec)									:: old_V, new_V, deltaLength
real(kind=prec)									:: deltaV, deltaE, deltaW
real(kind=prec)									:: box_min
integer									:: starten, slutet, i, j, sw_no !indeces
real(kind=prec)									:: randomno !random number
real(kind=prec)									:: new_x, new_y, new_z
real(kind=prec)									:: move_x, move_y, move_z
logical									:: acc
integer									:: longest , niter
real(kind=prec)									:: cubr_vol_ratio
real(kind=prec)									:: cm(1:3)
real(kind=prec)									::	old_EMorseD(max_qat)
real(kind=prec)									::	old_dMorse_i(3,max_qat)
real(kind=prec)									::	old_dMorse_j(3,max_qat)
integer									:: old_nbww_pair, old_nbpw_pair , old_nbpw_cgp_pair 
integer									:: old_nbpp_pair , old_nbpp_cgp_pair , old_nbqp_pair , old_nbqp_cgp_pair
integer									:: old_ww_max, old_pw_max, old_pp_max, old_qp_max

if (nodeid .eq. 0) then
write(*,8) 'Volume change', istep
write(*,*)
write(*,'(a)') '---------- Before move'
8 format('======================== ',A14,' at step ',i6,' ========================')
4 format(16X, 3A10)
6 format(A,T17, 3F10.3)
end if  !(nodeid .eq. 0)

!save the old energies,coordinates and forces
old_x(:) = x(:)
previous_E = E
old_boxl(:) = boxlength(:)
old_inv(:) = inv_boxl(:)

!qatom stuff
if(nqat.gt.0) then
        old_EMorseD = EMorseD
        old_dMorse_i = dMorse_i
        old_dMorse_j = dMorse_j
end if
!shake stuff
old_xx(:) = xx(:)



if (use_LRF) then
	old_nbww_pair = nbww_pair
	old_nbpw_pair = nbpw_pair
	old_nbpp_pair = nbpp_pair
	old_nbqp_pair = nbqp_pair

	old_nbpw_cgp_pair = nbpw_cgp_pair
	old_nbpp_cgp_pair = nbpp_cgp_pair
	old_nbqp_cgp_pair = nbqp_cgp_pair
!now only reallocate if the old size was too small
	if(calculation_assignment%ww%max.gt.SIZE(old_nbww,1)) then
		if (allocated(old_nbww)) deallocate(old_nbww)
		allocate(old_nbww(calculation_assignment%ww%max))
	endif
        if(calculation_assignment%pw%max.gt.SIZE(old_nbpw,1)) then
                if (allocated(old_nbpw)) deallocate(old_nbpw)
                allocate(old_nbpw(calculation_assignment%pw%max))
        endif
        if(calculation_assignment%pp%max.gt.SIZE(old_nbpp,1)) then
                if (allocated(old_nbpp)) deallocate(old_nbpp)
                allocate(old_nbpp(calculation_assignment%pp%max))
        endif
        if(calculation_assignment%qp%max.gt.SIZE(old_nbqp,1)) then
                if (allocated(old_nbqp)) deallocate(old_nbqp)
                allocate(old_nbqp(calculation_assignment%qp%max,nstates))
        endif
	old_nbww = nbww
	old_nbpw = nbpw
	old_nbpp = nbpp
	old_nbqp = nbqp
	old_nbpw_cgp = nbpw_cgp
	old_nbpp_cgp = nbpp_cgp
	old_nbqp_cgp = nbqp_cgp
	old_ww_max = calculation_assignment%ww%max
	old_pw_max = calculation_assignment%pw%max
	old_pp_max = calculation_assignment%pp%max
	old_qp_max = calculation_assignment%qp%max
	old_lrf(:) = lrf(:)
end if

#if defined(USE_MPI)
!Update modified coordinates  
call MPI_Bcast(x, natom*3, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast x')
#endif

call new_potential(previous_E)   !compute energies from previous md-step

if (nodeid .eq. 0 ) then
	old_E = E                !Update to fresh E before changing volume
        if (nqat.gt.0)	old_EQ(1:nstates) = EQ(1:nstates)
	old_V = old_boxl(1) * old_boxl(2) * old_boxl(3)


	!new volume randomized
	randomno = randm(pressure_seed) ! 0<=randomno<=1
	randomno = randomno*2 - 1    !-1 <= randomno <= 1
	deltaV = randomno * max_vol_displ
	new_V = deltaV + old_V
	cubr_vol_ratio = (new_V/old_V)**(one/3.0_prec)
	write(*,4) 'old', 'new', 'delta'
17	format('Volume',8x,f10.3,2x,f10.3,2x,f10.3)
	write(*,17) old_V, new_V, deltaV
	write(*,*)

	!compute new boxlenth and inv_boxl
	boxlength(1) = boxlength(1)*cubr_vol_ratio
	boxlength(2) = boxlength(2)*cubr_vol_ratio
	boxlength(3) = boxlength(3)*cubr_vol_ratio
	inv_boxl(:) = one/boxlength(:)
	write(*,10) old_boxl
	write(*,2) boxlength
	write(*,*)
	10 format('Old boxlength', 3f10.3)
	2 format('New boxlength ', 3f10.3)

	!compare cut-offs with new boxsize
	box_min = min( boxlength(1), boxlength(2), boxlength(3) )
	!Solute-solute cut-off radii
	if( .not. (box_min .gt. Rcpp*2) ) then
		write(*,*) 'Solute-solute cut-off radii too large', Rcpp
		call die('Solute-solute cut-off radii too large')
	!Solvent-solvent
	else if( .not. (box_min .gt. Rcww*2) ) then
		write(*,*) 'Solvent-solvent cut-off radii too large', Rcww
		call die('Solvent-solvent cut-off radii too large')
	!Solute-solvent
	else if( .not. (box_min .gt. Rcpw*2) ) then
		write(*,*) 'Solute-solvent cut-off radii too large', Rcpw
		call die('Solute-solvent cut-off radii too large')
	!Q-atom
	else if( .not. (box_min .gt. Rcq*2) ) then
		write(*,*) 'Q-atom cut-off radii too large', Rcq
		call die('Q-atom cut-off radii too large')
	!LRF
	else if( .not. (box_min .gt. RcLRF*2) ) then
		write(*,*) 'LRF cut-off radii too large', Rcq
		call die('LRF cut-off radii too large')
	end if




	if (.not. atom_based_scaling) then


		!compute new coordinates after molecules and centre of mass
		do i=1,nmol !looping over all molecules, can also do the last because of nmol +1 ...
			cm(:) =zero
			do j = istart_mol(i),istart_mol(i+1)-1 !loop over all atoms in molecule i
				cm(1) = cm(1) + x(j*3-2)*mass(j)
				cm(2) = cm(2) + x(j*3-1)*mass(j)
				cm(3) = cm(3) + x(j*3  )*mass(j)
			end do
!but mol_mass is actually 1/mol_mass
			cm(:) = cm(:) * mol_mass(i) !centre of mass of molecule i

			move_x = ( ( cm(1)-boxcentre(1) )*boxlength(1)/old_boxl(1) + boxcentre(1) ) - cm(1)
			move_y = ( ( cm(2)-boxcentre(2) )*boxlength(2)/old_boxl(2) + boxcentre(2) ) - cm(2)
			move_z = ( ( cm(3)-boxcentre(3) )*boxlength(3)/old_boxl(3) + boxcentre(3) ) - cm(3)

			do j= istart_mol(i),istart_mol(i+1)-1 !move the molecule
				x(j*3-2) = x(j*3-2) + move_x
				x(j*3-1) = x(j*3-1) + move_y
				x(j*3  ) = x(j*3  ) + move_z
			end do

		end do !over molecules

!we have nmol+1 to be able to iterate over all molecules

	else ! atom_based_scaling = .true.
	!move xx also if shake
		do j=1,natom
			x_move(j*3-2) = ( ( x(j*3-2)-boxcentre(1) )*boxlength(1)/old_boxl(1) + boxcentre(1) ) - x(j*3-2)
			x_move(j*3-1) = ( ( x(j*3-1)-boxcentre(2) )*boxlength(1)/old_boxl(2) + boxcentre(2) ) - x(j*3-1)
			x_move(j*3  ) = ( ( x(j*3  )-boxcentre(3) )*boxlength(1)/old_boxl(3) + boxcentre(3) ) - x(j*3  )
		end do
		!can be written: ?
		! x_move (1:3*nat_pro-2:3) = ( ( x(1:3*nat_pro-2:3) - boxcentre(1) )*boxlength(1)/old_boxl(1) + boxcentre(1) ) - x(1:3*nat_pro-2:3)
		! x_move (2:3*nat_pro-1:3) = ( ( x(2:3*nat_pro-1:3) - boxcentre(2) )*boxlength(2)/old_boxl(2) + boxcentre(2) ) - x(1:3*nat_pro-1:3)
		! x_move (3:3*nat_pro-0:3) = ( ( x(3:3*nat_pro-0:3) - boxcentre(3) )*boxlength(3)/old_boxl(3) + boxcentre(3) ) - x(1:3*nat_pro-0:3)

		x(:) = x(:) + x_move(:)
		
		! shake if necessary
		if(shake_constraints > 0) then
			xx(:) = xx(:) + x_move(:)
			niter=shake(xx, x)
		end if

	end if
end if

#if defined(USE_MPI)
!Update modified coordinates and boxlengths 
call MPI_Bcast(x, natom*3, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast x')
call MPI_Bcast(boxlength, 3, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
call MPI_Bcast(inv_boxl, 3, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
#endif

!Need to update entire LRF... sigh
if (use_LRF) then
	call cgp_centers
	if ( iuse_switch_atom == 1 ) then 
		call nbpplist_box_lrf
		call nbpwlist_box_lrf
		call nbqplist_box
	else
		call nbpplis2_box_lrf
		call nbpwlis2_box_lrf
		call nbqplis2_box
	endif
	call nbwwlist_box_lrf
	call nbqwlist_box
! we are gathering the LRF info here
#ifdef USE_MPI
call lrf_gather
#endif
end if !use_LRF



!compute the new potential, in parallel if possible
call new_potential( previous_E ) !do we need a broadcast before this if using LRF (i.e. if the nonbond lists have been updated)

if (nodeid .eq. 0) then
	!Jamfor nya med gamla
	deltaE = E%potential - old_E%potential
	deltaW = deltaE + pressure * deltaV - nmol*Boltz*Temp0*log(new_V/old_V)
	write(*,4) 'old', 'new', 'delta'
18	format('Potential',7x,f14.3,2x,f14.3,2x,f14.3)
	write(*,18) old_E%potential, E%potential, deltaE
	write(*,*)

	!accept or reject
	if( deltaW<=zero ) then
	acc = .true.
	else
	!slumpa tal mellan  0 coh 1
	randomno = randm(pressure_seed)
	if( randomno > exp(- deltaW / Boltz / Temp0) ) then
		acc = .false.
	else
		acc = .true.
	end if
	end if

	volume_try = volume_try + 1
	write(*,'(a)') '---------- After move' 

if( acc ) then
		write(*,'(a)') 'Volume change accepted'
		volume_acc = volume_acc + 1
else
		write(*,'(a)') 'Volume change rejected'

        !put stuff back to what they were before
        x(:) = old_x(:)
        E = previous_E
        EQ(1:nstates) = old_EQ(1:nstates)
        boxlength(:) = old_boxl(:)
        inv_boxl(:) = old_inv(:)
        if(nqat.gt.0) then
                EQ(1:nstates) = old_EQ(1:nstates)
                EMorseD = old_EMorseD
        	dMorse_i = old_dMorse_i
        	dMorse_j = old_dMorse_j
        end if
	xx(:) = old_xx(:)

	if (use_LRF) then
		nbww_pair = old_nbww_pair
		nbpw_pair = old_nbpw_pair
		nbpp_pair = old_nbpp_pair
		nbqp_pair = old_nbqp_pair

		nbpw_cgp_pair = old_nbpw_cgp_pair
		nbpp_cgp_pair = old_nbpp_cgp_pair
		nbqp_cgp_pair = old_nbqp_cgp_pair
		calculation_assignment%ww%max = old_ww_max
		calculation_assignment%pw%max = old_pw_max
		calculation_assignment%pp%max = old_pp_max
		calculation_assignment%qp%max = old_qp_max
        if(calculation_assignment%ww%max.gt.SIZE(nbww,1)) then
                call reallocate_nonbondlist_ww
        endif
        if(calculation_assignment%pw%max.gt.SIZE(nbpw,1)) then
                call reallocate_nonbondlist_pw
        endif
        if(calculation_assignment%pp%max.gt.SIZE(nbpp,1)) then
                call reallocate_nonbondlist_pp
        endif
        if(calculation_assignment%qp%max.gt.SIZE(nbqp,1)) then
                call reallocate_nonbondlist_qp
        endif
		nbww(1:old_ww_max) = old_nbww(1:old_ww_max)
		nbpw(1:old_pw_max) = old_nbpw(1:old_pw_max)
		nbpp(1:old_pp_max) = old_nbpp(1:old_pp_max)
		nbqp(1:old_qp_max,1:nstates) = old_nbqp(1:old_qp_max,1:nstates)
		nbpw_cgp(:) = old_nbpw_cgp(:)
		nbpp_cgp(:) = old_nbpp_cgp(:)
		nbqp_cgp(:) = old_nbqp_cgp(:)
		lrf(:) = old_lrf(:)
	end if
end if

	write(*,11) boxlength(1)*boxlength(2)*boxlength(3)
	write(*,12) boxlength
	write(*,13) sum(mass(:))/(boxlength(1)*boxlength(2)*boxlength(3)*1E-24_prec*6.02E23_prec)
	11 format('Final volume: ', f10.3)
	12 format('Final boxlength: ', 3f10.3)
	13 format('Final density (g/cm^3): ', f10.3)
	write(*,*)
	write(*,4) 'total', 'accepted', 'ratio'
	write(*,7) 'Attempts', volume_try, volume_acc, real(volume_acc, kind=prec)/volume_try
	write(*,*)
	7 format(A,T17, 2i10, f10.3)
	write(*,'(80a)') '==============================================================================='
end if !(nodeid .eq. 0)
#if defined(USE_MPI)
!Make slave nodes put things back if rejected
call MPI_Bcast(acc, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
if (.not. acc) then
  if (nodeid .ne. 0) then
    x(:) = old_x(:)
    boxlength(:) = old_boxl(:)
    inv_boxl(:) = old_inv(:)
  end if
end if
#endif
end subroutine MC_volume	

!-----------------------------------------------------------------------------------------------------
subroutine new_potential( old )

type(ENERGIES), intent(in)		:: old	
integer							:: istate,i, numrequest,ii,jj

!zero all energies
E%potential = zero
E%pp%el  = zero
E%pp%vdw = zero
E%pw%el  = zero
E%pw%vdw = zero
E%ww%el  = zero
E%ww%vdw = zero
E%qx%el    = zero
E%qx%vdw   = zero
E%restraint%protein = zero
E%LRF      = zero
E%p%bond =  zero
E%w%bond =  zero
E%p%angle =  zero
E%w%angle =  zero
E%p%angle =  zero
E%w%angle =  zero
E%p%torsion =  zero
E%w%torsion =  zero
E%p%improper =  zero
E%w%improper =  zero

do istate = 1, nstates
	EQ(istate)%qq(1:ene_header%arrays)%el = zero
	EQ(istate)%qq(1:ene_header%arrays)%vdw = zero
	EQ(istate)%qp(1:ene_header%arrays)%el = zero
	EQ(istate)%qp(1:ene_header%arrays)%vdw = zero
	EQ(istate)%qw(1:ene_header%arrays)%el = zero
	EQ(istate)%qw(1:ene_header%arrays)%vdw = zero
	EQ(istate)%restraint = zero
end do

!reset derivatives ---
d(:) = zero

#if defined (USE_MPI)
!First post recieves for gathering data from slaves
if (nodeid .eq. 0) then
	call gather_nonbond
end if
#endif


if (nodeid .eq. 0) then
        call pot_energy_bonds
        call p_restrain
end if
if(natom > nat_solute) then
        if((ivdw_rule.eq.VDW_GEOMETRIC).and. &
                (solvent_type == SOLVENT_SPC)) then
                call nonbond_qw_spc_box
                call nonbond_ww_spc_box
        else
                call nonbond_qw_box
                call nonbond_ww_box
        end if
        if(ntors>ntors_solute) then
                call nonbond_solvent_internal_box
        end if
        call nonbond_pw_box
end if
call nonbond_pp_box
call nonbond_qp_box
if (nodeid .eq. 0) then
        call nonbond_qqp
        call nonbond_qq
end if

if (use_LRF) then
        call lrf_taylor
end if

#if defined(USE_MPI)
if (nodeid .ne. 0) then  !Slave nodes
	call gather_nonbond
end if
#endif

if (nodeid .eq. 0) then 
#if (USE_MPI)
do i = 1, 3
    call MPI_WaitAll((numnodes-1),request_recv(1,i),mpi_status,ierr)
end do

	!Forces and energies are summarised
	do i=1,numnodes-1
	  d = d + d_recv(:,i)
	  E%pp%el   = E%pp%el  + E_recv(i)%pp%el
	  E%pp%vdw  = E%pp%vdw + E_recv(i)%pp%vdw
	  E%pw%el   = E%pw%el  + E_recv(i)%pw%el
	  E%pw%vdw  = E%pw%vdw + E_recv(i)%pw%vdw
	  E%ww%el   = E%ww%el  + E_recv(i)%ww%el
	  E%ww%vdw  = E%ww%vdw + E_recv(i)%ww%vdw
	  E%lrf     = E%lrf    + E_recv(i)%lrf
		do ii=1,nstates
		do jj=1,ene_header%arrays
	  EQ(ii)%qp(jj)%el  = EQ(ii)%qp(jj)%el  + EQ_recv(ii,jj,i)%qp%el
	  EQ(ii)%qp(jj)%vdw = EQ(ii)%qp(jj)%vdw + EQ_recv(ii,jj,i)%qp%vdw
	  EQ(ii)%qw(jj)%el  = EQ(ii)%qw(jj)%el  + EQ_recv(ii,jj,i)%qw%el
	  EQ(ii)%qw(jj)%vdw = EQ(ii)%qw(jj)%vdw + EQ_recv(ii,jj,i)%qw%vdw
		end do
		end do
	end do
#endif

	!summation of energies
	do istate = 1, nstates
			! update EQ
		do jj=1,ene_header%arrays
			EQ(istate)%qx(jj)%el  = EQ(istate)%qq(jj)%el &
			+ EQ(istate)%qp(jj)%el +EQ(istate)%qw(jj)%el
			EQ(istate)%qx(jj)%vdw = EQ(istate)%qq(jj)%vdw &
			+ EQ(istate)%qp(jj)%vdw+EQ(istate)%qw(jj)%vdw
		end do
		E%qx%el    = E%qx%el    + EQ(istate)%qx(1)%el   *EQ(istate)%lambda
		E%qx%vdw   = E%qx%vdw   + EQ(istate)%qx(1)%vdw  *EQ(istate)%lambda

		! update E%restraint%protein with an average of all states
		E%restraint%protein = E%restraint%protein + EQ(istate)%restraint*EQ(istate)%lambda
	end do


E%potential = old%p%bond + old%w%bond + old%p%angle + old%w%angle + old%p%torsion + &
old%p%improper + E%pp%el + E%pp%vdw + E%pw%el + E%pw%vdw + E%ww%el + &
E%ww%vdw + old%q%bond + old%q%angle + old%q%torsion + &
old%q%improper + E%qx%el + E%qx%vdw + E%restraint%protein + E%LRF


end if !(nodeid .eq. 0)
end subroutine new_potential

!----------------------------------------------------------------------------------------

#if defined (USE_MPI)
!***********************
!subroutine handling summation of nonbonded energies from slave nodes.
!***********************
! Use the global vars
!  request_recv, E_send,EQ_send,E_recv,EQ_Recv,d_recv
! Allocate  - status


subroutine gather_nonbond()

integer,parameter                       :: vars=3
integer,dimension(3,numnodes-1)         :: tag
integer,dimension(vars)	                :: blockcnt,ftype 
integer(kind=MPI_ADDRESS_KIND), dimension(vars)	:: fdisp, base
integer                                 :: mpitype_package,mpitype_send
integer                                 :: i,j,istate,ii

do i=1,numnodes-1
tag(1,i)=numnodes*100+i
tag(2,i)=numnodes*200+i
tag(3,i)=numnodes*300+i
end do

if (nodeid .eq. 0) then        !master

!reclength is defined in mpiglob and set in prep_sim
!used for the now not previously known length of the EQ array
!for sending same length distributed to nodes

! Post receives for each of the d/E/EQ_recv structures
! E/EQ_Recv should really be handled with MPI_Type_create_struct
! and d_recv's type should be handled correctly (it's KIND=wp8)
! should preferably use size(d_recv, 1) for count
do i = 1,numnodes-1
  call MPI_IRecv(d_recv(1,i), natom*3, MPI_REAL8, i, tag(1,i), MPI_COMM_WORLD, &
       request_recv(i,1),ierr)
  if (ierr .ne. 0) call die('gather_nonbond/MPI_IRecv d_recv')
  call MPI_IRecv(E_recv(i), 3*2+1, MPI_REAL8, i, tag(2,i), MPI_COMM_WORLD, &
       request_recv(i,2),ierr)
  if (ierr .ne. 0) call die('gather_nonbond/MPI_IRecv E_recv')
  call MPI_IRecv(EQ_recv(1,1,i), reclength, MPI_REAL8, i, tag(3,i), MPI_COMM_WORLD, &
	request_recv(i,3),ierr)
  if (ierr .ne. 0) call die('gather_nonbond/MPI_IRecv EQ_recv')
end do

else                  !slave nodes
E_send%pp%el  = E%pp%el
E_send%pp%vdw = E%pp%vdw
E_send%pw%el  = E%pw%el
E_send%pw%vdw = E%pw%vdw
E_send%ww%el  = E%ww%el
E_send%ww%vdw = E%ww%vdw
E_send%lrf    = E%lrf
do ii=1,nstates
do i=1,ene_header%arrays
EQ_send(ii,i)%qp%el  = EQ(ii)%qp(i)%el
EQ_send(ii,i)%qp%vdw = EQ(ii)%qp(i)%vdw
EQ_send(ii,i)%qw%el  = EQ(ii)%qw(i)%el
EQ_send(ii,i)%qw%vdw = EQ(ii)%qw(i)%vdw
end do
end do

! See comments above on the IRecv part
call MPI_Send(d, natom*3, MPI_REAL8, 0, tag(1,nodeid), MPI_COMM_WORLD,ierr) 
if (ierr .ne. 0) call die('gather_nonbond/Send d')
call MPI_Send(E_send, 3*2+1, MPI_REAL8, 0, tag(2,nodeid), MPI_COMM_WORLD,ierr) 
if (ierr .ne. 0) call die('gather_nonbond/Send E_send')
call MPI_Send(EQ_send, reclength, MPI_REAL8, 0, tag(3,nodeid), MPI_COMM_WORLD,ierr) 
if (ierr .ne. 0) call die('gather_nonbond/Send EQ_send')

end if
end subroutine gather_nonbond

#endif
!----------------------------------------------------------------------------------------
!*******************************************************
!Will find and return the xtop atom number from 
!  residue number and atom number in residue from
!  library sequence.
! Uses global variables: xtop,nres,res
!*******************************************************

integer function get_atom_from_resnum_atnum(aid)
!arguments
character(*), intent(in)	::	aid	!string=residue:atom
	
!locals
integer						::	separator_pos
character(len=20)			::	res_str
character(len=5)			::	atom_str
integer						::	filestat
integer						::	resnum, atnum

get_atom_from_resnum_atnum = 0

separator_pos = scan(aid, ':')
if(separator_pos < 2 .or. separator_pos == len_trim(aid)) return !no valid colon found
res_str = aid(1:separator_pos-1)
atom_str = aid(separator_pos+1:len_trim(aid))
read(res_str, *, iostat=filestat) resnum
read(atom_str, *, iostat=filestat) atnum
if(filestat > 0) return

!Residue must be in topology
if(resnum < 1 .or. resnum > nres) then                     
  return                                                 
end if

if(atnum .le. (res(resnum+1)%start - res(resnum)%start)) then
  get_atom_from_resnum_atnum = res(resnum)%start + atnum - 1
return
end if

!we have an error: 
write(*, 120) atnum, resnum
call die('error in finding atom number from resnum:atnum.')

120	format('>>>>> ERROR: There is no atom number ',i4,' in residue ',i4,'.')
end function get_atom_from_resnum_atnum

!----------------------------------------------------------------------------------------

end module MD




