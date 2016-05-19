! (C) 2014 Uppsala Molekylmekaniska HB, Uppsala, Sweden
! md.f90
! by Johan Åqvist, John Marelius, Anders Kaplan, Isabella Feierberg, Martin Nervall & Martin Almlöf
! molecular dynamics
!TODO: remove default real statment and change real(4) - in accordance with best practice

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
character(80)	                ::	MD_VERSION = ''
character(80)	                ::	MD_DATE    = ''
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

!-----------------------------------------------------------------------
!	profiling vars
!-----------------------------------------------------------------------
#if defined (PROFILING)


integer, parameter		:: num_profiling_times = 11

type profiling_var_type
	character(len=100)	:: name
	real(kind=prec)			:: time = zero
end type profiling_var_type

type(profiling_var_type)	:: profile(num_profiling_times)

#if defined (USE_MPI)
 !vectors for keeping track of node times
real(kind=prec),allocatable				:: all_node_times(:)
real(kind=prec),allocatable				:: node_times(:)
#endif

#endif

!-----------------------------------------------------------------------
!	temperature calculation variables
!-----------------------------------------------------------------------
integer						:: Ndegf,Ndegfree
integer						:: Ndegf_solute,Ndegfree_solute
integer						:: Ndegf_solvent,Ndegfree_solvent
logical						:: detail_temps			!controls whether or not solute and solvent temps are printed separately (true if solute and solvent degrees of freedom are both not zero)

!----END OF SHARED VARIABLES

!----START OF PUBLIC SUBROUTINES
contains


!----------------------------------------------------------------------


subroutine md_startup
! initialise used modules
call qatom_startup
call trj_startup


! initialise constants
pi = 4.0_prec*atan(one)
deg2rad = pi/180.0_prec

! read in version info
  MD_VERSION = trim(version_pass())
  MD_DATE    = trim(date_pass())


end subroutine md_startup


!----------------------------------------------------------------------

subroutine md_shutdown
! call used modules' shutdown subroutines
if (use_excluded_groups) then
call excluded_shutdown(ngroups_gc)
end if
call md_deallocate
call topo_deallocate
call qatom_shutdown
call index_shutdown
call trj_shutdown
end subroutine md_shutdown

!----------------------------------------------------------------------

subroutine die(cause)
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
        ! write final energies if run has started
        if (istep > 0) then
                if ( mod(istep,iout_cycle) .ne. 1 ) call write_out
        end if
		if(allocated(v)) then
	        !save restart file for diagnosing coordinate problems
		    write(*,*) 'restart file written at step', istep
			call write_xfin
		endif
        write (*,'(79a)') ('!',i=1,79)
        call close_output_files

        ! apologise
        write(*,'(a)') 'ABNORMAL TERMINATION of Qdyn5'
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
stop 'Qdyn5 terminated abnormally'
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

! --- Dynamics subroutines, alphabetically

real(kind=prec) function angle(istart, iend)
! *** arguments
integer						::	istart, iend

! *** local variables
integer						::	i,j,k,ia,ic,i3,j3,k3
real(kind=prec)						::	bjiinv, bjkinv, bji2inv, bjk2inv
real(kind=prec)						::	scp,angv,da,dv,f1
real(kind=prec)						::  rji(3),rjk(3),di(3),dk(3) 

#ifdef _OPENMP
integer :: quotient, remainder
#endif
! global variables used:
! ang, x, anglib, d

! calculate the total energy of all protein or water angles, depending
! updates d

! reset Eangle
angle = zero			!zero = 0.0_prec

!$omp parallel default(none) shared(threads_num, istart, iend, ang, x, anglib, d, angle) private(remainder, quotient,i,j,k,ia,ic,i3,j3,k3,bjiinv,bjkinv,bji2inv,bjk2inv,scp,angv,da,dv,f1,rji,rjk,di,dk)
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

do ia=mp_start,mp_end

#else
do ia=istart,iend
#endif
        ! for each angle in range:

        i  = ang(ia)%i
        j  = ang(ia)%j
        k  = ang(ia)%k
        ic = ang(ia)%cod
        ! calculate rji and rjk
i3=i*3-3
j3=j*3-3
k3=k*3-3
rji(1) = x(i3+1) - x(j3+1)
rji(2) = x(i3+2) - x(j3+2)
rji(3) = x(i3+3) - x(j3+3)
rjk(1) = x(k3+1) - x(j3+1)
rjk(2) = x(k3+2) - x(j3+2)
rjk(3) = x(k3+3) - x(j3+3)

        ! calculate bjiinv and bjkinv and their squares
bji2inv = one/(rji(1)**2 + rji(2)**2 + rji(3)**2 )
bjk2inv = one/(rjk(1)**2 + rjk(2)**2 + rjk(3)**2 )
        bjiinv = sqrt(bji2inv)
        bjkinv = sqrt(bjk2inv)

        ! calculate scp and angv
scp = ( rji(1)*rjk(1) + rji(2)*rjk(2) + rji(3)*rjk(3) )
scp = scp * bjiinv*bjkinv
if ( scp .gt.  one ) then
          scp =  one
        else if ( scp .lt. -one ) then
          scp = -one
        end if
angv = acos(scp)

        ! calculate da and dv
da = angv - anglib(ic)%ang0
#ifdef _OPENMP
    mp_real_tmp = mp_real_tmp + 0.5_prec*anglib(ic)%fk*da**2
#else 
    angle = angle + 0.5_prec*anglib(ic)%fk*da**2
#endif
dv = anglib(ic)%fk*da

        ! calculate f1
f1 = sin ( angv ) 
if ( abs(f1) .lt. 1.e-12_prec ) then	! E is for single precision adding _prec customise precision D for double but _prec addon could not work propielty.
          ! avoid division by zero
          f1 = -1.e12_prec
        else
  f1 =  -one / f1
        end if

        ! calculate di and dk
di(1) = f1 * ( rjk(1)*bjiinv*bjkinv - scp*rji(1)*bji2inv )
di(2) = f1 * ( rjk(2)*bjiinv*bjkinv - scp*rji(2)*bji2inv )
di(3) = f1 * ( rjk(3)*bjiinv*bjkinv - scp*rji(3)*bji2inv )
dk(1) = f1 * ( rji(1)*bjiinv*bjkinv - scp*rjk(1)*bjk2inv )
dk(2) = f1 * ( rji(2)*bjiinv*bjkinv - scp*rjk(2)*bjk2inv )
dk(3) = f1 * ( rji(3)*bjiinv*bjkinv - scp*rjk(3)*bjk2inv )

        ! update d
d(i3+1) = d(i3+1) + dv*di(1)
d(i3+2) = d(i3+2) + dv*di(2)
d(i3+3) = d(i3+3) + dv*di(3)
d(k3+1) = d(k3+1) + dv*dk(1)
d(k3+2) = d(k3+2) + dv*dk(2)
d(k3+3) = d(k3+3) + dv*dk(3)
d(j3+1) = d(j3+1) - dv*( di(1) + dk(1) )
d(j3+2) = d(j3+2) - dv*( di(2) + dk(2) )
d(j3+3) = d(j3+3) - dv*( di(3) + dk(3) )
end do
#ifdef _OPENMP
!$omp atomic update
angle = angle + mp_real_tmp
#endif
!$omp end parallel

end function angle

!-----------------------------------------------------------------------
real(kind=prec) function urey_bradley(istart, iend)
! *** arguments
integer						::	istart, iend

! *** local variables
integer						::	i,j,k,ia,ic,i3,j3,k3
real(kind=prec)						::	bjiinv, bjkinv, bji2inv, bjk2inv
real(kind=prec)						::	scp,angv,da,dv,f1
real(kind=prec)						::  rji(3),rjk(3),di(3),dk(3) 
real(kind=prec)						::	rik(3), dik, ru, du 
real(kind=prec)						::	Eurey

#ifdef _OPENMP
integer :: quotient, remainder
#endif
! global variables used:
! ang, x, anglib, d

! calculate the total energy of all protein or water angles, depending
! updates d
! reset energy
urey_bradley = zero                    !zero = 0.0_prec

!$omp parallel default(none) shared(threads_num, istart, iend, ang, x, anglib, d, urey_bradley) private(remainder, quotient,i,j,k,ia,ic,i3,j3,k3,bjiinv,bjkinv,bji2inv,bjk2inv,scp,angv,da,dv,f1,rji,rjk,di,dk,rik,dik,ru,du)
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

do ia=mp_start,mp_end

#else
do ia=istart,iend
#endif
        ! for each angle in range:

        i  = ang(ia)%i
        j  = ang(ia)%j
        k  = ang(ia)%k
        ic = ang(ia)%cod
        ! calculate rji and rjk
i3=i*3-3
j3=j*3-3
k3=k*3-3
        ! 1-3 distance for Urey-Bradley potential:
    if(anglib(ic)%ureyfk > zero) then
                rik(1) = x(k3+1) - x(i3+1)
                rik(2) = x(k3+2) - x(i3+2)
                rik(3) = x(k3+3) - x(i3+3)
                dik = sqrt(rik(1)*rik(1) + rik(2)*rik(2) + rik(3)*rik(3))
                ru = dik - anglib(ic)%ureyr0
#ifdef _OPENMP
    mp_real_tmp = mp_real_tmp + anglib(ic)%ureyfk*ru**2
#else 
                urey_bradley = urey_bradley + anglib(ic)%ureyfk*ru**2
#endif
                du = 2*anglib(ic)%ureyfk*ru/dik
                d(k3+1) = d(k3+1) + du*rik(1)
                d(k3+2) = d(k3+2) + du*rik(2)
                d(k3+3) = d(k3+3) + du*rik(3)
                d(i3+1) = d(i3+1) - du*rik(1)
                d(i3+2) = d(i3+2) - du*rik(2)
                d(i3+3) = d(i3+3) - du*rik(3)
        end if
end do
#ifdef _OPENMP
!$omp atomic update
urey_bradley = urey_bradley + mp_real_tmp
#endif
!$omp end parallel

end function urey_bradley

!-----------------------------------------------------------------------

real(kind=prec) function bond(istart, iend)
! *** arguments
integer						::	istart, iend

! *** local variables
integer						::	i,j,ib,ic,i3,j3
real(kind=prec)						::	b,db,dv
real(kind=prec)						::	rij(3)

#ifdef _OPENMP
integer :: quotient, remainder
#endif
! global variables used:
! bnd, x, bondlib, d

! reset Ebond
bond = zero

!$omp parallel default(none) shared(threads_num, istart, iend, bnd, x, bondlib, d, bond) private(remainder, quotient,i,j,ib,ic,i3,j3,b,db,dv,rij)
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

do ib=mp_start,mp_end

#else
do ib=istart,iend
#endif
        ! for each bond in range:

        i  = bnd(ib)%i
        j  = bnd(ib)%j
        ic = bnd(ib)%cod
        ! calculate rij
        i3=i*3-3
        j3=j*3-3
        rij(1) = x(j3+1) - x(i3+1)
        rij(2) = x(j3+2) - x(i3+2)
        rij(3) = x(j3+3) - x(i3+3)

        ! calculate b and db, update Ebond
        b = sqrt ( rij(1)**2 + rij(2)**2 + rij(3)**2 )
        db = b - bondlib(ic)%bnd0
#ifdef _OPENMP
    mp_real_tmp = mp_real_tmp + 0.5_prec*bondlib(ic)%fk*db**2
#else  
        bond = bond + 0.5_prec*bondlib(ic)%fk*db**2
#endif

        ! calculate dv and update d
        dv = bondlib(ic)%fk*db/b
        d(j3+1) = d(j3+1) + rij(1)*dv
        d(j3+2) = d(j3+2) + rij(2)*dv
        d(j3+3) = d(j3+3) + rij(3)*dv
        d(i3+1) = d(i3+1) - rij(1)*dv
        d(i3+2) = d(i3+2) - rij(2)*dv
        d(i3+3) = d(i3+3) - rij(3)*dv
end do
#ifdef _OPENMP
!$omp atomic update
bond = bond + mp_real_tmp
#endif
!$omp end parallel

end function bond
!-----------------------------------------------------------------------
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

subroutine cgp_centers
! *** local variables
integer						::	ig,i,i3
#ifdef _OPENMP
integer :: quotient, remainder


threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (ncgp)/threads_num
remainder = MOD(ncgp, threads_num)
mp_start = thread_id * quotient + MIN(thread_id, remainder) + 1
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif

do ig = mp_start,mp_end
#else
do ig = 1, ncgp
#endif
        lrf(ig)%cgp_cent(:) = zero
        lrf(ig)%phi0 = zero
lrf(ig)%phi1(:) = zero
lrf(ig)%phi2(:) = zero
lrf(ig)%phi3(:) = zero

        do i  = cgp(ig)%first, cgp(ig)%last
                lrf(ig)%cgp_cent(:) = lrf(ig)%cgp_cent(:) + x(cgpatom(i)*3-2:cgpatom(i)*3)
        end do

lrf(ig)%cgp_cent(:) = lrf(ig)%cgp_cent(:)/real(cgp(ig)%last - cgp(ig)%first +1, kind=prec)

end do

end subroutine cgp_centers
!-----------------------------------------------------------------------
subroutine make_nbqqlist
!locals
integer						::	is

call make_qconn
nbqq_max = nbqq_count()
allocate(nbqq(nbqq_max, nstates),nbqqp(nbqq_max, nstates),stat=alloc_status)
call check_alloc('Qatom-Qatom non-bond list')
do is =1, nstates
        write (*,200) nbqq_pair(is),is
end do
write (*,*)

200 format ('No. of Rcq indep. nb pairs involving q-atoms = ',i5, &
' in state :',i3)
!qqlist now done after precomputation
end subroutine make_nbqqlist

!-----------------------------------------------------------------------

subroutine distribute_nonbonds
!locals
integer					:: npp, npw, nqp, nww, nqw
type(NODE_ASSIGNMENT_TYPE),allocatable  :: node_assignment(:)
real					:: avgload, old_avgload
integer					:: i, last_cgp, last_pair
integer					:: mpitype_pair_assignment, mpitype_node_assignment
integer                                 :: average_pairs,inode,icgp,sum,less_than_sum
integer                                 :: n_bonded, n_nonbonded, master_assign
real                                    :: percent
integer                                 :: master_sum
!!!!Tmp vars för allokering
integer,parameter			:: vars = 5
integer     :: mpi_batch, i_loop
integer     :: blockcnt(vars),type(vars)
integer(8)  :: disp(vars)
!!!


! count the number of nonbonded interactions and distribute them among the nodes

if (nodeid .eq. 0) then
! nice header
call centered_heading('Distribution of charge groups','-')

!Allocate node_assignment
allocate(node_assignment(0:numnodes-1),stat=alloc_status)
call check_alloc('node_assignment')

!Allocate arrays that hold no. pairs per chargegroup.
call allocate_nbxx_per_cgp

!Count stuff for balancing nodes and allocating nonbond arrays nbxx()
nbpp_per_cgp = 0
nbpw_per_cgp = 0
nbqp_per_cgp = 0
nbqw_per_cgp = 0
nbww_per_cgp = 0

call nbpp_count(npp, nbpp_per_cgp)    !Only for switching atoms!!!?
call nbpw_count(npw, nbpw_per_cgp)
call nbqp_count(nqp, nbqp_per_cgp)
call nbqw_count(nqw, nbqw_per_cgp) 
call nbww_count(nww, nbww_per_cgp) 

!For keeping track of actual # of nonbonded pairs
totnbpp = npp
totnbpw = npw
totnbww = nww*solv_atom**2
totnbqp = nqp
totnbqw = nqw*solv_atom*nqat

if (numnodes .eq. 1) then
! only one node: no load balancing

! make the master node handle everything
calculation_assignment%pp%start = 1
calculation_assignment%pp%end = ncgp_solute
calculation_assignment%pw%start = 1
calculation_assignment%pw%end = ncgp_solute
calculation_assignment%qp%start = 1
calculation_assignment%qp%end = ncgp_solute
calculation_assignment%qw%start = 1
calculation_assignment%qw%end = nwat
calculation_assignment%ww%start = 1
calculation_assignment%ww%end = nwat

calculation_assignment%shake%start        = 1
calculation_assignment%shake%end          = shake_molecules
calculation_assignment%natom%start        = 1
calculation_assignment%natom%end          = natom
calculation_assignment%nat_solute%start   = 1
calculation_assignment%nat_solute%end     = nat_solute
calculation_assignment%nat_solvent%start  = nat_solute + 1
calculation_assignment%nat_solvent%end    = natom
calculation_assignment%nmol%start         = 1
calculation_assignment%nmol%end           = nmol
calculation_assignment%ncgp%start         = 1
calculation_assignment%ncgp%end           = ncgp

#if defined (USE_MPI)
else ! i.e. slave nodes exists


! A simple solution to avoid parallelising the bonded
! Calculate n_bonded and n_nonbonded 
! Approximate time of computing one bonded with one nonbonded
! The number of qq-interactions are neglected
n_bonded = nbonds + nangles + ntors + nimps
n_nonbonded = totnbpp + totnbpw + totnbww + totnbqw + totnbqp

! Compare to determine how many nonbonded master should get
! A bonded is faster, so this favours an early completion for master
master_assign =  n_nonbonded/numnodes - n_bonded * numnodes

! calculate the assignments ********************

!Calculate balanced assignment for p-p pairs
icgp=0
sum=0
 !First assign master a small part
node_assignment(0)%pp%start=icgp+1
percent=REAL(totnbpp, kind=prec)/n_nonbonded
less_than_sum = master_assign*percent  ! No. of pp-type to assign master
do while((icgp .lt. ncgp_solute) .and. (sum .lt. less_than_sum))
   icgp=icgp+1
   sum=sum + nbpp_per_cgp(icgp)
end do
node_assignment(0)%pp%end=icgp
master_sum=sum
 !Now assign slaves
average_pairs=(totnbpp-sum)/(numnodes-1)
do inode=1,numnodes-2
  node_assignment(inode)%pp%start=icgp+1
  less_than_sum=average_pairs*inode+master_sum
  do while (sum .lt. less_than_sum) 
     icgp=icgp+1
     sum=sum + nbpp_per_cgp(icgp)
  end do
  node_assignment(inode)%pp%end=icgp
end do
node_assignment(numnodes-1)%pp%start=icgp+1
node_assignment(numnodes-1)%pp%end=ncgp_solute

!Calculate balanced assignment for p-w pairs
icgp=0
sum=0
node_assignment(0)%pw%start=icgp+1
percent=REAL(totnbpw, kind=prec)/n_nonbonded
less_than_sum = master_assign*percent
do while((icgp .lt. ncgp_solute) .and. (sum .lt. less_than_sum))
   icgp=icgp+1
   sum=sum + nbpw_per_cgp(icgp)
end do
node_assignment(0)%pw%end=icgp
master_sum=sum
average_pairs=(totnbpw-sum)/(numnodes-1)
do inode=1,numnodes-2
  node_assignment(inode)%pw%start=icgp+1
  less_than_sum=average_pairs*inode+master_sum
  do while (sum .lt. less_than_sum)
     icgp=icgp+1
     sum=sum + nbpw_per_cgp(icgp)
  end do
  node_assignment(inode)%pw%end=icgp
end do
node_assignment(numnodes-1)%pw%start=icgp+1
node_assignment(numnodes-1)%pw%end=ncgp_solute

!Calculate balanced assignment for q-p pairs
icgp=0
sum=0
node_assignment(0)%qp%start=icgp+1
percent=REAL(totnbqp, kind=prec)/n_nonbonded
less_than_sum = master_assign*percent
do while((icgp .lt. ncgp_solute) .and. (sum .lt. less_than_sum))
   icgp=icgp+1
   sum=sum + nbqp_per_cgp(icgp)
end do
node_assignment(0)%qp%end=icgp
master_sum=sum
average_pairs=(totnbqp-sum)/(numnodes-1)
do inode=1,numnodes-2
  node_assignment(inode)%qp%start=icgp+1
  less_than_sum=average_pairs*inode+master_sum
  do while(sum .lt. less_than_sum)
     icgp=icgp+1
     sum=sum + nbqp_per_cgp(icgp)
  end do
  node_assignment(inode)%qp%end=icgp
end do
node_assignment(numnodes-1)%qp%start=icgp+1
node_assignment(numnodes-1)%qp%end=ncgp_solute

!Calculate balanced assignment for w-w pairs
icgp=0
sum=0
node_assignment(0)%ww%start=icgp+1
percent=REAL(totnbww, kind=prec)/n_nonbonded
less_than_sum = master_assign*percent
do while((icgp .lt. nwat) .and. (sum .lt. less_than_sum))
   icgp=icgp+1
   sum=sum + nbww_per_cgp(icgp)
end do
node_assignment(0)%ww%end=icgp
master_sum=sum
average_pairs=(totnbww-sum)/(numnodes-1)
do inode=1,numnodes-2
  node_assignment(inode)%ww%start=icgp+1
  less_than_sum=average_pairs*inode+master_sum
  do while(sum .lt. less_than_sum)
     icgp=icgp+1
     sum=sum + nbww_per_cgp(icgp)
  end do
  node_assignment(inode)%ww%end=icgp
end do
node_assignment(numnodes-1)%ww%start=icgp+1
node_assignment(numnodes-1)%ww%end=nwat

!Calculate balanced assignment for q-w pairs
icgp=0
sum=0
node_assignment(0)%qw%start=icgp+1
percent=REAL(totnbqw, kind=prec)/n_nonbonded
less_than_sum = master_assign*percent
do while((icgp .lt. nwat) .and. (sum .lt. less_than_sum))
   icgp=icgp+1
   sum=sum + nbqw_per_cgp(icgp)
end do
node_assignment(0)%qw%end=icgp
master_sum=sum
average_pairs=(totnbqw-sum)/(numnodes-1)
do inode=1,numnodes-2
  node_assignment(inode)%qw%start=icgp+1
  less_than_sum=average_pairs*inode+master_sum
  do while(sum .lt. less_than_sum)
     icgp=icgp+1
     sum=sum + nbqw_per_cgp(icgp)
  end do
  node_assignment(inode)%qw%end=icgp
end do
node_assignment(numnodes-1)%qw%start=icgp+1
node_assignment(numnodes-1)%qw%end=nwat

!now search the shake constrained molecules and give a balanced assignment for
!each node, based on the number of shake constraints in each searched molecule
icgp=0
sum=0
average_pairs=REAL((shake_constraints/numnodes),kind=prec)
do inode=0,numnodes-2
        node_assignment(inode)%shake%start=icgp+1
        less_than_sum = (inode+1)*average_pairs
        do while ((sum .lt. less_than_sum).and.(icgp.lt.shake_molecules))
             icgp=icgp+1
             sum=sum + shake_mol(icgp)%nconstraints
        end do
        node_assignment(inode)%shake%end=icgp
end do
node_assignment(numnodes-1)%shake%start=icgp+1
node_assignment(numnodes-1)%shake%end=shake_molecules

!assignment of atom numbers
sum=0
average_pairs=REAL((natom/numnodes),kind=prec)
do inode = 0,numnodes-2
        node_assignment(inode)%natom%start = sum + 1
        less_than_sum = (inode+1)*average_pairs
        do while (sum.lt.less_than_sum)
                sum = sum + average_pairs
        end do
        node_assignment(inode)%natom%end   = sum
end do
node_assignment(numnodes-1)%natom%start = node_assignment(numnodes-2)%natom%end + 1
node_assignment(numnodes-1)%natom%end   = natom

!assignment of solute atoms
sum=0
average_pairs=REAL((nat_solute/numnodes),kind=prec)
do inode = 0,numnodes-2
        node_assignment(inode)%nat_solute%start = sum + 1
        less_than_sum = (inode+1)*average_pairs
        do while (sum.lt.less_than_sum)
                sum = sum + average_pairs
        end do
        node_assignment(inode)%nat_solute%end   = sum
end do
node_assignment(numnodes-1)%nat_solute%start = & 
node_assignment(numnodes-2)%nat_solute%end + 1
node_assignment(numnodes-1)%nat_solute%end   = nat_solute

!assignment of solvent atoms
sum=nat_solute
average_pairs=REAL(((natom-nat_solute+1)/numnodes),kind=prec)
do inode = 0,numnodes-2
        node_assignment(inode)%nat_solvent%start = sum + 1
        less_than_sum = (inode+1)*average_pairs + nat_solute
        do while (sum.lt.less_than_sum)
                sum = sum + average_pairs
        end do
        node_assignment(inode)%nat_solvent%end   = sum
end do
node_assignment(numnodes-1)%nat_solvent%start = & 
node_assignment(numnodes-2)%nat_solvent%end + 1
node_assignment(numnodes-1)%nat_solvent%end   = natom

sum=0
average_pairs=REAL((nmol/numnodes),kind=prec)
do inode = 0, numnodes-2
        node_assignment(inode)%nmol%start = sum + 1
        less_than_sum = (inode+1)*average_pairs
        do while (sum.lt.less_than_sum)
                sum = sum + average_pairs
        end do
        node_assignment(inode)%nmol%end   = sum
end do
node_assignment(numnodes-1)%nmol%start = &   
node_assignment(numnodes-2)%nmol%end + 1
node_assignment(numnodes-1)%nmol%end   = nmol

sum=0
average_pairs=REAL((ncgp/numnodes),kind=prec)
do inode = 0, numnodes-2
        node_assignment(inode)%ncgp%start = sum + 1
        less_than_sum = (inode+1)*average_pairs
        do while (sum.lt.less_than_sum)
                sum = sum + average_pairs
        end do
        node_assignment(inode)%ncgp%end   = sum
end do
node_assignment(numnodes-1)%ncgp%start = &    
node_assignment(numnodes-2)%ncgp%end + 1
node_assignment(numnodes-1)%ncgp%end   = ncgp


#endif
end if    !if (numnodes .eq. 1)

! deallocate bookkeeping arrays
!deallocate(nppcgp, npwcgp, nqpcgp, nwwmol)

end if   !if (nodeid .eq. 0)

! distribute assignments to the nodes
#if defined (USE_MPI)
if (numnodes .gt. 1) then
    if (nodeid .ne. 0) then
	! Dummy allocation to avoid runtime errors when using pointer checking
	allocate(node_assignment(1),stat=alloc_status)
    endif 
! register data types
call MPI_Type_contiguous(3, MPI_INTEGER, mpitype_pair_assignment, ierr)
if (ierr .ne. 0) call die('failure while creating custom MPI data type')
call MPI_Type_commit(mpitype_pair_assignment, ierr)
if (ierr .ne. 0) call die('failure while creating custom MPI data type')

call MPI_Type_contiguous(11, mpitype_pair_assignment, mpitype_node_assignment, ierr)
if (ierr .ne. 0) call die('failure while creating custom MPI data type')
call MPI_Type_commit(mpitype_node_assignment, ierr)
if (ierr .ne. 0) call die('failure while creating custom MPI data type')

! distribute
call MPI_Scatter(node_assignment, 1, mpitype_node_assignment, &
    calculation_assignment, 1, mpitype_node_assignment, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('failure while sending node assignments')

! free data type
call MPI_Type_free(mpitype_node_assignment, ierr)
call MPI_Type_free(mpitype_pair_assignment, ierr)
    if (nodeid .ne. 0) then
	deallocate(node_assignment)
    endif 
end if
#endif

if (nodeid .eq. 0) then
 ! print a status report
 write(*,98) 'solute-solute', 'solute-water', 'water-water', 'Q-solute', 'Q-water', '#Shake', '#Atoms', '#Molecules'
 write(*,99) 'total',ncgp_solute,ncgp_solute,nwat,ncgp_solute,nwat,shake_molecules,natom,nmol
if (numnodes .gt. 1) then
 do i=0,numnodes-1
    write(*,100) i, 'assigned cgps', &
         node_assignment(i)%pp%end-node_assignment(i)%pp%start+1, &
         node_assignment(i)%pw%end-node_assignment(i)%pw%start+1, &
         node_assignment(i)%ww%end-node_assignment(i)%ww%start+1, &
         node_assignment(i)%qp%end-node_assignment(i)%qp%start+1, &
         node_assignment(i)%qw%end-node_assignment(i)%qw%start+1, &
         node_assignment(i)%shake%end-node_assignment(i)%shake%start+1, &
         node_assignment(i)%natom%end-node_assignment(i)%natom%start+1, &
         node_assignment(i)%nmol%end-node_assignment(i)%nmol%start+1
 end do
end if
end if

#if defined (USE_MPI)
blockcnt(:) = 1
type(:) = MPI_INTEGER
call MPI_Bcast(totnbpp, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
call MPI_Bcast(totnbpw, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
call MPI_Bcast(totnbqp, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
call MPI_Bcast(totnbqw, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
call MPI_Bcast(totnbww, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
#endif

! allocate
calculation_assignment%pp%max = totnbpp/numnodes + 0.20*totnbpp
allocate(nbpp(calculation_assignment%pp%max), stat=alloc_status)
call check_alloc('solute-solute non-bond list')

calculation_assignment%pw%max = (totnbpw+3)/numnodes + 0.20*totnbpw
allocate(nbpw(calculation_assignment%pw%max), stat=alloc_status)
call check_alloc('solute-solvent non-bond list')

calculation_assignment%ww%max = (totnbww+nwat)/numnodes + 0.20*totnbww
allocate(nbww(calculation_assignment%ww%max), stat=alloc_status)
call check_alloc('solvent-solvent non-bond list')

calculation_assignment%qp%max = totnbqp/numnodes + 0.20*totnbqp
allocate(nbqp(calculation_assignment%qp%max,nstates), stat=alloc_status)
call check_alloc('Qatom - solute non-bond list')

calculation_assignment%qw%max = totnbqw/numnodes + 0.20*totnbqw
allocate(nbqw(calculation_assignment%qw%max,nstates), stat=alloc_status)
call check_alloc('Qatom - water non-bond list')



  if (use_PBC) then
	!allocate array to keep track of chargegroups
	!approximate with one half of the number of atompairs
	allocate(nbpp_cgp(calculation_assignment%pp%max / 2), stat=alloc_status)
	call check_alloc('solute-solute non-bond charge group pair list')
	allocate(nbpw_cgp(calculation_assignment%pw%max / 2), stat=alloc_status)
	call check_alloc('solute-solvent non-bond charge group pair list')
	allocate(nbqp_cgp(calculation_assignment%qp%max / 2), stat=alloc_status)
	call check_alloc('qatom-solute non-bond charge group pair list')
!we now allocate the storage arrays for MC_volume update here, so we don't have to realloc every time 
!the function is called
        allocate(old_nbpp_cgp(calculation_assignment%pp%max / 2), stat=alloc_status)
        call check_alloc('solute-solute non-bond charge group pair list')
        allocate(old_nbpw_cgp(calculation_assignment%pw%max / 2), stat=alloc_status)
        call check_alloc('solute-solvent non-bond charge group pair list')
        allocate(old_nbqp_cgp(calculation_assignment%qp%max / 2), stat=alloc_status)
        call check_alloc('qatom-solute non-bond charge group pair list')

allocate(old_nbpp(calculation_assignment%pp%max), stat=alloc_status)
call check_alloc('solute-solute non-bond list')

allocate(old_nbpw(calculation_assignment%pw%max), stat=alloc_status)
call check_alloc('solute-solvent non-bond list')

allocate(old_nbww(calculation_assignment%ww%max), stat=alloc_status)
call check_alloc('solvent-solvent non-bond list')

allocate(old_nbqp(calculation_assignment%qp%max,nstates), stat=alloc_status)
call check_alloc('Qatom - solute non-bond list')
  end if	

!Kanske deallokera nbxx_per_cgp TODO

98 format('node value ',8a13)
99 format(a10,1x,8(1x,i12))
!99 format(a4,2x,a,t18,i13,3x,i13,3x,i13,3x,i13)
100 format(i4,1x,a5,1x,8(1x,i12))
!100 format(i4,2x,a,t18,i13,3x,i13,3x,i13,3x,i13)

if (nodeid .eq. 0)  call centered_heading('End of distribution', '-')

end subroutine distribute_nonbonds

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

!Restrain all excluded atoms plus heavy solute atoms in the inner shell.
subroutine fix_shell
! local variables
integer						::	i,i3
real(kind=prec)						::	fk,r2,erst
real(kind=prec)						::	dr(3)

! global variables used:
!  E, nat_pro, excl, shell, heavy, fk_fix, fk_pshell, x, xtop, d

do i = 1, nat_pro
if (excl(i) .or. shell(i)) then
  ! decide which fk to use
  i3 = i*3 - 3
if ( excl(i) ) then 
    fk = fk_fix
    if ( freeze ) then
       v(i3+1)=0
       v(i3+2)=0
       v(i3+3)=0
       x(i3+1)=xtop(i3+1)
       x(i3+2)=xtop(i3+2)
       x(i3+3)=xtop(i3+3)
       fk = 0
    endif
  else
fk = fk_pshell
  end if
!i3 = i*3-3

  ! calculate drift from topology
dr(1)   = x(i3+1) - xtop(i3+1)
dr(2)   = x(i3+2) - xtop(i3+2)
dr(3)   = x(i3+3) - xtop(i3+3)
r2      = dr(1)**2 + dr(2)**2 + dr(3)**2
erst    = 0.5_prec*fk*r2

  ! update restraint energies
if ( excl(i) ) E%restraint%fix   = E%restraint%fix + erst
if ( shell(i) ) E%restraint%shell = E%restraint%shell + erst 

  ! update forces
d(i3+1) = d(i3+1) + fk*dr(1)
d(i3+2) = d(i3+2) + fk*dr(2)
d(i3+3) = d(i3+3) + fk*dr(3)
end if
end do
end subroutine fix_shell

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

#ifdef USE_GRID
#ifdef USE_MPI
integer function grid_add(grid1,grid2,length,MPI_Datatype)
! arguments
        integer,INTENT (IN)      :: grid1(length)
        integer                  :: grid2(length)
        integer                         :: length,MPI_Datatype
! locals
        integer                 :: j,k

! we have to do something nasty to be able to combine the grid information form
! all the nodes, so here it is
! we search first for a zero element in the final grid, so we not overwrite any
! info that is already there
! as soons as we find the first zero element, we search the slave grid in the
! same row until we find a zero element, and copy all the groups from before
! into the master array
! quick and very dirty way
        j = 1
        do while (grid2(j) .ne. 0 )
                j = j + 1
        end do
        k = 1
        do while (grid1(k) .ne. 0 )
                grid2(j+k) = grid1(k)
        end do

        grid_add = -1

end function grid_add
#endif
#endif

real(kind=prec) function improper(istart, iend)
!arguments
integer						::	istart, iend

! evaluate harmonic impropers
! local variables
integer						::	ip
real(kind=prec)						::	scp,phi,dv,arg,f1
real(kind=prec)						::	bjinv, bkinv, bj2inv, bk2inv
real(kind=prec)						::	rji(3),rjk(3),rkl(3),rnj(3),rnk(3)
real(kind=prec)						::	rki(3),rlj(3),dp(12),di(3),dl(3)
type(TOR_TYPE), pointer		::	t
type(IMPLIB_TYPE), pointer	::	lib
#ifdef _OPENMP
integer :: quotient, remainder
#endif

! global variables used:
!  imp, implib, x, pi, d

improper = zero

!$omp parallel default(none) shared(threads_num, istart, iend, imp, implib, x, pi, d, improper) private(quotient, remainder,ip,scp,phi,dv,arg,f1,bjinv,bkinv,bj2inv,bk2inv,rji,rjk,rkl,rnj,rnk,rki,rlj,dp,di,dl,t,lib)
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
        t => imp(ip)
        lib => implib(t%cod)
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

bj2inv  = one/( rnj(1)**2 + rnj(2)**2 + rnj(3)**2)
bk2inv  = one/( rnk(1)**2 + rnk(2)**2 + rnk(3)**2)
        bjinv = sqrt(bj2inv)
        bkinv = sqrt(bk2inv)

scp = (rnj(1)*rnk(1)+rnj(2)*rnk(2)+rnj(3)*rnk(3))*(bjinv*bkinv)
if ( scp .gt.  one ) then
                scp =  one
else if ( scp .lt. -one ) then 
                scp = -one
        end if
phi = acos ( scp )
if(rjk(1)*(rnj(2)*rnk(3)-rnj(3)*rnk(2)) &
     +rjk(2)*(rnj(3)*rnk(1)-rnj(1)*rnk(3)) &
     +rjk(3)*(rnj(1)*rnk(2)-rnj(2)*rnk(1)) &
                 .lt. 0) then 
      phi = -phi
        end if

! ---       energy

arg = phi - lib%imp0
arg = arg - 2.0_prec*pi*nint(arg/(2.0_prec*pi))
dv  = lib%fk*arg
#ifdef _OPENMP
    mp_real_tmp = mp_real_tmp + 0.5_prec*dv*arg
#else 
    improper = improper + 0.5_prec*dv*arg
#endif

! ---       forces

f1 = sin ( phi ) 
if ( abs(f1) .lt. 1.e-12_prec ) f1 = 1.e-12_prec
f1 =  -one / f1
di(1) = f1 * ( rnk(1)*bjinv*bkinv - scp*rnj(1)*bj2inv )
di(2) = f1 * ( rnk(2)*bjinv*bkinv - scp*rnj(2)*bj2inv )
di(3) = f1 * ( rnk(3)*bjinv*bkinv - scp*rnj(3)*bj2inv )
dl(1) = f1 * ( rnj(1)*bjinv*bkinv - scp*rnk(1)*bk2inv )
dl(2) = f1 * ( rnj(2)*bjinv*bkinv - scp*rnk(2)*bk2inv )
dl(3) = f1 * ( rnj(3)*bjinv*bkinv - scp*rnk(3)*bk2inv )

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
improper = improper + mp_real_tmp
#endif
!$omp end parallel
end function improper

!-----------------------------------------------------------------------

real(kind=prec) function improper2(istart, iend)
!evaluate periodic impropers
!arguments
integer						::	istart, iend
! local variables
integer						::	ip
real(kind=prec)						::	scp,phi,dv,arg,f1
real(kind=prec)						::	bjinv, bkinv, bj2inv, bk2inv
real(kind=prec)						::	rji(3),rjk(3),rkl(3),rnj(3),rnk(3)
real(kind=prec)						::	rki(3),rlj(3),dp(12),di(3),dl(3)
type(TOR_TYPE), pointer		::	t
type(IMPLIB_TYPE), pointer	::	lib
#ifdef _OPENMP
integer :: quotient, remainder
#endif

! global variables used:
! imp, implib, x, pi, d

improper2 = zero

!$omp parallel default(none) shared(threads_num, istart, iend, imp, implib, x, pi, d, improper2) private(remainder, quotient,ip,scp,phi,dv,arg,f1,bjinv,bkinv,bj2inv,bk2inv,rji,rjk,rkl,rnj,rnk,rki,rlj,dp,di,dl,t,lib)
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
        t => imp(ip)
        lib => implib(t%cod)
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


bj2inv  = one/( rnj(1)**2 + rnj(2)**2 + rnj(3)**2)
bk2inv  = one/( rnk(1)**2 + rnk(2)**2 + rnk(3)**2)
        bjinv = sqrt(bj2inv)
        bkinv = sqrt(bk2inv)

scp = (rnj(1)*rnk(1)+rnj(2)*rnk(2)+rnj(3)*rnk(3))*(bjinv*bkinv)
if ( scp .gt.  one ) then
                scp =  one
else if ( scp .lt. -one ) then 
                scp = -one
        end if
phi = acos ( scp )
if(rjk(1)*(rnj(2)*rnk(3)-rnj(3)*rnk(2)) &
     +rjk(2)*(rnj(3)*rnk(1)-rnj(1)*rnk(3)) &
     +rjk(3)*(rnj(1)*rnk(2)-rnj(2)*rnk(1)) &
                 .lt. 0) then 
      phi = -phi
        end if

! ---       energy

arg = 2*phi - lib%imp0
#ifdef _OPENMP
    mp_real_tmp = mp_real_tmp + lib%fk * (1 + cos(arg))
#else
    improper2 = improper2 + lib%fk * (1 + cos(arg))
#endif
dv  = -2*lib%fk * sin(arg)

! ---       forces

f1 = sin ( phi ) 
if ( abs(f1) .lt. 1.e-12_prec ) f1 = 1.e-12_prec
f1 =  -one / f1
di(1) = f1 * ( rnk(1)*bjinv*bkinv - scp*rnj(1)*bj2inv )
di(2) = f1 * ( rnk(2)*bjinv*bkinv - scp*rnj(2)*bj2inv )
di(3) = f1 * ( rnk(3)*bjinv*bkinv - scp*rnj(3)*bj2inv )
dl(1) = f1 * ( rnj(1)*bjinv*bkinv - scp*rnk(1)*bk2inv )
dl(2) = f1 * ( rnj(2)*bjinv*bkinv - scp*rnk(2)*bk2inv )
dl(3) = f1 * ( rnj(3)*bjinv*bkinv - scp*rnk(3)*bk2inv )

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
improper2 = improper2 + mp_real_tmp
#endif
!$omp end parallel
end function improper2

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

subroutine initial_shaking
!
! initial shaking
!
integer						:: niter


xx(:)=x(:)
niter=shake(xx, x)	
write(*,100) 'x', niter
100	format('Initial ',a,'-shaking required',i4,&
' interations per molecule on average.')

xx(:)=x(:)-dt*v(:)
niter=shake(x, xx)	
write(*,100) 'v', niter

v(:)=(x(:)-xx(:))/dt


end subroutine initial_shaking

!-----------------------------------------------------------------------
logical function initialize()                  
! local variables
character					:: text*200
integer						:: i,j,length,ii
real(kind=prec)						:: stepsize
real(kind=prec)						:: lamda_tmp(max_states)
integer						:: fu, fstat
real(kind=prec)						::	rjunk
integer						::	ijunk

! local parameters
integer						:: num_args
character(200)				:: infilename
logical						::	yes
logical						::	need_restart
character(len=200)				::	instring
logical						::	inlog
integer						::	mask_rows, number
integer,parameter				:: maxmaskrows=30
character(len=200)				:: gc_mask_tmp(maxmaskrows),str,str2
logical                                         :: shake_all,shake_all_solvent,shake_all_solute
logical                                         :: shake_all_hydrogens,shake_all_heavy

! this subroutine will init:
!  nsteps, stepsize, dt
!  Temp0, tau_T, iseed, Tmaxw
!  use_LRF, NBcycle, Rcpp, Rcww, Rcpw, Rcq
!  shake_solute, shake_solvent, shake_hydrogens, shake_heavy
! fk_pshell
!  fk_wsphere=-1, wpol_restr, wpol_born
!  fkwpol=-1, Dwmz=-1 (values  ized to -1 will
!    be set in water_sphere, once target radius is known)
!  top_file
!  restart, [restart_file]
!  xfin_file
!  itrj_cycle, iene_cycle, iout_cycle, itemp_cycle, [trj_file], [ene_file]
!  fep_file
!  nstates, EQ (allocating memory for EQ)
!  implicit_rstr_from_file, [exrstr_file]
!  nrstr_seq, [rstseq] (allocating memory for rstseq)
!  nrstr_pos, [rstpos] (allocating memory for rstpos)
!  nrstr_dist, [rstdis] (allocating memory for rstdis)
!  nrstr_wall, [rstwal] (allocating memory for rstwal)

! external definition of iargc disabled for gfortran
!integer(4) iargc
!external iargc

! read name of input file from the command line
num_args = command_argument_count()
if (num_args .lt. 1) call die('no input file specified on the command line')
#if defined(CRAY)
call pxfgetarg(num_args, infilename, 200, i)
#elif defined(MPICH)
call getarg(1, infilename)
#else
call getarg(num_args, infilename)
#endif
text = 'Reading input from '//infilename
call centered_heading(trim(text), '-')

initialize = .true. 

if(.not. prm_open_section('PBC', infilename)) then
        box = .false.
        write(*,'(a)') 'Boundary: sphere'
else
        box = .true.
        write(*,'(a)') 'Boundary: periodic box'
        if( .not. prm_get_logical_by_key('rigid_box_centre', rigid_box_centre, .false. ) ) then
                write(*,'(a)') '>>> Error: rigid_box_centre must be on or off'
                initialize = .false.
        end if
        write(*,'(a,a3)') 'Rigid box centre ', onoff(rigid_box_centre)
        if( .not. prm_get_logical_by_key('constant_pressure', constant_pressure, .false.) ) then
                write(*,'(a)') '>>> Error: constant_pressure must be on or off'
                initialize = .false.
        end if

        if( constant_pressure ) then
                write(*,'(a)') 'NPT-ensemble'
                volume_try = 0
                volume_acc = 0
                if( .not. prm_get_real8_by_key('max_volume_displ', max_vol_displ) ) then
                        initialize = .false.
                        write(*,'(a)') '>>> ERROR: maximum volume displacement not specified (section PBC)'
                else
                        write(*,5) max_vol_displ
                end if
5	format ('Maximum volume displacemet = ', f10.3)

                if( .not. prm_get_integer_by_key('pressure_seed', pressure_seed)) then
                        pressure_seed = 3781
                end if

				write(*, '(a, i4 )' ) 'Pressure seed: ', pressure_seed

                if( .not. prm_get_real8_by_key('pressure', pressure) ) then
                        pressure = one
                end if
                write(*,9) pressure
9	format ('Pressure = ',f10.3,'  bar')
                !convert pressure to strange internal unit
                pressure = pressure * 1.43836e-5_prec
                yes = prm_get_logical_by_key('atom_based_scaling', atom_based_scaling, .false.)
                if (atom_based_scaling) then
                        write (*,'(a)') 'Coordinate scaling on volume changes:  Atom based'
                else
                        write (*,'(a)') 'Coordinate scaling on volume changes:  Molecule based'
                end if

        else
                write(*,'(a)') 'NVT-ensemble'
                if( prm_get_line_by_key('control_box', instring) ) then
                        read(instring, *) new_boxl(:)
                        control_box = .true.
                        write(*,'(a, 3f10.3)')'Boxsize will be changed to: ', new_boxl
                else
                        control_box = .false.
                end if


        end if !section constant_pressure

	yes = prm_get_logical_by_key('put_solvent_back_in_box', put_solvent_back_in_box)

	yes = prm_get_logical_by_key('put_solute_back_in_box', put_solute_back_in_box)


	if(put_solute_back_in_box .and. put_solvent_back_in_box) then
           write(*,'(a)') 'Solute and solvent molecules will be put back in box.'
	else
		if (put_solute_back_in_box) then
           	write(*,'(a)') 'Only solute molecules will be put back in box.'
		else
			if (put_solvent_back_in_box) then
				write(*,'(a)') 'Only solvent molecules will be put back in box.'
			else
				write(*,'(a)') 'No molecules will be put back in box.'				
			end if
		end if
	end if



end if !section PBC


if(.not. prm_open_section('md')) then
        call prm_close
        ! open input file
        fu = freefile()
        open(unit=fu, file=infilename, action='read', form='formatted', status='old', iostat=fstat)
        if (fstat .ne. 0) call die('error opening input file '//infilename)
                initialize = old_initialize(fu)
                close(fu)
        return
end if

need_restart = .false. !flag for restart file required
if(.not. prm_get_integer_by_key('steps', nsteps)) then
        write(*,*) '>>> ERROR: steps not specified (section MD)'
        initialize = .false.
end if
if(.not. prm_get_real8_by_key('stepsize', stepsize)) then
        write(*,*) '>>> ERROR: stepsize not specified (section MD)'
        initialize = .false.
end if
write (*,10) nsteps, stepsize
10	format ('Number of MD steps =',i10,'  Stepsize (fs)    =',f10.3)

! convert to internal time units once and for all.
dt=0.020462_prec*stepsize
dt2=0.5_prec*dt
! --- Temperature, Thermostat etc.
if(.not. prm_get_real8_by_key('temperature', Temp0)) then
        write(*,*) '>>> ERROR: temperature not specified (section MD)'
        initialize = .false.
end if

if(.not. prm_get_real8_by_key('bath_coupling', tau_T)) then
        write(*,*) 'Temperature bath relaxation time tau_T set to default'
        tau_T = tau_T_default
end if

write (*,15) Temp0,tau_T
tau_T=0.020462_prec*tau_T
if(Temp0 <= 0) then
        write(*,'(a)') &
                '>>> Error: No dynamics at zero temperature!'
        initialize = .false.
end if
if(tau_T < dt) then
        write(*,'(a)') '>>> Error: tau_t must be >= stepsize.'
        initialize = .false.
end if

if(.not. prm_get_string_by_key('thermostat', name_thermostat)) then
        write(*,*) 'No thermostat chosen. Berendsen thermostat will be used.'
        thermostat = BERENDSEN

else if( name_thermostat == 'langevin' ) then
	write(*,*) 'Thermostat chosen: ', name_thermostat
	thermostat = LANGEVIN
	if(.not. prm_get_real8_by_key('langevin_friction', friction)) then
		friction = 1/tau_T !***according to GROMACS manual, this is their default value. Need to check - A. Barrozo
!		gkT = 2*friction*Boltz*Temp0/dt !constant to be used to generate the random forces of the thermostat SET LATER	
		write(*,*) 'Langevin thermostat friction constant set to default: 1/tau_T'
	end if
        if(.not. prm_get_logical_by_key('langevin_random', langevin_gauss_rand )) then
        !we use uniform random noise by default for random numbers in the
        !langevin thermostat, Paul October 2014
                write(*,*) 'Using uniform random numbers for Langevin thermostat'
        else if (langevin_gauss_rand) then
                write(*,*) 'Using Gaussian random numbers for Langevin thermostat'
        else if (.not.langevin_gauss_rand) then
                 write(*,*) 'Using uniform random numbers for Langevin thermostat'
        end if
else if( name_thermostat == 'nose-hoover' ) then
        write(*,*) 'Thermostat chosen: Nose-Hoover'
		thermostat = NOSEHOOVER
        kbT = Boltz*Temp0

else if( name_thermostat /= 'berendsen' .and. name_thermostat /= 'langevin' .and. name_thermostat /= 'nose_hoover' ) then
	write(*,*) '>>> ERROR: this thermostat does not exist in Q.',name_thermostat
	initialize = .false.
else
		thermostat = BERENDSEN
		name_thermostat = 'berendsen'
        write(*,*) 'Thermostat chosen: ', name_thermostat
end if

if(.not. prm_get_integer_by_key('nhchains', numchain) .and. thermostat == NOSEHOOVER ) then
	numchain = 10
        write(*,*) 'Nose-Hoover thermostat chain number set to default: 10'
end if

if(.not. prm_get_real8_by_key('nose-hoover_mass', nhq) .and. thermostat == NOSEHOOVER ) then
	nhq = kbT/(tau_T*tau_T)
        write(*,*) 'Nose-Hoover thermostat mass set to default: kbT/tau_T^2' !Based on Martyna, Klein and Tuckerman J.Chem. Phys. 92
end if


if(.not. prm_get_string_by_key('integrator', name_integrator)) then
        write(*,*) 'Leap frog integrator used by default.'
        integrator = LEAPFROG

else if( name_integrator == 'velocity-verlet' ) then
	write(*,*) 'Integrator: ', name_integrator
	integrator = VELVERLET
else if( name_integrator /= 'leap-frog' .and. name_integrator /= 'velocity-verlet' ) then
	write(*,*) '>>> ERROR: no such integrator exists in Q.',name_integrator
	initialize = .false.
else
		name_integrator = 'leap-frog'
		integrator = LEAPFROG
        write(*,*) 'Integrator: ', name_integrator
end if
yes = prm_get_logical_by_key('separate_scaling', separate_scaling, .true.)
if(separate_scaling) then
	   write(*,'(a)') 'Solute and solvent atoms coupled separately to heat bath.'
else 
	   write(*,'(a)') 'Solute and solvent atoms coupled together to heat bath.'
end if

15	format ('Target temperature =',f10.2,'  T-relax time     =',f10.2)

yes = prm_get_integer_by_key('random_seed', iseed, 1) 
if(.not. prm_get_real8_by_key('initial_temperature', Tmaxw)) then
        iseed = 0 !set iseed = 0 if no initial temp
        need_restart = .true.
end if

if (iseed > 0) write (*,16) Tmaxw, iseed
16	format ('Initial velocities will be generated from Maxwell distribution:',&
                /,'Maxwell temperature=',f10.2,' Random number seed=',i10)

! --- shake, LRF
if(.not. prm_get_logical_by_key('shake_solvent', shake_solvent, .true.)) then
        write(*,'(a)') '>>> Error: shake_solvent must be on or off.'
        initialize = .false.
end if
!write(*,17) 'all solvent bonds', onoff(shake_solvent)
17	format('Shake constraints on ',a,' to ',a,' : ',a3)

if(.not. prm_get_logical_by_key('shake_solute', shake_solute, .false.)) then
        write(*,'(a)') '>>> Error: shake_solute must be on or off.'
        initialize = .false.
end if 
!write(*,17) 'all solute bonds', onoff(shake_solute)

if(.not. prm_get_logical_by_key('shake_hydrogens', shake_hydrogens, .true.)) then
        write(*,'(a)') '>>> Error: shake_hydrogens must be on or off.'
        initialize = .false.
end if 

if(.not. prm_get_logical_by_key('shake_heavy', shake_heavy, .false.) ) then
        write(*,'(a)') '>>> Error: shake_heavy must be on or off.'
        initialize = .false.
end if


if(prm_get_logical_by_key('shake_all',shake_all)) then
        if(shake_all) then
               shake_solute    = .true.
               shake_solvent   = .true.
               shake_hydrogens = .true.
               shake_heavy     = .true.
        else
               shake_solute    = .false.
               shake_solvent   = .false.
               shake_hydrogens = .false.
               shake_heavy     = .false.
        end if
end if

if(prm_get_logical_by_key('shake_all_solvent',shake_all_solvent)) then
        if(shake_all_solvent) then
               shake_solvent   = .true.
               shake_hydrogens = .true.
               shake_heavy     = .true.
        else
               shake_solvent   = .false.
        end if
end if

if(prm_get_logical_by_key('shake_all_solute',shake_all_solute)) then
        if(shake_all_solute) then
               shake_solute    = .true.
               shake_hydrogens = .true.
               shake_heavy     = .true.
        else
               shake_solute    = .false.
        end if
end if

if(prm_get_logical_by_key('shake_all_hydrogens',shake_all_hydrogens)) then
        if(shake_all_hydrogens) then
                shake_solute    = .true.
                shake_solvent   = .true.
                shake_hydrogens = .true.
        else
                shake_hydrogens = .false.
        end if
end if

if(prm_get_logical_by_key('shake_all_heavy',shake_all_heavy)) then
        if(shake_all_heavy) then
                 shake_solute  = .true.
                 shake_solvent = .true.
                 shake_heavy   = .true.
        else
                 shake_heavy   = .false.
        end if
end if

if(shake_solvent) then
        if(shake_heavy) write(*,17) 'solvent bonds','heavy atoms',onoff(shake_heavy)
        if(shake_hydrogens) write(*,17) 'solvent bonds','hydrogens',onoff(shake_hydrogens)
else
        write(*,17) 'solvent bonds','any atoms',onoff(shake_solvent)
endif
if(shake_solute) then
        if(shake_heavy) write(*,17) 'solute bonds','heavy atoms',onoff(shake_heavy)
        if(shake_hydrogens) write(*,17) 'solute bonds','hydrogens',onoff(shake_hydrogens)
else
        write(*,17) 'solute bonds','any atoms',onoff(shake_solute)
endif
!write(*,17) 'all bonds to hydrogen', onoff(shake_hydrogens)


       yes = prm_get_logical_by_key('lrf', use_LRF, .true.)
       if(use_LRF) then
               write(*,20) 'LRF Taylor expansion outside cut-off'
       else 
               write(*,20) 'standard cut-off'
       end if

20	format ('Nonbonded method   = ',a)

yes = prm_get_logical_by_key('force_rms', force_rms, .false.)
if(force_rms) then
        write(*,22) 
end if
22	format ('R.M.S. force will be calculated.')


! --- Rcpp, Rcww, Rcpw, Rcq, RcLRF
! change for PBC, q-atoms interact with everything for now, placeholder until we have 
! fixed it to calculate the interaction between the mirror images
if(.not. prm_open_section('cut-offs')) then
        write(*,'(a)') 'No cut-offs section, default cut-offs used'
        Rcpp = rcpp_default
        Rcww = rcww_default
        Rcpw = rcpw_default
	if (box) then
        Rcq = rcq_default_pbc
	else
	Rcq = rcq_default_sph
	end if
        RcLRF = rcLRF_default
else
        if(.not. prm_get_real8_by_key('solute_solute', Rcpp, rcpp_default)) then
                write(*,'(a)') 'solute-solute cut-off set to default'
        end if
        if(.not. prm_get_real8_by_key('solvent_solvent', Rcww, rcww_default)) then
                write(*,'(a)') 'solvent-solvent cut-off set to default'
        end if
        if(.not. prm_get_real8_by_key('solute_solvent', Rcpw, rcpw_default)) then
                write(*,'(a)') 'solute-solvent cut-off set to default'
        end if
	if (box) then
        if(.not. prm_get_real8_by_key('q_atom', Rcq, rcq_default_pbc)) then
                write(*,'(a)') 'q-atom cut-off set to default for PBC'
	end if
	else
	if(.not. prm_get_real8_by_key('q_atom', Rcq, rcq_default_sph)) then
                write(*,'(a)') 'q-atom cut-off set to default for sphere'
	end if
        end if
        if(use_LRF) then
                if(.not. prm_get_real8_by_key('lrf', rcLRF, rcLRF_default)) then
                        write(*,'(a)') 'LRF cut-off set to default'
                end if
                if(RcLRF < rcpp .or. RcLRF < rcpw .or. RcLRF < rcww) then
                        if (box) then
                              rcLRF = rcLRF_default_pbc
                              write(*,'(a)') 'LRF cut-off set to default for PBC'
                        else
                              write(*,'(a)') &
                                '>>> ERROR; LRF cut-off must not be smaller than solute or solvent cut-offs!'
                              initialize = .false.
                        end if
                end if
        end if
end if

write (*,25) Rcpp,Rcww,Rcpw,Rcq
if(use_LRF) write(*,26) RcLRF
25	format ('Cut-off radii for non-bonded interactions:',/, &
                'Solute-solute:    ',f6.2,/,&
                'Solvent-solvent:  ',f6.2,/,&
                'Solute-solvent:   ',f6.2,/,&
                'Q-atom-non-Q-atom:',f6.2)
26	format ('LRF:              ',f6.2)

30	format ('>>> WARNING: Ignoring obsolete keyword ',a,'.')
! --- simulation sphere

if( .not. box ) then
        if(.not. prm_open_section('sphere')) then
			fk_pshell = fk_pshell_default
			print*,'Radius of inner restrained shell set to 85% of exclusion shell radius.'
			rexcl_i = shell_default
			write(*,50) rexcl_i
        else
                if(prm_get_line_by_key('centre', instring)) then
                        write(*,30) 'centre'
                end if
                ! --- rexcl_o, rexcl_i, fk_pshell
                if(prm_get_real8_by_key('radius', rjunk)) then
                        write(*,30) 'radius'
                end if
                if(prm_get_real8_by_key('shell_radius', rexcl_i)) then  !inner radius of restrained shell
                    write(*,50) rexcl_i
                    if(rexcl_i < zero) then
                      call die('inner radius of restrained shell must be >= 0')
                    end if
                else
                    print*,'Radius of inner restrained shell set to 85% of exclusion shell radius.'
					rexcl_i = shell_default
                    write(*,50) rexcl_i
                end if
50 format('Radius of inner restrained shell       =    ',f8.3) 
                if(.not. prm_get_real8_by_key('shell_force', fk_pshell)) then
                        write(*,'(a)') 'Shell force constant set to default'
                        fk_pshell = fk_pshell_default
                end if
                if(fk_pshell > 0) then
                        write(*,47) fk_pshell
                end if
47		format('Shell restraint force constant         =',f8.2)
                if(.not. prm_get_real8_by_key('excluded_force', fk_fix   )) then
                        write(*,'(a)') 'Excluded atoms force constant set to default'
                        fk_fix = fk_fix_default
                end if
                if(fk_fix > 0) then
                        write(*,48) fk_fix
                else
                        fk_fix = fk_fix_default
                        write(*,'(a)')'Shell restraint force constant can not be less-equal zero, set to default'
                        write(*,48) fk_fix
                end if
48              format('Excluded atom restraint force constant         =',f8.2)

                yes = prm_get_logical_by_key('excluded_freeze', freeze, .false.)
                if(freeze) then
                        write(*,'(a)') &
                                'Excluded atoms will not move.'
                end if

                yes = prm_get_logical_by_key('exclude_bonded', exclude_bonded, .false.)
                if(exclude_bonded) then
                        write(*,'(a)') &
                                'Bonded interactions outside the sphere will be eliminated'
                end if
       end if

        ! --- solvent 
        inlog = prm_open_section('solvent')
        if(.not. inlog) inlog = prm_open_section('water') !try also the old name
        if(.not. inlog) then       !defaults
                fk_wsphere = -1
                Dwmz = -1
                awmz = -1
                wpol_restr = wpol_restr_default
                wpol_born = wpol_restr_default
                fkwpol = -1 
        else
                if(prm_get_real8_by_key('radius', rwat_in)) then
                        write(*,'(a,f8.2)') 'Target solvent radius =',rwat_in
                end if
                if(prm_get_line_by_key('centre', instring)) then
                        write(*,30) 'centre'
                end if
                if(prm_get_real8_by_key('pack', rjunk)) then
                        write(*,30) 'pack'
                end if


          if(.not. prm_get_real8_by_key('radial_force', fk_wsphere)) then
                write(*,'(a)') 'Solvent radial restraint force constant set to default'
                fk_wsphere = -1 ! this will be set in water_sphere, once target radius is known
          end if
          yes=prm_get_logical_by_key('polarisation', wpol_restr, wpol_restr_default)
          !default is on when pol. restr is on, otherwise off
          yes=prm_get_logical_by_key('charge_correction', wpol_born, wpol_restr)
          if(wpol_born .and. .not. wpol_restr) then
                write(*,'(a)') '>>> ERROR: charge_correction on requires polarisation on (section solvent)'
                initialize = .false.
          end if
          if(.not. prm_get_real8_by_key('polarisation_force', fkwpol)) then
                write(*,'(a)') 'Solvent polarisation force constant set to default'
                fkwpol = -1 ! this will be set in water_sphere, once target radius is known
          end if
          yes = prm_get_real8_by_key('morse_depth', Dwmz, -1._8)			
          yes = prm_get_real8_by_key('morse_width', awmz, -1._8)			
          if(prm_get_string_by_key('model', instring)) then
                write(*,30) 'model'
          end if
        end if !if (.not. inlog)
end if !if( .not. box )


if(.not. prm_open_section('intervals')) then
        write(*,'(a)') 'non-bond list update interval set to default.'
        NBcycle = NB_cycle_default
        write(*,'(a)') 'energy summary interval set to default.'
        iout_cycle = iout_cycle_default
        itemp_cycle = iout_cycle_default
        iene_cycle = 0 !no energy
        itrj_cycle = 0 !no trajectory

        ivolume_cycle = ivolume_cycle_default


else
        if(.not. prm_get_integer_by_key('non_bond', NBcycle)) then
                write(*,'(a)') 'non-bond list update interval set to default.'
                NBcycle = NB_cycle_default
        end if
        if(.not. prm_get_integer_by_key('output', iout_cycle)) then
                write(*,'(a)') 'energy summary interval set to default.'
                iout_cycle = iout_cycle_default
        end if
        if(.not. prm_get_integer_by_key('temperature', itemp_cycle)) then
                write(*,'(a)') 'temperature print-out interval set to default.'
                itemp_cycle = iout_cycle_default
        end if
        yes = prm_get_integer_by_key('energy', iene_cycle, 0)
        yes = prm_get_integer_by_key('trajectory', itrj_cycle, 0)

        if( constant_pressure ) then
                if( .not. prm_get_integer_by_key('volume_change', ivolume_cycle) ) then
                        write(*,'(a)') 'volume change intervall set to default'
                        ivolume_cycle = ivolume_cycle_default
                end if
        end if
end if

write(*,84) NBcycle
84	format('Non-bonded pair list update interval   =',i8)
86	format('Energy summary print-out interval      =',i8)
87	format('Temperature print-out interval         =',i8)
88	format('Trajectory write interval              =',i8)
89	format('Energy file write interval             =',i8)
83  format('Volume change interval                 =',i8)

if(iout_cycle > 0) then
        write (*,86) iout_cycle
else
        write(*,'(a)') 'No energy summaries written.'
        iout_cycle = -999999999 ! make sure mod(istep, iout_cycle) never = 0
end if
if(itemp_cycle > 0) then
        write (*,87) itemp_cycle
else
        write(*,'(a)') 'No temperatures written.'
        itemp_cycle = -999999999 ! make sure mod(istep, itemp_cycle) never = 0
end if
if(itrj_cycle > 0) then
        write (*,88) itrj_cycle
else
        itrj_cycle = -999999999 !no energy
        write(*,'(a)') 'No trajectory written.'
end if
if(iene_cycle > 0) then
        write (*,89) iene_cycle
else
        iene_cycle = -999999999 !no energy
        write(*,'(a)') 'No energy file written.'
end if
if( constant_pressure ) then
        write(*,83) ivolume_cycle
end if

!read trajectory atom mask
mask_rows = prm_count('trajectory_atoms')
if(itrj_cycle > 0) then
        if(mask_rows == 0) then
                write(*,'(a)') 'All atoms will be included in the trajectory.'
                yes = trj_store_mask('all')
        else
                do i=1,mask_rows
                        yes = prm_get_line(text)
                        yes = trj_store_mask(text)
                end do
        end if
elseif(mask_rows == 0) then
        write(*,'(a)') 'Ignoring section trajectory_atoms.'
end if

if(.not. prm_open_section('files')) then
        write(*,'(a)') '>>> ERROR: files section not found.'
        initialize = .false.
else
        if(.not. prm_get_string_by_key('topology', top_file)) then
                write(*,'(a)') '>>> ERROR: topology not specified (section files)'
                initialize = .false.
        end if
        write (*,60) trim(top_file)
60		format ('Topology file      = ',a)

        if(.not. prm_get_string_by_key('restart', restart_file)) then
                restart = .false.
                if(need_restart) then
                        write(*,'(a)') '>>> ERROR: Restart file required when initial temp. not given.'
                        initialize = .false.
                end if
        else
                restart = .true.
        end if

        if(restart) then
                write (*,65) trim(restart_file)
        else
                write (*,'(a)') 'Initial coordinates taken from topology.'
                if(iseed == 0) then
                        write(*,'(a)') &
                                '>>> ERROR: Need a random number seed to generate initial velocities, aborting.'
                        initialize = .false.
                end if
        end if
65		format ('Initial coord. file= ',a)

        if(.not. prm_get_string_by_key('final', xfin_file)) then
                write(*,'(a)') '>>> ERROR: final co-ordinate file not specified (section files, keyword final)'
                initialize = .false.
        end if
        write (*,80) trim(xfin_file)
80		format ('Final coord. file  = ',a)

        if(.not. prm_get_string_by_key('trajectory', trj_file)) then
                if(itrj_cycle > 0) then
                        write(*,'(a)') '>>> ERROR: Trajectory file name required to write trajectory!'
                        initialize = .false.
                end if
        else
                if(itrj_cycle < 0) then
                        write(*,*) '>>> Error: Trajectory file given but no output interval'
                        initialize = .false.
                end if
                if(itrj_cycle > 0) write (*,90) trim(trj_file)
        end if
90		format ('Trajectory file    = ',a)

        if(.not. prm_get_string_by_key('energy', ene_file)) then
                if(iene_cycle > 0) then
                        write(*,'(a)') '>>> ERROR: Energy file name required to write energies!'
                        initialize = .false.
                end if
        else

                if(iene_cycle < 0) then

                        write(*,'(a)') '>>> ERROR: Energy file given but no energy interval'

                        initialize=.false.

                end if
                if(iene_cycle > 0) write (*,94) trim(ene_file)
        end if
94		format ('Energy output file = ',a)

        if(.not. prm_get_string_by_key('fep', fep_file)) then
                write(*,'(a)') 'No FEP file.'
                !initialize = .false. !This condition IS OK.
                fep_file = ''
        else
                write (*,95) trim(fep_file)
95			format ('FEP input file     = ',a,/)
        end if
        if(.not. prm_get_string_by_key('restraint', exrstr_file)) then
                implicit_rstr_from_file = 0
        else
                implicit_rstr_from_file = 1
                write (*,104) trim(exrstr_file)
104			format ('External rstr file = ',a,/)
        end if
        if(prm_get_string_by_key('water', instring)) then
                write(*,30) 'water'
        end if
end if			

! --- states, EQ
nstates = 0
if(prm_open_section('lambdas')) then
        do while(prm_get_field(instring))
                nstates = nstates + 1
                read(instring, *, iostat=fstat) lamda_tmp(nstates)
                if(fstat /= 0) then
                        write(*,'(a)') '>>> ERROR: Invalid lambda value.'
                        initialize = .false.
                        exit
                end if
        end do
end if
if(nstates == 0 .and. fep_file /= '') then
        if(fep_file /= '') then
                write(*,'(a)') 'Defaulting to single FEP state.'
                nstates = 1
                lamda_tmp(1) = one
        end if
end if
if(nstates > 0 ) then
        if(fep_file == '') then
                write(*,'(a)') '>>> ERROR: FEP file required to use lambdas!'
                initialize = .false.
        else
                ! allocate memory for EQ
                allocate(EQ(nstates), stat=alloc_status)
                call check_alloc('Q-atom energy array')

                ! init EQ%lambda
                EQ(1:nstates)%lambda = lamda_tmp(1:nstates)
                write (*,98) (EQ(i)%lambda,i=1,nstates)
98			format ('lambda-values      = ',10f8.5)
        end if
end if
!Option to make additional calculation with atom groups excluded from the
!energy calculation to provide 'real' group contribution
!Added Paul Bauer 2014
use_excluded_groups = .false.
ngroups_gc = prm_count('group_contribution')
write (*,'(/,a)') 'List of excluded atom groups'
if ( ngroups_gc .gt. 0 ) then
!allocate new EQ arrays for each group
777 format ('Number of excluded group calculations =',i4)
	write (*,777) ngroups_gc
	allocate(ST_gc(ngroups_gc), stat=alloc_status)
	call check_alloc('gc param store')
778 format ('Groupnumber ',i4,'	contains')
	do i=1,ngroups_gc
		write(*,778,advance='no') i
		!read numbers from string
		number = 0
		do while (prm_get_field(instring))
			number = number + 1
			read(instring, *, iostat=fstat) gc_mask_tmp(number)
	                if(fstat /= 0) then
        	                write(*,'(a)') '>>> ERROR: Invalid mask.'
                	        initialize = .false.
                        	exit
                	end if
		end do

		mask_rows = number - 2
		ST_gc(i)%seltype = trim(gc_mask_tmp(1))
		ST_gc(i)%caltypen = trim(gc_mask_tmp(2))
		if ((trim(ST_gc(i)%seltype) /= 'atom' ).and. (trim(ST_gc(i)%seltype) /= 'residue')) then
779 format (/,'Could not understand name of the selection , ',a)
			write (*,779) trim(ST_gc(i)%seltype)
                	initialize = .false.
			exit
		end if
		if ((trim(ST_gc(i)%caltypen) /= 'full').and.(trim(ST_gc(i)%caltypen) /= 'electro') .and. &
			(trim(ST_gc(i)%caltypen) /= 'vdw').and.(trim(ST_gc(i)%caltypen) /= 'all')) then
781 format (/,'Could not understand name of the calculation type , ',a)
			write (*,781) trim(ST_GC(i)%caltypen)
			initialize = .false.
			exit
		end if
		if(iene_cycle > 0) then
        	if(mask_rows == 0) then
!You are not as funny as you think ...
                	write(*,780)
780 format ('no atoms')
                	ST_gc(i)%count=0
        	else
			ST_gc(i)%count=mask_rows
			allocate(ST_gc(i)%maskarray(mask_rows),stat=alloc_status)
			call check_alloc('gc maskarray')
			
	                do ii=3,mask_rows+2
                	        yes = gc_store_mask(ST_gc(i),gc_mask_tmp(ii))
	                end do
			use_excluded_groups = .true.
			write (str,'(i4)') i
782 format (i4,' atom groups')
                        write(*,782) mask_rows 
	        end if
		elseif(mask_rows == 0) then
		        write(*,'(a)') 'Ignoring section for excluding atoms.'
		end if
	end do

end if

!new section in *inp files to control QCP behaviour
!will only trigger if FEP file is in use -> if more than 0 states
if(nstates > 0 ) then
	if(.not. prm_open_section('QCP')) then
		write(*,'(a)') 'No QCP section found, will not try to use RPMD.'
		use_QCP = .false.
	else
		write(*,'(a)') 'Found QCP section, will use RPMD to describe atoms in Q region.'
		use_QCP = .true.
!chose which atoms should be treated as ring polymers
!important later when setting up NB list, RP will be treated different from classical
!this section can be overwritten in the FEP file
		if(.not. prm_get_string_by_key('selection', QCP_select)) then
			write(*,'(a)') 'Will default to Hydrogen atoms only treated as RP!'
			QCP_enum = QCP_HYDROGEN
		else
			QCP_select = upcase(QCP_select)
			 if ((QCP_select .eq. 'HYDROGEN') .or. &
				(QCP_select .eq. 'HYD') .or. &
				(QCP_select .eq. 'H')) then
				QCP_enum = QCP_HYDROGEN
				write(*,'(a)') 'Only treat Hydrogen atoms as RP!'
			elseif ((QCP_select .eq. 'ALL') .or. &
					(QCP_select .eq. 'QATOM') .or. &
					(QCP_select .eq. 'FEP')) then
				QCP_enum = QCP_ALLATOM
				write(*,'(a)') 'Treat all Q-atoms as RP!'
			elseif (QCP_select .eq. 'INDIVIDUAL') then
				QCP_enum = QCP_FEPATOM
				write(*,'(a)') 'Will use information in FEP file to select QCP atoms'
			else
				write(*,'(a)') ' >>> ERROR: No such QCP atom selection!'
				initialize = .false.
			end if
		end if
!how large should the RP be?
!can again be overwritten in FEP file for each atom itself
		if(.not. prm_get_string_by_key('size',QCP_size)) then
			write(*,'(a)') 'Will use default sizes for RP, xxx ring beads per atom!'
			QCP_size_enum = QCP_SIZE_DEFAULT
		else
			QCP_size = upcase(QCP_size)
			if(QCP_size .eq. 'DEFAULT') then
				write(*,'(a)') 'Will use default sizes for RP, xxx ring beads per atom!'
				QCP_size_enum = QCP_SIZE_DEFAULT
			else if (QCP_size .eq. 'SMALL') then
				write(*,'(a)') 'Will use small RP, xxx ring beads per atom!'
				QCP_size_enum = QCP_SIZE_SMALL
			else if (QCP_size .eq. 'INDIVIDUAL') then
				write(*,'(a)') 'Will use information in FEP file to build up ring polymer sizes!'
				QCP_size_enum = QCP_SIZE_INDIVIDUAL
			else if (QCP_size .eq. 'LARGE') then
				write(*,'(a)') 'Will use large RP, xxx beads per atom!'
				QCP_size_enum = QCP_SIZE_LARGE
			else
				write(*,'(a)') ' >>> ERROR: No such QCP size selection!'
				initialize = .false.
			end if
		end if
!number of PI steps at each calculation
		if(.not. prm_get_intg_by_key('steps',QCP_steps)) then
			write(*,'(a)') 'Will use default number of PI steps, n = 10!'
			QCP_steps = QCP_steps_default
		else
			if(QCP_steps .lt. 1) then
				write(*,'(a)') ' >>> ERROR: Can not use less than 1 PI step per classical step!'
				initialize = .false.
			end if
			write(*,'(a,i4,a)') 'Will use the follwoing number of PI steps, n = ',QCP_steps,' for each calculation!'
		end if
	end if
else
	write(*,'(a)') 'No RPMD in classical MD'
end if
!	--- restraints:
write (*,'(/,a)') 'Listing of restraining data:'

! --- nrstr_seq, [rstseq]
nrstr_seq = prm_count('sequence_restraints')
109 format (/,'No. of sequence restraints =',i10)
if ( nrstr_seq .gt. 0 ) then
        ! allocate memory for rstseq
        write (*,109) nrstr_seq
        allocate(rstseq(nrstr_seq), stat=alloc_status)
        call check_alloc('restraint list')
        write (*,110)
110		format ('  atom_i  atom_j      fc  H-flag to_centre')
        do i=1,nrstr_seq
                ! read rstseq(i)
                yes = prm_get_line(text)
                rstseq(i)%to_centre = 0 
                read(text,*, end=111, err=111) rstseq(i)
111			write(*,112) rstseq(i)
112			format (2i8,f8.2,i8,i10)
  end do
end if

! --- nrstr_pos, [rstpos]
nrstr_pos = prm_count('atom_restraints')
115 format (/,'No. of position restratints =',i10)
if ( nrstr_pos .gt. 0 ) then
        write (*,115) nrstr_pos
        ! allocate memory for rstpos
        allocate(rstpos(nrstr_pos), stat=alloc_status)
        call check_alloc('restraint list')
        write (*,120)
120		format ('atom_i      x0      y0      z0     fcx     fcy     fcz   state')
        do i=1,nrstr_pos ! read rstpos(i)
                yes = prm_get_line(text)
                read(text,*, iostat=fstat) rstpos(i)%i,(rstpos(i)%x(j),j=1,3), &
                        (rstpos(i)%fk(j),j=1,3), rstpos(i)%ipsi
                if(fstat /= 0) then
                        write(*,'(a)') '>>> ERROR: Invalid atom restraint data.'
                        initialize = .false.
                        exit
                end if
                write (*,122) rstpos(i)%i,(rstpos(i)%x(j),j=1,3), &
                        (rstpos(i)%fk(j),j=1,3), rstpos(i)%ipsi
        end do
122		format (i6,6f8.2,i8)
end if

! --- nrstr_dist, [rstdis]
nrstr_dist = prm_count('distance_restraints')
125	format (/,'No. of distance restraints =',i10)
if ( nrstr_dist .gt. 0 ) then
        write (*,125) nrstr_dist
        ! allocate memory for rstdis
        allocate(rstdis(nrstr_dist), stat=alloc_status)
        call check_alloc('restraint list')
        write (*,130)
130		format ('atom_i atom_j   dist1   dist2   fc        state') 
        do i=1,nrstr_dist
                yes=prm_get_line(text)
                ! read rstdis(i)
                if(scan(text, ':') > 0) then !got res:atnr
                  !Store in i&j as res:atnr and assign atom nr after topology is read (prep_coord)
                  read(text,*, iostat=fstat) rstdis(i)%itext,rstdis(i)%jtext,rstdis(i)%d1,& 
                     rstdis(i)%d2, rstdis(i)%fk, rstdis(i)%ipsi
                else !Plain numbers
                  read(text,*, iostat=fstat) rstdis(i)%i,rstdis(i)%j,rstdis(i)%d1,&
                        rstdis(i)%d2, rstdis(i)%fk, rstdis(i)%ipsi
                  rstdis(i)%itext = 'nil'
                  rstdis(i)%jtext = 'nil'
                end if
                if(fstat /= 0) then
                  write(*,'(a)') '>>> ERROR: Invalid distance restraint data.'
                  initialize = .false.
                  exit
                end if
                write (*,132) rstdis(i)%i,rstdis(i)%j,rstdis(i)%d1,rstdis(i)%d2,rstdis(i)%fk, &
                        rstdis(i)%ipsi
        end do
132		format (i6,1x,i6,3f8.2,i8)
end if

! --- nrstr_angl, [rstang]
nrstr_angl = prm_count('angle_restraints')
135     format (/,'No. of angle restraints =',i10)
if ( nrstr_angl .gt. 0 ) then
        write (*,135) nrstr_angl
        ! allocate memory for rstang
        allocate(rstang(nrstr_angl), stat=alloc_status)
        call check_alloc('restraint list')
        write (*,140)
140             format ('atom_i atom_j atom_k   angle   fc        state')
        do i=1,nrstr_angl
                yes=prm_get_line(text)
                ! read rstang(i)
                !if(scan(text, ':') > 0) then !got res:atnr
                  !Store in i&j as res:atnr and assign atom nr after topology is
                  !read (prep_coord)
                !  read(text,*, iostat=fstat) rstang(i)%itext,rstang(i)%jtext,rstang(i)%ktext,&
                !     rstang(i)%ang, rstang(i)%fk, rstang(i)%ipsi
                !else !Plain numbers
                  read(text,*, iostat=fstat) rstang(i)%i,rstang(i)%j,rstang(i)%k,&
                        rstang(i)%ang, rstang(i)%fk, rstang(i)%ipsi
                 ! rstang(i)%itext = 'nil'
                 ! rstang(i)%jtext = 'nil'
                 ! rstang(i)%ktext = 'nil'
                !end if
                if(fstat /= 0) then
                  write(*,'(a)') '>>> ERROR: Invalid angle restraint data.'
                  initialize = .false.
                  exit
                end if
                write (*,142) rstang(i)%i,rstang(i)%j,rstang(i)%k,rstang(i)%ang,rstang(i)%fk, &
                        rstang(i)%ipsi
        end do
142             format (i6,1x,i6,1x,i6,2f8.2,i8)
end if


if (.not. box )then
! --- nrstr_wall, [rstwal]
nrstr_wall = prm_count('wall_restraints')
145 format (/,'No. of wall sequence restraints=',i10)
if ( nrstr_wall .gt. 0) then
        write (*,145) nrstr_wall
        ! allocate memory for rstwal
        allocate(rstwal(nrstr_wall), stat=alloc_status)
        call check_alloc('restraint list')
        write (*,150)
150  format ('atom_i atom_j   dist.      fc  aMorse  dMorse  H-flag')
        do i=1,nrstr_wall
                ! read rstwal(:)
                yes = prm_get_line(text)
                read(text,*, iostat=fstat) rstwal(i)%i,rstwal(i)%j,rstwal(i)%d,rstwal(i)%fk, &
                        rstwal(i)%aMorse, rstwal(i)%dMorse, rstwal(i)%ih
                if(fstat /= 0) then
                        write(*,'(a)') '>>> ERROR: Invalid wall restraint data.'
                        initialize = .false.
                        exit
                end if
                write (*,152) rstwal(i)     
        end do
152  format (i6,1x,i6,4f8.2,i8)
end if
end if

call prm_close
end function initialize


!-----------------------------------------------------------------------------------


logical function old_initialize(fu)
!arguments
integer						::	fu
! local variables
integer						::	iuse_indip, shake_flag
character					:: text*200, watmodel*80
integer						:: i,j,length
integer						:: irestart
real(kind=prec)						:: stepsize
real(kind=prec)						:: lamda_tmp(max_states)
integer						:: fstat
integer						::	NBMethod
integer						::	iwpol_restr
real(kind=prec)						::	rjunk

!this is called by initialize to read old-style input file which is 
!alreadu open as unit fu

! this subroutine will init:
!  nsteps, stepsize, dt
!  Temp0, tau_T, iseed, Tmaxw
!  usr_LRF, NBcycle, Rcpp, Rcww, Rcpw, Rcq
!  shake_solvent, shake_solute, shake_hydrogens
!  fk_pshell


!  fk_wsphere=-1, wpol_restr, wpol_born fkwpol=-1, Dwmz=-1, awmz=-1
!    (values initialized to -1 will be set in water_sphere, 
!    once target radius is known)
!  top_file
!  restart, [restart_file]
!  xfin_file
!  itrj_cycle, iene_cycle, iout_cycle, itemp_cycle [trj_file], [ene_file]
!  fep_file
!  nstates, EQ (allocating memory for EQ)
!  implicit_rstr_from_file, [exrstr_file]
!  nrstr_seq, [rstseq] (allocating memory for rstseq)
!  nrstr_pos, [rstpos] (allocating memory for rstpos)
!  nrstr_dist, [rstdis] (allocating memory for rstdis)
!  nrstr_wall, [rstwal] (allocating memory for rstwal)


write(*,1)
1	format('>>> WARNING: Entering unsupported compatibility mode',/ &
           '             to read version 2 input file.',/,&
           '             New features unavailable.')

old_initialize = .true. 

!use default values for new features not in old kind of input.
RcLRF = 999.0_prec
exclude_bonded = .false.
force_rms = .false.
shake_hydrogens = .false.
itemp_cycle = iout_cycle_default
awmz = -1

! --- nsteps, stepsize
!read (fu,*, iostat=stat) nsteps,stepsize
if (.not. prm_get_int_real8(nsteps,stepsize)) then
        old_initialize = .false.
		call die("Wrong input format.")
end if

write (*,10) nsteps, stepsize
10	format ('Number of MD steps =',i10,'  Stepsize (fs)    =',f10.3)

! convert to internal time units once and for all.
dt=0.020462_prec*stepsize
dt2=0.5_prec*dt

! --- Temp0, tau_T, iseed, Tmaxw
read(fu,'(a80)') text !read line into buffer
!now read buffer (avoid reading more lines from inpunt in search for more values)
read(text,*, err=17, end=17) Temp0,tau_T, iseed,Tmaxw 
17	write (*,15) Temp0,tau_T
if(Temp0 <= 0) then
        write(*,'(a)') &
                '>>> Error: No dynamics at zero temperature! Aborting.'
        old_initialize = .false.
end if

15	format ('Target temperature =',f10.2,'  T-relax time     =',f10.2)
if (iseed > 0) write (*,16) Tmaxw, iseed
16	format ('Initial velocities will be generated from Maxwell distribution:',&
                /,'Maxwell temperature=',f10.2,' Random number seed=',i10)
tau_T=0.020462_prec*tau_T
if(tau_T < dt) then
        write(*,'(a)') '>>> Error: tau_t must be >= stepsize.'
        old_initialize = .false.
end if

! --- NBmethod, NBcycle, Rcpp, Rcww, Rcpw, Rcq
read (fu,*) NBmethod,NBcycle,Rcpp,Rcww,Rcpw,Rcq
if(NBMethod == 2) then
        use_LRF = .true.
else 
        use_LRF = .false.
end if
write (*,20) NBmethod,NBcycle
20	format ('Nonbonded method   =',i10,'  NB update cycle  =',i10,/)
write (*,25) Rcpp,Rcww,Rcpw,Rcq
25	format ('Cutoffs are: Rcpp  =',f6.2,'  Rcww =',f6.2,'  Rcpw =', &
f6.2,'  Rcqp =',f6.2,/)

! --- shake_flag
read (fu,*) shake_flag
shake_solvent = .false.
shake_solute = .false.
if(shake_flag >= 1) then
        shake_solvent = .true.
end if
if(shake_flag == 2) then
        shake_solute = .true.
end if

write (*,30) shake_flag
30	format ('Shake method       =',i10)

! --- iuse_indip
read (fu,*) 
write (*,35) 
35	format ('Ignoring induced dipole flag.')

! --- protein center: xpcent(:)
read (fu,*) 
write (*,40)
40	format ('Ignoring solute centre.')

! --- rexcl_o, rexcl_i, fk_pshell
read (fu,*) rjunk, rjunk, fk_pshell
write(*,44)
write (*,45) fk_pshell
44	format ('Ignoring exclusion and shell radii.')
45	format ('Restrained shell force const.          =',f8.2)

! --- water center: xwcent(:)
read (fu,*) 
write (*,50) 
50	format ('Ignoring solvent centre.')

! set default values before reading
! done this way because the SGI compiler initialises values to be read to zero

read(fu,'(a80)') text ! read line into buffer
! now read buffer (avoid reading more lines from input in search for more values)
read(text, fmt=*, err=58, end=58) rjunk, rjunk, fk_wsphere, iwpol_restr, fkwpol, Dwmz
goto 59

! set default values:
58	Dwmz = -1
if (fkwpol .eq. 0) then
  fkwpol = -1
        if (fk_wsphere .eq. 0) then
          fk_wsphere = -1
  end if
end if
59	if(iwpol_restr == 0) then
        wpol_restr = .false.
        wpol_born = .false.
elseif(iwpol_restr == 1) then
        wpol_restr = .true.
        wpol_born = .true.
elseif(iwpol_restr == 2) then
        wpol_restr = .true.
        wpol_born = .false.
else
        call die('unknown water polarisation restraining mode')
end if
write(*,57)
57	format('Ignoring solvent radius and min. packing distance.')
! --- top_file
read (fu,'(a80)') text
call get_fname (text,length,top_file)
write (*,60) top_file(1:length)
60	format ('Topology file      = ',a)

! --- restart, [restart_file]
read (fu,*) irestart
if ( irestart .eq. 1 ) then
        restart = .true.
        read (fu,'(a80)') text
        call get_fname (text,length,restart_file)
        write (*,65) restart_file(1:length)
else
        restart = .false.
        write (*,'(a)') 'Initial coordinates taken from topology.'
        if(iseed == 0) then
                write(*,'(a)') &
                        'Error: Need a random number seed to generate initial velocities, aborting.'
                call die('invalid data in input')
        end if
end if
65 format ('Initial coord. file= ',a)

! --- xfin_file
read (fu,'(a80)') text
call get_fname (text,length,xfin_file)
write (*,80) xfin_file(1:length)
80 format ('Final coord. file  = ',a,/)

! --- itrj_cycle, iene_cycle, iout_cycle, [trj_file], [ene_file]
read (fu,*) itrj_cycle, iene_cycle, iout_cycle
write (*,85) itrj_cycle, iene_cycle, iout_cycle
85 format ('Trajectory, Energy and Output cycles   =',3i8,/)

if ( itrj_cycle .gt. 0 ) then
read (fu,'(a80)') text
call get_fname (text,length,trj_file)
write (*,90) trj_file(1:length)
else
write (*,'(a)') 'No trajectory written.'
 itrj_cycle = -999999999 !make sure mod(istep, itrj_cycle) never = 0
end if
90 format ('Trajectory file    = ',a)

if ( iene_cycle .gt. 0 ) then
read (fu,'(a80)') text
call get_fname (text,length,ene_file)
write (*,94) ene_file(1:length)
else
write (*,'(a)') 'No energy file written'
 iene_cycle = -999999999 ! make sure mod(istep, iene_cycle) never = 0
end if
94 format ('Energy output file = ',a)

if(iout_cycle == 0) then
write(*,'(a)') 'No energy summaries written.'
iout_cycle = -999999999 ! make sure mod(istep, iout_cycle) never = 0
end if

! --- fep_file
read (fu,'(a80)') text
call get_fname (text,length,fep_file)
write (*,95) fep_file(1:length)
95 format ('FEP input file     = ',a,/)

! --- nstates, EQ
read (fu,*) nstates, (lamda_tmp(i),i=1,nstates)
if ( nstates .gt. 0 ) then
! allocate memory for EQ
allocate(EQ(nstates), stat=alloc_status)
call check_alloc('Q-atom energy array')

! init EQ%lambda
EQ(1:nstates)%lambda = lamda_tmp(1:nstates)
write (*,98) (EQ(i)%lambda,i=1,nstates)
98  format ('lambda-values      = ',10f8.5)
end if

!	--- restraints:
write (*,'(/,a)') 'Listing of restraining data:'

! --- implicit_rstr_from_file, [exrstr_file]
read (fu,*) implicit_rstr_from_file
write (*,101) implicit_rstr_from_file
101 format ('Read rstr file     =',i10)
if ( implicit_rstr_from_file .eq. 1 ) then
read (fu,'(a80)') text
call get_fname (text,length,exrstr_file)
write (*,104) exrstr_file(1:length)
else
write (*,105)
end if
104 format ('External rstr file = ',a,/)
105 format ('Implicit positional restraints from topology.',/)

! --- nrstr_seq, [rstseq]
read (fu,*) nrstr_seq
write (*,109) nrstr_seq
109 format (/,'No. sequence rstrs =',i10)
if ( nrstr_seq .gt. 0 ) then
! allocate memory for rstseq
allocate(rstseq(nrstr_seq), stat=alloc_status)
call check_alloc('restraint list')
write (*,110)
110 format (1x,'  atom_i  atom_j      fc  H-flag to_centre')
end if
do i=1,nrstr_seq
! read rstseq(i)
read (fu,'(a80)') text
rstseq(i)%to_centre = 0 
read(text,*, end=111, err=111) rstseq(i)
111	write(*,112) rstseq(i)
112 format (2i8,f8.2,i8,i10)
end do

! --- nrstr_pos, [rstpos]
read (fu,*) nrstr_pos
write (*,115) nrstr_pos
115 format (/,'No. position rstrs =',i10)
if ( nrstr_pos .gt. 0 ) then
! allocate memory for rstpos
allocate(rstpos(nrstr_pos), stat=alloc_status)
call check_alloc('restraint list')
write (*,120)
end if
120 format ('atom_i      x0      y0      z0     fcx     fcy     fcz  istate')
do i=1,nrstr_pos
! read rstpos(i)
read (fu,*) rstpos(i)%i,(rstpos(i)%x(j),j=1,3), &
  (rstpos(i)%fk(j),j=1,3), rstpos(i)%ipsi
write (*,122) rstpos(i)%i,(rstpos(i)%x(j),j=1,3), &
  (rstpos(i)%fk(j),j=1,3), rstpos(i)%ipsi
end do
122 format (i6,6f8.2,i8)

! --- nrstr_dist, [rstdis]
read (fu,*) nrstr_dist
write (*,125) nrstr_dist
125 format ('No. distance rstrs =',i10)
if ( nrstr_dist .gt. 0 ) then
! allocate memory for rstdis
allocate(rstdis(nrstr_dist), stat=alloc_status)
call check_alloc('restraint list')
write (*,130)
end if
130 format ('atom_i atom_j   dist.      fc  istate') 
do i=1,nrstr_dist
! read rstdis(i)
read (fu,*) rstdis(i)%i,rstdis(i)%j,rstdis(i)%d1,rstdis(i)%fk, &
  rstdis(i)%ipsi
        rstdis(i)%d2 = rstdis(i)%d1 !no flat-bottom
write (*,132) rstdis(i)%i,rstdis(i)%j,rstdis(i)%d1,rstdis(i)%fk, &
  rstdis(i)%ipsi
end do
132 format (i6,1x,i6,2f8.2,i8)
! --- nrstr_ang, [rstang]
read (fu,*) nrstr_angl
write (*,135) nrstr_angl
135 format ('No. angle rstrs =',i10)
if ( nrstr_angl .gt. 0 ) then
! allocate memory for rstang
allocate(rstang(nrstr_angl), stat=alloc_status)
call check_alloc('restraint list')
write (*,140)
end if
140 format ('atom_i atom_j atom_k   angle      fc  istate')
do i=1,nrstr_angl
! read rstang(i)
read (fu,*) rstang(i)%i,rstang(i)%j,rstang(i)%k,rstang(i)%ang, &
  rstang(i)%fk,rstang(i)%ipsi
write (*,142) rstang(i)%i,rstang(i)%j,rstang(i)%k,rstang(i)%ang, &
  rstang(i)%fk,rstang(i)%ipsi
end do
142 format (i6,1x,i6,1x,i6,2f8.2,i8)

! --- nrstr_wall, [rstwal]
read (fu,*) nrstr_wall
write (*,145) nrstr_wall
145 format ('No. wall seq. rstrs=',i10)
if ( nrstr_wall .gt. 0) then


! allocate memory for rstwal
allocate(rstwal(nrstr_wall), stat=alloc_status)
call check_alloc('restraint list')
write (*,150)
end if
150 format ('atom_i atom_j   dist.      fc  H-flag')
do i=1,nrstr_wall
! read rstwal(:)
read (fu,*) rstwal(i)%i,rstwal(i)%j,rstwal(i)%d,rstwal(i)%fk, &
  rstwal(i)%ih
write (*,152) rstwal(i)%i,rstwal(i)%j,rstwal(i)%d,rstwal(i)%fk, &
  rstwal(i)%ih     
end do
152 format (i6,1x,i6,2f8.2,i8)

read (fu,'(a80)') text
write (*,157) 
157 format ('Ignoring water file.')

! --- determine water model
read (fu,'(a80)') text
write (*,160) 
160 format ('Ignoring water model.')

end function old_initialize

!-----------------------------------------------------------------------
! function to gather lrf information from nodes and set it the values
! gathered there
#ifdef USE_MPI
integer function lrf_add(lrf1,lrf2,length,MPI_Datatype)
        type(LRF_TYPE),INTENT (IN)      :: lrf1(length)
        type(LRF_TYPE)                  :: lrf2(length)
        integer                         :: length,MPI_Datatype

        lrf2(1:length)%cgp_cent(1)        = lrf1(1:length)%cgp_cent(1)
        lrf2(1:length)%cgp_cent(2)        = lrf1(1:length)%cgp_cent(2)
        lrf2(1:length)%cgp_cent(3)        = lrf1(1:length)%cgp_cent(3)
        lrf2(1:length)%phi0               = lrf1(1:length)%phi0+lrf2(1:length)%phi0
        lrf2(1:length)%phi1(1)            = lrf1(1:length)%phi1(1 )+lrf2(1:length)%phi1(1 )
        lrf2(1:length)%phi1(2)            = lrf1(1:length)%phi1(2 )+lrf2(1:length)%phi1(2 )
        lrf2(1:length)%phi1(3)            = lrf1(1:length)%phi1(3 )+lrf2(1:length)%phi1(3 )
        lrf2(1:length)%phi2(1)            = lrf1(1:length)%phi2(1 )+lrf2(1:length)%phi2(1 )
        lrf2(1:length)%phi2(2)            = lrf1(1:length)%phi2(2 )+lrf2(1:length)%phi2(2 )
        lrf2(1:length)%phi2(3)            = lrf1(1:length)%phi2(3 )+lrf2(1:length)%phi2(3 )
        lrf2(1:length)%phi2(4)            = lrf1(1:length)%phi2(4 )+lrf2(1:length)%phi2(4 )
        lrf2(1:length)%phi2(5)            = lrf1(1:length)%phi2(5 )+lrf2(1:length)%phi2(5 )
        lrf2(1:length)%phi2(6)            = lrf1(1:length)%phi2(6 )+lrf2(1:length)%phi2(6 )
        lrf2(1:length)%phi2(7)            = lrf1(1:length)%phi2(7 )+lrf2(1:length)%phi2(7 )
        lrf2(1:length)%phi2(8)            = lrf1(1:length)%phi2(8 )+lrf2(1:length)%phi2(8 )
        lrf2(1:length)%phi2(9)            = lrf1(1:length)%phi2(9 )+lrf2(1:length)%phi2(9 )
        lrf2(1:length)%phi3(1)            = lrf1(1:length)%phi3(1 )+lrf2(1:length)%phi3(1 )
        lrf2(1:length)%phi3(2)            = lrf1(1:length)%phi3(2 )+lrf2(1:length)%phi3(2 )
        lrf2(1:length)%phi3(3)            = lrf1(1:length)%phi3(3 )+lrf2(1:length)%phi3(3 )
        lrf2(1:length)%phi3(4)            = lrf1(1:length)%phi3(4 )+lrf2(1:length)%phi3(4 )
        lrf2(1:length)%phi3(5)            = lrf1(1:length)%phi3(5 )+lrf2(1:length)%phi3(5 )
        lrf2(1:length)%phi3(6)            = lrf1(1:length)%phi3(6 )+lrf2(1:length)%phi3(6 )
        lrf2(1:length)%phi3(7)            = lrf1(1:length)%phi3(7 )+lrf2(1:length)%phi3(7 )
        lrf2(1:length)%phi3(8)            = lrf1(1:length)%phi3(8 )+lrf2(1:length)%phi3(8 )
        lrf2(1:length)%phi3(9)            = lrf1(1:length)%phi3(9 )+lrf2(1:length)%phi3(9 )
        lrf2(1:length)%phi3(10)           = lrf1(1:length)%phi3(10)+lrf2(1:length)%phi3(10)
        lrf2(1:length)%phi3(11)           = lrf1(1:length)%phi3(11)+lrf2(1:length)%phi3(11)
        lrf2(1:length)%phi3(12)           = lrf1(1:length)%phi3(12)+lrf2(1:length)%phi3(12)
        lrf2(1:length)%phi3(13)           = lrf1(1:length)%phi3(13)+lrf2(1:length)%phi3(13)
        lrf2(1:length)%phi3(14)           = lrf1(1:length)%phi3(14)+lrf2(1:length)%phi3(14)
        lrf2(1:length)%phi3(15)           = lrf1(1:length)%phi3(15)+lrf2(1:length)%phi3(15)
        lrf2(1:length)%phi3(16)           = lrf1(1:length)%phi3(16)+lrf2(1:length)%phi3(16)
        lrf2(1:length)%phi3(17)           = lrf1(1:length)%phi3(17)+lrf2(1:length)%phi3(17)
        lrf2(1:length)%phi3(18)           = lrf1(1:length)%phi3(18)+lrf2(1:length)%phi3(18)
        lrf2(1:length)%phi3(19)           = lrf1(1:length)%phi3(19)+lrf2(1:length)%phi3(19)
        lrf2(1:length)%phi3(20)           = lrf1(1:length)%phi3(20)+lrf2(1:length)%phi3(20)
        lrf2(1:length)%phi3(21)           = lrf1(1:length)%phi3(21)+lrf2(1:length)%phi3(21)
        lrf2(1:length)%phi3(22)           = lrf1(1:length)%phi3(22)+lrf2(1:length)%phi3(22)
        lrf2(1:length)%phi3(23)           = lrf1(1:length)%phi3(23)+lrf2(1:length)%phi3(23)
        lrf2(1:length)%phi3(24)           = lrf1(1:length)%phi3(24)+lrf2(1:length)%phi3(24)
        lrf2(1:length)%phi3(25)           = lrf1(1:length)%phi3(25)+lrf2(1:length)%phi3(25)
        lrf2(1:length)%phi3(26)           = lrf1(1:length)%phi3(26)+lrf2(1:length)%phi3(26)
        lrf2(1:length)%phi3(27)           = lrf1(1:length)%phi3(27)+lrf2(1:length)%phi3(27)

        lrf_add = -1

end function lrf_add

integer function lrf_cgp_rep(lrf1,lrf2,length,MPI_Datatype)
        type(LRF_TYPE),INTENT (IN)      :: lrf1(length)
        type(LRF_TYPE)                  :: lrf2(length)
        integer                         :: length,MPI_Datatype
        lrf2(1:length)%cgp_cent(1)        = lrf1(1:length)%cgp_cent(1)
        lrf2(1:length)%cgp_cent(2)        = lrf1(1:length)%cgp_cent(2)
        lrf2(1:length)%cgp_cent(3)        = lrf1(1:length)%cgp_cent(3)
        lrf_cgp_rep = -1

end function lrf_cgp_rep
#endif
! this one gets the lrf info from all nodes and puts it on master after a lrf
! update has been completed
! only used with MPI
#ifdef USE_MPI
subroutine lrf_gather
call MPI_Allreduce(MPI_IN_PLACE,LRF,ncgp,mpitype_batch_lrf,mpi_lrf_add,&
        MPI_COMM_WORLD,ierr)
if (ierr .ne. 0) call die('lrf_gather reduce')
!if (nodeid .eq.0 ) LRF = lrf_reduce
end subroutine lrf_gather
#endif

subroutine lrf_taylor
! *** local variables
integer						::	i,i3,ic
real(kind=prec)						::	Vij, q
real(kind=prec)						::	dr(3),df(3)

! global variables used:
!  E%LRF, natom, excl, iqatom, iwhich_cgp, lrf, x, crg, d
! change to use the calculation_assignment%natom variable
! to only calculate the stuff each node should do
do i = calculation_assignment%natom%start, calculation_assignment%natom%end
! for every atom on this node:

if ( ( use_PBC .and. (iqatom(i)==0) ) .or. ( (.not. excl(i) ) .and. (iqatom(i)==0) ) ) then
  ! unless excluded atom or q-atom:

  ! find the displacement dr from the center of the charge group
i3 = i*3-3
ic = iwhich_cgp(i)
dr(1) = lrf(ic)%cgp_cent(1) - x(i3+1)
dr(2) = lrf(ic)%cgp_cent(2) - x(i3+2)
dr(3) = lrf(ic)%cgp_cent(3) - x(i3+3)

! --- Electric potential
Vij=lrf(ic)%phi0 &
     +lrf(ic)%phi1(1)*dr(1)+lrf(ic)%phi1(2)*dr(2)+lrf(ic)%phi1(3)*dr(3) &
     +0.5_prec*(lrf(ic)%phi2(1)*dr(1)+lrf(ic)%phi2(2)*dr(2) &
     +lrf(ic)%phi2(3)*dr(3))*dr(1) &
     +0.5_prec*(lrf(ic)%phi2(4)*dr(1)+lrf(ic)%phi2(5)*dr(2) &
     +lrf(ic)%phi2(6)*dr(3))*dr(2) &
     +0.5_prec*(lrf(ic)%phi2(7)*dr(1)+lrf(ic)%phi2(8)*dr(2) &
     +lrf(ic)%phi2(9)*dr(3))*dr(3)

E%LRF = E%LRF + 0.5_prec * crg(i) * Vij

! --- Electric field
df(1)=lrf(ic)%phi1(1) &
     +lrf(ic)%phi2(1)*dr(1)+lrf(ic)%phi2(2)*dr(2)+lrf(ic)%phi2(3)*dr(3) &
     +0.5_prec*(lrf(ic)%phi3(1 )*dr(1)+lrf(ic)%phi3(2 )*dr(2) &
     +lrf(ic)%phi3(3 )*dr(3))*dr(1) &
     +0.5_prec*(lrf(ic)%phi3(4 )*dr(1)+lrf(ic)%phi3(5 )*dr(2) &
     +lrf(ic)%phi3(6 )*dr(3))*dr(2) &
     +0.5_prec*(lrf(ic)%phi3(7 )*dr(1)+lrf(ic)%phi3(8 )*dr(2) &
     +lrf(ic)%phi3(9 )*dr(3))*dr(3)
df(2)=lrf(ic)%phi1(2) &
     +lrf(ic)%phi2(4)*dr(1)+lrf(ic)%phi2(5)*dr(2)+lrf(ic)%phi2(6)*dr(3) &
     +0.5_prec*(lrf(ic)%phi3(10)*dr(1)+lrf(ic)%phi3(11)*dr(2) &
     +lrf(ic)%phi3(12)*dr(3))*dr(1) &
     +0.5_prec*(lrf(ic)%phi3(13)*dr(1)+lrf(ic)%phi3(14)*dr(2) &
     +lrf(ic)%phi3(15)*dr(3))*dr(2) &
     +0.5_prec*(lrf(ic)%phi3(16)*dr(1)+lrf(ic)%phi3(17)*dr(2) &
     +lrf(ic)%phi3(18)*dr(3))*dr(3)
df(3)=lrf(ic)%phi1(3) &
     +lrf(ic)%phi2(7)*dr(1)+lrf(ic)%phi2(8)*dr(2)+lrf(ic)%phi2(9)*dr(3) &
     +0.5_prec*(lrf(ic)%phi3(19)*dr(1)+lrf(ic)%phi3(20)*dr(2) &
     +lrf(ic)%phi3(21)*dr(3))*dr(1) &
     +0.5_prec*(lrf(ic)%phi3(22)*dr(1)+lrf(ic)%phi3(23)*dr(2) &
     +lrf(ic)%phi3(24)*dr(3))*dr(2) &
     +0.5_prec*(lrf(ic)%phi3(25)*dr(1)+lrf(ic)%phi3(26)*dr(2) &
     +lrf(ic)%phi3(27)*dr(3))*dr(3) 

  ! update d
d(i3+1)=d(i3+1)-crg(i)*df(1)
d(i3+2)=d(i3+2)-crg(i)*df(2)
d(i3+3)=d(i3+3)-crg(i)*df(3)
end if
end do
end subroutine lrf_taylor


!-----------------------------------------------------------------------
!DEC$ ATTRIBUTES FORCEINLINE :: lrf_update
subroutine lrf_update(group1,group2)
! --- input variables
integer				:: group1,group2

! --- local variables
integer				:: i,ia,i3,i_sw
real(kind=prec)			:: dr(3),r2,shift(3),field0,field1,field2

! get switch atoms for PBC
if (use_PBC) then
	i_sw = cgp(group1)%iswitch
	i3 = 3*i_sw-3

	shift(1) = x(i3+1) - lrf(group2)%cgp_cent(1)
	shift(2) = x(i3+2) - lrf(group2)%cgp_cent(2)
	shift(3) = x(i3+3) - lrf(group2)%cgp_cent(3)

	shift(1) = boxlength(1)*nint( shift(1)*inv_boxl(1) )
	shift(2) = boxlength(2)*nint( shift(2)*inv_boxl(2) )
	shift(3) = boxlength(3)*nint( shift(3)*inv_boxl(3) )
else
	shift(:) = 0
end if

iloop:        do ia = cgp(group1)%first, cgp(group1)%last
          ! skip if q-atom
          i = cgpatom(ia)
          if ( iqatom(i)/=0 ) cycle iloop

          i3 = i*3-3

          dr(1) = x(i3+1) - lrf(group2)%cgp_cent(1) - shift(1)
          dr(2) = x(i3+2) - lrf(group2)%cgp_cent(2) - shift(2)
          dr(3) = x(i3+3) - lrf(group2)%cgp_cent(3) - shift(3)
          r2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)

          field0=crg(i)/(r2*sqrt(r2))
          field1=3.0_prec*field0/r2
          field2=-field1/r2
!$omp critical
          lrf(group2)%phi0=lrf(group2)%phi0+field0*r2
          lrf(group2)%phi1(1)=lrf(group2)%phi1(1)-field0*dr(1)
          lrf(group2)%phi1(2)=lrf(group2)%phi1(2)-field0*dr(2)
          lrf(group2)%phi1(3)=lrf(group2)%phi1(3)-field0*dr(3)
          lrf(group2)%phi2(1)=lrf(group2)%phi2(1)+field1*dr(1)*dr(1)-field0
          lrf(group2)%phi2(2)=lrf(group2)%phi2(2)+field1*dr(1)*dr(2)
          lrf(group2)%phi2(3)=lrf(group2)%phi2(3)+field1*dr(1)*dr(3)
          lrf(group2)%phi2(4)=lrf(group2)%phi2(4)+field1*dr(2)*dr(1)
          lrf(group2)%phi2(5)=lrf(group2)%phi2(5)+field1*dr(2)*dr(2)-field0
          lrf(group2)%phi2(6)=lrf(group2)%phi2(6)+field1*dr(2)*dr(3)
          lrf(group2)%phi2(7)=lrf(group2)%phi2(7)+field1*dr(3)*dr(1)
          lrf(group2)%phi2(8)=lrf(group2)%phi2(8)+field1*dr(3)*dr(2)
          lrf(group2)%phi2(9)=lrf(group2)%phi2(9)+field1*dr(3)*dr(3)-field0
!$omp end critical
!$omp critical
          lrf(group2)%phi3(1 )=lrf(group2)%phi3(1 ) +field2*(5.0_prec*dr(1)*dr(1)*dr(1)-r2*3.0_prec*dr(1))
          lrf(group2)%phi3(2 )=lrf(group2)%phi3(2 ) +field2*(5.0_prec*dr(1)*dr(1)*dr(2)-r2*dr(2))
          lrf(group2)%phi3(3 )=lrf(group2)%phi3(3 ) +field2*(5.0_prec*dr(1)*dr(1)*dr(3)-r2*dr(3))
          lrf(group2)%phi3(4 )=lrf(group2)%phi3(4 ) +field2*(5.0_prec*dr(1)*dr(2)*dr(1)-r2*dr(2))
          lrf(group2)%phi3(5 )=lrf(group2)%phi3(5 ) +field2*(5.0_prec*dr(1)*dr(2)*dr(2)-r2*dr(1))
          lrf(group2)%phi3(6 )=lrf(group2)%phi3(6 ) +field2*(5.0_prec*dr(1)*dr(2)*dr(3))
          lrf(group2)%phi3(7 )=lrf(group2)%phi3(7 ) +field2*(5.0_prec*dr(1)*dr(3)*dr(1)-r2*dr(3))
          lrf(group2)%phi3(8 )=lrf(group2)%phi3(8 ) +field2*(5.0_prec*dr(1)*dr(3)*dr(2))
          lrf(group2)%phi3(9 )=lrf(group2)%phi3(9 ) +field2*(5.0_prec*dr(1)*dr(3)*dr(3)-r2*dr(1))
          lrf(group2)%phi3(10)=lrf(group2)%phi3(10) +field2*(5.0_prec*dr(2)*dr(1)*dr(1)-r2*dr(2))
          lrf(group2)%phi3(11)=lrf(group2)%phi3(11) +field2*(5.0_prec*dr(2)*dr(1)*dr(2)-r2*dr(1))
          lrf(group2)%phi3(12)=lrf(group2)%phi3(12) +field2*(5.0_prec*dr(2)*dr(1)*dr(3))
          lrf(group2)%phi3(13)=lrf(group2)%phi3(13) +field2*(5.0_prec*dr(2)*dr(2)*dr(1)-r2*dr(1))
          lrf(group2)%phi3(14)=lrf(group2)%phi3(14) +field2*(5.0_prec*dr(2)*dr(2)*dr(2)-r2*3.0_prec*dr(2))
          lrf(group2)%phi3(15)=lrf(group2)%phi3(15) +field2*(5.0_prec*dr(2)*dr(2)*dr(3)-r2*dr(3))
          lrf(group2)%phi3(16)=lrf(group2)%phi3(16) +field2*(5.0_prec*dr(2)*dr(3)*dr(1))
          lrf(group2)%phi3(17)=lrf(group2)%phi3(17) +field2*(5.0_prec*dr(2)*dr(3)*dr(2)-r2*dr(3))
          lrf(group2)%phi3(18)=lrf(group2)%phi3(18) +field2*(5.0_prec*dr(2)*dr(3)*dr(3)-r2*dr(2))
          lrf(group2)%phi3(19)=lrf(group2)%phi3(19) +field2*(5.0_prec*dr(3)*dr(1)*dr(1)-r2*dr(3))
          lrf(group2)%phi3(20)=lrf(group2)%phi3(20) +field2*(5.0_prec*dr(3)*dr(1)*dr(2))
          lrf(group2)%phi3(21)=lrf(group2)%phi3(21) +field2*(5.0_prec*dr(3)*dr(1)*dr(3)-r2*dr(1))
          lrf(group2)%phi3(22)=lrf(group2)%phi3(22) +field2*(5.0_prec*dr(3)*dr(2)*dr(1))
          lrf(group2)%phi3(23)=lrf(group2)%phi3(23) +field2*(5.0_prec*dr(3)*dr(2)*dr(2)-r2*dr(3))
          lrf(group2)%phi3(24)=lrf(group2)%phi3(24) +field2*(5.0_prec*dr(3)*dr(2)*dr(3)-r2*dr(2))
          lrf(group2)%phi3(25)=lrf(group2)%phi3(25) +field2*(5.0_prec*dr(3)*dr(3)*dr(1)-r2*dr(1))
          lrf(group2)%phi3(26)=lrf(group2)%phi3(26) +field2*(5.0_prec*dr(3)*dr(3)*dr(2)-r2*dr(2))
          lrf(group2)%phi3(27)=lrf(group2)%phi3(27) +field2*(5.0_prec*dr(3)*dr(3)*dr(3)-r2*3.0_prec*dr(3))
!$omp end critical
        end do iloop

end subroutine lrf_update

!----------------------------------------------------------

subroutine make_pair_lists
#if defined (PROFILING)
real(kind=prec)                                         :: start_loop_time
start_loop_time = rtime()
#endif
! start with populating the grid
! each node populates its assigned ncgp_solute
! and nwat and then intercommunication does the rest

! we make this whole thing now omp parallel
! start opening region here so we don't have to pay the penalty
! to do this in every function

!$omp parallel default(none) shared(use_PBC,use_LRF,iuse_switch_atom) 

#ifdef USE_GRID
call populate_grids
#endif
if( use_PBC ) then
        if (.not. use_LRF)then
		!cutoff
           if(iuse_switch_atom == 1) then
                call nbpplist_box
                call nbpwlist_box
                call nbqplist_box
            else
                call nbpplis2_box
                call nbpwlis2_box
                call nbqplis2_box
            end if
            call nbwwlist_box
        else
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
        endif
		call nbqwlist_box
else !spherical case
        if(.not. use_LRF) then
                ! cutoff
                if( iuse_switch_atom .eq. 1 ) then
                        call nbpplist
                        call nbpwlist
                        call nbqplist
                else
                        call nbpplis2
                        call nbpwlis2
                        call nbqplis2
                end if
                call nbwwlist
        else 
                ! cutoff with lrf
                call cgp_centers ! *** måste anropas av alla noder (nollställer lrf)
                if( iuse_switch_atom .eq. 1 ) then
                        call nbpplist_lrf   
                        call nbpwlist_lrf   
                        call nbqplist       
                else
                        call nbpplis2_lrf
                        call nbpwlis2_lrf
                        call nbqplis2
                end if
                call nbwwlist_lrf   
        end if

        call nbqwlist   
end if

!$omp end parallel

! we are gathering the LRF info here
#ifdef USE_MPI
if (use_LRF) call lrf_gather
#endif
#if defined (PROFILING)
profile(1)%time = profile(1)%time + rtime() - start_loop_time
#endif

end subroutine make_pair_lists

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



subroutine temperature(Tscale_solute,Tscale_solvent,Ekinmax)
! calculate the temperature
!arguments
real(kind=prec)						:: Tscale_solute,Tscale_solvent,Ekinmax

!locals
integer						::	i, i3
real(kind=prec)						::	Ekin

Temp = zero
Temp_solute = zero
Tfree_solute = zero
Texcl_solute = zero

!get kinetic energies for solute atoms
do i=1,nat_solute
        i3=i*3-3
        Ekin = 0.5_prec*iaclib(iac(i))%mass*(v(i3+1)**2+v(i3+2)**2+v(i3+3)**2)
        Temp_solute = Temp_solute + Ekin

        !******PWadded if
        if( use_PBC .or. ( (.not. use_PBC) .and. (.not. excl(i)) ) ) then
                Tfree_solute = Tfree_solute +Ekin
        else
        	Texcl_solute = Texcl_solute +Ekin
        end if
        !if ( .not. excl(i)) Tfree = Tfree + Ekin
        if ( Ekin .gt. Ekinmax ) then
                ! hot atom warning
                write (*,180) i,2.0_prec*Ekin/Boltz/3.0_prec
        end if
end do

Tfree_solvent = zero
Temp_solvent = zero
Texcl_solvent = zero
Ekin = zero


!get kinetic energies for solvent atoms
do i=nat_solute+1,natom
        i3=i*3-3
        Ekin = 0.5_prec*iaclib(iac(i))%mass*(v(i3+1)**2+v(i3+2)**2+v(i3+3)**2)
        Temp_solvent = Temp_solvent + Ekin

        !******PWadded if
        if( use_PBC .or. ( (.not. use_PBC) .and. (.not. excl(i)) ) ) then
                Tfree_solvent = Tfree_solvent +Ekin
        else
        	Texcl_solvent = Texcl_solvent +Ekin
        end if
        !if ( .not. excl(i)) Tfree = Tfree + Ekin
        if ( Ekin .gt. Ekinmax ) then
                ! hot atom warning
                write (*,180) i,2.0_prec*Ekin/Boltz/3.0_prec
        end if
end do

Tfree = Tfree_solvent + Tfree_solute
Temp = Temp_solute + Temp_solvent

E%kinetic = Temp

Temp  = 2.0_prec*Temp/Boltz/real(Ndegf, kind=prec)
Tfree = 2.0_prec*Tfree/Boltz/real(Ndegfree, kind=prec)

if (detail_temps) then
	Temp_solute  = 2.0_prec*Temp_solute /Boltz/real(Ndegf_solute, kind=prec)
	Tfree_solute = 2.0_prec*Tfree_solute/Boltz/real(Ndegfree_solute, kind=prec)
	if ( Ndegf_solute .ne. Ndegfree_solute) Texcl_solute = 2.0_prec*Texcl_solute/Boltz/real(Ndegf_solute - Ndegfree_solute, kind=prec)

	Temp_solvent  = 2.0_prec*Temp_solvent /Boltz/real(Ndegf_solvent, kind =prec)
	Tfree_solvent = 2.0_prec*Tfree_solvent/Boltz/real(Ndegfree_solvent, kind=prec)
	if ( Ndegf_solvent .ne. Ndegfree_solvent) Texcl_solvent = 2.0_prec*Texcl_solvent/Boltz/real(Ndegf_solvent - Ndegfree_solvent,&
			 kind=prec)
end if

if (thermostat == BERENDSEN) then
	if (separate_scaling) then
		if ( Tfree_solvent .ne. 0 ) Tscale_solvent = Temp0/Tfree_solvent - one
		Tscale_solvent = sqrt ( 1 + dt/tau_T * Tscale_solvent )
		if ( Tfree_solute .ne. 0 ) Tscale_solute = Temp0/Tfree_solute - one
		Tscale_solute = sqrt ( 1 + dt/tau_T * Tscale_solute )
	else
		if ( Tfree .ne. 0 ) Tscale_solvent = Temp0/Tfree - one
		Tscale_solvent = sqrt ( 1 + dt/tau_T * Tscale_solvent )
		Tscale_solute = Tscale_solvent
	end if
end if

180 format ('>>> WARNING: hot atom, i =',i10,' Temp(i)=',f10.2)

end subroutine temperature

!----------------------------------------------------------------------
! Subroutines for velocity update, one for each thermostat
subroutine newvel_ber(step,tempscale,startatom,endatom)!dt_mod,Tscale_solute,i,nat_solute)
real(kind=prec)				:: step,tempscale
integer				:: startatom, endatom

!locals
integer				:: i,i3

do i=startatom,endatom
	i3=i*3-3
	v(i3+1)= ( v(i3+1)-d(i3+1)*winv(i)*step ) * tempscale
	xx(i3+1) = x(i3+1)
	x(i3+1) = x(i3+1) + v(i3+1)*dt

	v(i3+2)= ( v(i3+2)-d(i3+2)*winv(i)*step ) * tempscale
	xx(i3+2) = x(i3+2)
	x(i3+2) = x(i3+2) + v(i3+2)*dt

	v(i3+3)= ( v(i3+3)-d(i3+3)*winv(i)*step ) * tempscale
	xx(i3+3) = x(i3+3)
	x(i3+3) = x(i3+3) + v(i3+3)*dt

	if (integrator == VELVERLET) then
		v(i3+1)= ( v(i3+1)-d(i3+1)*winv(i)*step ) * tempscale
		v(i3+2)= ( v(i3+2)-d(i3+2)*winv(i)*step ) * tempscale
		v(i3+3)= ( v(i3+3)-d(i3+3)*winv(i)*step ) * tempscale
	end if

end do

end subroutine newvel_ber

subroutine newvel_lan(step,lfriction,rands,rande)
real(kind=prec)					:: step,lfriction,randnum
real(kind=prec)					:: rands(:),rande(:,:)
!locals
integer				:: i,j,i3


if (.not.langevin_gauss_rand) then
!using uniform random numbers, default
	do i=1,natom
	do j=1,3
	randnum = randm(iseed)
	rande(i,j) =  rands(i) * ( randnum  - 0.5_prec )
	end do 
	end do
else
!using gaussian random numbers, not default
	do i=1,natom
	do j=1,3
	call gauss(zero,rands(i),rande(i,j),iseed)
	end do 
	end do
	end if  

do i=1,natom
	i3=i*3-3
	v(i3+1)  = (v(i3+1)*lfriction) -(d(i3+1)-rande(i,1))*winv(i)*step
	xx(i3+1) = x(i3+1)
	x(i3+1)  = x(i3+1) + v(i3+1)*dt

	v(i3+2)  = (v(i3+2)*lfriction) -(d(i3+2)-rande(i,2))*winv(i)*step
	xx(i3+2) = x(i3+2)
	x(i3+2)  = x(i3+2) + v(i3+2)*dt

	v(i3+3)  = (v(i3+3)*lfriction) -(d(i3+3)-rande(i,3))*winv(i)*step
	xx(i3+3) = x(i3+3)
	x(i3+3)  = x(i3+3) + v(i3+3)*dt

	if (integrator == VELVERLET) then
		v(i3+1)= (v(i3+1)*lfriction) -(d(i3+1)-rande(i,1))*winv(i)*step
		v(i3+2)= (v(i3+2)*lfriction) -(d(i3+2)-rande(i,2))*winv(i)*step
		v(i3+3)= (v(i3+3)*lfriction) -(d(i3+3)-rande(i,3))*winv(i)*step
	end if
end do


end subroutine newvel_lan

subroutine newvel_nos(step)
real(kind=prec)				:: step

!locals
integer				:: i,i3

do i=1,natom
        i3=i*3-3
        v(i3+1)= ( v(i3+1)-d(i3+1)*winv(i)*step )
        xx(i3+1) = x(i3+1)
        x(i3+1) = x(i3+1) + v(i3+1)*dt

        v(i3+2)= ( v(i3+2)-d(i3+2)*winv(i)*step )
        xx(i3+2) = x(i3+2)
        x(i3+2) = x(i3+2) + v(i3+2)*dt

        v(i3+3)= ( v(i3+3)-d(i3+3)*winv(i)*step )
        xx(i3+3) = x(i3+3)
        x(i3+3) = x(i3+3) + v(i3+3)*dt

        if (integrator == VELVERLET) then
                v(i3+1)= ( v(i3+1)-d(i3+1)*winv(i)*step )
                v(i3+2)= ( v(i3+2)-d(i3+2)*winv(i)*step )
                v(i3+3)= ( v(i3+3)-d(i3+3)*winv(i)*step )
        end if

end do

end subroutine newvel_nos
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
!******PWchanged 2002-10-01
subroutine md_run

! local variables
integer				:: i,j,k,niter,iii

integer				:: i3
real(kind=prec)                         :: Tlast
real(kind=prec)                         ::Ekinmax
real(kind=prec)				::Tscale_solute,Tscale_solvent
!new for temperature control
real(kind=prec)				::Tscale,dv_mod,dv_friction,dv_mod2,randnum
integer				:: n_max = -1
real(kind=prec),allocatable	::randva(:),randfa(:,:)
!Random variable temperature control array
real(kind=prec)				:: time0, time1, time_per_step, startloop
integer(4)				:: time_completion

!Local variables for file writing to energy files
integer				:: loc_arrays
#if defined(PROFILING)
real(kind=prec)                                         :: start_loop_time1, start_loop_time2

profile(1)%name = 'NB_update'
profile(2)%name = '   nbwwlist_time'
profile(3)%name = '   nbpplist_time'
profile(4)%name = '   nbpwlist_time'
profile(5)%name = '   nbqplist_time'
profile(6)%name = '   nbqwlist_time'
profile(7)%name = 'SHAKE'
profile(8)%name = 'Bonded Terms'
profile(9)%name = 'Restraints'
profile(10)%name = 'Nonbonded Terms'
profile(11)%name = 'Update vel. & coords.'


#endif

#if defined(PROFILING)
#if defined(USE_MPI)
if (nodeid .eq. 0) then
	allocate(all_node_times(num_profiling_times*numnodes), stat=alloc_status) !vector for storing all node's node_times, used by mpi_gather at end of md_run
	call check_alloc('MPI profiling')
end if
allocate(node_times(num_profiling_times), stat=alloc_status) !each node's profiling times, used at end of md_run by mpi_gather
call check_alloc('MPI profiling')

all_node_times(:) = zero
node_times(:) = zero

#endif
#endif

!Define array length for EQ array in local usage
loc_arrays=ene_header%arrays

!Define number of coord to send/recieve
nat3=natom*3

! calculate maximum temperature
!**MN-> Only master calc. temp for now.
if (nodeid .eq. 0) then
Ekinmax = 1000.0_prec*Ndegf*Boltz*Temp0/2.0_prec/real(natom, kind=prec)

call temperature(Tscale_solute,Tscale_solvent,Ekinmax)
!store old Temp
Tlast = Temp
end if

if (nodeid .eq. 0) then
        ! master node only: print initial temperatures
        write (*,*)
        write (*,120) 'Initial', Temp, Tfree
		if ( detail_temps ) then
			write (*,120) 'Solvent', Temp_solvent, Tfree_solvent
			write (*,120) 'Solute', Temp_solute, Tfree_solute
!			write (*,120) 'Excl solute, solvent', Texcl_solute, Texcl_solvent
		end if
120		format(a7,' temperatures are : Ttot =',f10.2,' Tfree =',f10.2)
        write (*,*)

        ! init timer
        time0 = rtime()

		! Init timer of total loop time
		startloop = rtime()
end if

!set thermostat variables to the values needed in the time steps
!moved here to have definition before we do every step
if (nodeid.eq.0) then
!	allocate temp control arrays
	if ( thermostat == LANGEVIN ) then
	allocate(randva(natom))
	call check_alloc('Random variable temperature control array')
	allocate(randfa(natom,3))
	call check_alloc('Random variable temperature control array 2')
	end if
! if velocity verlet is used, take half steps
if ( integrator == VELVERLET ) then
        dv_mod = dt2
else
        dv_mod = dt
end if
if ( thermostat == LANGEVIN ) then
	dv_friction = (1.0_prec-friction*dv_mod)
	if (.not.langevin_gauss_rand) then
!we use the fast random number generation from unifrom random numbers
!this is the default! -- Paul and Alex October 2014
	gkT = 24.0_prec*friction*Boltz*Temp0/dv_mod ! this 24 comes from the fact that we are using uniform random numbers
	else
!we use gaussian random numbers
    gkT = 2.0_prec*friction*Boltz*Temp0/dv_mod
	end if
	randva(:)= sqrt (gkT/winv(:))
end if
end if !nodeid .eq. 0

!***********************************************************************
!	begin MAIN DYNAMICS LOOP (Verlet leap-frog algorithm)
!***********************************************************************

! No loop (only calc. energies) if compiling with the DUM flag
#ifndef DUM
do istep = 0, nsteps-1
#endif
!change volume

        if ( mod(istep, NBcycle) .eq. 0 ) then


                ! every NBcycle steps:

				!Put molecules back in box for nice visualisation, needs to be here to prevent problems with LRF
				!Update cgp_centers for LRF
				!only call put_back_in_box if using PBC and either solute or solvent should be put back in box
				if( use_PBC .and. (put_solute_back_in_box .or. put_solvent_back_in_box) ) then 
				  call put_back_in_box() 
				end if
!Time estimate removed, can be activated again by passing -DTIME to the compiler
!Paul Bauer October 2014
#ifdef TIME
                if ((nodeid .eq. 0) .and. (istep > 0)) then
                        ! print timing info
                        call centered_heading('Timing', '-')
                        time1 = rtime()
                        time_per_step = 1000*(time1-time0)/NBcycle
                        time_completion = int(time_per_step*(nsteps-istep)/60000)
                        time0 = time1
                        write(*,222) time_per_step, time_completion
222    format('Milliseconds per step (wall-clock): ',f5.2,&
                                ' Estimated completion in',i6,' minutes')
                end if
#endif
                ! update lists of nonbonded interaction pairs
                if (nodeid .eq. 0) then
                   call centered_heading('Nonbonded pair list generation', '-')
                end if
                call make_pair_lists
#if defined(DUMP)
				write(*,332) 'solute-solute', 'solute-water', 'water-water', 'Q-solute', 'Q-water'
				write(*,333) nodeid, 'count', nbpp_pair, nbpw_pair, &
					 &  nbww_pair, nbqp_pair, 3*nqat*nbqw_pair

#if defined(USE_MPI)
				!reduce totnxx, i.e. collect # pairs found by slave nodes
				nbxx(1)=nbpp_pair
				nbxx(2)=nbpw_pair
				nbxx(3)=nbww_pair
				nbxx(4)=nbqp_pair
				nbxx(5)=3*nqat*nbqw_pair

				call MPI_Reduce(nbxx,nbxx_tot,5,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr) 
				if (ierr .ne. 0) call die('run/Reduce')

				if (nodeid .eq. 0) then
				totnbpp=nbxx_tot(1)
				totnbpw=nbxx_tot(2)
				totnbww=nbxx_tot(3)
				totnbqp=nbxx_tot(4)
				totnbqw=nbxx_tot(5)
				write(*,99) 'total', totnbpp,totnbpw,totnbww,totnbqp,totnbqw
				end if
				99  format(a10,1x,5(1x,i12))
#endif
				332	format('node value ',5a13)
				333	format(i4,1x,a5,1x,5(1x,i12))
#endif



		! do MC_volume step here, after NBupdate for consistency with LRFs and such.
		! Implies a limitation in the volume_update interval, i.e. mod(volume_update,nb_update)=0
				if( use_PBC .and. constant_pressure) then
				   if( mod(istep, ivolume_cycle)==0 .and. istep>0 ) then
					  call MC_volume
				   end if
				end if

        end if ! every NBcycle steps




        ! --- start of time step ---



        ! get potential energy and derivatives from FF
        call pot_energy





if(nodeid .eq. 0) then
        if ( mod(istep,iout_cycle) == 0 .and. monitor_group_pairs > 0) then
           call nonbond_monitor  
        end if

        ! off-diagonals
        if ( noffd .gt. 0 ) call offdiag

#ifndef DUM
        ! update velocities from accelerations,
        ! scale velocities & update positions from velocities
#if defined (PROFILING)
start_loop_time1 = rtime()
#endif

! new selection and subfunction for thermostats
!if Berendsen thermostat was chosen, call two times for diff scaling
if( thermostat == BERENDSEN ) then
	call newvel_ber(dv_mod,Tscale_solute,1,nat_solute)
	call newvel_ber(dv_mod,Tscale_solvent,nat_solute+1,natom)
else if (thermostat == LANGEVIN ) then
	call newvel_lan(dv_mod,dv_friction,randva,randfa)
else if (thermostat == NOSEHOOVER ) then
	call nh_prop
	call newvel_nos(dv_mod)
end if
! --- end of thermostat section ---

#if defined (PROFILING)
profile(11)%time = profile(11)%time + rtime() - start_loop_time1
#endif

        ! shake if necessary
        if(shake_constraints > 0) then
                niter=shake(xx, x)
                v(:) = (x(:) - xx(:)) / dt
		if ( thermostat == NOSEHOOVER ) then
			call nh_prop !scaling velocities after applying SHAKE
		end if
        end if

        ! --- end of time step ---
#if defined (PROFILING)
start_loop_time2 = rtime()
#endif

        ! calculate temperature and scaling factor
        call temperature(Tscale_solute,Tscale_solvent,Ekinmax)
#if defined (PROFILING)
profile(12)%time = profile(12)%time + rtime() - start_loop_time2
#endif

end if !if(nodeid .eq. 0)


#if defined(USE_MPI)
call MPI_Bcast(x, nat3, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('MD Bcast x')
#endif

        ! print [intermediate] results (master node only)
        if (nodeid .eq. 0) then
                ! trajectory, energy data, output and backup restart file
                if ( mod(istep,itrj_cycle) == 0 .and. istep > 0) then
                ! write_trj: write x to the trajectory file
                        call write_trj
                end if
                ! energies
                if ( mod(istep, iene_cycle) == 0 .and. istep > 0) then
                ! nrgy_put_ene(unit, e2, OFFD): print 'e2'=EQ and OFFD to unit 'unit'=11
                        call put_ene(11, EQ, OFFD,loc_arrays,nstates)
                end if
                ! end-of-line, then call write_out, which will print a report on E and EQ
                if ( mod(istep,iout_cycle) == 0 ) then
                        call write_out
                end if
                ! backup file of coordinates and velocities
                if ( mod(istep,1000) .eq. 0 ) then
                        call write_xfin
                end if
                if ( abs(Temp-Tlast)/Temp > TEMP_PRINT_THRESHOLD .or. &
                        (mod(istep, itemp_cycle) == 0 .and. istep > 0)) then
                        ! temperatures
                        Tlast = Temp
                        write(*,201) istep, Temp, Tfree
						if (detail_temps) then
							write(*,2020) Tfree_solute, Tfree_solvent
!							write(*,2030) Texcl_solute, Texcl_solvent

						end if

                        ! test for NaN
                        if (Temp.ne.Temp) then 
                            call die('a detected NaN.')
                        end if

                end if


        end if ! print results


end do ! time step
201	 format('Temperature at step',i8,':         T_tot=',f10.1,'         T_free=',f10.1)
2020 format('                             T_free_solute=',f10.1,' T_free_solvent=',f10.1)
2030 format('                             T_excl_solute=',f10.1,' T_excl_solvent=',f10.1)


!***********************************************************************
!	end MAIN DYNAMICS LOOP
!***********************************************************************

! end of Qdum exclusion
#else
end if !(nodeid .eq. 0) from far above
#endif

        ! write final trajectory image when istep = nsteps
#ifndef DUM
if (nodeid .eq. 0) then
    if ( mod(istep,itrj_cycle) == 0) call write_trj
end if
#endif

        ! write output for final step and final coords
call make_pair_lists
call pot_energy
if (nodeid .eq. 0) then
        write(*,*)
        call write_out
        call write_xfin
end if



if (nodeid .eq. 0) then
!Deallocate temperature control arrays
        if (allocated(randfa)) deallocate(randfa)
        if (allocated(randva)) deallocate(randva)

	time1 = rtime()
	write (*,202) time1 - startloop
	202 format('Total time of main loop:                     ', f15.1,'(s)')
end if
#if defined(PROFILING)
!Print more profiling info

#if defined(USE_MPI)
do i=1,num_profiling_times
	node_times(i) = profile(i)%time
end do
call MPI_GATHER(node_times,num_profiling_times,MPI_REAL8,all_node_times,num_profiling_times,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
if (ierr .ne. 0) call die('md_run/MPI_GATHER profiling times')

if (nodeid .eq. 0) then
write (*,210,advance='no')
do j=0,numnodes-1
	write (*,209,advance='no') j
end do
write(*,*)

do i=1,num_profiling_times
	write (*,207,advance='no') profile(i)%name
	do j=0,numnodes-1
		write (*,208,advance='no') all_node_times(i+j*num_profiling_times)
	end do
	write (*,*) ' (s)'
end do

207	format('Total time of ',A25,T40,': ')
208	format(f10.1,' ')
209	format(I11)
210 format(T30,'node:     ')

write (*,*)

end if

#else

do i=1,num_profiling_times
	write (*,207) profile(i)%name,profile(i)%time
end do
207	format('Total time of ',A25,T40,': ',f15.1,' (s).')
#endif
#endif


end subroutine md_run

!-----------------------------------------------------------------------
subroutine nbpp_count(npp, nppcgp)
! arguments
integer						:: npp
integer						:: nppcgp(:)

! local variables
integer						:: i,j,ig,jg,ia,ja,i3,j3,nl
real(kind=prec)						:: rcut2,r2
integer						:: LJ_code

real(kind=prec)						:: dx, dy, dz
! This routine counts non-bonded solute-solute atom pairs 
! excluding any Q-atoms.

! uses the global variables:
!  Rcpp, ncgp, cgp, excl, x, cgpatom, iqatom, ljcod, crg, 
!  iaclib, max_nbr_range, listex, nexlong, listexlong


npp = 0
rcut2 = Rcpp*Rcpp

igloop: do ig = 1, ncgp_solute
nppcgp(ig) = 0

ia = cgp(ig)%iswitch
 
! skip if excluded group
if ( .not. use_PBC .and. excl(ia) ) cycle igloop

i3 = 3*ia-3

jgloop:	do jg = 1, ncgp_solute
  ja = cgp(jg)%iswitch
 
   ! skip if excluded group
  if ( .not. use_PBC .and. excl(ja) ) cycle jgloop

  ! count each charge group pair once only
  if ( ((ig .gt. jg) .and. (mod(ig+jg,2) .eq. 0)) .or. &
           ((ig .lt. jg) .and. (mod(ig+jg,2) .eq. 1)) ) &
           cycle jgloop

  j3 = 3*ja-3

  !******PWadded if-statement 2001-10-01

  if( .not. use_PBC ) then
  r2 = ( x(i3+1) -x(j3+1) )**2 &
        +( x(i3+2) -x(j3+2) )**2 &
        +( x(i3+3) -x(j3+3) )**2
  else
        dx = x(i3+1) -x(j3+1)
        dy = x(i3+2) -x(j3+2)
        dz = x(i3+3) -x(j3+3)
        dx = dx - boxlength(1)*nint(dx*inv_boxl(1))
        dy = dy - boxlength(2)*nint(dx*inv_boxl(2))
        dz = dz - boxlength(3)*nint(dx*inv_boxl(3))
        r2 = dx**2 + dy**2 + dz**2
  end if


  ! skip if outside cutoff
  if ( r2 .gt. rcut2 ) cycle jgloop

ialoop: do ia = cgp(ig)%first, cgp(ig)%last
        i = cgpatom(ia)

        ! skip if q-atom
        if ( iqatom(i)/=0 ) cycle ialoop

jaloop:	do ja = cgp(jg)%first, cgp(jg)%last
          j = cgpatom(ja)

          ! skip if q-atom
          if ( iqatom(j)/=0 ) cycle jaloop

          if ( ig .eq. jg .and. i .ge. j ) cycle jaloop

          ! skip if all interactions zero
          LJ_code = ljcod(iac(i),iac(j))
          if((crg(i) * crg(j) == zero) &
                .and. &
                (iaclib(iac(i))%avdw(LJ_code)*iaclib(iac(j))%avdw(LJ_code) == zero) &
                .and. &
                (iaclib(iac(i))%bvdw(LJ_code)*iaclib(iac(j))%bvdw(LJ_code) == zero)) &
                cycle jaloop

          ! check bonded exclusions and 1-4 nbors
          if ( abs(j-i) .le. max_nbr_range ) then
                if ( i .lt. j ) then
                  if ( listex(j-i,i) ) cycle jaloop
                else
                  if ( listex(i-j,j) ) cycle jaloop
                end if
          else
                do nl = 1, nexlong
                  if ( (listexlong(1,nl) .eq. i .and. &
                        listexlong(2,nl) .eq. j      ) .or. &
                        (listexlong(1,nl) .eq. j .and. &
                        listexlong(2,nl) .eq. i      ) ) cycle jaloop
                end do
          end if

          ! passed all tests -- count the pair
          npp = npp + 1
          nppcgp(ig) = nppcgp(ig) + 1

        end do jaloop
  end do ialoop
end do jgloop
end do igloop

end subroutine nbpp_count

!-----------------------------------------------------------------------


subroutine nbpplis2
! local variables
integer						:: i,j,ig,jg,ia,ja,i3,j3,nl,inside,jgr
real(kind=prec)						:: rcut2,r2

integer						::iagrid, igrid, jgrid, kgrid, gridnum

! for spherical boundary 
!	This routine makes a list of non-bonded solute-solute atom pairs 
!	excluding any Q-atoms.
! now rewritten to take into account that we already know all the stuff from the
! precomputation

#ifdef _OPENMP
integer :: quotient, remainder
#endif
#if defined (PROFILING)
real(kind=prec)                                         :: start_loop_time
start_loop_time = rtime()
#endif
!$omp single 
nbpp_pair = 0
!$omp end single 
!$omp barrier
rcut2 = Rcpp*Rcpp

! new stuff because we are now having a grid 
! check first if the atoms are in one the adjecent grids
! if not, cycle them
! grid information taken from new costum structure, updated every nb steps
! after getting the first group, we cycle through the grids to find those that interact
! and only search further in them

! interaction matrix is saved in grid_pp_int array 

#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (calculation_assignment%pp%end - calculation_assignment%pp%start + 1)/threads_num
remainder = MOD(calculation_assignment%pp%end - calculation_assignment%pp%start + 1, threads_num)
mp_start = thread_id * quotient + calculation_assignment%pp%start + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif

igloop: do ig = mp_start,mp_end
#else
igloop:  do ig = calculation_assignment%pp%start, calculation_assignment%pp%end
#endif

ia = cgp(ig)%iswitch
if ( excl(ia) ) cycle igloop
#ifdef USE_GRID
iagrid = pp_igrid(ig)
!start at the first grid
gridnum = 0
igridloop:  do igrid = 1, pp_ndim
jgridloop:   do jgrid = 1, pp_ndim
kgridloop:    do kgrid = 1, pp_ndim
		gridnum = gridnum + 1
		if (.not. (grid_pp_int(iagrid,igrid,jgrid,kgrid))) cycle kgridloop
! else, we have interacting grids
! now search for the charge groups in there
! and cycle over all of them
! save the real number in jgr
jgloop: do jg = 1, grid_pp_ngrp(gridnum)
                jgr = grid_pp_grp(gridnum,jg)
#else
jgloop: do jgr = 1, ncgp_solute
#endif
  ! count each charge group pair once only
                if ( ((ig .gt. jgr) .and. (mod(ig+jgr,2) .eq. 0)) .or. &
                        ((ig .lt. jgr) .and. (mod(ig+jgr,2) .eq. 1)) ) &
                        cycle jgloop
                ja = cgp(jgr)%iswitch
                if ( excl(ja) ) cycle jgloop

!	      --- outside cutoff ? ---
                inside = 0
                ia = cgp(ig)%first
                do while ((ia .le. cgp(ig)%last) .and. (inside .eq. 0))
                        i = cgpatom(ia)
                        i3 = 3*i-3
                        ja = cgp(jgr)%first
                        do while ((ja .le. cgp(jgr)%last) .and. (inside .eq. 0))
                                j = cgpatom(ja)
                                j3 = 3*j-3

                                r2 = ( x(i3+1) -x(j3+1) )**2 &
                                        +( x(i3+2) -x(j3+2) )**2 &
                                        +( x(i3+3) -x(j3+3) )**2

                                if ( r2 .le. rcut2 ) then
! one atom pair is within cutoff: set inside
                                        inside = 1
                                end if
                                ja = ja + 1
                        end do
                        ia = ia + 1
                end do
                if (inside .eq. 0) cycle jgloop
ialoop:         do ia = cgp(ig)%first, cgp(ig)%last
                i = cgpatom(ia)
!	         --- q-atom ? ---
                if ( iqatom(i).ne.0 ) cycle ialoop

jaloop:                 do ja = cgp(jgr)%first, cgp(jgr)%last
                                j = cgpatom(ja)
!                --- q-atom ? ---
                                if ( iqatom(j).ne.0 ) cycle jaloop
! make sure each pair is only counted once
                                if ( ig .eq. jgr .and. i .ge. j ) cycle jaloop
! do they interact?
                                if(.not.pp_precomp(i,j)%set) cycle jaloop
! if out of space then make more space
!$omp critical
                                if (nbpp_pair .eq. calculation_assignment%pp%max) call reallocate_nonbondlist_pp
                                nbpp_pair = nbpp_pair + 1
                                nbpp(nbpp_pair)%i = i
                                nbpp(nbpp_pair)%j = j
                                nbpp(nbpp_pair)%vdWA = pp_precomp(i,j)%vdWA
                                nbpp(nbpp_pair)%vdWB = pp_precomp(i,j)%vdWB
                                nbpp(nbpp_pair)%elec = pp_precomp(i,j)%elec
!$omp end critical

                        end do jaloop
                end do ialoop
        end do jgloop
#ifdef USE_GRID
	end do kgridloop
	end do jgridloop
	end do igridloop
#endif
end do igloop
#if defined (PROFILING) 
profile(3)%time = profile(3)%time + rtime() - start_loop_time
#endif 

end subroutine nbpplis2


!--------------------------------------------------------------------

subroutine nbpplis2_box
! local variables
integer						:: i,j,ig,jg,ia,ja,i3,j3,nl,inside,jgr
real(kind=prec)					:: rcut2,r2
real(kind=prec)					:: dx, dy, dz

integer                                         ::iagrid, igrid, jgrid, kgrid, gridnum
#ifdef _OPENMP
integer :: quotient, remainder
#endif
#if defined (PROFILING)
real(kind=prec)                                         :: start_loop_time
start_loop_time = rtime()
#endif
 
  ! for periodic boundary conditions
  !	This routine makes a list of non-bonded solute-solute atom pairs 
  !	excluding any Q-atoms.
!$omp single
  nbpp_pair = 0
  nbpp_cgp_pair = 0
!$omp end single
!$omp barrier
  rcut2 = Rcpp*Rcpp
#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (calculation_assignment%pp%end - calculation_assignment%pp%start + 1)/threads_num
remainder = MOD(calculation_assignment%pp%end - calculation_assignment%pp%start + 1, threads_num)
mp_start = thread_id * quotient + calculation_assignment%pp%start + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif

igloop: do ig = mp_start,mp_end
#else
igloop: do ig = calculation_assignment%pp%start, calculation_assignment%pp%end
#endif
#ifdef USE_GRID
        iagrid = pp_igrid(ig)
        gridnum = 0
igridloop: do igrid = 1, pp_ndim
jgridloop:  do jgrid = 1, pp_ndim
kgridloop:   do kgrid = 1, pp_ndim
		gridnum = gridnum + 1
		if (.not. (grid_pp_int(iagrid,igrid,jgrid,kgrid))) cycle kgridloop
! else, we have interacting grids
! now search for the charge groups in there
! and cycle over all of them
! save the real number in jgr
jgloop:         do jg = 1, grid_pp_ngrp(gridnum)
                        jgr = grid_pp_grp(gridnum,jg)
#else
jgloop:         do jgr = 1, ncgp_solute
#endif
! count each charge group pair once only
                        if ( ((ig .gt. jgr) .and. (mod(ig+jgr,2) .eq. 0)) .or. &
                                ((ig .lt. jgr) .and. (mod(ig+jgr,2) .eq. 1)) ) &
                                cycle jgloop

!      --- outside cutoff ? ---
                        inside = 0
                        ia = cgp(ig)%first
                        do while ((ia .le. cgp(ig)%last) .and. (inside .eq. 0))
                                i = cgpatom(ia)
                                i3 = 3*i-3
                                ja = cgp(jgr)%first
                                do while ((ja .le. cgp(jgr)%last) .and. (inside .eq. 0))
                                        j = cgpatom(ja)
                                        j3 = 3*j-3
                                        dx = x(i3+1) - x(j3+1)
                                        dy = x(i3+2) - x(j3+2)
                                        dz = x(i3+3) - x(j3+3)
                                        dx = dx - boxlength(1)*nint( dx*inv_boxl(1) )
                                        dy = dy - boxlength(2)*nint( dy*inv_boxl(2) )
                                        dz = dz - boxlength(3)*nint( dz*inv_boxl(3) )
                                        r2 = dx**2 + dy**2 + dz**2

                                        if ( r2 .le. rcut2 ) then
! one atom pair is within cutoff: set inside
                                                inside = 1
!$omp critical
                                                if (nbpp_cgp_pair .eq. size(nbpp_cgp, 1) )  call reallocate_nbpp_cgp
                                                nbpp_cgp_pair = nbpp_cgp_pair + 1
                                                nbpp_cgp(nbpp_cgp_pair)%i = i
                                                nbpp_cgp(nbpp_cgp_pair)%j = j
!$omp end critical
                                        end if
                                        ja = ja + 1
                                end do
                                ia = ia + 1
                        end do
                        if (inside .eq. 0) cycle jgloop

ialoop:                 do ia = cgp(ig)%first, cgp(ig)%last
                                i = cgpatom(ia)
!	         --- q-atom ? ---
                                if ( iqatom(i).ne.0 ) cycle ialoop
jaloop:                         do ja = cgp(jgr)%first, cgp(jgr)%last
                                        j = cgpatom(ja)
!                --- q-atom ? ---
                                        if ( iqatom(j).ne.0 ) cycle jaloop
! make sure each pair is only counted once
                                        if ( ig .eq. jgr .and. i .ge. j ) cycle jaloop
! do they interact?
                                        if(.not.pp_precomp(i,j)%set) cycle jaloop
! if out of space then make more space
!$omp critical
                                        if (nbpp_pair .eq. calculation_assignment%pp%max) call reallocate_nonbondlist_pp
                                        nbpp_pair = nbpp_pair + 1
                                        nbpp(nbpp_pair)%i = i
                                        nbpp(nbpp_pair)%j = j
                                        nbpp(nbpp_pair)%vdWA = pp_precomp(i,j)%vdWA
                                        nbpp(nbpp_pair)%vdWB = pp_precomp(i,j)%vdWB
                                        nbpp(nbpp_pair)%elec = pp_precomp(i,j)%elec
                                        nbpp(nbpp_pair)%cgp_pair = nbpp_cgp_pair
!$omp end critical
                                end do jaloop
                        end do ialoop
                end do jgloop
#ifdef USE_GRID
	end do kgridloop
	end do jgridloop
	end do igridloop
#endif
end do igloop
#if defined (PROFILING) 
profile(3)%time = profile(3)%time + rtime() - start_loop_time
#endif 

end subroutine nbpplis2_box
!-----------------------------------------------------------------------------
subroutine nbpplis2_box_lrf
  ! local variables
  integer						:: i,j,ig,jg,ia,ja,i3,j3,nl,inside,jgr
  real(kind=prec)						:: rcut2,r2
  real(kind=prec)						:: dx, dy, dz
 
  real(kind=prec)						::RcLRF2
  integer						::inside_LRF, is3
  integer						::iagrid, igrid, jgrid, kgrid, gridnum
  ! for periodic boundary conditions
  !	This routine makes a list of non-bonded solute-solute atom pairs 
  !	excluding any Q-atoms.
#ifdef _OPENMP
integer :: quotient, remainder
#endif
#if defined (PROFILING)
real(kind=prec)                                         :: start_loop_time
start_loop_time = rtime()
#endif
!$omp single
  nbpp_pair = 0
  nbpp_cgp_pair = 0
!$omp end single
!$omp barrier
  rcut2 = Rcpp*Rcpp
  RcLRF2 = RcLRF*RcLRF
#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (calculation_assignment%pp%end - calculation_assignment%pp%start + 1)/threads_num
remainder = MOD(calculation_assignment%pp%end - calculation_assignment%pp%start + 1, threads_num)
mp_start = thread_id * quotient + calculation_assignment%pp%start + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif

igloop: do ig = mp_start,mp_end
#else
igloop: do ig = calculation_assignment%pp%start, calculation_assignment%pp%end
#endif
#ifdef USE_GRID
         iagrid = pp_igrid(ig)
	  gridnum = 0
igridloop:  do igrid = 1, pp_ndim
jgridloop:   do jgrid = 1, pp_ndim
kgridloop:    do kgrid = 1, pp_ndim
		gridnum = gridnum + 1
		if (.not. (grid_pp_int(iagrid,igrid,jgrid,kgrid))) then
ggloop:		 do jg = 1, grid_pp_ngrp(gridnum)
		  jgr = grid_pp_grp(gridnum,jg)
		   if ( ((ig .gt. jgr) .and. (mod(ig+jgr,2) .eq. 0)) .or. &
                   ((ig .lt. jgr) .and. (mod(ig+jgr,2) .eq. 1)) .or.&
                   (ig.eq.jgr)) &
                   cycle ggloop
		    call lrf_update(ig,jgr)
		    call lrf_update(jgr,ig)
		   end do ggloop
		 else 
! we have real interacting grids , do the rest of the calculations
jgloop: do jg = 1, grid_pp_ngrp(gridnum)
                jgr = grid_pp_grp(gridnum,jg)
#else
jgloop: do jgr = 1, ncgp_solute
#endif
! count each charge group pair once only
                if ( ((ig .gt. jgr) .and. (mod(ig+jgr,2) .eq. 0)) .or. &
                        ((ig .lt. jgr) .and. (mod(ig+jgr,2) .eq. 1)) ) &
                        cycle jgloop
!      --- outside cutoff ? ---
                inside = 0
                inside_LRF = 0
                ia = cgp(ig)%first
                do while ((ia .le. cgp(ig)%last) .and. (inside .eq. 0))
                        i = cgpatom(ia)
                        i3 = 3*i-3
                        ja = cgp(jgr)%first
                        do while ((ja .le. cgp(jgr)%last) .and. (inside .eq. 0))
                                j = cgpatom(ja)
                                j3 = 3*j-3
                                dx = x(i3+1) - x(j3+1)
                                dy = x(i3+2) - x(j3+2)
                                dz = x(i3+3) - x(j3+3)
                                dx = dx - boxlength(1)*nint( dx*inv_boxl(1) )
                                dy = dy - boxlength(2)*nint( dy*inv_boxl(2) )
                                dz = dz - boxlength(3)*nint( dz*inv_boxl(3) )
                                r2 = dx**2 + dy**2 + dz**2
                                if ( r2 .le. rcut2 ) then
! one atom pair is within cutoff: set inside
                                        inside = 1
!$omp critical
                                        if (nbpp_cgp_pair .eq. size(nbpp_cgp, 1) )  call reallocate_nbpp_cgp
                                        nbpp_cgp_pair = nbpp_cgp_pair + 1
                                        nbpp_cgp(nbpp_cgp_pair)%i = i
                                        nbpp_cgp(nbpp_cgp_pair)%j = j
!$omp end critical
                                elseif ((r2 .le. RcLRF2) .or. (RcLRF .eq. -one)) then
                                        inside_LRF = 1
                                end if
                                ja = ja + 1
                        end do ! ia .le. cgp(ig)%last
                        ia = ia + 1
                end do ! ja .le. cgp(jgr)%last
		
                if (inside .eq. 1) then
ialoop:                 do ia = cgp(ig)%first, cgp(ig)%last
                                i = cgpatom(ia)
!	         --- q-atom ? ---
                                if ( iqatom(i).ne.0 ) cycle ialoop
jaloop:                         do ja = cgp(jgr)%first, cgp(jgr)%last
                                        j = cgpatom(ja)
!                --- q-atom ? ---
                                        if ( iqatom(j).ne.0 ) cycle jaloop
! make sure each pair is only counted once
                                        if ( ig .eq. jgr .and. i .ge. j ) cycle jaloop
! do they interact?
                                        if(.not.pp_precomp(i,j)%set) cycle jaloop
! if out of space then make more space
!$omp critical
                                        if (nbpp_pair .eq. calculation_assignment%pp%max) call reallocate_nonbondlist_pp
                                        nbpp_pair = nbpp_pair + 1
                                        nbpp(nbpp_pair)%i = i
                                        nbpp(nbpp_pair)%j = j
                                        nbpp(nbpp_pair)%vdWA = pp_precomp(i,j)%vdWA
                                        nbpp(nbpp_pair)%vdWB = pp_precomp(i,j)%vdWB
                                        nbpp(nbpp_pair)%elec = pp_precomp(i,j)%elec
                                        nbpp(nbpp_pair)%cgp_pair = nbpp_cgp_pair
!$omp end critical
                                end do jaloop
                        end do ialoop
                elseif((inside_LRF .eq. 1) .and. (inside .eq. 0)) then
! outside pp-cutoff but inside LRF cut-off use LRF
!ig : jgr calculation
                        call lrf_update(ig,jgr)
!jgr : ig calculations
                        call lrf_update(jgr,ig)
                end if ! outside cutoff
        end do jgloop
#ifdef USE_GRID
	end if ! interacting boxes
	end do kgridloop
	end do jgridloop
	end do igridloop
#endif
end do igloop
#if defined (PROFILING) 
profile(3)%time = profile(3)%time + rtime() - start_loop_time
#endif 

end subroutine nbpplis2_box_lrf
!-----------------------------------------------------------------------
subroutine nbpplis2_lrf
! local variables
integer						:: i,j,ig,jg,ia,ja,i3,j3,nl,is,jgr
logical						::	inside
real(kind=prec)						:: rcut2,r2
real(kind=prec)						::	RcLRF2

integer							::iagrid, igrid, jgrid, kgrid, gridnum
!	This routine makes a list of non-bonded solute-solute atom pairs 
!	excluding any Q-atoms.
#ifdef _OPENMP
integer :: quotient, remainder
#endif
#if defined (PROFILING)
real(kind=prec)                                         :: start_loop_time
start_loop_time = rtime()
#endif
!$omp single
nbpp_pair = 0
!$omp end single
!$omp barrier
rcut2 = Rcpp*Rcpp
RcLRF2 = RcLRF*RcLRF

#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (calculation_assignment%pp%end - calculation_assignment%pp%start + 1)/threads_num
remainder = MOD(calculation_assignment%pp%end - calculation_assignment%pp%start + 1, threads_num)
mp_start = thread_id * quotient + calculation_assignment%pp%start + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif

igloop: do ig = mp_start,mp_end
#else
igloop: do ig = calculation_assignment%pp%start, calculation_assignment%pp%end
#endif

        ! skip if excluded group
        is = cgp(ig)%iswitch
        if ( excl(is) ) cycle igloop
#ifdef USE_GRID
        gridnum = 0
	iagrid = pp_igrid(ig)
igridloop:  do igrid = 1, pp_ndim
jgridloop:   do jgrid = 1, pp_ndim
kgridloop:    do kgrid = 1, pp_ndim
		gridnum = gridnum + 1
                if (.not. (grid_pp_int(iagrid,igrid,jgrid,kgrid))) then
ggloop:          do jg = 1, grid_pp_ngrp(gridnum)
                  jgr = grid_pp_grp(gridnum,jg)
		  ja = cgp(jgr)%iswitch
                  if ( excl(ja) ) cycle ggloop
                   if ( ((ig .gt. jgr) .and. (mod(ig+jgr,2) .eq. 0)) .or. &
                   ((ig .lt. jgr) .and. (mod(ig+jgr,2) .eq. 1)) .or. &
                   (ig.eq.jgr)) &
                   cycle ggloop
                    call lrf_update(ig,jgr)
                    call lrf_update(jgr,ig)
                   end do ggloop
                 else
! we have real interacting grids , do the rest of the calculations
jgloop: do jg = 1, grid_pp_ngrp(gridnum)
                jgr = grid_pp_grp(gridnum,jg)
#else
jgloop: do jgr = 1, ncgp_solute
#endif
! count each charge group pair once only
                if( ((ig .gt. jgr) .and. (mod(ig+jgr,2) .eq. 0)) .or. &
                        ((ig .lt. jgr) .and. (mod(ig+jgr,2) .eq. 1)) ) &
                        cycle jgloop

!	      --- excluded group ? ---
                ja = cgp(jgr)%iswitch
                if ( excl(ja) ) cycle jgloop

!	      --- outside cutoff ? ---
                inside = .false.
                ia = cgp(ig)%first
pairloop:	do while( (ia .le. cgp(ig)%last) .and. (inside.eqv..false.) )
                        i = cgpatom(ia)
                        i3 = 3*i-3
			ja = cgp(jgr)%first
                        do while (( ja.le.cgp(jgr)%last) .and. (inside.eqv..false.))
                                j = cgpatom(ja)
                                j3 = 3*j-3

                                r2 = ( x(i3+1) -x(j3+1) )**2 &
                                   +( x(i3+2) -x(j3+2) )**2 &
                                   +( x(i3+3) -x(j3+3) )**2

                                if ( r2 <= rcut2 ) then
                                        ! one atom pair is within cutoff: set inside
                                        inside = .true.
                                end if
                        ja = ja + 1
                        end do
                        ia = ia + 1
                end do pairloop
!	      --- inside cutoff ? ---
                if (inside) then
ialoop:                 do ia = cgp(ig)%first, cgp(ig)%last
                                i = cgpatom(ia)

                                !	             --- q-atom ? ---
                                if ( iqatom(i).ne.0 ) cycle ialoop
jaloop:				do ja = cgp(jgr)%first, cgp(jgr)%last
                                        j = cgpatom(ja)
!                --- q-atom ? ---
                                        if ( iqatom(j).ne.0 ) cycle jaloop
! make sure each pair is only counted once
                                        if ( ig .eq. jgr .and. i .ge. j ) cycle jaloop
! do they interact?
                                        if(.not.pp_precomp(i,j)%set) cycle jaloop

                                        ! if out of space then make more space
!$omp critical
                                        if (nbpp_pair == calculation_assignment%pp%max) call reallocate_nonbondlist_pp

                                        nbpp_pair = nbpp_pair + 1
                                        nbpp(nbpp_pair)%i = i
                                        nbpp(nbpp_pair)%j = j
                                        nbpp(nbpp_pair)%vdWA = pp_precomp(i,j)%vdWA
                                        nbpp(nbpp_pair)%vdWB = pp_precomp(i,j)%vdWB
                                        nbpp(nbpp_pair)%elec = pp_precomp(i,j)%elec
!$omp end critical
                                end do jaloop
                        end do ialoop
                elseif(r2 .le. RcLRF2) then
                        ! outside pp-cutoff but inside LRF cut-off: use LRF
			call lrf_update(ig,jgr)
			call lrf_update(jgr,ig)
                end if
        end do jgloop
#ifdef USE_GRID
	end if ! interacting boxes
	end do kgridloop
	end do jgridloop
	end do igridloop
#endif
end do igloop
#if defined (PROFILING) 
profile(3)%time = profile(3)%time + rtime() - start_loop_time
#endif 

end subroutine nbpplis2_lrf

!-----------------------------------------------------------------------

subroutine nbpplist
! local variables
integer						:: i,j,ig,jg,ia,ja,i3,j3,nl,jgr
real(kind=prec)						:: rcut2,r2
integer						:: LJ_code
integer                                                 ::iagrid, igrid, jgrid, kgrid, gridnum

! For use with spherical boundary   
!	This routine makes a list of non-bonded solute-solute atom pairs 
!	excluding any Q-atoms.

! uses the global variables:
!  Rcpp, ncgp, cgp, excl, x, cgpatom, iqatom, ljcod, crg, iaclib, max_nbr_range, listex
!  nexlong, listexlong, calculation_assignment%pp%max, alloc_status, list14, n14long, list14long

! reset #pairs
#ifdef _OPENMP
integer :: quotient, remainder
#endif
#if defined (PROFILING)
real(kind=prec)                                         :: start_loop_time
start_loop_time = rtime()
#endif
!$omp single
nbpp_pair = 0
!$omp end single
!$omp barrier
rcut2 = Rcpp*Rcpp
#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (calculation_assignment%pp%end - calculation_assignment%pp%start + 1)/threads_num
remainder = MOD(calculation_assignment%pp%end - calculation_assignment%pp%start + 1, threads_num)
mp_start = thread_id * quotient + calculation_assignment%pp%start + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif

igloop: do ig = mp_start,mp_end
#else
igloop: do ig = calculation_assignment%pp%start, calculation_assignment%pp%end
#endif
! for every assigned charge group:

! skip if excluded group
        ia = cgp(ig)%iswitch
        if ( excl(ia) ) cycle igloop
        i3 = 3*ia-3
#ifdef USE_GRID
iagrid = pp_igrid(ig)
gridnum = 0

igridloop:	do igrid = 1, pp_ndim
jgridloop:	 do jgrid = 1, pp_ndim
kgridloop:	  do kgrid = 1, pp_ndim
			gridnum = gridnum + 1
			if (.not. (grid_pp_int(iagrid,igrid,jgrid,kgrid))) cycle kgridloop


jgloop: do jg = 1, grid_pp_ngrp(gridnum) 
                jgr = grid_pp_grp(gridnum,jg)
#else
jgloop: do jgr = 1, ncgp_solute
#endif
! for every charge group:
! count each charge group pair once only
                if ( ((ig .gt. jgr) .and. (mod(ig+jgr,2) .eq. 0)) .or. &
                        ((ig .lt. jgr) .and. (mod(ig+jgr,2) .eq. 1)) ) &
                        cycle jgloop
! skip if excluded group
                ja = cgp(jgr)%iswitch
                if ( excl(ja) ) cycle jgloop
                j3 = 3*ja-3
                r2 = ( x(i3+1) - x(j3+1) )**2 &
                        +( x(i3+2) - x(j3+2) )**2 &
                        +( x(i3+3) - x(j3+3) )**2

! skip if outside cutoff
                if ( r2 .gt. rcut2 ) cycle jgloop
ialoop:         do ia = cgp(ig)%first, cgp(ig)%last
! for every atom in the charge group (of the outermost loop):
                        i = cgpatom(ia)
! skip if q-atom
                        if ( iqatom(i).ne.0 ) cycle ialoop
jaloop:                 do ja = cgp(jgr)%first, cgp(jgr)%last
! for every atom in the charge group (innermost loop)
                                j = cgpatom(ja)
!                --- q-atom ? ---
                                if ( iqatom(j).ne.0 ) cycle jaloop
! make sure each pair is only counted once
                                if ( ig .eq. jgr .and. i .ge. j ) cycle jaloop
! do they interact?
                                if(.not.pp_precomp(i,j)%set) cycle jaloop
! if out of space then make more space
!$omp critical
                                if (nbpp_pair .eq. calculation_assignment%pp%max) call reallocate_nonbondlist_pp
! all tests passed, add the pair
                                nbpp_pair = nbpp_pair + 1
                                nbpp(nbpp_pair)%i = i
                                nbpp(nbpp_pair)%j = j 
                                nbpp(nbpp_pair)%vdWA = pp_precomp(i,j)%vdWA
                                nbpp(nbpp_pair)%vdWB = pp_precomp(i,j)%vdWB
                                nbpp(nbpp_pair)%elec = pp_precomp(i,j)%elec
!$omp end critical

                        end do jaloop
                end do ialoop
        end do jgloop
#ifdef USE_GRID
end do kgridloop
end do jgridloop
end do igridloop
#endif
end do igloop
#if defined (PROFILING) 
profile(3)%time = profile(3)%time + rtime() - start_loop_time
#endif 

end subroutine nbpplist

!-----------------------------------------------------------------------

subroutine nbpplist_box
  ! local variables
  integer						:: i,j,ig,jg,ia,ja,i3,j3,nl,ig_sw, jg_sw,jgr
  real(kind=prec)						:: rcut2,r2
  integer						:: LJ_code
  real(kind=prec)						:: dx, dy, dz
  integer                                                 ::iagrid, igrid, jgrid, kgrid, gridnum
  
  ! For use with periodic boundary conditions
  !	This routine makes a list of non-bonded solute-solute atom pairs 
  !	excluding any Q-atoms.

  ! uses the global variables:
  !  Rcpp, ncgp, cgp, excl, x, cgpatom, iqatom, ljcod, crg, iaclib, max_nbr_range, listex
  !  nexlong, listexlong, calculation_assignment%pp%max, alloc_status, list14, n14long, list14long

  ! reset #pairs
#ifdef _OPENMP
integer :: quotient, remainder
#endif
#if defined (PROFILING)
real(kind=prec)                                         :: start_loop_time
start_loop_time = rtime()
#endif
!$omp single
  nbpp_pair = 0 !atom pairs
  nbpp_cgp_pair = 0 !chargegroup pairs
!$omp end single
!$omp barrier
  rcut2 = Rcpp*Rcpp

#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (calculation_assignment%pp%end - calculation_assignment%pp%start + 1)/threads_num
remainder = MOD(calculation_assignment%pp%end - calculation_assignment%pp%start + 1, threads_num)
mp_start = thread_id * quotient + calculation_assignment%pp%start + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif

igloop: do ig = mp_start,mp_end
#else
igloop: do ig = calculation_assignment%pp%start, calculation_assignment%pp%end
#endif
	! for every assigned charge group:
	ig_sw = cgp(ig)%iswitch !switching atom in charge group ig
	i3 = 3*ig_sw-3
#ifdef USE_GRID
        iagrid = pp_igrid(ig)
	gridnum = 0
igridloop: do igrid = 1, pp_ndim
jgridloop:  do jgrid = 1, pp_ndim
kgridloop:   do kgrid = 1, pp_ndim
		gridnum = gridnum + 1
		if (.not. (grid_pp_int(iagrid,igrid,jgrid,kgrid))) cycle kgridloop

jgloop: do jg = 1, grid_pp_ngrp(gridnum)
		jgr = grid_pp_grp(gridnum,jg)
#else
jgloop: do jgr = 1, ncgp_solute
#endif
	  ! for every charge group:

	  ! count each charge group pair once only
	  if ( ((ig .gt. jgr) .and. (mod(ig+jgr,2) .eq. 0)) .or. &
		   ((ig .lt. jgr) .and. (mod(ig+jgr,2) .eq. 1)) ) &
		   cycle jgloop

	  jg_sw = cgp(jgr)%iswitch !switching atom in charge group jg
	  j3 = 3*jg_sw-3
	  dx = x(i3+1) - x(j3+1)
	  dy = x(i3+2) - x(j3+2)
	  dz = x(i3+3) - x(j3+3)
	  dx = dx - boxlength(1)*nint( dx*inv_boxl(1) )
	  dy = dy - boxlength(2)*nint( dy*inv_boxl(2) )
	  dz = dz - boxlength(3)*nint( dz*inv_boxl(3) )
	  r2 = dx**2 + dy**2 + dz**2
	  ! skip if outside cutoff
	  if ( r2 .gt. rcut2 ) cycle jgloop

	  !inside cutoff

	  !check if more memory is needed
!$omp critical
	  if(nbpp_cgp_pair .eq. size(nbpp_cgp, 1)) call reallocate_nbpp_cgp
	  !add the charge group pair
	  nbpp_cgp_pair = nbpp_cgp_pair + 1
	  nbpp_cgp(nbpp_cgp_pair)%i = ig_sw !the switching atoms of the charge groups in the pair
	  nbpp_cgp(nbpp_cgp_pair)%j = jg_sw
!$omp end critical

ialoop:         do ia = cgp(ig)%first, cgp(ig)%last
! for every atom in the charge group ig (of the outermost loop):
                        i = cgpatom(ia)
! skip if q-atom
                        if ( iqatom(i).ne.0 ) cycle ialoop
jaloop:                 do ja = cgp(jgr)%first, cgp(jgr)%last
! for every atom in the charge group jg (innermost loop)
                                j = cgpatom(ja)
!                --- q-atom ? ---
                                if ( iqatom(j).ne.0 ) cycle jaloop
! make sure each pair is only counted once
                                if ( ig .eq. jgr .and. i .ge. j ) cycle jaloop
! do they interact?
                                if(.not.pp_precomp(i,j)%set) cycle jaloop
! if out of space then make more space
!$omp critical
                                if (nbpp_pair .eq. calculation_assignment%pp%max) call reallocate_nonbondlist_pp
! all tests passed, add the pair
                                nbpp_pair = nbpp_pair + 1
                                nbpp(nbpp_pair)%i = i
                                nbpp(nbpp_pair)%j = j 
                                nbpp(nbpp_pair)%vdWA = pp_precomp(i,j)%vdWA
                                nbpp(nbpp_pair)%vdWB = pp_precomp(i,j)%vdWB
                                nbpp(nbpp_pair)%elec = pp_precomp(i,j)%elec
                                nbpp(nbpp_pair)%cgp_pair = nbpp_cgp_pair !which pair of charge groups the atom pair belongs to
!$omp end critical		
                        end do jaloop
                end do ialoop
	end do jgloop
#ifdef USE_GRID
	end do kgridloop
	end do jgridloop
	end do igridloop
#endif
end do igloop
#if defined (PROFILING) 
profile(3)%time = profile(3)%time + rtime() - start_loop_time
#endif 

end subroutine nbpplist_box

!--------------------------------------------------------------------------------------

subroutine nbpplist_lrf
! local variables
integer						:: i,j,ig,jg,ia,ja,i3,j3,nl,is,is3,jgr
real(kind=prec)						:: rcut2,r2
integer						:: LJ_code
real(kind=prec)						::	RcLRF2
integer                                                 ::iagrid, igrid, jgrid, kgrid, gridnum
#ifdef _OPENMP
integer :: quotient, remainder
#endif
#if defined (PROFILING)
real(kind=prec)						:: start_loop_time
start_loop_time = rtime()
#endif


!	This routine makes a list of non-bonded solute-solute atom pairs 
!	excluding any Q-atoms.

! uses the global variables:
!  nbpp_pair, Rcpp, RcLRF, cgp, excl, ncgp, x, cgpatom, iqatom, ljcod, crg, iaclib, max_nbr_range,
!  listex, nexlong, listexlong, nbpp, alloc_status, list14, n14long, list14long, lrf

!$omp single
nbpp_pair = 0
!$omp end single
!$omp barrier
rcut2 = Rcpp*Rcpp
RcLRF2 = RcLRF*RcLRF
#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (calculation_assignment%pp%end - calculation_assignment%pp%start + 1)/threads_num
remainder = MOD(calculation_assignment%pp%end - calculation_assignment%pp%start + 1, threads_num)
mp_start = thread_id * quotient + calculation_assignment%pp%start + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif

igloop: do ig = mp_start,mp_end
#else
igloop: do ig = calculation_assignment%pp%start, calculation_assignment%pp%end
#endif
! for every assigned charge group:

! skip if excluded group
        is = cgp(ig)%iswitch
        if ( excl(is) ) cycle igloop
        is3 = 3*is-3
#ifdef USE_GRID
iagrid = pp_igrid(ig)
gridnum = 0
igridloop:	do igrid = 1, pp_ndim
jgridloop:	 do jgrid = 1, pp_ndim
kgridloop:	  do kgrid = 1, pp_ndim
			gridnum = gridnum + 1
			if (.not. (grid_pp_int(iagrid,igrid,jgrid,kgrid))) then
ggloop:			 do jg = 1, grid_pp_ngrp(gridnum)
			  jgr = grid_pp_grp(gridnum,jg)
                          ja = cgp(jgr)%iswitch
                          if ( excl(ja) ) cycle ggloop
			  if ( ((ig .gt. jgr) .and. (mod(ig+jgr,2) .eq. 0)) .or. &
			  ((ig .lt. jgr) .and. (mod(ig+jgr,2) .eq. 1)) .or. &
                          (ig.eq.jgr)) &
			  cycle ggloop
			  call lrf_update(ig,jgr)
			  call lrf_update(jgr,ig)
			 end do ggloop
			else
jgloop: do jg = 1, grid_pp_ngrp(gridnum)
                jgr = grid_pp_grp(gridnum,jg)
#else
jgloop: do jgr = 1, ncgp_solute
#endif
! count each charge group pair once only
                if ( ((ig .gt. jgr) .and. (mod(ig+jgr,2) .eq. 0)) .or. &
                        ((ig .lt. jgr) .and. (mod(ig+jgr,2) .eq. 1)) ) &
                        cycle jgloop
! skip if excluded group
                ja = cgp(jgr)%iswitch
                if ( excl(ja) ) cycle jgloop
                j3 = 3*ja-3
                r2 = ( x(is3+1) -x(j3+1) )**2 &
                        +( x(is3+2) -x(j3+2) )**2 &
                        +( x(is3+3) -x(j3+3) )**2

! inside cutoff?
                if ( r2 .le. rcut2 ) then
ialoop:                 do ia = cgp(ig)%first, cgp(ig)%last
! skip if q-atom
                                i = cgpatom(ia)
                                if ( iqatom(i).ne.0 ) cycle ialoop
jaloop:                         do ja = cgp(jgr)%first, cgp(jgr)%last
                                        j = cgpatom(ja)
! skip if q-atom
                                        if ( iqatom(j).ne.0 ) cycle jaloop
! make sure each pair is only counted once
                                        if ( ig .eq. jgr .and. i .ge. j ) cycle jaloop
! do they interact?
                                        if(.not.pp_precomp(i,j)%set) cycle jaloop
! if out of space then make more space
!$omp critical
                                        if (nbpp_pair .eq. calculation_assignment%pp%max) call reallocate_nonbondlist_pp
! all tests passed, add the pair
                                        nbpp_pair = nbpp_pair + 1
                                        nbpp(nbpp_pair)%i = i
                                        nbpp(nbpp_pair)%j = j
                                        nbpp(nbpp_pair)%vdWA = pp_precomp(i,j)%vdWA
                                        nbpp(nbpp_pair)%vdWB = pp_precomp(i,j)%vdWB
                                        nbpp(nbpp_pair)%elec = pp_precomp(i,j)%elec
!$omp end critical
                                end do jaloop
                        end do ialoop
                elseif(r2 .le. RcLRF2) then
! outside pp-cutoff but inside LRF cut-off use LRF
                        call lrf_update(ig,jgr)
                        call lrf_update(jgr,ig)
                end if ! LRF cutoff
        end do jgloop
#ifdef USE_GRID
end if ! interaction
end do kgridloop
end do jgridloop
end do igridloop
#endif
end do igloop
#if defined (PROFILING)
profile(3)%time = profile(3)%time + rtime() - start_loop_time
#endif

end subroutine nbpplist_lrf
!----------------LRF version of PW PBC-----------------------
subroutine nbpplist_box_lrf
  ! local variables
  integer						:: i,j,ig,jg,ia,ja,i3,j3,nl,ig_sw, jg_sw, is3, jgr
  real(kind=prec)						:: rcut2,r2
  integer						:: LJ_code
  real(kind=prec)						:: dx, dy, dz
  
  real(kind=prec)						::RcLRF2

  integer                                                 ::iagrid, igrid, jgrid, kgrid, gridnum
#ifdef _OPENMP
integer :: quotient, remainder
#endif
#if defined (PROFILING)
real(kind=prec)                                         :: start_loop_time
start_loop_time = rtime()
#endif

  ! For use with periodic boundary conditions
  !	This routine makes a list of non-bonded solute-solute atom pairs 
  !	excluding any Q-atoms.

  ! uses the global variables:
  !  Rcpp, ncgp, cgp, excl, x, cgpatom, iqatom, ljcod, crg, iaclib, max_nbr_range, listex
  !  nexlong, listexlong, calculation_assignment%pp%max, alloc_status, list14, n14long, list14long

  ! reset #pairs
!$omp single
  nbpp_pair = 0 !atom pairs
  nbpp_cgp_pair = 0 !chargegroup pairs
!$omp end single
!$omp barrier
  rcut2 = Rcpp*Rcpp
  RcLRF2 = RcLRF*RcLRF
#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (calculation_assignment%pp%end - calculation_assignment%pp%start + 1)/threads_num
remainder = MOD(calculation_assignment%pp%end - calculation_assignment%pp%start + 1, threads_num)
mp_start = thread_id * quotient + calculation_assignment%pp%start + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif

igloop: do ig = mp_start,mp_end
#else
igloop: do ig = calculation_assignment%pp%start, calculation_assignment%pp%end
#endif
	! for every assigned charge group:
	ig_sw = cgp(ig)%iswitch !switching atom in charge group ig
	i3 = 3*ig_sw-3
#ifdef USE_GRID
	iagrid = pp_igrid(ig)
	gridnum = 0
igridloop:	do igrid = 1, pp_ndim
jgridloop:	 do jgrid = 1, pp_ndim
kgridloop:	  do kgrid = 1, pp_ndim
			gridnum = gridnum + 1
			if (.not. (grid_pp_int(iagrid,igrid,jgrid,kgrid))) then
ggloop:			 do jg = 1, grid_pp_ngrp(gridnum)
			  jgr = grid_pp_grp(gridnum,jg)
			  if ( ((ig .gt. jgr) .and. (mod(ig+jgr,2) .eq. 0)) .or. &
	                   ((ig .lt. jgr) .and. (mod(ig+jgr,2) .eq. 1)) .or. &
                           (ig.eq.jgr)) &
        	           cycle ggloop
			  call lrf_update(ig,jgr)
			  call lrf_update(jgr,ig)
			 end do ggloop
			 else
jgloop:	do jg = 1, grid_pp_ngrp(gridnum)
	        jgr = grid_pp_grp(gridnum,jg)
#else
jgloop: do jgr = 1, ncgp_solute
#endif
! for every charge group:
! count each charge group pair once only
                if ( ((ig .gt. jgr) .and. (mod(ig+jgr,2) .eq. 0)) .or. &
                        ((ig .lt. jgr) .and. (mod(ig+jgr,2) .eq. 1)) ) &
                        cycle jgloop
                jg_sw = cgp(jgr)%iswitch !switching atom in charge group jg
                j3 = 3*jg_sw-3
                dx = x(i3+1) - x(j3+1)
                dy = x(i3+2) - x(j3+2)
                dz = x(i3+3) - x(j3+3)
                dx = dx - boxlength(1)*nint( dx*inv_boxl(1) )
                dy = dy - boxlength(2)*nint( dy*inv_boxl(2) )
                dz = dz - boxlength(3)*nint( dz*inv_boxl(3) )
                r2 = dx**2 + dy**2 + dz**2
	
! skip if outside cutoff
                if ( r2 .le. rcut2 ) then 
!inside cutoff
!check if more memory is needed
!$omp critical
                        if(nbpp_cgp_pair .eq. size(nbpp_cgp, 1)) call reallocate_nbpp_cgp
!add the charge group pair
                        nbpp_cgp_pair = nbpp_cgp_pair + 1
                        nbpp_cgp(nbpp_cgp_pair)%i = ig_sw !the switching atoms of the charge groups in the pair
                        nbpp_cgp(nbpp_cgp_pair)%j = jg_sw
!$omp end critical
ialoop:                 do ia = cgp(ig)%first, cgp(ig)%last
! for every atom in the charge group ig (of the outermost loop):
                                i = cgpatom(ia)
!skip if q-atom
                                if ( iqatom(i).ne.0 ) cycle ialoop

jaloop:                         do ja = cgp(jgr)%first, cgp(jgr)%last
! for every atom in the charge group jg (innermost loop)
                                        j = cgpatom(ja)
!                --- q-atom ? ---
                                        if ( iqatom(j).ne.0 ) cycle jaloop
! make sure each pair is only counted once
                                        if ( ig .eq. jgr .and. i .ge. j ) cycle jaloop
! do they interact?
                                        if(.not.pp_precomp(i,j)%set) cycle jaloop
! if out of space then make more space
!$omp critical
                                        if (nbpp_pair .eq. calculation_assignment%pp%max) call reallocate_nonbondlist_pp
! all tests passed, add the pair
                                        nbpp_pair = nbpp_pair + 1
                                        nbpp(nbpp_pair)%i = i
                                        nbpp(nbpp_pair)%j = j 
                                        nbpp(nbpp_pair)%vdWA = pp_precomp(i,j)%vdWA
                                        nbpp(nbpp_pair)%vdWB = pp_precomp(i,j)%vdWB
                                        nbpp(nbpp_pair)%elec = pp_precomp(i,j)%elec
                                        nbpp(nbpp_pair)%cgp_pair = nbpp_cgp_pair !which pair of charge groups the atom pair belongs to
!$omp end critical
                                end do jaloop
                        end do ialoop
                elseif ((r2 .le. RcLRF2) .or. (rcLRF .eq. -one)) then
! outside pp-cutoff but inside LRF cut-off use LRF
!ig : jg calculation
                        call lrf_update(ig,jgr)
!jg : ig calculations
                        call lrf_update(jgr,ig)
                end if ! outside cutoff

        end do jgloop
#ifdef USE_GRID
	end if ! interaction
	end do kgridloop
	end do jgridloop
	end do igridloop
#endif
end do igloop
#if defined (PROFILING)
profile(3)%time = profile(3)%time + rtime() - start_loop_time
#endif

end subroutine nbpplist_box_lrf
!-----------------------------------------------------------------------
!******PWchanged 2001-10-01
subroutine nbpw_count(npw, npwcgp)
! arguments
integer						:: npw
integer						:: npwcgp(:)

! local variables
integer						:: i,ig,jg,ia,ja,i3,j3
real(kind=prec)						:: rcut2,r2
integer						:: LJ_code

!******PWadded variables 2001-10-01

real(kind=prec)						:: dx, dy, dz

!	This routine makes a list of non-bonded solute-solvent atom pairs
!	excluding Q-atoms.

! uses the global variables:
!  Rcpw, ncgp, cgp, excl, nwat, nat_solute, x, cgpatom, iqatom, ljcod, crg, iaclib

npw = 0
rcut2 = Rcpw*Rcpw

igloop: do ig = 1, ncgp_solute
! for each charge group of the protein:
npwcgp(ig) = 0

! skip if excluded charge group
ia = cgp(ig)%iswitch
if ( .not.use_PBC .and. excl(ia) ) cycle igloop

i3 = 3*ia-3

jgloop: do jg = 1, nwat
  ! for each water molecule:
ja = nat_solute + solv_atom*jg-(solv_atom-1)
if(.not. use_PBC .and. excl(ja) ) cycle jgloop ! skip excluded waters

j3 = 3*ja-3

  if( .not. use_PBC ) then
        r2 = ( x(i3+1) -x(j3+1) )**2 &
                +( x(i3+2) -x(j3+2) )**2 &
                +( x(i3+3) -x(j3+3) )**2
  else
        dx = x(i3+1) - x(j3+1)
        dy = x(i3+2) - x(j3+2)
        dz = x(i3+3) - x(j3+3)
        dx = dx - boxlength(1)*nint( dx*inv_boxl(1) )
        dy = dy - boxlength(2)*nint( dy*inv_boxl(2) )
        dz = dz - boxlength(3)*nint( dz*inv_boxl(3) )
        r2 = dx**2 + dy**2 + dz**2
  end if

! skip if outside cutoff
if ( r2 .gt. rcut2 ) cycle jgloop

ialoop:  do ia = cgp(ig)%first, cgp(ig)%last
! for each atom in the protein charge group:
i = cgpatom(ia)

! skip if q-atom
if ( iqatom(i)/=0 ) cycle ialoop

jaloop:	do ja = nat_solute + solv_atom*jg-(solv_atom-1), nat_solute + solv_atom*jg
          ! for every atom of the water molecule:

          ! calculate LJ_code for the pair
          LJ_code = ljcod(iac(i),iac(ja))

          ! skip pairs with zero interaction
          if((crg(i) * crg(ja) == zero) &
                  .and. &
                  (iaclib(iac(i))%avdw(LJ_code)*iaclib(iac(ja))%avdw(LJ_code) == zero) &
                  .and. &
                  (iaclib(iac(i))%bvdw(LJ_code)*iaclib(iac(ja))%bvdw(LJ_code) == zero)) &
                  cycle jaloop

          ! count the pair
          npw = npw + 1
          npwcgp(ig) = npwcgp(ig) + 1

        end do jaloop
end do ialoop
end do jgloop
end do igloop

end subroutine nbpw_count

!-----------------------------------------------------------------------

subroutine nbpwlis2
! local variables
integer						:: i,ig,jg,ia,ja,i3,j3, inside,jgr,j
real(kind=prec)						:: rcut2,r2
integer                                                 ::iagrid, igrid, jgrid, kgrid, gridnum, LJ_code
#ifdef _OPENMP
integer :: quotient, remainder
#endif
#if defined (PROFILING)
real(kind=prec)						:: start_loop_time
start_loop_time = rtime()
#endif

!	This routine makes a list of non-bonded solute-solvent atom pairs
!	excluding Q-atoms.

!$omp single
nbpw_pair = 0
!$omp end single
!$omp barrier
rcut2 = Rcpw*Rcpw

#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (calculation_assignment%pw%end - calculation_assignment%pw%start + 1)/threads_num
remainder = MOD(calculation_assignment%pw%end - calculation_assignment%pw%start + 1, threads_num)
mp_start = thread_id * quotient + calculation_assignment%pw%start + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif

igloop: do ig = mp_start,mp_end
#else
igloop: do ig = calculation_assignment%pw%start, calculation_assignment%pw%end
#endif
!	   --- excluded group ? ---
        ia = cgp(ig)%iswitch
        if ( excl(ia) ) cycle igloop
#ifdef USE_GRID
        iagrid = pp_igrid(ig)
igridloop:	do igrid = 1, pw_ndim
jgridloop:	 do jgrid = 1, pw_ndim
kgridloop:	  do kgrid = 1, pw_ndim
			gridnum = gridnum + 1
			if (.not. (grid_pw_int(iagrid,igrid,jgrid,kgrid))) cycle kgridloop
jgloop: do jg = 1, grid_pw_ngrp(gridnum)
                jgr = grid_pw_grp(gridnum,jg)
#else
jgloop: do jgr = 1, nwat
#endif
                ja = nat_solute + solv_atom*jgr-(solv_atom-1)
                if ( excl(ja) ) cycle jgloop ! skip excluded waters
                j3 = 3*ja-3

!	      --- outside cutoff ? ---
                inside = 0
                ia = cgp(ig)%first
                do while ((ia .le. cgp(ig)%last) .and. (inside .eq. 0))
                        i = cgpatom(ia)
                        i3 = 3*i-3
                        r2 = ( x(i3+1) -x(j3+1) )**2 &
                                +( x(i3+2) -x(j3+2) )**2 &
                                +( x(i3+3) -x(j3+3) )**2

                        if ( r2 .le. rcut2 ) then
! inside cutoff, raise the flag
                                inside = 1
                        end if
                        ia = ia + 1
                end do
                if (inside .eq. 0) cycle jgloop

ialoop:         do ia = cgp(ig)%first, cgp(ig)%last
!	         --- q-atom ? ---
                        i = cgpatom(ia)
                        if ( iqatom(i)/=0 ) cycle ialoop
! if out of space then make more space
!$omp critical
                        if (nbpw_pair .gt. calculation_assignment%pw%max-solv_atom) call reallocate_nonbondlist_pw
jaloop:                 do j = 1, solv_atom
                                ja = nat_solute + (solv_atom*jgr) - solv_atom + j
! add the pair
                                nbpw_pair = nbpw_pair + 1
                                nbpw(nbpw_pair)%i = i
                                nbpw(nbpw_pair)%j = ja
                                nbpw(nbpw_pair)%vdWA = pw_precomp(i,j)%vdWA
                                nbpw(nbpw_pair)%vdWB = pw_precomp(i,j)%vdWB
                                nbpw(nbpw_pair)%elec = pw_precomp(i,j)%elec
                                nbpw(nbpw_pair)%cgp_pair = nbpw_cgp_pair
                        end do jaloop
!$omp end critical
                end do ialoop
        end do jgloop
#ifdef USE_GRID
end do kgridloop
end do jgridloop
end do igridloop
#endif
end do igloop
#if defined (PROFILING)
profile(4)%time = profile(4)%time + rtime() - start_loop_time
#endif

end subroutine nbpwlis2

!------------------------------------------------------------------------------

!******PWadded 2001-10-18
subroutine nbpwlis2_box
  ! local variables
  integer						:: i,ig,jg,ia,ja,i3,j3,ig_atom, inside, jgr, LJ_code,j
  real(kind=prec)						:: rcut2,r2
  real(kind=prec)						:: dx, dy,dz
  integer                                                 ::iagrid, igrid, jgrid, kgrid, gridnum
  ! for periodic boundary conditions
  !	This routine makes a list of non-bonded solute-solvent atom pairs
  !	excluding Q-atoms.
#ifdef _OPENMP
integer :: quotient, remainder
#endif
#if defined (PROFILING)
real(kind=prec)						:: start_loop_time
start_loop_time = rtime()
#endif
!$omp single
nbpw_pair = 0
nbpw_cgp_pair = 0
!$omp end single
!$omp barrier
rcut2 = Rcpw*Rcpw

#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (calculation_assignment%pw%end - calculation_assignment%pw%start + 1)/threads_num
remainder = MOD(calculation_assignment%pw%end - calculation_assignment%pw%start + 1, threads_num)
mp_start = thread_id * quotient + calculation_assignment%pw%start + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif

igloop: do ig = mp_start,mp_end
#else
igloop: do ig = calculation_assignment%pw%start, calculation_assignment%pw%end
#endif
#ifdef USE_GRID
	iagrid = pw_igrid(ig)
	gridnum = 0
igridloop: do igrid = 1, pw_ndim
jgridloop:  do jgrid = 1, pw_ndim
kgridloop:   do kgrid = 1, pw_ndim
		gridnum = gridnum + 1
		if (.not. (grid_pw_int(iagrid,igrid,jgrid,kgrid))) cycle kgridloop

jgloop:	do jg = 1, grid_pw_ngrp(gridnum)
                jgr = grid_pw_grp(gridnum,jg)
#else
jgloop: do jgr = 1, nwat
#endif
                ja = nat_solute + solv_atom*jgr-(solv_atom-1)
                j3 = 3*ja-3
!	      --- outside cutoff ? ---
		inside = 0
		ig_atom = cgp(ig)%first
		do while ((ig_atom .le. cgp(ig)%last) .and. (inside .eq. 0))
                   	i = cgpatom(ig_atom)
                  	i3 = 3*i-3
        		dx = x(i3+1) - x(j3+1)
        		dy = x(i3+2) - x(j3+2)
        		dz = x(i3+3) - x(j3+3)
        		dx = dx - boxlength(1)*nint( dx*inv_boxl(1) )
        		dy = dy - boxlength(2)*nint( dy*inv_boxl(2) )
        		dz = dz - boxlength(3)*nint( dz*inv_boxl(3) )
        		r2 = dx**2 + dy**2 + dz**2
		   
                        if ( r2 .le. rcut2 ) then
! inside cutoff, raise the flag
                                inside = 1
!$omp critical
                                if (nbpw_cgp_pair .eq. size(nbpw_cgp, 1)) call reallocate_nbpw_cgp
                                nbpw_cgp_pair = nbpw_cgp_pair + 1
                                nbpw_cgp(nbpw_cgp_pair)%i = i
                                nbpw_cgp(nbpw_cgp_pair)%j = ja
!$omp end critical
                        end if
                        ig_atom = ig_atom + 1 !ia = ia + 1
                end do
		if (inside .eq. 0) cycle jgloop
ialoop:         do ia = cgp(ig)%first, cgp(ig)%last
!	         --- q-atom ? ---
                        i = cgpatom(ia)
                        if ( iqatom(i)/=0 ) cycle ialoop
! if out of space then make more space
!$omp critical
                        if (nbpw_pair .gt. calculation_assignment%pw%max-solv_atom) call reallocate_nonbondlist_pw
jaloop:                 do j = 1, solv_atom
                                ja = nat_solute + (solv_atom*jgr) - solv_atom + j
! for every atom of the water molecule:
! add the pair
                                nbpw_pair = nbpw_pair + 1
                                nbpw(nbpw_pair)%i = i
                                nbpw(nbpw_pair)%j = ja
                                nbpw(nbpw_pair)%vdWA = pw_precomp(i,j)%vdWA
                                nbpw(nbpw_pair)%vdWB = pw_precomp(i,j)%vdWB
                                nbpw(nbpw_pair)%elec = pw_precomp(i,j)%elec
                                nbpw(nbpw_pair)%cgp_pair = nbpw_cgp_pair
                        end do jaloop
!$omp end critical
                end do ialoop
        end do jgloop
#ifdef USE_GRID
	end do kgridloop
	end do jgridloop
	end do igridloop
#endif
end do igloop
#if defined (PROFILING)
profile(4)%time = profile(4)%time + rtime() - start_loop_time
#endif

end subroutine nbpwlis2_box
!-----------------------------------------------------------------------
subroutine nbpwlis2_box_lrf
  ! local variables
  integer						:: i,ig,jg,ia,ja,i3,j3,ig_atom, inside, jgr, LJ_code, j
  real(kind=prec)						:: rcut2,r2
  real(kind=prec)						:: dx, dy,dz
  !LRF
  real(kind=prec)						:: RcLRF2, field0, field1, field2
  integer						:: jg_cgp, inside_LRF, is3

  integer                                                 ::iagrid, igrid, jgrid, kgrid, gridnum
  
  ! for periodic boundary conditions
  !	This routine makes a list of non-bonded solute-solvent atom pairs
  !	excluding Q-atoms.
#ifdef _OPENMP
integer :: quotient, remainder
#endif
#if defined (PROFILING)
real(kind=prec)						:: start_loop_time
start_loop_time = rtime()
#endif
!$omp single
  nbpw_pair = 0
  nbpw_cgp_pair = 0
!$omp end single
!$omp barrier
  rcut2 = Rcpw*Rcpw
  RcLRF2 = RcLRF*RcLRF

#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (calculation_assignment%pw%end - calculation_assignment%pw%start + 1)/threads_num
remainder = MOD(calculation_assignment%pw%end - calculation_assignment%pw%start + 1, threads_num)
mp_start = thread_id * quotient + calculation_assignment%pw%start + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif

igloop: do ig = mp_start,mp_end
#else
igloop: do ig = calculation_assignment%pw%start, calculation_assignment%pw%end
#endif
#ifdef USE_GRID
	iagrid = pw_igrid(ig)
igridloop:	do igrid = 1, pw_ndim
jgridloop:	 do jgrid = 1, pw_ndim
kgridloop:	  do kgrid = 1, pw_ndim
			gridnum = gridnum + 1
			if (.not. (grid_pw_int(iagrid,igrid,jgrid,kgrid))) then
ggloop:			 do jg = 1, grid_pw_ngrp(gridnum)
			  jgr = grid_pw_grp(gridnum,jg)
			  ja = nat_solute + solv_atom*jgr-(solv_atom-1)
			  jg_cgp = iwhich_cgp(ja)
			  call lrf_update(ig,jg_cgp)
			  call lrf_update(jg_cgp,ig)
			 end do ggloop
			else
jgloop:	do jg = 1, grid_pw_ngrp(gridnum)
                jgr = grid_pw_grp(gridnum,jg)
#else
jgloop: do jgr = 1, nwat
#endif
                ja = nat_solute + solv_atom*jgr-(solv_atom-1)
                jg_cgp = iwhich_cgp(ja)
                j3 = 3*ja-3
!	      --- outside cutoff ? ---
                inside = 0
                inside_LRF = 0
                ig_atom = cgp(ig)%first
                do while ((ig_atom .le. cgp(ig)%last) .and. (inside .eq. 0))
                        i = cgpatom(ig_atom)
                        i3 = 3*i-3
                        dx = x(i3+1) - x(j3+1)
                        dy = x(i3+2) - x(j3+2)
                        dz = x(i3+3) - x(j3+3)
                        dx = dx - boxlength(1)*nint( dx*inv_boxl(1) )
                        dy = dy - boxlength(2)*nint( dy*inv_boxl(2) )
                        dz = dz - boxlength(3)*nint( dz*inv_boxl(3) )
                        r2 = dx**2 + dy**2 + dz**2
                        if ( r2 .le. rcut2 ) then
! inside cutoff, raise the flag
                                inside = 1
!$omp critical
                                if (nbpw_cgp_pair .eq. size(nbpw_cgp, 1)) call reallocate_nbpw_cgp
                                nbpw_cgp_pair = nbpw_cgp_pair + 1
                                nbpw_cgp(nbpw_cgp_pair)%i = i
                                nbpw_cgp(nbpw_cgp_pair)%j = ja
!$omp end critical
                        elseif ((r2 .le. RcLRF2) .or. (RcLRF .eq. -one)) then
                                inside_LRF = 1
                        end if
                        ig_atom = ig_atom + 1 !ia = ia + 1
                end do
                if (inside .eq. 1) then
ialoop:                 do ia = cgp(ig)%first, cgp(ig)%last
!	         --- q-atom ? ---
                        i = cgpatom(ia)
                        if ( iqatom(i)/=0 ) cycle ialoop
! if out of space then make more space
!$omp critical
                                if (nbpw_pair .gt. calculation_assignment%pw%max - 3) call reallocate_nonbondlist_pw
jaloop:                         do j = 1, solv_atom
                                        ja = nat_solute + (solv_atom*jgr) - solv_atom + j 
! for every atom of the water molecule:
! add the pair
                                        nbpw_pair = nbpw_pair + 1
                                        nbpw(nbpw_pair)%i = i
                                        nbpw(nbpw_pair)%j = ja 
                                        nbpw(nbpw_pair)%vdWA = pw_precomp(i,j)%vdWA
                                        nbpw(nbpw_pair)%vdWB = pw_precomp(i,j)%vdWB
                                        nbpw(nbpw_pair)%elec = pw_precomp(i,j)%elec
                                        nbpw(nbpw_pair)%cgp_pair = nbpw_cgp_pair
                                end do jaloop
!$omp end critical

                        end do ialoop
                elseif((inside_LRF .eq. 1) .and. (inside .eq. 0)) then   
! outside pw-cutoff but inside LRF cut-off: use LRF
!solut : solvent
                        call lrf_update(ig,jg_cgp)
!solvent : solut	
                        call lrf_update(jg_cgp,ig)
                end if
        end do jgloop
#ifdef USE_GRID
	end if ! interaction
end do kgridloop
end do jgridloop
end do igridloop
#endif
end do igloop
#if defined (PROFILING)
profile(4)%time = profile(4)%time + rtime() - start_loop_time
#endif

end subroutine nbpwlis2_box_lrf
!-----------------------------------------------------------------------
subroutine nbpwlis2_lrf
! local variables
integer						:: i,j,ig,jg,jg_cgp,ia,ja,i3,j3,inside,is,is3,jgr
real(kind=prec)						:: rcut2,r2
real(kind=prec)						::	RcLRF2
integer                                                 ::iagrid, igrid, jgrid, kgrid, gridnum, LJ_code
#ifdef _OPENMP
integer :: quotient, remainder
#endif

!	This routine makes a list of non-bonded solute-solvent atom pairs
!	excluding Q-atoms.
#if defined (PROFILING)
real(kind=prec)						:: start_loop_time
start_loop_time = rtime()
#endif
!$omp single
nbpw_pair = 0
!$omp end single
!$omp barrier
rcut2 = Rcpw*Rcpw
RcLRF2 = RcLRF*RcLRF
#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (calculation_assignment%pw%end - calculation_assignment%pw%start + 1)/threads_num
remainder = MOD(calculation_assignment%pw%end - calculation_assignment%pw%start + 1, threads_num)
mp_start = thread_id * quotient + calculation_assignment%pw%start + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif

igloop: do ig = mp_start,mp_end
#else
igloop: do ig = calculation_assignment%pw%start, calculation_assignment%pw%end
#endif

!	   --- excluded group ? ---
        is = cgp(ig)%iswitch
        if ( excl(is) ) cycle igloop
        is3 = 3*is-3
#ifdef USE_GRID
        iagrid = pw_igrid(ig)
        gridnum = 0
igridloop:	do igrid = 1, pw_ndim
jgridloop:	 do jgrid = 1, pw_ndim
kgridloop:	  do kgrid = 1, pw_ndim
			gridnum = gridnum + 1
			if (.not. (grid_pw_int(iagrid,igrid,jgrid,kgrid))) then
ggloop:			do jg = 1, grid_pw_ngrp(gridnum)
			 jgr = grid_pw_grp(gridnum,jg)
			 ja = nat_solute + solv_atom*jgr-(solv_atom-1)
			 if(excl(ja)) cycle ggloop ! skip excluded waters
			 jg_cgp = iwhich_cgp(ja)
			 call lrf_update(ig,jg_cgp)
			 call lrf_update(jg_cgp,ig)
			end do ggloop
			else
jgloop:	do jg = 1, grid_pw_ngrp(gridnum)
	        jgr = grid_pw_grp(gridnum,jg)
#else
jgloop: do jgr = 1, nwat
#endif
                ja = nat_solute + solv_atom*jgr-(solv_atom-1)
                if(excl(ja)) cycle jgloop
                jg_cgp = iwhich_cgp(ja)
                j3 = 3*ja-3

!	      --- outside cutoff ? ---
                inside = 0
                ia = cgp(ig)%first
                do while ((ia .le. cgp(ig)%last) .and. (inside .eq. 0 ))
                        i = cgpatom(ia)
                        i3 = 3*i-3
                        r2 = ( x(i3+1) -x(j3+1) )**2 &
                                +( x(i3+2) -x(j3+2) )**2 &
                                +( x(i3+3) -x(j3+3) )**2
                        if ( r2 .le. rcut2 ) then
! inside cutoff, raise the flag
                                inside = 1
                        end if
                        ia = ia + 1
                end do

                if ( inside .eq. 1 ) then
ialoop:                 do ia = cgp(ig)%first, cgp(ig)%last
!	            --- q-atom ? ---
                                i = cgpatom(ia)
                                if ( iqatom(i)/=0 ) cycle ialoop
! if out of space then make more space
!$omp critical
                                if (nbpw_pair .gt. calculation_assignment%pw%max-solv_atom) call reallocate_nonbondlist_pw
jaloop:                         do j = 1, solv_atom 
                                        ja = nat_solute + (solv_atom*jgr) - solv_atom + j
! for every atom of the water molecule:
! add the pair
                                        nbpw_pair = nbpw_pair + 1
                                        nbpw(nbpw_pair)%i = i
                                        nbpw(nbpw_pair)%j = ja
                                        nbpw(nbpw_pair)%vdWA = pw_precomp(i,j)%vdWA
                                        nbpw(nbpw_pair)%vdWB = pw_precomp(i,j)%vdWB
                                        nbpw(nbpw_pair)%elec = pw_precomp(i,j)%elec
                                        nbpw(nbpw_pair)%cgp_pair = nbpw_cgp_pair
                                end do jaloop
!$omp end critical
                        end do ialoop

                elseif(r2 .le. RcLRF2) then   
! outside pw-cutoff but inside LRF cut-off: use LRF
                	call lrf_update(ig,jg_cgp)
	                call lrf_update(jg_cgp,ig)
                end if
        end do jgloop
#ifdef USE_GRID
end if ! interaction
end do kgridloop
end do jgridloop
end do igridloop
#endif
end do igloop
#if defined (PROFILING)
profile(4)%time = profile(4)%time + rtime() - start_loop_time
#endif

end subroutine nbpwlis2_lrf

!-----------------------------------------------------------------------

subroutine nbpwlist
! local variables
integer						:: i,ig,jg,ia,ja,i3,j3,jgr,j
real(kind=prec)						:: rcut2,r2
integer						:: LJ_code
integer                                                 ::iagrid, igrid, jgrid, kgrid, gridnum

! for use with spherical boundary

!	This routine makes a list of non-bonded solute-solvent atom pairs
!	excluding Q-atoms.
#ifdef _OPENMP
integer :: quotient, remainder
#endif
#if defined (PROFILING)
real(kind=prec)						:: start_loop_time
start_loop_time = rtime()
#endif

! reset nbpw_pair
!$omp single
nbpw_pair = 0
!$omp end single
!$omp barrier
rcut2 = Rcpw*Rcpw
#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (calculation_assignment%pw%end - calculation_assignment%pw%start + 1)/threads_num
remainder = MOD(calculation_assignment%pw%end - calculation_assignment%pw%start + 1, threads_num)
mp_start = thread_id * quotient + calculation_assignment%pw%start + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif

igloop: do ig = mp_start,mp_end
#else
igloop: do ig = calculation_assignment%pw%start, calculation_assignment%pw%end
#endif
! for every charge group:

! skip excluded groups
        ia = cgp(ig)%iswitch
        if ( excl(ia) ) cycle igloop
        i3 = 3*ia-3
#ifdef USE_GRID
        iagrid = pw_igrid(ig)
        gridnum = 0
igridloop:	do igrid = 1, pw_ndim
jgridloop:	 do jgrid = 1, pw_ndim
kgridloop:	  do kgrid = 1, pw_ndim
			gridnum = gridnum + 1
			if (.not. (grid_pw_int(iagrid,igrid,jgrid,kgrid))) cycle kgridloop
jgloop:	do jg = 1, grid_pw_ngrp(gridnum)
	        jgr = grid_pw_grp(gridnum,jg)
#else
jgloop: do jgr = 1, nwat
#endif
! for every water molecule in grid:
                ja = nat_solute + solv_atom*jgr-(solv_atom-1)
                if(excl(ja)) cycle jgloop ! skip excluded waters

                j3 = 3*ja-3
                r2 = ( x(i3+1) -x(j3+1) )**2 &
                        +( x(i3+2) -x(j3+2) )**2 &
                        +( x(i3+3) -x(j3+3) )**2

! skip water outside cutoff
                if ( r2 .gt. rcut2 ) cycle jgloop

ialoop:         do ia = cgp(ig)%first, cgp(ig)%last
! for every atom in the charge group:

! find the atom index of the atom in the charge group
                        i = cgpatom(ia)
!	skip q-atoms
                        if ( iqatom(i)/=0 ) cycle ialoop
!$omp critical
! if out of space then make more space
                        if (nbpw_pair .gt. calculation_assignment%pw%max - solv_atom) call reallocate_nonbondlist_pw
jaloop:                 do j = 1, solv_atom
                                ja = nat_solute + (solv_atom*jgr) - solv_atom + j
! for every atom of the water molecule:
! add the pair
                                nbpw_pair = nbpw_pair + 1
                                nbpw(nbpw_pair)%i = i
                                nbpw(nbpw_pair)%j = ja 
                                nbpw(nbpw_pair)%vdWA = pw_precomp(i,j)%vdWA
                                nbpw(nbpw_pair)%vdWB = pw_precomp(i,j)%vdWB
                                nbpw(nbpw_pair)%elec = pw_precomp(i,j)%elec
                        end do jaloop
!$omp end critical
                end do ialoop
        end do jgloop
#ifdef USE_GRID
end do kgridloop
end do jgridloop
end do igridloop
#endif
end do igloop
#if defined (PROFILING)
profile(4)%time = profile(4)%time + rtime() - start_loop_time
#endif

end subroutine nbpwlist

!-----------------------------------------------------------------------
!******PWadded 2001-10-18

subroutine nbpwlist_box
  ! local variables
  integer						:: i,ig,jg,ia,ja,i3,j3,ig_sw,jg_sw,jgr,j
  real(kind=prec)						:: rcut2,r2
  integer						:: LJ_code
  real(kind=prec)						:: dx, dy, dz
  integer                                                 ::iagrid, igrid, jgrid, kgrid, gridnum
  ! For use with periodic boundary conditions
  !	This routine makes a list of non-bonded solute-solvent atom pairs
  !	excluding Q-atoms.
#ifdef _OPENMP
integer :: quotient, remainder
#endif
#if defined (PROFILING)
real(kind=prec)						:: start_loop_time
start_loop_time = rtime()
#endif

! reset nbpw_pair
!$omp single
nbpw_pair = 0
nbpw_cgp_pair = 0
!$omp end single
!$omp barrier
rcut2 = Rcpw*Rcpw
#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (calculation_assignment%pw%end - calculation_assignment%pw%start + 1)/threads_num
remainder = MOD(calculation_assignment%pw%end - calculation_assignment%pw%start + 1, threads_num)
mp_start = thread_id * quotient + calculation_assignment%pw%start + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif

igloop: do ig = mp_start,mp_end
#else
igloop: do ig = calculation_assignment%pw%start, calculation_assignment%pw%end
#endif
! for every charge group:
        ig_sw = cgp(ig)%iswitch
        i3 = 3*ig_sw-3
#ifdef USE_GRID
        iagrid = pw_igrid(ig)
        gridnum = 0
igridloop:	do igrid = 1, pw_ndim
jgridloop:	 do jgrid = 1, pw_ndim
kgridloop:	  do kgrid = 1, pw_ndim
			gridnum = gridnum + 1
			if (.not. (grid_pw_int(iagrid,igrid,jgrid,kgrid))) cycle kgridloop
jgloop: do jg = 1, grid_pw_ngrp(gridnum)
	        jgr = grid_pw_grp(gridnum,jg)
#else
jgloop: do jgr = 1, nwat
#endif
! for every water molecule:
                jg_sw = nat_solute + solv_atom*jgr-(solv_atom-1)
                j3 = 3*jg_sw-3
                dx = x(i3+1) - x(j3+1)
                dy = x(i3+2) - x(j3+2)
                dz = x(i3+3) - x(j3+3)
                dx = dx - boxlength(1)*nint( dx*inv_boxl(1) )
                dy = dy - boxlength(2)*nint( dy*inv_boxl(2) )
                dz = dz - boxlength(3)*nint( dz*inv_boxl(3) )
                r2 = dx**2 + dy**2 + dz**2

! skip water outside cutoff
                if ( r2 .gt. rcut2 ) cycle jgloop

	  !inside cut-off
!check if charge group pair list is big enough
!$omp critical
                if(nbpw_cgp_pair .eq. size(nbpw_cgp, 1) ) call reallocate_nbpw_cgp
                nbpw_cgp_pair = nbpw_cgp_pair + 1
                nbpw_cgp(nbpw_cgp_pair)%i = ig_sw  !solute
                nbpw_cgp(nbpw_cgp_pair)%j = jg_sw  !water
!$omp end critical
ialoop:         do ia = cgp(ig)%first, cgp(ig)%last
! for every atom in the charge group:

! find the atom index of the atom in the charge group
                        i = cgpatom(ia)
!	skip q-atoms
                        if ( iqatom(i)/=0 ) cycle ialoop
! if out of space then make more space
!$omp critical
                        if (nbpw_pair .gt. calculation_assignment%pw%max - solv_atom) call reallocate_nonbondlist_pw
jaloop:                 do j = 1, solv_atom
                                ja = nat_solute + (solv_atom*jgr) - solv_atom + j
! for every atom of the water molecule:
! add the pair
                                nbpw_pair = nbpw_pair + 1
                                nbpw(nbpw_pair)%i = i
                                nbpw(nbpw_pair)%j = ja 
                                nbpw(nbpw_pair)%vdWA = pw_precomp(i,j)%vdWA
                                nbpw(nbpw_pair)%vdWB = pw_precomp(i,j)%vdWB
                                nbpw(nbpw_pair)%elec = pw_precomp(i,j)%elec
                                nbpw(nbpw_pair)%cgp_pair = nbpw_cgp_pair
                        end do jaloop
!$omp end critical
                end do ialoop
        end do jgloop
#ifdef USE_GRID
end do kgridloop
end do jgridloop
end do igridloop
#endif
  end do igloop
#if defined (PROFILING)
profile(4)%time = profile(4)%time + rtime() - start_loop_time
#endif

end subroutine nbpwlist_box
!-------------------------------------------------------------------------------------

subroutine nbpwlist_lrf
! local variables
integer						:: i,j,ig,jg,jg_cgp,ia,ja,i3,j3,is,is3,jgr
real(kind=prec)						:: rcut2,r2
integer						:: LJ_code
real(kind=prec)						::	RcLRF2
integer                                                 ::iagrid, igrid, jgrid, kgrid, gridnum
#ifdef _OPENMP
integer :: quotient, remainder
#endif
#if defined (PROFILING)
real(kind=prec)						:: start_loop_time
start_loop_time = rtime()
#endif

!	This routine makes a list of non-bonded solute-solvent atom pairs
!	excluding Q-atoms.

! reset nbpw_pair
!$omp single
nbpw_pair = 0
!$omp end single
!$omp barrier
rcut2 = Rcpw*Rcpw
RcLRF2 = RcLRF*RcLRF
#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (calculation_assignment%pw%end - calculation_assignment%pw%start + 1)/threads_num
remainder = MOD(calculation_assignment%pw%end - calculation_assignment%pw%start + 1, threads_num)
mp_start = thread_id * quotient + calculation_assignment%pw%start + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif

igloop: do ig = mp_start,mp_end
#else
igloop: do ig = calculation_assignment%pw%start, calculation_assignment%pw%end
#endif
! for every charge group:

! skip excluded groups
        is = cgp(ig)%iswitch
        if ( excl(is) ) cycle igloop
        is3 = 3*is-3
#ifdef USE_GRID
        iagrid = pw_igrid(ig)
        gridnum = 0
igridloop:	do igrid = 1, pw_ndim
jgridloop:	 do jgrid = 1, pw_ndim
kgridloop:	  do kgrid = 1, pw_ndim
			gridnum = gridnum + 1
			if (.not. (grid_pw_int(iagrid,igrid,jgrid,kgrid))) then
ggloop:			do jg = 1, grid_pw_ngrp(gridnum)
			 jgr = grid_pw_grp(gridnum,jg)
			 ja = nat_solute + solv_atom*jgr-(solv_atom-1)
			 if(excl(ja)) cycle ggloop ! skip excluded waters
			 jg_cgp = iwhich_cgp(ja)
			 call lrf_update(ig,jg_cgp)
			 call lrf_update(jg_cgp,ig)
			end do ggloop
			else
jgloop:	do jg = 1, grid_pw_ngrp(gridnum)
	        jgr = grid_pw_grp(gridnum,jg)
#else
jgloop: do jgr = 1,nwat
#endif
                ja = nat_solute + solv_atom*jgr-(solv_atom-1)
                jg_cgp = iwhich_cgp(ja)
                if(excl(ja)) cycle jgloop
                j3 = 3*ja-3
                r2 = ( x(is3+1) -x(j3+1) )**2 &
                     +( x(is3+2) -x(j3+2) )**2 &
                     +( x(is3+3) -x(j3+3) )**2

                if ( r2 .le. rcut2 ) then
! within the cutoff radix:

ialoop:                 do ia = cgp(ig)%first, cgp(ig)%last
                                i = cgpatom(ia)
! skip q-atoms
                                if ( iqatom(i)/=0 ) cycle ialoop
!$omp critical
! if out of space then make more space
                                if (nbpw_pair .gt. calculation_assignment%pw%max - solv_atom) call reallocate_nonbondlist_pw
jaloop:                         do j = 1, solv_atom
                                        ja = nat_solute + (solv_atom*jgr) - solv_atom + j
! add the pair
                                        nbpw_pair = nbpw_pair + 1
                                        nbpw(nbpw_pair)%i = i
                                        nbpw(nbpw_pair)%j = ja 
                                        nbpw(nbpw_pair)%vdWA = pw_precomp(i,j)%vdWA
                                        nbpw(nbpw_pair)%vdWB = pw_precomp(i,j)%vdWB
                                        nbpw(nbpw_pair)%elec = pw_precomp(i,j)%elec
                                end do jaloop
!$omp end critical
                        end do ialoop
                elseif(r2 .le. RcLRF2) then   
! outside pw-cutoff but inside LRF cut-off: use LRF
                        call lrf_update(ig,jg_cgp)
                        call lrf_update(jg_cgp,ig)
                end if
        end do jgloop
#ifdef USE_GRID
end if ! interaction
end do kgridloop
end do jgridloop
end do igridloop
#endif
end do igloop
#if defined (PROFILING)
profile(4)%time = profile(4)%time + rtime() - start_loop_time
#endif

end subroutine nbpwlist_lrf
!---------------LRF version of PW PBC-----------------------
subroutine nbpwlist_box_lrf
  ! local variables
  integer						:: i,ig,jg,ia,ja,i3,j3,ig_sw,jg_sw,jgr
  real(kind=prec)						:: rcut2,r2
  integer						:: LJ_code
  real(kind=prec)						:: dx, dy, dz
  ! LRF
  real(kind=prec)						:: RcLRF2
  integer						:: jg_cgp, j, is3
  integer                                                 ::iagrid, igrid, jgrid, kgrid, gridnum
#ifdef _OPENMP
integer :: quotient, remainder
#endif
#if defined (PROFILING)
real(kind=prec)                                         :: start_loop_time
start_loop_time = rtime()
#endif

  ! For use with periodic boundary conditions
  !	This routine makes a list of non-bonded solute-solvent atom pairs
  !	excluding Q-atoms.

  ! reset nbpw_pair
!$omp single
  nbpw_pair = 0
  nbpw_cgp_pair = 0
!$omp end single
!$omp barrier
  rcut2 = Rcpw*Rcpw
  RcLRF2 = RcLRF*RcLRF
#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (calculation_assignment%pw%end - calculation_assignment%pw%start + 1)/threads_num
remainder = MOD(calculation_assignment%pw%end - calculation_assignment%pw%start + 1, threads_num)
mp_start = thread_id * quotient + calculation_assignment%pw%start + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif

igloop: do ig = mp_start,mp_end
#else
igloop: do ig = calculation_assignment%pw%start, calculation_assignment%pw%end
#endif
	! for every charge group:
        ig_sw = cgp(ig)%iswitch
        i3 = 3*ig_sw-3
#ifdef USE_GRID
        iagrid = pw_igrid(ig)
        gridnum = 0
igridloop:	do igrid = 1, pw_ndim
jgridloop:	 do jgrid = 1, pw_ndim
kgridloop:	  do kgrid = 1, pw_ndim
			gridnum = gridnum + 1
			if (.not. (grid_pw_int(iagrid,igrid,jgrid,kgrid))) then
ggloop:			do jg = 1, grid_pw_ngrp(gridnum)
			 jgr = grid_pw_grp(gridnum,jg)
			 jg_sw = nat_solute + (solv_atom*jgr-(solv_atom-1))
			 jg_cgp = iwhich_cgp(jg_sw)
			 call lrf_update(ig,jg_cgp)
			 call lrf_update(jg_cgp,ig)
			end do ggloop
			else
jgloop: do jg = 1, grid_pw_ngrp(gridnum)
	        jgr = grid_pw_grp(gridnum,jg)
#else
jgloop: do jgr = 1, nwat
#endif
                jg_sw = nat_solute + (solv_atom*jgr-(solv_atom-1))
                jg_cgp = iwhich_cgp(jg_sw)
                j3 = 3*jg_sw-3
                dx = x(i3+1) - x(j3+1)
                dy = x(i3+2) - x(j3+2)
                dz = x(i3+3) - x(j3+3)
                dx = dx - boxlength(1)*nint( dx*inv_boxl(1) )
                dy = dy - boxlength(2)*nint( dy*inv_boxl(2) )
                dz = dz - boxlength(3)*nint( dz*inv_boxl(3) )
                r2 = dx**2 + dy**2 + dz**2

! skip water outside cutoff
                if ( r2 .le. rcut2 ) then
!inside cut-off
!check if charge group pair list is big enough
!$omp critical
                        if(nbpw_cgp_pair .eq. size(nbpw_cgp, 1) ) call reallocate_nbpw_cgp
                        nbpw_cgp_pair = nbpw_cgp_pair + 1
                        nbpw_cgp(nbpw_cgp_pair)%i = ig_sw  !solute
                        nbpw_cgp(nbpw_cgp_pair)%j = jg_sw  !water
!$omp end critical
ialoop:                 do ia = cgp(ig)%first, cgp(ig)%last
! for every atom in the charge group:
! find the atom index of the atom in the charge group
                                i = cgpatom(ia)
!	skip q-atoms
                                if ( iqatom(i)/=0 ) cycle ialoop
! if out of space then make more space
!$omp critical
                                if (nbpw_pair .gt. calculation_assignment%pw%max - solv_atom) call reallocate_nonbondlist_pw
jaloop:                         do j = 1, solv_atom
                                        ja = nat_solute + (solv_atom*jgr) - solv_atom + j
! for every atom of the solvent molecule:
! add the pair
                                        nbpw_pair = nbpw_pair + 1
                                        nbpw(nbpw_pair)%i = i
                                        nbpw(nbpw_pair)%j = ja 
                                        nbpw(nbpw_pair)%vdWA = pw_precomp(i,j)%vdWA
                                        nbpw(nbpw_pair)%vdWB = pw_precomp(i,j)%vdWB
                                        nbpw(nbpw_pair)%elec = pw_precomp(i,j)%elec
                                        nbpw(nbpw_pair)%cgp_pair = nbpw_cgp_pair
                                end do jaloop
!$omp end critical
                        end do ialoop
                elseif((r2 .le. RcLRF2) .or. (RcLRF .eq. -one)) then   
! outside pw-cutoff but inside LRF cut-off: use LRF
!solut : solvent
                        call lrf_update(ig,jg_cgp)
!solvent : solut	
                        call lrf_update(jg_cgp,ig)
                end if
        end do jgloop
#ifdef USE_GRID
end if !interaction
end do kgridloop
end do jgridloop
end do igridloop
#endif
  end do igloop
#if defined (PROFILING)
profile(4)%time = profile(4)%time + rtime() - start_loop_time
#endif    

end subroutine nbpwlist_box_lrf
!-------------------------------------------------------------------------------------
subroutine make_qconn

integer						::	 i, iq, is, ia

allocate(qconn(nstates,nat_solute, nqat))
qconn(:,:,:) = 9

do iq = 1, nqat
        qconn(:, iqseq(iq), iq) = 1
end do

do iq = 1, nqat
        do is = 1, nstates
                i = iqseq(iq)
                call find_bonded(origin=i, current=i, level=1, state=is)
        end do
end do

!modify matrix to take special exclusions into account
do i = 1, nexspec
        iq = iqatom(exspec(i)%i)
        if(iq > 0) then
                do is = 1, nstates
                        if(exspec(i)%flag(is)) then
                                qconn(is, exspec(i)%j, iq) = 0 !exclude by setting to 0
                        end if
                end do
        end if
        iq = iqatom(exspec(i)%j)
        if(iq > 0) then
                do is = 1, nstates
                        if(exspec(i)%flag(is)) then
                                qconn(is, exspec(i)%i, iq) = 0 !exclude by setting to 0
                        end if
                end do
        end if
end do

end subroutine make_qconn


!------------------------------------------------------------------------------


recursive subroutine find_bonded(origin, current, level, state)
!args
integer, intent(in)			::	origin, current, level, state
!locals
integer						::	b, newcurrent, newlevel

!find q-atom connectivity using the bond list and the q-bond list
!shaken bonds (code -1) must be taken into account, but not 
!redefined bonds in the topology
do b = 1, nbonds_solute
        if(bnd(b)%cod == 0) cycle !skip redefined (but not shaken)
        if(bnd(b)%i  == current) then
                newlevel = level + 1 
                newcurrent = bnd(b)%j
                if(qconn(state, newcurrent, iqatom(origin)) > newlevel)  then
                        qconn(state, newcurrent, iqatom(origin)) = newlevel
                        if(newlevel < 4) then
                                call find_bonded(origin, newcurrent, newlevel, state)
                        end if
                end if
        elseif(bnd(b)%j  == current) then
                newlevel = level + 1 
                newcurrent = bnd(b)%i
                if(qconn(state, newcurrent, iqatom(origin)) > newlevel)  then
                        qconn(state, newcurrent, iqatom(origin)) = newlevel
                        if(newlevel < 4) then
                                call find_bonded(origin, newcurrent, newlevel, state)
                        end if
                end if
        end if
end do
do b = 1, nqbond
        if(qbnd(b)%cod(state) > 0) then
                if(qbnd(b)%i  == current) then
                        newlevel = level + 1 
                        newcurrent = qbnd(b)%j
                        if(qconn(state, newcurrent, iqatom(origin)) > newlevel)  then
                                qconn(state, newcurrent, iqatom(origin)) = newlevel
                                if(newlevel < 4) then
                                        call find_bonded(origin, newcurrent, newlevel, state)
                                end if
                        end if
                elseif(qbnd(b)%j  == current) then
                        newlevel = level + 1 
                        newcurrent = qbnd(b)%i
                        if(qconn(state, newcurrent, iqatom(origin)) > newlevel)  then
                                qconn(state, newcurrent, iqatom(origin)) = newlevel
                                if(newlevel < 4) then
                                        call find_bonded(origin, newcurrent, newlevel, state)
                                end if
                        end if
                end if
        end if
end do

end subroutine find_bonded


!------------------------------------------------------------------------------


integer function nbqq_count()

integer						::	iq, j, jq, is

nbqq_pair(:) = 0

!count Q-Q
do iq = 1, nqat - 1
        do jq = iq + 1, nqat
                do is = 1, nstates
                        if(qconn(is, iqseq(jq), iq) > 3) then
                                nbqq_pair(is) = nbqq_pair(is)+1
                        end if
                end do
        end do
end do
!count Q-non-Q
do j = 1, nat_solute
        if(iqatom(j) > 0) cycle
        if(any(qconn(:,j,:) <= 3)) then 
                !bonded or angled to at least one Q-atom in any state
                do iq = 1, nqat
                        do is = 1, nstates
                                if(qconn(is, j, iq) >= 4) then
                                        nbqq_pair(is) = nbqq_pair(is)+1
                                end if
                        end do
                end do
        end if
end do

nbqq_count = maxval(nbqq_pair(:))
end function nbqq_count

!---------------------------------------------------------------------------------

subroutine nbqqlist
integer						::	iq, j, jq, is, i, k,l, ia , ja
real(kind=prec)                     :: el_scale
logical                     :: set

nbqq_pair(:) = 0

nbqqp_pair(:) = 0
!list Q-Q
do iq = 1, nqat - 1
        ia = iqseq(iq)
        do jq = iq + 1, nqat
                ja = iqseq(jq)
                do is = 1, nstates
                        if(.not.qq_precomp(iq,jq,is)%set) cycle
                        nbqq_pair(is) = nbqq_pair(is)+1
                        nbqq(nbqq_pair(is),is)%iq = iq
                        nbqq(nbqq_pair(is),is)%jq = jq
                        nbqq(nbqq_pair(is),is)%vdWA = qq_precomp(iq,jq,is)%vdWA
                        nbqq(nbqq_pair(is),is)%vdWB = qq_precomp(iq,jq,is)%vdWB
                        nbqq(nbqq_pair(is),is)%elec = qq_precomp(iq,jq,is)%elec
                        nbqq(nbqq_pair(is),is)%score = qq_precomp(iq,jq,is)%score
                        nbqq(nbqq_pair(is),is)%soft = qq_precomp(iq,jq,is)%soft
                end do
        end do
end do

!list variable Q-P atoms
do ja = 1, nat_solute
        if(iqatom(ja) .ne. 0) cycle
        if(any(qconn(:,ja,:) <= 3)) then
                !bonded or angled to at least one Q-atom
                do iq = 1, nqat
                        do is = 1, nstates
                                if(.not.qp_precomp(ja,iq,is)%set) cycle
! the list style for the qp interactions has now been changed to be state
! dependent, so we can use a single list to calculate all qp interactions
! as this one here will be static, the minimum number of qp pairs will be kept
! the same always and is saved under the new variable nbqqp_pair(:)
! the lists will now be allocated based on the maximum size needed that is
! determined after running the precomputation
                                nbqqp_pair(is) = nbqqp_pair(is)+1
                                nbqqp(nbqqp_pair(is),is)%iq = iq
                                nbqqp(nbqqp_pair(is),is)%jq = ja
                                nbqqp(nbqqp_pair(is),is)%vdWA = qp_precomp(ja,iq,is)%vdWA
                                nbqqp(nbqqp_pair(is),is)%vdWB = qp_precomp(ja,iq,is)%vdWB
                                nbqqp(nbqqp_pair(is),is)%elec = qp_precomp(ja,iq,is)%elec
                                nbqqp(nbqqp_pair(is),is)%score = qp_precomp(ja,iq,is)%score
                        end do
                end do
        end if
end do

end subroutine nbqqlist

!-----------------------------------------------------------------------
!******PWchanged 2001-10-01
subroutine nbqp_count(nqp, nqpcgp)
! arguments
integer						:: nqp
integer						:: nqpcgp(:)

! local variables
integer						:: ig,ia,i,j,iq,i3
real(kind=prec)						:: rcut2,r2

!******PWadded variables 2001-10-01

real(kind=prec)						:: dx, dy, dz

!	This routine counts non-bonded atom pairs involving
!	*one* Q-atom and *one* non-Q-atom, where the latter is *not connected*
!	(meaning not bonded or angled) to any Q-atom.
!
!	( i , j ) pairs correspond to
!	( iq, j ) with first index being a Q-atom with the *Q-atom numbering*, &
!	          and the second index is the non-Q-atom.

!	Note that for PBC the default is to count all atoms as interacting with each other
!	if the user did not specify a q - q cutoff of > 0

nqp = 0
rcut2 = Rcq*Rcq



if(nqat .eq. 0) return

! --- solute - Q-atoms

igloop: do ig = 1, ncgp_solute

nqpcgp(ig) = 0

! skip if excluded group
ia = cgp(ig)%iswitch
if ( .not. use_PBC .and. excl(ia) ) cycle igloop
i3 = 3*ia-3

!******PWadded if 2001-10-01

if( .not. use_PBC ) then
r2 = ( x(i3+1) - xpcent(1) )**2 &
        +( x(i3+2) - xpcent(2) )**2 &
        +( x(i3+3) - xpcent(3) )**2
else
	if (Rcq .gt. zero ) then
        dx = x(i3+1) - x(3*qswitch-2)
        dy = x(i3+2) - x(3*qswitch-1)
        dz = x(i3+3) - x(3*qswitch)
        dx = dx - boxlength(1)*nint( dx*inv_boxl(1) )
        dy = dy - boxlength(2)*nint( dy*inv_boxl(2) )
        dz = dz - boxlength(3)*nint( dz*inv_boxl(3) )
        r2 = dx**2 + dy**2 + dz**2
end if
end if

! skip if outside cutoff
if ( ( r2 .gt. rcut2 ) .or. (use_PBC.and.(Rcq.gt.zero))) cycle igloop

ialoop: do ia = cgp(ig)%first, cgp(ig)%last
  i = cgpatom(ia)

  if(iqatom(i) .ne. 0) cycle
  ! check if already on qq list
  if(any(qconn(:,i,:) <= 3)) cycle ialoop

  ! count the pairs
  nqp = nqp + nqat
  nqpcgp(ig) = nqpcgp(ig) + nqat

end do ialoop
end do igloop

end subroutine nbqp_count

!-----------------------------------------------------------------------
!******PWchanged 2001-10-01
subroutine nbqw_count(nqw, nqwmol)
! arguments
integer						:: nqw
integer						:: nqwmol(:)

! local variables
integer						:: ig,ia,i,j,iq,i3
real(kind=prec)						:: rcut2,r2

!******PWadded variables

real(kind=prec)						:: dx, dy, dz

!	This routine counts water molecules that interact with q-atoms
!	Note that in PBC the default is no cutoff (rcq < 0), so we just count all of them
nqw = 0
rcut2 = Rcq*Rcq


if(nqat .eq. 0) return

! --- solvent - Q-atoms

iwloop: do ig = 1, nwat
nqwmol(ig) = 0
ia = nat_solute + solv_atom*ig-(solv_atom-1)
if(.not. use_PBC .and. excl(ia)) cycle iwloop ! skip excluded waters
i3 = 3*ia-3

!******PWadded if-statement 2001-10-01
if( .not. use_PBC ) then
        r2 = ( x(i3+1) - xpcent(1) )**2 &
                +( x(i3+2) - xpcent(2) )**2 &
                +( x(i3+3) - xpcent(3) )**2
else
	if (rcq .gt. zero) then
        dx = x(i3+1) - x(3*qswitch-2)
        dy = x(i3+2) - x(3*qswitch-1)
        dz = x(i3+3) - x(3*qswitch)
        dx = dx - boxlength(1)*nint( dx*inv_boxl(1) )
        dy = dy - boxlength(2)*nint( dy*inv_boxl(2) )
        dz = dz - boxlength(3)*nint( dz*inv_boxl(3) )
        r2 = dx**2 + dy**2 + dz**2
end if
end if

! skip if outside cutoff
if ( (r2 .lt. rcut2 ) .or. (use_PBC.and.(rcq.lt.zero))) then
        nqw = nqw + 1
        nqwmol(ig) = solv_atom*nqat
end if
end do iwloop

end subroutine nbqw_count

!-----------------------------------------------------------------------

subroutine nbqplis2
! local variables
integer						:: ig,ia,i,j,iq,i3,nl,inside,is
real(kind=prec)						:: rcut2,r2
integer						:: xspec
logical, save					:: list_done
#ifdef _OPENMP
integer :: quotient, remainder
#endif

!	This routine makes a list of non-bonded atom pairs involving
!	*one* Q-atom and *one* non-Q-atom, where the latter is *not connected*
!	(meaning not bonded or angled) to any Q-atom.
!
!	( i , j ) pairs correspond to
!	( iq, j ) with first index being a Q-atom with the *Q-atom numbering*, &
!	          and the second index is the non-Q-atom.

! uses the global variables:
!  nbqs_pair, Rcq, cgp, excl, cgpatom, x, xpcent, nqat, iqseq, 
!  qconn, calculation_assignment%qs%max, nbqs, ljcod, nwat, nat_solute

!don't remake pair list if q-atoms interact with all solute atoms (rcq>rexcl_o)
!and list already made
#if defined (PROFILING)
real(kind=prec)                                         :: start_loop_time
start_loop_time = rtime()
#endif



! we don't remake the q - p list in PBC either now, by just including all freaking atoms into it
! if the q-q cutoff is set to < 0
! this is made default for the standart inputs, you can still set it to something else if you want to
if(list_done .and. (((.not.use_PBC).and.(Rcq > rexcl_o)))) return
!$omp single
nbqp_pair = 0
!$omp end single
!$omp barrier
rcut2 = Rcq*Rcq


if(nqat .eq. 0) return


! --- solute - Q-atoms

#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (calculation_assignment%qp%end - calculation_assignment%qp%start + 1)/threads_num
remainder = MOD(calculation_assignment%qp%end - calculation_assignment%qp%start + 1, threads_num)
mp_start = thread_id * quotient + calculation_assignment%qp%start + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif
igloop: do ig = mp_start,mp_end
#else
igloop: do ig = calculation_assignment%qp%start, calculation_assignment%qp%end
#endif
! for every assigned charge group:

! skip if excluded group
        ia = cgp(ig)%iswitch
        if ( excl(ia) ) cycle igloop

! check cutoff
        inside = 0
        ia = cgp(ig)%first
        do while ((ia .le. cgp(ig)%last) .and. (inside .eq. 0))
                i = cgpatom(ia)
                i3 = 3*i-3


                r2 = ( x(i3+1) - xpcent(1) )**2 &
                        +( x(i3+2) - xpcent(2) )**2 &
                        +( x(i3+3) - xpcent(3) )**2
                if ( r2 .le. rcut2 ) then
                        inside = 1
                end if
                ia = ia + 1
        end do
        if (inside .eq. 0) cycle igloop

ialoop: do ia = cgp(ig)%first, cgp(ig)%last
                i = cgpatom(ia)
! see if it is already on the list after nbqq
                if(any(qconn(:,i,:) <= 3)) cycle ialoop
!$omp critical
                ! if out of space then make more space
                if (nbqp_pair .ge. calculation_assignment%qp%max-nqat) call reallocate_nonbondlist_qp
qaloop:         do iq = 1, nqat
                        nbqp_pair = nbqp_pair + 1
                        do is = 1, nstates
! store the pair
                                nbqp(nbqp_pair,is)%i = iq
                                nbqp(nbqp_pair,is)%j = i
                                nbqp(nbqp_pair,is)%vdWA = qp_precomp(i,iq,is)%vdWA
                                nbqp(nbqp_pair,is)%vdWB = qp_precomp(i,iq,is)%vdWB
                                nbqp(nbqp_pair,is)%elec = qp_precomp(i,iq,is)%elec
                                nbqp(nbqp_pair,is)%score = qp_precomp(i,iq,is)%score
                        end do
                end do qaloop
!$omp end critical
        end do ialoop
end do igloop
!$omp single
list_done = .true. !save this value
!$omp end single
#if defined (PROFILING)
profile(5)%time = profile(5)%time + rtime() - start_loop_time
#endif

end subroutine nbqplis2



!----------------------------------------------------------------------------

subroutine nbqplis2_box
!  ! local variables
  integer						:: ig,ia,i,j,iq,i3,nl,inside,ig_atom,is
  real(kind=prec)						:: rcut2,r2
  integer						:: xspec
  real(kind=prec)						:: dx, dy, dz
  logical,save						:: list_done = .false.
 
  ! for periodic boundary conditions
  !	This routine makes a list of non-bonded atom pairs involving
  !	*one* Q-atom and *one* non-Q-atom, where the latter is *not connected*
  !	(meaning not bonded or angled) to any Q-atom.
  !
  !	( i , j ) pairs correspond to
  !	( iq, j ) with first index being a Q-atom with the *Q-atom numbering*, &
  !	          and the second index is the non-Q-atom.

  ! uses the global variables:
  !  nbqs_pair, Rcq, cgp, excl, cgpatom, x, xpcent, nqat, iqseq, 
  !  qconn, calculation_assignment%qs%max, nbqs, ljcod, nwat, nat_solute

! now with more parallel
#ifdef _OPENMP
integer :: quotient, remainder
#endif
#if defined (PROFILING)
real(kind=prec)                                         :: start_loop_time
start_loop_time = rtime()
#endif

! we don't remake the q - p list in PBC either now, by just including all freaking atoms into it
! if the q-q cutoff is set to < 0
! this is made default for the standart inputs, you can still set it to something else if you want to
if(list_done .and. (((use_PBC).and.(Rcq.lt.zero)))) return
if(nqat .eq. 0) return

!$omp single
nbqp_pair = 0
nbqp_cgp_pair = 0
!$omp end single
!$omp barrier
rcut2 = Rcq*Rcq

 ! --- solute - Q-atoms
#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (calculation_assignment%qp%end - calculation_assignment%qp%start + 1)/threads_num
remainder = MOD(calculation_assignment%qp%end - calculation_assignment%qp%start + 1, threads_num)
mp_start = thread_id * quotient + calculation_assignment%qp%start + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif

igloop: do ig = mp_start,mp_end
#else
igloop: do ig = calculation_assignment%qp%start, calculation_assignment%qp%end
#endif
    ! for every assigned charge group:

	! check cutoff
	! if Rcq < 0 -> no cutoff, skip distance calculation here
	if (((use_PBC).and.(Rcq.lt.zero))) then
	        inside = 1
	else
	        inside = 0
	end if
	ig_atom = cgp(ig)%first
	do while ((ig_atom .le. cgp(ig)%last) .and. (inside .eq. 0))
                i = cgpatom(ig_atom)
                i3 = 3*i-3
                dx = x(i3+1) - x(3*qswitch-2)
                dy = x(i3+2) - x(3*qswitch-1)
                dz = x(i3+3) - x(3*qswitch)
                dx = dx - boxlength(1)*nint( dx*inv_boxl(1) )
                dy = dy - boxlength(2)*nint( dy*inv_boxl(2) )
                dz = dz - boxlength(3)*nint( dz*inv_boxl(3) )
                r2 = dx**2 + dy**2 + dz**2

                if ( r2 .le. rcut2 ) then
		        inside = 1
!$omp critical
                        if(nbqp_cgp_pair .eq. size(nbqp_cgp, 1) ) call reallocate_nbqp_cgp
                        nbqp_cgp_pair = nbqp_cgp_pair + 1
                        nbqp_cgp(nbqp_cgp_pair)%i = i !leave %j empty, equals qswitch
!$omp end critical
                end if

                ig_atom = ig_atom + 1 !ia = ia + 1
        end do
	if (inside .eq. 0) cycle igloop

ialoop: do ia = cgp(ig)%first, cgp(ig)%last
                i = cgpatom(ia)

! check if already on qq list
                if(any(qconn(:,i,:) <= 3)) cycle ialoop

!$omp critical
                ! if out of space then make more space
                if (nbqp_pair .ge. calculation_assignment%qp%max-nqat) call reallocate_nonbondlist_qp
qaloop:         do iq = 1, nqat
                        nbqp_pair = nbqp_pair + 1
                        do is = 1, nstates
! store the pair
                                nbqp(nbqp_pair,is)%i = iq
                                nbqp(nbqp_pair,is)%j = i
                                nbqp(nbqp_pair,is)%vdWA = qp_precomp(i,iq,is)%vdWA
                                nbqp(nbqp_pair,is)%vdWB = qp_precomp(i,iq,is)%vdWB
                                nbqp(nbqp_pair,is)%elec = qp_precomp(i,iq,is)%elec
                                nbqp(nbqp_pair,is)%score = qp_precomp(i,iq,is)%score
                                nbqp(nbqp_pair,is)%cgp_pair = nbqp_cgp_pair
                        end do
                end do qaloop
!$omp end critical
        end do ialoop
end do igloop
#if defined (PROFILING)
profile(5)%time = profile(5)%time + rtime() - start_loop_time
#endif
!$omp single
list_done = .true. !save this value
!$omp end single

end subroutine nbqplis2_box
!-----------------------------------------------------------------------

subroutine nbqplist
! local variables
integer						:: ig,ia,i,j,iq,i3,nl,is
real(kind=prec)						:: rcut2,r2
integer						::	xspec
logical, save				::	list_done = .false.
#ifdef _OPENMP
integer :: quotient, remainder
#endif
#if defined (PROFILING)
real(kind=prec)						:: start_loop_time
start_loop_time = rtime()
#endif


!For spherical boundary	
!	This routine makes a list of non-bonded atom pairs involving
!	*one* Q-atom and *one* non-Q-atom, where the latter is *not connected*
!	(meaning not bonded or angled) to any Q-atom.
!
!	( i , j ) pairs correspond to
!	( iq, j ) with first index being a Q-atom with the *Q-atom numbering*, &
!	          and the second index is the non-Q-atom.

! global variables used:
!  nbqs_pair, Rcq, cgp, excl, x, xpcent, cgpatom, nqat, iqseq,
!  calculation_assignment%qp%max, nbqs, ljcod, nwat, nat_solute

!don't remake pair list if q-atoms interact with all solute atoms (rcq>rexcl_o)
!and list already made
if(list_done .and. (Rcq .gt. rexcl_o)) return
if(nqat .eq. 0) return
!$omp single
nbqp_pair = 0
!$omp end single
!$omp barrier
rcut2 = Rcq*Rcq


! --- solute - Q-atoms

#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (calculation_assignment%qp%end - calculation_assignment%qp%start + 1)/threads_num
remainder = MOD(calculation_assignment%qp%end - calculation_assignment%qp%start + 1, threads_num)
mp_start = thread_id * quotient + calculation_assignment%qp%start + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif
igloop: do ig = mp_start,mp_end
#else
igloop: do ig = calculation_assignment%qp%start, calculation_assignment%qp%end
#endif

! skip if excluded group
        ia = cgp(ig)%iswitch
        if ( excl(ia) ) cycle igloop
        i3 = 3*ia-3
        r2 = ( x(i3+1) - xpcent(1) )**2 &
                +( x(i3+2) - xpcent(2) )**2 &
                +( x(i3+3) - xpcent(3) )**2

! skip if outside cutoff
        if ( r2 .gt. rcut2 ) cycle igloop

ialoop: do ia = cgp(ig)%first, cgp(ig)%last

                i = cgpatom(ia)

  ! check if already on qq list
                if(any(qconn(:,i,:) <= 3)) cycle ialoop
! if out of space then make more space
!$omp critical
                if (nbqp_pair .ge. calculation_assignment%qp%max-nqat) call reallocate_nonbondlist_qp
qaloop:         do iq = 1, nqat
                        nbqp_pair = nbqp_pair + 1
                        do is = 1, nstates
! store the pair
                                nbqp(nbqp_pair,is)%i = iq
                                nbqp(nbqp_pair,is)%j = i
                                nbqp(nbqp_pair,is)%vdWA = qp_precomp(i,iq,is)%vdWA
                                nbqp(nbqp_pair,is)%vdWB = qp_precomp(i,iq,is)%vdWB
                                nbqp(nbqp_pair,is)%elec = qp_precomp(i,iq,is)%elec
                                nbqp(nbqp_pair,is)%score = qp_precomp(i,iq,is)%score
                        end do
                end do qaloop
!$omp end critical
        end do ialoop
end do igloop
#if defined (PROFILING)
profile(5)%time = profile(5)%time + rtime() - start_loop_time
#endif
!$omp single
list_done = .true. !save this value
!$omp end single
end subroutine nbqplist

!-----------------------------------------------------------------------

!******PWadded 2001-10-18
subroutine nbqplist_box
  ! local variables
  integer						:: ig,ia,i,j,iq,i3,nl,ig_atom,is
  real(kind=prec)						:: rcut2,r2
  integer						::	xspec
  real(kind=prec)						:: dx, dy, dz
  integer						:: ga, gb, inside
  logical,save						:: list_done = .false.
#ifdef _OPENMP
integer :: quotient, remainder
#endif
#if defined (PROFILING)
real(kind=prec)                                         :: start_loop_time
start_loop_time = rtime()
#endif

  !For periodic boundary conditions	
  !	This routine makes a list of non-bonded atom pairs involving
  !	*one* Q-atom and *one* non-Q-atom, where the latter is *not connected*
  !	(meaning not bonded or angled) to any Q-atom.
  !
  !	( i , j ) pairs correspond to
  !	( iq, j ) with first index being a Q-atom with the *Q-atom numbering*, &
  !	          and the second index is the non-Q-atom.

  ! global variables used:
  !  nbqs_pair, Rcq, cgp, x, xpcent, cgpatom, nqat, iqseq,
  !  calculation_assignment%qs%max, nbqs, ljcod, nwat, nat_solute

! we don't remake the q - p list in PBC either now, by just including all freaking atoms into it
! if the q-q cutoff is set to < 0
! this is made default for the standart inputs, you can still set it to something else if you want to
if(list_done .and. (((use_PBC).and.(Rcq.lt.zero)))) return

if(nqat .eq. 0) return
!$omp single
nbqp_pair = 0
!$omp end single
!$omp barrier
rcut2 = Rcq*Rcq
! --- solute - Q-atoms
#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (calculation_assignment%qp%end - calculation_assignment%qp%start + 1)/threads_num
remainder = MOD(calculation_assignment%qp%end - calculation_assignment%qp%start + 1, threads_num)
mp_start = thread_id * quotient + calculation_assignment%qp%start + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif

igloop: do ig = mp_start,mp_end
#else
igloop: do ig = calculation_assignment%qp%start, calculation_assignment%qp%end
#endif
        if (((use_PBC).and.(Rcq.lt.zero))) then
                inside = 1
        else
                inside = 0
        end if
        ig_atom = cgp(ig)%first
        do while ((ig_atom .le. cgp(ig)%last) .and. (inside .eq. 0))
                i = cgpatom(ig_atom)
                i3 = 3*i-3
                dx = x(i3+1) - x(3*qswitch-2)
                dy = x(i3+2) - x(3*qswitch-1)
                dz = x(i3+3) - x(3*qswitch)
                dx = dx - boxlength(1)*nint( dx*inv_boxl(1) )
                dy = dy - boxlength(2)*nint( dy*inv_boxl(2) )
                dz = dz - boxlength(3)*nint( dz*inv_boxl(3) )
                r2 = dx**2 + dy**2 + dz**2

                if ( r2 .le. rcut2 ) then
                        inside = 1

!$omp critical
                        if( nbqp_cgp_pair .eq. size(nbqp_cgp, 1) ) call reallocate_nbqp_cgp
                        nbqp_cgp_pair = nbqp_cgp_pair + 1
                        nbqp_cgp(nbqp_cgp_pair)%i = i !leave %j empty, equals qswitch
!$omp end critical
                end if
                ig_atom = ig_atom + 1 !ia = ia + 1
        end do
        if (inside .eq. 0) cycle igloop

ialoop: do ia = cgp(ig)%first, cgp(ig)%last

                i = cgpatom(ia)
! check if already on qq list
                if(any(qconn(:,i,:) <= 3)) cycle ialoop
! if out of space then make more space
!$omp critical
                ! if out of space then make more space
                if (nbqp_pair .ge. calculation_assignment%qp%max-nqat) call reallocate_nonbondlist_qp
qaloop:         do iq = 1, nqat
                        nbqp_pair = nbqp_pair + 1
                        do is = 1, nstates
! store the pair
                                nbqp(nbqp_pair,is)%i = iq
                                nbqp(nbqp_pair,is)%j = i
                                nbqp(nbqp_pair,is)%vdWA = qp_precomp(i,iq,is)%vdWA
                                nbqp(nbqp_pair,is)%vdWB = qp_precomp(i,iq,is)%vdWB
                                nbqp(nbqp_pair,is)%elec = qp_precomp(i,iq,is)%elec
                                nbqp(nbqp_pair,is)%score = qp_precomp(i,iq,is)%score
                        end do
                end do qaloop
!$omp end critical
        end do ialoop
end do igloop

#if defined (PROFILING)
profile(5)%time = profile(5)%time + rtime() - start_loop_time
#endif
!$omp single
list_done = .true.
!$omp end single

end subroutine nbqplist_box
!------------------------------------------------------------------------------------

subroutine nbqwlist

! local variables
integer						:: ig,ia,i,i3,is,iq,istate
real(kind=prec)						:: rcut2,r2
logical,save					:: list_done = .false.

! we do not know the number of molecules in the solvent, so we need to get them first before allocating
! this is done in prep_sim and stored in solv_atom

!	This routine makes a list of water molecules within rcq from xpcent,
! i.e. the q-atom - water non-bond lists which implicitly includes all
! q-atoms with all atoms of the listed water
! waters may not have bonded interactions with q-atoms !
#ifdef _OPENMP
integer :: quotient, remainder
#endif
#if defined (PROFILING)
real(kind=prec)                                         :: start_loop_time
start_loop_time = rtime()
#endif



!We don't have to remake the list if q-atoms interact with all waters
!(rcq > rexcl_o) and list already made (nbwq_pair >0)
if( list_done .and.  (Rcq .gt. rexcl_o)) return

if(nqat .eq. 0) return
!$omp single
nbqw_pair = 0
!$omp end single
!$omp barrier
rcut2 = Rcq*Rcq



#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (calculation_assignment%qw%end - calculation_assignment%qw%start + 1)/threads_num
remainder = MOD(calculation_assignment%qw%end - calculation_assignment%qw%start + 1, threads_num)
mp_start = thread_id * quotient + calculation_assignment%qw%start + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif

iwloop: do ig = mp_start,mp_end
#else
iwloop: do ig = calculation_assignment%qw%start, calculation_assignment%qw%end
#endif
        ia = nat_solute + solv_atom*ig-(solv_atom-1)
        if( excl(ia) ) cycle iwloop ! skip excluded waters
        i3 = 3*ia-3
        r2 = ( x(i3+1) - xpcent(1) )**2 &
                +( x(i3+2) - xpcent(2) )**2 &
                +( x(i3+3) - xpcent(3) )**2


! store if inside cutoff
        if ( r2 <= rcut2 ) then
!$omp critical
qaloop:         do iq = 1, nqat
                nbqw_pair = nbqw_pair + solv_atom
solvloop:               do is = 1, solv_atom
                                nbqw(nbqw_pair-(solv_atom-is),1:nstates)%i=iq
                                nbqw(nbqw_pair-(solv_atom-is),1:nstates)%j=ia +(is-1)
                                nbqw(nbqw_pair-(solv_atom-is),1:nstates)%elec  = qw_precomp(iq,is,1:nstates)%elec 
                                nbqw(nbqw_pair-(solv_atom-is),1:nstates)%vdWA  = qw_precomp(iq,is,1:nstates)%vdWA
                                nbqw(nbqw_pair-(solv_atom-is),1:nstates)%vdWB  = qw_precomp(iq,is,1:nstates)%vdWB
                                nbqw(nbqw_pair-(solv_atom-is),1:nstates)%score = qw_precomp(iq,is,1:nstates)%score
                        end do solvloop
                end do qaloop
!$omp end critical
        end if
end do iwloop
#if defined (PROFILING)
profile(6)%time = profile(6)%time + rtime() - start_loop_time
#endif
! save this to prevent remaking of the list
!$omp single
list_done = .true.
!$omp end single

end subroutine nbqwlist

!-----------------------------------------------------------------------

!******PWadded 2001-10-18
subroutine nbqwlist_box

  ! local variables
  integer						:: ig,ia,i,i3,is,iq
  real(kind=prec)						:: rcut2,r2
  real(kind=prec)						:: dx, dy, dz
  logical,save						:: list_done = .false.
  ! we do not know the number of molecules in the solvent, so we need to get them first before allocating
  ! this is done in prep_sim and stored in solv_atom
  !	This routine makes a list of water molecules within rcq from xswitch,
  ! i.e. the q-atom - water non-bond lists which implicitly includes all
  ! q-atoms with all atoms of the listed water
  ! waters may not have bonded interactions with q-atoms !
!now with more parallel
#ifdef _OPENMP
integer :: quotient, remainder
#endif
#if defined (PROFILING)
real(kind=prec)                                         :: start_loop_time
start_loop_time = rtime()
#endif

if(list_done .and. (Rcq.lt.zero)) return
if(nqat .eq. 0) return
!$omp single
nbqw_pair = 0
!$omp end single
!$omp barrier
rcut2 = Rcq*Rcq

#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (calculation_assignment%qw%end - calculation_assignment%qw%start + 1)/threads_num
remainder = MOD(calculation_assignment%qw%end - calculation_assignment%qw%start + 1, threads_num)
mp_start = thread_id * quotient + calculation_assignment%qw%start + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif

iwloop: do ig = mp_start,mp_end
#else
iwloop: do ig = calculation_assignment%qw%start, calculation_assignment%qw%end
#endif
	ia = nat_solute + solv_atom*ig-(solv_atom-1)
	i3 = 3*ia-3

	if (Rcq.gt.zero) then
	dx = x(i3+1) - x(3*qswitch-2)
	dy = x(i3+2) - x(3*qswitch-1)
	dz = x(i3+3) - x(3*qswitch)
	dx = dx - boxlength(1)*nint( dx*inv_boxl(1) )
	dy = dy - boxlength(2)*nint( dy*inv_boxl(2) )
	dz = dz - boxlength(3)*nint( dz*inv_boxl(3) )
	r2 = dx**2 + dy**2 + dz**2
	end if
	! store if inside cutoff
	if ( ( r2 .le. rcut2 ).or. (Rcq.lt.zero)) then
!$omp critical
qaloop:	        do iq = 1, nqat
		        nbqw_pair = nbqw_pair + solv_atom
                        nbqw(nbqw_pair-(solv_atom-is),1:nstates)%i=iq
                        nbqw(nbqw_pair-(solv_atom-is),1:nstates)%j=ia +(is-1)
solvloop:	        do is = 1, solv_atom
                                nbqw(nbqw_pair-(solv_atom-is),1:nstates)%elec  = qw_precomp(iq,is,1:nstates)%elec 
                                nbqw(nbqw_pair-(solv_atom-is),1:nstates)%vdWA  = qw_precomp(iq,is,1:nstates)%vdWA
                                nbqw(nbqw_pair-(solv_atom-is),1:nstates)%vdWB  = qw_precomp(iq,is,1:nstates)%vdWB
                                nbqw(nbqw_pair-(solv_atom-is),1:nstates)%score = qw_precomp(iq,is,1:nstates)%score
		        end do solvloop
		end do qaloop
!$omp end critical
        end if
end do iwloop
#if defined (PROFILING)
profile(6)%time = profile(6)%time + rtime() - start_loop_time
#endif
!$omp single
list_done = .true.
!$omp end single
end subroutine nbqwlist_box
!-------------------------------------------------------------------------------------
!******PWchanged 2001-10-01
subroutine nbww_count(nww, nwwmol)
! arguments
integer						:: nww
integer						:: nwwmol(:)

! local variables
integer						:: iw,jw,ia,ja,i3,j3
real(kind=prec)						:: rcut2,r2

!******PWadded variables		

real(kind=prec)						:: dx, dy, dz

! This routine counts non-bonded solvent-solvent atom pairs.

nww = 0
rcut2 = Rcww*Rcww

iwloop: do iw = 1, nwat
nwwmol(iw) = 0

ia = nat_solute + solv_atom*iw-(solv_atom-1)
if(.not. use_PBC .and. excl(ia)) cycle iwloop ! skip excluded waters

i3 = 3*ia-3

jwloop: do jw = 1, nwat
  ja = nat_solute + solv_atom*jw-(solv_atom-1)
  if(.not. use_PBC .and. excl(ja)) cycle jwloop ! skip excluded waters

  j3 = 3*ja-3

  ! count each w-w pair once only
  if ( ((iw .gt. jw) .and. (mod(iw+jw,2) .eq. 0)) .or. &
           ((iw .lt. jw) .and. (mod(iw+jw,2) .eq. 1)) .or. &
           (iw .eq. jw)) &
           cycle jwloop
  !******PWadded if 2001-10-01

  if( use_PBC ) then
        dx = x(i3+1) - x(j3+1)
        dy = x(i3+2) - x(j3+2)
        dz = x(i3+3) - x(j3+3)
        dx = dx - boxlength(1)*nint( dx*inv_boxl(1) )
        dy = dy - boxlength(2)*nint( dy*inv_boxl(2) )
        dz = dz - boxlength(3)*nint( dz*inv_boxl(3) )
        r2 = dx**2 + dy**2 + dz**2
  else
        r2 = ( x(i3+1) - x(j3+1) )**2 &
           + ( x(i3+2) - x(j3+2) )**2 &
           + ( x(i3+3) - x(j3+3) )**2
  end if
  ! count the pair if inside cutoff
  if ( r2 .le. rcut2 ) then
        nww = nww + 1
        nwwmol(iw) = nwwmol(iw) + solv_atom**2
  end if

end do jwloop
end do iwloop
end subroutine nbww_count

!-----------------------------------------------------------------------
subroutine nbwwlist
! local variables
integer						:: iw,jw,ia,ja,i3,j3,la,ka,jwr
real(kind=prec)						:: rcut2,r2
integer                                                 ::iagrid, igrid, jgrid, kgrid, gridnum
! This routine makes a list of non-bonded solvent-solvent atom pairs
! The number of atoms in each solvent molecule is stored in the global solv_atom
! uses the global variables:
!  nbww_pair, Rcww, nat_solute, excl, nwat, x, nbww, calculation_assignment%ww%max
#ifdef _OPENMP
integer :: quotient, remainder
#endif
#if defined (PROFILING)
real(kind=prec)                                         :: start_loop_time
start_loop_time = rtime()
#endif

!$omp single 
nbww_pair = 0
!$omp end single
!$omp barrier
rcut2 = Rcww*Rcww
#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (calculation_assignment%ww%end - calculation_assignment%ww%start + 1)/threads_num
remainder = MOD(calculation_assignment%ww%end - calculation_assignment%ww%start + 1, threads_num)
mp_start = thread_id * quotient + calculation_assignment%ww%start + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif

iwloop: do iw = mp_start,mp_end
#else
iwloop: do iw = calculation_assignment%ww%start, calculation_assignment%ww%end
#endif

	ia = nat_solute + solv_atom*iw-(solv_atom-1)
        if (excl(ia)) cycle iwloop
        i3 = 3*ia-3
#ifdef USE_GRID
        iagrid = ww_igrid(iw)
	gridnum = 0
igridloop:	do igrid = 1, ww_ndim
jgridloop:	 do jgrid = 1, ww_ndim
kgridloop:	  do kgrid = 1, ww_ndim
			gridnum = gridnum + 1
			if (.not. (grid_ww_int(iagrid,igrid,jgrid,kgrid))) cycle kgridloop
jwloop: do jw=1, grid_ww_ngrp(gridnum)
	        jwr = grid_ww_grp(gridnum,jw)
#else
jwloop:	do jwr =1, nwat
#endif
                ja = nat_solute + solv_atom*jwr-(solv_atom-1)
                if (excl(ja)) cycle jwloop ! skip excluded waters
                j3 = 3*ja-3
! count each w-w pair once only
                if ( ((iw .gt. jwr) .and. (mod(iw+jwr,2) .eq. 0)) .or. &
                     ((iw .lt. jwr) .and. (mod(iw+jwr,2) .eq. 1)) .or. (iw .eq. jwr))  cycle jwloop

                r2    = ( x(i3+1) -x(j3+1) )**2 &
                      + ( x(i3+2) -x(j3+2) )**2 &
                      + ( x(i3+3) -x(j3+3) )**2

                if ( r2 .le. rcut2 ) then
! inside cutoff: add the pair

! check if there is enough space
!$omp critical
			if (nbww_pair .ge. (calculation_assignment%ww%max - solv_atom**2)) call reallocate_nonbondlist_ww
laloop:			do la = 1, solv_atom
kaloop:			        do ka = 1, solv_atom
                                        nbww_pair = nbww_pair + 1
                                        nbww(nbww_pair)%i = nat_solute +solv_atom*iw-(solv_atom-la)
        				nbww(nbww_pair)%j = nat_solute +solv_atom*jwr-(solv_atom-ka)
			        end do kaloop
			end do laloop
!$omp end critical
                end if
        end do jwloop
#ifdef USE_GRID
end do kgridloop
end do jgridloop
end do igridloop
#endif
end do iwloop
#if defined (PROFILING)
profile(2)%time = profile(2)%time + rtime() - start_loop_time
#endif
end subroutine nbwwlist
!------------PWadded 2001-10-18------------------------------------------
subroutine nbwwlist_box

! local variables
integer						:: iw,jw,ia,ja,i3,j3,la,ka,jwr
real(kind=prec)						:: rcut2,r2
real(kind=prec)						:: dx, dy, dz
integer                                                 ::iagrid, igrid, jgrid, kgrid, gridnum
! for periodic boundary conditions
! This routine makes a list of non-bonded solvent-solvent atom pairs
! uses the global variables:
!  nbww_pair, Rcww, nat_solute, excl, nwat, x, nbww, calculation_assignment%ww%max

! now with more parallel
#ifdef _OPENMP
integer :: quotient, remainder
#endif
#if defined (PROFILING)
real(kind=prec)                                         :: start_loop_time
start_loop_time = rtime()
#endif


!$omp single
nbww_pair = 0
!$omp end single
!$omp barrier
rcut2 = Rcww*Rcww

#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (calculation_assignment%ww%end - calculation_assignment%ww%start + 1)/threads_num
remainder = MOD(calculation_assignment%ww%end - calculation_assignment%ww%start + 1, threads_num)
mp_start = thread_id * quotient + calculation_assignment%ww%start + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif
iwloop: do iw = mp_start,mp_end
#else
iwloop: do iw = calculation_assignment%ww%start, calculation_assignment%ww%end
#endif
        ia = nat_solute + solv_atom*iw-(solv_atom-1)
        i3 = 3*ia-3
#ifdef USE_GRID
        iagrid = ww_igrid(iw)
	gridnum = 0
igridloop:	do igrid = 1, ww_ndim
jgridloop:	 do jgrid = 1, ww_ndim
kgridloop:	  do kgrid = 1, ww_ndim
			gridnum = gridnum + 1
			if (.not. (grid_ww_int(iagrid,igrid,jgrid,kgrid))) cycle kgridloop
jwloop: do jw=1, grid_ww_ngrp(gridnum)
                jwr = grid_ww_grp(gridnum,jw)
#else
jwloop:	do jwr= 1, nwat
#endif
                ja = nat_solute + solv_atom*jwr-(solv_atom-1)
                j3 = 3*ja-3
! count each w-w pair once only
                if ( ((iw .gt. jwr) .and. (mod(iw+jwr,2) .eq. 0)) .or. &
                        ((iw .lt. jwr) .and. (mod(iw+jwr,2) .eq. 1)) .or. &
                        (iw .eq. jwr)) &
                        cycle jwloop

                dx = x(i3+1) - x(j3+1)
                dy = x(i3+2) - x(j3+2)
                dz = x(i3+3) - x(j3+3)
                dx = dx - boxlength(1)*nint( dx*inv_boxl(1) )
                dy = dy - boxlength(2)*nint( dy*inv_boxl(2) )
                dz = dz - boxlength(3)*nint( dz*inv_boxl(3) )
                r2 = dx**2 + dy**2 + dz**2
                if ( r2 .le. rcut2 ) then
! inside cutoff: add the pair
! check if there is enough space
!$omp critical
                        if (nbww_pair .ge. (calculation_assignment%ww%max - solv_atom**2)) call reallocate_nonbondlist_ww
laloop:                 do la = 1, solv_atom
kaloop:                         do ka = 1, solv_atom
                                        nbww_pair = nbww_pair + 1
                                        nbww(nbww_pair)%i = nat_solute+solv_atom*iw-(solv_atom-la)
                                        nbww(nbww_pair)%j = nat_solute+solv_atom*jwr-(solv_atom-ka)
                                end do kaloop
                        end do laloop
!$omp end critical
                end if
        end do jwloop
#ifdef USE_GRID
end do kgridloop
end do jgridloop
end do igridloop
#endif
end do iwloop
#if defined (PROFILING)
profile(2)%time = profile(2)%time + rtime() - start_loop_time
#endif

end subroutine nbwwlist_box

!---------------------------------------------------------------------------

subroutine nbwwlist_lrf
! local variables
integer						:: i,j,ig,jg,iw,jw,ia,ja,i3,j3,is,is3,la,ka,jwr
real(kind=prec)						:: rcut2,r2
real(kind=prec)						:: dr(3)
real(kind=prec)						::	RcLRF2
integer                                                 ::iagrid, igrid, jgrid, kgrid, gridnum
! now with more parallel -> OMP added
#ifdef _OPENMP
integer :: quotient, remainder
#endif
#if defined (PROFILING)
real(kind=prec)						:: start_loop_time
start_loop_time = rtime()
#endif
!	This routine makes a list of non-bonded solvent-solvent atom pairs.

! uses the global variables:
!  nbww_pair, Rcww, nwat, nat_solute, excl, x, nbww, ncgp, lrf, crg, calculation_assignment%ww%max



!$omp single
nbww_pair = 0
!$omp end single
!$omp barrier
rcut2 = Rcww*Rcww
RcLRF2 = RcLRF*RcLRF

#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (calculation_assignment%ww%end - calculation_assignment%ww%start + 1)/threads_num
remainder = MOD(calculation_assignment%ww%end - calculation_assignment%ww%start + 1, threads_num)
mp_start = thread_id * quotient + calculation_assignment%ww%start + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif

iwloop: do iw = mp_start,mp_end
#else
iwloop: do iw = calculation_assignment%ww%start, calculation_assignment%ww%end
#endif
        is  = nat_solute + solv_atom*iw-(solv_atom-1)
        if( excl(is)) cycle iwloop
        is3 = 3*is-3
        ig = iwhich_cgp(is)
#ifdef USE_GRID
        iagrid = ww_igrid(iw)
        gridnum = 0
igridloop:	do igrid = 1, ww_ndim
jgridloop:	 do jgrid = 1, ww_ndim
kgridloop:	  do kgrid = 1, ww_ndim
			gridnum = gridnum + 1
			if (.not. (grid_ww_int(iagrid,igrid,jgrid,kgrid))) then
ggloop:			 do jw = 1, grid_ww_ngrp(gridnum)
			  jwr = grid_ww_grp(gridnum,jw)
			  ja = nat_solute + solv_atom*jwr-(solv_atom-1)
			  if(excl(ja)) cycle ggloop
			  jg = iwhich_cgp(ja)
                          if ( ((iw .gt. jwr) .and. (mod(iw+jwr,2) .eq. 0)) .or. &
                                  ((iw .lt. jwr) .and. (mod(iw+jwr,2) .eq. 1)).or. &
                                  (iw .eq. jwr) ) &
                                  cycle ggloop
			  call lrf_update(ig,jg)
			  call lrf_update(jg,ig)
			 end do ggloop
			else
jwloop:	do jw = 1, grid_ww_ngrp(gridnum)
	        jwr = grid_ww_grp(gridnum,jw)
#else
jwloop: do jwr = 1, nwat
#endif
                ja = nat_solute + solv_atom*jwr-(solv_atom-1)
                if(excl(ja)) cycle jwloop
                jg = iwhich_cgp(ja)
                j3 = 3*ja-3

! count each w-w pair once only
                if ( ((iw .gt. jwr) .and. (mod(iw+jwr,2) .eq. 0)) .or. &
                        ((iw .lt. jwr) .and. (mod(iw+jwr,2) .eq. 1)) .or. &
                        (iw .eq. jwr) ) &
                        cycle jwloop

                r2 = ( x(is3+1) -x(j3+1) )**2 &
                        +( x(is3+2) -x(j3+2) )**2 &
                        +( x(is3+3) -x(j3+3) )**2

                if ( r2 .le. rcut2 ) then
! inside cutoff: add the pair
! this now has to be OMP critical so the memory is not getting corrupted
!$omp critical 
! check if there is enough space
                        if (nbww_pair .ge. (calculation_assignment%ww%max - solv_atom**2)) call reallocate_nonbondlist_ww
laloop:                 do la = 1, solv_atom
kaloop:                         do ka = 1, solv_atom
                                        nbww_pair = nbww_pair + 1
                                        nbww(nbww_pair)%i = nat_solute+solv_atom*iw-(solv_atom-la)
                                        nbww(nbww_pair)%j = nat_solute+solv_atom*jwr-(solv_atom-ka)
                                end do kaloop
                        end do laloop
!$omp end critical

                elseif (r2 .le. RcLRF2) then
! outside ww-cutoff but inside LRF cut-off: use LRF
                        call lrf_update(ig,jg)
                        call lrf_update(jg,ig)
                end if
        end do jwloop
#ifdef USE_GRID
end if ! interaction
end do kgridloop
end do jgridloop
end do igridloop
#endif
end do iwloop
#if defined (PROFILING)
profile(2)%time = profile(2)%time + rtime() - start_loop_time
#endif

end subroutine nbwwlist_lrf
!--------------LRF version of PW PBC-----------------------------------
subroutine nbwwlist_box_lrf
! local variables
integer						:: i,j,ig,jg,iw,jw,ia,ja,i3,j3,is,is3,la,ka,jwr
real(kind=prec)						:: rcut2,r2
real(kind=prec)						::	RcLRF2
real(kind=prec)						::  dx, dy, dz
integer                                                 ::iagrid, igrid, jgrid, kgrid, gridnum

#ifdef _OPENMP
integer :: quotient, remainder
#endif

#if defined (PROFILING)
real(kind=prec)                                         :: start_loop_time
start_loop_time = rtime()
#endif


!	This routine makes a list of non-bonded solvent-solvent atom pairs.

! uses the global variables:
!  nbww_pair, Rcww, nwat, nat_solute, excl, x, nbww, ncgp, lrf, crg, calculation_assignment%ww%max


!$omp single
nbww_pair = 0
!$omp end single
!$omp barrier
rcut2 = Rcww*Rcww
RcLRF2 = RcLRF*RcLRF

#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (calculation_assignment%ww%end - calculation_assignment%ww%start + 1)/threads_num
remainder = MOD(calculation_assignment%ww%end - calculation_assignment%ww%start + 1, threads_num)
mp_start = thread_id * quotient + calculation_assignment%ww%start + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif

iwloop: do iw = mp_start,mp_end
#else
iwloop: do iw = calculation_assignment%ww%start, calculation_assignment%ww%end
#endif
        is  = nat_solute + solv_atom*iw-(solv_atom-1)
        is3 = 3*is-3
        ig = iwhich_cgp(is)
#ifdef USE_GRID
        iagrid = ww_igrid(iw)
        gridnum = 0
igridloop:      do igrid = 1, ww_ndim
jgridloop:       do jgrid = 1, ww_ndim
kgridloop:        do kgrid = 1, ww_ndim
                        gridnum = gridnum + 1
                        if (.not. (grid_ww_int(iagrid,igrid,jgrid,kgrid))) then
ggloop:                  do jw = 1, grid_ww_ngrp(gridnum)
                          jwr = grid_ww_grp(gridnum,jw)
                          ja = nat_solute + solv_atom*jwr-(solv_atom-1)
                          jg = iwhich_cgp(ja)
                          call lrf_update(ig,jg)
                          call lrf_update(jg,ig)
                         end do ggloop
                        else
jwloop: do jw = 1, grid_ww_ngrp(gridnum)
                jwr = grid_ww_grp(gridnum,jw)
#else
jwloop: do jwr = 1, nwat
#endif
                ja = nat_solute + solv_atom*jwr-(solv_atom-1)
                j3 = 3*ja-3

! count each w-w pair once only
                if ( ((iw .gt. jwr) .and. (mod(iw+jwr,2) .eq. 0)) .or. &
                        ((iw .lt. jwr) .and. (mod(iw+jwr,2) .eq. 1)) .or. &
                        (iw .eq. jwr) ) &
                        cycle jwloop

                jg = iwhich_cgp(ja)
                dx = x(is3+1) -x(j3+1)
                dy = x(is3+2) -x(j3+2)
                dz = x(is3+3) -x(j3+3)
                dx = dx - boxlength(1)*nint( dx*inv_boxl(1) )
                dy = dy - boxlength(2)*nint( dy*inv_boxl(2) )
                dz = dz - boxlength(3)*nint( dz*inv_boxl(3) )
                r2 = dx*dx + dy*dy + dz*dz

                if ( r2 .le. rcut2 ) then
! inside cutoff: add the pair
! check if there is enough space
!$omp critical
                        if (nbww_pair .ge. (calculation_assignment%ww%max - solv_atom**2)) call reallocate_nonbondlist_ww
laloop:                 do la = 1, solv_atom
kaloop:                         do ka = 1, solv_atom
                                        nbww_pair = nbww_pair + 1
                                        nbww(nbww_pair)%i = nat_solute+solv_atom*iw-(solv_atom-la)
                                        nbww(nbww_pair)%j = nat_solute+solv_atom*jwr-(solv_atom-ka)
                                end do kaloop
                        end do laloop
!$omp end critical

                elseif((r2 .le. RcLRF2) .or. (RcLRF .eq. -one)) then
! outside ww-cutoff but inside LRF cut-off: use LRF
!iw interaction     
                        call lrf_update(ig,jg)	
!jw interaction
                        call lrf_update(jg,ig)
                end if
        end do jwloop
#ifdef USE_GRID
end if ! interaction
end do kgridloop
end do jgridloop
end do igridloop
#endif
end do iwloop
#if defined (PROFILING)
profile(2)%time = profile(2)%time + rtime() - start_loop_time
#endif
end subroutine nbwwlist_box_lrf
!---------------------------------------------------------------------------
subroutine nbmonitorlist
!precalculate interactions for later use here, too
!old stuff was not usefull any more
! local variables
integer         :: i,j,ig,jg,ia,ja,i3,j3,nl,istate,LJ_code,maxingroups,par, atomnri
integer         :: max_int
integer,allocatable :: num_int(:)
integer         :: grpi,grpj,atomi,atomj,qq_pair,aLJ,bLJ


if (monitor_group_pairs == 0) return

!check the size of the largest group
maxingroups=maxval(monitor_atom_group(:)%n)
max_int = sum(monitor_atom_group(:)%n)

allocate(monitor_group_int(max_int,nstates))
monitor_group_int(:,:)%score = zero
monitor_group_int(:,:)%soft  = .false.
monitor_group_int(:,:)%i     = -1
monitor_group_int(:,:)%j     = -1
monitor_group_int(:,:)%elec  = -1E35_prec
monitor_group_int(:,:)%vdWA  = -1E35_prec
monitor_group_int(:,:)%vdWB  = -1E35_prec
allocate(num_int(nstates))
num_int(:) = 0

do par=1,monitor_group_pairs
        monitor_group_pair(par)%lstart(:) = num_int(:) + 1
        grpi=monitor_group_pair(par)%i      
        grpj=monitor_group_pair(par)%j 
        do i=1,monitor_atom_group(grpi)%n
                atomi=monitor_atom_group(grpi)%atom(i)
                do j=1,monitor_atom_group(grpj)%n
                        atomj=monitor_atom_group(grpj)%atom(j)    
!assert that atoms are different
                        if(atomi == atomj) then 
                                call die('two paired monitor atom groups contain the same atom')
                        end if
!check which list we need to check
                        if((iqatom(atomi).ne.0) .and. (iqatom(atomj).ne.0)) then
!both q_atoms, are on qq list
!check if those bastards are allowed to interact
!we now do the real checking to see if the atoms actually interact
!to make sure we get the interaction energy right 
                                do istate = 1 , nstates
                                if (.not.(qq_precomp(iqatom(atomi),iqatom(atomj),istate)%set)) cycle
                                        num_int(istate) = num_int(istate) + 1
                                        monitor_group_int(num_int(istate),istate)%i= atomi
                                        monitor_group_int(num_int(istate),istate)%j= atomj
                                        monitor_group_int(num_int(istate),istate)%elec = &
                                                qq_precomp(iqatom(atomi),iqatom(atomj),istate)%elec
                                        monitor_group_int(num_int(istate),istate)%vdWA = &
                                                qq_precomp(iqatom(atomi),iqatom(atomj),istate)%vdWA
                                        monitor_group_int(num_int(istate),istate)%vdWB = &
                                                qq_precomp(iqatom(atomi),iqatom(atomj),istate)%vdWB
! for q-q atoms we also need the info on soft pairs and soft core stuff
! so the group_int stuff needs to be able to store this, too
                                        monitor_group_int(num_int(istate),istate)%soft = &
                                                qq_precomp(iqatom(atomi),iqatom(atomj),istate)%soft
                                        monitor_group_int(num_int(istate),istate)%score = &
                                                qq_precomp(iqatom(atomi),iqatom(atomj),istate)%score
                                end do

                        else if(iqatom(atomi).ne.0) then
!only atom i is Q atom, atom j will be on q-p list then
!so we check the q-p precompute table for them, to also account for qqp atoms
                                do istate = 1 , nstates
                                if(.not.(qp_precomp(iqatom(atomi),atomj,istate)%set)) cycle
                                        num_int(istate) = num_int(istate) + 1
                                        monitor_group_int(num_int(istate),istate)%i= atomi
                                        monitor_group_int(num_int(istate),istate)%j= atomj
                                        monitor_group_int(num_int(istate),istate)%elec = &
                                                qp_precomp(iqatom(atomi),atomj,istate)%elec
                                        monitor_group_int(num_int(istate),istate)%vdWA = &
                                                qp_precomp(iqatom(atomi),atomj,istate)%vdWA
                                        monitor_group_int(num_int(istate),istate)%vdWB = &
                                                qp_precomp(iqatom(atomi),atomj,istate)%vdWB
! for q-q atoms we also need the info on soft pairs and soft core stuff
! so the group_int stuff needs to be able to store this, too
                                        monitor_group_int(num_int(istate),istate)%score = &
                                                qp_precomp(iqatom(atomi),atomj,istate)%score
                                end do
!other way around, j is Q atom, same code but inverted
                        else if(iqatom(atomj).ne.0) then
                                do istate = 1 , nstates
                                if(.not.(qp_precomp(iqatom(atomj),atomi,istate)%set)) cycle
                                        num_int(istate) = num_int(istate) + 1
                                        monitor_group_int(num_int(istate),istate)%i= atomi
                                        monitor_group_int(num_int(istate),istate)%j= atomj
                                        monitor_group_int(num_int(istate),istate)%elec = &
                                                qp_precomp(iqatom(atomj),atomi,istate)%elec
                                        monitor_group_int(num_int(istate),istate)%vdWA = &
                                                qp_precomp(iqatom(atomj),atomi,istate)%vdWA
                                        monitor_group_int(num_int(istate),istate)%vdWB = &
                                                qp_precomp(iqatom(atomj),atomi,istate)%vdWB
! for q-q atoms we also need the info on soft pairs and soft core stuff
! so the group_int stuff needs to be able to store this, too
                                        monitor_group_int(num_int(istate),istate)%score = &
                                                qp_precomp(iqatom(atomj),atomi,istate)%score
                                end do
! and if neither is a qatom, we use the pp-list
                        else
                                if(.not.(pp_precomp(i,j)%set)) cycle
                                num_int(:) = num_int(:) + 1
                                monitor_group_int(num_int(:),:)%i= atomi
                                monitor_group_int(num_int(:),:)%j= atomj
                                monitor_group_int(num_int(:),:)%elec = &
                                        pp_precomp(i,j)%elec
                                monitor_group_int(num_int(:),:)%vdWA = &
                                        pp_precomp(i,j)%vdWA
                                monitor_group_int(num_int(:),:)%vdWB = &
                                        pp_precomp(i,j)%vdWB
                        end if
                end do
        end do
        monitor_group_pair(par)%lend(:) = num_int(:)
end do
! cleaning up again
deallocate(num_int)                                
end subroutine nbmonitorlist

!---------------------------------------------------------------------------------------------------------		   

subroutine nonbond_monitor
!monitor nonbonded energies between selected groups of atoms
!rewritten tu use precomputed interaction energies

real(kind=prec)  :: dx1,dx2,dx3,r,r2,r6,r12,r6_hc
integer	 :: i,j,istate,par
integer	 :: atomi,atomj,cstart,cend,calcgroup
real(kind=prec)  :: V_a,V_b,Vel,Vvdw,Vwel,Vwvdw,Vwsum


do par=1,monitor_group_pairs
        monitor_group_pair(par)%Vel(:)=zero
        monitor_group_pair(par)%Vlj(:)=zero
        monitor_group_pair(par)%Vwel = zero
        monitor_group_pair(par)%Vwlj = zero
        monitor_group_pair(par)%Vwsum= zero  

        do istate = 1, nstates   
! loop over all our nice pairs for each state
! defined under lstart and lend
! that then point to the right number in monitor_int
        cstart = monitor_group_pair(par)%lstart(istate)
        cend   = monitor_group_pair(par)%lend(istate)

                do calcgroup = cstart, cend
                        atomi = monitor_group_int(calcgroup,istate)%i
                        atomj = monitor_group_int(calcgroup,istate)%j
        
                        dx1  = x(3*atomi-2)-x(3*atomj-2)      ! calculate the distance
                        dx2  = x(3*atomi-1)-x(3*atomj-1)
                        dx3  = x(3*atomi)-x(3*atomj)
                        if (use_PBC) then
                        	dx1 = dx1 - boxlength(1)*nint( dx1*inv_boxl(1) )
                        	dx2 = dx2 - boxlength(2)*nint( dx2*inv_boxl(2) )
                        	dx3 = dx3 - boxlength(3)*nint( dx3*inv_boxl(3) )
                        end if

                        r2   = dx1*dx1 + dx2*dx2 + dx3*dx3
                        r6_hc = r2*r2*r2  !for softcore
                        r6   = r6_hc + monitor_group_int(calcgroup,istate)%score !softcore, default = zero
                        r6   = one/r6
                        r2   = one/r2
                        r    = sqrt ( r2 )
                        r12  = r6*r6
        
                        if (monitor_group_int(calcgroup,istate)%soft) then
                                V_b = monitor_group_int(calcgroup,istate)%vdWB
                                V_a = monitor_group_int(calcgroup,istate)%vdWA*exp(-V_b/r)
                        else
                                V_a  = monitor_group_int(calcgroup,istate)%vdWA *r12
                                V_b  = monitor_group_int(calcgroup,istate)%vdWB *r6
                        end if
                        Vel  = monitor_group_int(calcgroup,istate)%elec * r
                        Vvdw = V_a - V_b
                        monitor_group_pair(par)%Vel(istate) = &
                                monitor_group_pair(par)%Vel(istate)+Vel
                        monitor_group_pair(par)%Vlj(istate) = &
                                monitor_group_pair(par)%Vlj(istate)+Vvdw
                end do ! groups
        end do ! nstates
        !calc lambda-weighted sum
        monitor_group_pair(par)%Vwel=dot_product(monitor_group_pair(par)%Vel(1:nstates),EQ(1:nstates)%lambda)
        monitor_group_pair(par)%Vwlj=dot_product(monitor_group_pair(par)%Vlj(1:nstates),EQ(1:nstates)%lambda)
        monitor_group_pair(par)%Vwsum= monitor_group_pair(par)%Vwlj+monitor_group_pair(par)%Vwel
end do !par

end subroutine nonbond_monitor

!-------------------------------------------------------------------------
! nonbond routines are now consolidated to just sphere/box versions
! because combination rules are applied earlier

subroutine nonbond_pp
! local variables
integer						:: ip,i,j,i3,j3
real(kind=prec)						:: dx1,dx2,dx3,r2,r,r6,r12
real(kind=prec)						:: Vel,V_a,V_b,dv

#ifdef _OPENMP
integer :: quotient, remainder
real(kind=prec) :: Vel_omp, Vvdw_omp
#endif

! global variables used:
! x,d,E


#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (nbpp_pair)/threads_num
remainder = MOD(nbpp_pair, threads_num)
mp_start = thread_id * quotient + 1 + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif
Vel_omp = zero
Vvdw_omp = zero
do ip = mp_start,mp_end
#else
do ip = 1, nbpp_pair 
#endif

! do this for every pair and use the compiler to unroll the loop

! get atom indicies
i = nbpp(ip)%i
j = nbpp(ip)%j

i3   = i*3-3
j3   = j*3-3
! calculate dx and r
dx1  = x(j3+1) - x(i3+1)
dx2  = x(j3+2) - x(i3+2)
dx3  = x(j3+3) - x(i3+3)

r2    = dx1*dx1 + dx2*dx2 + dx3*dx3
r2    = one/r2
r     = sqrt ( r2 )
r6    = r2*r2*r2
r12   = r6*r6

! calculate Vel and dv
Vel  = nbpp(ip)%elec * r
V_a  = nbpp(ip)%vdWA *r12 
V_b  = nbpp(ip)%vdWB *r6
dv   = r2*( -Vel -12.0_prec*V_a +6.0_prec*V_b )
! update d
d(i3+1) = d(i3+1) - dv*dx1
d(i3+2) = d(i3+2) - dv*dx2
d(i3+3) = d(i3+3) - dv*dx3
d(j3+1) = d(j3+1) + dv*dx1
d(j3+2) = d(j3+2) + dv*dx2
d(j3+3) = d(j3+3) + dv*dx3

! update energies
#ifdef _OPENMP
Vel_omp  = Vel_omp + Vel
Vvdw_omp = Vvdw_omp +  V_a - V_b
#else
E%pp%el  = E%pp%el + Vel 
E%pp%vdw = E%pp%vdw + V_a - V_b
#endif
end do

#ifdef _OPENMP
!$omp atomic update
E%pp%el = E%pp%el + Vel_omp
!$omp atomic update
E%pp%vdw = E%pp%vdw + Vvdw_omp
#endif


end subroutine nonbond_pp

!------------------------------------------------------------------------
subroutine nonbond_pp_box
  ! local variables
  integer						:: ip, i, j, ga, gb, group
  real(kind=prec)						:: dx1,dx2,dx3,r2,r,r6
  real(kind=prec)						:: Vel,V_a,V_b,dv

#ifdef _OPENMP
integer :: quotient, remainder
real(kind=prec) :: Vel_omp, Vvdw_omp
#endif
! global variables used:
! x,d,E

#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (nbpp_cgp_pair)/threads_num
remainder = MOD(nbpp_cgp_pair, threads_num)
mp_start = thread_id * quotient + 1 + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif

do group = mp_start,mp_end
#else
do group = 1, nbpp_cgp_pair
#endif
ga = nbpp_cgp(group)%i !atom index for the two switching atoms
gb = nbpp_cgp(group)%j

!the distance between the two switching atoms
dx1 = x(3*gb-2) - x(3*ga-2)
dx2 = x(3*gb-1) - x(3*ga-1)
dx3 = x(3*gb  ) - x(3*ga  )

nbpp_cgp(group)%x = boxlength(1)*nint( dx1*inv_boxl(1) )
nbpp_cgp(group)%y = boxlength(2)*nint( dx2*inv_boxl(2) )	
nbpp_cgp(group)%z = boxlength(3)*nint( dx3*inv_boxl(3) )

end do
!$omp barrier
#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (nbpp_pair)/threads_num
remainder = MOD(nbpp_pair, threads_num)
mp_start = thread_id * quotient + 1 + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif
Vel_omp = zero
Vvdw_omp = zero
do ip = mp_start,mp_end
#else
do ip = 1, nbpp_pair !- 1, 2
#endif
! for every pair, not every second pair
! get indicies
ga = nbpp(ip)%cgp_pair
i  = nbpp(ip)%i
j  = nbpp(ip)%j

! calculate dx, r and r2
dx1  = x(j*3-2) - x(i*3-2)
dx2  = x(j*3-1) - x(i*3-1)
dx3  = x(j*3-0) - x(i*3-0)
dx1 = dx1 - nbpp_cgp(ga)%x         
dx2 = dx2 - nbpp_cgp(ga)%y    
dx3 = dx3 - nbpp_cgp(ga)%z  

r2   = one/(dx1*dx1 + dx2*dx2 + dx3*dx3)
r = sqrt(r2)
r6   = r2*r2*r2

Vel  = nbpp(ip)%elec * r
V_a  = nbpp(ip)%vdWA *r6*r6 
V_b  = nbpp(ip)%vdWB *r6
dv   = r2*( -Vel -12.0_prec*V_a +6.0_prec*V_b )

! update d
d(i*3-2) = d(i*3-2) - dv*dx1
d(i*3-1) = d(i*3-1) - dv*dx2
d(i*3-0) = d(i*3-0) - dv*dx3
d(j*3-2) = d(j*3-2) + dv*dx1
d(j*3-1) = d(j*3-1) + dv*dx2
d(j*3-0) = d(j*3-0) + dv*dx3

! update energies
#ifdef _OPENMP
Vel_omp = Vel_omp + Vel
Vvdw_omp = Vvdw_omp + V_a - V_b
#else
E%pp%el  = E%pp%el + Vel
E%pp%vdw = E%pp%vdw + V_a - V_b
#endif
end do

#ifdef _OPENMP
!$omp atomic update
E%pp%el  = E%pp%el + Vel_omp
!$omp atomic update
E%pp%vdw = E%pp%vdw + Vvdw_omp
#endif


end subroutine nonbond_pp_box
!----------------------------------------------------------------------

subroutine nonbond_pw
! local variables
integer						:: ip,i,j,i3,j3
real(kind=prec)						:: dx1,dx2,dx3,r2,r,r6,r12
real(kind=prec)						:: Vel,V_a,V_b,dv
#ifdef _OPENMP
integer :: quotient, remainder
real(kind=prec) :: Vel_omp, Vvdw_omp
#endif
#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (nbpw_pair)/threads_num
remainder = MOD(nbpw_pair, threads_num)
mp_start = thread_id * quotient + 1 + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif
Vel_omp = zero
Vvdw_omp = zero
do ip = mp_start,mp_end
#else
! global variables used:
! x, d, E

do ip = 1, nbpw_pair
#endif
! for every assigned pair:
i    = nbpw(ip)%i
j    = nbpw(ip)%j
i3   = i*3-3
j3   = j*3-3

! calculate dx and r
dx1  = x(j3+1) - x(i3+1)
dx2  = x(j3+2) - x(i3+2)
dx3  = x(j3+3) - x(i3+3)

r2   = dx1*dx1 + dx2*dx2 + dx3*dx3
r2   = one/r2
r    = sqrt ( r2 ) 
r6   = r2*r2*r2
r12  = r6*r6

! calculate Vel and dv
Vel  = nbpw(ip)%elec *r
V_a  = nbpw(ip)%vdWA *r12 
V_b  = nbpw(ip)%vdWB *r6
dv   = r2*( -Vel -12.0_prec*V_a +6.0_prec*V_b )

! update forces
d(i3+1) = d(i3+1) - dv*dx1
d(i3+2) = d(i3+2) - dv*dx2
d(i3+3) = d(i3+3) - dv*dx3
d(j3+1) = d(j3+1) + dv*dx1
d(j3+2) = d(j3+2) + dv*dx2
d(j3+3) = d(j3+3) + dv*dx3

! update energies
#ifdef _OPENMP
Vel_omp =Vel_omp + Vel
Vvdw_omp = Vvdw_omp + V_a - V_b
#else
E%pw%el  = E%pw%el + Vel       
E%pw%vdw = E%pw%vdw + V_a - V_b
#endif
end do
#ifdef _OPENMP
!$omp atomic update
E%pw%el = E%pw%el + Vel_omp
!$omp atomic update
E%pw%vdw = E%pw%vdw + Vvdw_omp
#endif
end subroutine nonbond_pw

!-----------------------------------------------------------------------

subroutine nonbond_pw_box
  ! local variables
  integer						:: ip,i,j,i3,j3
  real(kind=prec)						:: dx1,dx2,dx3,r2,r,r6,r12
  real(kind=prec)						:: Vel,V_a,V_b,dv
  integer						:: group, ga, gb
#ifdef _OPENMP
integer :: quotient, remainder
real(kind=prec) :: Vel_omp, Vvdw_omp
#endif
  ! global variables used:
  !  iac, crg, iaclib, x, d, E

#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (nbpw_cgp_pair)/threads_num
remainder = MOD(nbpw_cgp_pair, threads_num)
mp_start = thread_id * quotient + 1 + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
mp_end = mp_end + 1
endif

do group = mp_start,mp_end
#else
!compute the peridocal shift for every charge group pair
do group = 1, nbpw_cgp_pair
#endif
ga = nbpw_cgp(group)%i !atom index for solute switching atom
gb = nbpw_cgp(group)%j  !atom index for the solvent switching atom

!the distance between the two switching atoms
dx1 = x(3*gb-2) - x(3*ga-2)
dx2 = x(3*gb-1) - x(3*ga-1)
dx3 = x(3*gb  ) - x(3*ga  )
nbpw_cgp(group)%x = boxlength(1)*nint( dx1*inv_boxl(1) )
nbpw_cgp(group)%y = boxlength(2)*nint( dx2*inv_boxl(2) )	
nbpw_cgp(group)%z = boxlength(3)*nint( dx3*inv_boxl(3) )

end do
!$omp barrier
#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (nbpw_pair)/threads_num
remainder = MOD(nbpw_pair, threads_num)
mp_start = thread_id * quotient + 1 + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
mp_end = mp_end + 1
endif
Vel_omp = zero
Vvdw_omp = zero
do ip = mp_start,mp_end
#else
do ip = 1, nbpw_pair
#endif
! for every assigned pair:
i    = nbpw(ip)%i !solute atom 
j    = nbpw(ip)%j !solvent atom
group = nbpw(ip)%cgp_pair
i3   = i*3-3
j3   = j*3-3

! calculate dx and r
dx1  = x(j3+1) - x(i3+1)
dx2  = x(j3+2) - x(i3+2)
dx3  = x(j3+3) - x(i3+3)
dx1 = dx1 - nbpw_cgp(group)%x  
dx2 = dx2 - nbpw_cgp(group)%y  
dx3 = dx3 - nbpw_cgp(group)%z  


r2   = dx1*dx1 + dx2*dx2 + dx3*dx3
r2   = one/r2
r    = sqrt ( r2 ) 
r6   = r2*r2*r2
r12  = r6*r6

! calculate Vel and dv
Vel  = nbpw(ip)%elec *r
V_a  = nbpw(ip)%vdWA *r12 
V_b  = nbpw(ip)%vdWB *r6
dv   = r2*( -Vel -12.0_prec*V_a +6.0_prec*V_b )
! update forces
d(i3+1) = d(i3+1) - dv*dx1
d(i3+2) = d(i3+2) - dv*dx2
d(i3+3) = d(i3+3) - dv*dx3
d(j3+1) = d(j3+1) + dv*dx1
d(j3+2) = d(j3+2) + dv*dx2
d(j3+3) = d(j3+3) + dv*dx3

! update energies
#ifdef _OPENMP
Vel_omp =Vel_omp + Vel
Vvdw_omp = Vvdw_omp + V_a - V_b
#else
E%pw%el  = E%pw%el + Vel       
E%pw%vdw = E%pw%vdw + V_a - V_b
#endif
end do
#ifdef _OPENMP
!$omp atomic update
E%pw%el = E%pw%el + Vel_omp
!$omp atomic update
E%pw%vdw = E%pw%vdw + Vvdw_omp
#endif

end subroutine nonbond_pw_box

!-----------------------------------------------------------------------

subroutine nonbond_qq
! local variables
integer						:: istate,jj
integer						:: ip,iq,jq,i,j,k,i3,j3
real(kind=prec)						:: qi,qj,dx1,dx2,dx3,r2,r,r6,r12,r6_hc
real(kind=prec)						:: Vel,V_a,V_b,dv
#ifdef _OPENMP
integer :: quotient, remainder
#endif
#ifdef _OPENMP
qomp_elec  = zero
qomp_vdw   = zero
#endif
do istate = 1, nstates
! for every state:
#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (nbqq_pair(istate))/threads_num
remainder = MOD(nbqq_pair(istate), threads_num)
mp_start = thread_id * quotient + 1 + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
mp_end = mp_end + 1
endif

do ip = mp_start,mp_end
#else
do ip = 1, nbqq_pair(istate)
#endif
  ! for every pair:
iq   = nbqq(ip,istate)%iq
i    = iqseq(iq)
jq   = nbqq(ip,istate)%jq
j    = iqseq(jq)
i3   = i*3-3
j3   = j*3-3

! calculate dx and r
dx1  = x(j3+1) - x(i3+1)
dx2  = x(j3+2) - x(i3+2)
dx3  = x(j3+3) - x(i3+3)

r2    = dx1*dx1 + dx2*dx2 + dx3*dx3
r6    = r2*r2*r2
r6_hc = r6
r2    = one/r2
r     = sqrt ( r2 )
r6    = r6_hc + nbqq(ip,istate)%score !softcore
r6    = one/r6
r12   = r6*r6

! calculate Vel, V_a, V_b and dv
Vel  = nbqq(ip,istate)%elec *r
if (nbqq(ip,istate)%soft) then
V_b = zero
V_a = nbqq(ip,istate)%vdWA*exp(-nbqq(ip,istate)%vdWB/r)
dv  = r2*( -Vel -V_b*V_a/r )*EQ(istate)%lambda
else
V_a  = nbqq(ip,istate)%vdWA *r12 
V_b  = nbqq(ip,istate)%vdWB *r6
dv  = r2*( -Vel -(12.0_prec*V_a -6.0_prec*V_b)*r6*r6_hc )*EQ(istate)%lambda
endif

! update forces
d(i3+1) = d(i3+1) - dv*dx1
d(i3+2) = d(i3+2) - dv*dx2
d(i3+3) = d(i3+3) - dv*dx3
d(j3+1) = d(j3+1) + dv*dx1
d(j3+2) = d(j3+2) + dv*dx2
d(j3+3) = d(j3+3) + dv*dx3
! update energies
#ifdef _OPENMP
qomp_elec(istate,1)  = qomp_elec(istate,1) + Vel
qomp_vdw(istate,1)   = qomp_vdw(istate,1) + V_a - V_b
#else
EQ(istate)%qq(1)%el  = EQ(istate)%qq(1)%el + Vel
EQ(istate)%qq(1)%vdw = EQ(istate)%qq(1)%vdw + V_a - V_b
#endif
if (use_excluded_groups) then
do jj = 2, ene_header%arrays
if (.not.((ST_gc(ene_header%gcnum(jj))%gcmask%mask(i)).or. &
ST_gc(ene_header%gcnum(jj))%gcmask%mask(j))) then
#ifdef _OPENMP
qomp_elec(istate,jj)  = qomp_elec(istate,jj) + Vel
qomp_vdw(istate,jj)   = qomp_vdw(istate,jj) + V_a - V_b
#else
EQ(istate)%qq(jj)%el  = EQ(istate)%qq(jj)%el + Vel
EQ(istate)%qq(jj)%vdw = EQ(istate)%qq(jj)%vdw + V_a - V_b
#endif
else
select case (ST_gc(ene_header%gcnum(jj))%caltype)
case (FULL)
!do nothing!
case (ELECTRO)
#ifdef _OPENMP
qomp_vdw(istate,jj)   = qomp_vdw(istate,jj) + V_a - V_b
#else
EQ(istate)%qq(jj)%vdw = EQ(istate)%qq(jj)%vdw + V_a -V_b
#endif
case (VDW)
#ifdef _OPENMP
qomp_elec(istate,jj)  = qomp_elec(istate,jj) + Vel
#else
EQ(istate)%qq(jj)%el  = EQ(istate)%qq(jj)%el + Vel
#endif
end select
end if
end do
end if
end do ! ip

end do ! istate
#ifdef _OPENMP
do istate=1, nstates
do jj=1, ene_header%arrays
!$omp atomic update
EQ(istate)%qq(jj)%el  = EQ(istate)%qq(jj)%el + qomp_elec(istate,jj)
!$omp atomic update
EQ(istate)%qq(jj)%vdw = EQ(istate)%qq(jj)%vdw + qomp_vdw(istate,jj)
end do
end do
#endif
end subroutine nonbond_qq

!-----------------------------------------------------------------------
subroutine nonbond_qqp
! local variables
integer                                         :: istate,jj
integer                                         :: ip,iq,jq,i,j,k,i3,j3
real(kind=prec)                                         :: qi,qj,dx1,dx2,dx3,r2,r,r6,r12,r6_hc
real(kind=prec)                                         :: Vel,V_a,V_b,dv
#ifdef _OPENMP
integer :: quotient, remainder
#endif
#ifdef _OPENMP
qomp_elec  = zero
qomp_vdw   = zero
#endif
do istate = 1, nstates
! for every state:
#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (nbqqp_pair(istate))/threads_num
remainder = MOD(nbqqp_pair(istate), threads_num)
mp_start = thread_id * quotient + 1 + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
mp_end = mp_end + 1
endif

do ip = mp_start,mp_end
#else
do ip = 1, nbqqp_pair(istate)
#endif
  ! for every pair:
iq   = nbqqp(ip,istate)%iq
i    = iqseq(iq)
j    = nbqqp(ip,istate)%jq
i3   = i*3-3
j3   = j*3-3
! calculate dx and r
dx1  = x(j3+1) - x(i3+1)
dx2  = x(j3+2) - x(i3+2)
dx3  = x(j3+3) - x(i3+3)

r2    = dx1*dx1 + dx2*dx2 + dx3*dx3
r6    = r2*r2*r2
r6_hc = r6
r2    = one/r2
r     = sqrt ( r2 )
r6    = r6_hc + nbqqp(ip,istate)%score !softcore
r6    = one/r6
r12   = r6*r6

! calculate Vel, V_a, V_b and dv
Vel  = nbqqp(ip,istate)%elec *r
V_a  = nbqqp(ip,istate)%vdWA *r12
V_b  = nbqqp(ip,istate)%vdWB *r6
dv  = r2*( -Vel -(12.0_prec*V_a -6.0_prec*V_b)*r6*r6_hc )*EQ(istate)%lambda

! update forces
d(i3+1) = d(i3+1) - dv*dx1
d(i3+2) = d(i3+2) - dv*dx2
d(i3+3) = d(i3+3) - dv*dx3
d(j3+1) = d(j3+1) + dv*dx1
d(j3+2) = d(j3+2) + dv*dx2
d(j3+3) = d(j3+3) + dv*dx3
! update energies
#ifdef _OPENMP
qomp_elec(istate,1) = qomp_elec(istate,1) + Vel
qomp_vdw(istate,1) = qomp_vdw(istate,1) + V_a - V_b
#else
EQ(istate)%qp(1)%el  = EQ(istate)%qp(1)%el + Vel
EQ(istate)%qp(1)%vdw = EQ(istate)%qp(1)%vdw + V_a - V_b
#endif
if (use_excluded_groups) then
do jj = 2, ene_header%arrays
if (.not.((ST_gc(ene_header%gcnum(jj))%gcmask%mask(i)).or. &
ST_gc(ene_header%gcnum(jj))%gcmask%mask(j))) then
#ifdef _OPENMP
qomp_elec(istate,jj) = qomp_elec(istate,jj) + Vel
qomp_vdw(istate,jj)  = qomp_vdw(istate,jj) + V_a - V_b
#else
EQ(istate)%qp(jj)%el  = EQ(istate)%qp(jj)%el + Vel
EQ(istate)%qp(jj)%vdw = EQ(istate)%qp(jj)%vdw + V_a - V_b
#endif
else
select case (ST_gc(ene_header%gcnum(jj))%caltype)
case (FULL)
!do nothing!
case (ELECTRO)
#ifdef _OPENMP
qomp_vdw(istate,jj) = qomp_vdw(istate,jj) + V_a - V_b
#else
EQ(istate)%qp(jj)%vdw = EQ(istate)%qp(jj)%vdw + V_a -V_b
#endif
case (VDW)
#ifdef _OPENMP
qomp_elec(istate,jj) = qomp_elec(istate,jj) + Vel
#else
EQ(istate)%qp(jj)%el  = EQ(istate)%qp(jj)%el + Vel
#endif
end select
end if
end do
end if
end do ! ip

end do ! istate
#ifdef _OPENMP
do istate=1, nstates
do jj=1, ene_header%arrays
!$omp atomic update
EQ(istate)%qp(jj)%el  = EQ(istate)%qp(jj)%el + qomp_elec(istate,jj)
!$omp atomic update
EQ(istate)%qp(jj)%vdw = EQ(istate)%qp(jj)%vdw + qomp_vdw(istate,jj)
end do
end do
#endif
end subroutine nonbond_qqp

!-----------------------------------------------------------------------

subroutine nonbond_qp
! local variables
integer						:: ip,iq,i,j,i3,j3
integer						:: istate,jj
real(kind=prec)						:: dx1,dx2,dx3,r2,r,r6,r6_hc,r12
real(kind=prec)						:: Vel,V_a,V_b,dv
#ifdef _OPENMP
integer :: quotient, remainder
#endif
! global variables used:
!  iqseq, iac, crg, x, nstates, qvdw_flag, iaclib, qiac, qavdw, qbvdw, qcrg, el14_scale, EQ, d, nat_solute
#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (nbqp_pair)/threads_num
remainder = MOD(nbqp_pair, threads_num)
mp_start = thread_id * quotient + 1 + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif

qomp_elec  = zero
qomp_vdw   = zero
do ip = mp_start,mp_end
#else
do ip = 1, nbqp_pair
#endif
! for every assigned q-s pair:

! init state-invariant variables:
iq   = nbqp(ip,1)%i
i    = iqseq(iq)
j    = nbqp(ip,1)%j
i3   = i*3-3
j3   = j*3-3


dx1  = x(j3+1) - x(i3+1)
dx2  = x(j3+2) - x(i3+2)
dx3  = x(j3+3) - x(i3+3)


r2    = dx1*dx1 + dx2*dx2 + dx3*dx3
r6    = r2*r2*r2
r6_hc = r6
r2    = one/r2
r     = sqrt ( r2 )



do istate = 1, nstates
r6    = r6_hc + nbqp(ip,istate)%score !softcore
r6    = one/r6
r12   = r6*r6
! calculate qi, Vel, V_a, V_b and dv
Vel  = nbqp(ip,istate)%elec *r
V_a  = nbqp(ip,istate)%vdWA *r12
V_b  = nbqp(ip,istate)%vdWB *r6
dv   = r2*( -Vel -(12.0_prec*V_a -6.0_prec*V_b)*r6*r6_hc )*EQ(istate)%lambda   !softcore r6*r6_hc is (r^6/(r^6+alpha))

! update forces
d(i3+1) = d(i3+1) - dv*dx1
d(i3+2) = d(i3+2) - dv*dx2
d(i3+3) = d(i3+3) - dv*dx3
d(j3+1) = d(j3+1) + dv*dx1
d(j3+2) = d(j3+2) + dv*dx2
d(j3+3) = d(j3+3) + dv*dx3

  ! update q-protein or q-water energies
#ifdef _OPENMP
qomp_elec(istate,1) = qomp_elec(istate,1) + Vel
qomp_vdw(istate,1)  = qomp_vdw(istate,1) + V_a - V_b
#else
EQ(istate)%qp(1)%el  = EQ(istate)%qp(1)%el + Vel
EQ(istate)%qp(1)%vdw = EQ(istate)%qp(1)%vdw + V_a - V_b
#endif
if (use_excluded_groups) then
do jj = 2, ene_header%arrays
if (.not.((ST_gc(ene_header%gcnum(jj))%gcmask%mask(i)).or. &
        ST_gc(ene_header%gcnum(jj))%gcmask%mask(j))) then
#ifdef _OPENMP
qomp_elec(istate,jj) = qomp_elec(istate,jj) + Vel
qomp_vdw(istate,jj)  = qomp_vdw(istate,jj) + V_a - V_b
#else
EQ(istate)%qp(jj)%el  = EQ(istate)%qp(jj)%el + Vel
EQ(istate)%qp(jj)%vdw = EQ(istate)%qp(jj)%vdw + V_a - V_b
#endif
else
select case (ST_gc(ene_header%gcnum(jj))%caltype)
case (FULL)
!do nothing!
case (ELECTRO)
#ifdef _OPENMP
qomp_vdw(istate,jj)  = qomp_vdw(istate,jj) + V_a - V_b
#else
EQ(istate)%qp(jj)%vdw = EQ(istate)%qp(jj)%vdw + V_a -V_b
#endif
case (VDW)
#ifdef _OPENMP
qomp_elec(istate,jj) = qomp_elec(istate,jj) + Vel
#else
EQ(istate)%qp(jj)%el  = EQ(istate)%qp(jj)%el + Vel
#endif
end select
end if
end do 
end if
end do ! istate

end do

#ifdef _OPENMP
do istate=1, nstates
do jj=1, ene_header%arrays
!$omp atomic update
EQ(istate)%qp(jj)%el  = EQ(istate)%qp(jj)%el + qomp_vdw(istate,jj)
!$omp atomic update
EQ(istate)%qp(jj)%vdw = EQ(istate)%qp(jj)%vdw + qomp_elec(istate,jj)
end do
end do
#endif
end subroutine nonbond_qp
!-----------------------------------------------------------------------

!******PWadded 2001-10-23
subroutine nonbond_qp_box
! local variables
integer						:: ip,iq,i,j,i3,j3
integer						:: istate,jj
real(kind=prec)						:: dx1,dx2,dx3,r2,r,r6,r6_hc,r12
real(kind=prec)						:: Vel,V_a,V_b,dv
integer						:: group, gr, ia
#ifdef _OPENMP
integer :: quotient, remainder
#endif
  ! global variables used:
  !  iqseq, iac, crg, x, nstates, qvdw_flag, iaclib, qiac, qavdw, qbvdw, qcrg, el14_scale, EQ, d, nat_solute

#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (nbqp_cgp_pair)/threads_num
remainder = MOD(nbqp_cgp_pair, threads_num)
mp_start = thread_id * quotient + 1 + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
mp_end = mp_end + 1
endif

do gr = mp_start,mp_end
#else
!compute the peridocal shift for every charge group pair
do gr = 1, nbqp_cgp_pair
#endif
ia = nbqp_cgp(gr)%i !atom index for the atom

!the distance between the two switching atoms
dx1 = x(3*ia-2) - x(3*qswitch-2)
dx2 = x(3*ia-1) - x(3*qswitch-1)
dx3 = x(3*ia  ) - x(3*qswitch  )

nbqp_cgp(gr)%x = boxlength(1)*nint( dx1*inv_boxl(1) )
nbqp_cgp(gr)%y = boxlength(2)*nint( dx2*inv_boxl(2) )	
nbqp_cgp(gr)%z = boxlength(3)*nint( dx3*inv_boxl(3) )
end do
!$omp barrier
#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (nbqp_pair)/threads_num
remainder = MOD(nbqp_pair, threads_num)
mp_start = thread_id * quotient + 1 + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
mp_end = mp_end + 1
endif
qomp_elec = zero
qomp_vdw  = zero
do ip = mp_start,mp_end
#else
do ip = 1, nbqp_pair
#endif
! for every assigned q-s pair:

! init state-invariant variables:
iq   = nbqp(ip,1)%i
i    = iqseq(iq)
j    = nbqp(ip,1)%j
i3   = i*3-3
j3   = j*3-3
group = nbqp(ip,1)%cgp_pair

dx1  = x(j3+1) - x(i3+1)
dx2  = x(j3+2) - x(i3+2)
dx3  = x(j3+3) - x(i3+3)
dx1 = dx1 - nbqp_cgp(group)%x  
dx2 = dx2 - nbqp_cgp(group)%y  
dx3 = dx3 - nbqp_cgp(group)%z  

r2    = dx1*dx1 + dx2*dx2 + dx3*dx3
r6    = r2*r2*r2
r6_hc = r6
r2    = one/r2
r     = sqrt ( r2 )

do istate = 1, nstates
r6    = r6_hc + nbqp(ip,istate)%score !softcore
r6    = one/r6
r12   = r6*r6

! calculate qi, Vel, V_a, V_b and dv
Vel  = nbqp(ip,istate)%elec *r
V_a  = nbqp(ip,istate)%vdWA *r12
V_b  = nbqp(ip,istate)%vdWB *r6
dv   = r2*( -Vel -(12.0_prec*V_a -6.0_prec*V_b)*r6*r6_hc )*EQ(istate)%lambda

! update forces
d(i3+1) = d(i3+1) - dv*dx1
d(i3+2) = d(i3+2) - dv*dx2
d(i3+3) = d(i3+3) - dv*dx3
d(j3+1) = d(j3+1) + dv*dx1
d(j3+2) = d(j3+2) + dv*dx2
d(j3+3) = d(j3+3) + dv*dx3

! update q-protein or q-water energies
#ifdef _OPENMP
qomp_elec(istate,1) = qomp_elec(istate,1) + Vel
qomp_vdw(istate,1)  = qomp_vdw(istate,1) + V_a - V_b
#else
EQ(istate)%qp(1)%el  = EQ(istate)%qp(1)%el + Vel
EQ(istate)%qp(1)%vdw = EQ(istate)%qp(1)%vdw + V_a - V_b
#endif
if (use_excluded_groups) then
do jj = 2, ene_header%arrays
if (.not.((ST_gc(ene_header%gcnum(jj))%gcmask%mask(i)).or. &
        ST_gc(ene_header%gcnum(jj))%gcmask%mask(j))) then
#ifdef _OPENMP
qomp_elec(istate,jj) = qomp_elec(istate,jj) + Vel
qomp_vdw(istate,jj)  = qomp_vdw(istate,jj) + V_a - V_b
#else
EQ(istate)%qp(jj)%el  = EQ(istate)%qp(jj)%el + Vel
EQ(istate)%qp(jj)%vdw = EQ(istate)%qp(jj)%vdw + V_a - V_b
#endif
else
select case (ST_gc(ene_header%gcnum(jj))%caltype)
case (FULL)
!do nothing!
case (ELECTRO)
#ifdef _OPENMP
qomp_vdw(istate,jj)   = qomp_vdw(istate,jj) + V_a - V_b
#else
EQ(istate)%qp(jj)%vdw = EQ(istate)%qp(jj)%vdw + V_a -V_b
#endif
case (VDW)
#ifdef _OPENMP
qomp_elec(istate,jj)  = qomp_elec(istate,jj) + Vel
#else
EQ(istate)%qp(jj)%el  = EQ(istate)%qp(jj)%el + Vel
#endif
end select
end if
end do
end if
end do ! istate

end do
#ifdef _OPENMP
do istate=1, nstates
do jj=1, ene_header%arrays
!$omp atomic update
EQ(istate)%qp(jj)%el  = EQ(istate)%qp(jj)%el + qomp_elec(istate,jj)
!$omp atomic update
EQ(istate)%qp(jj)%vdw = EQ(istate)%qp(jj)%vdw + qomp_vdw(istate,jj)
end do
end do
#endif
end subroutine nonbond_qp_box

!-----------------------------------------------------------------------

subroutine nonbond_qw

! local variables
integer						:: jw,iq,i,j,jj,i3,j3,iw
integer						:: istate
real(kind=prec)					:: dx1,dx2,dx3,r2,r,r6,r6_hc,dv,Vel,V_a,V_b,r12
#ifdef _OPENMP
integer :: quotient, remainder
real(kind=prec),allocatable :: Vel_omp(:),Vvdw_omp(:)
#endif
! global variables used:
!  iqseq, iac, crg, x, nstates, qvdw_flag, iaclib, qiac, qavdw, qbvdw, qcrg, el14_scale, EQ, d, nat_solute


#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (nbww_pair/solv_atom)/threads_num
remainder = MOD(nbww_pair/solv_atom, threads_num)
mp_start = thread_id * quotient + 1 + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
mp_start = (mp_start*solv_atom) - solv_atom + 1

if (remainder .gt. thread_id) then
mp_end = mp_end + 1
endif
mp_end = (mp_end*solv_atom) - solv_atom + 1
allocate(Vel_omp(istate),Vvdw_omp(istate))
qomp_elec = zero
qomp_vdw  = zero

do iw = mp_start,mp_end,solv_atom
#else
do iw = 1, nbqw_pair, solv_atom
#endif
! init state-invariant variables:
iq   = nbqw(iw,1)%i
i = iqseq(iq)
i3   = i*3-3


do jw = 0, solv_atom-1
! for every assigned q-w pair:

! init state-invariant variables:
j    = nbqw(iw+jw,1)%j
j3   = j*3-3


dx1  = x(j3+1) - x(i3+1)
dx2  = x(j3+2) - x(i3+2)
dx3  = x(j3+3) - x(i3+3)



r2    = dx1*dx1 + dx2*dx2 + dx3*dx3
r6    = r2*r2*r2
r6_hc = r6
r2    = one/r2
r     = sqrt ( r2 )

do istate = 1, nstates
r6    = r6_hc + nbqw(iw+jw,istate)%score !softcore
r6    = one/r6
r12   = r6*r6
! calculate qi, Vel, V_a, V_b and dv
Vel  = nbqw(iw+jw,istate)%elec *r
V_a  = nbqw(iw+jw,istate)%vdWA *r12
V_b  = nbqw(iw+jw,istate)%vdWB *r6
dv   = r2*( -Vel -(12.0_prec*V_a -6.0_prec*V_b)*r6*r6_hc )*EQ(istate)%lambda   !softcore r6*r6_hc is (r^6/(r^6+alpha))

! update forces
d(i3+1) = d(i3+1) - dv*dx1
d(i3+2) = d(i3+2) - dv*dx2
d(i3+3) = d(i3+3) - dv*dx3
d(j3+1) = d(j3+1) + dv*dx1
d(j3+2) = d(j3+2) + dv*dx2
d(j3+3) = d(j3+3) + dv*dx3

do jj = 1, ene_header%arrays
#ifdef _OPENMP
qomp_elec(istate,jj) = qomp_elec(istate,jj) + Vel
qomp_vdw(istate,jj)  = qomp_vdw(istate,jj) + V_a - V_b
#else
EQ(istate)%qw(jj)%el  = EQ(istate)%qw(jj)%el + Vel
EQ(istate)%qw(jj)%vdw = EQ(istate)%qw(jj)%vdw + V_a - V_b 
#endif
end do

end do ! nstates
end do !jw

end do ! nbqw
#ifdef _OPENMP
do istate=1, nstates
do jj=1, ene_header%arrays
!$omp atomic update
EQ(istate)%qp(jj)%el  = EQ(istate)%qp(jj)%el + qomp_elec(istate,jj)
!$omp atomic update
EQ(istate)%qp(jj)%vdw = EQ(istate)%qp(jj)%vdw + qomp_vdw(istate,jj)
end do
end do
#endif
end subroutine nonbond_qw

!-----------------------------------------------------------------------
!******PWadded 2001-10-23
subroutine nonbond_qw_box
	! local variables
integer                                         :: jw,iq,i,j,jj,i3,j3,iw
integer						:: istate
real(kind=prec)                                 :: dx1,dx2,dx3,r2,r,r6,r6_hc,dv,Vel,V_a,V_b,r12
real(kind=prec)					:: boxshiftx,boxshifty,boxshiftz
#ifdef _OPENMP
integer :: quotient, remainder
#endif
	! global variables used:
	!  iqseq, iac, crg, x, nstates, qvdw_flag, iaclib, qiac, qavdw, qbvdw, qcrg, el14_scale, EQ, d, nat_solute
  
#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (nbww_pair/solv_atom)/threads_num
remainder = MOD(nbww_pair/solv_atom, threads_num)
mp_start = thread_id * quotient + 1 + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
mp_start = (mp_start*solv_atom) - solv_atom + 1

if (remainder .gt. thread_id) then
mp_end = mp_end + 1
endif
mp_end = (mp_end*solv_atom) - solv_atom + 1

qomp_elec = zero
qomp_vdw  = zero

do iw = mp_start,mp_end,solv_atom
#else
do iw = 1, nbqw_pair, solv_atom
#endif
! for every assigned q-s pair:
! init state-invariant variables:
iq   = nbqw(iw,1)%i
i = iqseq(iq)
i3   = i*3-3
do jw = 0, solv_atom-1

! init state-invariant variables:
j    = nbqw(iw+jw,1)%j
j3   = j*3-3

!compute the periodical shift
dx1 = x(3*j-2) - x(3*qswitch-2)
dx2 = x(3*j-1) - x(3*qswitch-1)
dx3 = x(3*j  ) - x(3*qswitch  )
boxshiftx = boxlength(1)*nint(dx1*inv_boxl(1))
boxshifty = boxlength(2)*nint(dx2*inv_boxl(2))
boxshiftz = boxlength(3)*nint(dx3*inv_boxl(3))

dx1  = x(j3+1) - x(i3+1)
dx2  = x(j3+2) - x(i3+2)
dx3  = x(j3+3) - x(i3+3)
dx1 = dx1 - boxshiftx
dx2 = dx2 - boxshifty
dx3 = dx3 - boxshiftz

r2    = dx1*dx1 + dx2*dx2 + dx3*dx3
r6    = r2*r2*r2
r6_hc = r6
r2    = one/r2
r     = sqrt ( r2 )


do istate = 1, nstates
r6    = r6_hc + nbqw(iw+jw,istate)%score !softcore
r6    = one/r6
r12   = r6*r6
! calculate qi, Vel, V_a, V_b and dv
Vel  = nbqw(iw+jw,istate)%elec *r
V_a  = nbqw(iw+jw,istate)%vdWA *r12
V_b  = nbqw(iw+jw,istate)%vdWB *r6
dv   = r2*( -Vel -(12.0_prec*V_a -6.0_prec*V_b)*r6*r6_hc )*EQ(istate)%lambda   !softcore r6*r6_hc is (r^6/(r^6+alpha))

! update forces
d(i3+1) = d(i3+1) - dv*dx1
d(i3+2) = d(i3+2) - dv*dx2
d(i3+3) = d(i3+3) - dv*dx3
d(j3+1) = d(j3+1) + dv*dx1
d(j3+2) = d(j3+2) + dv*dx2
d(j3+3) = d(j3+3) + dv*dx3

do jj = 1, ene_header%arrays
#ifdef _OPENMP
qomp_elec(istate,jj) = qomp_elec(istate,jj) + Vel
qomp_vdw(istate,jj)  = qomp_vdw(istate,jj) + V_a - V_b
#else
EQ(istate)%qw(jj)%el  = EQ(istate)%qw(jj)%el + Vel
EQ(istate)%qw(jj)%vdw = EQ(istate)%qw(jj)%vdw + V_a - V_b
#endif
end do

end do ! nstates
end do ! jw

end do ! nbqw
#ifdef _OPENMP
do istate=1, nstates
do jj=1, ene_header%arrays
!$omp atomic update
EQ(istate)%qp(jj)%el  = EQ(istate)%qp(jj)%el + qomp_elec(istate,jj)
!$omp atomic update
EQ(istate)%qp(jj)%vdw = EQ(istate)%qp(jj)%vdw + qomp_vdw(istate,jj)
end do
end do
#endif
end subroutine nonbond_qw_box

!----------------------------------------------------------------------------------------------------

subroutine nonbond_ww
! local variables
integer						:: iw,ip,i,j,i3,j3,ia,ichg
integer						:: ja
integer						:: ipstart
real(kind=prec)						:: dx1,dx2,dx3,r2,r,r6,r12
real(kind=prec)						:: Vel,V_a,V_b,dv
#ifdef _OPENMP
integer :: quotient, remainder
real(kind=prec) :: Vel_omp,Vvdw_omp
#endif
! global variables used:
!  nat_solute, iac, crg, ljcod, iaclib, x, d, E

! totally rewritten now
! we are now using the nbww_pair index and nbww true list
! this function before was garbage
ichg = 1
#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (nbww_pair/solv_atom)/threads_num
remainder = MOD(nbww_pair/solv_atom, threads_num)
mp_start = thread_id * quotient + 1 + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
mp_start = (mp_start*solv_atom) - solv_atom + 1
if (remainder .gt. thread_id) then
mp_end = mp_end + 1
endif
mp_end = (mp_end*solv_atom) - solv_atom + 1

ichg = ceiling(REAL(MOD(mp_start,solv_atom**2),kind=prec)/solv_atom)
Vel_omp=zero
Vvdw_omp=zero
do iw = mp_start,mp_end,solv_atom
#else
do iw = 1, nbww_pair, solv_atom
#endif
if (ichg .gt. solv_atom) ichg=1
i    = nbww(iw)%i
i3   = i*3-3
do ip = 0, solv_atom - 1
j    = nbww(iw+ip)%j
j3   = j*3-3

dx1  = x(j3+1) - x(i3+1)
dx2  = x(j3+2) - x(i3+2)
dx3  = x(j3+3) - x(i3+3)



r2   = dx1*dx1 + dx2*dx2 + dx3*dx3
r2   = one/r2
r    = sqrt ( r2 ) 
r6   = r2*r2*r2
r12  = r6*r6
Vel  = ww_precomp(ichg,ip+1)%elec *r
V_a  = ww_precomp(ichg,ip+1)%vdWA *r12 
V_b  = ww_precomp(ichg,ip+1)%vdWB *r6
dv   = r2*( -Vel -12.0_prec*V_a +6.0_prec*V_b )
d(i3+1) = d(i3+1) - dv*dx1
d(i3+2) = d(i3+2) - dv*dx2
d(i3+3) = d(i3+3) - dv*dx3
d(j3+1) = d(j3+1) + dv*dx1
d(j3+2) = d(j3+2) + dv*dx2
d(j3+3) = d(j3+3) + dv*dx3
#ifdef _OPENMP
Vel_omp = Vel_omp + Vel
Vvdw_omp = Vvdw_omp + V_a - V_b
#else
E%ww%el  = E%ww%el + Vel       
E%ww%vdw = E%ww%vdw + V_a - V_b
#endif

end do ! ip
ichg = ichg + 1
end do ! iw
#ifdef _OPENMP
!$omp atomic update
E%ww%el  = E%ww%el + Vel_omp
!$omp atomic update
E%ww%vdw = E%ww%vdw + Vvdw_omp
#endif

end subroutine nonbond_ww

!----------------------------------------------------------------------------------------------------
subroutine nonbond_ww_box
  ! local variables
  integer						:: iw,ip,i,j,i3,j3,ia,ichg
  integer						:: ja
  integer						:: ipstart
  real(kind=prec)						:: dx1,dx2,dx3,r2,r,r6,r12
  real(kind=prec)						:: Vel,V_a,V_b,dv
  real(kind=prec)						:: do1, do2, do3
#ifdef _OPENMP
integer :: quotient, remainder
real(kind=prec) :: Vel_omp,Vvdw_omp
#endif
  ! global variables used:
  !  nat_solute, iac, crg, ljcod, iaclib, x, d, E

! totally rewritten now
! we are now using the nbww_pair index and nbww true list
! this function before was garbage
ichg = 1
#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (nbww_pair/solv_atom)/threads_num
remainder = MOD(nbww_pair/solv_atom, threads_num)
mp_start = thread_id * quotient + 1 + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
mp_start = (mp_start*solv_atom) - solv_atom + 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif
mp_end = (mp_end*solv_atom) - solv_atom + 1

ichg = ceiling(REAL(MOD(mp_start,solv_atom**2),kind=prec)/solv_atom)
Vel_omp=zero
Vvdw_omp=zero
do iw = mp_start,mp_end,solv_atom
#else
do iw = 1, nbww_pair, solv_atom
#endif
if (ichg .gt. solv_atom) ichg=1
i    = nbww(iw)%i
i3   = i*3-3
j    = nbww(iw)%j
j3   = j*3-3

!distance between this oxygen atom and the oxygen atom of the above watermolecule, iw
!this is always the first interaction :P
!only need to calculate this for ichg = 1
if (ichg .eq.1) then
do1 = x(j3+1) - x( i3+1 )
do2 = x(j3+2) - x( i3+2 )
do3 = x(j3+3) - x( i3+3 )
!the peridic shift
do1 = boxlength(1)*nint( do1*inv_boxl(1) )
do2 = boxlength(2)*nint( do2*inv_boxl(2) )
do3 = boxlength(3)*nint( do3*inv_boxl(3) )
end if

do ip = 0, solv_atom - 1
j    = nbww(iw+ip)%j
j3   = j*3-3

dx1  = x(j3+1) - x(i3+1)
dx2  = x(j3+2) - x(i3+2)
dx3  = x(j3+3) - x(i3+3)
dx1 = dx1 - do1
dx2 = dx2 - do2
dx3 = dx3 - do3
r2   = dx1*dx1 + dx2*dx2 + dx3*dx3
r2   = one/r2
r    = sqrt ( r2 ) 
r6   = r2*r2*r2
r12  = r6*r6

Vel  = ww_precomp(ichg,ip+1)%elec *r
V_a  = ww_precomp(ichg,ip+1)%vdWA *r12
V_b  = ww_precomp(ichg,ip+1)%vdWB *r6
dv   = r2*( -Vel -12.0_prec*V_a +6.0_prec*V_b )

d(i3+1) = d(i3+1) - dv*dx1
d(i3+2) = d(i3+2) - dv*dx2
d(i3+3) = d(i3+3) - dv*dx3
d(j3+1) = d(j3+1) + dv*dx1
d(j3+2) = d(j3+2) + dv*dx2
d(j3+3) = d(j3+3) + dv*dx3

#ifdef _OPENMP
Vel_omp = Vel_omp + Vel
Vvdw_omp = Vvdw_omp + V_a - V_b
#else
E%ww%el  = E%ww%el + Vel       
E%ww%vdw = E%ww%vdw + V_a - V_b
#endif
end do ! ip
ichg = ichg + 1
end do ! iw
#ifdef _OPENMP
!$omp atomic update
E%ww%el  = E%ww%el + Vel_omp
!$omp atomic update
E%ww%vdw = E%ww%vdw + Vvdw_omp
#endif

end subroutine nonbond_ww_box

!-----------------------------------------------------------------------

subroutine nonbond_qw_spc
!calculate non-bonded interactions between Q-atoms and SPC water molecules
!(optimisations rely on LJ params = 0 for water H) using geometric comb. rule

! local variables
integer						:: iw,jw,iq,i,j,jj,ip,ja,i3,j3
integer						:: istate
real(kind=prec)						::	dx1,dx2,dx3
real(kind=prec)						::	r, r2, r6, r12
real(kind=prec)						::	Vel, dv
real(kind=prec)						::  V_a, V_b, r6_hc
#ifdef _OPENMP
integer :: quotient, remainder
#endif
! has been extensivly rewritten to use n atom spc solvent
! and the new nbqw list
! old function made Paul very sad

! global variables used:
!  iqseq, iac, crg, x, nstates, qvdw_flag, iaclib, qiac, qavdw, qbvdw, qcrg, el14_scale, EQ, d, nat_solute

#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (nbqw_pair/solv_atom)/threads_num
remainder = MOD(nbqw_pair/solv_atom, threads_num)
mp_start = thread_id * quotient + 1 + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
mp_start = (mp_start*solv_atom) - solv_atom + 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif
mp_end = (mp_end*solv_atom) 

qomp_elec = zero
qomp_vdw  = zero

do iw = mp_start,mp_end,solv_atom
#else
do iw = 1, nbqw_pair, solv_atom
#endif
! init state-invariant variables:
iq   = nbqw(iw,1)%i
i    = iqseq(iq)
i3   = i*3-3
j    = nbqw(iw,1)%j
j3   = j*3-3

! only q - O distance first, only one with vdW
dx1  = x(j3+1) - x(i3+1)
dx2  = x(j3+2) - x(i3+2)
dx3  = x(j3+3) - x(i3+3)

r2    = dx1*dx1 + dx2*dx2 + dx3*dx3
r6    = r2*r2*r2
r6_hc = r6
r2    = one/r2
r     = sqrt ( r2 )

dv = zero
do istate = 1, nstates
r6    = r6_hc + nbqw(iw,istate)%score !softcore
r6    = one/r6
r12   = r6*r6
! calculate qi, Vel, V_a, V_b and dv
V_a = nbqw(iw,istate)%vdWA *r6*r6
V_b = nbqw(iw,istate)%vdWB *r6
Vel = nbqw(iw,istate)%elec *r
dv  = dv  + r2*( -Vel -( (12.0_prec*V_a - 6.0_prec*V_b)*r6*r6_hc ))*EQ(istate)%lambda
!softcore r6*r6_hc is (r^6/(r^6+alpha))
! update q-water energies
#ifdef _OPENMP
do jj=1,ene_header%arrays
qomp_elec(istate,jj) = qomp_elec(istate,jj) + Vel
qomp_vdw(istate,jj)  = qomp_vdw(istate,jj) + V_a - V_b
end do ! ene_header%arrays
#else
do jj=1,ene_header%arrays
EQ(istate)%qw(jj)%el  = EQ(istate)%qw(jj)%el + Vel
EQ(istate)%qw(jj)%vdw = EQ(istate)%qw(jj)%vdw + V_a - V_b
end do ! ene_header%arrays
#endif
end do !istate
! update forces on q atom
d(i3+1) = d(i3+1) - dv*dx1
d(i3+2) = d(i3+2) - dv*dx2
d(i3+3) = d(i3+3) - dv*dx3
                ! update forces on water
d(j3+1) = d(j3+1) + dv*dx1
d(j3+2) = d(j3+2) + dv*dx2
d(j3+3) = d(j3+3) + dv*dx3

!now calculate only charge-charge interaction for hydrogens
do ip = 1, solv_atom-1

ja   = j + ip
j3   = ja*3-3

dx1  = x(j3+1) - x(i3+1)
dx2  = x(j3+2) - x(i3+2)
dx3  = x(j3+3) - x(i3+3)

r2   = one/(dx1*dx1 + dx2*dx2 + dx3*dx3)
r    = sqrt(r2)

dv = 0
do istate = 1, nstates
Vel = nbqw(iw+ip,istate)%elec *r
#ifdef _OPENMP
do jj=1,ene_header%arrays
qomp_elec(istate,jj) = qomp_elec(istate,jj) + Vel
end do ! ene_header
#else
do jj=1,ene_header%arrays
EQ(istate)%qw(jj)%el  = EQ(istate)%qw(jj)%el + Vel
end do ! ene_header
#endif
dv = dv - r2*Vel*EQ(istate)%lambda
end do ! istate
! update forces on q atom
d(i3+1) = d(i3+1) - dv*dx1
d(i3+2) = d(i3+2) - dv*dx2
d(i3+3) = d(i3+3) - dv*dx3
! update forces on water
d(j3+1) = d(j3+1) + dv*dx1
d(j3+2) = d(j3+2) + dv*dx2
d(j3+3) = d(j3+3) + dv*dx3
end do ! ip

end do ! nbqw
#ifdef _OPENMP
do istate = 1, nstates
do jj=1,ene_header%arrays
!$omp atomic update
EQ(istate)%qw(jj)%vdw = EQ(istate)%qw(jj)%vdw + qomp_vdw(istate,jj)
!$omp atomic update
EQ(istate)%qw(jj)%el  = EQ(istate)%qw(jj)%el + qomp_elec(istate,jj)
end do ! ene_header
end do ! nstates
#endif
end subroutine nonbond_qw_spc

!-----------------------------------------------------------------------
!******PWadded 2001-10-23
subroutine nonbond_qw_spc_box
	!calculate non-bonded interactions between Q-atoms and SPC water molecules
	!(optimisations rely on LJ params = 0 for water H) using geometric comb. rule

	! local variables
integer                                         :: iw,iq,i,j,iLJ,jj,ip,ja,i3,j3
integer						:: istate
real(kind=prec)                                         ::      dx1,dx2,dx3
real(kind=prec)                                         ::      r, r2, r6, r12
real(kind=prec)                                         ::      Vel, dv, boxshiftx,boxshifty,boxshiftz
real(kind=prec)                                         ::  V_a, V_b, r6_hc
#ifdef _OPENMP
integer :: quotient, remainder
#endif

! has been extensivly rewritten to use n atom spc solvent
! and the new nbqw list

#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (nbqw_pair/solv_atom)/threads_num
remainder = MOD(nbqw_pair/solv_atom, threads_num)
mp_start = thread_id * quotient + 1 + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
mp_start = (mp_start*solv_atom) - solv_atom + 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif
mp_end = (mp_end*solv_atom) 
qomp_elec = zero
qomp_vdw  = zero

do iw = mp_start,mp_end,solv_atom
#else

do iw = 1, nbqw_pair, solv_atom
#endif
! init state-invariant variables:
iq   = nbqw(iw,1)%i
i    = iqseq(iq)
i3   = i*3-3
j    = nbqw(iw,1)%j
j3   = j*3-3

! first get the periodic shift done
dx1  = x(j3+1) - x(3*qswitch-2)
dx2  = x(j3+2) - x(3*qswitch-1)
dx3  = x(j3+3) - x(3*qswitch  )

!compute the periodical shift
boxshiftx = boxlength(1)*nint(dx1*inv_boxl(1))
boxshifty = boxlength(2)*nint(dx2*inv_boxl(2))
boxshiftz = boxlength(3)*nint(dx3*inv_boxl(3))

! only q - O distance first, only one with vdW
dx1  = x(j3+1) - x(i3+1)
dx2  = x(j3+2) - x(i3+2)
dx3  = x(j3+3) - x(i3+3)

dx1 = dx1 - boxshiftx
dx2 = dx2 - boxshifty
dx3 = dx3 - boxshiftz

r2    = dx1*dx1 + dx2*dx2 + dx3*dx3
r6    = r2*r2*r2
r6_hc = r6
r2    = one/r2
r     = sqrt ( r2 )


dv = zero
do istate = 1, nstates
r6    = r6_hc + nbqw(iw,istate)%score !softcore
r6    = one/r6
r12   = r6*r6

! calculate qi, Vel, V_a, V_b and dv
V_a = nbqw(iw,istate)%vdWA *r6*r6
V_b = nbqw(iw,istate)%vdWB *r6
Vel = nbqw(iw,istate)%elec *r
dv  = dv  + r2*( -Vel -( (12.0_prec*V_a - 6.0_prec*V_b)*r6*r6_hc))*EQ(istate)%lambda
!softcore r6*r6_hc is (r^6/(r^6+alpha))
! update q-water energies
#ifdef _OPENMP
do jj=1,ene_header%arrays
qomp_elec(istate,jj) = qomp_elec(istate,jj) + Vel
qomp_vdw(istate,jj)  = qomp_vdw(istate,jj) + V_a - V_b
end do ! ene_header%arrays
#else
do jj=1,ene_header%arrays
EQ(istate)%qw(jj)%el  = EQ(istate)%qw(jj)%el + Vel
EQ(istate)%qw(jj)%vdw = EQ(istate)%qw(jj)%vdw + V_a - V_b
end do ! ene_header%arrays
#endif
end do !istate			
! update forces on q atom
d(i3+1) = d(i3+1) - dv*dx1
d(i3+2) = d(i3+2) - dv*dx2
d(i3+3) = d(i3+3) - dv*dx3
! update forces on water
d(j3+1) = d(j3+1) + dv*dx1
d(j3+2) = d(j3+2) + dv*dx2
d(j3+3) = d(j3+3) + dv*dx3

!now calculate only charge-charge interaction for hydrogens
do ip = 1, solv_atom-1

ja   = j + ip
j3   = ja*3-3

dx1  = x(j3+1) - x(i3+1) - boxshiftx
dx2  = x(j3+2) - x(i3+2) - boxshifty
dx3  = x(j3+3) - x(i3+3) - boxshiftz

r2   = one/(dx1*dx1 + dx2*dx2 + dx3*dx3)
r    = sqrt(r2)

dv = 0
do istate = 1, nstates
Vel = nbqw(iw+ip,istate)%elec *r
#ifdef _OPENMP
do jj=1,ene_header%arrays
qomp_elec(istate,jj) = qomp_elec(istate,jj) + Vel
end do ! ene_header%arrays
#else
do jj=1,ene_header%arrays
EQ(istate)%qw(jj)%el  = EQ(istate)%qw(jj)%el + Vel
end do ! ene_header
#endif
dv = dv - r2*Vel*EQ(istate)%lambda
end do ! istate
! update forces on q atom
d(i3+1) = d(i3+1) - dv*dx1
d(i3+2) = d(i3+2) - dv*dx2
d(i3+3) = d(i3+3) - dv*dx3
! update forces on water
d(j3+1) = d(j3+1) + dv*dx1
d(j3+2) = d(j3+2) + dv*dx2
d(j3+3) = d(j3+3) + dv*dx3
end do ! ip

end do ! nbqw
#ifdef _OPENMP
do istate = 1, nstates
do jj=1,ene_header%arrays
!$omp atomic update
EQ(istate)%qw(jj)%vdw = EQ(istate)%qw(jj)%vdw + qomp_vdw(istate,jj)
!$omp atomic update
EQ(istate)%qw(jj)%el  = EQ(istate)%qw(jj)%el + qomp_elec(istate,jj)
end do ! ene_header
end do ! nstates
#endif
end subroutine nonbond_qw_spc_box


!-----------------------------------------------------------------------

subroutine nonbond_ww_spc
! function has been totally rewritten to be used for n - atom spc solvent
! only vdW interactions between first two atoms are calculated -> those are the heavy atoms
! rest of the atoms only interact via coloumb interactions
! old function gave me nightmares

! local variables
integer                                         :: iw,ip,i,j,i3,j3,ia,ichg,ja
integer						:: ipstart
real(kind=prec)                                 :: dx1,dx2,dx3,r2,r,r6
real(kind=prec)					:: Vel,V_a,V_b,dv
#ifdef _OPENMP
integer :: quotient, remainder
real(kind=prec) :: Vel_omp, Vvdw_omp
#endif

! the loop is now organised in the following way (as before for the arithmetic case)
! we always jump solv_atom positions, to calculate atomAn-atomB(1-solv_atom)
! and iterate a separate counter to get the current position of the atomA(1-solv_atom) position
! vdW is only calculated if atomA1 - atomB1

ichg = 1
#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (nbww_pair/solv_atom)/threads_num
remainder = MOD(nbww_pair/solv_atom, threads_num)
mp_start = thread_id * quotient + 1 + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
mp_start = (mp_start*solv_atom) - solv_atom + 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif
mp_end = (mp_end*solv_atom) 

ichg = ceiling(REAL(MOD(mp_start,solv_atom**2),kind=prec)/solv_atom)
Vel_omp = zero
Vvdw_omp = zero
do iw = mp_start,mp_end,solv_atom
#else
do iw = 1, nbww_pair, solv_atom
#endif
if (ichg .gt. solv_atom) ichg=1
i    = nbww(iw)%i
i3   = i*3-3
j    = nbww(iw)%j
do ip = 0, solv_atom - 1
ja   = j + ip
j3   = ja*3-3
dx1  = x(j3+1) - x(i3+1)
dx2  = x(j3+2) - x(i3+2)
dx3  = x(j3+3) - x(i3+3)
r2   = dx1*dx1 + dx2*dx2 + dx3*dx3 
r2   = one/r2
r    = sqrt ( r2 )
Vel  = ww_precomp(ichg,ip+1)%elec *r
#ifdef _OPENMP
Vel_omp = Vel_omp + Vel
#else
E%ww%el = E%ww%el + Vel
#endif
dv   = r2*( -Vel)
if ((ichg .eq. 1).and.(ip .eq. 0)) then
r6   = r2*r2*r2
V_a = ww_precomp(ichg,ip+1)%vdWA * r6 * r6
V_b = ww_precomp(ichg,ip+1)%vdWB * r6
#ifdef _OPENMP
Vvdw_omp = Vvdw_omp + V_a - V_b
#else
E%ww%vdw = E%ww%vdw + V_a - V_b
#endif
dv = dv + r2*(-12.0_prec*V_a +6.0_prec*V_b )
end if
d(i3+1) = d(i3+1) - dv*dx1
d(i3+2) = d(i3+2) - dv*dx2
d(i3+3) = d(i3+3) - dv*dx3
d(j3+1) = d(j3+1) + dv*dx1
d(j3+2) = d(j3+2) + dv*dx2
d(j3+3) = d(j3+3) + dv*dx3
end do ! ip
ichg = ichg + 1
end do ! iw
#ifdef _OPENMP
!$omp atomic update
E%ww%el = E%ww%el + Vel_omp
!$omp atomic update
E%ww%vdw = E%ww%vdw + Vvdw_omp
#endif
end subroutine nonbond_ww_spc

!-----------------------------------------------------------------------
!******PWadded 2001-10-23
subroutine nonbond_ww_spc_box

! function has been totally rewritten to be used for n - atom spc solvent
! only vdW interactions between first two atoms are calculated -> those are the heavy atoms
! rest of the atoms only interact via coloumb interactions
! old function gave me nightmares

  ! local variables
integer                                         :: iw,ip,i,j,i3,j3,ia,ichg,ja
integer						:: ipstart
real(kind=prec)                                         :: dx1,dx2,dx3,r2,r,r6
real(kind=prec)                                         :: Vel,V_a,V_b,dv,boxshiftx,boxshifty,boxshiftz
#ifdef _OPENMP
integer :: quotient, remainder
real(kind=prec) :: Vel_omp, Vvdw_omp
#endif

! the loop is now organised in the following way (as before for the arithmetic case)
! we always jump solv_atom positions, to calculate atomAn-atomB(1-solv_atom)
! and iterate a separate counter to get the current position of the atomA(1-solv_atom) position
! vdW is only calculated if atomA1 - atomB1

ichg = 1

#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (nbww_pair/solv_atom)/threads_num
remainder = MOD(nbww_pair/solv_atom, threads_num)
mp_start = thread_id * quotient + 1 + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
mp_start = (mp_start*solv_atom) - solv_atom + 1
if (remainder .gt. thread_id) then
mp_end = mp_end + 1
endif
mp_end = (mp_end*solv_atom) 

ichg = ceiling(REAL(MOD(mp_start,solv_atom**2),kind=prec)/solv_atom)

Vel_omp = zero
Vvdw_omp = zero

do iw = mp_start,mp_end,solv_atom
#else
do iw = 1, nbww_pair, solv_atom
#endif
if (ichg .gt. solv_atom) ichg=1
i    = nbww(iw)%i
i3   = i*3-3
j    = nbww(iw)%j
do ip = 0, solv_atom - 1
ja   = j + ip
j3   = ja*3-3
dx1  = x(j3+1) - x(i3+1)
dx2  = x(j3+2) - x(i3+2)
dx3  = x(j3+3) - x(i3+3)
if ((ichg .eq. 1).and.(ip .eq. 0)) then
boxshiftx = boxlength(1)*nint( dx1*inv_boxl(1) )
boxshifty = boxlength(2)*nint( dx2*inv_boxl(2) )
boxshiftz = boxlength(3)*nint( dx3*inv_boxl(3) )
end if
dx1 = dx1 -boxshiftx
dx2 = dx2 -boxshifty
dx3 = dx3 -boxshiftz
r2   = dx1*dx1 + dx2*dx2 + dx3*dx3
r2   = one/r2
r    = sqrt ( r2 )
Vel  = ww_precomp(ichg,ip+1)%elec *r
#ifdef _OPENMP
Vel_omp = Vel_omp + Vel
#else
E%ww%el = E%ww%el + Vel
#endif
dv   = r2*( -Vel)
if ((ichg .eq. 1).and.(ip .eq. 0)) then
r6   = r2*r2*r2
V_a = ww_precomp(ichg,ip+1)%vdWA * r6 * r6
V_b = ww_precomp(ichg,ip+1)%vdWB * r6
#ifdef _OPENMP
Vvdw_omp = Vvdw_omp + V_a - V_b
#else
E%ww%vdw = E%ww%vdw + V_a - V_b
#endif
dv = dv + r2*(-12.0_prec*V_a +6.0_prec*V_b )
end if
d(i3+1) = d(i3+1) - dv*dx1
d(i3+2) = d(i3+2) - dv*dx2
d(i3+3) = d(i3+3) - dv*dx3
d(j3+1) = d(j3+1) + dv*dx1
d(j3+2) = d(j3+2) + dv*dx2
d(j3+3) = d(j3+3) + dv*dx3

end do ! ip
ichg = ichg + 1
end do ! iw
#ifdef _OPENMP
!$omp atomic update
E%ww%el = E%ww%el + Vel_omp
!$omp atomic update
E%ww%vdw = E%ww%vdw + Vvdw_omp
#endif
end subroutine nonbond_ww_spc_box


!-----------------------------------------------------------------------

subroutine nonbond_solvent_internal
! new function that handles non bonded interactions within a solvent molecule
! only called for solvents with several heavy atoms and if the solvent has torsions
! all interactions are precomputed, only need distances
! each node only works on its own solvent atoms
! local variables
integer                                         :: iw,is,ip,ia,ib,ia3,ib3
real(kind=prec)					:: dx1,dx2,dx3,r2,r6,r
real(kind=prec)					:: Vel,V_a,V_b,dv
#ifdef _OPENMP
integer :: quotient, remainder
real(kind=prec) :: Vel_omp, Vvdw_omp
#endif

#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (calculation_assignment%ww%end - calculation_assignment%ww%start + 1)/threads_num
remainder = MOD(calculation_assignment%ww%end - calculation_assignment%ww%start + 1, threads_num)
mp_start = thread_id * quotient + calculation_assignment%ww%start + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
mp_start = (mp_start*solv_atom) - solv_atom + 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif

Vel_omp = zero
Vvdw_omp = zero

do iw = mp_start,mp_end
#else
do iw = calculation_assignment%ww%start, calculation_assignment%ww%end
#endif
is  = nat_solute + solv_atom*iw-(solv_atom)

! loop over all interactions within this solvent molecule
do ip = 1 , num_solv_int 
ia   = is + nonbnd_solv_int(ip)%i
ib   = is + nonbnd_solv_int(ip)%j

ia3  = ia*3-3
ib3  = ib*3-3

dx1  = x(ib3+1) - x(ia3+1)
dx2  = x(ib3+2) - x(ia3+2)
dx3  = x(ib3+3) - x(ia3+3)

r2   = dx1*dx1 + dx2*dx2 + dx3*dx3
r6   = r2*r2*r2
r2   = one/r2
r    = sqrt ( r2 )
r6   = one/r6

Vel  = nonbnd_solv_int(ip)%elec * r
V_a  = nonbnd_solv_int(ip)%vdWA * r6 *r6
V_b  = nonbnd_solv_int(ip)%vdWB * r6

dv   = r2*( -Vel) + r2*(-12.0_prec*V_a +6.0_prec*V_b )
#ifdef _OPENMP
Vel_omp  = Vel_omp + Vel
Vvdw_omp = Vvdw_omp + V_a - V_b
#else
E%ww%el  = E%ww%el + Vel
E%ww%vdw = E%ww%vdw + V_a - V_b
#endif
d(ia3+1) = d(ia3+1) - dv*dx1
d(ia3+2) = d(ia3+2) - dv*dx2
d(ia3+3) = d(ia3+3) - dv*dx3
d(ib3+1) = d(ib3+1) + dv*dx1
d(ib3+2) = d(ib3+2) + dv*dx2
d(ib3+3) = d(ib3+3) + dv*dx3
end do
end do
#ifdef _OPENMP
!$omp atomic update
E%ww%el = E%ww%el + Vel_omp
!$omp atomic update
E%ww%vdw = E%ww%vdw + Vvdw_omp
#endif

end subroutine nonbond_solvent_internal

!-----------------------------------------------------------------------

subroutine nonbond_solvent_internal_box
! new function that handles non bonded interactions within a solvent molecule
! only called for solvents with several heavy atoms and if the solvent has torsions
! all interactions are precomputed, only need distances
! each node only works on its own solvent atoms
! local variables
integer                                         :: iw,is,ip,ia,ib,ia3,ib3
real(kind=prec)                                 :: dx1,dx2,dx3,r2,r6,r
real(kind=prec)					:: boxshiftx,boxshifty,boxshiftz
real(kind=prec)                                 :: Vel,V_a,V_b,dv
#ifdef _OPENMP
integer :: quotient, remainder
real(kind=prec) :: Vel_omp, Vvdw_omp
#endif

#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (calculation_assignment%ww%end - calculation_assignment%ww%start + 1)/threads_num
remainder = MOD(calculation_assignment%ww%end - calculation_assignment%ww%start + 1, threads_num)
mp_start = thread_id * quotient + calculation_assignment%ww%start + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
mp_start = (mp_start*solv_atom) - solv_atom + 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif

Vel_omp = zero
Vvdw_omp = zero

do iw = mp_start,mp_end
#else
do iw = calculation_assignment%ww%start, calculation_assignment%ww%end
#endif
is  = nat_solute + solv_atom*iw-(solv_atom)

! loop over all interactions within this solvent molecule
do ip = 1 , num_solv_int
ia   = is + nonbnd_solv_int(ip)%i
ib   = is + nonbnd_solv_int(ip)%j

ia3  = ia*3-3
ib3  = ib*3-3

dx1  = x(ib3+1) - x(ia3+1)
dx2  = x(ib3+2) - x(ia3+2)
dx3  = x(ib3+3) - x(ia3+3)

boxshiftx = boxlength(1)*nint( dx1*inv_boxl(1) )
boxshifty = boxlength(2)*nint( dx2*inv_boxl(2) )
boxshiftz = boxlength(3)*nint( dx3*inv_boxl(3) )

dx1 = dx1 -boxshiftx
dx2 = dx2 -boxshifty
dx3 = dx3 -boxshiftz

r2   = dx1*dx1 + dx2*dx2 + dx3*dx3
r6   = r2*r2*r2
r2   = one/r2
r    = sqrt ( r2 )
r6   = one/r6

Vel  = nonbnd_solv_int(ip)%elec * r
V_a  = nonbnd_solv_int(ip)%vdWA * r6 *r6
V_b  = nonbnd_solv_int(ip)%vdWB * r6

dv   = r2*( -Vel) + r2*(-12.0_prec*V_a +6.0_prec*V_b )
#ifdef _OPENMP
Vel_omp  = Vel_omp + Vel
Vvdw_omp = Vvdw_omp + V_a - V_b
#else
E%ww%el  = E%ww%el + Vel
E%ww%vdw = E%ww%vdw + V_a - V_b
#endif
d(ia3+1) = d(ia3+1) - dv*dx1
d(ia3+2) = d(ia3+2) - dv*dx2
d(ia3+3) = d(ia3+3) - dv*dx3
d(ib3+1) = d(ib3+1) + dv*dx1
d(ib3+2) = d(ib3+2) + dv*dx2
d(ib3+3) = d(ib3+3) + dv*dx3
end do
end do
#ifdef _OPENMP
!$omp atomic update
E%ww%el = E%ww%el + Vel_omp
!$omp atomic update
E%ww%vdw = E%ww%vdw + Vvdw_omp
#endif


end subroutine nonbond_solvent_internal_box

!-----------------------------------------------------------------------

subroutine offdiag
! local variables
integer						:: io,i,j,k,l,k3,l3
real(kind=prec)						:: r

! global variables used:
!  offd, noffd, iqseq, x, Hij, offd2

do io = 1, noffd
! for every offd:

i  = offd(io)%i
j  = offd(io)%j
k  = iqseq(offd2(io)%k)
l  = iqseq(offd2(io)%l)
k3 = k*3-3
l3 = l*3-3

r = sqrt ( (x(l3+1)-x(k3+1))**2 + &
  (x(l3+2)-x(k3+2))**2 + &
  (x(l3+3)-x(k3+3))**2 )

Hij(i,j) = offd2(io)%A * exp(-offd2(io)%mu*r)
offd(io)%Hij = Hij(i,j)	! store for save
offd(io)%rkl = r
end do
end subroutine offdiag

!-----------------------------------------------------------------------
subroutine p_restrain
! *** Local variables
integer :: ir,i,j,k,i3,j3,k3,istate, n_ctr
real(kind=prec) :: fk,r2,erst,Edum,x2,y2,z2,wgt,b,db,dv, totmass, theta,rij,r2ij,rjk,r2jk, scp,f1
real(kind=prec) :: dr(3), dr2(3), ctr(3), di(3), dk(3)
real(kind=prec)						::	fexp

! global variables used:
!  E, nstates, EQ, nrstr_seq, rstseq, heavy, x, xtop, d, nrstr_pos, rstpos, nrstr_dist, 
!  rstdis, nrstr_wall, rstang, nrstr_ang, rstwal, xwcent, excl, freeze

! sequence restraints (independent of Q-state)
do ir = 1, nrstr_seq
fk = rstseq(ir)%fk

if(rstseq(ir)%to_centre == 1) then     ! Put == 1, then equal to 2
  ! restrain to geometrical centre

  ! reset dr & atom counter
  dr(:) = zero
  n_ctr = 0

  ! calculate deviation from center
  do i = rstseq(ir)%i, rstseq(ir)%j
        if ( heavy(i) .or. rstseq(ir)%ih .eq. 1 .and.(.not.(excl(i).and.freeze))) then
          n_ctr = n_ctr + 1
          dr(:) = dr(:) + x(i*3-2:i*3) - xtop(i*3-2:i*3)
        end if
  end do

  if(n_ctr > 0) then 
    ! only if atoms were found:

        ! form average
        dr(:) = dr(:) / n_ctr 
        r2      = dr(1)**2 + dr(2)**2 + dr(3)**2
        erst    = 0.5_prec*fk*r2
        E%restraint%protein  = E%restraint%protein + erst

        ! apply same force to all atoms
        do i = rstseq(ir)%i, rstseq(ir)%j
        if ( heavy(i) .or. rstseq(ir)%ih .eq. 1 .and.(.not.(excl(i).and.freeze))) then
                d(i*3-2:i*3) = d(i*3-2:i*3) + fk*dr(:)*iaclib(iac(i))%mass/12.010_prec
          end if
        end do
  end if
 
else if(rstseq(ir)%to_centre == 2) then     ! Put == 1, then equal to 2
  ! restrain to mass centre
  ! reset dr & variable to put masses
  dr(:) = zero
  totmass = zero
  
! calculate deviation from mass center
  do i = rstseq(ir)%i, rstseq(ir)%j
        if ( heavy(i) .or. rstseq(ir)%ih .eq. 1 .and.(.not.(excl(i).and.freeze))) then
          totmass = totmass + iaclib(iac(i))%mass                              ! Add masses
          dr(:) = dr(:) + (x(i*3-2:i*3) - xtop(i*3-2:i*3))*iaclib(iac(i))%mass ! Massweight distances
		end if
  end do

 if(totmass > 0) then 
    ! only if atoms were found: (i.e has a total mass)

        ! form average
        dr(:) = dr(:)/totmass                                  ! divide by total mass
        r2      = dr(1)**2 + dr(2)**2 + dr(3)**2
        erst    = 0.5_prec*fk*r2
        E%restraint%protein  = E%restraint%protein + erst

        ! apply same force to all atoms
        do i = rstseq(ir)%i, rstseq(ir)%j
        if ( heavy(i) .or. rstseq(ir)%ih .eq. 1 .and.(.not.(excl(i).and.freeze))) then
                d(i*3-2:i*3) = d(i*3-2:i*3) + fk*dr(:)
          end if
        end do
  end if

else 
  ! restrain each atom to its topology co-ordinate
  do i = rstseq(ir)%i, rstseq(ir)%j
        if ( heavy(i) .or. rstseq(ir)%ih .eq. 1 .and.(.not.(excl(i).and.freeze))) then
          i3 = i*3-3

          dr(1)   = x(i3+1) - xtop(i3+1)
          dr(2)   = x(i3+2) - xtop(i3+2)
          dr(3)   = x(i3+3) - xtop(i3+3)
       !use the periodically minimal distance:
       if( use_PBC ) then
               dr(1) = dr(1) - boxlength(1)*nint( dr(1)*inv_boxl(1) )
               dr(2) = dr(2) - boxlength(2)*nint( dr(2)*inv_boxl(2) )
               dr(3) = dr(3) - boxlength(3)*nint( dr(3)*inv_boxl(3) )
       end if
          r2      = dr(1)**2 + dr(2)**2 + dr(3)**2

          erst    = 0.5_prec*fk*r2
          E%restraint%protein  = E%restraint%protein + erst

          d(i3+1) = d(i3+1) + fk*dr(1)
          d(i3+2) = d(i3+2) + fk*dr(2)
          d(i3+3) = d(i3+3) + fk*dr(3)
        end if
  end do
end if
end do

! extra positional restraints (Q-state dependent)
do ir = 1, nrstr_pos
istate = rstpos(ir)%ipsi
i      = rstpos(ir)%i
i3     = i*3-3

dr(1)  = x(i3+1) - rstpos(ir)%x(1)
dr(2)  = x(i3+2) - rstpos(ir)%x(2)
dr(3)  = x(i3+3) - rstpos(ir)%x(3)

if ( istate .ne. 0 ) then
wgt = EQ(istate)%lambda
else
wgt = one
end if

x2 = dr(1)**2
y2 = dr(2)**2
z2 = dr(3)**2

Edum = 0.5_prec*rstpos(ir)%fk(1)*x2 + &
   0.5_prec*rstpos(ir)%fk(2)*y2 + &
   0.5_prec*rstpos(ir)%fk(3)*z2

d(i3+1) = d(i3+1) + rstpos(ir)%fk(1)*dr(1)*wgt
d(i3+2) = d(i3+2) + rstpos(ir)%fk(2)*dr(2)*wgt
d(i3+3) = d(i3+3) + rstpos(ir)%fk(3)*dr(3)*wgt

if ( istate .eq. 0 ) then
do k = 1, nstates
EQ(k)%restraint = EQ(k)%restraint + Edum
end do
if ( nstates .eq. 0 ) E%restraint%protein = E%restraint%protein + Edum
else
EQ(istate)%restraint = EQ(istate)%restraint + Edum
end if
end do

! atom-atom distance restraints (Q-state dependent)
do ir = 1, nrstr_dist
istate = rstdis(ir)%ipsi
i      = rstdis(ir)%i
j      = rstdis(ir)%j
i3     = i*3-3
j3     = j*3-3

dr(1)  = x(j3+1) - x(i3+1)
dr(2)  = x(j3+2) - x(i3+2)
dr(3)  = x(j3+3) - x(i3+3)

! if PBC then adjust lengths according to periodicity - MA
if( use_PBC ) then
	dr(1) = dr(1) - boxlength(1)*nint( dr(1)*inv_boxl(1) )
	dr(2) = dr(2) - boxlength(2)*nint( dr(2)*inv_boxl(2) )
	dr(3) = dr(3) - boxlength(3)*nint( dr(3)*inv_boxl(3) )
end if

if ( istate .ne. 0 ) then
wgt = EQ(istate)%lambda
else
wgt = one
end if

b      = sqrt ( dr(1)**2 + dr(2)**2 + dr(3)**2 )
if(b < rstdis(ir)%d1) then !shorter than d1
        db     = b - rstdis(ir)%d1
elseif(b > rstdis(ir)%d2) then !longer than d2
        db     = b - rstdis(ir)%d2
else
        db = 0
        cycle !skip zero force calculation
endif

Edum   = 0.5_prec*rstdis(ir)%fk*db**2
dv     = wgt*rstdis(ir)%fk*db/b

d(j3+1) = d(j3+1) + dr(1)*dv
d(j3+2) = d(j3+2) + dr(2)*dv
d(j3+3) = d(j3+3) + dr(3)*dv
d(i3+1) = d(i3+1) - dr(1)*dv
d(i3+2) = d(i3+2) - dr(2)*dv
d(i3+3) = d(i3+3) - dr(3)*dv

if ( istate .eq. 0 ) then
do k = 1, nstates
EQ(k)%restraint = EQ(k)%restraint + Edum
end do
if ( nstates .eq. 0 ) E%restraint%protein = E%restraint%protein + Edum
else
EQ(istate)%restraint = EQ(istate)%restraint + Edum
end if
end do

! atom-atom-atom angle restraints (Q-state dependent)
do ir = 1, nrstr_angl

istate = rstang(ir)%ipsi
i      = rstang(ir)%i
j      = rstang(ir)%j
k      = rstang(ir)%k
i3     = i*3-3
j3     = j*3-3
k3     = k*3-3

! distance from atom i to atom j
dr(1)  = x(i3+1) - x(j3+1)
dr(2)  = x(i3+2) - x(j3+2)
dr(3)  = x(i3+3) - x(j3+3)

! distance from atom k to atom j
dr2(1) = x(k3+1) - x(j3+1)
dr2(2) = x(k3+2) - x(j3+2)
dr2(3) = x(k3+3) - x(j3+3)

! if PBC then adjust lengths according to periodicity - MA
if( use_PBC ) then
        dr(1) = dr(1) - boxlength(1)*nint( dr(1)*inv_boxl(1) )
        dr(2) = dr(2) - boxlength(2)*nint( dr(2)*inv_boxl(2) )
        dr(3) = dr(3) - boxlength(3)*nint( dr(3)*inv_boxl(3) )

        dr2(1) = dr2(1) - boxlength(1)*nint( dr2(1)*inv_boxl(1) )
        dr2(2) = dr2(2) - boxlength(2)*nint( dr2(2)*inv_boxl(2) )
        dr2(3) = dr2(3) - boxlength(3)*nint( dr2(3)*inv_boxl(3) )

end if

if ( istate .ne. 0 ) then
wgt = EQ(istate)%lambda
else
wgt = one
end if

! square distances from the triangle formed by the atoms i, j and k
r2ij = dr(1)**2 + dr(2)**2 + dr(3)**2
r2jk = dr2(1)**2 + dr2(2)**2 + dr2(3)**2

rij  = sqrt ( r2ij )
rjk  = sqrt ( r2jk )

! calculate the scalar product (scp) and the angle theta from it
scp = ( dr(1)*dr2(1) + dr(2)*dr2(2) + dr(3)*dr2(3) )
scp = scp / ( rij*rjk )

! criteria inserted on the real angle force calculations to ensure no weird scp.
if ( scp .gt.  one ) scp =  one
if ( scp .lt. -one ) scp = -one

theta = acos(scp)
db    = theta - (rstang(ir)%ang)*deg2rad

! dv is the force to be added in module
Edum   = 0.5_prec*rstang(ir)%fk*db**2
dv     = wgt*rstang(ir)%fk*db

! calculate sin(theta) to use in forces
f1 = sin ( theta )
if ( abs(f1) .lt. 1.e-12_prec ) then
        ! avoid division by zero
        f1 = -1.e12_prec
else
        f1 =  -one / f1
end if

        ! calculate di and dk
di(1) = f1 * ( dr2(1) / ( rij * rjk ) - scp * dr(1) / r2ij )
di(2) = f1 * ( dr2(2) / ( rij * rjk ) - scp * dr(2) / r2ij )
di(3) = f1 * ( dr2(3) / ( rij * rjk ) - scp * dr(3) / r2ij )
dk(1) = f1 * ( dr(1) / ( rij * rjk ) - scp * dr2(1) / r2jk )
dk(2) = f1 * ( dr(2) / ( rij * rjk ) - scp * dr2(2) / r2jk )
dk(3) = f1 * ( dr(3) / ( rij * rjk ) - scp * dr2(3) / r2jk )

        ! update d
d(i3+1) = d(i3+1) + dv*di(1)
d(i3+2) = d(i3+2) + dv*di(2)
d(i3+3) = d(i3+3) + dv*di(3)
d(k3+1) = d(k3+1) + dv*dk(1)
d(k3+2) = d(k3+2) + dv*dk(2)
d(k3+3) = d(k3+3) + dv*dk(3)
d(j3+1) = d(j3+1) - dv*( di(1) + dk(1) )
d(j3+2) = d(j3+2) - dv*( di(2) + dk(2) )
d(j3+3) = d(j3+3) - dv*( di(3) + dk(3) )


if ( istate .eq. 0 ) then
do k = 1, nstates
EQ(k)%restraint = EQ(k)%restraint + Edum
end do
if ( nstates .eq. 0 ) E%restraint%protein = E%restraint%protein + Edum
else
EQ(istate)%restraint = EQ(istate)%restraint + Edum
end if
end do
if( .not. use_PBC ) then
! extra half-harmonic wall restraints
do ir = 1, nrstr_wall
        fk = rstwal(ir)%fk
        do i = rstwal(ir)%i, rstwal(ir)%j
                if ( heavy(i) .or. rstwal(ir)%ih .eq. 1 ) then
                        i3 = i*3-3

                        dr(1) = x(i3+1) - xwcent(1)
                        dr(2) = x(i3+2) - xwcent(2)
                        dr(3) = x(i3+3) - xwcent(3)

                        b     = sqrt ( dr(1)**2 + dr(2)**2 + dr(3)**2 )
                        db = b - rstwal(ir)%d

                        if(db > zero) then
                                erst =  0.5_prec * fk * db**2 - rstwal(ir)%Dmorse
                                dv = fk*db/b
                        else
                                fexp = exp(rstwal(ir)%aMorse*db)
                                erst = rstwal(ir)%dMorse*(fexp*fexp-2.0_prec*fexp)
                                dv=-2.0_prec*rstwal(ir)%dMorse*rstwal(ir)%aMorse*(fexp-fexp*fexp)/b
                        end if
                        E%restraint%protein = E%restraint%protein + erst

                        d(i3+1) = d(i3+1) + dv*dr(1)
                        d(i3+2) = d(i3+2) + dv*dr(2)
                        d(i3+3) = d(i3+3) + dv*dr(3)
                end if
        end do
end do

end if

end subroutine p_restrain

!-----------------------------------------------------------------------
#ifdef USE_GRID
subroutine populate_grids
! this function searches the number of charge groups and number of water molecules
! and places them into the new md grids
! this is done every NB steps in make_pair_lists as the first thing
! in the end master broadcasts the new grids to the slave nodes
! needs some trickery so that we don't overflow any arrays

integer					:: ichg,igrid
integer					:: i,i3
real(kind=prec)				:: boxmin(1:3),boxmax(1:3),xc(1:3)
logical                                 :: set = .false.
#ifdef USE_MPI
integer, parameter                      :: vars = 40    !increment this var when adding data to broadcast in batch 1
integer                                 :: blockcnt(vars), ftype(vars)
integer(kind=MPI_ADDRESS_KIND)                          :: fdisp(vars)
#endif

!each node now does its share

if (use_PBC) then
boxmin(1:3) = boxcentre(1:3)-boxlength/2
boxmax(1:3) = boxcentre(1:3)+boxlength/2
end if

! reset the group information for all grids to zero
! and null all arrays
grid_pp_ngrp(:)=0
grid_pw_ngrp(:)=0
grid_ww_ngrp(:)=0
ww_igrid(:)=0
pw_igrid(:)=0
pp_igrid(:)=0
grid_pp_grp(:,:)=0
grid_pw_grp(:,:)=0
grid_ww_grp(:,:)=0

do ichg=calculation_assignment%pp%start,calculation_assignment%pp%end
! we decide the grid by the position of the switch atom of each charge group
i   = cgp(ichg)%iswitch
i3  = i*3-3
xc(1) = x(i3+1)
xc(2) = x(i3+2)
xc(3) = x(i3+3)
! this is to account for the possibility of atoms to move out of the PBC box if they are not put back in
! in this case we just shift the atom implicitly by the distance
! this is just a bookkeeping issue
do while (use_PBC.and.((xc(1).gt.(boxmax(1))).or.(xc(1).lt.(boxmin(1))).or.(xc(2).gt.(boxmax(2))).or.(xc(2).lt.(boxmin(2))).or.(xc(3).gt.(boxmax(3))).or.(xc(3).lt.(boxmin(3)))))
if (xc(1).gt.boxmax(1)) then
xc(1) = boxmin(1) + (xc(1)-boxmax(1))
elseif (xc(1).lt.boxmin(1)) then
xc(1) = boxmax(1) - (xc(1)-boxmin(1))
end if
if (xc(2).gt.boxmax(2)) then
xc(2) = boxmin(2) + (xc(2)-boxmax(2))
elseif (xc(2).lt.boxmin(2)) then
xc(2) = boxmax(2) - (xc(2)-boxmin(2))
end if
if (xc(3).gt.boxmax(3)) then
xc(3) = boxmin(3) + (xc(3)-boxmax(3))
elseif (xc(3).lt.boxmin(3)) then
xc(3) = boxmax(3) - (xc(3)-boxmin(3))
end if
end do

set   = .false.
igrid = 1
! now loop first over the pp_grid and then over the pw_grid for ncgp_solute
pploop: do while (( igrid .le. pp_gridnum).and.(.not.set))
		if ((xc(1).ge.grid_pp(igrid)%x).and.(xc(1).lt.grid_pp(igrid)%xend)) then
			if ((xc(2).ge.grid_pp(igrid)%y).and.(xc(2).lt.grid_pp(igrid)%yend)) then
				if ((xc(3).ge.grid_pp(igrid)%z).and.(xc(3).lt.grid_pp(igrid)%zend)) then  
! yeah, we found our grid, now place the atom in here and set all the bookkeeping stuff
! then cycle the grid and continue with the next one
! needs temp storage because of the OMP part!!!!!
				pp_igrid(ichg) = igrid
				grid_pp_ngrp(igrid) = grid_pp_ngrp(igrid) + 1
				grid_pp_grp(igrid,grid_pp_ngrp(igrid)) = ichg
                                set = .true.
				cycle pploop
! this should (!!!!) ensure we place all groups, but needs to be tested
				end if
			end if
		end if
                igrid = igrid + 1
	end do pploop
! now the pw_grid, same procedure
set   = .false.
igrid = 1
pwloop1: do while (( igrid .le. pw_gridnum).and.(.not.set))
		if ((xc(1).ge.grid_pw(igrid)%x).and.(xc(1).lt.grid_pw(igrid)%xend)) then
			if ((xc(2).ge.grid_pw(igrid)%y).and.(xc(2).lt.grid_pw(igrid)%yend)) then
				if ((xc(3).ge.grid_pw(igrid)%z).and.(xc(3).lt.grid_pw(igrid)%zend)) then
				pw_igrid(ichg) = igrid
                                set = .true.
				cycle pwloop1
				end if
			end if
		end if
                igrid = igrid + 1
	end do pwloop1
end do ! ncgp_solute

! now do the same for nwat
do ichg=calculation_assignment%ww%start,calculation_assignment%ww%end
! we decide the grid by the position of the first atom of each solvent molecule
i   = nat_solute+(ichg*solv_atom)-(solv_atom-1)
i3  = i*3-3
xc(1) = x(i3+1)
xc(2) = x(i3+2)
xc(3) = x(i3+3)
! this is to account for the possibility of atoms to move out of the PBC box if they are not put back in
! in this case we just shift the atom implicitly by the distance
! this is just a bookkeeping issue
do while (use_PBC.and.((xc(1).gt.(boxmax(1))).or.(xc(1).lt.(boxmin(1))).or.(xc(2).gt.(boxmax(2))).or.(xc(2).lt.(boxmin(2))).or.(xc(3).gt.(boxmax(3))).or.(xc(3).lt.(boxmin(3)))))
if (xc(1).gt.boxmax(1)) then
xc(1) = boxmin(1) + (xc(1)-boxmax(1))
elseif (xc(1).lt.boxmin(1)) then
xc(1) = boxmax(1) - (xc(1)-boxmin(1))
end if
if (xc(2).gt.boxmax(2)) then
xc(2) = boxmin(2) + (xc(2)-boxmax(2))
elseif (xc(2).lt.boxmin(2)) then
xc(2) = boxmax(2) - (xc(2)-boxmin(2))
end if
if (xc(3).gt.boxmax(3)) then
xc(3) = boxmin(3) + (xc(3)-boxmax(3))
elseif (xc(3).lt.boxmin(3)) then
xc(3) = boxmax(3) - (xc(3)-boxmin(3))
end if
end do

set   = .false.
igrid = 1
pwloop2: do while (( igrid .le. pw_gridnum).and.(.not.set)) 
                if ((xc(1).ge.grid_pw(igrid)%x).and.(xc(1).lt.grid_pw(igrid)%xend)) then
                        if ((xc(2).ge.grid_pw(igrid)%y).and.(xc(2).lt.grid_pw(igrid)%yend)) then
                                if ((xc(3).ge.grid_pw(igrid)%z).and.(xc(3).lt.grid_pw(igrid)%zend)) then
				grid_pw_ngrp(igrid) = grid_pw_ngrp(igrid) + 1
				grid_pw_grp(igrid,grid_pw_ngrp(igrid)) = ichg
                                set = .true.
                                cycle pwloop2
                                end if
                        end if
                end if
                igrid = igrid + 1
        end do pwloop2
set   = .false.
igrid = 1
wwloop: do while (( igrid .le. ww_gridnum).and.(.not.set))
		if ((xc(1).ge.grid_ww(igrid)%x).and.(xc(1).lt.grid_ww(igrid)%xend)) then
			if ((xc(2).ge.grid_ww(igrid)%y).and.(xc(2).lt.grid_ww(igrid)%yend)) then
				if ((xc(3).ge.grid_ww(igrid)%z).and.(xc(3).lt.grid_ww(igrid)%zend)) then
				grid_ww_ngrp(igrid) = grid_ww_ngrp(igrid) + 1
				grid_ww_grp(igrid,grid_ww_ngrp(igrid)) = ichg
				ww_igrid(ichg) = igrid
                                set = .true.
				cycle wwloop
				end if
			end if
		end if
                igrid = igrid + 1
	end do wwloop
end do ! nwat


! now collect and distribute all the info to all nodes
#ifdef USE_MPI

if (ncgp_solute .ne. 0 ) then

call MPI_Allreduce(MPI_IN_PLACE,grid_pp_ngrp,pp_gridnum,MPI_INTEGER,MPI_SUM, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('populate grids MPI_Allreduce grid_pp_ngrp')
call MPI_Allreduce(MPI_IN_PLACE,pp_igrid,ncgp_solute,MPI_INTEGER, MPI_SUM,MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('populate grids MPI_Allreduce pp_igrid')

do igrid = 1, pp_gridnum
call MPI_Allreduce(MPI_IN_PLACE,grid_pp_grp(igrid,:),gridstor_pp,MPI_INTEGER,&
                        mpi_grid_add, MPI_COMM_WORLD, ierr)
end do
if (ierr .ne. 0) call die('populate grids MPI_Allreduce grid_pp_grp')

end if
if (nwat .gt. 0 ) then
call MPI_Allreduce(MPI_IN_PLACE,grid_pw_ngrp,pw_gridnum,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('populate grids MPI_Allreduce grid_pw_ngrp')
call MPI_Allreduce(MPI_IN_PLACE,pw_igrid,ncgp_solute,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('populate grids MPI_Allreduce pw_igrid')
call MPI_Allreduce(MPI_IN_PLACE,grid_ww_ngrp,ww_gridnum,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('populate grids MPI_Allreduce grid_ww_ngrp')
call MPI_Allreduce(MPI_IN_PLACE,ww_igrid,nwat,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('populate grids MPI_Allreduce ww_igrid')

do igrid = 1, pw_gridnum
call MPI_Allreduce(MPI_IN_PLACE,grid_pw_grp(igrid,:),gridstor_pw,MPI_INTEGER,&
                        mpi_grid_add, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('populate grids MPI_Allreduce grid_pw_grp')
end do
do igrid = 1, ww_gridnum
call MPI_Allreduce(MPI_IN_PLACE,grid_ww_grp(igrid,:),gridstor_ww,MPI_INTEGER,&
                        mpi_grid_add, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('populate grids MPI_Allreduce grid_ww_grp')
end do
end if
#endif
 

end subroutine populate_grids
#endif

subroutine pot_energy
! local variables


integer					:: istate, i, nat3
integer					:: is, j,jj
#if defined (PROFILING)
real(kind=prec)					:: start_loop_time1
real(kind=prec)					:: start_loop_time2
real(kind=prec)					:: start_loop_time3
#endif

! --- reset all energies

E%potential = zero
!E%kinetic = zero	! no need to reset because it will be assigned its final value at once
E%LRF = zero
E%p%bond  = zero
E%p%angle = zero
E%p%torsion = zero
E%p%improper = zero
E%w%bond  = zero
E%w%angle = zero
E%w%torsion = zero
E%w%improper = zero
E%q%bond  = zero
E%q%angle = zero
E%q%torsion   = zero
E%q%improper   = zero
E%pp%el  = zero
E%pp%vdw = zero
E%pw%el  = zero
E%pw%vdw = zero
E%ww%el  = zero
E%ww%vdw = zero
E%qx%el    = zero
E%qx%vdw   = zero
!E%restraint%total = zero	! will be assigned its final value later
E%restraint%fix = zero
E%restraint%shell = zero
E%restraint%protein = zero
E%restraint%solvent_radial = zero
E%restraint%water_pol = zero
do istate = 1, nstates
!EQ(istate)%lambda set by initialize
!EQ(istate)%total assigned its final value later
EQ(istate)%q%bond = zero
EQ(istate)%q%angle = zero
EQ(istate)%q%torsion = zero
EQ(istate)%q%improper = zero
!EQ(istate)%qx%el = zero	! assigned its final value later
!EQ(istate)%qx%vdw = zero	! assigned its final value later
do jj=1,ene_header%arrays
EQ(istate)%qq(jj)%el = zero
EQ(istate)%qq(jj)%vdw = zero
EQ(istate)%qp(jj)%el = zero
EQ(istate)%qp(jj)%vdw = zero
EQ(istate)%qw(jj)%el = zero
EQ(istate)%qw(jj)%vdw = zero
end do
EQ(istate)%restraint = zero
end do

!reset derivatives ---
d(:) = zero

! --- calculate the potential energy and derivatives ---
! *** nonbonds distribueras

#if defined (USE_MPI)
if (nodeid .eq. 0) then
!First post recieves for gathering data from slaves
call gather_nonbond
end if
#endif

#if defined (PROFILING)
start_loop_time2 = rtime()
#endif

! classical nonbonds
call pot_energy_nonbonds
#if defined (PROFILING)
profile(10)%time = profile(10)%time + rtime() - start_loop_time2
#endif

if (nodeid .eq. 0) then


! classical bond interactions (master node only)
#if defined (PROFILING)
start_loop_time1 = rtime()
#endif
call pot_energy_bonds
#if defined (PROFILING)
profile(8)%time = profile(8)%time + rtime() - start_loop_time1
#endif

! various restraints
if( .not. use_PBC ) then
   call fix_shell     !Restrain all excluded atoms plus heavy solute atoms in the inner shell.
end if

call p_restrain       !Seq. restraints, dist. restaints, etc

if( .not. use_PBC ) then
        if(nwat > 0) then
          call restrain_solvent 
        if (wpol_restr) call watpol
        end if
end if

! q-q nonbonded interactions
call nonbond_qq
call nonbond_qqp

! q-atom bonded interactions: loop over q-atom states
do istate = 1, nstates
  ! bonds, angles, torsions and impropers
  call qbond (istate)
  call qangle (istate)
  if(ff_type == FF_CHARMM) call qurey_bradley(istate)
  call qtorsion (istate)
  call qimproper (istate)
end do
#if defined (PROFILING)
profile(9)%time = profile(9)%time + rtime() - start_loop_time1 - profile(8)%time
#endif
#if defined(USE_MPI)
else  !Slave nodes
call gather_nonbond
#endif
end if
if (nodeid .eq. 0) then 
#if (USE_MPI)
do i = 1, 3
    call MPI_WaitAll(numnodes-1,request_recv(1,i),mpi_status,ierr)
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
	do istate=1,nstates
	do jj=1,ene_header%arrays
  EQ(istate)%qp(jj)%el  = EQ(istate)%qp(jj)%el  + EQ_recv(istate,jj,i)%qp%el
  EQ(istate)%qp(jj)%vdw = EQ(istate)%qp(jj)%vdw + EQ_recv(istate,jj,i)%qp%vdw
  EQ(istate)%qw(jj)%el  = EQ(istate)%qw(jj)%el  + EQ_recv(istate,jj,i)%qw%el
  EQ(istate)%qw(jj)%vdw = EQ(istate)%qw(jj)%vdw + EQ_recv(istate,jj,i)%qw%vdw
	end do
	end do
end do
#endif
! q-atom energy summary
do istate = 1, nstates
! update EQ
do jj=1,ene_header%arrays
EQ(istate)%qx(jj)%el  = EQ(istate)%qq(jj)%el +EQ(istate)%qp(jj)%el +EQ(istate)%qw(jj)%el
EQ(istate)%qx(jj)%vdw = EQ(istate)%qq(jj)%vdw+EQ(istate)%qp(jj)%vdw+EQ(istate)%qw(jj)%vdw

EQ(istate)%total(jj) =  EQ(istate)%q%bond + EQ(istate)%q%angle   &
  + EQ(istate)%q%torsion  + EQ(istate)%q%improper + EQ(istate)%qx(jj)%el &
  + EQ(istate)%qx(jj)%vdw  + EQ(istate)%restraint
end do

! update E with an average of all states
E%q%bond  = E%q%bond  + EQ(istate)%q%bond *EQ(istate)%lambda
E%q%angle = E%q%angle + EQ(istate)%q%angle*EQ(istate)%lambda
E%q%torsion   = E%q%torsion   + EQ(istate)%q%torsion  *EQ(istate)%lambda
E%q%improper   = E%q%improper   + EQ(istate)%q%improper  *EQ(istate)%lambda
! only use full system to update total system energy -> what if we change this
! to get the effects on the total trajectory? (entropic stuff and so)???
E%qx%el    = E%qx%el    + EQ(istate)%qx(1)%el   *EQ(istate)%lambda
E%qx%vdw   = E%qx%vdw   + EQ(istate)%qx(1)%vdw  *EQ(istate)%lambda

! update E%restraint%protein with an average of all states
E%restraint%protein = E%restraint%protein + EQ(istate)%restraint*EQ(istate)%lambda
end do

! total energy summary
E%restraint%total = E%restraint%fix + E%restraint%shell + &
E%restraint%protein + E%restraint%solvent_radial + E%restraint%water_pol

E%potential = E%p%bond + E%w%bond + E%p%angle + E%w%angle + E%p%torsion + &
E%p%improper + E%pp%el + E%pp%vdw + E%pw%el + E%pw%vdw + E%ww%el + &
E%ww%vdw + E%q%bond + E%q%angle + E%q%torsion + &
E%q%improper + E%qx%el + E%qx%vdw + E%restraint%total + E%LRF
end if

end subroutine pot_energy

!-----------------------------------------------------------------------

subroutine pot_energy_bonds
! bond, angle, torsion and improper potential energy
select case(ff_type)
case(FF_GROMOS)
        E%p%bond = bond(1, nbonds_solute)
        E%w%bond = bond(nbonds_solute+1, nbonds)
        E%p%angle = angle(1, nangles_solute)
        E%w%angle = angle(nangles_solute+1, nangles)
        E%p%torsion = torsion(1, ntors_solute)
        E%w%torsion = torsion(ntors_solute+1, ntors)
        E%p%improper = improper(1, nimps_solute)
        E%w%improper = improper(nimps_solute+1, nimps)
case(FF_AMBER)
        E%p%bond = bond(1, nbonds_solute)
        E%w%bond = bond(nbonds_solute+1, nbonds)
        E%p%angle = angle(1, nangles_solute)
        E%w%angle = angle(nangles_solute+1, nangles)
        E%p%torsion = torsion(1, ntors_solute)
        E%w%torsion = torsion(ntors_solute+1, ntors)
        E%p%improper = improper2(1, nimps_solute)
        E%w%improper = improper2(nimps_solute+1, nimps)
case(FF_CHARMM)
        E%p%bond = bond(1, nbonds_solute)
        E%w%bond = bond(nbonds_solute+1, nbonds)
        E%p%angle = angle(1, nangles_solute)
        E%w%angle = angle(nangles_solute+1, nangles)
        E%p%angle = E%p%angle + urey_bradley(1, nangles_solute)
        E%w%angle = E%w%angle + urey_bradley(nangles_solute+1, nangles)
        E%p%torsion = torsion(1, ntors_solute)
        E%w%torsion = torsion(ntors_solute+1, ntors)
        E%p%improper = improper(1, nimps_solute)
        E%w%improper = improper(nimps_solute+1, nimps)
end select
end subroutine pot_energy_bonds

!-----------------------------------------------------------------------
subroutine pot_energy_nonbonds

!nonbonded interactions
!$omp parallel default(none) shared(use_PBC,ivdw_rule,natom,nat_solute,solvent_type,use_LRF,ntors,ntors_solute) &
!$omp reduction(+:d)
if( use_PBC ) then !periodic box
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
else
        if(natom > nat_solute) then
                if((ivdw_rule.eq.VDW_GEOMETRIC).and. &
                        (solvent_type == SOLVENT_SPC)) then
                        call nonbond_qw_spc
                        call nonbond_ww_spc
                else
                        call nonbond_qw
                        call nonbond_ww
                end if
! now we get started on the funky stuff with solvent torsions
! only run this if we have torsions in the solvent
! because it means we have at least 1-4 interactions
                if(ntors>ntors_solute) then
                        call nonbond_solvent_internal
                end if
                call nonbond_pw
        end if
        call nonbond_pp
        call nonbond_qp
end if

!LRF
!$omp single
if (use_LRF) then
        call lrf_taylor
end if
!$omp end single
!$omp end parallel
end subroutine pot_energy_nonbonds

!-----------------------------------------------------------------------

subroutine prep_coord


! local variables
integer(4)          :: i,nat3
logical             :: old_restart = .false.
real(8),allocatable :: x_old(:),v_old(:)
real(8)             :: old_boxlength(3), old_boxcentre(3)
integer             :: headercheck,myprec
!new variables for differen tprecision restart files
real(kind=singleprecision),allocatable :: x_single(:),v_single(:)
real(kind=singleprecision)             :: boxl_single(3),boxc_single(3)
real(kind=doubleprecision),allocatable :: x_double(:),v_double(:)
real(kind=doubleprecision)             :: boxl_double(3),boxc_double(3)
#ifndef PGI
real(kind=quadprecision),allocatable   :: x_quad(:),v_quad(:)
real(kind=quadprecision)               :: boxl_quad(3),boxc_quad(3)
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
         allocate(x_single(nat3))
         read (12,err=112,end=112) nat3, (x_single(i),i=1,nat_pro*3)
         xtop(1:nat_pro*3) = x_single(1:nat_pro*3)
         deallocate(x_single)
       else if (headercheck .eq. -1337) then
         allocate(x_double(nat3))
         read (12,err=112,end=112) nat3, (x_double(i),i=1,nat_pro*3)
         xtop(1:nat_pro*3) = x_double(1:nat_pro*3)
         deallocate(x_double)
       else if (headercheck .eq. -13337) then
#ifndef PGI
         allocate(x_quad(nat3))
         read (12,err=112,end=112) nat3, (x_quad(i),i=1,nat_pro*3)
         xtop(1:nat_pro*3) = x_quad(1:nat_pro*3)
         deallocate(x_quad)
#else
         call die('Quadruple precision not supported in PGI')
#endif
       end if
     else
       read (12,err=112,end=112) nat3, (xtop(i),i=1,nat_pro*3)
     end if
  else
     allocate(x_old(nat3))
     read (12,err=112,end=112) nat3, (x_old(i),i=1,nat_pro*3)
     xtop(1:nat_pro*3) = x_old(1:nat_pro*3)
     deallocate(x_old)
  end if
!read (12) nat3,(xtop(i),i=1,nat_pro*3)
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
              allocate(x_single(nat3),v_single(nat3))
              read (2,err=112,end=112) nat3, (x_single(i),i=1,nat3)
              read (2,err=112,end=112) nat3, (v_single(i),i=1,nat3)
              x(1:nat3) = x_single(1:nat3)
              v(1:nat3) = v_single(1:nat3)
              deallocate(x_single,v_single)
              if( use_PBC) then
                 read(2,err=112,end=112) boxl_single(:)
                 read(2,err=112,end=112) boxc_single(:)
                 boxlength(1:3) = boxl_single(1:3)
                 boxcentre(1:3) = boxc_single(1:3)
              end if
            else if (headercheck .eq. -1337) then
              allocate(x_double(nat3),v_double(nat3))
              read (2,err=112,end=112) nat3, (x_double(i),i=1,nat3)
              read (2,err=112,end=112) nat3, (v_double(i),i=1,nat3)
              x(1:nat3) = x_double(1:nat3)
              v(1:nat3) = v_double(1:nat3)
              deallocate(x_double,v_double)
              if( use_PBC) then
                 read(2,err=112,end=112) boxl_double(:)
                 read(2,err=112,end=112) boxc_double(:)
                 boxlength(1:3) = boxl_double(1:3)
                 boxcentre(1:3) = boxc_double(1:3)
              end if
            else if (headercheck .eq. -13337) then
#ifndef PGI
              allocate(x_quad(nat3),v_quad(nat3))
              read (2,err=112,end=112) nat3, (x_quad(i),i=1,nat3)
              read (2,err=112,end=112) nat3, (v_quad(i),i=1,nat3)
              x(1:nat3) = x_quad(1:nat3)
              v(1:nat3) = v_quad(1:nat3)
              deallocate(x_quad,v_quad)
              if( use_PBC) then
                 read(2,err=112,end=112) boxl_quad(:)
                 read(2,err=112,end=112) boxc_quad(:)
                 boxlength(1:3) = boxl_quad(1:3)
                 boxcentre(1:3) = boxc_quad(1:3)
              end if
#else
              call die('Quadruple precision not supported in PGI')
#endif
            end if
          else
              read (2,err=112,end=112) nat3, (x(i),i=1,nat3)
              read (2,err=112,end=112) nat3, (v(i),i=1,nat3)
              if( use_PBC) then
                read(2,err=112,end=112) boxlength(:)
                read(2,err=112,end=112) boxcentre(:)
              end if
          end if
        else
          allocate(x_old(nat3),v_old(nat3))
          read (2,err=112,end=112) nat3, (x_old(i),i=1,nat3)
          read (2,err=112,end=112) nat3, (v_old(i),i=1,nat3)
          write(*,*) 'Read coordinates and velocities from previous version of qdyn'
          x(1:nat3) = x_old(1:nat3)
          v(1:nat3) = v_old(1:nat3)
          deallocate(x_old,v_old)
          if( use_PBC) then
             read(2,err=112,end=112) old_boxlength(:)
             read(2,err=112,end=112) old_boxcentre(:)
             write(*,*) 'Read boxlength and center from previous version of qdyn'
             boxlength(1:3) = old_boxlength(1:3)
             boxcentre(1:3) = old_boxcentre(1:3)
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
        x(1:nat_pro*3) = xtop(1:nat_pro*3)
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




