! (C) 2014 Uppsala Molekylmekaniska HB, Uppsala, Sweden
! globals.f90, based on
! md.f90
! by Johan Åqvist, John Marelius, Anders Kaplan, Isabella Feierberg, Martin Nervall & Martin Almlöf
! global variables and structures used by routines
! by Paul Bauer

module GLOBALS

! used modules
!use PROFILING
use MPIGLOB
use ATOM_MASK
use QMATH
use SIZES

implicit none

!-----------------------------------------------------------------------
!	shared variables
!-----------------------------------------------------------------------
!	Constants
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
logical         :: temp_spam = .false. ! controlling if temperatures should be spammed to the user


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
TYPE(qr_vec)					:: new_boxl 

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
	ENUMERATOR	:: LEAPFROG,VELVERLET,MIN_STEEP,MIN_CG
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
#ifdef HAVEQUAD
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
#ifdef HAVEQUAD
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
type EXC_PARAMS
        character*10            ::seltype
        integer                 ::count=0,caltype=NOGC
        integer                 ::curcount=0
        character*80,allocatable::maskarray(:)
        character*30,allocatable::sendtomask(:)
        type(MASK_TYPE)         ::gcmask
        character*10            ::caltypen
end type EXC_PARAMS
type(EXC_PARAMS),allocatable :: ST_gc(:)

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
        TYPE(qr_vec)					::	fk
        TYPE(qr_vec)					::	x
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

! New ! for internal interactions of solvent atoms
! will be precomputed in prep sim
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


type(PRECOMPUTE_INT),allocatable                :: pp_precomp(:)
type(PRECOMPUTE_INT),allocatable                :: pw_precomp(:,:)
type(PRECOMPUTE_INT),allocatable                :: ww_precomp(:,:)
type(PRECOMPUTE_INT),allocatable                :: qp_precomp(:,:,:)
type(PRECOMPUTE_INT),allocatable                :: qw_precomp(:,:,:)
type(PRECOMPUTE_INT_QQ),allocatable             :: qq_precomp(:,:,:)
integer(AI),allocatable                         :: pp_map(:,:)
integer                                         :: pp_low,pw_low

!-----------------------------------------------------------------------
!	Coordinates, velocities, forces
!-----------------------------------------------------------------------
TYPE(qr_vec),allocatable		::	d(:)
TYPE(qr_vec),allocatable		::	x(:)
TYPE(qr_vec),allocatable		::	xx(:) !for shake
TYPE(qr_vec),allocatable		::	v(:)
real(kind=prec), allocatable			::	winv(:)
real(kind=prec)						::	grms !RMS force


!-----------------------------------------------------------------------
!	Energies , EQ is defined in qatom.f90
!-----------------------------------------------------------------------
type(ENERGIES)				::	E,old_E,previous_E
!added the variables for the MC_volume function as globals
real(kind=prec)						::	Tfree, Tfree_solvent, Tfree_solute, Temp
real(kind=prec)						::	Temp_solvent, Temp_solute, Texcl_solute, Texcl_solvent
real(kind=prec)                                         :: energy_cutoff

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
        TYPE(qr_vec)				:: shift !periodical shifts
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

! return type for nobonded energy calculations
type ENERET_TYPE
        real(kind=prec)                         :: Vel, V_a, V_b, dv
        TYPE(qr_vec)                            :: vec
end type ENERET_TYPE

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
        TYPE(qr_vec)			:: gstart           ! grid starting coordinates
        TYPE(qr_vec)			:: gend             ! grid end coordinates
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
integer,allocatable		:: nbqq_pair(:),nbqqp_pair(:)
type(NBQ_TYPE), allocatable ::nbqq(:,:)
type(NBQP_TYPE), allocatable ::nbqqp(:,:)

integer						::	nbqp_pair !current no of qatom-solute pairs
type(NBQP_TYPE), allocatable, target::nbqp(:,:)

type(NBQP_TYPE), allocatable, target::old_nbqp(:,:)
integer						::	nbqw_pair !current no of q-atom-water mol. pairs
type(NBQP_TYPE), allocatable	::	nbqw(:,:)

type(NBQ_TYPE),allocatable      :: monitor_group_int(:,:)

TYPE EXC_NBQ_TYPE
        integer                         :: inter
        TYPE(NBQ_TYPE),allocatable      :: list(:)
end TYPE EXC_NBQ_TYPE

TYPE EXC_NBQP_TYPE
        integer                         :: inter
        TYPE(NBQP_TYPE),allocatable     :: list(:)
end TYPE EXC_NBQP_TYPE

type(EXC_NBQ_TYPE),allocatable          :: exc_nbqq(:,:)
type(EXC_NBQP_TYPE),allocatable         :: exc_nbqqp(:,:)
type(EXC_NBQP_TYPE),allocatable         :: exc_nbqp(:,:)
! for internal solvent stuff
type(NB_TYPE),allocatable			:: nonbnd_solv_int(:)

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
        TYPE(qr_vec)					::	cgp_cent
        real(kind=prec)					::	phi0
        TYPE(qr_vec)					::	phi1
        TYPE(qr_vec)					::	phi2(3)
        TYPE(qr_vec)					::	phi3(9)
end type LRF_TYPE

type(LRF_TYPE), allocatable	::	lrf(:)
type(LRF_TYPE), allocatable	::	old_lrf(:)  !for constant pressure: MC_volume routine

type(node_assignment_type)	:: calculation_assignment


!constraint types & variables
!convergence criterion (fraction of distance)
real(kind=prec), parameter                      ::      CONST_TOL = 0.0001_prec
integer, parameter                              ::      CONST_MAX_ITER = 1000


type CONST_BOND_TYPE
        integer(AI)                             ::      i,j
        ! new for LINCS, number of constraints per constraint
        integer(AI)                             ::      ncc
        real(kind=prec)                         ::      dist2
        logical                                 ::      ready
end type CONST_BOND_TYPE

! additional variables needed for LINCS, mainly arrays of lengths and constants
type LINCS_MOL_TYPE
        integer(AI),allocatable                 ::      con(:,:)
        TYPE(qr_vec),allocatable                ::      B(:)
        real(kind=prec),allocatable             ::      A(:,:),S(:),coef(:,:)
        real(kind=prec),allocatable             ::      rhs(:,:),sol(:),length(:),length2(:)
end type LINCS_MOL_TYPE

type CONST_MOL_TYPE
        integer                                 ::      nconstraints
        type(LINCS_MOL_TYPE)                    ::      lin
        type(CONST_BOND_TYPE), allocatable      ::      bond(:)
end type CONST_MOL_TYPE


integer                                         ::      constraints, const_molecules
integer                                         ::      const_method = SHAKE_CONST
integer                                         ::      lincs_recursion = 4
type(CONST_MOL_TYPE), allocatable               ::      const_mol(:)


! variables for QCP calculations, concerning parameter read in and usage
integer                         :: QCP_N = QCP_OFF
integer                         :: QCP_steps_default = 10
integer                         :: QCP_size_default = 32
integer                         :: QCP_size_small = 16
integer                         :: QCP_size_verysmall = 4
integer                         :: QCP_size_large = 64
integer                         :: qcp_enum,qcp_size
integer                         :: qcp_level = -1
integer                         :: qcp_filen
character(80)                   :: qcp_pdb_name
! one for equilibration, second for sampling
integer                         :: qcp_steps(2),qcp_atnum,qcp_pos,qcp_pos2
logical                         :: use_qcp,qcp_veryverbose=.false.,qcp_verbose=.false.,use_qcp_mass=.false.
logical                         :: qcp_write = .true.
integer,allocatable             :: qcp_atom(:)
real(kind=prec),allocatable     :: qcp_mass(:)
TYPE(ENERGIES)                  :: qcp_E
TYPE(OQ_ENERGIES),allocatable   :: qcp_EQ(:)
TYPE(qr_vec),allocatable        :: x_save(:),d_save(:)
real(kind=prec),parameter       :: angstrom = 1E-10_prec
real(kind=prec),parameter       :: amu = 1.6605655E-27_prec
real(kind=prec),parameter       :: avogadro = 6.022E23_prec
real(kind=prec),parameter       :: convert = 4186.05_prec / avogadro
real(kind=prec),parameter       :: cboltz = 1.380662E-23_prec
real(kind=prec),parameter       :: planck = 6.626176E-34_prec
real(kind=prec),parameter       :: ddiff  = 1E-6_prec ! for checking equality

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

 !vectors for keeping track of node times
real(kind=prec),allocatable				:: all_node_times(:)
real(kind=prec),allocatable				:: node_times(:)

#endif

!-----------------------------------------------------------------------
!	temperature calculation variables
!-----------------------------------------------------------------------
integer						:: Ndegf,Ndegfree
logical						:: detail_temps			!controls whether or not solute and solvent temps are printed separately (true if solute and solvent degrees of freedom are both not zero)
integer                                         :: ntgroups,ntgroups_kind
TYPE TGROUP_TYPE
        real(kind=prec)                         :: sfact,exclshk,temp,tfree,texcl
        integer                                 :: starta,enda
        integer                                 :: Ndegf,Ndegfree,shake,nexcl
        character(20)                           :: tname
end TYPE TGROUP_TYPE
TYPE(TGROUP_TYPE),allocatable                   :: tscale(:)

!----END OF SHARED VARIABLES
end module GLOBALS
