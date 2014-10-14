! (C) 2014 Uppsala Molekylmekaniska HB, Uppsala, Sweden
! md.f90
! by Johan Åqvist, John Marelius, Anders Kaplan, Isabella Feierberg, Martin Nervall & Martin Almlöf
! molecular dynamics

module MD

! used modules
!use PROFILING
use SIZES
use TRJ
use MPIGLOB
use QATOM
use EXC
#if defined (_DF_VERSION_)
use DFPORT
use DFLIB
#endif
implicit none
#if defined (USE_MPI)
include "mpif.h"
#endif



!-----------------------------------------------------------------------
!	shared variables
!-----------------------------------------------------------------------
!	Constants
real(8)			::	pi, deg2rad	!set in sub startup
real(8),parameter		::	zero = 0.0
character*(*), parameter	::	MD_VERSION = '5.06'
character*(*), parameter	::	MD_DATE    = '2014-04-21'
real, parameter		::  rho_wat = 0.0335  ! molecules / A**3
real, parameter		::  Boltz = 0.001986

!Read status
integer                     :: stat

!print temperature if it changed more than 2% in one time step
real, parameter				::	TEMP_PRINT_THRESHOLD=0.02


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
real(8)						::	crg_ow, crg_hw, mu_w


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
real(8)						:: pressure, max_vol_displ
integer						:: pressure_seed
real(8), allocatable		:: mol_mass(:)
real(8), allocatable		:: mass(:)



!variabels used when user change boxsize in inputfile
logical						:: control_box
real(8)						:: new_boxl(3) 

!-----------------------------------------------------------------------
!	Dynamics control information     
!-----------------------------------------------------------------------


! --- MD parameters
! friction, gkT, randf and randv are variables to be used for the Langevin thermostat - A. Barrozo
! numchain and nhchain are variables to be used for the Nosé-Hoover thermostat
integer						::	nsteps, istep
integer						::	iseed, numchain
logical						::	restart
real(8)						::	dt, dt2, Temp0, Tmaxw, tau_T, friction, gkT, randf, randv, kbT, nhq
logical						::	shake_solvent, shake_solute
logical						::	shake_hydrogens
logical						::	separate_scaling
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
real(8), allocatable				::	xnh(:), vnh(:), qnh(:), Gnh(:) 

! --- Non-bonded strategy
logical						::	use_LRF
integer						::	NBcycle
real(8)						::	Rcpp, Rcww, Rcpw, Rcq, RcLRF
integer, parameter			                ::	max_atyp = 255
integer						::	ljcod(max_atyp,max_atyp)


! --- Output control parameters
integer						::	itrj_cycle, iene_cycle
integer						::	itemp_cycle, iout_cycle
logical						::	force_rms
!******PWadded 2001-10-23
integer						::	ivolume_cycle

! --- Protein boundary
logical						::	exclude_bonded
real(8)						::	fk_pshell
real(8)      :: fk_fix    = 200.0

! --- Water sphere
real(8)						::	Dwmz, awmz, rwat_in
real(8)						::	fkwpol
logical						::	wpol_restr, wpol_Born
logical      :: freeze
real(8)						::	fk_wsphere, crgtot, crgQtot
integer(AI), allocatable	::	list_sh(:,:), nsort(:,:)
real(8),  allocatable		::	theta(:), theta0(:), tdum(:)
integer						::	nwpolr_shell, n_max_insh


type SHELL_TYPE
        real					::	rout, dr, cstb
        real					::	avtheta, avn_insh, theta_corr
        integer					::	n_insh
end type SHELL_TYPE


type(SHELL_TYPE), allocatable::	wshell(:)


! constants & default values
integer, parameter			::	itdis_update		= 100
real,  parameter			::	wpolr_layer			= 3.0001
real, parameter				::	drout				= 0.5
real(8), parameter			::	tau_T_default		= 100.
real(8), parameter			::	rcpp_default		= 10.
real(8), parameter			::	rcww_default		= 10.
real(8), parameter			::	rcpw_default		= 10.
real(8), parameter			::	rcq_default			= 99.
real(8), parameter			::	rcLRF_default		= 99.
!recxl_i is set to rexcl_o * shell_default
real(8), parameter   :: fk_fix_default = 200.
real(8), parameter			::	shell_default		= 0.85
real(8), parameter			::	fk_pshell_default	= 10.
integer, parameter			::	itrj_cycle_default	= 100
integer, parameter			::	iene_cycle_default	= 10
integer, parameter			::	iout_cycle_default	= 10
integer, parameter			::	nb_cycle_default	= 10
real(8), parameter			::	fkwpol_default		= 20.
real(8), parameter			::	fk_wsphere_default	= 60.
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

! --- Restraints
integer						::	implicit_rstr_from_file
integer      :: nrstr_seq, nrstr_pos, nrstr_dist, nrstr_angl, nrstr_wall


type RSTRSEQ_TYPE
        integer(AI)				::	i,j
        real(8)					::	fk
        integer(TINY)			::	ih
        integer					::	to_centre !flag for restraining to geom. or mass centre
end type RSTRSEQ_TYPE


type RSTRPOS_TYPE
        integer(AI)				::	i
        integer(TINY)			::	ipsi
        real(8)					::	fk(3)
        real(8)					::	x(3)
end type RSTRPOS_TYPE


type RSTRDIS_TYPE
        integer(AI)				::	i,j
        integer(TINY)			::	ipsi
        real(8)					::	fk
        real(8)					::	d1, d2
        character(len=20)       ::  itext,jtext
end type RSTRDIS_TYPE

type RSTRANG_TYPE 
        integer(AI)    :: i,j,k
        integer(TINY)  :: ipsi
        real(8)        :: fk
        real(8)        :: ang
!       character(len=20)       ::  itext,jtext,ktext
end type RSTRANG_TYPE                                                                                                

type RSTRWAL_TYPE

        integer(AI)				::	i,j
        real(8)					::	d, fk, aMorse, dMorse
        integer(TINY)			::	ih
end type RSTRWAL_TYPE


type(RSTRSEQ_TYPE), allocatable::	rstseq(:)
type(RSTRPOS_TYPE), allocatable::	rstpos(:)
type(RSTRDIS_TYPE), allocatable::	rstdis(:)
type(RSTRANG_TYPE), allocatable:: rstang(:)
type(RSTRWAL_TYPE), allocatable::	rstwal(:)


!-----------------------------------------------------------------------
!	Coordinates, velocities, forces
!-----------------------------------------------------------------------
real(8), allocatable		::	d(:)
real(8), allocatable		::	x(:)
real(8), allocatable		::	xx(:) !for shake
real(8), allocatable		::	v(:)
real(8), allocatable			::	winv(:)
real(8)						::	grms !RMS force


!-----------------------------------------------------------------------
!	Energies , EQ is defined in qatom.f90
!-----------------------------------------------------------------------
type(ENERGIES)				::	E
real(8)						::	Tfree, Tfree_solvent, Tfree_solute, Temp
real(8)						::	Temp_solvent, Temp_solute, Texcl_solute, Texcl_solvent


!-----------------------------------------------------------------------
!	Nonbonded pair information
!-----------------------------------------------------------------------
type NB_TYPE
        integer(AI)				::	i, j
        integer(TINY)			::	LJcod
        integer					:: cgp_pair ! cgp_pair only used with periodic conditions
end type NB_TYPE

type CGP_PAIR_TYPE
        integer(AI)				:: i, j !switching atoms (or equal in case of no switching atoms) of the chargegroups
        real(8)					:: x, y, z !periodical shifts
end type CGP_PAIR_TYPE

type NBQP_TYPE
        integer(AI)				::	i, j
        integer(TINY)			::	LJcod, qLJcod
        integer					::	cgp_pair !this variable only used with periodic conditions
end type NBQP_TYPE


type NBQ_TYPE
        integer(AI)				::	j       !atom number
        integer(AI)				::	iq, jq  !q-atom numbers
        integer(TINY)			::	LJcod
        real(8)                 ::  el_scale !scale factor for electostatic interactions in qq-pairs
end type NBQ_TYPE

integer						::	nbpp_pair !current no solute-solute pairs
type(NB_TYPE), allocatable, target::nbpp(:)

integer	::	nbww_pair,nbww_true_pair !current no solvent-solvent pairs, implicit and explicit
integer(AI), allocatable, target::nbww(:)

integer						::	nbpw_pair !current no solute-solvent pairs
type(NB_TYPE), allocatable, target::nbpw(:)

integer						::	nbqq_max !max number of q-q pairs in any state
integer(TINY), allocatable	::	qconn(:,:,:) !Q-atom connectivity list
integer						::	nbqq_pair(max_states)
type(NBQ_TYPE), allocatable ::nbqq(:,:)

integer						::	nbqp_pair !current no of qatom-solute pairs
type(NBQP_TYPE), allocatable, target::nbqp(:)

integer						::	nbqw_pair !current no of q-atom-water mol. pairs
integer(AI), allocatable	::	nbqw(:)


!these three used only under periodic conditions
integer						:: nbpp_cgp_pair !number of solute-solute chargegroups interacting
type(CGP_PAIR_TYPE), allocatable:: nbpp_cgp(:)
integer						:: nbpw_cgp_pair
type(CGP_PAIR_TYPE), allocatable :: nbpw_cgp(:)
integer						::	nbqp_cgp_pair
type(CGP_PAIR_TYPE), allocatable :: nbqp_cgp(:)

!special monitoring of pairs
integer (TINY),allocatable  :: special_LJcod(:,:,:,:)

! LRF related variables
integer(AI), allocatable	::	iwhich_cgp(:)


type LRF_TYPE
        real(8)					::	cgp_cent(3)
        real(8)					::	phi0
        real(8)					::	phi1(3)
        real(8)					::	phi2(9)
        real(8)					::	phi3(27)
end type LRF_TYPE

type(LRF_TYPE), allocatable	::	lrf(:)
type(LRF_TYPE), allocatable	::	old_lrf(:)  !for constant pressure: MC_volume routine


type(node_assignment_type)	:: calculation_assignment


!shake types & variables
!convergence criterion (fraction of distance)
real(8), parameter			::	SHAKE_TOL = 0.0001
integer, parameter			::	SHAKE_MAX_ITER = 1000


type SHAKE_BOND_TYPE
        integer(AI)				::	i,j
        real(8)					::	dist2
        logical					::	ready
end type SHAKE_BOND_TYPE


type SHAKE_MOL_TYPE
        integer					::	nconstraints
        type(SHAKE_BOND_TYPE), pointer :: bond(:)
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
	real(8)			:: time = 0.0
end type profiling_var_type

type(profiling_var_type)	:: profile(num_profiling_times)

#if defined (USE_MPI)
 !vectors for keeping track of node times
real(8),allocatable				:: all_node_times(:)
real(8),allocatable				:: node_times(:)
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
pi = 4.0*atan(1.0)
deg2rad = pi/180.0

end subroutine md_startup


!----------------------------------------------------------------------

subroutine md_shutdown
! call used modules' shutdown subroutines
call md_deallocate
call topo_deallocate
call qatom_shutdown
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

subroutine allocate_natom_arrays

allocate(x(natom*3), &
                xx(natom*3), &
                v(natom*3), &
                d(natom*3), &
                winv(natom), &
                iqatom(natom), &
                stat=alloc_status)
call check_alloc('atom data arrays')
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
         request_recv(numnodes-1,4), &
         d_recv(natom*3,numnodes-1), &
         E_recv(numnodes-1), &
         EQ_recv(nstates,numnodes-1), &
         stat=alloc_status)
 call check_alloc('MPI data arrays')
 if (use_excluded_groups) then
	allocate(EQ_gc_recv(ngroups_gc,nstates,numnodes-1), &
        stat=alloc_status)
	call check_alloc('MPI data array excluded groups')
 end if
else
 allocate(E_send(1), &
          EQ_send(nstates), &
          stat=alloc_status)
 call check_alloc('MPI data arrays')
 if (use_excluded_groups) then
	allocate(EQ_gc_send(ngroups_gc,nstates), &
	stat=alloc_status)
	call check_alloc('MPI data array excluded groups')
 end if
end if
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

! use status to avoid errors if not allocated

! atom arrays
deallocate (x, stat=alloc_status)
deallocate (xx, stat=alloc_status)
deallocate (v, stat=alloc_status)
deallocate (d, stat=alloc_status)
deallocate (winv, stat=alloc_status)
deallocate (iqatom, stat=alloc_status)

! Nosé-Hoover array
if ( thermostat == NOSEHOOVER) then
deallocate (xnh, stat=alloc_status)
deallocate (vnh, stat=alloc_status)
deallocate (qnh, stat=alloc_status)
deallocate (Gnh, stat=alloc_status)
end if
! nonbond lists
deallocate(nbpp, nbpw, nbww, nbqq, nbqp, nbqw, qconn, stat=alloc_status)

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
   if (allocated(EQ_gc_recv)) then
	deallocate(EQ_gc_recv,stat=alloc_status)
   end if
else
   deallocate(E_send, stat=alloc_status)
   deallocate(EQ_send, stat=alloc_status)
   if (allocated(EQ_gc_send)) then
	deallocate(EQ_gc_send,stat=alloc_status)
   end if
end if
#endif

end subroutine md_deallocate

!-----------------------------------------------------------------------

subroutine reallocate_nonbondlist_pp

! variables
type(NB_TYPE), allocatable	:: old_nbxx(:)
integer						:: old_max,temp_max

! copy
!make it safer and add some padding
!memory is cheap
old_max = calculation_assignment%pp%max
temp_max = int(calculation_assignment%pp%max * 1.05) + 200
allocate(old_nbxx(temp_max), stat=alloc_status)
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
!making some padding already here
new_max = int( old_max*1.05) + 200
allocate(old_nbxx(new_max), stat=alloc_status)
call check_alloc('reallocating non-bonded charge group list')
old_nbxx(1:old_max) = nbpp_cgp(1:old_max)

! deallocate and copy back
deallocate(nbpp_cgp)
!new_max = int( old_max*1.05) + 200
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
!again, giving some padding
new_max = int( old_max*1.05) + 200
allocate(old_nbxx(new_max), stat=alloc_status)
call check_alloc('reallocating non-bonded charge group list')
old_nbxx(1:old_max) = nbpw_cgp(1:old_max)

! deallocate and copy back
deallocate(nbpw_cgp)
!new_max = int( old_max*1.05) + 200
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
new_max = int( old_max*1.05) + 200
allocate(old_nbxx(new_max), stat=alloc_status)
call check_alloc('reallocating non-bonded charge group list')
old_nbxx(1:old_max) = nbqp_cgp(1:old_max)

! deallocate and copy back
deallocate(nbqp_cgp)
!new_max = int( old_max*1.05) + 200
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
integer						:: old_max,temp_max

! copy
old_max = calculation_assignment%pw%max
temp_max = int(calculation_assignment%pw%max * 1.05) + 200
allocate(old_nbxx(temp_max), stat=alloc_status)
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
type(NBQP_TYPE), allocatable	:: old_nbxx(:)
integer						:: old_max,temp_max

! copy
old_max = calculation_assignment%qp%max
temp_max = int(old_max * 1.05 ) +200
allocate(old_nbxx(temp_max), stat=alloc_status)
call check_alloc('reallocating non-bonded interaction list')
old_nbxx(1:old_max) = nbqp(1:old_max)

! deallocate and copy back
deallocate(nbqp)
calculation_assignment%qp%max = int(calculation_assignment%qp%max * 1.05) + 200
allocate(nbqp(calculation_assignment%qp%max), stat = alloc_status)
call check_alloc('reallocating non-bonded interaction list')
nbqp(1:old_max) = old_nbxx(1:old_max)

! deallocate copy
deallocate(old_nbxx)

! tell the world
write (*,100) calculation_assignment%qp%max
100 format('>>> reallocating q-s pair list, new max is ', i8)

end subroutine reallocate_nonbondlist_qp

!----------------------------------------------------------------------

subroutine reallocate_nonbondlist_ww
! variables
integer(AI), allocatable		:: old_nbxx(:)
integer						:: old_max,temp_max

! copy
old_max = calculation_assignment%ww%max
temp_max = int(old_max * 1.05 ) + 200
allocate(old_nbxx(temp_max), stat=alloc_status)
call check_alloc('reallocating non-bonded interaction list')
old_nbxx(1:old_max) = nbww(1:old_max)
!old_max = calculation_assignment%ww%max

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

real(8) function angle(istart, iend)
! *** arguments
integer						::	istart, iend

! *** local variables
integer						::	i,j,k,ia,ic,i3,j3,k3
real(8)						::	bjiinv, bjkinv, bji2inv, bjk2inv
real(8)						::	scp,angv,da,dv,f1
real(8)						::  rji(3),rjk(3),di(3),dk(3) 

! global variables used:
! ang, x, anglib, d

! calculate the total energy of all protein or water angles, depending
! updates d

! reset Eangle
angle = 0.

do ia=istart,iend
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
bji2inv = 1./(rji(1)**2 + rji(2)**2 + rji(3)**2 )
bjk2inv = 1./(rjk(1)**2 + rjk(2)**2 + rjk(3)**2 )
        bjiinv = sqrt(bji2inv)
        bjkinv = sqrt(bjk2inv)

        ! calculate scp and angv
scp = ( rji(1)*rjk(1) + rji(2)*rjk(2) + rji(3)*rjk(3) )
scp = scp * bjiinv*bjkinv
if ( scp .gt.  1.0 ) then
          scp =  1.0
        else if ( scp .lt. -1.0 ) then
          scp = -1.0
        end if
angv = acos(scp)

        ! calculate da and dv
da = angv - anglib(ic)%ang0
angle = angle + 0.5*anglib(ic)%fk*da**2
dv = anglib(ic)%fk*da

        ! calculate f1
f1 = sin ( angv ) 
if ( abs(f1) .lt. 1.e-12 ) then
          ! avoid division by zero
          f1 = -1.e12
        else
  f1 =  -1.0 / f1
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

end function angle

!-----------------------------------------------------------------------
real(8) function urey_bradley(istart, iend)
! *** arguments
integer						::	istart, iend

! *** local variables
integer						::	i,j,k,ia,ic,i3,j3,k3
real(8)						::	bjiinv, bjkinv, bji2inv, bjk2inv
real(8)						::	scp,angv,da,dv,f1
real(8)						::  rji(3),rjk(3),di(3),dk(3) 
real(8)						::	rik(3), dik, ru, du 
real(8)						::	Eurey

! global variables used:
! ang, x, anglib, d

! reset energy
urey_bradley = 0.

do ia=istart,iend
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
    if(anglib(ic)%ureyfk > 0.) then
                rik(1) = x(k3+1) - x(i3+1)
                rik(2) = x(k3+2) - x(i3+2)
                rik(3) = x(k3+3) - x(i3+3)
                dik = sqrt(rik(1)*rik(1) + rik(2)*rik(2) + rik(3)*rik(3))
                ru = dik - anglib(ic)%ureyr0
                urey_bradley = urey_bradley + anglib(ic)%ureyfk*ru**2
                du = 2*anglib(ic)%ureyfk*ru/dik
                d(k3+1) = d(k3+1) + du*rik(1)
                d(k3+2) = d(k3+2) + du*rik(2)
                d(k3+3) = d(k3+3) + du*rik(3)
                d(i3+1) = d(i3+1) - du*rik(1)
                d(i3+2) = d(i3+2) - du*rik(2)
                d(i3+3) = d(i3+3) - du*rik(3)
        end if
end do

end function urey_bradley

!-----------------------------------------------------------------------

real(8) function bond(istart, iend)
! *** arguments
integer						::	istart, iend

! *** local variables
integer						::	i,j,ib,ic,i3,j3
real(8)						::	b,db,dv
real(8)						::	rij(3)

! global variables used:
! bnd, x, bondlib, d

! reset Ebond
bond = 0

do ib=istart,iend
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
        bond = bond + 0.5*bondlib(ic)%fk*db**2

        ! calculate dv and update d
        dv = bondlib(ic)%fk*db/b
        d(j3+1) = d(j3+1) + rij(1)*dv
        d(j3+2) = d(j3+2) + rij(2)*dv
        d(j3+3) = d(j3+3) + rij(3)*dv
        d(i3+1) = d(i3+1) - rij(1)*dv
        d(i3+2) = d(i3+2) - rij(2)*dv
        d(i3+3) = d(i3+3) - rij(3)*dv
end do

end function bond
!-----------------------------------------------------------------------
subroutine cgp_centers
! *** local variables
integer						::	ig,i,i3

do ig = 1, ncgp
        lrf(ig)%cgp_cent(:) = 0.
        lrf(ig)%phi0 = 0.
lrf(ig)%phi1(:) = 0.
lrf(ig)%phi2(:) = 0.
lrf(ig)%phi3(:) = 0.

        do i  = cgp(ig)%first, cgp(ig)%last
                lrf(ig)%cgp_cent(:) = lrf(ig)%cgp_cent(:) + x(cgpatom(i)*3-2:cgpatom(i)*3)
        end do

lrf(ig)%cgp_cent(:) = lrf(ig)%cgp_cent(:)/real(cgp(ig)%last - cgp(ig)%first +1)

end do

end subroutine cgp_centers
!-----------------------------------------------------------------------
subroutine make_nbqqlist
!locals
integer						::	is

call make_qconn
nbqq_max = nbqq_count()
allocate(nbqq(nbqq_max, nstates), stat=alloc_status)
call check_alloc('Qatom-Qatom non-bond list')
! prepare q-atom nonbond lists that do not need updating
call nbqqlist

do is =1, nstates
        write (*,200) nbqq_pair(is),is
end do
write (*,*)

200 format ('No. of Rcq indep. nb pairs involving q-atoms = ',i5, &
' in state :',i3)
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
totnbww = nww*9
totnbqp = nqp
totnbqw = nqw*3*nqat

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
percent=REAL(totnbpp)/n_nonbonded
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
percent=REAL(totnbpw)/n_nonbonded
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
percent=REAL(totnbqp)/n_nonbonded
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
percent=REAL(totnbww)/n_nonbonded
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
percent=REAL(totnbqw)/n_nonbonded
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

call MPI_Type_contiguous(5, mpitype_pair_assignment, mpitype_node_assignment, ierr)
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
 write(*,98) 'solute-solute', 'solute-water', 'water-water', 'Q-solute', 'Q-water'
 write(*,99) 'total', ncgp_solute,ncgp_solute,nwat,ncgp_solute,nwat
if (numnodes .gt. 1) then
 do i=0,numnodes-1
    write(*,100) i, 'assigned cgps', &
         node_assignment(i)%pp%end-node_assignment(i)%pp%start+1, &
         node_assignment(i)%pw%end-node_assignment(i)%pw%start+1, &
         node_assignment(i)%ww%end-node_assignment(i)%ww%start+1, &
         node_assignment(i)%qp%end-node_assignment(i)%qp%start+1, &
         node_assignment(i)%qw%end-node_assignment(i)%qw%start+1
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
allocate(nbqp(calculation_assignment%qp%max), stat=alloc_status)
call check_alloc('Qatom - solute non-bond list')

calculation_assignment%qw%max = nwat
allocate(nbqw(calculation_assignment%qw%max), stat=alloc_status)
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
  end if	

!Kanske deallokera nbxx_per_cgp TODO

98 format('node value ',5a13)
99 format(a10,1x,5(1x,i12))
!99 format(a4,2x,a,t18,i13,3x,i13,3x,i13,3x,i13)
100 format(i4,1x,a5,1x,5(1x,i12))
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
integer :: iii
close (3)
if ( itrj_cycle .gt. 0 ) close (10)
if ( iene_cycle .gt. 0 ) close (11)

if (use_excluded_groups) then
        do iii=1,ngroups_gc
        close (ST_gc(iii)%fileunit)
        end do
end if
end subroutine close_output_files

!-----------------------------------------------------------------------

subroutine open_files
	integer :: iii
! --> restart file (2)
if(restart) then
open (unit=2, file=restart_file, status='old', form='unformatted', action='read', err=2)
end if

! --> final coords (3)
open (unit=3, file=xfin_file, status='unknown', form='unformatted', action='write', err=3)

! --> energy output file (11)
if ( iene_cycle .gt. 0 ) then
open (unit=11, file=ene_file, status='unknown', form='unformatted', action='write', err=11)
	if (use_excluded_groups) then
		do iii=1,ngroups_gc
                ST_gc(iii)%fileunit=freefile()
		open (unit=ST_gc(iii)%fileunit,file=ST_gc(iii)%filename,status='unknown',form='unformatted',action='write',err=11)
		end do
	end if
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
real(8)						::	fk,r2,erst
real(8)						::	dr(3)

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
erst    = 0.5*fk*r2

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
real(8)					::	am,sd,v
integer					::	ig

! local variables
integer					::	i
real(8)					::	a,y

a=0.0
do i=1,12
y=randm(ig)
a=a+y
end do
v=(a-6.0)*sd+am
end subroutine gauss

!-----------------------------------------------------------------------
subroutine get_fep
! local variables
character					::	libtext*80,qaname*2
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
sc_lookup(:,:,:)=0.0

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
!			if(( (tor(i)%i.eq.iqtor(j) .and. tor(i)%j.eq.jqtor(j) .and. &
!			        tor(i)%k.eq.kqtor(j) .and. tor(i)%l.eq.lqtor(j)) .or. &
!                               (tor(i)%i.eq.lqtor(j) .and. tor(i)%j.eq.kqtor(j) .and. &
!				tor(i)%k.eq.jqtor(j) .and. tor(i)%l.eq.iqtor(j)) ) .and. &
!                                tor(i)%cod /= 0) then
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
!                                if(((imp(i)%i .eq. iqimp(j)) .or. &
!                                        (imp(i)%i .eq. lqimp(j)) .or. &
!                                        (imp(i)%l .eq. iqimp(j)) .or. &
!                                        (imp(i)%l .eq. lqimp(j))) .and. &
!                                        ((imp(i)%j .eq. iqimp(j)) .or. &
!                                        (imp(i)%j .eq. jqimp(j))  .or. &
!                                        (imp(i)%j .eq. kqimp(j))  .or. &
!                                        (imp(i)%j .eq. lqimp(j))) .and. &
!                                        ((imp(i)%k .eq. iqimp(j)) .or. &
!                                        (imp(i)%k .eq. jqimp(j)) .or. &
!                                        (imp(i)%k .eq. kqimp(j)) .or. &
!                                        (imp(i)%k .eq. lqimp(j))) .and. &
!                                        imp(i)%cod /= 0) then
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
!                        if(((imp(i)%j.eq.jqimp(j) .and. imp(i)%k.eq.kqimp(j)) .or. &
!                                (imp(i)%j.eq.kqimp(j) .and. imp(i)%k.eq.jqimp(j))) .and. &
!                                imp(i)%cod ne. 0) then
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
character					::	text*80,filnam*80
integer					::	length
! local variables


integer					::	i

length=80
do i=1,80
if ( text(i:i) .eq. ' ' ) then
length=i-1
goto 10
end if
end do
10 filnam(1:length)=text(1:length)

end subroutine get_fname

!-----------------------------------------------------------------------

real(8) function improper(istart, iend)
!arguments
integer						::	istart, iend

! evaluate harmonic impropers
! local variables
integer						::	ip
real(8)						::	scp,phi,dv,arg,f1
real(8)						::	bjinv, bkinv, bj2inv, bk2inv
real(8)						::	rji(3),rjk(3),rkl(3),rnj(3),rnk(3)
real(8)						::	rki(3),rlj(3),dp(12),di(3),dl(3)
type(TOR_TYPE), pointer		::	t
type(IMPLIB_TYPE), pointer	::	lib

! global variables used:
!  imp, implib, x, pi, d

improper = 0.

do ip = iend, istart,-1
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

bj2inv  = 1./( rnj(1)**2 + rnj(2)**2 + rnj(3)**2)
bk2inv  = 1./( rnk(1)**2 + rnk(2)**2 + rnk(3)**2)
        bjinv = sqrt(bj2inv)
        bkinv = sqrt(bk2inv)

scp = (rnj(1)*rnk(1)+rnj(2)*rnk(2)+rnj(3)*rnk(3))*(bjinv*bkinv)
if ( scp .gt.  1.0 ) then
                scp =  1.0
else if ( scp .lt. -1.0 ) then 
                scp = -1.0
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
arg = arg - 2.*pi*nint(arg/(2.*pi))
dv  = lib%fk*arg
improper = improper + 0.5*dv*arg

! ---       forces

f1 = sin ( phi ) 
if ( abs(f1) .lt. 1.e-12 ) f1 = 1.e-12
f1 =  -1.0 / f1
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
end function improper

!-----------------------------------------------------------------------

real(8) function improper2(istart, iend)
!evaluate periodic impropers
!arguments
integer						::	istart, iend
! local variables
integer						::	ip
real(8)						::	scp,phi,dv,arg,f1
real(8)						::	bjinv, bkinv, bj2inv, bk2inv
real(8)						::	rji(3),rjk(3),rkl(3),rnj(3),rnk(3)
real(8)						::	rki(3),rlj(3),dp(12),di(3),dl(3)
type(TOR_TYPE), pointer		::	t
type(IMPLIB_TYPE), pointer	::	lib

! global variables used:
! imp, implib, x, pi, d

improper2 = 0.

do ip = iend, istart,-1
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


bj2inv  = 1./( rnj(1)**2 + rnj(2)**2 + rnj(3)**2)
bk2inv  = 1./( rnk(1)**2 + rnk(2)**2 + rnk(3)**2)
        bjinv = sqrt(bj2inv)
        bkinv = sqrt(bk2inv)

scp = (rnj(1)*rnk(1)+rnj(2)*rnk(2)+rnj(3)*rnk(3))*(bjinv*bkinv)
if ( scp .gt.  1.0 ) then
                scp =  1.0
else if ( scp .lt. -1.0 ) then 
                scp = -1.0
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
improper2 = improper2 + lib%fk * (1 + cos(arg))
dv  = -2*lib%fk * sin(arg)

! ---       forces

f1 = sin ( phi ) 
if ( abs(f1) .lt. 1.e-12 ) f1 = 1.e-12
f1 =  -1.0 / f1
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
integer					:: nat3,j,jj
real(kind=wp8), allocatable		:: temp_lambda(:)
integer, parameter                      ::maxint=2147483647
real(kind=wp8), parameter                        ::maxreal=1E35
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



! water parameters: crg_ow, crg_hw (used by nonbond_ww)
call MPI_Bcast(crg_ow, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast crg_ow')
call MPI_Bcast(crg_hw, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast crg_hw')

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

!**MN-> Behövs om shake ska parallelliseras
! shake/temperature parameters
!call MPI_Bcast(shake_constraints, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)  !bara i init_shake & md_run
!call MPI_Bcast(shake_molecules, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)    !bara i div init_
!call MPI_Bcast(Ndegf, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)    !bara i div init_
!call MPI_Bcast(Ndegfree, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)    !bara i div init_

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
call MPI_Bcast(nat_pro, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast nat_pro')


!vars from QATOM
call MPI_Bcast(nstates, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast nstates')
call MPI_Bcast(use_excluded_groups,1,MPI_LOGICAL, 0,MPI_COMM_WORLD,ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast use_excluded_groups')

!first time force syncing so that everyone knows and has all parameters to this
!point, should give no hits to performance
!Paul Bauer 2014
call MPI_BARRIER(MPI_COMM_WORLD, ierr)



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
                tempmask(:,jj)=ST_gc(jj)%gcmask%mask(:)
                end do
        end if
        call MPI_Bcast(tempmask,nat_pro*ngroups_gc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        if (ierr .ne. 0) call die('init_nodes/MPI_Bcast gc atom masks')
        if (nodeid.ne.0) then
                do jj=1,ngroups_gc
                call mask_initialize(ST_gc(jj)%gcmask)
                ST_gc(jj)%gcmask%mask(:)=tempmask(:,jj)
                end do
        end if
        deallocate(tempmask,stat=alloc_status)
end if


!Setting all vars not sent to slaves to 2147483647. To avoid hidden bugs.
if (nodeid .ne. 0) then 
shake_constraints=maxint 
shake_molecules=maxint 
Ndegf=maxint 
Ndegfree=maxint
xwcent(:)=maxreal 
end if

!first part completed, syncing again
call MPI_BARRIER(MPI_COMM_WORLD, ierr)

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

!Setting all vars not sent to slaves to 2147483647. To avoid conflicts.
!winv is an array of type real8 -> set to maxreal
if (nodeid .ne. 0) then 
winv(:)=maxreal
end if

!Broadcast iqatom
call MPI_Bcast(iqatom, natom, MPI_INTEGER2, 0, MPI_COMM_WORLD, ierr) !(TINY)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast iqatom')

!Broadcast ljcod
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




! --- shake data ---

!if (shake_solute .or. shake_solvent .or. shake_hydrogens) then
! shake stuff

!if (nodeid .eq. 0) write (*,'(80a)') 'shake data'

!if (nodeid .ne. 0) then
! allocate shake arrays
! ADD CODE HERE to allocate shake array!
!end if

! bake all shake data into a big packet & bcast
! ADD CODE HERE to broadcast shake data
!end if

!and syncing again
call MPI_BARRIER(MPI_COMM_WORLD, ierr)
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
ftype(:) = MPI_REAL8
blockcnt(1) = 3					! real(8) cgp_cent(3)
fdisp(1) = 0
blockcnt(2) = 1					! real(8) phi0
fdisp(2) = 3*8
blockcnt(3) = 3					! real(8) phi1(3)
fdisp(3) = 3*8 + 8
blockcnt(4) = 9					! real(8) phi2(9)
fdisp(4) = 3*8 + 8 + 3*8
blockcnt(5) = 27				! real(8) phi3(27)
fdisp(5) = 3*8 + 8 + 3*8 + 9*8
call MPI_Type_create_struct(5, blockcnt, fdisp, ftype, mpitype_batch, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Type_create_struct')
call MPI_Type_commit(mpitype_batch, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Type_commit')
call MPI_Bcast(lrf, ncgp, mpitype_batch, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast shake parameters')
call MPI_Type_free(mpitype_batch, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Type_free')
end if !(use_LRF)

!sync to make sure everything is there
call MPI_BARRIER(MPI_COMM_WORLD, ierr)
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
blockcnt(1) = 1					! real(8) mass
fdisp(1) = 0
blockcnt(2) = nljtyp				! real(8) avdw(nljtyp)
fdisp(2) = 8
blockcnt(3) = nljtyp				! real(8) bvdw(nljtyp)
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
!end of part -> sync
call MPI_BARRIER(MPI_COMM_WORLD, ierr)

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
call check_alloc('Q-atom arrays')
	if (use_excluded_groups) then
		allocate(EQ_gc(ngroups_gc,nstates),stat=alloc_status)
		call check_alloc('Q-energy arrays excluded groups')
	end if
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
! qavdw and qbvdw share the same format: real(8) qxvdw(nqlib,nljtyp)
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

if (nodeid .eq. 0) then 
call centered_heading('End of initiation', '-')
print *
end if

!and we sync again
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
integer						::	mol, b, ia, ja, constr, angle
real(8)						:: exclshk
integer						::	src, trg
integer						::	solute_shake_constraints

!allocate molecule list
allocate(shake_mol(nmol), stat=alloc_status)
call check_alloc('shake molecule array')

shake_mol(:)%nconstraints = 0
mol = 0
exclshk = 0.

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
        if((shake_hydrogens .and. (.not. heavy(ia) .or. .not. heavy(ja))) .or. &
           (shake_solute .and. ia <= nat_solute) .or. &
           (shake_solvent .and. ia > nat_solute)) then
           shake_mol(mol)%nconstraints = shake_mol(mol)%nconstraints + 1

           if( .not. use_PBC ) then
                if(excl(ia)) exclshk = exclshk + 0.5
                if(excl(ja)) exclshk = exclshk + 0.5
           end if

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
        ia = bnd(b)%i
        ja = bnd(b)%j
        do while(ia >= istart_mol(mol+1)) 
                !new molecule
                mol = mol +1
                shake_mol(mol)%nconstraints = 0
        end do
        !skip redefined bonds
        if(bnd(b)%cod == 0) cycle
if((shake_hydrogens .and. (.not. heavy(ia) .or. .not. heavy(ja))) .or.&
           (shake_solute .and. ia <= nat_solute) .or. &
           (shake_solvent .and. ia > nat_solute)) then
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
                nullify(shake_mol(src)%bond) 
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
character					:: text*80
integer						:: i,j,length,ii
real(8)						:: stepsize
real(8)						:: lamda_tmp(max_states)
integer						:: fu, fstat
real(8)						::	rjunk
integer						::	ijunk

! local parameters
integer						:: num_args
character(200)				:: infilename
logical						::	yes
logical						::	need_restart
character(len=80)			::	instring
logical						::	inlog
integer						::	mask_rows, number
integer,parameter				:: maxmaskrows=30
character*80					:: gc_mask_tmp(maxmaskrows),str,str2

! this subroutine will init:
!  nsteps, stepsize, dt
!  Temp0, tau_T, iseed, Tmaxw
!  use_LRF, NBcycle, Rcpp, Rcww, Rcpw, Rcq
!  shake_solute, shake_solvent, shake_hydrogens
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
                        pressure = 1.0
                end if
                write(*,9) pressure
9	format ('Pressure = ',f10.3,'  bar')
                !convert pressure to strange internal unit
                pressure = pressure * 1.43836e-5
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
dt=0.020462*stepsize
dt2=0.5*dt
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
tau_T=0.020462*tau_T
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
		gkT = 2*friction*Boltz*Temp0/dt !constant to be used to generate the random forces of the thermostat
		write(*,*) 'Langevin thermostat friction constant set to default: 1/tau_T'
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
write(*,17) 'all solvent bonds', onoff(shake_solvent)
17	format('Shake constaints on ',a,t42,': ',a3)

if(.not. prm_get_logical_by_key('shake_solute', shake_solute, .false.)) then
        write(*,'(a)') '>>> Error: shake_solute must be on or off.'
        initialize = .false.
end if 
write(*,17) 'all solute bonds', onoff(shake_solute)

if(.not. prm_get_logical_by_key('shake_hydrogens', shake_hydrogens, .false.)) then
        write(*,'(a)') '>>> Error: shake_hydrogens must be on or off.'
        initialize = .false.
end if 
write(*,17) 'all bonds to hydrogen', onoff(shake_hydrogens)


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
if(.not. prm_open_section('cut-offs')) then
        write(*,'(a)') 'No cut-offs section, default cut-offs used'
        rcpp = rcpp_default
        rcww = rcww_default
        rcpw = rcpw_default
        rcq = rcq_default
        rcLRF = rcLRF_default
else
        if(.not. prm_get_real8_by_key('solute_solute', rcpp, rcpp_default)) then
                write(*,'(a)') 'solute-solute cut-off set to default'
        end if
        if(.not. prm_get_real8_by_key('solvent_solvent', rcww, rcww_default)) then
                write(*,'(a)') 'solvent-solvent cut-off set to default'
        end if
        if(.not. prm_get_real8_by_key('solute_solvent', rcpw, rcpw_default)) then
                write(*,'(a)') 'solute-solvent cut-off set to default'
        end if
        if(.not. prm_get_real8_by_key('q_atom', rcq, rcq_default)) then
                write(*,'(a)') 'q-atom cut-off set to default'
        end if
        if(use_LRF) then
                if(.not. prm_get_real8_by_key('lrf', rcLRF, rcLRF_default)) then
                        write(*,'(a)') 'LRF cut-off set to default'
                end if
                if(RcLRF < rcpp .or. RcLRF < rcpw .or. RcLRF < rcww) then
                        write(*,'(a)') &
                                '>>> ERROR; LRF cut-off must not be smaller than solute or solvent cut-offs!'
                        initialize = .false.
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

30	format ('>>> WARNING: Ingnoring obsolete keyword ',a,'.')
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
                    if(rexcl_i < 0.) then
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
                lamda_tmp(1) = 1.
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
ngroups_gc = prm_count('exclude_groups')
write (*,'(/,a)') 'List of excluded atom groups'
if ( ngroups_gc .gt. 0 ) then
!allocate new EQ arrays for each group
777 format (/,'No of excluded group calculations =',i10)
	write (*,777) ngroups_gc
	allocate(ST_gc(ngroups_gc), stat=alloc_status)
	call check_alloc('gc param store')
	allocate(EQ_gc(ngroups_gc,nstates), stat=alloc_status)
	call check_alloc('gc excluded list')
778 format ('groupno	contains')
	write (*,778)
	do i=1,ngroups_gc
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

		mask_rows = number - 1
		ST_gc(i)%seltype = trim(gc_mask_tmp(1))
		if ((trim(ST_gc(i)%seltype) /= 'atom' ).and. (trim(ST_gc(i)%seltype) /= 'residue')) then
779 format (/,'Could not understand name of the selection , ',a)
			write (*,779) text
                	initialize = .false.
		end if
		if(iene_cycle > 0) then
        	if(mask_rows == 0) then
!You are not as funny as you think ...
                	write(*,780) i
780 format (/,'No atoms excluded for group ',i10)
                	ST_gc(i)%count=0
        	else
!			call mask_initialize(ST_gc(i)%gcmask)
			ST_gc(i)%count=mask_rows
			allocate(ST_gc(i)%maskarray(mask_rows),stat=alloc_status)
			call check_alloc('gc maskarray')
			
	                do ii=2,mask_rows+1
                	        yes = gc_store_mask(ST_gc(i),gc_mask_tmp(ii))
	                end do
			use_excluded_groups = .true.
!			ST_gc(i)%fileunit=freefile()
			write (str,'(i10)') i
			ST_gc(i)%filename=trim(adjustl(str))//'_'
                        str2=ene_file
                        str=trim(ST_gc(i)%filename)//trim(str2)
			ST_gc(i)%filename=trim(str)
781 format (/,i10,i10,' atom groups')
                        write(*,781) i,mask_rows 
	        end if
		elseif(mask_rows == 0) then
		        write(*,'(a)') 'Ignoring section for excluding atoms.'
		end if
	end do

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
character					:: text*80, watmodel*80
integer						:: i,j,length
integer						:: irestart
real(8)						:: stepsize
real(8)						:: lamda_tmp(max_states)
integer						:: fstat
integer						::	NBMethod
integer						::	iwpol_restr
real(8)						::	rjunk

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
RcLRF = 999.
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
dt=0.020462*stepsize
dt2=0.5*dt

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
tau_T=0.020462*tau_T
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

subroutine lrf_taylor
! *** local variables
integer						::	i,i3,ic
real(8)						::	Vij, q
real(8)						::	dr(3),df(3)

! global variables used:
!  E%LRF, natom, excl, iqatom, iwhich_cgp, lrf, x, crg, d

do i = 1, natom
! for every atom:

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
     +0.5*(lrf(ic)%phi2(1)*dr(1)+lrf(ic)%phi2(2)*dr(2) &
     +lrf(ic)%phi2(3)*dr(3))*dr(1) &
     +0.5*(lrf(ic)%phi2(4)*dr(1)+lrf(ic)%phi2(5)*dr(2) &
     +lrf(ic)%phi2(6)*dr(3))*dr(2) &
     +0.5*(lrf(ic)%phi2(7)*dr(1)+lrf(ic)%phi2(8)*dr(2) &
     +lrf(ic)%phi2(9)*dr(3))*dr(3)

E%LRF = E%LRF + .5 * crg(i) * Vij

! --- Electric field
df(1)=lrf(ic)%phi1(1) &
     +lrf(ic)%phi2(1)*dr(1)+lrf(ic)%phi2(2)*dr(2)+lrf(ic)%phi2(3)*dr(3) &
     +0.5*(lrf(ic)%phi3(1 )*dr(1)+lrf(ic)%phi3(2 )*dr(2) &
     +lrf(ic)%phi3(3 )*dr(3))*dr(1) &
     +0.5*(lrf(ic)%phi3(4 )*dr(1)+lrf(ic)%phi3(5 )*dr(2) &
     +lrf(ic)%phi3(6 )*dr(3))*dr(2) &
     +0.5*(lrf(ic)%phi3(7 )*dr(1)+lrf(ic)%phi3(8 )*dr(2) &
     +lrf(ic)%phi3(9 )*dr(3))*dr(3)
df(2)=lrf(ic)%phi1(2) &
     +lrf(ic)%phi2(4)*dr(1)+lrf(ic)%phi2(5)*dr(2)+lrf(ic)%phi2(6)*dr(3) &
     +0.5*(lrf(ic)%phi3(10)*dr(1)+lrf(ic)%phi3(11)*dr(2) &
     +lrf(ic)%phi3(12)*dr(3))*dr(1) &
     +0.5*(lrf(ic)%phi3(13)*dr(1)+lrf(ic)%phi3(14)*dr(2) &
     +lrf(ic)%phi3(15)*dr(3))*dr(2) &
     +0.5*(lrf(ic)%phi3(16)*dr(1)+lrf(ic)%phi3(17)*dr(2) &
     +lrf(ic)%phi3(18)*dr(3))*dr(3)
df(3)=lrf(ic)%phi1(3) &
     +lrf(ic)%phi2(7)*dr(1)+lrf(ic)%phi2(8)*dr(2)+lrf(ic)%phi2(9)*dr(3) &
     +0.5*(lrf(ic)%phi3(19)*dr(1)+lrf(ic)%phi3(20)*dr(2) &
     +lrf(ic)%phi3(21)*dr(3))*dr(1) &
     +0.5*(lrf(ic)%phi3(22)*dr(1)+lrf(ic)%phi3(23)*dr(2) &
     +lrf(ic)%phi3(24)*dr(3))*dr(2) &
     +0.5*(lrf(ic)%phi3(25)*dr(1)+lrf(ic)%phi3(26)*dr(2) &
     +lrf(ic)%phi3(27)*dr(3))*dr(3) 

  ! update d
d(i3+1)=d(i3+1)-crg(i)*df(1)
d(i3+2)=d(i3+2)-crg(i)*df(2)
d(i3+3)=d(i3+3)-crg(i)*df(3)
end if
end do
end subroutine lrf_taylor


!-----------------------------------------------------------------------

subroutine make_pair_lists
#if defined (PROFILING)
real(8)                                         :: start_loop_time
start_loop_time = rtime()
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

#if defined (PROFILING)
profile(1)%time = profile(1)%time + rtime() - start_loop_time
#endif

end subroutine make_pair_lists

!-----------------------------------------------------------------------

subroutine maxwell
! *** local variables
integer						:: i,j,k
real(8)						:: sd,vg,kT

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
real(8)						:: Tscale_solute,Tscale_solvent,Ekinmax

!locals
integer						::	i, i3
real(8)						::	Ekin

Temp = 0.
Temp_solute = 0.
Tfree_solute = 0.
Texcl_solute = 0.

!get kinetic energies for solute atoms
do i=1,nat_solute
        i3=i*3-3
        Ekin = 0.5*iaclib(iac(i))%mass*(v(i3+1)**2+v(i3+2)**2+v(i3+3)**2)
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
                write (*,180) i,2.*Ekin/Boltz/3.
        end if
end do

Tfree_solvent = 0.
Temp_solvent = 0.
Texcl_solvent = 0.
Ekin = 0


!get kinetic energies for solvent atoms
do i=nat_solute+1,natom
        i3=i*3-3
        Ekin = 0.5*iaclib(iac(i))%mass*(v(i3+1)**2+v(i3+2)**2+v(i3+3)**2)
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
                write (*,180) i,2.*Ekin/Boltz/3.
        end if
end do

Tfree = Tfree_solvent + Tfree_solute
Temp = Temp_solute + Temp_solvent

E%kinetic = Temp

Temp  = 2.0*Temp/Boltz/real(Ndegf)
Tfree = 2.0*Tfree/Boltz/real(Ndegfree)

if (detail_temps) then
	Temp_solute  = 2.0*Temp_solute /Boltz/real(Ndegf_solute)
	Tfree_solute = 2.0*Tfree_solute/Boltz/real(Ndegfree_solute)
	if ( Ndegf_solute .ne. Ndegfree_solute) Texcl_solute = 2.0*Texcl_solute/Boltz/real(Ndegf_solute - Ndegfree_solute)

	Temp_solvent  = 2.0*Temp_solvent /Boltz/real(Ndegf_solvent)
	Tfree_solvent = 2.0*Tfree_solvent/Boltz/real(Ndegfree_solvent)
	if ( Ndegf_solvent .ne. Ndegfree_solvent) Texcl_solvent = 2.0*Texcl_solvent/Boltz/real(Ndegf_solvent - Ndegfree_solvent)
end if

if (thermostat == BERENDSEN) then
	if (separate_scaling) then
		if ( Tfree_solvent .ne. 0 ) Tscale_solvent = Temp0/Tfree_solvent - 1.0
		Tscale_solvent = sqrt ( 1 + dt/tau_T * Tscale_solvent )
		if ( Tfree_solute .ne. 0 ) Tscale_solute = Temp0/Tfree_solute - 1.0
		Tscale_solute = sqrt ( 1 + dt/tau_T * Tscale_solute )
	else
		if ( Tfree .ne. 0 ) Tscale_solvent = Temp0/Tfree - 1.0
		Tscale_solvent = sqrt ( 1 + dt/tau_T * Tscale_solvent )
		Tscale_solute = Tscale_solvent
	end if
end if

180 format ('>>> WARNING: hot atom, i =',i10,' Temp(i)=',f10.2)

end subroutine temperature

!----------------------------------------------------------------------
! Subroutine for the Nosé-Hoover chain propagation
subroutine nh_prop

real(8)				:: dt4, dt8, expf, s
integer				:: i, i3

!initializing all the constants to be used
dt4=0.5*dt2
dt8=0.5*dt4

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
real(8)                         :: Tlast
real(8)                         ::Ekinmax
real(8)				::Tscale_solute,Tscale_solvent
!new for temperature control
real(8)				::Tscale,dv_mod,dv_friction,dv_mod2
integer				:: n_max = -1
real(8),allocatable	::randva(:),randfa(:,:),dv(:)
!Random variable temperature control array
real(8)				:: time0, time1, time_per_step, startloop
integer(4)				:: time_completion

#if defined(PROFILING)
real(8)                                         :: start_loop_time1, start_loop_time2

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

all_node_times(:) = 0.0
node_times(:) = 0.0

#endif
#endif


!Define number of coord to send/recieve
nat3=natom*3

! calculate maximum temperature
!**MN-> Only master calc. temp for now.
if (nodeid .eq. 0) then
Ekinmax = 1000.0*Ndegf*Boltz*Temp0/2.0/real(natom)

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
	allocate(randva(natom))
	call check_alloc('Random variable temperature control array')
	allocate(randfa(natom,3))
	call check_alloc('Random variable temperature control array 2')
        allocate(dv(nat3))
        call check_alloc('Delta v array for temperature control')
! if velocity verlet is used, take half steps
if ( integrator == VELVERLET ) then
        dv_mod = dt2
else
        dv_mod = dt
end if
if ( thermostat == BERENDSEN ) then
	dv_friction = 1
	n_max = nat_solute
	randfa(:,:) = 0
else if ( thermostat == NOSEHOOVER )then
	dv_friction = 1
	randfa(:,:) = 0
	n_max = natom
else if ( thermostat == LANGEVIN ) then
	dv_friction = (1-friction*dv_mod)
        gkT = 2*friction*Boltz*Temp0/dv_mod
	randva(:)= sqrt (gkT/winv(:))
	n_max = natom
else
	write(*,*) 'No such thermostat'
	stop
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
					 &  nbww_true_pair, nbqp_pair, 3*nqat*nbqw_pair

#if defined(USE_MPI)
				!reduce totnxx, i.e. collect # pairs found by slave nodes
				nbxx(1)=nbpp_pair
				nbxx(2)=nbpw_pair
				nbxx(3)=nbww_true_pair
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

! if Berendsen thermostat was chosen, applied here
!because Tscale changes all the time
if( thermostat == BERENDSEN ) then
		Tscale = Tscale_solute
else
		Tscale = 1
end if
if (thermostat == LANGEVIN ) then
        do i=1,natom
                do j=1,3
	        call gauss(zero,randva(i),randfa(i,j),iseed)
                end do 
        end do
end if
if (thermostat == NOSEHOOVER ) then
	call nh_prop
end if
!Now general equation for all thermostats
!to reduce code duplication
!Paul Bauer October 2014
	dv_mod2 = dv_mod * Tscale
        do i=1,n_max
                i3=i*3-3
		dv(i3+1) = (d(i3+1)-randfa(i,1))*winv(i)*dv_mod2
		dv(i3+2) = (d(i3+2)-randfa(i,2))*winv(i)*dv_mod2
		dv(i3+3) = (d(i3+3)-randfa(i,3))*winv(i)*dv_mod2

                v(i3+1)  = (v(i3+1)*Tscale*dv_friction) -dv(i3+1)
                xx(i3+1) = x(i3+1)
                x(i3+1)  = x(i3+1) + v(i3+1)*dt

                v(i3+2)  = (v(i3+2)*Tscale*dv_friction) -dv(i3+2)
                xx(i3+2) = x(i3+2)
                x(i3+2)  = x(i3+2) + v(i3+2)*dt

                v(i3+3)  = (v(i3+3)*Tscale*dv_friction) -dv(i3+3)
                xx(i3+3) = x(i3+3)
                x(i3+3)  = x(i3+3) + v(i3+3)*dt

		if (integrator == VELVERLET) then
			v(i3+1)= (v(i3+1)*Tscale*dv_friction) -dv(i3+1)
			v(i3+2)= (v(i3+2)*Tscale*dv_friction) -dv(i3+2)
			v(i3+3)= (v(i3+3)*Tscale*dv_friction) -dv(i3+3)
		end if
        end do
if( thermostat == BERENDSEN ) then
	Tscale = Tscale_solvent
	dv_mod2 = dv_mod * Tscale
        do i=n_max+1,natom
                i3=i*3-3
		dv(i3+1) = d(i3+1)*winv(i)*dv_mod2
		dv(i3+2) = d(i3+2)*winv(i)*dv_mod2
		dv(i3+3) = d(i3+3)*winv(i)*dv_mod2

                v(i3+1)= (v(i3+1)*Tscale) -dv(i3+1)
                xx(i3+1) = x(i3+1)
                x(i3+1) = x(i3+1) + v(i3+1)*dt

                v(i3+2)= (v(i3+2)*Tscale) -dv(i3+2)
                xx(i3+2) = x(i3+2)
                x(i3+2) = x(i3+2) + v(i3+2)*dt

                v(i3+3)= (v(i3+3)*Tscale) -dv(i3+3)
                xx(i3+3) = x(i3+3)
                x(i3+3) = x(i3+3) + v(i3+3)*dt

		if (integrator == VELVERLET) then
			v(i3+1)= (v(i3+1)*Tscale) -dv(i3+1)
			v(i3+2)= (v(i3+2)*Tscale) -dv(i3+2)
			v(i3+3)= (v(i3+3)*Tscale) -dv(i3+3)
		end if
        end do
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
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast x')
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
                        call put_ene(11, EQ, OFFD)
		! for new excluded groups, write more files for each group
			if(use_excluded_groups) then
				do iii=1,ngroups_gc
				call put_ene(ST_gc(iii)%fileunit,EQ_gc(iii,:),OFFD)
				end do
			end if
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
        deallocate(randfa)
        deallocate(randva)
        deallocate(dv)

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
real(8)						:: rcut2,r2
integer						:: LJ_code

real(8)						:: dx, dy, dz
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
          if((crg(i) * crg(j) == 0.) &
                .and. &
                (iaclib(iac(i))%avdw(LJ_code)*iaclib(iac(j))%avdw(LJ_code) == 0.) &
                .and. &
                (iaclib(iac(i))%bvdw(LJ_code)*iaclib(iac(j))%bvdw(LJ_code) == 0.)) &
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
integer						:: i,j,ig,jg,ia,ja,i3,j3,nl,inside
real(8)						:: rcut2,r2


! for spherical boundary 
!	This routine makes a list of non-bonded solute-solute atom pairs 
!	excluding any Q-atoms.

#if defined (PROFILING)
real(8)                                         :: start_loop_time
start_loop_time = rtime()
#endif

nbpp_pair = 0
rcut2 = Rcpp*Rcpp


igloop:  do ig = calculation_assignment%pp%start, calculation_assignment%pp%end

ia = cgp(ig)%iswitch
if ( excl(ia) ) cycle igloop

jgloop:     do jg = 1, ncgp_solute

  ! count each charge group pair once only
  if ( ((ig .gt. jg) .and. (mod(ig+jg,2) .eq. 0)) .or. &
           ((ig .lt. jg) .and. (mod(ig+jg,2) .eq. 1)) ) &
           cycle jgloop

ja = cgp(jg)%iswitch
         if ( excl(ja) ) cycle jgloop

!	      --- outside cutoff ? ---
        inside = 0
        ia = cgp(ig)%first
        do while ((ia .le. cgp(ig)%last) .and. (inside .eq. 0))
   i = cgpatom(ia)
   i3 = 3*i-3

           ja = cgp(jg)%first
           do while ((ja .le. cgp(jg)%last) .and. (inside .eq. 0))
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

ialoop:    do ia = cgp(ig)%first, cgp(ig)%last
   i = cgpatom(ia)

   !	         --- q-atom ? ---
   if ( iqatom(i)/=0 ) cycle ialoop

jaloop:       do ja = cgp(jg)%first, cgp(jg)%last
      j = cgpatom(ja)

      !	            --- q-atom ? ---
      if ( iqatom(j)/=0 ) cycle jaloop

                        ! count once
      if ( ig .eq. jg .and. i .ge. j ) cycle jaloop

      !	            --- check bonded exclusions and 1-4 nbors ---

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

                  ! if out of space then make more space
      if (nbpp_pair .eq. calculation_assignment%pp%max) call reallocate_nonbondlist_pp

      nbpp_pair = nbpp_pair + 1
      nbpp(nbpp_pair)%i = i
      nbpp(nbpp_pair)%j = j 
      nbpp(nbpp_pair)%LJcod = ljcod(iac(i),iac(j))

      if ( abs(j-i) .le. max_nbr_range ) then
         if ( i .lt. j ) then
            if ( list14(j-i,i) ) nbpp(nbpp_pair)%LJcod = 3
         else
            if ( list14(i-j,j) ) nbpp(nbpp_pair)%LJcod = 3
         end if
      else
         do nl = 1, n14long
            if ( (list14long(1,nl) .eq. i .and. &
                 list14long(2,nl) .eq. j      ) .or. &
                 (list14long(1,nl) .eq. j .and. &
                 list14long(2,nl) .eq. i      ) ) then
                 nbpp(nbpp_pair)%LJcod = 3
                                end if
         end do
      end if

                  end do jaloop
           end do ialoop
        end do jgloop
end do igloop
#if defined (PROFILING) 
profile(3)%time = profile(3)%time + rtime() - start_loop_time
#endif 

end subroutine nbpplis2


!--------------------------------------------------------------------

subroutine nbpplis2_box
  ! local variables
  integer						:: i,j,ig,jg,ia,ja,i3,j3,nl,inside
  real(8)						:: rcut2,r2
  real(8)						:: dx, dy, dz

#if defined (PROFILING)
real(8)                                         :: start_loop_time
start_loop_time = rtime()
#endif
 
  ! for periodic boundary conditions
  !	This routine makes a list of non-bonded solute-solute atom pairs 
  !	excluding any Q-atoms.

  nbpp_pair = 0
  nbpp_cgp_pair = 0
  rcut2 = Rcpp*Rcpp

igloop:  do ig = calculation_assignment%pp%start, calculation_assignment%pp%end
jgloop:     do jg = 1, ncgp_solute

	  ! count each charge group pair once only
	  if ( ((ig .gt. jg) .and. (mod(ig+jg,2) .eq. 0)) .or. &
		   ((ig .lt. jg) .and. (mod(ig+jg,2) .eq. 1)) ) &
		   cycle jgloop
		
        !	      --- outside cutoff ? ---
		inside = 0
		ia = cgp(ig)%first
		do while ((ia .le. cgp(ig)%last) .and. (inside .eq. 0))
           i = cgpatom(ia)
           i3 = 3*i-3

		   ja = cgp(jg)%first
		   do while ((ja .le. cgp(jg)%last) .and. (inside .eq. 0))
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

				if (nbpp_cgp_pair .eq. size(nbpp_cgp, 1) )  call reallocate_nbpp_cgp

				nbpp_cgp_pair = nbpp_cgp_pair + 1
				nbpp_cgp(nbpp_cgp_pair)%i = i
				nbpp_cgp(nbpp_cgp_pair)%j = j

			  end if

			  ja = ja + 1
           end do

		   ia = ia + 1
        end do
		if (inside .eq. 0) cycle jgloop

ialoop:    do ia = cgp(ig)%first, cgp(ig)%last
           i = cgpatom(ia)

           !	         --- q-atom ? ---
           if ( iqatom(i)/=0 ) cycle ialoop

jaloop:       do ja = cgp(jg)%first, cgp(jg)%last
              j = cgpatom(ja)

              !	            --- q-atom ? ---
              if ( iqatom(j)/=0 ) cycle jaloop

				! count once
              if ( ig .eq. jg .and. i .ge. j ) cycle jaloop

              !	            --- check bonded exclusions and 1-4 nbors ---

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

			  ! if out of space then make more space
              if (nbpp_pair .eq. calculation_assignment%pp%max) call reallocate_nonbondlist_pp

              nbpp_pair = nbpp_pair + 1
              nbpp(nbpp_pair)%i = i
              nbpp(nbpp_pair)%j = j 
              nbpp(nbpp_pair)%LJcod = ljcod(iac(i),iac(j))
			  nbpp(nbpp_pair)%cgp_pair = nbpp_cgp_pair

              if ( abs(j-i) .le. max_nbr_range ) then
                 if ( i .lt. j ) then
                    if ( list14(j-i,i) ) nbpp(nbpp_pair)%LJcod = 3
                 else
                    if ( list14(i-j,j) ) nbpp(nbpp_pair)%LJcod = 3
                 end if
              else
                 do nl = 1, n14long
                    if ( (list14long(1,nl) .eq. i .and. &
                         list14long(2,nl) .eq. j      ) .or. &
                         (list14long(1,nl) .eq. j .and. &
                         list14long(2,nl) .eq. i      ) ) then
                         nbpp(nbpp_pair)%LJcod = 3
					end if
                 end do
              end if

			  end do jaloop
		   end do ialoop
		end do jgloop
	end do igloop
#if defined (PROFILING) 
profile(3)%time = profile(3)%time + rtime() - start_loop_time
#endif 

end subroutine nbpplis2_box
!-----------------------------------------------------------------------------
subroutine nbpplis2_box_lrf
  ! local variables
  integer						:: i,j,ig,jg,ia,ja,i3,j3,nl,inside
  real(8)						:: rcut2,r2
  real(8)						:: dx, dy, dz
 
  real(8)						::RcLRF2,field0, field1, field2
  real(8)						::dr(3)
  real(8)						::boxshiftx, boxshifty, boxshiftz	
  integer						::inside_LRF, is3
  ! for periodic boundary conditions
  !	This routine makes a list of non-bonded solute-solute atom pairs 
  !	excluding any Q-atoms.
#if defined (PROFILING)
real(8)                                         :: start_loop_time
start_loop_time = rtime()
#endif

  nbpp_pair = 0
  nbpp_cgp_pair = 0
  rcut2 = Rcpp*Rcpp
  RcLRF2 = RcLRF*RcLRF

igloop:  do ig = calculation_assignment%pp%start, calculation_assignment%pp%end
jgloop:     do jg = 1, ncgp_solute

	  ! count each charge group pair once only
	  if ( ((ig .gt. jg) .and. (mod(ig+jg,2) .eq. 0)) .or. &
		   ((ig .lt. jg) .and. (mod(ig+jg,2) .eq. 1)) ) &
		   cycle jgloop
		
        !	      --- outside cutoff ? ---
		inside = 0
		inside_LRF = 0
		ia = cgp(ig)%first
		do while ((ia .le. cgp(ig)%last) .and. (inside .eq. 0))
           i = cgpatom(ia)
           i3 = 3*i-3

		   ja = cgp(jg)%first
		   do while ((ja .le. cgp(jg)%last) .and. (inside .eq. 0))
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

				if (nbpp_cgp_pair .eq. size(nbpp_cgp, 1) )  call reallocate_nbpp_cgp

				nbpp_cgp_pair = nbpp_cgp_pair + 1
				nbpp_cgp(nbpp_cgp_pair)%i = i
				nbpp_cgp(nbpp_cgp_pair)%j = j
              elseif (r2 <= RcLRF2) then
			    inside_LRF = 1
              end if
			  ja = ja + 1
           end do
		   ia = ia + 1
        end do
		
		if (inside .eq. 1) then
ialoop:    do ia = cgp(ig)%first, cgp(ig)%last
           i = cgpatom(ia)

           !	         --- q-atom ? ---
           if ( iqatom(i)/=0 ) cycle ialoop

jaloop:       do ja = cgp(jg)%first, cgp(jg)%last
              j = cgpatom(ja)

              !	            --- q-atom ? ---
              if ( iqatom(j)/=0 ) cycle jaloop

				! count once
              if ( ig .eq. jg .and. i .ge. j ) cycle jaloop

              !	            --- check bonded exclusions and 1-4 nbors ---

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

			  ! if out of space then make more space
              if (nbpp_pair .eq. calculation_assignment%pp%max) call reallocate_nonbondlist_pp

              nbpp_pair = nbpp_pair + 1
              nbpp(nbpp_pair)%i = i
              nbpp(nbpp_pair)%j = j 
              nbpp(nbpp_pair)%LJcod = ljcod(iac(i),iac(j))
			  nbpp(nbpp_pair)%cgp_pair = nbpp_cgp_pair

              if ( abs(j-i) .le. max_nbr_range ) then
                 if ( i .lt. j ) then
                    if ( list14(j-i,i) ) nbpp(nbpp_pair)%LJcod = 3
                 else
                    if ( list14(i-j,j) ) nbpp(nbpp_pair)%LJcod = 3
                 end if
              else
                 do nl = 1, n14long
                    if ( (list14long(1,nl) .eq. i .and. &
                         list14long(2,nl) .eq. j      ) .or. &
                         (list14long(1,nl) .eq. j .and. &
                         list14long(2,nl) .eq. i      ) ) then
                         nbpp(nbpp_pair)%LJcod = 3
					end if
                 end do
              end if

			  end do jaloop
		   end do ialoop
    	elseif((inside_LRF ==1) .and. (inside == 0)) then
        ! outside pp-cutoff but inside LRF cut-off use LRF
		
		!ig : jg calculation
		boxshiftx = x(i3+1) - lrf(jg)%cgp_cent(1)
		boxshifty = x(i3+2) - lrf(jg)%cgp_cent(2)
		boxshiftz = x(i3+3) - lrf(jg)%cgp_cent(3)

		boxshiftx = boxlength(1)*nint( boxshiftx*inv_boxl(1) )
		boxshifty = boxlength(2)*nint( boxshifty*inv_boxl(2) )
		boxshiftz = boxlength(3)*nint( boxshiftz*inv_boxl(3) )

        do ia = cgp(ig)%first, cgp(ig)%last

          ! skip if q-atom
          i = cgpatom(ia)
          if ( iqatom(i)/=0 ) cycle

          is3 = i*3-3

          dr(1) = x(is3+1) - lrf(jg)%cgp_cent(1) - boxshiftx
          dr(2) = x(is3+2) - lrf(jg)%cgp_cent(2) - boxshifty
          dr(3) = x(is3+3) - lrf(jg)%cgp_cent(3) - boxshiftz
       	  r2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)

          field0=crg(i)/(r2*sqrt(r2))
          lrf(jg)%phi0=lrf(jg)%phi0+field0*r2
          lrf(jg)%phi1(1)=lrf(jg)%phi1(1)-field0*dr(1)
          lrf(jg)%phi1(2)=lrf(jg)%phi1(2)-field0*dr(2)
          lrf(jg)%phi1(3)=lrf(jg)%phi1(3)-field0*dr(3)
          field1=3.*field0/r2
          lrf(jg)%phi2(1)=lrf(jg)%phi2(1)+field1*dr(1)*dr(1)-field0
          lrf(jg)%phi2(2)=lrf(jg)%phi2(2)+field1*dr(1)*dr(2)
          lrf(jg)%phi2(3)=lrf(jg)%phi2(3)+field1*dr(1)*dr(3)
          lrf(jg)%phi2(4)=lrf(jg)%phi2(4)+field1*dr(2)*dr(1)
          lrf(jg)%phi2(5)=lrf(jg)%phi2(5)+field1*dr(2)*dr(2)-field0
          lrf(jg)%phi2(6)=lrf(jg)%phi2(6)+field1*dr(2)*dr(3)
          lrf(jg)%phi2(7)=lrf(jg)%phi2(7)+field1*dr(3)*dr(1)
          lrf(jg)%phi2(8)=lrf(jg)%phi2(8)+field1*dr(3)*dr(2)
          lrf(jg)%phi2(9)=lrf(jg)%phi2(9)+field1*dr(3)*dr(3)-field0
          field2=-field1/r2
          lrf(jg)%phi3(1 )=lrf(jg)%phi3(1 ) +field2*(5.*dr(1)*dr(1)*dr(1)-r2*3.*dr(1))
          lrf(jg)%phi3(2 )=lrf(jg)%phi3(2 ) +field2*(5.*dr(1)*dr(1)*dr(2)-r2*dr(2))
          lrf(jg)%phi3(3 )=lrf(jg)%phi3(3 ) +field2*(5.*dr(1)*dr(1)*dr(3)-r2*dr(3))
          lrf(jg)%phi3(4 )=lrf(jg)%phi3(4 ) +field2*(5.*dr(1)*dr(2)*dr(1)-r2*dr(2))
          lrf(jg)%phi3(5 )=lrf(jg)%phi3(5 ) +field2*(5.*dr(1)*dr(2)*dr(2)-r2*dr(1))
          lrf(jg)%phi3(6 )=lrf(jg)%phi3(6 ) +field2*(5.*dr(1)*dr(2)*dr(3))
          lrf(jg)%phi3(7 )=lrf(jg)%phi3(7 ) +field2*(5.*dr(1)*dr(3)*dr(1)-r2*dr(3))
          lrf(jg)%phi3(8 )=lrf(jg)%phi3(8 ) +field2*(5.*dr(1)*dr(3)*dr(2))
          lrf(jg)%phi3(9 )=lrf(jg)%phi3(9 ) +field2*(5.*dr(1)*dr(3)*dr(3)-r2*dr(1))
          lrf(jg)%phi3(10)=lrf(jg)%phi3(10) +field2*(5.*dr(2)*dr(1)*dr(1)-r2*dr(2))
          lrf(jg)%phi3(11)=lrf(jg)%phi3(11) +field2*(5.*dr(2)*dr(1)*dr(2)-r2*dr(1))
          lrf(jg)%phi3(12)=lrf(jg)%phi3(12) +field2*(5.*dr(2)*dr(1)*dr(3))
          lrf(jg)%phi3(13)=lrf(jg)%phi3(13) +field2*(5.*dr(2)*dr(2)*dr(1)-r2*dr(1))
          lrf(jg)%phi3(14)=lrf(jg)%phi3(14) +field2*(5.*dr(2)*dr(2)*dr(2)-r2*3.*dr(2))
          lrf(jg)%phi3(15)=lrf(jg)%phi3(15) +field2*(5.*dr(2)*dr(2)*dr(3)-r2*dr(3))
          lrf(jg)%phi3(16)=lrf(jg)%phi3(16) +field2*(5.*dr(2)*dr(3)*dr(1))
          lrf(jg)%phi3(17)=lrf(jg)%phi3(17) +field2*(5.*dr(2)*dr(3)*dr(2)-r2*dr(3))
          lrf(jg)%phi3(18)=lrf(jg)%phi3(18) +field2*(5.*dr(2)*dr(3)*dr(3)-r2*dr(2))
          lrf(jg)%phi3(19)=lrf(jg)%phi3(19) +field2*(5.*dr(3)*dr(1)*dr(1)-r2*dr(3))
          lrf(jg)%phi3(20)=lrf(jg)%phi3(20) +field2*(5.*dr(3)*dr(1)*dr(2))
          lrf(jg)%phi3(21)=lrf(jg)%phi3(21) +field2*(5.*dr(3)*dr(1)*dr(3)-r2*dr(1))
          lrf(jg)%phi3(22)=lrf(jg)%phi3(22) +field2*(5.*dr(3)*dr(2)*dr(1))
          lrf(jg)%phi3(23)=lrf(jg)%phi3(23) +field2*(5.*dr(3)*dr(2)*dr(2)-r2*dr(3))
          lrf(jg)%phi3(24)=lrf(jg)%phi3(24) +field2*(5.*dr(3)*dr(2)*dr(3)-r2*dr(2))
          lrf(jg)%phi3(25)=lrf(jg)%phi3(25) +field2*(5.*dr(3)*dr(3)*dr(1)-r2*dr(1))
          lrf(jg)%phi3(26)=lrf(jg)%phi3(26) +field2*(5.*dr(3)*dr(3)*dr(2)-r2*dr(2))
          lrf(jg)%phi3(27)=lrf(jg)%phi3(27) +field2*(5.*dr(3)*dr(3)*dr(3)-r2*3.*dr(3))
        end do

		!jg : ig calculations
		boxshiftx = x(j3+1) - lrf(ig)%cgp_cent(1)
		boxshifty = x(j3+2) - lrf(ig)%cgp_cent(2)
		boxshiftz = x(j3+3) - lrf(ig)%cgp_cent(3)
	
		boxshiftx = boxlength(1)*nint( boxshiftx*inv_boxl(1) )
		boxshifty = boxlength(2)*nint( boxshifty*inv_boxl(2) )
		boxshiftz = boxlength(3)*nint( boxshiftz*inv_boxl(3) )

	    do ja = cgp(jg)%first, cgp(jg)%last

          ! skip if q-atom
          j = cgpatom(ja)
          if ( iqatom(j)/=0 ) cycle

          j3 = j*3-3

          dr(1) = x(j3+1) - lrf(ig)%cgp_cent(1) - boxshiftx
          dr(2) = x(j3+2) - lrf(ig)%cgp_cent(2) - boxshifty
          dr(3) = x(j3+3) - lrf(ig)%cgp_cent(3) - boxshiftz
          r2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)

          field0=crg(j)/(r2*sqrt(r2))
          lrf(ig)%phi0=lrf(ig)%phi0+field0*r2
          lrf(ig)%phi1(1)=lrf(ig)%phi1(1)-field0*dr(1)
          lrf(ig)%phi1(2)=lrf(ig)%phi1(2)-field0*dr(2)
          lrf(ig)%phi1(3)=lrf(ig)%phi1(3)-field0*dr(3)
          field1=3.*field0/r2
          lrf(ig)%phi2(1)=lrf(ig)%phi2(1)+field1*dr(1)*dr(1)-field0
          lrf(ig)%phi2(2)=lrf(ig)%phi2(2)+field1*dr(1)*dr(2)
          lrf(ig)%phi2(3)=lrf(ig)%phi2(3)+field1*dr(1)*dr(3)
          lrf(ig)%phi2(4)=lrf(ig)%phi2(4)+field1*dr(2)*dr(1)
          lrf(ig)%phi2(5)=lrf(ig)%phi2(5)+field1*dr(2)*dr(2)-field0
          lrf(ig)%phi2(6)=lrf(ig)%phi2(6)+field1*dr(2)*dr(3)
          lrf(ig)%phi2(7)=lrf(ig)%phi2(7)+field1*dr(3)*dr(1)
          lrf(ig)%phi2(8)=lrf(ig)%phi2(8)+field1*dr(3)*dr(2)
          lrf(ig)%phi2(9)=lrf(ig)%phi2(9)+field1*dr(3)*dr(3)-field0
          field2=-field1/r2
          lrf(ig)%phi3(1 )=lrf(ig)%phi3(1 ) +field2*(5.*dr(1)*dr(1)*dr(1)-r2*3.*dr(1))
          lrf(ig)%phi3(2 )=lrf(ig)%phi3(2 ) +field2*(5.*dr(1)*dr(1)*dr(2)-r2*dr(2))
          lrf(ig)%phi3(3 )=lrf(ig)%phi3(3 ) +field2*(5.*dr(1)*dr(1)*dr(3)-r2*dr(3))
          lrf(ig)%phi3(4 )=lrf(ig)%phi3(4 ) +field2*(5.*dr(1)*dr(2)*dr(1)-r2*dr(2))
          lrf(ig)%phi3(5 )=lrf(ig)%phi3(5 ) +field2*(5.*dr(1)*dr(2)*dr(2)-r2*dr(1))
          lrf(ig)%phi3(6 )=lrf(ig)%phi3(6 ) +field2*(5.*dr(1)*dr(2)*dr(3))
          lrf(ig)%phi3(7 )=lrf(ig)%phi3(7 ) +field2*(5.*dr(1)*dr(3)*dr(1)-r2*dr(3))
          lrf(ig)%phi3(8 )=lrf(ig)%phi3(8 ) +field2*(5.*dr(1)*dr(3)*dr(2))
          lrf(ig)%phi3(9 )=lrf(ig)%phi3(9 ) +field2*(5.*dr(1)*dr(3)*dr(3)-r2*dr(1))
          lrf(ig)%phi3(10)=lrf(ig)%phi3(10) +field2*(5.*dr(2)*dr(1)*dr(1)-r2*dr(2))
          lrf(ig)%phi3(11)=lrf(ig)%phi3(11) +field2*(5.*dr(2)*dr(1)*dr(2)-r2*dr(1))
          lrf(ig)%phi3(12)=lrf(ig)%phi3(12) +field2*(5.*dr(2)*dr(1)*dr(3))
          lrf(ig)%phi3(13)=lrf(ig)%phi3(13) +field2*(5.*dr(2)*dr(2)*dr(1)-r2*dr(1))
          lrf(ig)%phi3(14)=lrf(ig)%phi3(14) +field2*(5.*dr(2)*dr(2)*dr(2)-r2*3.*dr(2))
          lrf(ig)%phi3(15)=lrf(ig)%phi3(15) +field2*(5.*dr(2)*dr(2)*dr(3)-r2*dr(3))
          lrf(ig)%phi3(16)=lrf(ig)%phi3(16) +field2*(5.*dr(2)*dr(3)*dr(1))
          lrf(ig)%phi3(17)=lrf(ig)%phi3(17) +field2*(5.*dr(2)*dr(3)*dr(2)-r2*dr(3))
          lrf(ig)%phi3(18)=lrf(ig)%phi3(18) +field2*(5.*dr(2)*dr(3)*dr(3)-r2*dr(2))
          lrf(ig)%phi3(19)=lrf(ig)%phi3(19) +field2*(5.*dr(3)*dr(1)*dr(1)-r2*dr(3))
          lrf(ig)%phi3(20)=lrf(ig)%phi3(20) +field2*(5.*dr(3)*dr(1)*dr(2))
          lrf(ig)%phi3(21)=lrf(ig)%phi3(21) +field2*(5.*dr(3)*dr(1)*dr(3)-r2*dr(1))
          lrf(ig)%phi3(22)=lrf(ig)%phi3(22) +field2*(5.*dr(3)*dr(2)*dr(1))
          lrf(ig)%phi3(23)=lrf(ig)%phi3(23) +field2*(5.*dr(3)*dr(2)*dr(2)-r2*dr(3))
          lrf(ig)%phi3(24)=lrf(ig)%phi3(24) +field2*(5.*dr(3)*dr(2)*dr(3)-r2*dr(2))
          lrf(ig)%phi3(25)=lrf(ig)%phi3(25) +field2*(5.*dr(3)*dr(3)*dr(1)-r2*dr(1))
          lrf(ig)%phi3(26)=lrf(ig)%phi3(26) +field2*(5.*dr(3)*dr(3)*dr(2)-r2*dr(2))
          lrf(ig)%phi3(27)=lrf(ig)%phi3(27) +field2*(5.*dr(3)*dr(3)*dr(3)-r2*3.*dr(3))
        enddo

      end if ! outside cutoff

	end do jgloop
  end do igloop
#if defined (PROFILING) 
profile(3)%time = profile(3)%time + rtime() - start_loop_time
#endif 

end subroutine nbpplis2_box_lrf
!-----------------------------------------------------------------------
subroutine nbpplis2_lrf
! local variables
integer						:: i,j,ig,jg,ia,ja,i3,j3,nl,is 
logical						::	inside
real(8)						:: rcut2,r2,field0,field1,field2
real(8)						:: dr(3)
real(8)						::	RcLRF2

!	This routine makes a list of non-bonded solute-solute atom pairs 
!	excluding any Q-atoms.
#if defined (PROFILING)
real(8)                                         :: start_loop_time
start_loop_time = rtime()
#endif

nbpp_pair = 0
rcut2 = Rcpp*Rcpp
RcLRF2 = RcLRF*RcLRF


igloop: do ig = calculation_assignment%pp%start, calculation_assignment%pp%end

        ! skip if excluded group
        is = cgp(ig)%iswitch
        if ( excl(is) ) cycle igloop

jgloop:	do jg = 1, ncgp_solute

                ! count each charge group pair once only
                if( ((ig .gt. jg) .and. (mod(ig+jg,2) .eq. 0)) .or. &
                        ((ig .lt. jg) .and. (mod(ig+jg,2) .eq. 1)) ) &
                        cycle jgloop

                !	      --- excluded group ? ---
                ja = cgp(jg)%iswitch
                if ( excl(ja) ) cycle jgloop

                !	      --- outside cutoff ? ---
                inside = .false.
pairloop:	do ia=cgp(ig)%first, cgp(ig)%last
                        i = cgpatom(ia)
                        i3 = 3*i-3

                        do ja=cgp(jg)%first, cgp(jg)%last
                                j = cgpatom(ja)
                                j3 = 3*j-3

                                r2 = ( x(i3+1) -x(j3+1) )**2 &
                                   +( x(i3+2) -x(j3+2) )**2 &
                                   +( x(i3+3) -x(j3+3) )**2

                                if ( r2 <= rcut2 ) then
                                        ! one atom pair is within cutoff: set inside
                                        inside = .true.
                                        exit pairloop
                                end if
                        end do
            end do pairloop

                !	      --- inside cutoff ? ---
                if (inside) then
ialoop:			do ia = cgp(ig)%first, cgp(ig)%last
                                i = cgpatom(ia)

                                !	             --- q-atom ? ---
                                if ( iqatom(i)/=0 ) cycle ialoop
jaloop:				do ja = cgp(jg)%first, cgp(jg)%last
                                        j = cgpatom(ja)
                                        !	               --- q-atom ? ---
                                        if ( iqatom(j)/=0 ) cycle jaloop
                                        ! count once
                                        if ( ig .eq. jg .and. i .ge. j ) cycle jaloop
                                        !	               --- check bonded exclusions and 1-4 nbors ---
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

                                        ! if out of space then make more space
                                        if (nbpp_pair == calculation_assignment%pp%max) call reallocate_nonbondlist_pp

                                        nbpp_pair = nbpp_pair + 1
                                        nbpp(nbpp_pair)%i = i
                                        nbpp(nbpp_pair)%j = j 
                                        nbpp(nbpp_pair)%LJcod = ljcod(iac(i),iac(j))

                                        if ( abs(j-i) .le. max_nbr_range ) then
                                                if ( i .lt. j ) then
                                                        if ( list14(j-i,i) ) nbpp(nbpp_pair)%LJcod = 3
                                                else
                                                        if ( list14(i-j,j) ) nbpp(nbpp_pair)%LJcod = 3
                                                end if
                                        else
                                                do nl = 1, n14long
                                                        if ( (list14long(1,nl) .eq. i .and. &
                                                                list14long(2,nl) .eq. j      ) .or. &
                                                                (list14long(1,nl) .eq. j .and. &
                                                                list14long(2,nl) .eq. i      ) ) &
                                                                nbpp(nbpp_pair)%LJcod = 3
                                                end do
                                        end if




                        end do jaloop
                        end do ialoop
                elseif(r2 <= RcLRF2) then
                        ! outside pp-cutoff but inside LRF cut-off: use LRF

ialoop2:		do ia = cgp(ig)%first, cgp(ig)%last

                                !	             --- q-atom ? ---
                                i = cgpatom(ia)
                                if ( iqatom(i)/=0 ) cycle ialoop2

                                i3 = i*3-3

                                dr(1) = x(i3+1) - lrf(jg)%cgp_cent(1)
                                dr(2) = x(i3+2) - lrf(jg)%cgp_cent(2)
                                dr(3) = x(i3+3) - lrf(jg)%cgp_cent(3)
                                r2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)

                                field0=crg(i)/(r2*sqrt(r2))
                                lrf(jg)%phi0=lrf(jg)%phi0+field0*r2
                                lrf(jg)%phi1(1)=lrf(jg)%phi1(1)-field0*dr(1)
                                lrf(jg)%phi1(2)=lrf(jg)%phi1(2)-field0*dr(2)
                                lrf(jg)%phi1(3)=lrf(jg)%phi1(3)-field0*dr(3)
                                field1=3.*field0/r2
                                lrf(jg)%phi2(1)=lrf(jg)%phi2(1)+field1*dr(1)*dr(1)-field0
                                lrf(jg)%phi2(2)=lrf(jg)%phi2(2)+field1*dr(1)*dr(2)
                                lrf(jg)%phi2(3)=lrf(jg)%phi2(3)+field1*dr(1)*dr(3)
                                lrf(jg)%phi2(4)=lrf(jg)%phi2(4)+field1*dr(2)*dr(1)
                                lrf(jg)%phi2(5)=lrf(jg)%phi2(5)+field1*dr(2)*dr(2)-field0
                                lrf(jg)%phi2(6)=lrf(jg)%phi2(6)+field1*dr(2)*dr(3)
                                lrf(jg)%phi2(7)=lrf(jg)%phi2(7)+field1*dr(3)*dr(1)
                                lrf(jg)%phi2(8)=lrf(jg)%phi2(8)+field1*dr(3)*dr(2)
                                lrf(jg)%phi2(9)=lrf(jg)%phi2(9)+field1*dr(3)*dr(3)-field0
                                field2=-field1/r2
                                lrf(jg)%phi3(1 )=lrf(jg)%phi3(1 ) +field2*(5.*dr(1)*dr(1)*dr(1)-r2*3.*dr(1))
                                lrf(jg)%phi3(2 )=lrf(jg)%phi3(2 ) +field2*(5.*dr(1)*dr(1)*dr(2)-r2*dr(2))
                                lrf(jg)%phi3(3 )=lrf(jg)%phi3(3 ) +field2*(5.*dr(1)*dr(1)*dr(3)-r2*dr(3))
                                lrf(jg)%phi3(4 )=lrf(jg)%phi3(4 ) +field2*(5.*dr(1)*dr(2)*dr(1)-r2*dr(2))
                                lrf(jg)%phi3(5 )=lrf(jg)%phi3(5 ) +field2*(5.*dr(1)*dr(2)*dr(2)-r2*dr(1))
                                lrf(jg)%phi3(6 )=lrf(jg)%phi3(6 ) +field2*(5.*dr(1)*dr(2)*dr(3))
                                lrf(jg)%phi3(7 )=lrf(jg)%phi3(7 ) +field2*(5.*dr(1)*dr(3)*dr(1)-r2*dr(3))
                                lrf(jg)%phi3(8 )=lrf(jg)%phi3(8 ) +field2*(5.*dr(1)*dr(3)*dr(2))
                                lrf(jg)%phi3(9 )=lrf(jg)%phi3(9 ) +field2*(5.*dr(1)*dr(3)*dr(3)-r2*dr(1))
                                lrf(jg)%phi3(10)=lrf(jg)%phi3(10) +field2*(5.*dr(2)*dr(1)*dr(1)-r2*dr(2))
                                lrf(jg)%phi3(11)=lrf(jg)%phi3(11) +field2*(5.*dr(2)*dr(1)*dr(2)-r2*dr(1))
                                lrf(jg)%phi3(12)=lrf(jg)%phi3(12) +field2*(5.*dr(2)*dr(1)*dr(3))
                                lrf(jg)%phi3(13)=lrf(jg)%phi3(13) +field2*(5.*dr(2)*dr(2)*dr(1)-r2*dr(1))
                                lrf(jg)%phi3(14)=lrf(jg)%phi3(14) +field2*(5.*dr(2)*dr(2)*dr(2)-r2*3.*dr(2))
                                lrf(jg)%phi3(15)=lrf(jg)%phi3(15) +field2*(5.*dr(2)*dr(2)*dr(3)-r2*dr(3))
                                lrf(jg)%phi3(16)=lrf(jg)%phi3(16) +field2*(5.*dr(2)*dr(3)*dr(1))
                                lrf(jg)%phi3(17)=lrf(jg)%phi3(17) +field2*(5.*dr(2)*dr(3)*dr(2)-r2*dr(3))
                                lrf(jg)%phi3(18)=lrf(jg)%phi3(18) +field2*(5.*dr(2)*dr(3)*dr(3)-r2*dr(2))
                                lrf(jg)%phi3(19)=lrf(jg)%phi3(19) +field2*(5.*dr(3)*dr(1)*dr(1)-r2*dr(3))
                                lrf(jg)%phi3(20)=lrf(jg)%phi3(20) +field2*(5.*dr(3)*dr(1)*dr(2))
                                lrf(jg)%phi3(21)=lrf(jg)%phi3(21) +field2*(5.*dr(3)*dr(1)*dr(3)-r2*dr(1))
                                lrf(jg)%phi3(22)=lrf(jg)%phi3(22) +field2*(5.*dr(3)*dr(2)*dr(1))
                                lrf(jg)%phi3(23)=lrf(jg)%phi3(23) +field2*(5.*dr(3)*dr(2)*dr(2)-r2*dr(3))
                                lrf(jg)%phi3(24)=lrf(jg)%phi3(24) +field2*(5.*dr(3)*dr(2)*dr(3)-r2*dr(2))
                                lrf(jg)%phi3(25)=lrf(jg)%phi3(25) +field2*(5.*dr(3)*dr(3)*dr(1)-r2*dr(1))
                                lrf(jg)%phi3(26)=lrf(jg)%phi3(26) +field2*(5.*dr(3)*dr(3)*dr(2)-r2*dr(2))
                                lrf(jg)%phi3(27)=lrf(jg)%phi3(27) +field2*(5.*dr(3)*dr(3)*dr(3)-r2*3.*dr(3))
                        end do ialoop2

jaloop2:		do ja = cgp(jg)%first, cgp(jg)%last

                                !	             --- q-atom ? ---
                                j = cgpatom(ja)
                                if ( iqatom(j)/=0 ) cycle jaloop2

                                j3 = j*3-3

                                dr(1) = x(j3+1) - lrf(ig)%cgp_cent(1)
                                dr(2) = x(j3+2) - lrf(ig)%cgp_cent(2)
                                dr(3) = x(j3+3) - lrf(ig)%cgp_cent(3)
                                r2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)

                                field0=crg(j)/(r2*sqrt(r2))
                                lrf(ig)%phi0=lrf(ig)%phi0+field0*r2
                                lrf(ig)%phi1(1)=lrf(ig)%phi1(1)-field0*dr(1)
                                lrf(ig)%phi1(2)=lrf(ig)%phi1(2)-field0*dr(2)
                                lrf(ig)%phi1(3)=lrf(ig)%phi1(3)-field0*dr(3)
                                field1=3.*field0/r2
                                lrf(ig)%phi2(1)=lrf(ig)%phi2(1)+field1*dr(1)*dr(1)-field0
                                lrf(ig)%phi2(2)=lrf(ig)%phi2(2)+field1*dr(1)*dr(2)


                        lrf(ig)%phi2(3)=lrf(ig)%phi2(3)+field1*dr(1)*dr(3)
                                lrf(ig)%phi2(4)=lrf(ig)%phi2(4)+field1*dr(2)*dr(1)
                                lrf(ig)%phi2(5)=lrf(ig)%phi2(5)+field1*dr(2)*dr(2)-field0
                                lrf(ig)%phi2(6)=lrf(ig)%phi2(6)+field1*dr(2)*dr(3)
                                lrf(ig)%phi2(7)=lrf(ig)%phi2(7)+field1*dr(3)*dr(1)
                                lrf(ig)%phi2(8)=lrf(ig)%phi2(8)+field1*dr(3)*dr(2)
                                lrf(ig)%phi2(9)=lrf(ig)%phi2(9)+field1*dr(3)*dr(3)-field0
                                field2=-field1/r2
                                lrf(ig)%phi3(1 )=lrf(ig)%phi3(1 ) +field2*(5.*dr(1)*dr(1)*dr(1)-r2*3.*dr(1))
                                lrf(ig)%phi3(2 )=lrf(ig)%phi3(2 ) +field2*(5.*dr(1)*dr(1)*dr(2)-r2*dr(2))
                                lrf(ig)%phi3(3 )=lrf(ig)%phi3(3 ) +field2*(5.*dr(1)*dr(1)*dr(3)-r2*dr(3))
                                lrf(ig)%phi3(4 )=lrf(ig)%phi3(4 ) +field2*(5.*dr(1)*dr(2)*dr(1)-r2*dr(2))
                                lrf(ig)%phi3(5 )=lrf(ig)%phi3(5 ) +field2*(5.*dr(1)*dr(2)*dr(2)-r2*dr(1))
                                lrf(ig)%phi3(6 )=lrf(ig)%phi3(6 ) +field2*(5.*dr(1)*dr(2)*dr(3))
                                lrf(ig)%phi3(7 )=lrf(ig)%phi3(7 ) +field2*(5.*dr(1)*dr(3)*dr(1)-r2*dr(3))
                                lrf(ig)%phi3(8 )=lrf(ig)%phi3(8 ) +field2*(5.*dr(1)*dr(3)*dr(2))
                                lrf(ig)%phi3(9 )=lrf(ig)%phi3(9 ) +field2*(5.*dr(1)*dr(3)*dr(3)-r2*dr(1))
                                lrf(ig)%phi3(10)=lrf(ig)%phi3(10) +field2*(5.*dr(2)*dr(1)*dr(1)-r2*dr(2))
                                lrf(ig)%phi3(11)=lrf(ig)%phi3(11) +field2*(5.*dr(2)*dr(1)*dr(2)-r2*dr(1))
                                lrf(ig)%phi3(12)=lrf(ig)%phi3(12) +field2*(5.*dr(2)*dr(1)*dr(3))
                                lrf(ig)%phi3(13)=lrf(ig)%phi3(13) +field2*(5.*dr(2)*dr(2)*dr(1)-r2*dr(1))
                                lrf(ig)%phi3(14)=lrf(ig)%phi3(14) +field2*(5.*dr(2)*dr(2)*dr(2)-r2*3.*dr(2))
                                lrf(ig)%phi3(15)=lrf(ig)%phi3(15) +field2*(5.*dr(2)*dr(2)*dr(3)-r2*dr(3))
                                lrf(ig)%phi3(16)=lrf(ig)%phi3(16) +field2*(5.*dr(2)*dr(3)*dr(1))
                                lrf(ig)%phi3(17)=lrf(ig)%phi3(17) +field2*(5.*dr(2)*dr(3)*dr(2)-r2*dr(3))
                                lrf(ig)%phi3(18)=lrf(ig)%phi3(18) +field2*(5.*dr(2)*dr(3)*dr(3)-r2*dr(2))
                                lrf(ig)%phi3(19)=lrf(ig)%phi3(19) +field2*(5.*dr(3)*dr(1)*dr(1)-r2*dr(3))
                                lrf(ig)%phi3(20)=lrf(ig)%phi3(20) +field2*(5.*dr(3)*dr(1)*dr(2))
                                lrf(ig)%phi3(21)=lrf(ig)%phi3(21) +field2*(5.*dr(3)*dr(1)*dr(3)-r2*dr(1))
                                lrf(ig)%phi3(22)=lrf(ig)%phi3(22) +field2*(5.*dr(3)*dr(2)*dr(1))
                                lrf(ig)%phi3(23)=lrf(ig)%phi3(23) +field2*(5.*dr(3)*dr(2)*dr(2)-r2*dr(3))
                                lrf(ig)%phi3(24)=lrf(ig)%phi3(24) +field2*(5.*dr(3)*dr(2)*dr(3)-r2*dr(2))
                                lrf(ig)%phi3(25)=lrf(ig)%phi3(25) +field2*(5.*dr(3)*dr(3)*dr(1)-r2*dr(1))
                                lrf(ig)%phi3(26)=lrf(ig)%phi3(26) +field2*(5.*dr(3)*dr(3)*dr(2)-r2*dr(2))
                                lrf(ig)%phi3(27)=lrf(ig)%phi3(27) +field2*(5.*dr(3)*dr(3)*dr(3)-r2*3.*dr(3))

                        end do jaloop2
                end if
        end do jgloop
end do igloop
#if defined (PROFILING) 
profile(3)%time = profile(3)%time + rtime() - start_loop_time
#endif 

end subroutine nbpplis2_lrf

!-----------------------------------------------------------------------

subroutine nbpplist
! local variables
integer						:: i,j,ig,jg,ia,ja,i3,j3,nl
real(8)						:: rcut2,r2
integer						:: LJ_code

! For use with spherical boundary   
!	This routine makes a list of non-bonded solute-solute atom pairs 
!	excluding any Q-atoms.

! uses the global variables:
!  Rcpp, ncgp, cgp, excl, x, cgpatom, iqatom, ljcod, crg, iaclib, max_nbr_range, listex
!  nexlong, listexlong, calculation_assignment%pp%max, alloc_status, list14, n14long, list14long

! reset #pairs

#if defined (PROFILING)
real(8)                                         :: start_loop_time
start_loop_time = rtime()
#endif

nbpp_pair = 0
rcut2 = Rcpp*Rcpp

igloop: do ig = calculation_assignment%pp%start, calculation_assignment%pp%end
! for every assigned charge group:

! skip if excluded group
ia = cgp(ig)%iswitch
if ( excl(ia) ) cycle igloop

i3 = 3*ia-3

jgloop: do jg = 1, ncgp_solute
  ! for every charge group:

  ! count each charge group pair once only
  if ( ((ig .gt. jg) .and. (mod(ig+jg,2) .eq. 0)) .or. &
           ((ig .lt. jg) .and. (mod(ig+jg,2) .eq. 1)) ) &
           cycle jgloop

  ! skip if excluded group
  ja = cgp(jg)%iswitch
  if ( excl(ja) ) cycle jgloop

  j3 = 3*ja-3
  r2 = ( x(i3+1) - x(j3+1) )**2 &
          +( x(i3+2) - x(j3+2) )**2 &
          +( x(i3+3) - x(j3+3) )**2

  ! skip if outside cutoff
  if ( r2 .gt. rcut2 ) cycle jgloop

ialoop: do ia = cgp(ig)%first, cgp(ig)%last
        ! for every atom in the charge group (of the outermost loop):
        i = cgpatom(ia)

        ! skip if q-atom
        if ( iqatom(i)/=0 ) cycle ialoop

jaloop: do ja = cgp(jg)%first, cgp(jg)%last
          ! for every atom in the charge group (innermost loop)
          j = cgpatom(ja)

          ! make sure each pair is only counted once
          if ( ig .eq. jg .and. i .ge. j ) cycle jaloop

          ! skip if q-atom
          if ( iqatom(j)/=0 ) cycle jaloop

          LJ_code = ljcod(iac(i),iac(j))

          ! skip if all interactions zero
          if((crg(i) * crg(j) == 0.) &
                .and. &
                (iaclib(iac(i))%avdw(LJ_code)*iaclib(iac(j))%avdw(LJ_code) == 0.) &
                .and. &
                (iaclib(iac(i))%bvdw(LJ_code)*iaclib(iac(j))%bvdw(LJ_code) == 0.)) &
                cycle jaloop

          ! skip bonded exclusions and 1-4 nbors
          if ( abs(j-i) .le. max_nbr_range ) then
                if ( i .lt. j ) then
                  if ( listex(j-i, i) ) cycle jaloop
                else
                  if ( listex(i-j, j) ) cycle jaloop
                end if
          else
                do nl = 1, nexlong
                  if ((listexlong(1, nl) .eq. i .and. listexlong(2, nl) .eq. j) .or. &
                        (listexlong(1, nl) .eq. j .and. listexlong(2, nl) .eq. i) ) &
                        cycle jaloop
                end do
          end if

          ! if out of space then make more space
          if (nbpp_pair .eq. calculation_assignment%pp%max) call reallocate_nonbondlist_pp

          ! all tests passed, add the pair
          nbpp_pair = nbpp_pair + 1
          nbpp(nbpp_pair)%i = i
          nbpp(nbpp_pair)%j = j 
          nbpp(nbpp_pair)%LJcod = LJ_code

          ! set LJcod of the pair to 3 if the atoms have bonded interactions
          if ( abs(j-i) .le. max_nbr_range ) then
                if ( i .lt. j ) then
                  if ( list14(j-i, i) ) nbpp(nbpp_pair)%LJcod = 3
                else
                  if ( list14(i-j, j) ) nbpp(nbpp_pair)%LJcod = 3
                end if
          else


        do nl = 1, n14long
                  if ((list14long(1, nl) .eq. i .and. list14long(2, nl) .eq. j) .or. &
                        (list14long(1, nl) .eq. j .and. list14long(2, nl) .eq. i)) &
                        nbpp(nbpp_pair)%LJcod = 3
                end do
          end if
        end do jaloop
  end do ialoop
end do jgloop
end do igloop
#if defined (PROFILING) 
profile(3)%time = profile(3)%time + rtime() - start_loop_time
#endif 

end subroutine nbpplist

!-----------------------------------------------------------------------

subroutine nbpplist_box
  ! local variables
  integer						:: i,j,ig,jg,ia,ja,i3,j3,nl,ig_sw, jg_sw
  real(8)						:: rcut2,r2
  integer						:: LJ_code
  real(8)						:: dx, dy, dz
  
  ! For use with periodic boundary conditions
  !	This routine makes a list of non-bonded solute-solute atom pairs 
  !	excluding any Q-atoms.

  ! uses the global variables:
  !  Rcpp, ncgp, cgp, excl, x, cgpatom, iqatom, ljcod, crg, iaclib, max_nbr_range, listex
  !  nexlong, listexlong, calculation_assignment%pp%max, alloc_status, list14, n14long, list14long

  ! reset #pairs

#if defined (PROFILING)
real(8)                                         :: start_loop_time
start_loop_time = rtime()
#endif

  nbpp_pair = 0 !atom pairs
  nbpp_cgp_pair = 0 !chargegroup pairs
  rcut2 = Rcpp*Rcpp

igloop: do ig = calculation_assignment%pp%start, calculation_assignment%pp%end
	! for every assigned charge group:
	ig_sw = cgp(ig)%iswitch !switching atom in charge group ig
	i3 = 3*ig_sw-3

jgloop: do jg = 1, ncgp_solute
	  ! for every charge group:

	  ! count each charge group pair once only
	  if ( ((ig .gt. jg) .and. (mod(ig+jg,2) .eq. 0)) .or. &
		   ((ig .lt. jg) .and. (mod(ig+jg,2) .eq. 1)) ) &
		   cycle jgloop

	  jg_sw = cgp(jg)%iswitch !switching atom in charge group jg
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
	  if(nbpp_cgp_pair .eq. size(nbpp_cgp, 1)) call reallocate_nbpp_cgp
	  !add the charge group pair
	  nbpp_cgp_pair = nbpp_cgp_pair + 1
	  nbpp_cgp(nbpp_cgp_pair)%i = ig_sw !the switching atoms of the charge groups in the pair
	  nbpp_cgp(nbpp_cgp_pair)%j = jg_sw

ialoop: do ia = cgp(ig)%first, cgp(ig)%last
		! for every atom in the charge group ig (of the outermost loop):
		i = cgpatom(ia)

		! skip if q-atom
		if ( iqatom(i)/=0 ) cycle ialoop

jaloop: do ja = cgp(jg)%first, cgp(jg)%last
		  ! for every atom in the charge group jg (innermost loop)
		  j = cgpatom(ja)

		  ! make sure each pair is only counted once
		  if ( ig .eq. jg .and. i .ge. j ) cycle jaloop

		  ! skip if q-atom
		  if ( iqatom(j)/=0 ) cycle jaloop

		  LJ_code = ljcod(iac(i),iac(j))

		  ! skip if all interactions zero
		  if((crg(i) * crg(j) == 0.) &
			.and. &
			(iaclib(iac(i))%avdw(LJ_code)*iaclib(iac(j))%avdw(LJ_code) == 0.) &
			.and. &
			(iaclib(iac(i))%bvdw(LJ_code)*iaclib(iac(j))%bvdw(LJ_code) == 0.)) &
			cycle jaloop

		  ! skip bonded exclusions and 1-4 nbors
		  if ( abs(j-i) .le. max_nbr_range ) then
			if ( i .lt. j ) then
			  if ( listex(j-i, i) ) cycle jaloop
			else
			  if ( listex(i-j, j) ) cycle jaloop
			end if
		  else
			do nl = 1, nexlong
			  if ((listexlong(1, nl) .eq. i .and. listexlong(2, nl) .eq. j) .or. &
				(listexlong(1, nl) .eq. j .and. listexlong(2, nl) .eq. i) ) &
				cycle jaloop
			end do
		  end if

		  ! if out of space then make more space
		  if (nbpp_pair .eq. calculation_assignment%pp%max) call reallocate_nonbondlist_pp

		  ! all tests passed, add the pair
		  nbpp_pair = nbpp_pair + 1
		  nbpp(nbpp_pair)%i = i
		  nbpp(nbpp_pair)%j = j 
		  nbpp(nbpp_pair)%LJcod = LJ_code
		  nbpp(nbpp_pair)%cgp_pair = nbpp_cgp_pair !which pair of charge groups the atom pair belongs to
		
		  ! set LJcod of the pair to 3 if the atoms have bonded interactions
		  if ( abs(j-i) .le. max_nbr_range ) then
			if ( i .lt. j ) then
			  if ( list14(j-i, i) ) nbpp(nbpp_pair)%LJcod = 3
			else
			  if ( list14(i-j, j) ) nbpp(nbpp_pair)%LJcod = 3
			end if
		  else

		do nl = 1, n14long
			  if ((list14long(1, nl) .eq. i .and. list14long(2, nl) .eq. j) .or. &
				(list14long(1, nl) .eq. j .and. list14long(2, nl) .eq. i)) &
				nbpp(nbpp_pair)%LJcod = 3
			end do
		  end if
		end do jaloop
	  end do ialoop
	end do jgloop
  end do igloop
#if defined (PROFILING) 
profile(3)%time = profile(3)%time + rtime() - start_loop_time
#endif 

end subroutine nbpplist_box

!--------------------------------------------------------------------------------------

subroutine nbpplist_lrf
! local variables
integer						:: i,j,ig,jg,ia,ja,i3,j3,nl,is,is3
real(8)						:: rcut2,r2,field0,field1,field2
real(8)						:: dr(3)
integer						:: LJ_code
real(8)						::	RcLRF2

#if defined (PROFILING)
real(8)						:: start_loop_time
start_loop_time = rtime()
#endif


!	This routine makes a list of non-bonded solute-solute atom pairs 
!	excluding any Q-atoms.

! uses the global variables:
!  nbpp_pair, Rcpp, RcLRF, cgp, excl, ncgp, x, cgpatom, iqatom, ljcod, crg, iaclib, max_nbr_range,
!  listex, nexlong, listexlong, nbpp, alloc_status, list14, n14long, list14long, lrf


nbpp_pair = 0
rcut2 = Rcpp*Rcpp
RcLRF2 = RcLRF*RcLRF
igloop: do ig = calculation_assignment%pp%start, calculation_assignment%pp%end
! for every assigned charge group:

! skip if excluded group
is = cgp(ig)%iswitch
if ( excl(is) ) cycle igloop

is3 = 3*is-3

jgloop: do jg = 1, ncgp_solute
  ! for every charge group:


  ! count each charge group pair once only
  if ( ((ig .gt. jg) .and. (mod(ig+jg,2) .eq. 0)) .or. &
           ((ig .lt. jg) .and. (mod(ig+jg,2) .eq. 1)) ) &
           cycle jgloop

  ! skip if excluded group
  ja = cgp(jg)%iswitch
  if ( excl(ja) ) cycle jgloop

  j3 = 3*ja-3
  r2 = ( x(is3+1) -x(j3+1) )**2 &
     +( x(is3+2) -x(j3+2) )**2 &
     +( x(is3+3) -x(j3+3) )**2

  ! inside cutoff?
  if ( r2 .le. rcut2 ) then
ialoop:	do ia = cgp(ig)%first, cgp(ig)%last

          ! skip if q-atom
          i = cgpatom(ia)
          if ( iqatom(i)/=0 ) cycle ialoop

jaloop:	  do ja = cgp(jg)%first, cgp(jg)%last
                j = cgpatom(ja)

                ! make sure each pair is only counted once
                if ( ig .eq. jg .and. i .ge. j ) cycle jaloop

                ! skip if q-atom
                if ( iqatom(j)/=0 ) cycle jaloop

                LJ_code = ljcod(iac(i),iac(j))

                ! skip if all interactions zero
                if ((crg(i) * crg(j) == 0.) &
                  .and. &
                  (iaclib(iac(i))%avdw(LJ_code)*iaclib(iac(j))%avdw(LJ_code) == 0.) &
                  .and. &
                  (iaclib(iac(i))%bvdw(LJ_code)*iaclib(iac(j))%bvdw(LJ_code) == 0.)) &
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
                          listexlong(2,nl) .eq. i ) ) &
                          cycle jaloop
                  end do
                end if

                ! if out of space then make more space
                if (nbpp_pair .eq. calculation_assignment%pp%max) call reallocate_nonbondlist_pp

                ! all tests passed, add the pair
                nbpp_pair = nbpp_pair + 1
                nbpp(nbpp_pair)%i = i
                nbpp(nbpp_pair)%j = j 
                nbpp(nbpp_pair)%LJcod = ljcod(iac(i),iac(j))
!TMP                nbpp_per_cgp(ig)=nbpp_per_cgp(ig)+1

                if ( abs(j-i) .le. max_nbr_range ) then
                  if ( i .lt. j ) then
                        if ( list14(j-i,i) ) nbpp(nbpp_pair)%LJcod = 3
                  else
                        if ( list14(i-j,j) ) nbpp(nbpp_pair)%LJcod = 3
                  end if
                else
                  do nl = 1, n14long
                        if ( (list14long(1,nl) .eq. i .and. &
                          list14long(2,nl) .eq. j      ) .or. &
                          (list14long(1,nl) .eq. j .and. &
                          list14long(2,nl) .eq. i      ) ) &
                          nbpp(nbpp_pair)%LJcod = 3
                  end do
                end if
          end do jaloop
        end do ialoop

  elseif(r2 <= RcLRF2) then
        ! outside pp-cutoff but inside LRF cut-off use LRF

        do ia = cgp(ig)%first, cgp(ig)%last

          ! skip if q-atom
          i = cgpatom(ia)
          if ( iqatom(i)/=0 ) cycle

          i3 = i*3-3

          dr(1) = x(i3+1) - lrf(jg)%cgp_cent(1)
          dr(2) = x(i3+2) - lrf(jg)%cgp_cent(2)
          dr(3) = x(i3+3) - lrf(jg)%cgp_cent(3)
          r2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)

          field0=crg(i)/(r2*sqrt(r2))
          lrf(jg)%phi0=lrf(jg)%phi0+field0*r2
          lrf(jg)%phi1(1)=lrf(jg)%phi1(1)-field0*dr(1)
          lrf(jg)%phi1(2)=lrf(jg)%phi1(2)-field0*dr(2)
          lrf(jg)%phi1(3)=lrf(jg)%phi1(3)-field0*dr(3)
          field1=3.*field0/r2
          lrf(jg)%phi2(1)=lrf(jg)%phi2(1)+field1*dr(1)*dr(1)-field0
          lrf(jg)%phi2(2)=lrf(jg)%phi2(2)+field1*dr(1)*dr(2)
          lrf(jg)%phi2(3)=lrf(jg)%phi2(3)+field1*dr(1)*dr(3)
          lrf(jg)%phi2(4)=lrf(jg)%phi2(4)+field1*dr(2)*dr(1)
          lrf(jg)%phi2(5)=lrf(jg)%phi2(5)+field1*dr(2)*dr(2)-field0
          lrf(jg)%phi2(6)=lrf(jg)%phi2(6)+field1*dr(2)*dr(3)
          lrf(jg)%phi2(7)=lrf(jg)%phi2(7)+field1*dr(3)*dr(1)
          lrf(jg)%phi2(8)=lrf(jg)%phi2(8)+field1*dr(3)*dr(2)
          lrf(jg)%phi2(9)=lrf(jg)%phi2(9)+field1*dr(3)*dr(3)-field0
          field2=-field1/r2
          lrf(jg)%phi3(1 )=lrf(jg)%phi3(1 ) +field2*(5.*dr(1)*dr(1)*dr(1)-r2*3.*dr(1))
          lrf(jg)%phi3(2 )=lrf(jg)%phi3(2 ) +field2*(5.*dr(1)*dr(1)*dr(2)-r2*dr(2))
          lrf(jg)%phi3(3 )=lrf(jg)%phi3(3 ) +field2*(5.*dr(1)*dr(1)*dr(3)-r2*dr(3))
          lrf(jg)%phi3(4 )=lrf(jg)%phi3(4 ) +field2*(5.*dr(1)*dr(2)*dr(1)-r2*dr(2))
          lrf(jg)%phi3(5 )=lrf(jg)%phi3(5 ) +field2*(5.*dr(1)*dr(2)*dr(2)-r2*dr(1))
          lrf(jg)%phi3(6 )=lrf(jg)%phi3(6 ) +field2*(5.*dr(1)*dr(2)*dr(3))
          lrf(jg)%phi3(7 )=lrf(jg)%phi3(7 ) +field2*(5.*dr(1)*dr(3)*dr(1)-r2*dr(3))
          lrf(jg)%phi3(8 )=lrf(jg)%phi3(8 ) +field2*(5.*dr(1)*dr(3)*dr(2))
          lrf(jg)%phi3(9 )=lrf(jg)%phi3(9 ) +field2*(5.*dr(1)*dr(3)*dr(3)-r2*dr(1))
          lrf(jg)%phi3(10)=lrf(jg)%phi3(10) +field2*(5.*dr(2)*dr(1)*dr(1)-r2*dr(2))
          lrf(jg)%phi3(11)=lrf(jg)%phi3(11) +field2*(5.*dr(2)*dr(1)*dr(2)-r2*dr(1))
          lrf(jg)%phi3(12)=lrf(jg)%phi3(12) +field2*(5.*dr(2)*dr(1)*dr(3))
          lrf(jg)%phi3(13)=lrf(jg)%phi3(13) +field2*(5.*dr(2)*dr(2)*dr(1)-r2*dr(1))
          lrf(jg)%phi3(14)=lrf(jg)%phi3(14) +field2*(5.*dr(2)*dr(2)*dr(2)-r2*3.*dr(2))
          lrf(jg)%phi3(15)=lrf(jg)%phi3(15) +field2*(5.*dr(2)*dr(2)*dr(3)-r2*dr(3))
          lrf(jg)%phi3(16)=lrf(jg)%phi3(16) +field2*(5.*dr(2)*dr(3)*dr(1))
          lrf(jg)%phi3(17)=lrf(jg)%phi3(17) +field2*(5.*dr(2)*dr(3)*dr(2)-r2*dr(3))
          lrf(jg)%phi3(18)=lrf(jg)%phi3(18) +field2*(5.*dr(2)*dr(3)*dr(3)-r2*dr(2))
          lrf(jg)%phi3(19)=lrf(jg)%phi3(19) +field2*(5.*dr(3)*dr(1)*dr(1)-r2*dr(3))
          lrf(jg)%phi3(20)=lrf(jg)%phi3(20) +field2*(5.*dr(3)*dr(1)*dr(2))
          lrf(jg)%phi3(21)=lrf(jg)%phi3(21) +field2*(5.*dr(3)*dr(1)*dr(3)-r2*dr(1))
          lrf(jg)%phi3(22)=lrf(jg)%phi3(22) +field2*(5.*dr(3)*dr(2)*dr(1))
          lrf(jg)%phi3(23)=lrf(jg)%phi3(23) +field2*(5.*dr(3)*dr(2)*dr(2)-r2*dr(3))
          lrf(jg)%phi3(24)=lrf(jg)%phi3(24) +field2*(5.*dr(3)*dr(2)*dr(3)-r2*dr(2))
          lrf(jg)%phi3(25)=lrf(jg)%phi3(25) +field2*(5.*dr(3)*dr(3)*dr(1)-r2*dr(1))
          lrf(jg)%phi3(26)=lrf(jg)%phi3(26) +field2*(5.*dr(3)*dr(3)*dr(2)-r2*dr(2))
          lrf(jg)%phi3(27)=lrf(jg)%phi3(27) +field2*(5.*dr(3)*dr(3)*dr(3)-r2*3.*dr(3))
        end do

        do ja = cgp(jg)%first, cgp(jg)%last

          ! skip if q-atom
          j = cgpatom(ja)
          if ( iqatom(j)/=0 ) cycle

          j3 = j*3-3

          dr(1) = x(j3+1) - lrf(ig)%cgp_cent(1)
          dr(2) = x(j3+2) - lrf(ig)%cgp_cent(2)
          dr(3) = x(j3+3) - lrf(ig)%cgp_cent(3)
          r2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)

          field0=crg(j)/(r2*sqrt(r2))
          lrf(ig)%phi0=lrf(ig)%phi0+field0*r2
          lrf(ig)%phi1(1)=lrf(ig)%phi1(1)-field0*dr(1)
          lrf(ig)%phi1(2)=lrf(ig)%phi1(2)-field0*dr(2)
          lrf(ig)%phi1(3)=lrf(ig)%phi1(3)-field0*dr(3)
          field1=3.*field0/r2
          lrf(ig)%phi2(1)=lrf(ig)%phi2(1)+field1*dr(1)*dr(1)-field0
          lrf(ig)%phi2(2)=lrf(ig)%phi2(2)+field1*dr(1)*dr(2)
          lrf(ig)%phi2(3)=lrf(ig)%phi2(3)+field1*dr(1)*dr(3)
          lrf(ig)%phi2(4)=lrf(ig)%phi2(4)+field1*dr(2)*dr(1)
          lrf(ig)%phi2(5)=lrf(ig)%phi2(5)+field1*dr(2)*dr(2)-field0
          lrf(ig)%phi2(6)=lrf(ig)%phi2(6)+field1*dr(2)*dr(3)
          lrf(ig)%phi2(7)=lrf(ig)%phi2(7)+field1*dr(3)*dr(1)
          lrf(ig)%phi2(8)=lrf(ig)%phi2(8)+field1*dr(3)*dr(2)
          lrf(ig)%phi2(9)=lrf(ig)%phi2(9)+field1*dr(3)*dr(3)-field0
          field2=-field1/r2
          lrf(ig)%phi3(1 )=lrf(ig)%phi3(1 ) +field2*(5.*dr(1)*dr(1)*dr(1)-r2*3.*dr(1))
          lrf(ig)%phi3(2 )=lrf(ig)%phi3(2 ) +field2*(5.*dr(1)*dr(1)*dr(2)-r2*dr(2))
          lrf(ig)%phi3(3 )=lrf(ig)%phi3(3 ) +field2*(5.*dr(1)*dr(1)*dr(3)-r2*dr(3))
          lrf(ig)%phi3(4 )=lrf(ig)%phi3(4 ) +field2*(5.*dr(1)*dr(2)*dr(1)-r2*dr(2))
          lrf(ig)%phi3(5 )=lrf(ig)%phi3(5 ) +field2*(5.*dr(1)*dr(2)*dr(2)-r2*dr(1))
          lrf(ig)%phi3(6 )=lrf(ig)%phi3(6 ) +field2*(5.*dr(1)*dr(2)*dr(3))
          lrf(ig)%phi3(7 )=lrf(ig)%phi3(7 ) +field2*(5.*dr(1)*dr(3)*dr(1)-r2*dr(3))
          lrf(ig)%phi3(8 )=lrf(ig)%phi3(8 ) +field2*(5.*dr(1)*dr(3)*dr(2))
          lrf(ig)%phi3(9 )=lrf(ig)%phi3(9 ) +field2*(5.*dr(1)*dr(3)*dr(3)-r2*dr(1))
          lrf(ig)%phi3(10)=lrf(ig)%phi3(10) +field2*(5.*dr(2)*dr(1)*dr(1)-r2*dr(2))
          lrf(ig)%phi3(11)=lrf(ig)%phi3(11) +field2*(5.*dr(2)*dr(1)*dr(2)-r2*dr(1))
          lrf(ig)%phi3(12)=lrf(ig)%phi3(12) +field2*(5.*dr(2)*dr(1)*dr(3))
          lrf(ig)%phi3(13)=lrf(ig)%phi3(13) +field2*(5.*dr(2)*dr(2)*dr(1)-r2*dr(1))
          lrf(ig)%phi3(14)=lrf(ig)%phi3(14) +field2*(5.*dr(2)*dr(2)*dr(2)-r2*3.*dr(2))
          lrf(ig)%phi3(15)=lrf(ig)%phi3(15) +field2*(5.*dr(2)*dr(2)*dr(3)-r2*dr(3))
          lrf(ig)%phi3(16)=lrf(ig)%phi3(16) +field2*(5.*dr(2)*dr(3)*dr(1))
          lrf(ig)%phi3(17)=lrf(ig)%phi3(17) +field2*(5.*dr(2)*dr(3)*dr(2)-r2*dr(3))
          lrf(ig)%phi3(18)=lrf(ig)%phi3(18) +field2*(5.*dr(2)*dr(3)*dr(3)-r2*dr(2))
          lrf(ig)%phi3(19)=lrf(ig)%phi3(19) +field2*(5.*dr(3)*dr(1)*dr(1)-r2*dr(3))
          lrf(ig)%phi3(20)=lrf(ig)%phi3(20) +field2*(5.*dr(3)*dr(1)*dr(2))
          lrf(ig)%phi3(21)=lrf(ig)%phi3(21) +field2*(5.*dr(3)*dr(1)*dr(3)-r2*dr(1))
          lrf(ig)%phi3(22)=lrf(ig)%phi3(22) +field2*(5.*dr(3)*dr(2)*dr(1))
          lrf(ig)%phi3(23)=lrf(ig)%phi3(23) +field2*(5.*dr(3)*dr(2)*dr(2)-r2*dr(3))
          lrf(ig)%phi3(24)=lrf(ig)%phi3(24) +field2*(5.*dr(3)*dr(2)*dr(3)-r2*dr(2))
          lrf(ig)%phi3(25)=lrf(ig)%phi3(25) +field2*(5.*dr(3)*dr(3)*dr(1)-r2*dr(1))
          lrf(ig)%phi3(26)=lrf(ig)%phi3(26) +field2*(5.*dr(3)*dr(3)*dr(2)-r2*dr(2))
          lrf(ig)%phi3(27)=lrf(ig)%phi3(27) +field2*(5.*dr(3)*dr(3)*dr(3)-r2*3.*dr(3))
        enddo

  end if ! outside cutoff

end do jgloop
end do igloop

#if defined (PROFILING)
profile(3)%time = profile(3)%time + rtime() - start_loop_time
#endif

end subroutine nbpplist_lrf
!----------------LRF version of PW PBC-----------------------
subroutine nbpplist_box_lrf
  ! local variables
  integer						:: i,j,ig,jg,ia,ja,i3,j3,nl,ig_sw, jg_sw, is3
  real(8)						:: rcut2,r2
  integer						:: LJ_code
  real(8)						:: dx, dy, dz
  
  real(8)						::RcLRF2,field0, field1, field2
  real(8)						::dr(3)
  real(8)						::boxshiftx, boxshifty, boxshiftz	

#if defined (PROFILING)
real(8)                                         :: start_loop_time
start_loop_time = rtime()
#endif

  ! For use with periodic boundary conditions
  !	This routine makes a list of non-bonded solute-solute atom pairs 
  !	excluding any Q-atoms.

  ! uses the global variables:
  !  Rcpp, ncgp, cgp, excl, x, cgpatom, iqatom, ljcod, crg, iaclib, max_nbr_range, listex
  !  nexlong, listexlong, calculation_assignment%pp%max, alloc_status, list14, n14long, list14long

  ! reset #pairs
  nbpp_pair = 0 !atom pairs
  nbpp_cgp_pair = 0 !chargegroup pairs
  rcut2 = Rcpp*Rcpp
  RcLRF2 = RcLRF*RcLRF

igloop: do ig = calculation_assignment%pp%start, calculation_assignment%pp%end
	! for every assigned charge group:
	ig_sw = cgp(ig)%iswitch !switching atom in charge group ig
	i3 = 3*ig_sw-3

jgloop: do jg = 1, ncgp_solute
	  ! for every charge group:

	  ! count each charge group pair once only
	  if ( ((ig .gt. jg) .and. (mod(ig+jg,2) .eq. 0)) .or. &
		   ((ig .lt. jg) .and. (mod(ig+jg,2) .eq. 1)) ) &
		   cycle jgloop

	  jg_sw = cgp(jg)%iswitch !switching atom in charge group jg
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

        !innside cutoff

	    !check if more memory is needed
	    if(nbpp_cgp_pair .eq. size(nbpp_cgp, 1)) call reallocate_nbpp_cgp
	    !add the charge group pair
	    nbpp_cgp_pair = nbpp_cgp_pair + 1
	    nbpp_cgp(nbpp_cgp_pair)%i = ig_sw !the switching atoms of the charge groups in the pair
	    nbpp_cgp(nbpp_cgp_pair)%j = jg_sw

ialoop: do ia = cgp(ig)%first, cgp(ig)%last
	    ! for every atom in the charge group ig (of the outermost loop):
		i = cgpatom(ia)

		! skip if q-atom
		if ( iqatom(i)/=0 ) cycle ialoop

jaloop:   do ja = cgp(jg)%first, cgp(jg)%last
		    ! for every atom in the charge group jg (innermost loop)
		    j = cgpatom(ja)

		    ! make sure each pair is only counted once
		    if ( ig .eq. jg .and. i .ge. j ) cycle jaloop

		    ! skip if q-atom
		    if ( iqatom(j)/=0 ) cycle jaloop

		    LJ_code = ljcod(iac(i),iac(j))

		    ! skip if all interactions zero
		    if((crg(i) * crg(j) == 0.) &
			  .and. &
			  (iaclib(iac(i))%avdw(LJ_code)*iaclib(iac(j))%avdw(LJ_code) == 0.) &
			  .and. &
			  (iaclib(iac(i))%bvdw(LJ_code)*iaclib(iac(j))%bvdw(LJ_code) == 0.)) &
			  cycle jaloop

		    ! skip bonded exclusions and 1-4 nbors
		    if ( abs(j-i) .le. max_nbr_range ) then
			  if ( i .lt. j ) then
			    if ( listex(j-i, i) ) cycle jaloop
			  else
			    if ( listex(i-j, j) ) cycle jaloop
			  end if
		    else
			  do nl = 1, nexlong
			    if ((listexlong(1, nl) .eq. i .and. listexlong(2, nl) .eq. j) .or. &
				  (listexlong(1, nl) .eq. j .and. listexlong(2, nl) .eq. i) ) &
				  cycle jaloop
			  end do
		    end if

		    ! if out of space then make more space
		    if (nbpp_pair .eq. calculation_assignment%pp%max) call reallocate_nonbondlist_pp

		    ! all tests passed, add the pair
		    nbpp_pair = nbpp_pair + 1
		    nbpp(nbpp_pair)%i = i
		    nbpp(nbpp_pair)%j = j 
		    nbpp(nbpp_pair)%LJcod = LJ_code
		    nbpp(nbpp_pair)%cgp_pair = nbpp_cgp_pair !which pair of charge groups the atom pair belongs to
		
		    ! set LJcod of the pair to 3 if the atoms have bonded interactions
		    if ( abs(j-i) .le. max_nbr_range ) then
			  if ( i .lt. j ) then
			    if ( list14(j-i, i) ) nbpp(nbpp_pair)%LJcod = 3
			  else
			    if ( list14(i-j, j) ) nbpp(nbpp_pair)%LJcod = 3
			  end if
		    else

		      do nl = 1, n14long
		       if ((list14long(1, nl) .eq. i .and. list14long(2, nl) .eq. j) .or. &
			      (list14long(1, nl) .eq. j .and. list14long(2, nl) .eq. i)) &
				  nbpp(nbpp_pair)%LJcod = 3
			  end do
		    end if
	
		  end do jaloop
	    end do ialoop
      
	  elseif(r2 <= RcLRF2) then
        ! outside pp-cutoff but inside LRF cut-off use LRF
		
		!ig : jg calculation
		boxshiftx = x(i3+1) - lrf(jg)%cgp_cent(1)
		boxshifty = x(i3+2) - lrf(jg)%cgp_cent(2)
		boxshiftz = x(i3+3) - lrf(jg)%cgp_cent(3)

		boxshiftx = boxlength(1)*nint( boxshiftx*inv_boxl(1) )
		boxshifty = boxlength(2)*nint( boxshifty*inv_boxl(2) )
		boxshiftz = boxlength(3)*nint( boxshiftz*inv_boxl(3) )

        do ia = cgp(ig)%first, cgp(ig)%last

          ! skip if q-atom
          i = cgpatom(ia)
          if ( iqatom(i)/=0 ) cycle

          is3 = i*3-3

          dr(1) = x(is3+1) - lrf(jg)%cgp_cent(1) - boxshiftx
          dr(2) = x(is3+2) - lrf(jg)%cgp_cent(2) - boxshifty
          dr(3) = x(is3+3) - lrf(jg)%cgp_cent(3) - boxshiftz
       	  r2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)

          field0=crg(i)/(r2*sqrt(r2))
          lrf(jg)%phi0=lrf(jg)%phi0+field0*r2
          lrf(jg)%phi1(1)=lrf(jg)%phi1(1)-field0*dr(1)
          lrf(jg)%phi1(2)=lrf(jg)%phi1(2)-field0*dr(2)
          lrf(jg)%phi1(3)=lrf(jg)%phi1(3)-field0*dr(3)
          field1=3.*field0/r2
          lrf(jg)%phi2(1)=lrf(jg)%phi2(1)+field1*dr(1)*dr(1)-field0
          lrf(jg)%phi2(2)=lrf(jg)%phi2(2)+field1*dr(1)*dr(2)
          lrf(jg)%phi2(3)=lrf(jg)%phi2(3)+field1*dr(1)*dr(3)
          lrf(jg)%phi2(4)=lrf(jg)%phi2(4)+field1*dr(2)*dr(1)
          lrf(jg)%phi2(5)=lrf(jg)%phi2(5)+field1*dr(2)*dr(2)-field0
          lrf(jg)%phi2(6)=lrf(jg)%phi2(6)+field1*dr(2)*dr(3)
          lrf(jg)%phi2(7)=lrf(jg)%phi2(7)+field1*dr(3)*dr(1)
          lrf(jg)%phi2(8)=lrf(jg)%phi2(8)+field1*dr(3)*dr(2)
          lrf(jg)%phi2(9)=lrf(jg)%phi2(9)+field1*dr(3)*dr(3)-field0
          field2=-field1/r2
          lrf(jg)%phi3(1 )=lrf(jg)%phi3(1 ) +field2*(5.*dr(1)*dr(1)*dr(1)-r2*3.*dr(1))
          lrf(jg)%phi3(2 )=lrf(jg)%phi3(2 ) +field2*(5.*dr(1)*dr(1)*dr(2)-r2*dr(2))
          lrf(jg)%phi3(3 )=lrf(jg)%phi3(3 ) +field2*(5.*dr(1)*dr(1)*dr(3)-r2*dr(3))
          lrf(jg)%phi3(4 )=lrf(jg)%phi3(4 ) +field2*(5.*dr(1)*dr(2)*dr(1)-r2*dr(2))
          lrf(jg)%phi3(5 )=lrf(jg)%phi3(5 ) +field2*(5.*dr(1)*dr(2)*dr(2)-r2*dr(1))
          lrf(jg)%phi3(6 )=lrf(jg)%phi3(6 ) +field2*(5.*dr(1)*dr(2)*dr(3))
          lrf(jg)%phi3(7 )=lrf(jg)%phi3(7 ) +field2*(5.*dr(1)*dr(3)*dr(1)-r2*dr(3))
          lrf(jg)%phi3(8 )=lrf(jg)%phi3(8 ) +field2*(5.*dr(1)*dr(3)*dr(2))
          lrf(jg)%phi3(9 )=lrf(jg)%phi3(9 ) +field2*(5.*dr(1)*dr(3)*dr(3)-r2*dr(1))
          lrf(jg)%phi3(10)=lrf(jg)%phi3(10) +field2*(5.*dr(2)*dr(1)*dr(1)-r2*dr(2))
          lrf(jg)%phi3(11)=lrf(jg)%phi3(11) +field2*(5.*dr(2)*dr(1)*dr(2)-r2*dr(1))
          lrf(jg)%phi3(12)=lrf(jg)%phi3(12) +field2*(5.*dr(2)*dr(1)*dr(3))
          lrf(jg)%phi3(13)=lrf(jg)%phi3(13) +field2*(5.*dr(2)*dr(2)*dr(1)-r2*dr(1))
          lrf(jg)%phi3(14)=lrf(jg)%phi3(14) +field2*(5.*dr(2)*dr(2)*dr(2)-r2*3.*dr(2))
          lrf(jg)%phi3(15)=lrf(jg)%phi3(15) +field2*(5.*dr(2)*dr(2)*dr(3)-r2*dr(3))
          lrf(jg)%phi3(16)=lrf(jg)%phi3(16) +field2*(5.*dr(2)*dr(3)*dr(1))
          lrf(jg)%phi3(17)=lrf(jg)%phi3(17) +field2*(5.*dr(2)*dr(3)*dr(2)-r2*dr(3))
          lrf(jg)%phi3(18)=lrf(jg)%phi3(18) +field2*(5.*dr(2)*dr(3)*dr(3)-r2*dr(2))
          lrf(jg)%phi3(19)=lrf(jg)%phi3(19) +field2*(5.*dr(3)*dr(1)*dr(1)-r2*dr(3))
          lrf(jg)%phi3(20)=lrf(jg)%phi3(20) +field2*(5.*dr(3)*dr(1)*dr(2))
          lrf(jg)%phi3(21)=lrf(jg)%phi3(21) +field2*(5.*dr(3)*dr(1)*dr(3)-r2*dr(1))
          lrf(jg)%phi3(22)=lrf(jg)%phi3(22) +field2*(5.*dr(3)*dr(2)*dr(1))
          lrf(jg)%phi3(23)=lrf(jg)%phi3(23) +field2*(5.*dr(3)*dr(2)*dr(2)-r2*dr(3))
          lrf(jg)%phi3(24)=lrf(jg)%phi3(24) +field2*(5.*dr(3)*dr(2)*dr(3)-r2*dr(2))
          lrf(jg)%phi3(25)=lrf(jg)%phi3(25) +field2*(5.*dr(3)*dr(3)*dr(1)-r2*dr(1))
          lrf(jg)%phi3(26)=lrf(jg)%phi3(26) +field2*(5.*dr(3)*dr(3)*dr(2)-r2*dr(2))
          lrf(jg)%phi3(27)=lrf(jg)%phi3(27) +field2*(5.*dr(3)*dr(3)*dr(3)-r2*3.*dr(3))
        end do

		!jg : ig calculations
		boxshiftx = x(j3+1) - lrf(ig)%cgp_cent(1)
		boxshifty = x(j3+2) - lrf(ig)%cgp_cent(2)
		boxshiftz = x(j3+3) - lrf(ig)%cgp_cent(3)
	
		boxshiftx = boxlength(1)*nint( boxshiftx*inv_boxl(1) )
		boxshifty = boxlength(2)*nint( boxshifty*inv_boxl(2) )
		boxshiftz = boxlength(3)*nint( boxshiftz*inv_boxl(3) )

	    do ja = cgp(jg)%first, cgp(jg)%last

          ! skip if q-atom
          j = cgpatom(ja)
          if ( iqatom(j)/=0 ) cycle

          j3 = j*3-3

          dr(1) = x(j3+1) - lrf(ig)%cgp_cent(1) - boxshiftx
          dr(2) = x(j3+2) - lrf(ig)%cgp_cent(2) - boxshifty
          dr(3) = x(j3+3) - lrf(ig)%cgp_cent(3) - boxshiftz
          r2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)

          field0=crg(j)/(r2*sqrt(r2))
          lrf(ig)%phi0=lrf(ig)%phi0+field0*r2
          lrf(ig)%phi1(1)=lrf(ig)%phi1(1)-field0*dr(1)
          lrf(ig)%phi1(2)=lrf(ig)%phi1(2)-field0*dr(2)
          lrf(ig)%phi1(3)=lrf(ig)%phi1(3)-field0*dr(3)
          field1=3.*field0/r2
          lrf(ig)%phi2(1)=lrf(ig)%phi2(1)+field1*dr(1)*dr(1)-field0
          lrf(ig)%phi2(2)=lrf(ig)%phi2(2)+field1*dr(1)*dr(2)
          lrf(ig)%phi2(3)=lrf(ig)%phi2(3)+field1*dr(1)*dr(3)
          lrf(ig)%phi2(4)=lrf(ig)%phi2(4)+field1*dr(2)*dr(1)
          lrf(ig)%phi2(5)=lrf(ig)%phi2(5)+field1*dr(2)*dr(2)-field0
          lrf(ig)%phi2(6)=lrf(ig)%phi2(6)+field1*dr(2)*dr(3)
          lrf(ig)%phi2(7)=lrf(ig)%phi2(7)+field1*dr(3)*dr(1)
          lrf(ig)%phi2(8)=lrf(ig)%phi2(8)+field1*dr(3)*dr(2)
          lrf(ig)%phi2(9)=lrf(ig)%phi2(9)+field1*dr(3)*dr(3)-field0
          field2=-field1/r2
          lrf(ig)%phi3(1 )=lrf(ig)%phi3(1 ) +field2*(5.*dr(1)*dr(1)*dr(1)-r2*3.*dr(1))
          lrf(ig)%phi3(2 )=lrf(ig)%phi3(2 ) +field2*(5.*dr(1)*dr(1)*dr(2)-r2*dr(2))
          lrf(ig)%phi3(3 )=lrf(ig)%phi3(3 ) +field2*(5.*dr(1)*dr(1)*dr(3)-r2*dr(3))
          lrf(ig)%phi3(4 )=lrf(ig)%phi3(4 ) +field2*(5.*dr(1)*dr(2)*dr(1)-r2*dr(2))
          lrf(ig)%phi3(5 )=lrf(ig)%phi3(5 ) +field2*(5.*dr(1)*dr(2)*dr(2)-r2*dr(1))
          lrf(ig)%phi3(6 )=lrf(ig)%phi3(6 ) +field2*(5.*dr(1)*dr(2)*dr(3))
          lrf(ig)%phi3(7 )=lrf(ig)%phi3(7 ) +field2*(5.*dr(1)*dr(3)*dr(1)-r2*dr(3))
          lrf(ig)%phi3(8 )=lrf(ig)%phi3(8 ) +field2*(5.*dr(1)*dr(3)*dr(2))
          lrf(ig)%phi3(9 )=lrf(ig)%phi3(9 ) +field2*(5.*dr(1)*dr(3)*dr(3)-r2*dr(1))
          lrf(ig)%phi3(10)=lrf(ig)%phi3(10) +field2*(5.*dr(2)*dr(1)*dr(1)-r2*dr(2))
          lrf(ig)%phi3(11)=lrf(ig)%phi3(11) +field2*(5.*dr(2)*dr(1)*dr(2)-r2*dr(1))
          lrf(ig)%phi3(12)=lrf(ig)%phi3(12) +field2*(5.*dr(2)*dr(1)*dr(3))
          lrf(ig)%phi3(13)=lrf(ig)%phi3(13) +field2*(5.*dr(2)*dr(2)*dr(1)-r2*dr(1))
          lrf(ig)%phi3(14)=lrf(ig)%phi3(14) +field2*(5.*dr(2)*dr(2)*dr(2)-r2*3.*dr(2))
          lrf(ig)%phi3(15)=lrf(ig)%phi3(15) +field2*(5.*dr(2)*dr(2)*dr(3)-r2*dr(3))
          lrf(ig)%phi3(16)=lrf(ig)%phi3(16) +field2*(5.*dr(2)*dr(3)*dr(1))
          lrf(ig)%phi3(17)=lrf(ig)%phi3(17) +field2*(5.*dr(2)*dr(3)*dr(2)-r2*dr(3))
          lrf(ig)%phi3(18)=lrf(ig)%phi3(18) +field2*(5.*dr(2)*dr(3)*dr(3)-r2*dr(2))
          lrf(ig)%phi3(19)=lrf(ig)%phi3(19) +field2*(5.*dr(3)*dr(1)*dr(1)-r2*dr(3))
          lrf(ig)%phi3(20)=lrf(ig)%phi3(20) +field2*(5.*dr(3)*dr(1)*dr(2))
          lrf(ig)%phi3(21)=lrf(ig)%phi3(21) +field2*(5.*dr(3)*dr(1)*dr(3)-r2*dr(1))
          lrf(ig)%phi3(22)=lrf(ig)%phi3(22) +field2*(5.*dr(3)*dr(2)*dr(1))
          lrf(ig)%phi3(23)=lrf(ig)%phi3(23) +field2*(5.*dr(3)*dr(2)*dr(2)-r2*dr(3))
          lrf(ig)%phi3(24)=lrf(ig)%phi3(24) +field2*(5.*dr(3)*dr(2)*dr(3)-r2*dr(2))
          lrf(ig)%phi3(25)=lrf(ig)%phi3(25) +field2*(5.*dr(3)*dr(3)*dr(1)-r2*dr(1))
          lrf(ig)%phi3(26)=lrf(ig)%phi3(26) +field2*(5.*dr(3)*dr(3)*dr(2)-r2*dr(2))
          lrf(ig)%phi3(27)=lrf(ig)%phi3(27) +field2*(5.*dr(3)*dr(3)*dr(3)-r2*3.*dr(3))
        enddo

      end if ! outside cutoff

	end do jgloop
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
real(8)						:: rcut2,r2
integer						:: LJ_code

!******PWadded variables 2001-10-01

real(8)						:: dx, dy, dz

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
ja = nat_solute + 3*jg-2
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

jaloop:	do ja = nat_solute + 3*jg-2, nat_solute + 3*jg
          ! for every atom of the water molecule:

          ! calculate LJ_code for the pair
          LJ_code = ljcod(iac(i),iac(ja))

          ! skip pairs with zero interaction
          if((crg(i) * crg(ja) == 0.) &
                  .and. &
                  (iaclib(iac(i))%avdw(LJ_code)*iaclib(iac(ja))%avdw(LJ_code) == 0.) &
                  .and. &
                  (iaclib(iac(i))%bvdw(LJ_code)*iaclib(iac(ja))%bvdw(LJ_code) == 0.)) &
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
integer						:: i,ig,jg,ia,ja,i3,j3, inside
real(8)						:: rcut2,r2
#if defined (PROFILING)
real(8)						:: start_loop_time
start_loop_time = rtime()
#endif

!	This routine makes a list of non-bonded solute-solvent atom pairs
!	excluding Q-atoms.


nbpw_pair = 0
rcut2 = Rcpw*Rcpw


igloop:  do ig = calculation_assignment%pw%start, calculation_assignment%pw%end

!	   --- excluded group ? ---
ia = cgp(ig)%iswitch
if ( excl(ia) ) cycle igloop

jgloop:     do jg = 1, nwat
ja = nat_solute + 3*jg-2
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

ialoop:    do ia = cgp(ig)%first, cgp(ig)%last

   !	         --- q-atom ? ---
   i = cgpatom(ia)
   if ( iqatom(i)/=0 ) cycle ialoop

           ! if out of space then make more space
           if (nbpw_pair > calculation_assignment%pw%max-3) call reallocate_nonbondlist_pw

   nbpw_pair = nbpw_pair + 3
   nbpw(nbpw_pair-2)%i = i
   nbpw(nbpw_pair-2)%j = ja 
   nbpw(nbpw_pair-2)%LJcod = ljcod(iac(i),iac(ja))
   nbpw(nbpw_pair-1)%i = i
   nbpw(nbpw_pair-1)%j = ja+1
   nbpw(nbpw_pair-1)%LJcod = ljcod(iac(i),iac(ja+1))
   nbpw(nbpw_pair  )%i = i
   nbpw(nbpw_pair  )%j = ja+2
   nbpw(nbpw_pair  )%LJcod = ljcod(iac(i),iac(ja+2))

end do ialoop
end do jgloop
end do igloop
#if defined (PROFILING)
profile(4)%time = profile(4)%time + rtime() - start_loop_time
#endif

end subroutine nbpwlis2

!------------------------------------------------------------------------------

!******PWadded 2001-10-18
subroutine nbpwlis2_box
  ! local variables
  integer						:: i,ig,jg,ia,ja,i3,j3,ig_atom, inside
  real(8)						:: rcut2,r2
  real(8)						:: dx, dy,dz
  
  ! for periodic boundary conditions
  !	This routine makes a list of non-bonded solute-solvent atom pairs
  !	excluding Q-atoms.
#if defined (PROFILING)
real(8)						:: start_loop_time
start_loop_time = rtime()
#endif

  nbpw_pair = 0
  nbpw_cgp_pair = 0
  rcut2 = Rcpw*Rcpw

igloop:  do ig = calculation_assignment%pw%start, calculation_assignment%pw%end

jgloop:     do jg = 1, nwat

        	ja = nat_solute + 3*jg-2
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
			 if (nbpw_cgp_pair .eq. size(nbpw_cgp, 1)) call reallocate_nbpw_cgp
			 nbpw_cgp_pair = nbpw_cgp_pair + 1
			 nbpw_cgp(nbpw_cgp_pair)%i = i
			 nbpw_cgp(nbpw_cgp_pair)%j = ja
		   end if

		  ig_atom = ig_atom + 1 !ia = ia + 1
        end do
		if (inside .eq. 0) cycle jgloop

ialoop:    do ia = cgp(ig)%first, cgp(ig)%last

           !	         --- q-atom ? ---
           i = cgpatom(ia)
           if ( iqatom(i)/=0 ) cycle ialoop

		   ! if out of space then make more space
		   if (nbpw_pair > calculation_assignment%pw%max-3) call reallocate_nonbondlist_pw

           nbpw_pair = nbpw_pair + 3
           nbpw(nbpw_pair-2)%i = i
           nbpw(nbpw_pair-2)%j = ja 
           nbpw(nbpw_pair-2)%LJcod = ljcod(iac(i),iac(ja))
		   nbpw(nbpw_pair-2)%cgp_pair = nbpw_cgp_pair
           nbpw(nbpw_pair-1)%i = i
           nbpw(nbpw_pair-1)%j = ja+1
           nbpw(nbpw_pair-1)%LJcod = ljcod(iac(i),iac(ja+1))
		   nbpw(nbpw_pair-1)%cgp_pair = nbpw_cgp_pair
           nbpw(nbpw_pair  )%i = i
           nbpw(nbpw_pair  )%j = ja+2
           nbpw(nbpw_pair  )%LJcod = ljcod(iac(i),iac(ja+2))
		   nbpw(nbpw_pair  )%cgp_pair = nbpw_cgp_pair

        end do ialoop
     end do jgloop
  end do igloop
#if defined (PROFILING)
profile(4)%time = profile(4)%time + rtime() - start_loop_time
#endif

end subroutine nbpwlis2_box
!-----------------------------------------------------------------------
subroutine nbpwlis2_box_lrf
  ! local variables
  integer						:: i,ig,jg,ia,ja,i3,j3,ig_atom, inside
  real(8)						:: rcut2,r2
  real(8)						:: dx, dy,dz
  !LRF
  real(8)						:: RcLRF2, field0, field1, field2
  real(8)						:: dr(3)
  integer						:: jg_cgp, j, inside_LRF, is3
  real(8)						:: boxshiftx, boxshifty, boxshiftz
  
  ! for periodic boundary conditions
  !	This routine makes a list of non-bonded solute-solvent atom pairs
  !	excluding Q-atoms.
#if defined (PROFILING)
real(8)						:: start_loop_time
start_loop_time = rtime()
#endif

  nbpw_pair = 0
  nbpw_cgp_pair = 0
  rcut2 = Rcpw*Rcpw
  RcLRF2 = RcLRF*RcLRF

igloop:  do ig = calculation_assignment%pw%start, calculation_assignment%pw%end
jgloop:    do jg = 1, nwat
        		ja = nat_solute + 3*jg-2
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
					   if (nbpw_cgp_pair .eq. size(nbpw_cgp, 1)) call reallocate_nbpw_cgp
					   nbpw_cgp_pair = nbpw_cgp_pair + 1
					   nbpw_cgp(nbpw_cgp_pair)%i = i
					   nbpw_cgp(nbpw_cgp_pair)%j = ja
					 elseif (r2 <= RcLRF2) then
    			     inside_LRF = 1
			         end if

		             ig_atom = ig_atom + 1 !ia = ia + 1
                     end do

	  if (inside .eq. 1) then

ialoop: do ia = cgp(ig)%first, cgp(ig)%last

           !	         --- q-atom ? ---
           i = cgpatom(ia)
           if ( iqatom(i)/=0 ) cycle ialoop

		   ! if out of space then make more space
		   if (nbpw_pair > calculation_assignment%pw%max-3) call reallocate_nonbondlist_pw

           nbpw_pair = nbpw_pair + 3
           nbpw(nbpw_pair-2)%i = i
           nbpw(nbpw_pair-2)%j = ja 
           nbpw(nbpw_pair-2)%LJcod = ljcod(iac(i),iac(ja))
		   nbpw(nbpw_pair-2)%cgp_pair = nbpw_cgp_pair
           nbpw(nbpw_pair-1)%i = i
           nbpw(nbpw_pair-1)%j = ja+1
           nbpw(nbpw_pair-1)%LJcod = ljcod(iac(i),iac(ja+1))
		   nbpw(nbpw_pair-1)%cgp_pair = nbpw_cgp_pair
           nbpw(nbpw_pair  )%i = i
           nbpw(nbpw_pair  )%j = ja+2
           nbpw(nbpw_pair  )%LJcod = ljcod(iac(i),iac(ja+2))
		   nbpw(nbpw_pair  )%cgp_pair = nbpw_cgp_pair

        end do ialoop
      elseif((inside_LRF ==1) .and. (inside == 0)) then   
                ! outside pw-cutoff but inside LRF cut-off: use LRF
		!solut : solvent
		jg_cgp = iwhich_cgp(ja)
		
		boxshiftx = x(i3+1) - lrf(jg_cgp)%cgp_cent(1)
		boxshifty = x(i3+2) - lrf(jg_cgp)%cgp_cent(2)
		boxshiftz = x(i3+3) - lrf(jg_cgp)%cgp_cent(3)

		boxshiftx = boxlength(1)*nint( boxshiftx*inv_boxl(1) )
		boxshifty = boxlength(2)*nint( boxshifty*inv_boxl(2) )
		boxshiftz = boxlength(3)*nint( boxshiftz*inv_boxl(3) )

ialoop2: do ia = cgp(ig)%first, cgp(ig)%last
          i = cgpatom(ia)

          ! skip q-atoms
          if ( iqatom(i)/=0 ) cycle
		  is3 = i*3-3
		  
		  !jg = ncgp + iw
   		  ! calculate dr and (d)r2
		  dr(1) = x(is3+1) - lrf(jg_cgp)%cgp_cent(1) - boxshiftx
		  dr(2) = x(is3+2) - lrf(jg_cgp)%cgp_cent(2) - boxshifty
		  dr(3) = x(is3+3) - lrf(jg_cgp)%cgp_cent(3) - boxshiftz
	   	  		  	  	
		  r2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)

          ! calculate lrf parameters for the charge group
		  field0=crg(i)/(r2*sqrt(r2))
		  lrf(jg_cgp)%phi0=lrf(jg_cgp)%phi0+field0*r2
		  lrf(jg_cgp)%phi1(1)=lrf(jg_cgp)%phi1(1)-field0*dr(1)
		  lrf(jg_cgp)%phi1(2)=lrf(jg_cgp)%phi1(2)-field0*dr(2)
          lrf(jg_cgp)%phi1(3)=lrf(jg_cgp)%phi1(3)-field0*dr(3)
		  field1=3.*field0/r2
		  lrf(jg_cgp)%phi2(1)=lrf(jg_cgp)%phi2(1)+field1*dr(1)*dr(1)-field0
		  lrf(jg_cgp)%phi2(2)=lrf(jg_cgp)%phi2(2)+field1*dr(1)*dr(2)
		  lrf(jg_cgp)%phi2(3)=lrf(jg_cgp)%phi2(3)+field1*dr(1)*dr(3)
		  lrf(jg_cgp)%phi2(4)=lrf(jg_cgp)%phi2(4)+field1*dr(2)*dr(1)
		  lrf(jg_cgp)%phi2(5)=lrf(jg_cgp)%phi2(5)+field1*dr(2)*dr(2)-field0
		  lrf(jg_cgp)%phi2(6)=lrf(jg_cgp)%phi2(6)+field1*dr(2)*dr(3)
		  lrf(jg_cgp)%phi2(7)=lrf(jg_cgp)%phi2(7)+field1*dr(3)*dr(1)
		  lrf(jg_cgp)%phi2(8)=lrf(jg_cgp)%phi2(8)+field1*dr(3)*dr(2)
		  lrf(jg_cgp)%phi2(9)=lrf(jg_cgp)%phi2(9)+field1*dr(3)*dr(3)-field0
		  field2=-field1/r2
		  lrf(jg_cgp)%phi3(1 )=lrf(jg_cgp)%phi3(1 )+field2*(5.*dr(1)*dr(1)*dr(1)-r2*3.*dr(1))
		  lrf(jg_cgp)%phi3(2 )=lrf(jg_cgp)%phi3(2 )+field2*(5.*dr(1)*dr(1)*dr(2)-r2*dr(2))
		  lrf(jg_cgp)%phi3(3 )=lrf(jg_cgp)%phi3(3 )+field2*(5.*dr(1)*dr(1)*dr(3)-r2*dr(3))
		  lrf(jg_cgp)%phi3(4 )=lrf(jg_cgp)%phi3(4 )+field2*(5.*dr(1)*dr(2)*dr(1)-r2*dr(2))
		  lrf(jg_cgp)%phi3(5 )=lrf(jg_cgp)%phi3(5 )+field2*(5.*dr(1)*dr(2)*dr(2)-r2*dr(1))
		  lrf(jg_cgp)%phi3(6 )=lrf(jg_cgp)%phi3(6 )+field2*(5.*dr(1)*dr(2)*dr(3))
		  lrf(jg_cgp)%phi3(7 )=lrf(jg_cgp)%phi3(7 )+field2*(5.*dr(1)*dr(3)*dr(1)-r2*dr(3))
		  lrf(jg_cgp)%phi3(8 )=lrf(jg_cgp)%phi3(8 )+field2*(5.*dr(1)*dr(3)*dr(2))
		  lrf(jg_cgp)%phi3(9 )=lrf(jg_cgp)%phi3(9 )+field2*(5.*dr(1)*dr(3)*dr(3)-r2*dr(1))
		  lrf(jg_cgp)%phi3(10)=lrf(jg_cgp)%phi3(10)+field2*(5.*dr(2)*dr(1)*dr(1)-r2*dr(2))
		  lrf(jg_cgp)%phi3(11)=lrf(jg_cgp)%phi3(11)+field2*(5.*dr(2)*dr(1)*dr(2)-r2*dr(1))
		  lrf(jg_cgp)%phi3(12)=lrf(jg_cgp)%phi3(12)+field2*(5.*dr(2)*dr(1)*dr(3))
		  lrf(jg_cgp)%phi3(13)=lrf(jg_cgp)%phi3(13)+field2*(5.*dr(2)*dr(2)*dr(1)-r2*dr(1))
		  lrf(jg_cgp)%phi3(14)=lrf(jg_cgp)%phi3(14)+field2*(5.*dr(2)*dr(2)*dr(2)-r2*3.*dr(2))
		  lrf(jg_cgp)%phi3(15)=lrf(jg_cgp)%phi3(15)+field2*(5.*dr(2)*dr(2)*dr(3)-r2*dr(3))
		  lrf(jg_cgp)%phi3(16)=lrf(jg_cgp)%phi3(16)+field2*(5.*dr(2)*dr(3)*dr(1))
		  lrf(jg_cgp)%phi3(17)=lrf(jg_cgp)%phi3(17)+field2*(5.*dr(2)*dr(3)*dr(2)-r2*dr(3))
		  lrf(jg_cgp)%phi3(18)=lrf(jg_cgp)%phi3(18)+field2*(5.*dr(2)*dr(3)*dr(3)-r2*dr(2))
		  lrf(jg_cgp)%phi3(19)=lrf(jg_cgp)%phi3(19)+field2*(5.*dr(3)*dr(1)*dr(1)-r2*dr(3))
		  lrf(jg_cgp)%phi3(20)=lrf(jg_cgp)%phi3(20)+field2*(5.*dr(3)*dr(1)*dr(2))
		  lrf(jg_cgp)%phi3(21)=lrf(jg_cgp)%phi3(21)+field2*(5.*dr(3)*dr(1)*dr(3)-r2*dr(1))
		  lrf(jg_cgp)%phi3(22)=lrf(jg_cgp)%phi3(22)+field2*(5.*dr(3)*dr(2)*dr(1))
		  lrf(jg_cgp)%phi3(23)=lrf(jg_cgp)%phi3(23)+field2*(5.*dr(3)*dr(2)*dr(2)-r2*dr(3))
		  lrf(jg_cgp)%phi3(24)=lrf(jg_cgp)%phi3(24)+field2*(5.*dr(3)*dr(2)*dr(3)-r2*dr(2))
		  lrf(jg_cgp)%phi3(25)=lrf(jg_cgp)%phi3(25)+field2*(5.*dr(3)*dr(3)*dr(1)-r2*dr(1))
		  lrf(jg_cgp)%phi3(26)=lrf(jg_cgp)%phi3(26)+field2*(5.*dr(3)*dr(3)*dr(2)-r2*dr(2))
		  lrf(jg_cgp)%phi3(27)=lrf(jg_cgp)%phi3(27)+field2*(5.*dr(3)*dr(3)*dr(3)-r2*3.*dr(3))
		 enddo ialoop2

		!solvent : solut	
		boxshiftx = x(j3+1) - lrf(ig)%cgp_cent(1)
		boxshifty = x(j3+2) - lrf(ig)%cgp_cent(2)
		boxshiftz = x(j3+3) - lrf(ig)%cgp_cent(3)

		boxshiftx = boxlength(1)*nint( boxshiftx*inv_boxl(1) )
		boxshifty = boxlength(2)*nint( boxshifty*inv_boxl(2) )
		boxshiftz = boxlength(3)*nint( boxshiftz*inv_boxl(3) )

jaloop2: do ja = 1, 3
		  j = nat_solute + 3*jg-3 + ja
		  j3 = j*3-3

		  ! calculate dr and (d)r2
		  dr(1) = x(j3+1) - lrf(ig)%cgp_cent(1) - boxshiftx
		  dr(2) = x(j3+2) - lrf(ig)%cgp_cent(2) - boxshifty
		  dr(3) = x(j3+3) - lrf(ig)%cgp_cent(3) - boxshiftz
		  
	      r2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)

		  ! calculate lrf for the water molecule
		  field0=crg(j)/(r2*sqrt(r2))
		  lrf(ig)%phi0=lrf(ig)%phi0+field0*r2
		  lrf(ig)%phi1(1)=lrf(ig)%phi1(1)-field0*dr(1)
		  lrf(ig)%phi1(2)=lrf(ig)%phi1(2)-field0*dr(2)
		  lrf(ig)%phi1(3)=lrf(ig)%phi1(3)-field0*dr(3)
		  field1=3.*field0/r2
		  lrf(ig)%phi2(1)=lrf(ig)%phi2(1)+field1*dr(1)*dr(1)-field0
		  lrf(ig)%phi2(2)=lrf(ig)%phi2(2)+field1*dr(1)*dr(2)
		  lrf(ig)%phi2(3)=lrf(ig)%phi2(3)+field1*dr(1)*dr(3)
		  lrf(ig)%phi2(4)=lrf(ig)%phi2(4)+field1*dr(2)*dr(1)
		  lrf(ig)%phi2(5)=lrf(ig)%phi2(5)+field1*dr(2)*dr(2)-field0
		  lrf(ig)%phi2(6)=lrf(ig)%phi2(6)+field1*dr(2)*dr(3)
		  lrf(ig)%phi2(7)=lrf(ig)%phi2(7)+field1*dr(3)*dr(1)
		  lrf(ig)%phi2(8)=lrf(ig)%phi2(8)+field1*dr(3)*dr(2)
		  lrf(ig)%phi2(9)=lrf(ig)%phi2(9)+field1*dr(3)*dr(3)-field0
		  field2=-field1/r2
		  lrf(ig)%phi3(1 )=lrf(ig)%phi3(1 )+field2*(5.*dr(1)*dr(1)*dr(1)-r2*3.*dr(1))
		  lrf(ig)%phi3(2 )=lrf(ig)%phi3(2 )+field2*(5.*dr(1)*dr(1)*dr(2)-r2*dr(2))
		  lrf(ig)%phi3(3 )=lrf(ig)%phi3(3 )+field2*(5.*dr(1)*dr(1)*dr(3)-r2*dr(3))
		  lrf(ig)%phi3(4 )=lrf(ig)%phi3(4 )+field2*(5.*dr(1)*dr(2)*dr(1)-r2*dr(2))
		  lrf(ig)%phi3(5 )=lrf(ig)%phi3(5 )+field2*(5.*dr(1)*dr(2)*dr(2)-r2*dr(1))
		  lrf(ig)%phi3(6 )=lrf(ig)%phi3(6 )+field2*(5.*dr(1)*dr(2)*dr(3))
		  lrf(ig)%phi3(7 )=lrf(ig)%phi3(7 )+field2*(5.*dr(1)*dr(3)*dr(1)-r2*dr(3))
		  lrf(ig)%phi3(8 )=lrf(ig)%phi3(8 )+field2*(5.*dr(1)*dr(3)*dr(2))
		  lrf(ig)%phi3(9 )=lrf(ig)%phi3(9 )+field2*(5.*dr(1)*dr(3)*dr(3)-r2*dr(1))
		  lrf(ig)%phi3(10)=lrf(ig)%phi3(10)+field2*(5.*dr(2)*dr(1)*dr(1)-r2*dr(2))
		  lrf(ig)%phi3(11)=lrf(ig)%phi3(11)+field2*(5.*dr(2)*dr(1)*dr(2)-r2*dr(1))
		  lrf(ig)%phi3(12)=lrf(ig)%phi3(12)+field2*(5.*dr(2)*dr(1)*dr(3))
		  lrf(ig)%phi3(13)=lrf(ig)%phi3(13)+field2*(5.*dr(2)*dr(2)*dr(1)-r2*dr(1))
		  lrf(ig)%phi3(14)=lrf(ig)%phi3(14)+field2*(5.*dr(2)*dr(2)*dr(2)-r2*3.*dr(2))
		  lrf(ig)%phi3(15)=lrf(ig)%phi3(15)+field2*(5.*dr(2)*dr(2)*dr(3)-r2*dr(3))
		  lrf(ig)%phi3(16)=lrf(ig)%phi3(16)+field2*(5.*dr(2)*dr(3)*dr(1))
		  lrf(ig)%phi3(17)=lrf(ig)%phi3(17)+field2*(5.*dr(2)*dr(3)*dr(2)-r2*dr(3))
		  lrf(ig)%phi3(18)=lrf(ig)%phi3(18)+field2*(5.*dr(2)*dr(3)*dr(3)-r2*dr(2))
		  lrf(ig)%phi3(19)=lrf(ig)%phi3(19)+field2*(5.*dr(3)*dr(1)*dr(1)-r2*dr(3))
		  lrf(ig)%phi3(20)=lrf(ig)%phi3(20)+field2*(5.*dr(3)*dr(1)*dr(2))
		  lrf(ig)%phi3(21)=lrf(ig)%phi3(21)+field2*(5.*dr(3)*dr(1)*dr(3)-r2*dr(1))
		  lrf(ig)%phi3(22)=lrf(ig)%phi3(22)+field2*(5.*dr(3)*dr(2)*dr(1))
		  lrf(ig)%phi3(23)=lrf(ig)%phi3(23)+field2*(5.*dr(3)*dr(2)*dr(2)-r2*dr(3))
		  lrf(ig)%phi3(24)=lrf(ig)%phi3(24)+field2*(5.*dr(3)*dr(2)*dr(3)-r2*dr(2))
		  lrf(ig)%phi3(25)=lrf(ig)%phi3(25)+field2*(5.*dr(3)*dr(3)*dr(1)-r2*dr(1))
		  lrf(ig)%phi3(26)=lrf(ig)%phi3(26)+field2*(5.*dr(3)*dr(3)*dr(2)-r2*dr(2))
		  lrf(ig)%phi3(27)=lrf(ig)%phi3(27)+field2*(5.*dr(3)*dr(3)*dr(3)-r2*3.*dr(3))
         enddo jaloop2
      end if

     end do jgloop
  end do igloop
#if defined (PROFILING)
profile(4)%time = profile(4)%time + rtime() - start_loop_time
#endif

end subroutine nbpwlis2_box_lrf
!-----------------------------------------------------------------------
subroutine nbpwlis2_lrf
! local variables
integer						:: i,j,ig,iw,jg,ia,ja,i3,j3,inside,is,is3
real(8)						:: rcut2,r2,field0,field1,field2
real(8)						:: dr(3)
real(8)						::	RcLRF2
!	This routine makes a list of non-bonded solute-solvent atom pairs
!	excluding Q-atoms.
#if defined (PROFILING)
real(8)						:: start_loop_time
start_loop_time = rtime()
#endif

nbpw_pair = 0
rcut2 = Rcpw*Rcpw
RcLRF2 = RcLRF*RcLRF

igloop:  do ig = calculation_assignment%pw%start, calculation_assignment%pw%end

!	   --- excluded group ? ---
is = cgp(ig)%iswitch
if ( excl(is) ) cycle igloop
is3 = 3*is-3

iwloop:     do iw = 1, nwat

ja = nat_solute + 3*iw-2
if(excl(ja)) cycle iwloop ! skip excluded waters
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

if ( inside .eq. 1 ) then

ialoop:       do ia = cgp(ig)%first, cgp(ig)%last

      !	            --- q-atom ? ---
      i = cgpatom(ia)
      if ( iqatom(i)/=0 ) cycle ialoop

                  ! if out of space then make more space
                  if (nbpw_pair > calculation_assignment%pw%max-3) call reallocate_nonbondlist_pw

      nbpw_pair = nbpw_pair + 3
      nbpw(nbpw_pair-2)%i = i
      nbpw(nbpw_pair-2)%j = ja 
      nbpw(nbpw_pair-2)%LJcod = ljcod(iac(i),iac(ja))
      nbpw(nbpw_pair-1)%i = i
      nbpw(nbpw_pair-1)%j = ja+1
      nbpw(nbpw_pair-1)%LJcod = ljcod(iac(i),iac(ja+1))
      nbpw(nbpw_pair  )%i = i
      nbpw(nbpw_pair  )%j = ja+2
      nbpw(nbpw_pair  )%LJcod = ljcod(iac(i),iac(ja+2))

   end do ialoop

elseif(r2 <= RcLRF2) then   
                ! outside pw-cutoff but inside LRF cut-off: use LRF

ialoop2:      do ia = cgp(ig)%first, cgp(ig)%last

      !	             --- q-atom ? ---
      i = cgpatom(ia)
      if ( iqatom(i)/=0 ) cycle ialoop2

      i3 = i*3-3
      !jg = ncgp + iw
                  jg = iwhich_cgp(ja)

      dr(1) = x(i3+1) - lrf(jg)%cgp_cent(1)
      dr(2) = x(i3+2) - lrf(jg)%cgp_cent(2)
      dr(3) = x(i3+3) - lrf(jg)%cgp_cent(3)
      r2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)

      field0=crg(i)/(r2*sqrt(r2))
      lrf(jg)%phi0=lrf(jg)%phi0+field0*r2
      lrf(jg)%phi1(1)=lrf(jg)%phi1(1)-field0*dr(1)
      lrf(jg)%phi1(2)=lrf(jg)%phi1(2)-field0*dr(2)
      lrf(jg)%phi1(3)=lrf(jg)%phi1(3)-field0*dr(3)
      field1=3.*field0/r2
      lrf(jg)%phi2(1)=lrf(jg)%phi2(1)+field1*dr(1)*dr(1)-field0
      lrf(jg)%phi2(2)=lrf(jg)%phi2(2)+field1*dr(1)*dr(2)
      lrf(jg)%phi2(3)=lrf(jg)%phi2(3)+field1*dr(1)*dr(3)
      lrf(jg)%phi2(4)=lrf(jg)%phi2(4)+field1*dr(2)*dr(1)
      lrf(jg)%phi2(5)=lrf(jg)%phi2(5)+field1*dr(2)*dr(2)-field0
      lrf(jg)%phi2(6)=lrf(jg)%phi2(6)+field1*dr(2)*dr(3)
      lrf(jg)%phi2(7)=lrf(jg)%phi2(7)+field1*dr(3)*dr(1)
      lrf(jg)%phi2(8)=lrf(jg)%phi2(8)+field1*dr(3)*dr(2)
      lrf(jg)%phi2(9)=lrf(jg)%phi2(9)+field1*dr(3)*dr(3)-field0
      field2=-field1/r2
      lrf(jg)%phi3(1 )=lrf(jg)%phi3(1 ) +field2*(5.*dr(1)*dr(1)*dr(1)-r2*3.*dr(1))
      lrf(jg)%phi3(2 )=lrf(jg)%phi3(2 ) +field2*(5.*dr(1)*dr(1)*dr(2)-r2*dr(2))
      lrf(jg)%phi3(3 )=lrf(jg)%phi3(3 ) +field2*(5.*dr(1)*dr(1)*dr(3)-r2*dr(3))
      lrf(jg)%phi3(4 )=lrf(jg)%phi3(4 ) +field2*(5.*dr(1)*dr(2)*dr(1)-r2*dr(2))
      lrf(jg)%phi3(5 )=lrf(jg)%phi3(5 ) +field2*(5.*dr(1)*dr(2)*dr(2)-r2*dr(1))
      lrf(jg)%phi3(6 )=lrf(jg)%phi3(6 ) +field2*(5.*dr(1)*dr(2)*dr(3))
      lrf(jg)%phi3(7 )=lrf(jg)%phi3(7 ) +field2*(5.*dr(1)*dr(3)*dr(1)-r2*dr(3))
      lrf(jg)%phi3(8 )=lrf(jg)%phi3(8 ) +field2*(5.*dr(1)*dr(3)*dr(2))
      lrf(jg)%phi3(9 )=lrf(jg)%phi3(9 ) +field2*(5.*dr(1)*dr(3)*dr(3)-r2*dr(1))
      lrf(jg)%phi3(10)=lrf(jg)%phi3(10) +field2*(5.*dr(2)*dr(1)*dr(1)-r2*dr(2))
      lrf(jg)%phi3(11)=lrf(jg)%phi3(11) +field2*(5.*dr(2)*dr(1)*dr(2)-r2*dr(1))
      lrf(jg)%phi3(12)=lrf(jg)%phi3(12) +field2*(5.*dr(2)*dr(1)*dr(3))
      lrf(jg)%phi3(13)=lrf(jg)%phi3(13) +field2*(5.*dr(2)*dr(2)*dr(1)-r2*dr(1))
      lrf(jg)%phi3(14)=lrf(jg)%phi3(14) +field2*(5.*dr(2)*dr(2)*dr(2)-r2*3.*dr(2))
      lrf(jg)%phi3(15)=lrf(jg)%phi3(15) +field2*(5.*dr(2)*dr(2)*dr(3)-r2*dr(3))
      lrf(jg)%phi3(16)=lrf(jg)%phi3(16) +field2*(5.*dr(2)*dr(3)*dr(1))
      lrf(jg)%phi3(17)=lrf(jg)%phi3(17) +field2*(5.*dr(2)*dr(3)*dr(2)-r2*dr(3))
      lrf(jg)%phi3(18)=lrf(jg)%phi3(18) +field2*(5.*dr(2)*dr(3)*dr(3)-r2*dr(2))
      lrf(jg)%phi3(19)=lrf(jg)%phi3(19) +field2*(5.*dr(3)*dr(1)*dr(1)-r2*dr(3))
      lrf(jg)%phi3(20)=lrf(jg)%phi3(20) +field2*(5.*dr(3)*dr(1)*dr(2))
      lrf(jg)%phi3(21)=lrf(jg)%phi3(21) +field2*(5.*dr(3)*dr(1)*dr(3)-r2*dr(1))
      lrf(jg)%phi3(22)=lrf(jg)%phi3(22) +field2*(5.*dr(3)*dr(2)*dr(1))
      lrf(jg)%phi3(23)=lrf(jg)%phi3(23) +field2*(5.*dr(3)*dr(2)*dr(2)-r2*dr(3))
      lrf(jg)%phi3(24)=lrf(jg)%phi3(24) +field2*(5.*dr(3)*dr(2)*dr(3)-r2*dr(2))
      lrf(jg)%phi3(25)=lrf(jg)%phi3(25) +field2*(5.*dr(3)*dr(3)*dr(1)-r2*dr(1))
      lrf(jg)%phi3(26)=lrf(jg)%phi3(26) +field2*(5.*dr(3)*dr(3)*dr(2)-r2*dr(2))
      lrf(jg)%phi3(27)=lrf(jg)%phi3(27) 	+field2*(5.*dr(3)*dr(3)*dr(3)-r2*3.*dr(3))

   end do ialoop2

jaloop2:   do ja = 1, 3

      j = nat_solute + 3*iw-3 + ja

      j3 = j*3-3

      dr(1) = x(j3+1) - lrf(ig)%cgp_cent(1)
      dr(2) = x(j3+2) - lrf(ig)%cgp_cent(2)
      dr(3) = x(j3+3) - lrf(ig)%cgp_cent(3)
      r2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)

      field0=crg(j)/(r2*sqrt(r2))
      lrf(ig)%phi0=lrf(ig)%phi0+field0*r2
      lrf(ig)%phi1(1)=lrf(ig)%phi1(1)-field0*dr(1)
      lrf(ig)%phi1(2)=lrf(ig)%phi1(2)-field0*dr(2)
      lrf(ig)%phi1(3)=lrf(ig)%phi1(3)-field0*dr(3)
      field1=3.*field0/r2
      lrf(ig)%phi2(1)=lrf(ig)%phi2(1)+field1*dr(1)*dr(1)-field0
      lrf(ig)%phi2(2)=lrf(ig)%phi2(2)+field1*dr(1)*dr(2)
      lrf(ig)%phi2(3)=lrf(ig)%phi2(3)+field1*dr(1)*dr(3)
      lrf(ig)%phi2(4)=lrf(ig)%phi2(4)+field1*dr(2)*dr(1)
      lrf(ig)%phi2(5)=lrf(ig)%phi2(5)+field1*dr(2)*dr(2)-field0
      lrf(ig)%phi2(6)=lrf(ig)%phi2(6)+field1*dr(2)*dr(3)
      lrf(ig)%phi2(7)=lrf(ig)%phi2(7)+field1*dr(3)*dr(1)
      lrf(ig)%phi2(8)=lrf(ig)%phi2(8)+field1*dr(3)*dr(2)
      lrf(ig)%phi2(9)=lrf(ig)%phi2(9)+field1*dr(3)*dr(3)-field0
      field2=-field1/r2
      lrf(ig)%phi3(1 )=lrf(ig)%phi3(1 ) +field2*(5.*dr(1)*dr(1)*dr(1)-r2*3.*dr(1))
      lrf(ig)%phi3(2 )=lrf(ig)%phi3(2 ) +field2*(5.*dr(1)*dr(1)*dr(2)-r2*dr(2))
      lrf(ig)%phi3(3 )=lrf(ig)%phi3(3 ) +field2*(5.*dr(1)*dr(1)*dr(3)-r2*dr(3))
      lrf(ig)%phi3(4 )=lrf(ig)%phi3(4 ) +field2*(5.*dr(1)*dr(2)*dr(1)-r2*dr(2))
      lrf(ig)%phi3(5 )=lrf(ig)%phi3(5 ) +field2*(5.*dr(1)*dr(2)*dr(2)-r2*dr(1))
      lrf(ig)%phi3(6 )=lrf(ig)%phi3(6 ) +field2*(5.*dr(1)*dr(2)*dr(3))
      lrf(ig)%phi3(7 )=lrf(ig)%phi3(7 ) +field2*(5.*dr(1)*dr(3)*dr(1)-r2*dr(3))
      lrf(ig)%phi3(8 )=lrf(ig)%phi3(8 ) +field2*(5.*dr(1)*dr(3)*dr(2))
      lrf(ig)%phi3(9 )=lrf(ig)%phi3(9 ) +field2*(5.*dr(1)*dr(3)*dr(3)-r2*dr(1))
      lrf(ig)%phi3(10)=lrf(ig)%phi3(10) +field2*(5.*dr(2)*dr(1)*dr(1)-r2*dr(2))
      lrf(ig)%phi3(11)=lrf(ig)%phi3(11) +field2*(5.*dr(2)*dr(1)*dr(2)-r2*dr(1))
      lrf(ig)%phi3(12)=lrf(ig)%phi3(12) +field2*(5.*dr(2)*dr(1)*dr(3))
      lrf(ig)%phi3(13)=lrf(ig)%phi3(13) +field2*(5.*dr(2)*dr(2)*dr(1)-r2*dr(1))
      lrf(ig)%phi3(14)=lrf(ig)%phi3(14) +field2*(5.*dr(2)*dr(2)*dr(2)-r2*3.*dr(2))
      lrf(ig)%phi3(15)=lrf(ig)%phi3(15) +field2*(5.*dr(2)*dr(2)*dr(3)-r2*dr(3))
      lrf(ig)%phi3(16)=lrf(ig)%phi3(16) +field2*(5.*dr(2)*dr(3)*dr(1))
      lrf(ig)%phi3(17)=lrf(ig)%phi3(17) +field2*(5.*dr(2)*dr(3)*dr(2)-r2*dr(3))
      lrf(ig)%phi3(18)=lrf(ig)%phi3(18) +field2*(5.*dr(2)*dr(3)*dr(3)-r2*dr(2))
      lrf(ig)%phi3(19)=lrf(ig)%phi3(19) +field2*(5.*dr(3)*dr(1)*dr(1)-r2*dr(3))
      lrf(ig)%phi3(20)=lrf(ig)%phi3(20) +field2*(5.*dr(3)*dr(1)*dr(2))
      lrf(ig)%phi3(21)=lrf(ig)%phi3(21) +field2*(5.*dr(3)*dr(1)*dr(3)-r2*dr(1))
      lrf(ig)%phi3(22)=lrf(ig)%phi3(22) +field2*(5.*dr(3)*dr(2)*dr(1))
      lrf(ig)%phi3(23)=lrf(ig)%phi3(23) +field2*(5.*dr(3)*dr(2)*dr(2)-r2*dr(3))
      lrf(ig)%phi3(24)=lrf(ig)%phi3(24) +field2*(5.*dr(3)*dr(2)*dr(3)-r2*dr(2))
      lrf(ig)%phi3(25)=lrf(ig)%phi3(25) +field2*(5.*dr(3)*dr(3)*dr(1)-r2*dr(1))
      lrf(ig)%phi3(26)=lrf(ig)%phi3(26) +field2*(5.*dr(3)*dr(3)*dr(2)-r2*dr(2))
      lrf(ig)%phi3(27)=lrf(ig)%phi3(27) +field2*(5.*dr(3)*dr(3)*dr(3)-r2*3.*dr(3))

   end do jaloop2

end if

end do iwloop

end do igloop
#if defined (PROFILING)
profile(4)%time = profile(4)%time + rtime() - start_loop_time
#endif

end subroutine nbpwlis2_lrf

!-----------------------------------------------------------------------

subroutine nbpwlist
! local variables
integer						:: i,ig,jg,ia,ja,i3,j3
real(8)						:: rcut2,r2
integer						:: LJ_code


! for use with spherical boundary

!	This routine makes a list of non-bonded solute-solvent atom pairs
!	excluding Q-atoms.
#if defined (PROFILING)
real(8)						:: start_loop_time
start_loop_time = rtime()
#endif

! reset nbpw_pair
nbpw_pair = 0
rcut2 = Rcpw*Rcpw

igloop: do ig = calculation_assignment%pw%start, calculation_assignment%pw%end
! for every charge group:

! skip excluded groups
ia = cgp(ig)%iswitch
if ( excl(ia) ) cycle igloop
i3 = 3*ia-3

jgloop: do jg = 1, nwat
  ! for every water molecule:

ja = nat_solute + 3*jg-2
if(excl(ja)) cycle jgloop ! skip excluded waters

  j3 = 3*ja-3
  r2 = ( x(i3+1) -x(j3+1) )**2 &
  +( x(i3+2) -x(j3+2) )**2 &
  +( x(i3+3) -x(j3+3) )**2

  ! skip water outside cutoff
if ( r2 .gt. rcut2 ) cycle jgloop

ialoop: do ia = cgp(ig)%first, cgp(ig)%last
    ! for every atom in the charge group:

        ! find the atom index of the atom in the charge group
i = cgpatom(ia)

!	skip q-atoms
        if ( iqatom(i)/=0 ) cycle ialoop

        ! if out of space then make more space
        if (nbpw_pair .gt. calculation_assignment%pw%max - 3) call reallocate_nonbondlist_pw

jaloop:	do ja = nat_solute + 3*jg-2, nat_solute + 3*jg
          ! for every atom of the water molecule:

          ! calculate LJ_code for the pair
          LJ_code = ljcod(iac(i),iac(ja))

          ! skip pairs with zero interaction
          if((crg(i) * crg(ja) == 0.) &
                  .and. &
                  (iaclib(iac(i))%avdw(LJ_code)*iaclib(iac(ja))%avdw(LJ_code) == 0.) &
                  .and. &
                  (iaclib(iac(i))%bvdw(LJ_code)*iaclib(iac(ja))%bvdw(LJ_code) == 0.)) &
                  cycle jaloop

          ! add the pair
          nbpw_pair = nbpw_pair + 1
          nbpw(nbpw_pair)%i = i
          nbpw(nbpw_pair)%j = ja 
          nbpw(nbpw_pair)%LJcod = LJ_code
        end do jaloop
end do ialoop
end do jgloop
end do igloop
#if defined (PROFILING)
profile(4)%time = profile(4)%time + rtime() - start_loop_time
#endif

end subroutine nbpwlist

!-----------------------------------------------------------------------
!******PWadded 2001-10-18

subroutine nbpwlist_box
  ! local variables
  integer						:: i,ig,jg,ia,ja,i3,j3,ig_sw,jg_sw
  real(8)						:: rcut2,r2
  integer						:: LJ_code
  real(8)						:: dx, dy, dz

  ! For use with periodic boundary conditions
  !	This routine makes a list of non-bonded solute-solvent atom pairs
  !	excluding Q-atoms.
#if defined (PROFILING)
real(8)						:: start_loop_time
start_loop_time = rtime()
#endif

  ! reset nbpw_pair
  nbpw_pair = 0
  nbpw_cgp_pair = 0
  rcut2 = Rcpw*Rcpw

igloop: do ig = calculation_assignment%pw%start, calculation_assignment%pw%end
	! for every charge group:
    ig_sw = cgp(ig)%iswitch
    i3 = 3*ig_sw-3

jgloop: do jg = 1, nwat
	  ! for every water molecule:

      jg_sw = nat_solute + 3*jg-2

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
	  if(nbpw_cgp_pair .eq. size(nbpw_cgp, 1) ) call reallocate_nbpw_cgp

	  nbpw_cgp_pair = nbpw_cgp_pair + 1
	  nbpw_cgp(nbpw_cgp_pair)%i = ig_sw  !solute
	  nbpw_cgp(nbpw_cgp_pair)%j = jg_sw  !water

ialoop: do ia = cgp(ig)%first, cgp(ig)%last
	    ! for every atom in the charge group:

		! find the atom index of the atom in the charge group
        i = cgpatom(ia)

       !	skip q-atoms
		if ( iqatom(i)/=0 ) cycle ialoop

		! if out of space then make more space
		if (nbpw_pair .gt. calculation_assignment%pw%max - 3) call reallocate_nonbondlist_pw
                   
jaloop:	do ja = nat_solute + 3*jg-2, nat_solute + 3*jg
  		  ! for every atom of the water molecule:

		  ! calculate LJ_code for the pair
		  LJ_code = ljcod(iac(i),iac(ja))

		  ! skip pairs with zero interaction
		  if((crg(i) * crg(ja) == 0.) &
  			  .and. &
			  (iaclib(iac(i))%avdw(LJ_code)*iaclib(iac(ja))%avdw(LJ_code) == 0.) &
			  .and. &
			  (iaclib(iac(i))%bvdw(LJ_code)*iaclib(iac(ja))%bvdw(LJ_code) == 0.)) &
			  cycle jaloop

		  ! add the pair
		  nbpw_pair = nbpw_pair + 1
		  nbpw(nbpw_pair)%i = i
		  nbpw(nbpw_pair)%j = ja 
		  nbpw(nbpw_pair)%LJcod = LJ_code
		  nbpw(nbpw_pair)%cgp_pair = nbpw_cgp_pair
		  !write(*,'(a, i8, a, i7)') 'pw-pair ', nbpw_pair, ' belongs to charge group ' , nbpw(nbpw_pair)%cgp_pair
		end do jaloop
      end do ialoop
    end do jgloop
  end do igloop
#if defined (PROFILING)
profile(4)%time = profile(4)%time + rtime() - start_loop_time
#endif

end subroutine nbpwlist_box
!-------------------------------------------------------------------------------------

subroutine nbpwlist_lrf
! local variables
integer						:: i,j,ig,iw,jg,ia,ja,i3,j3,is,is3
real(8)						:: rcut2,r2,field0,field1,field2
real(8)						:: dr(3)
integer						:: LJ_code
real(8)						::	RcLRF2

#if defined (PROFILING)
real(8)						:: start_loop_time
start_loop_time = rtime()
#endif

!	This routine makes a list of non-bonded solute-solvent atom pairs
!	excluding Q-atoms.

! reset nbpw_pair
nbpw_pair = 0
rcut2 = Rcpw*Rcpw
RcLRF2 = RcLRF*RcLRF

igloop: do ig = calculation_assignment%pw%start, calculation_assignment%pw%end
! for every charge group:

! skip excluded groups
is = cgp(ig)%iswitch
if ( excl(is) ) cycle igloop
is3 = 3*is-3

iwloop: do iw = 1, nwat
  ! for every water molecule:

ja = nat_solute + 3*iw-2
if(excl(ja)) cycle iwloop ! skip excluded waters
j3 = 3*ja-3

r2 = ( x(is3+1) -x(j3+1) )**2 &
     +( x(is3+2) -x(j3+2) )**2 &
     +( x(is3+3) -x(j3+3) )**2

if ( r2 .le. rcut2 ) then
    ! within the cutoff radix:

ialoop: do ia = cgp(ig)%first, cgp(ig)%last
  i = cgpatom(ia)

  ! skip q-atoms
  if ( iqatom(i)/=0 ) cycle ialoop

          ! if out of space then make more space
          if (nbpw_pair .gt. calculation_assignment%pw%max - 3) call reallocate_nonbondlist_pw

jaloop:	  do ja = nat_solute + 3*iw-2, nat_solute + 3*iw
                ! calculate the LJ_code of the pair
                LJ_code = ljcod(iac(i),iac(ja))

                ! skip pairs with zero interaction
                if((crg(i) * crg(ja) == 0.) &
                        .and. &
                        (iaclib(iac(i))%avdw(LJ_code)*iaclib(iac(ja))%avdw(LJ_code) == 0.) &
                        .and. &
                        (iaclib(iac(i))%bvdw(LJ_code)*iaclib(iac(ja))%bvdw(LJ_code) == 0.)) &
                        cycle jaloop

                ! add the pair
                nbpw_pair = nbpw_pair + 1
                nbpw(nbpw_pair)%i = i
                nbpw(nbpw_pair)%j = ja 
                nbpw(nbpw_pair)%LJcod = LJ_code
!TMP                nbpw_per_cgp(ig) = nbpw_per_cgp(ig) + 1
          end do jaloop
end do ialoop

elseif(r2 <= RcLRF2) then   
                ! outside pw-cutoff but inside LRF cut-off: use LRF

ialoop2: do ia = cgp(ig)%first, cgp(ig)%last
  i = cgpatom(ia)

  ! skip q-atoms
  if ( iqatom(i)/=0 ) cycle

  i3 = i*3-3
  !jg = ncgp + iw
          jg = iwhich_cgp(ja)

          ! calculate dr and (d)r2
  dr(1) = x(i3+1) - lrf(jg)%cgp_cent(1)
  dr(2) = x(i3+2) - lrf(jg)%cgp_cent(2)
  dr(3) = x(i3+3) - lrf(jg)%cgp_cent(3)
  r2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)

          ! calculate lrf parameters for the charge group
  field0=crg(i)/(r2*sqrt(r2))
  lrf(jg)%phi0=lrf(jg)%phi0+field0*r2
  lrf(jg)%phi1(1)=lrf(jg)%phi1(1)-field0*dr(1)
  lrf(jg)%phi1(2)=lrf(jg)%phi1(2)-field0*dr(2)
  lrf(jg)%phi1(3)=lrf(jg)%phi1(3)-field0*dr(3)
  field1=3.*field0/r2
  lrf(jg)%phi2(1)=lrf(jg)%phi2(1)+field1*dr(1)*dr(1)-field0
  lrf(jg)%phi2(2)=lrf(jg)%phi2(2)+field1*dr(1)*dr(2)
  lrf(jg)%phi2(3)=lrf(jg)%phi2(3)+field1*dr(1)*dr(3)
  lrf(jg)%phi2(4)=lrf(jg)%phi2(4)+field1*dr(2)*dr(1)
  lrf(jg)%phi2(5)=lrf(jg)%phi2(5)+field1*dr(2)*dr(2)-field0
  lrf(jg)%phi2(6)=lrf(jg)%phi2(6)+field1*dr(2)*dr(3)
  lrf(jg)%phi2(7)=lrf(jg)%phi2(7)+field1*dr(3)*dr(1)
  lrf(jg)%phi2(8)=lrf(jg)%phi2(8)+field1*dr(3)*dr(2)
  lrf(jg)%phi2(9)=lrf(jg)%phi2(9)+field1*dr(3)*dr(3)-field0
  field2=-field1/r2
  lrf(jg)%phi3(1 )=lrf(jg)%phi3(1 )+field2*(5.*dr(1)*dr(1)*dr(1)-r2*3.*dr(1))
  lrf(jg)%phi3(2 )=lrf(jg)%phi3(2 )+field2*(5.*dr(1)*dr(1)*dr(2)-r2*dr(2))
  lrf(jg)%phi3(3 )=lrf(jg)%phi3(3 )+field2*(5.*dr(1)*dr(1)*dr(3)-r2*dr(3))
  lrf(jg)%phi3(4 )=lrf(jg)%phi3(4 )+field2*(5.*dr(1)*dr(2)*dr(1)-r2*dr(2))
  lrf(jg)%phi3(5 )=lrf(jg)%phi3(5 )+field2*(5.*dr(1)*dr(2)*dr(2)-r2*dr(1))
  lrf(jg)%phi3(6 )=lrf(jg)%phi3(6 )+field2*(5.*dr(1)*dr(2)*dr(3))
  lrf(jg)%phi3(7 )=lrf(jg)%phi3(7 )+field2*(5.*dr(1)*dr(3)*dr(1)-r2*dr(3))
  lrf(jg)%phi3(8 )=lrf(jg)%phi3(8 )+field2*(5.*dr(1)*dr(3)*dr(2))
  lrf(jg)%phi3(9 )=lrf(jg)%phi3(9 )+field2*(5.*dr(1)*dr(3)*dr(3)-r2*dr(1))
  lrf(jg)%phi3(10)=lrf(jg)%phi3(10)+field2*(5.*dr(2)*dr(1)*dr(1)-r2*dr(2))
  lrf(jg)%phi3(11)=lrf(jg)%phi3(11)+field2*(5.*dr(2)*dr(1)*dr(2)-r2*dr(1))
  lrf(jg)%phi3(12)=lrf(jg)%phi3(12)+field2*(5.*dr(2)*dr(1)*dr(3))
  lrf(jg)%phi3(13)=lrf(jg)%phi3(13)+field2*(5.*dr(2)*dr(2)*dr(1)-r2*dr(1))
  lrf(jg)%phi3(14)=lrf(jg)%phi3(14)+field2*(5.*dr(2)*dr(2)*dr(2)-r2*3.*dr(2))
  lrf(jg)%phi3(15)=lrf(jg)%phi3(15)+field2*(5.*dr(2)*dr(2)*dr(3)-r2*dr(3))
  lrf(jg)%phi3(16)=lrf(jg)%phi3(16)+field2*(5.*dr(2)*dr(3)*dr(1))
  lrf(jg)%phi3(17)=lrf(jg)%phi3(17)+field2*(5.*dr(2)*dr(3)*dr(2)-r2*dr(3))
  lrf(jg)%phi3(18)=lrf(jg)%phi3(18)+field2*(5.*dr(2)*dr(3)*dr(3)-r2*dr(2))
  lrf(jg)%phi3(19)=lrf(jg)%phi3(19)+field2*(5.*dr(3)*dr(1)*dr(1)-r2*dr(3))
  lrf(jg)%phi3(20)=lrf(jg)%phi3(20)+field2*(5.*dr(3)*dr(1)*dr(2))
  lrf(jg)%phi3(21)=lrf(jg)%phi3(21)+field2*(5.*dr(3)*dr(1)*dr(3)-r2*dr(1))
  lrf(jg)%phi3(22)=lrf(jg)%phi3(22)+field2*(5.*dr(3)*dr(2)*dr(1))
  lrf(jg)%phi3(23)=lrf(jg)%phi3(23)+field2*(5.*dr(3)*dr(2)*dr(2)-r2*dr(3))
  lrf(jg)%phi3(24)=lrf(jg)%phi3(24)+field2*(5.*dr(3)*dr(2)*dr(3)-r2*dr(2))
  lrf(jg)%phi3(25)=lrf(jg)%phi3(25)+field2*(5.*dr(3)*dr(3)*dr(1)-r2*dr(1))
  lrf(jg)%phi3(26)=lrf(jg)%phi3(26)+field2*(5.*dr(3)*dr(3)*dr(2)-r2*dr(2))
  lrf(jg)%phi3(27)=lrf(jg)%phi3(27)+field2*(5.*dr(3)*dr(3)*dr(3)-r2*3.*dr(3))
enddo ialoop2

jaloop2: do ja = 1, 3
  j = nat_solute + 3*iw-3 + ja
  j3 = j*3-3

          ! calculate dr and (d)r2
  dr(1) = x(j3+1) - lrf(ig)%cgp_cent(1)
  dr(2) = x(j3+2) - lrf(ig)%cgp_cent(2)
  dr(3) = x(j3+3) - lrf(ig)%cgp_cent(3)
  
  r2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)

          ! calculate lrf for the water molecule
  field0=crg(j)/(r2*sqrt(r2))
  lrf(ig)%phi0=lrf(ig)%phi0+field0*r2
  lrf(ig)%phi1(1)=lrf(ig)%phi1(1)-field0*dr(1)
  lrf(ig)%phi1(2)=lrf(ig)%phi1(2)-field0*dr(2)
  lrf(ig)%phi1(3)=lrf(ig)%phi1(3)-field0*dr(3)
  field1=3.*field0/r2
  lrf(ig)%phi2(1)=lrf(ig)%phi2(1)+field1*dr(1)*dr(1)-field0
  lrf(ig)%phi2(2)=lrf(ig)%phi2(2)+field1*dr(1)*dr(2)
  lrf(ig)%phi2(3)=lrf(ig)%phi2(3)+field1*dr(1)*dr(3)
  lrf(ig)%phi2(4)=lrf(ig)%phi2(4)+field1*dr(2)*dr(1)
  lrf(ig)%phi2(5)=lrf(ig)%phi2(5)+field1*dr(2)*dr(2)-field0
  lrf(ig)%phi2(6)=lrf(ig)%phi2(6)+field1*dr(2)*dr(3)
  lrf(ig)%phi2(7)=lrf(ig)%phi2(7)+field1*dr(3)*dr(1)
  lrf(ig)%phi2(8)=lrf(ig)%phi2(8)+field1*dr(3)*dr(2)
  lrf(ig)%phi2(9)=lrf(ig)%phi2(9)+field1*dr(3)*dr(3)-field0
  field2=-field1/r2
  lrf(ig)%phi3(1 )=lrf(ig)%phi3(1 )+field2*(5.*dr(1)*dr(1)*dr(1)-r2*3.*dr(1))
  lrf(ig)%phi3(2 )=lrf(ig)%phi3(2 )+field2*(5.*dr(1)*dr(1)*dr(2)-r2*dr(2))
  lrf(ig)%phi3(3 )=lrf(ig)%phi3(3 )+field2*(5.*dr(1)*dr(1)*dr(3)-r2*dr(3))
  lrf(ig)%phi3(4 )=lrf(ig)%phi3(4 )+field2*(5.*dr(1)*dr(2)*dr(1)-r2*dr(2))
  lrf(ig)%phi3(5 )=lrf(ig)%phi3(5 )+field2*(5.*dr(1)*dr(2)*dr(2)-r2*dr(1))
  lrf(ig)%phi3(6 )=lrf(ig)%phi3(6 )+field2*(5.*dr(1)*dr(2)*dr(3))
  lrf(ig)%phi3(7 )=lrf(ig)%phi3(7 )+field2*(5.*dr(1)*dr(3)*dr(1)-r2*dr(3))
  lrf(ig)%phi3(8 )=lrf(ig)%phi3(8 )+field2*(5.*dr(1)*dr(3)*dr(2))
  lrf(ig)%phi3(9 )=lrf(ig)%phi3(9 )+field2*(5.*dr(1)*dr(3)*dr(3)-r2*dr(1))
  lrf(ig)%phi3(10)=lrf(ig)%phi3(10)+field2*(5.*dr(2)*dr(1)*dr(1)-r2*dr(2))
  lrf(ig)%phi3(11)=lrf(ig)%phi3(11)+field2*(5.*dr(2)*dr(1)*dr(2)-r2*dr(1))
  lrf(ig)%phi3(12)=lrf(ig)%phi3(12)+field2*(5.*dr(2)*dr(1)*dr(3))
  lrf(ig)%phi3(13)=lrf(ig)%phi3(13)+field2*(5.*dr(2)*dr(2)*dr(1)-r2*dr(1))
  lrf(ig)%phi3(14)=lrf(ig)%phi3(14)+field2*(5.*dr(2)*dr(2)*dr(2)-r2*3.*dr(2))
  lrf(ig)%phi3(15)=lrf(ig)%phi3(15)+field2*(5.*dr(2)*dr(2)*dr(3)-r2*dr(3))
  lrf(ig)%phi3(16)=lrf(ig)%phi3(16)+field2*(5.*dr(2)*dr(3)*dr(1))
  lrf(ig)%phi3(17)=lrf(ig)%phi3(17)+field2*(5.*dr(2)*dr(3)*dr(2)-r2*dr(3))
  lrf(ig)%phi3(18)=lrf(ig)%phi3(18)+field2*(5.*dr(2)*dr(3)*dr(3)-r2*dr(2))
  lrf(ig)%phi3(19)=lrf(ig)%phi3(19)+field2*(5.*dr(3)*dr(1)*dr(1)-r2*dr(3))
  lrf(ig)%phi3(20)=lrf(ig)%phi3(20)+field2*(5.*dr(3)*dr(1)*dr(2))
  lrf(ig)%phi3(21)=lrf(ig)%phi3(21)+field2*(5.*dr(3)*dr(1)*dr(3)-r2*dr(1))
  lrf(ig)%phi3(22)=lrf(ig)%phi3(22)+field2*(5.*dr(3)*dr(2)*dr(1))
  lrf(ig)%phi3(23)=lrf(ig)%phi3(23)+field2*(5.*dr(3)*dr(2)*dr(2)-r2*dr(3))
  lrf(ig)%phi3(24)=lrf(ig)%phi3(24)+field2*(5.*dr(3)*dr(2)*dr(3)-r2*dr(2))
  lrf(ig)%phi3(25)=lrf(ig)%phi3(25)+field2*(5.*dr(3)*dr(3)*dr(1)-r2*dr(1))
  lrf(ig)%phi3(26)=lrf(ig)%phi3(26)+field2*(5.*dr(3)*dr(3)*dr(2)-r2*dr(2))
  lrf(ig)%phi3(27)=lrf(ig)%phi3(27)+field2*(5.*dr(3)*dr(3)*dr(3)-r2*3.*dr(3))
enddo jaloop2

end if

end do iwloop 
end do igloop

#if defined (PROFILING)
profile(4)%time = profile(4)%time + rtime() - start_loop_time
#endif

end subroutine nbpwlist_lrf
!---------------LRF version of PW PBC-----------------------
subroutine nbpwlist_box_lrf
  ! local variables
  integer						:: i,ig,jg,ia,ja,i3,j3,ig_sw,jg_sw
  real(8)						:: rcut2,r2
  integer						:: LJ_code
  real(8)						:: dx, dy, dz
  ! LRF
  real(8)						:: RcLRF2, field0, field1, field2
  real(8)						:: dr(3)
  integer						:: jg_cgp, j, is3
  real(8)						:: boxshiftx, boxshifty, boxshiftz

#if defined (PROFILING)
real(8)                                         :: start_loop_time
start_loop_time = rtime()
#endif

  ! For use with periodic boundary conditions
  !	This routine makes a list of non-bonded solute-solvent atom pairs
  !	excluding Q-atoms.

  ! reset nbpw_pair
  nbpw_pair = 0
  nbpw_cgp_pair = 0
  rcut2 = Rcpw*Rcpw
  RcLRF2 = RcLRF*RcLRF

igloop: do ig = calculation_assignment%pw%start, calculation_assignment%pw%end
	! for every charge group:
    ig_sw = cgp(ig)%iswitch
    i3 = 3*ig_sw-3

jgloop: do jg = 1, nwat
	  ! for every water molecule:

      jg_sw = nat_solute + 3*jg-2

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
	  if(nbpw_cgp_pair .eq. size(nbpw_cgp, 1) ) call reallocate_nbpw_cgp

	  nbpw_cgp_pair = nbpw_cgp_pair + 1
	  nbpw_cgp(nbpw_cgp_pair)%i = ig_sw  !solute
	  nbpw_cgp(nbpw_cgp_pair)%j = jg_sw  !water

ialoop: do ia = cgp(ig)%first, cgp(ig)%last
	      ! for every atom in the charge group:

		  ! find the atom index of the atom in the charge group
          i = cgpatom(ia)

          !	skip q-atoms
		  if ( iqatom(i)/=0 ) cycle ialoop

		  ! if out of space then make more space
		  if (nbpw_pair .gt. calculation_assignment%pw%max - 3) call reallocate_nonbondlist_pw
                   
jaloop:	  do ja = nat_solute + 3*jg-2, nat_solute + 3*jg
  		    ! for every atom of the water molecule:

		    ! calculate LJ_code for the pair
		    LJ_code = ljcod(iac(i),iac(ja))

		    ! skip pairs with zero interaction
		    if((crg(i) * crg(ja) == 0.) &
  			    .and. &
			    (iaclib(iac(i))%avdw(LJ_code)*iaclib(iac(ja))%avdw(LJ_code) == 0.) &
			    .and. &
			    (iaclib(iac(i))%bvdw(LJ_code)*iaclib(iac(ja))%bvdw(LJ_code) == 0.)) &
			    cycle jaloop

		    ! add the pair
		    nbpw_pair = nbpw_pair + 1
		    nbpw(nbpw_pair)%i = i
		    nbpw(nbpw_pair)%j = ja 
		    nbpw(nbpw_pair)%LJcod = LJ_code
		    nbpw(nbpw_pair)%cgp_pair = nbpw_cgp_pair
		    !write(*,'(a, i8, a, i7)') 'pw-pair ', nbpw_pair, ' belongs to charge group ' , nbpw(nbpw_pair)%cgp_pair
		  end do jaloop
        end do ialoop

      elseif(r2 <= RcLRF2) then   
                ! outside pw-cutoff but inside LRF cut-off: use LRF
		!solut : solvent
		jg_cgp = iwhich_cgp(jg_sw)
		
		boxshiftx = x(i3+1) - lrf(jg_cgp)%cgp_cent(1)
		boxshifty = x(i3+2) - lrf(jg_cgp)%cgp_cent(2)
		boxshiftz = x(i3+3) - lrf(jg_cgp)%cgp_cent(3)

		boxshiftx = boxlength(1)*nint( boxshiftx*inv_boxl(1) )
		boxshifty = boxlength(2)*nint( boxshifty*inv_boxl(2) )
		boxshiftz = boxlength(3)*nint( boxshiftz*inv_boxl(3) )

ialoop2: do ia = cgp(ig)%first, cgp(ig)%last
          i = cgpatom(ia)

          ! skip q-atoms
          if ( iqatom(i)/=0 ) cycle
		  is3 = i*3-3
		  
		  !jg = ncgp + iw
          
		  ! calculate dr and (d)r2
		  dr(1) = x(is3+1) - lrf(jg_cgp)%cgp_cent(1) - boxshiftx
		  dr(2) = x(is3+2) - lrf(jg_cgp)%cgp_cent(2) - boxshifty
		  dr(3) = x(is3+3) - lrf(jg_cgp)%cgp_cent(3) - boxshiftz
	   	  		  	  	
		  r2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)

          ! calculate lrf parameters for the charge group
		  field0=crg(i)/(r2*sqrt(r2))
		  lrf(jg_cgp)%phi0=lrf(jg_cgp)%phi0+field0*r2
		  lrf(jg_cgp)%phi1(1)=lrf(jg_cgp)%phi1(1)-field0*dr(1)
		  lrf(jg_cgp)%phi1(2)=lrf(jg_cgp)%phi1(2)-field0*dr(2)
          lrf(jg_cgp)%phi1(3)=lrf(jg_cgp)%phi1(3)-field0*dr(3)
		  field1=3.*field0/r2
		  lrf(jg_cgp)%phi2(1)=lrf(jg_cgp)%phi2(1)+field1*dr(1)*dr(1)-field0
		  lrf(jg_cgp)%phi2(2)=lrf(jg_cgp)%phi2(2)+field1*dr(1)*dr(2)
		  lrf(jg_cgp)%phi2(3)=lrf(jg_cgp)%phi2(3)+field1*dr(1)*dr(3)
		  lrf(jg_cgp)%phi2(4)=lrf(jg_cgp)%phi2(4)+field1*dr(2)*dr(1)
		  lrf(jg_cgp)%phi2(5)=lrf(jg_cgp)%phi2(5)+field1*dr(2)*dr(2)-field0
		  lrf(jg_cgp)%phi2(6)=lrf(jg_cgp)%phi2(6)+field1*dr(2)*dr(3)
		  lrf(jg_cgp)%phi2(7)=lrf(jg_cgp)%phi2(7)+field1*dr(3)*dr(1)
		  lrf(jg_cgp)%phi2(8)=lrf(jg_cgp)%phi2(8)+field1*dr(3)*dr(2)
		  lrf(jg_cgp)%phi2(9)=lrf(jg_cgp)%phi2(9)+field1*dr(3)*dr(3)-field0
		  field2=-field1/r2
		  lrf(jg_cgp)%phi3(1 )=lrf(jg_cgp)%phi3(1 )+field2*(5.*dr(1)*dr(1)*dr(1)-r2*3.*dr(1))
		  lrf(jg_cgp)%phi3(2 )=lrf(jg_cgp)%phi3(2 )+field2*(5.*dr(1)*dr(1)*dr(2)-r2*dr(2))
		  lrf(jg_cgp)%phi3(3 )=lrf(jg_cgp)%phi3(3 )+field2*(5.*dr(1)*dr(1)*dr(3)-r2*dr(3))
		  lrf(jg_cgp)%phi3(4 )=lrf(jg_cgp)%phi3(4 )+field2*(5.*dr(1)*dr(2)*dr(1)-r2*dr(2))
		  lrf(jg_cgp)%phi3(5 )=lrf(jg_cgp)%phi3(5 )+field2*(5.*dr(1)*dr(2)*dr(2)-r2*dr(1))
		  lrf(jg_cgp)%phi3(6 )=lrf(jg_cgp)%phi3(6 )+field2*(5.*dr(1)*dr(2)*dr(3))
		  lrf(jg_cgp)%phi3(7 )=lrf(jg_cgp)%phi3(7 )+field2*(5.*dr(1)*dr(3)*dr(1)-r2*dr(3))
		  lrf(jg_cgp)%phi3(8 )=lrf(jg_cgp)%phi3(8 )+field2*(5.*dr(1)*dr(3)*dr(2))
		  lrf(jg_cgp)%phi3(9 )=lrf(jg_cgp)%phi3(9 )+field2*(5.*dr(1)*dr(3)*dr(3)-r2*dr(1))
		  lrf(jg_cgp)%phi3(10)=lrf(jg_cgp)%phi3(10)+field2*(5.*dr(2)*dr(1)*dr(1)-r2*dr(2))
		  lrf(jg_cgp)%phi3(11)=lrf(jg_cgp)%phi3(11)+field2*(5.*dr(2)*dr(1)*dr(2)-r2*dr(1))
		  lrf(jg_cgp)%phi3(12)=lrf(jg_cgp)%phi3(12)+field2*(5.*dr(2)*dr(1)*dr(3))
		  lrf(jg_cgp)%phi3(13)=lrf(jg_cgp)%phi3(13)+field2*(5.*dr(2)*dr(2)*dr(1)-r2*dr(1))
		  lrf(jg_cgp)%phi3(14)=lrf(jg_cgp)%phi3(14)+field2*(5.*dr(2)*dr(2)*dr(2)-r2*3.*dr(2))
		  lrf(jg_cgp)%phi3(15)=lrf(jg_cgp)%phi3(15)+field2*(5.*dr(2)*dr(2)*dr(3)-r2*dr(3))
		  lrf(jg_cgp)%phi3(16)=lrf(jg_cgp)%phi3(16)+field2*(5.*dr(2)*dr(3)*dr(1))
		  lrf(jg_cgp)%phi3(17)=lrf(jg_cgp)%phi3(17)+field2*(5.*dr(2)*dr(3)*dr(2)-r2*dr(3))
		  lrf(jg_cgp)%phi3(18)=lrf(jg_cgp)%phi3(18)+field2*(5.*dr(2)*dr(3)*dr(3)-r2*dr(2))
		  lrf(jg_cgp)%phi3(19)=lrf(jg_cgp)%phi3(19)+field2*(5.*dr(3)*dr(1)*dr(1)-r2*dr(3))
		  lrf(jg_cgp)%phi3(20)=lrf(jg_cgp)%phi3(20)+field2*(5.*dr(3)*dr(1)*dr(2))
		  lrf(jg_cgp)%phi3(21)=lrf(jg_cgp)%phi3(21)+field2*(5.*dr(3)*dr(1)*dr(3)-r2*dr(1))
		  lrf(jg_cgp)%phi3(22)=lrf(jg_cgp)%phi3(22)+field2*(5.*dr(3)*dr(2)*dr(1))
		  lrf(jg_cgp)%phi3(23)=lrf(jg_cgp)%phi3(23)+field2*(5.*dr(3)*dr(2)*dr(2)-r2*dr(3))
		  lrf(jg_cgp)%phi3(24)=lrf(jg_cgp)%phi3(24)+field2*(5.*dr(3)*dr(2)*dr(3)-r2*dr(2))
		  lrf(jg_cgp)%phi3(25)=lrf(jg_cgp)%phi3(25)+field2*(5.*dr(3)*dr(3)*dr(1)-r2*dr(1))
		  lrf(jg_cgp)%phi3(26)=lrf(jg_cgp)%phi3(26)+field2*(5.*dr(3)*dr(3)*dr(2)-r2*dr(2))
		  lrf(jg_cgp)%phi3(27)=lrf(jg_cgp)%phi3(27)+field2*(5.*dr(3)*dr(3)*dr(3)-r2*3.*dr(3))
		 enddo ialoop2

		!solvent : solut	
		boxshiftx = x(j3+1) - lrf(ig)%cgp_cent(1)
		boxshifty = x(j3+2) - lrf(ig)%cgp_cent(2)
		boxshiftz = x(j3+3) - lrf(ig)%cgp_cent(3)

		boxshiftx = boxlength(1)*nint( boxshiftx*inv_boxl(1) )
		boxshifty = boxlength(2)*nint( boxshifty*inv_boxl(2) )
		boxshiftz = boxlength(3)*nint( boxshiftz*inv_boxl(3) )

jaloop2: do ja = 1, 3
		  j = nat_solute + 3*jg-3 + ja
		  j3 = j*3-3

				  ! calculate dr and (d)r2
		  dr(1) = x(j3+1) - lrf(ig)%cgp_cent(1) - boxshiftx
		  dr(2) = x(j3+2) - lrf(ig)%cgp_cent(2) - boxshifty
		  dr(3) = x(j3+3) - lrf(ig)%cgp_cent(3) - boxshiftz
		  
	      r2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)

				  ! calculate lrf for the water molecule
		  field0=crg(j)/(r2*sqrt(r2))
		  lrf(ig)%phi0=lrf(ig)%phi0+field0*r2
		  lrf(ig)%phi1(1)=lrf(ig)%phi1(1)-field0*dr(1)
		  lrf(ig)%phi1(2)=lrf(ig)%phi1(2)-field0*dr(2)
		  lrf(ig)%phi1(3)=lrf(ig)%phi1(3)-field0*dr(3)
		  field1=3.*field0/r2
		  lrf(ig)%phi2(1)=lrf(ig)%phi2(1)+field1*dr(1)*dr(1)-field0
		  lrf(ig)%phi2(2)=lrf(ig)%phi2(2)+field1*dr(1)*dr(2)
		  lrf(ig)%phi2(3)=lrf(ig)%phi2(3)+field1*dr(1)*dr(3)
		  lrf(ig)%phi2(4)=lrf(ig)%phi2(4)+field1*dr(2)*dr(1)
		  lrf(ig)%phi2(5)=lrf(ig)%phi2(5)+field1*dr(2)*dr(2)-field0
		  lrf(ig)%phi2(6)=lrf(ig)%phi2(6)+field1*dr(2)*dr(3)
		  lrf(ig)%phi2(7)=lrf(ig)%phi2(7)+field1*dr(3)*dr(1)
		  lrf(ig)%phi2(8)=lrf(ig)%phi2(8)+field1*dr(3)*dr(2)
		  lrf(ig)%phi2(9)=lrf(ig)%phi2(9)+field1*dr(3)*dr(3)-field0
		  field2=-field1/r2
		  lrf(ig)%phi3(1 )=lrf(ig)%phi3(1 )+field2*(5.*dr(1)*dr(1)*dr(1)-r2*3.*dr(1))
		  lrf(ig)%phi3(2 )=lrf(ig)%phi3(2 )+field2*(5.*dr(1)*dr(1)*dr(2)-r2*dr(2))
		  lrf(ig)%phi3(3 )=lrf(ig)%phi3(3 )+field2*(5.*dr(1)*dr(1)*dr(3)-r2*dr(3))
		  lrf(ig)%phi3(4 )=lrf(ig)%phi3(4 )+field2*(5.*dr(1)*dr(2)*dr(1)-r2*dr(2))
		  lrf(ig)%phi3(5 )=lrf(ig)%phi3(5 )+field2*(5.*dr(1)*dr(2)*dr(2)-r2*dr(1))
		  lrf(ig)%phi3(6 )=lrf(ig)%phi3(6 )+field2*(5.*dr(1)*dr(2)*dr(3))
		  lrf(ig)%phi3(7 )=lrf(ig)%phi3(7 )+field2*(5.*dr(1)*dr(3)*dr(1)-r2*dr(3))
		  lrf(ig)%phi3(8 )=lrf(ig)%phi3(8 )+field2*(5.*dr(1)*dr(3)*dr(2))
		  lrf(ig)%phi3(9 )=lrf(ig)%phi3(9 )+field2*(5.*dr(1)*dr(3)*dr(3)-r2*dr(1))
		  lrf(ig)%phi3(10)=lrf(ig)%phi3(10)+field2*(5.*dr(2)*dr(1)*dr(1)-r2*dr(2))
		  lrf(ig)%phi3(11)=lrf(ig)%phi3(11)+field2*(5.*dr(2)*dr(1)*dr(2)-r2*dr(1))
		  lrf(ig)%phi3(12)=lrf(ig)%phi3(12)+field2*(5.*dr(2)*dr(1)*dr(3))
		  lrf(ig)%phi3(13)=lrf(ig)%phi3(13)+field2*(5.*dr(2)*dr(2)*dr(1)-r2*dr(1))
		  lrf(ig)%phi3(14)=lrf(ig)%phi3(14)+field2*(5.*dr(2)*dr(2)*dr(2)-r2*3.*dr(2))
		  lrf(ig)%phi3(15)=lrf(ig)%phi3(15)+field2*(5.*dr(2)*dr(2)*dr(3)-r2*dr(3))
		  lrf(ig)%phi3(16)=lrf(ig)%phi3(16)+field2*(5.*dr(2)*dr(3)*dr(1))
		  lrf(ig)%phi3(17)=lrf(ig)%phi3(17)+field2*(5.*dr(2)*dr(3)*dr(2)-r2*dr(3))
		  lrf(ig)%phi3(18)=lrf(ig)%phi3(18)+field2*(5.*dr(2)*dr(3)*dr(3)-r2*dr(2))
		  lrf(ig)%phi3(19)=lrf(ig)%phi3(19)+field2*(5.*dr(3)*dr(1)*dr(1)-r2*dr(3))
		  lrf(ig)%phi3(20)=lrf(ig)%phi3(20)+field2*(5.*dr(3)*dr(1)*dr(2))
		  lrf(ig)%phi3(21)=lrf(ig)%phi3(21)+field2*(5.*dr(3)*dr(1)*dr(3)-r2*dr(1))
		  lrf(ig)%phi3(22)=lrf(ig)%phi3(22)+field2*(5.*dr(3)*dr(2)*dr(1))
		  lrf(ig)%phi3(23)=lrf(ig)%phi3(23)+field2*(5.*dr(3)*dr(2)*dr(2)-r2*dr(3))
		  lrf(ig)%phi3(24)=lrf(ig)%phi3(24)+field2*(5.*dr(3)*dr(2)*dr(3)-r2*dr(2))
		  lrf(ig)%phi3(25)=lrf(ig)%phi3(25)+field2*(5.*dr(3)*dr(3)*dr(1)-r2*dr(1))
		  lrf(ig)%phi3(26)=lrf(ig)%phi3(26)+field2*(5.*dr(3)*dr(3)*dr(2)-r2*dr(2))
		  lrf(ig)%phi3(27)=lrf(ig)%phi3(27)+field2*(5.*dr(3)*dr(3)*dr(3)-r2*3.*dr(3))
         enddo jaloop2
      end if
    end do jgloop
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
integer						::	iq, j, jq, is, i, k,l
real(8)                     :: el_scale
logical                     :: set

nbqq_pair(:) = 0

!list Q-Q
do iq = 1, nqat - 1
   do jq = iq + 1, nqat
       j = iqseq(jq)
       do is = 1, nstates
           if(qconn(is, j, iq) > 3) then
                   nbqq_pair(is) = nbqq_pair(is)+1
                   nbqq(nbqq_pair(is),is)%iq = iq
                   nbqq(nbqq_pair(is),is)%j = j
                   nbqq(nbqq_pair(is),is)%jq = jq

                   if(qconn(is, j, iq) == 4) then
                         nbqq(nbqq_pair(is),is)%LJcod = 3
                   elseif(.not. qvdw_flag) then
                         nbqq(nbqq_pair(is),is)%LJcod = ljcod(iac(iqseq(iq)),iac(j))
                   else
                         nbqq(nbqq_pair(is),is)%LJcod = 1
                         do i = 1, nqexpnb
                           if ((iq == iqexpnb(i) .and.  &
                             jq == jqexpnb(i))  .or. &
                             ( jq == iqexpnb(i) .and. &
                             iq == jqexpnb(i)))  then
                                  nbqq(nbqq_pair(is),is)%LJcod = 2
                                  exit
                           end if
                         end do
                   end if !if (qconn = 4)
				   if (nel_scale .eq. 0) then
                     nbqq(nbqq_pair(is),is)%el_scale = 1.0
				   else
  				     set=.false.
				     do i=1,nel_scale
                       k=qq_el_scale(i)%iqat
                       l=qq_el_scale(i)%jqat
                       if ((iq == k .and. jq == l) .or. &
                         (iq == l .and. jq == k)) then
                           nbqq(nbqq_pair(is),is)%el_scale = qq_el_scale(i)%el_scale(is)!masoud
						   set=.true.
						   exit
                       end if
                       if (.not. set) nbqq(nbqq_pair(is),is)%el_scale = 1.0
				     end do
				   end if
           end if !if (qconn > 3)
       end do
   end do
end do

!list Q-non-Q
do j = 1, nat_solute
        if(iqatom(j) > 0) cycle
        if(any(qconn(:,j,:) <= 3)) then 
                !bonded or angled to at least one Q-atom
                do iq = 1, nqat
                        do is = 1, nstates
                                if(qconn(is, j, iq) >= 4) then
                                        nbqq_pair(is) = nbqq_pair(is)+1
                                        nbqq(nbqq_pair(is),is)%iq = iq
                                        nbqq(nbqq_pair(is),is)%j = j
                                        nbqq(nbqq_pair(is),is)%jq = 0
										nbqq(nbqq_pair(is),is)%el_scale=1.0
                                        if(qconn(is, j, iq) == 4) then
                                                nbqq(nbqq_pair(is),is)%LJcod = 3
                                        elseif(qvdw_flag) then
                                                nbqq(nbqq_pair(is),is)%LJcod = 1
                                        else
                                                nbqq(nbqq_pair(is),is)%LJcod = ljcod(iac(iqseq(iq)),iac(j))
                                        end if
                                end if
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
real(8)						:: rcut2,r2

!******PWadded variables 2001-10-01

real(8)						:: dx, dy, dz

!	This routine counts non-bonded atom pairs involving
!	*one* Q-atom and *one* non-Q-atom, where the latter is *not connected*
!	(meaning not bonded or angled) to any Q-atom.
!
!	( i , j ) pairs correspond to
!	( iq, j ) with first index being a Q-atom with the *Q-atom numbering*, &
!	          and the second index is the non-Q-atom.


nqp = 0
rcut2 = Rcq*Rcq



if(nqat==0) return

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
        dx = x(i3+1) - x(3*qswitch-2)
        dy = x(i3+2) - x(3*qswitch-1)
        dz = x(i3+3) - x(3*qswitch)
        dx = dx - boxlength(1)*nint( dx*inv_boxl(1) )
        dy = dy - boxlength(2)*nint( dy*inv_boxl(2) )
        dz = dz - boxlength(3)*nint( dz*inv_boxl(3) )
        r2 = dx**2 + dy**2 + dz**2
end if

! skip if outside cutoff
if ( r2 .gt. rcut2 ) cycle igloop

ialoop: do ia = cgp(ig)%first, cgp(ig)%last
  i = cgpatom(ia)

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
real(8)						:: rcut2,r2

!******PWadded variables

real(8)						:: dx, dy, dz

!	This routine counts water molecules that interact with q-atoms
nqw = 0
rcut2 = Rcq*Rcq


if(nqat==0) return

! --- solvent - Q-atoms

iwloop: do ig = 1, nwat
nqwmol(ig) = 0
ia = nat_solute + 3*ig-2
if(.not. use_PBC .and. excl(ia)) cycle iwloop ! skip excluded waters
i3 = 3*ia-3

!******PWadded if-statement 2001-10-01
if( .not. use_PBC ) then
        r2 = ( x(i3+1) - xpcent(1) )**2 &
                +( x(i3+2) - xpcent(2) )**2 &
                +( x(i3+3) - xpcent(3) )**2
else
        dx = x(i3+1) - x(3*qswitch-2)
        dy = x(i3+2) - x(3*qswitch-1)
        dz = x(i3+3) - x(3*qswitch)
        dx = dx - boxlength(1)*nint( dx*inv_boxl(1) )
        dy = dy - boxlength(2)*nint( dy*inv_boxl(2) )
        dz = dz - boxlength(3)*nint( dz*inv_boxl(3) )
        r2 = dx**2 + dy**2 + dz**2
end if

! skip if outside cutoff
if ( r2 <= rcut2 ) then
        nqw = nqw + 1
        nqwmol(ig) = 3*nqat
end if
end do iwloop

end subroutine nbqw_count

!-----------------------------------------------------------------------

subroutine nbqplis2
! local variables
integer						:: ig,ia,i,j,iq,i3,nl,inside
real(8)						:: rcut2,r2
integer						:: xspec
logical, save					:: list_done


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
real(8)                                         :: start_loop_time
start_loop_time = rtime()
#endif




if(list_done .and. rcq > rexcl_o) return

nbqp_pair = 0
rcut2 = Rcq*Rcq


if(nqat == 0) return


! --- solute - Q-atoms

igloop: do ig = calculation_assignment%qp%start, calculation_assignment%qp%end
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

  ! check if already on qq list
  if(any(qconn(:,i,:) <= 3)) cycle ialoop

        ! if out of space then make more space
        if (nbqp_pair .ge. calculation_assignment%qp%max-nqat) call reallocate_nonbondlist_qp
qaloop: do iq = 1, nqat

        !check special exclusions
        !is already done in make_qconn

        ! store the pair
        nbqp_pair = nbqp_pair + 1
        nbqp(nbqp_pair)%i = iq
        nbqp(nbqp_pair)%j = i
        nbqp(nbqp_pair)%LJcod = ljcod(iac(i),iac(iqseq(iq)))

        ! adjust LJcod for neighbors
        ! check only first state, should be same for all in this list
        if(qconn(1, i, iq) == 4) nbqp(nbqp_pair)%LJcod = 3

        !store special q-LJ code = 1,1,3 for normal code = 1,2,3 resp.
        nbqp(nbqp_pair)%qLJcod = nbqp(nbqp_pair)%LJcod
        if(nbqp(nbqp_pair)%qLJcod == 2) nbqp(nbqp_pair)%qLJcod = 1

  end do qaloop
end do ialoop
end do igloop

list_done = .true. !save this value
#if defined (PROFILING)
profile(5)%time = profile(5)%time + rtime() - start_loop_time
#endif

end subroutine nbqplis2



!----------------------------------------------------------------------------

subroutine nbqplis2_box
!  ! local variables
  integer						:: ig,ia,i,j,iq,i3,nl,inside,ig_atom
  real(8)						:: rcut2,r2
  integer						:: xspec
  real(8)						:: dx, dy, dz
 
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
#if defined (PROFILING)
real(8)                                         :: start_loop_time
start_loop_time = rtime()
#endif


  nbqp_pair = 0
  nbqp_cgp_pair = 0
  rcut2 = Rcq*Rcq
	
  if(nqat == 0) return

 ! --- solute - Q-atoms

igloop: do ig = calculation_assignment%qp%start, calculation_assignment%qp%end
    ! for every assigned charge group:

	! check cutoff
	inside = 0
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

		if(nbqp_cgp_pair .eq. size(nbqp_cgp, 1) ) call reallocate_nbqp_cgp
		nbqp_cgp_pair = nbqp_cgp_pair + 1
		nbqp_cgp(nbqp_cgp_pair)%i = i !leave %j empty, equals qswitch
	  end if

	  ig_atom = ig_atom + 1 !ia = ia + 1
	end do
	if (inside .eq. 0) cycle igloop

ialoop: do ia = cgp(ig)%first, cgp(ig)%last
	  i = cgpatom(ia)

	  ! check if already on qq list
	  if(any(qconn(:,i,:) <= 3)) cycle ialoop

		! if out of space then make more space
		if (nbqp_pair .ge. calculation_assignment%qp%max-nqat) call reallocate_nonbondlist_qp
qaloop: do iq = 1, nqat

		!check special exclusions
		!is already done in make_qconn

		! store the pair
		nbqp_pair = nbqp_pair + 1
		nbqp(nbqp_pair)%i = iq
		nbqp(nbqp_pair)%j = i
		nbqp(nbqp_pair)%LJcod = ljcod(iac(i),iac(iqseq(iq)))
		nbqp(nbqp_pair)%cgp_pair = nbqp_cgp_pair

		! adjust LJcod for neighbors
		! check only first state, should be same for all in this list
		if(qconn(1, i, iq) == 4) nbqp(nbqp_pair)%LJcod = 3
		
		!store special q-LJ code = 1,1,3 for normal code = 1,2,3 resp.
		nbqp(nbqp_pair)%qLJcod = nbqp(nbqp_pair)%LJcod
		if(nbqp(nbqp_pair)%qLJcod == 2) nbqp(nbqp_pair)%qLJcod = 1

	  end do qaloop
	end do ialoop
  end do igloop
#if defined (PROFILING)
profile(5)%time = profile(5)%time + rtime() - start_loop_time
#endif

end subroutine nbqplis2_box
!-----------------------------------------------------------------------

subroutine nbqplist
! local variables
integer						:: ig,ia,i,j,iq,i3,nl
real(8)						:: rcut2,r2
integer						::	xspec
logical, save				::	list_done = .false.

#if defined (PROFILING)
real(8)						:: start_loop_time
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
if(list_done .and. rcq > rexcl_o) return

nbqp_pair = 0
rcut2 = Rcq*Rcq

if(nqat==0) return

! --- solute - Q-atoms

igloop: do ig = calculation_assignment%qp%start, calculation_assignment%qp%end

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
        if (nbqp_pair .ge. calculation_assignment%qp%max-nqat) call reallocate_nonbondlist_qp
qaloop: do iq = 1, nqat

        !check special exclusions
        !is already done in make_qconn

        ! add the pair
        nbqp_pair = nbqp_pair + 1
        nbqp(nbqp_pair)%i = iq
        nbqp(nbqp_pair)%j = i
        nbqp(nbqp_pair)%LJcod = ljcod(iac(i),iac(iqseq(iq)))

        ! adjust LJcod for neighbors
        ! check only first state, should be same for all in this list
        if(qconn(1, i, iq) == 4) nbqp(nbqp_pair)%LJcod = 3

        !store special q-LJ code = 1,1,3 for normal code = 1,2,3 resp.
        nbqp(nbqp_pair)%qLJcod = nbqp(nbqp_pair)%LJcod
        if(nbqp(nbqp_pair)%qLJcod == 2) nbqp(nbqp_pair)%qLJcod = 1
  end do qaloop
!TMP        nbqp_per_cgp(ig) = nbqp_per_cgp(ig) + nqat
end do ialoop
end do igloop

list_done = .true. !save this value

#if defined (PROFILING)
profile(5)%time = profile(5)%time + rtime() - start_loop_time
#endif

end subroutine nbqplist

!-----------------------------------------------------------------------

!******PWadded 2001-10-18
subroutine nbqplist_box
  ! local variables
  integer						:: ig,ia,i,j,iq,i3,nl,ig_sw
  real(8)						:: rcut2,r2
  integer						::	xspec
  real(8)						:: dx, dy, dz
  integer						:: ga, gb

#if defined (PROFILING)
real(8)                                         :: start_loop_time
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


  nbqp_pair = 0
  nbqp_cgp_pair = 0
  rcut2 = Rcq*Rcq

  if(nqat==0) return

  ! --- solute - Q-atoms
igloop: do ig = calculation_assignment%qp%start, calculation_assignment%qp%end

	ig_sw = cgp(ig)%iswitch
	
	i3 = 3*ig_sw-3
	dx = x(i3+1) - x(3*qswitch-2)
	dy = x(i3+2) - x(3*qswitch-1)
	dz = x(i3+3) - x(3*qswitch)
	dx = dx - boxlength(1)*nint( dx*inv_boxl(1) )
	dy = dy - boxlength(2)*nint( dy*inv_boxl(2) )
	dz = dz - boxlength(3)*nint( dz*inv_boxl(3) )
	r2 = dx**2 + dy**2 + dz**2

	! skip if outside cutoff
	if ( r2 .gt. rcut2 ) cycle igloop

	if( nbqp_cgp_pair .eq. size(nbqp_cgp, 1) ) call reallocate_nbqp_cgp
	
	nbqp_cgp_pair = nbqp_cgp_pair + 1
	nbqp_cgp(nbqp_cgp_pair)%i = ig_sw
	

ialoop: do ia = cgp(ig)%first, cgp(ig)%last

	  i = cgpatom(ia)

	  ! check if already on qq list
	  if(any(qconn(:,i,:) <= 3)) cycle ialoop

          ! if out of space then make more space
          if (nbqp_pair .ge. calculation_assignment%qp%max-nqat) call reallocate_nonbondlist_qp
qaloop: do iq = 1, nqat

		
		!check special exclusions
		!is already done in make_qconn

		! add the pair
		nbqp_pair = nbqp_pair + 1
		nbqp(nbqp_pair)%i = iq
		nbqp(nbqp_pair)%j = i
		nbqp(nbqp_pair)%LJcod = ljcod(iac(i),iac(iqseq(iq)))
		nbqp(nbqp_pair)%cgp_pair = nbqp_cgp_pair
		
		! adjust LJcod for neighbors
		! check only first state, should be same for all in this list
		if(qconn(1, i, iq) == 4) nbqp(nbqp_pair)%LJcod = 3
		
		!store special q-LJ code = 1,1,3 for normal code = 1,2,3 resp.
		nbqp(nbqp_pair)%qLJcod = nbqp(nbqp_pair)%LJcod
		if(nbqp(nbqp_pair)%qLJcod == 2) nbqp(nbqp_pair)%qLJcod = 1

	  end do qaloop
	end do ialoop
   end do igloop

#if defined (PROFILING)
profile(5)%time = profile(5)%time + rtime() - start_loop_time
#endif

end subroutine nbqplist_box
!------------------------------------------------------------------------------------

subroutine nbqwlist

! local variables
integer						:: ig,ia,i,i3
real(8)						:: rcut2,r2


!	This routine makes a list of water molecules within rcq from xpcent,
! i.e. the q-atom - water non-bond lists which implicitly includes all
! q-atoms with all atoms of the listed water
! waters may not have bonded interactions with q-atoms !
#if defined (PROFILING)
real(8)                                         :: start_loop_time
start_loop_time = rtime()
#endif



!We don't have to remake the list if q-atoms interact with all waters
!(rcq > rexcl_o) and list already made (nbwq_pair >0)
if(rcq > rexcl_o .and. nbqw_pair > 0) return


nbqw_pair = 0
rcut2 = Rcq*Rcq


if(nqat==0) return


iwloop: do ig = calculation_assignment%qw%start, calculation_assignment%qw%end
ia = nat_solute + 3*ig-2
if( excl(ia) ) cycle iwloop ! skip excluded waters

i3 = 3*ia-3


r2 = ( x(i3+1) - xpcent(1) )**2 &
        +( x(i3+2) - xpcent(2) )**2 &
        +( x(i3+3) - xpcent(3) )**2


! store if inside cutoff
if ( r2 <= rcut2 ) then
        nbqw_pair = nbqw_pair + 1
        nbqw(nbqw_pair) = nat_solute + 3*ig-2
end if
end do iwloop
#if defined (PROFILING)
profile(6)%time = profile(6)%time + rtime() - start_loop_time
#endif

end subroutine nbqwlist

!-----------------------------------------------------------------------

!******PWadded 2001-10-18
subroutine nbqwlist_box

  ! local variables
  integer						:: ig,ia,i,i3
  real(8)						:: rcut2,r2
  real(8)						:: dx, dy, dz

  !	This routine makes a list of water molecules within rcq from xswitch,
  ! i.e. the q-atom - water non-bond lists which implicitly includes all
  ! q-atoms with all atoms of the listed water
  ! waters may not have bonded interactions with q-atoms !
#if defined (PROFILING)
real(8)                                         :: start_loop_time
start_loop_time = rtime()
#endif


	nbqw_pair = 0
	rcut2 = Rcq*Rcq

	if(nqat==0) return

iwloop: do ig = calculation_assignment%qw%start, calculation_assignment%qw%end
	ia = nat_solute + 3*ig-2
	i3 = 3*ia-3

	dx = x(i3+1) - x(3*qswitch-2)
	dy = x(i3+2) - x(3*qswitch-1)
	dz = x(i3+3) - x(3*qswitch)
	dx = dx - boxlength(1)*nint( dx*inv_boxl(1) )
	dy = dy - boxlength(2)*nint( dy*inv_boxl(2) )
	dz = dz - boxlength(3)*nint( dz*inv_boxl(3) )
	r2 = dx**2 + dy**2 + dz**2

	! store if inside cutoff
	if ( r2 <= rcut2 ) then
		nbqw_pair = nbqw_pair + 1
		nbqw(nbqw_pair) = nat_solute + 3*ig-2
	end if
  end do iwloop
#if defined (PROFILING)
profile(6)%time = profile(6)%time + rtime() - start_loop_time
#endif

end subroutine nbqwlist_box
!-------------------------------------------------------------------------------------
!******PWchanged 2001-10-01
subroutine nbww_count(nww, nwwmol)
! arguments
integer						:: nww
integer						:: nwwmol(:)

! local variables
integer						:: iw,jw,ia,ja,i3,j3
real(8)						:: rcut2,r2

!******PWadded variables		

real(8)						:: dx, dy, dz

! This routine counts non-bonded solvent-solvent atom pairs.

nww = 0
rcut2 = Rcww*Rcww

iwloop: do iw = 1, nwat
nwwmol(iw) = 0

ia = nat_solute + 3*iw-2
if(.not. use_PBC .and. excl(ia)) cycle iwloop ! skip excluded waters

i3 = 3*ia-3

jwloop: do jw = 1, nwat
  ja = nat_solute + 3*jw-2
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
        nwwmol(iw) = nwwmol(iw) + 9
  end if

end do jwloop
end do iwloop
end subroutine nbww_count

!-----------------------------------------------------------------------
subroutine nbwwlist
! local variables
integer						:: iw,jw,ia,ja,i3,j3
real(8)						:: rcut2,r2

! This routine makes a list of non-bonded solvent-solvent atom pairs

! uses the global variables:
!  nbww_pair, Rcww, nat_solute, excl, nwat, x, nbww, calculation_assignment%ww%max,nbww_true_pair
#if defined (PROFILING)
real(8)                                         :: start_loop_time
start_loop_time = rtime()
#endif


nbww_pair = 0
nbww_true_pair=0
rcut2 = Rcww*Rcww

iwloop: do iw = calculation_assignment%ww%start, calculation_assignment%ww%end

ia = nat_solute + 3*iw-2

        if (.not.excl(ia)) then
        i3 = 3*ia-3

jwloop:		do jw=1, nwat

                  ja = nat_solute + 3*jw-2
                  if (excl(ja)) cycle jwloop ! skip excluded waters
                  j3 = 3*ja-3


                 ! count each w-w pair once only
                if ( ((iw .gt. jw) .and. (mod(iw+jw,2) .eq. 0)) .or. &
                     ((iw .lt. jw) .and. (mod(iw+jw,2) .eq. 1)) .or. (iw .eq. jw))  cycle jwloop

                r2    = ( x(i3+1) -x(j3+1) )**2 &
                      + ( x(i3+2) -x(j3+2) )**2 &
                      + ( x(i3+3) -x(j3+3) )**2

                if ( r2 .le. rcut2 ) then
                                ! inside cutoff: add the pair
                                nbww_true_pair = nbww_true_pair + 9
                                nbww_pair = nbww_pair + 1
                                nbww(nbww_pair) = ja
                end if
                        ! if out of space then make more space	
                        if (nbww_pair .ge. calculation_assignment%ww%max) then
                                call reallocate_nonbondlist_ww
                        end if
                end do jwloop

        end if !if ia not excluded
        ! now mark the end of the list of molecules interacting with molecule iw
        ! by means of a zero element in the list
        nbww_pair = nbww_pair + 1
        nbww(nbww_pair) = 0	
end do iwloop
#if defined (PROFILING)
profile(2)%time = profile(2)%time + rtime() - start_loop_time
#endif

end subroutine nbwwlist
!------------PWadded 2001-10-18------------------------------------------
subroutine nbwwlist_box

! local variables
integer						:: iw,jw,ia,ja,i3,j3
real(8)						:: rcut2,r2
real(8)						:: dx, dy, dz

! for periodic boundary conditions
! This routine makes a list of non-bonded solvent-solvent atom pairs
! uses the global variables:
!  nbww_pair, Rcww, nat_solute, excl, nwat, x, nbww, calculation_assignment%ww%max,nbww_true_pair
#if defined (PROFILING)
real(8)                                         :: start_loop_time
start_loop_time = rtime()
#endif


nbww_true_pair = 0
nbww_pair = 0
rcut2 = Rcww*Rcww

iwloop: do iw = calculation_assignment%ww%start, calculation_assignment%ww%end

        ia = nat_solute + 3*iw-2
        i3 = 3*ia-3

jwloop:		do jw=1, nwat

                        ja = nat_solute + 3*jw-2
                        j3 = 3*ja-3

                        ! count each w-w pair once only
                        if ( ((iw .gt. jw) .and. (mod(iw+jw,2) .eq. 0)) .or. &
                           ((iw .lt. jw) .and. (mod(iw+jw,2) .eq. 1)) .or. &
                           (iw .eq. jw)) &
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
                                nbww_true_pair = nbww_true_pair + 9
                                nbww_pair = nbww_pair + 1
                                nbww(nbww_pair) = ja
                        end if

                        ! if out of space then make more space	
                        if (nbww_pair .ge. calculation_assignment%ww%max) then
                                call reallocate_nonbondlist_ww
                        end if

                end do jwloop

        ! now mark the end of the list of molecules interacting with molecule iw
        ! by means of a zero element in the list
        nbww_pair = nbww_pair + 1
        nbww(nbww_pair) = 0	
end do iwloop
#if defined (PROFILING)
profile(2)%time = profile(2)%time + rtime() - start_loop_time
#endif

end subroutine nbwwlist_box

!---------------------------------------------------------------------------

subroutine nbwwlist_lrf
! local variables
integer						:: i,j,ig,jg,iw,jw,ia,ja,i3,j3,is,is3
real(8)						:: rcut2,r2,field0,field1,field2
real(8)						:: dr(3)
real(8)						::	RcLRF2

#if defined (PROFILING)
real(8)						:: start_loop_time
start_loop_time = rtime()
#endif

!	This routine makes a list of non-bonded solvent-solvent atom pairs.

! uses the global variables:
!  nbww_pair, Rcww, nwat, nat_solute, excl, x, nbww, ncgp, lrf, crg, calculation_assignment%ww%max,nbww_true_pair






nbww_true_pair=0
nbww_pair = 0
rcut2 = Rcww*Rcww
RcLRF2 = RcLRF*RcLRF

iwloop: do iw = calculation_assignment%ww%start, calculation_assignment%ww%end
is  = nat_solute + 3*iw-2
if(.not. excl(is)) then
is3 = 3*is-3

jwloop: do jw = 1, nwat
  ja = nat_solute + 3*jw-2
  if(excl(ja)) cycle jwloop ! skip excluded waters
  j3 = 3*ja-3

  ! count each w-w pair once only
  if ( ((iw .gt. jw) .and. (mod(iw+jw,2) .eq. 0)) .or. &
           ((iw .lt. jw) .and. (mod(iw+jw,2) .eq. 1)) .or. &
           (iw .eq. jw) ) &
           cycle jwloop

r2 = ( x(is3+1) -x(j3+1) )**2 &
        +( x(is3+2) -x(j3+2) )**2 &
        +( x(is3+3) -x(j3+3) )**2

  if ( r2 .le. rcut2 ) then
        ! inside cutoff: add the pair
        nbww_pair = nbww_pair + 1
        nbww_true_pair = nbww_true_pair + 9  !To get explicit no. of interactions
        nbww(nbww_pair) = ja 
  elseif(r2 <= RcLRF2) then
        ! outside ww-cutoff but inside LRF cut-off: use LRF
        do ia=1,3
          i = nat_solute+iw*3-3+ia
          i3 = i*3-3
          !jg = ncgp + jw
          jg = iwhich_cgp(ja)

          dr(1) = x(i3+1) - lrf(jg)%cgp_cent(1)
          dr(2) = x(i3+2) - lrf(jg)%cgp_cent(2)
          dr(3) = x(i3+3) - lrf(jg)%cgp_cent(3)
          r2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)



          field0=crg(i)/(r2*sqrt(r2))
          lrf(jg)%phi0=lrf(jg)%phi0+field0*r2
          lrf(jg)%phi1(1)=lrf(jg)%phi1(1)-field0*dr(1)
          lrf(jg)%phi1(2)=lrf(jg)%phi1(2)-field0*dr(2)
          lrf(jg)%phi1(3)=lrf(jg)%phi1(3)-field0*dr(3)
          field1=3.*field0/r2
          lrf(jg)%phi2(1)=lrf(jg)%phi2(1)+field1*dr(1)*dr(1)-field0
          lrf(jg)%phi2(2)=lrf(jg)%phi2(2)+field1*dr(1)*dr(2)
          lrf(jg)%phi2(3)=lrf(jg)%phi2(3)+field1*dr(1)*dr(3)
          lrf(jg)%phi2(4)=lrf(jg)%phi2(4)+field1*dr(2)*dr(1)
          lrf(jg)%phi2(5)=lrf(jg)%phi2(5)+field1*dr(2)*dr(2)-field0
          lrf(jg)%phi2(6)=lrf(jg)%phi2(6)+field1*dr(2)*dr(3)
          lrf(jg)%phi2(7)=lrf(jg)%phi2(7)+field1*dr(3)*dr(1)
          lrf(jg)%phi2(8)=lrf(jg)%phi2(8)+field1*dr(3)*dr(2)
          lrf(jg)%phi2(9)=lrf(jg)%phi2(9)+field1*dr(3)*dr(3)-field0
          field2=-field1/r2
          lrf(jg)%phi3(1 )=lrf(jg)%phi3(1 )+field2*(5.*dr(1)*dr(1)*dr(1)-r2*3.*dr(1))
          lrf(jg)%phi3(2 )=lrf(jg)%phi3(2 )+field2*(5.*dr(1)*dr(1)*dr(2)-r2*dr(2))
          lrf(jg)%phi3(3 )=lrf(jg)%phi3(3 )+field2*(5.*dr(1)*dr(1)*dr(3)-r2*dr(3))
          lrf(jg)%phi3(4 )=lrf(jg)%phi3(4 )+field2*(5.*dr(1)*dr(2)*dr(1)-r2*dr(2))
          lrf(jg)%phi3(5 )=lrf(jg)%phi3(5 )+field2*(5.*dr(1)*dr(2)*dr(2)-r2*dr(1))
          lrf(jg)%phi3(6 )=lrf(jg)%phi3(6 )+field2*(5.*dr(1)*dr(2)*dr(3))
          lrf(jg)%phi3(7 )=lrf(jg)%phi3(7 )+field2*(5.*dr(1)*dr(3)*dr(1)-r2*dr(3))
          lrf(jg)%phi3(8 )=lrf(jg)%phi3(8 )+field2*(5.*dr(1)*dr(3)*dr(2))
          lrf(jg)%phi3(9 )=lrf(jg)%phi3(9 )+field2*(5.*dr(1)*dr(3)*dr(3)-r2*dr(1))
          lrf(jg)%phi3(10)=lrf(jg)%phi3(10)+field2*(5.*dr(2)*dr(1)*dr(1)-r2*dr(2))
          lrf(jg)%phi3(11)=lrf(jg)%phi3(11)+field2*(5.*dr(2)*dr(1)*dr(2)-r2*dr(1))
          lrf(jg)%phi3(12)=lrf(jg)%phi3(12)+field2*(5.*dr(2)*dr(1)*dr(3))
          lrf(jg)%phi3(13)=lrf(jg)%phi3(13)+field2*(5.*dr(2)*dr(2)*dr(1)-r2*dr(1))
          lrf(jg)%phi3(14)=lrf(jg)%phi3(14)+field2*(5.*dr(2)*dr(2)*dr(2)-r2*3.*dr(2))
          lrf(jg)%phi3(15)=lrf(jg)%phi3(15)+field2*(5.*dr(2)*dr(2)*dr(3)-r2*dr(3))
          lrf(jg)%phi3(16)=lrf(jg)%phi3(16)+field2*(5.*dr(2)*dr(3)*dr(1))
          lrf(jg)%phi3(17)=lrf(jg)%phi3(17)+field2*(5.*dr(2)*dr(3)*dr(2)-r2*dr(3))
          lrf(jg)%phi3(18)=lrf(jg)%phi3(18)+field2*(5.*dr(2)*dr(3)*dr(3)-r2*dr(2))
          lrf(jg)%phi3(19)=lrf(jg)%phi3(19)+field2*(5.*dr(3)*dr(1)*dr(1)-r2*dr(3))
          lrf(jg)%phi3(20)=lrf(jg)%phi3(20)+field2*(5.*dr(3)*dr(1)*dr(2))
          lrf(jg)%phi3(21)=lrf(jg)%phi3(21)+field2*(5.*dr(3)*dr(1)*dr(3)-r2*dr(1))
          lrf(jg)%phi3(22)=lrf(jg)%phi3(22)+field2*(5.*dr(3)*dr(2)*dr(1))
          lrf(jg)%phi3(23)=lrf(jg)%phi3(23)+field2*(5.*dr(3)*dr(2)*dr(2)-r2*dr(3))
          lrf(jg)%phi3(24)=lrf(jg)%phi3(24)+field2*(5.*dr(3)*dr(2)*dr(3)-r2*dr(2))
          lrf(jg)%phi3(25)=lrf(jg)%phi3(25)+field2*(5.*dr(3)*dr(3)*dr(1)-r2*dr(1))
          lrf(jg)%phi3(26)=lrf(jg)%phi3(26)+field2*(5.*dr(3)*dr(3)*dr(2)-r2*dr(2))
          lrf(jg)%phi3(27)=lrf(jg)%phi3(27)+field2*(5.*dr(3)*dr(3)*dr(3)-r2*3.*dr(3))
        enddo

        do ja = 1, 3
          j = nat_solute + 3*jw-3 + ja
          j3 = j*3-3
          !ig = ncgp + iw
          ig = iwhich_cgp(is)


          dr(1) = x(j3+1) - lrf(ig)%cgp_cent(1)
          dr(2) = x(j3+2) - lrf(ig)%cgp_cent(2)
          dr(3) = x(j3+3) - lrf(ig)%cgp_cent(3)
          r2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)

          field0=crg(j)/(r2*sqrt(r2))
          lrf(ig)%phi0=lrf(ig)%phi0+field0*r2
          lrf(ig)%phi1(1)=lrf(ig)%phi1(1)-field0*dr(1)
          lrf(ig)%phi1(2)=lrf(ig)%phi1(2)-field0*dr(2)
          lrf(ig)%phi1(3)=lrf(ig)%phi1(3)-field0*dr(3)
          field1=3.*field0/r2
          lrf(ig)%phi2(1)=lrf(ig)%phi2(1)+field1*dr(1)*dr(1)-field0
          lrf(ig)%phi2(2)=lrf(ig)%phi2(2)+field1*dr(1)*dr(2)
          lrf(ig)%phi2(3)=lrf(ig)%phi2(3)+field1*dr(1)*dr(3)
          lrf(ig)%phi2(4)=lrf(ig)%phi2(4)+field1*dr(2)*dr(1)
          lrf(ig)%phi2(5)=lrf(ig)%phi2(5)+field1*dr(2)*dr(2)-field0
          lrf(ig)%phi2(6)=lrf(ig)%phi2(6)+field1*dr(2)*dr(3)
          lrf(ig)%phi2(7)=lrf(ig)%phi2(7)+field1*dr(3)*dr(1)
          lrf(ig)%phi2(8)=lrf(ig)%phi2(8)+field1*dr(3)*dr(2)
          lrf(ig)%phi2(9)=lrf(ig)%phi2(9)+field1*dr(3)*dr(3)-field0
          field2=-field1/r2
          lrf(ig)%phi3(1 )=lrf(ig)%phi3(1 )+field2*(5.*dr(1)*dr(1)*dr(1)-r2*3.*dr(1))
          lrf(ig)%phi3(2 )=lrf(ig)%phi3(2 )+field2*(5.*dr(1)*dr(1)*dr(2)-r2*dr(2))
          lrf(ig)%phi3(3 )=lrf(ig)%phi3(3 )+field2*(5.*dr(1)*dr(1)*dr(3)-r2*dr(3))
          lrf(ig)%phi3(4 )=lrf(ig)%phi3(4 )+field2*(5.*dr(1)*dr(2)*dr(1)-r2*dr(2))
          lrf(ig)%phi3(5 )=lrf(ig)%phi3(5 )+field2*(5.*dr(1)*dr(2)*dr(2)-r2*dr(1))
          lrf(ig)%phi3(6 )=lrf(ig)%phi3(6 )+field2*(5.*dr(1)*dr(2)*dr(3))
          lrf(ig)%phi3(7 )=lrf(ig)%phi3(7 )+field2*(5.*dr(1)*dr(3)*dr(1)-r2*dr(3))
          lrf(ig)%phi3(8 )=lrf(ig)%phi3(8 )+field2*(5.*dr(1)*dr(3)*dr(2))
          lrf(ig)%phi3(9 )=lrf(ig)%phi3(9 )+field2*(5.*dr(1)*dr(3)*dr(3)-r2*dr(1))
          lrf(ig)%phi3(10)=lrf(ig)%phi3(10)+field2*(5.*dr(2)*dr(1)*dr(1)-r2*dr(2))
          lrf(ig)%phi3(11)=lrf(ig)%phi3(11)+field2*(5.*dr(2)*dr(1)*dr(2)-r2*dr(1))
          lrf(ig)%phi3(12)=lrf(ig)%phi3(12)+field2*(5.*dr(2)*dr(1)*dr(3))
          lrf(ig)%phi3(13)=lrf(ig)%phi3(13)+field2*(5.*dr(2)*dr(2)*dr(1)-r2*dr(1))
          lrf(ig)%phi3(14)=lrf(ig)%phi3(14)+field2*(5.*dr(2)*dr(2)*dr(2)-r2*3.*dr(2))
          lrf(ig)%phi3(15)=lrf(ig)%phi3(15)+field2*(5.*dr(2)*dr(2)*dr(3)-r2*dr(3))
          lrf(ig)%phi3(16)=lrf(ig)%phi3(16)+field2*(5.*dr(2)*dr(3)*dr(1))
          lrf(ig)%phi3(17)=lrf(ig)%phi3(17)+field2*(5.*dr(2)*dr(3)*dr(2)-r2*dr(3))
          lrf(ig)%phi3(18)=lrf(ig)%phi3(18)+field2*(5.*dr(2)*dr(3)*dr(3)-r2*dr(2))
          lrf(ig)%phi3(19)=lrf(ig)%phi3(19)+field2*(5.*dr(3)*dr(1)*dr(1)-r2*dr(3))
          lrf(ig)%phi3(20)=lrf(ig)%phi3(20)+field2*(5.*dr(3)*dr(1)*dr(2))
          lrf(ig)%phi3(21)=lrf(ig)%phi3(21)+field2*(5.*dr(3)*dr(1)*dr(3)-r2*dr(1))
          lrf(ig)%phi3(22)=lrf(ig)%phi3(22)+field2*(5.*dr(3)*dr(2)*dr(1))
          lrf(ig)%phi3(23)=lrf(ig)%phi3(23)+field2*(5.*dr(3)*dr(2)*dr(2)-r2*dr(3))
          lrf(ig)%phi3(24)=lrf(ig)%phi3(24)+field2*(5.*dr(3)*dr(2)*dr(3)-r2*dr(2))
          lrf(ig)%phi3(25)=lrf(ig)%phi3(25)+field2*(5.*dr(3)*dr(3)*dr(1)-r2*dr(1))
          lrf(ig)%phi3(26)=lrf(ig)%phi3(26)+field2*(5.*dr(3)*dr(3)*dr(2)-r2*dr(2))
          lrf(ig)%phi3(27)=lrf(ig)%phi3(27)+field2*(5.*dr(3)*dr(3)*dr(3)-r2*3.*dr(3))
        enddo
  end if

  ! if out of space then make more space
  if (nbww_pair .ge. calculation_assignment%ww%max) call reallocate_nonbondlist_ww

end do jwloop
end if !if ia not excluded
 !now mark the end of the list of molecules interacting with molecule iw
 !by means of a zero element in the list
 nbww_pair = nbww_pair + 1
 nbww(nbww_pair) = 0	

end do iwloop

#if defined (PROFILING)
profile(2)%time = profile(2)%time + rtime() - start_loop_time
#endif

end subroutine nbwwlist_lrf
!--------------LRF version of PW PBC-----------------------------------
subroutine nbwwlist_box_lrf
! local variables
integer						:: i,j,ig,jg,iw,jw,ia,ja,i3,j3,is,is3
real(8)						:: rcut2,r2,field0,field1,field2
real(8)						:: dr(3)
real(8)						::	RcLRF2
real(8)						::  dx, dy, dz
real(8)                     :: boxshiftx, boxshifty, boxshiftz  

#if defined (PROFILING)
real(8)                                         :: start_loop_time
start_loop_time = rtime()
#endif


!	This routine makes a list of non-bonded solvent-solvent atom pairs.

! uses the global variables:
!  nbww_pair, Rcww, nwat, nat_solute, excl, x, nbww, ncgp, lrf, crg, calculation_assignment%ww%max,nbww_true_pair

nbww_true_pair=0
nbww_pair = 0
rcut2 = Rcww*Rcww
RcLRF2 = RcLRF*RcLRF

iwloop: do iw = calculation_assignment%ww%start, calculation_assignment%ww%end
is  = nat_solute + 3*iw-2
is3 = 3*is-3

jwloop: do jw = 1, nwat
  ja = nat_solute + 3*jw-2
  j3 = 3*ja-3

  ! count each w-w pair once only
  if ( ((iw .gt. jw) .and. (mod(iw+jw,2) .eq. 0)) .or. &
           ((iw .lt. jw) .and. (mod(iw+jw,2) .eq. 1)) .or. &
           (iw .eq. jw) ) &
           cycle jwloop

dx = x(is3+1) -x(j3+1)
dy = x(is3+2) -x(j3+2)
dz = x(is3+3) -x(j3+3)

dx = dx - boxlength(1)*nint( dx*inv_boxl(1) )
dy = dy - boxlength(2)*nint( dy*inv_boxl(2) )
dz = dz - boxlength(3)*nint( dz*inv_boxl(3) )
 
r2 = dx*dx + dy*dy + dz*dz

  if ( r2 .le. rcut2 ) then
        ! inside cutoff: add the pair
        nbww_pair = nbww_pair + 1
        nbww_true_pair = nbww_true_pair + 9  !To get explicit no. of interactions
        nbww(nbww_pair) = ja 
  elseif(r2 <= RcLRF2) then
        ! outside ww-cutoff but inside LRF cut-off: use LRF
        
		!iw interaction     
	    jg = iwhich_cgp(ja)
		
		boxshiftx = x(is3+1) -lrf(jg)%cgp_cent(1)
		boxshifty = x(is3+2) -lrf(jg)%cgp_cent(2)
		boxshiftz = x(is3+3) -lrf(jg)%cgp_cent(3)

		boxshiftx = boxlength(1)*nint( boxshiftx*inv_boxl(1) )
		boxshifty = boxlength(2)*nint( boxshifty*inv_boxl(2) )
		boxshiftz = boxlength(3)*nint( boxshiftz*inv_boxl(3) )

		do ia=1,3
          i = nat_solute+iw*3-3+ia
          i3 = i*3-3
          !jg = ncgp + jw
          

          dr(1) = x(i3+1) - lrf(jg)%cgp_cent(1) - boxshiftx
          dr(2) = x(i3+2) - lrf(jg)%cgp_cent(2) - boxshifty
          dr(3) = x(i3+3) - lrf(jg)%cgp_cent(3) - boxshiftz
          
		  r2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)

          field0=crg(i)/(r2*sqrt(r2))
          lrf(jg)%phi0=lrf(jg)%phi0+field0*r2
          lrf(jg)%phi1(1)=lrf(jg)%phi1(1)-field0*dr(1)
          lrf(jg)%phi1(2)=lrf(jg)%phi1(2)-field0*dr(2)
          lrf(jg)%phi1(3)=lrf(jg)%phi1(3)-field0*dr(3)
          field1=3.*field0/r2
          lrf(jg)%phi2(1)=lrf(jg)%phi2(1)+field1*dr(1)*dr(1)-field0
          lrf(jg)%phi2(2)=lrf(jg)%phi2(2)+field1*dr(1)*dr(2)
          lrf(jg)%phi2(3)=lrf(jg)%phi2(3)+field1*dr(1)*dr(3)
          lrf(jg)%phi2(4)=lrf(jg)%phi2(4)+field1*dr(2)*dr(1)
          lrf(jg)%phi2(5)=lrf(jg)%phi2(5)+field1*dr(2)*dr(2)-field0
          lrf(jg)%phi2(6)=lrf(jg)%phi2(6)+field1*dr(2)*dr(3)
          lrf(jg)%phi2(7)=lrf(jg)%phi2(7)+field1*dr(3)*dr(1)
          lrf(jg)%phi2(8)=lrf(jg)%phi2(8)+field1*dr(3)*dr(2)
          lrf(jg)%phi2(9)=lrf(jg)%phi2(9)+field1*dr(3)*dr(3)-field0
          field2=-field1/r2
          lrf(jg)%phi3(1 )=lrf(jg)%phi3(1 )+field2*(5.*dr(1)*dr(1)*dr(1)-r2*3.*dr(1))
          lrf(jg)%phi3(2 )=lrf(jg)%phi3(2 )+field2*(5.*dr(1)*dr(1)*dr(2)-r2*dr(2))
          lrf(jg)%phi3(3 )=lrf(jg)%phi3(3 )+field2*(5.*dr(1)*dr(1)*dr(3)-r2*dr(3))


  lrf(jg)%phi3(4 )=lrf(jg)%phi3(4 )+field2*(5.*dr(1)*dr(2)*dr(1)-r2*dr(2))
          lrf(jg)%phi3(5 )=lrf(jg)%phi3(5 )+field2*(5.*dr(1)*dr(2)*dr(2)-r2*dr(1))
          lrf(jg)%phi3(6 )=lrf(jg)%phi3(6 )+field2*(5.*dr(1)*dr(2)*dr(3))
          lrf(jg)%phi3(7 )=lrf(jg)%phi3(7 )+field2*(5.*dr(1)*dr(3)*dr(1)-r2*dr(3))
          lrf(jg)%phi3(8 )=lrf(jg)%phi3(8 )+field2*(5.*dr(1)*dr(3)*dr(2))
          lrf(jg)%phi3(9 )=lrf(jg)%phi3(9 )+field2*(5.*dr(1)*dr(3)*dr(3)-r2*dr(1))
          lrf(jg)%phi3(10)=lrf(jg)%phi3(10)+field2*(5.*dr(2)*dr(1)*dr(1)-r2*dr(2))
          lrf(jg)%phi3(11)=lrf(jg)%phi3(11)+field2*(5.*dr(2)*dr(1)*dr(2)-r2*dr(1))
          lrf(jg)%phi3(12)=lrf(jg)%phi3(12)+field2*(5.*dr(2)*dr(1)*dr(3))
          lrf(jg)%phi3(13)=lrf(jg)%phi3(13)+field2*(5.*dr(2)*dr(2)*dr(1)-r2*dr(1))
          lrf(jg)%phi3(14)=lrf(jg)%phi3(14)+field2*(5.*dr(2)*dr(2)*dr(2)-r2*3.*dr(2))
          lrf(jg)%phi3(15)=lrf(jg)%phi3(15)+field2*(5.*dr(2)*dr(2)*dr(3)-r2*dr(3))
          lrf(jg)%phi3(16)=lrf(jg)%phi3(16)+field2*(5.*dr(2)*dr(3)*dr(1))
          lrf(jg)%phi3(17)=lrf(jg)%phi3(17)+field2*(5.*dr(2)*dr(3)*dr(2)-r2*dr(3))
          lrf(jg)%phi3(18)=lrf(jg)%phi3(18)+field2*(5.*dr(2)*dr(3)*dr(3)-r2*dr(2))
          lrf(jg)%phi3(19)=lrf(jg)%phi3(19)+field2*(5.*dr(3)*dr(1)*dr(1)-r2*dr(3))
          lrf(jg)%phi3(20)=lrf(jg)%phi3(20)+field2*(5.*dr(3)*dr(1)*dr(2))
          lrf(jg)%phi3(21)=lrf(jg)%phi3(21)+field2*(5.*dr(3)*dr(1)*dr(3)-r2*dr(1))
          lrf(jg)%phi3(22)=lrf(jg)%phi3(22)+field2*(5.*dr(3)*dr(2)*dr(1))
          lrf(jg)%phi3(23)=lrf(jg)%phi3(23)+field2*(5.*dr(3)*dr(2)*dr(2)-r2*dr(3))
          lrf(jg)%phi3(24)=lrf(jg)%phi3(24)+field2*(5.*dr(3)*dr(2)*dr(3)-r2*dr(2))
          lrf(jg)%phi3(25)=lrf(jg)%phi3(25)+field2*(5.*dr(3)*dr(3)*dr(1)-r2*dr(1))
          lrf(jg)%phi3(26)=lrf(jg)%phi3(26)+field2*(5.*dr(3)*dr(3)*dr(2)-r2*dr(2))
          lrf(jg)%phi3(27)=lrf(jg)%phi3(27)+field2*(5.*dr(3)*dr(3)*dr(3)-r2*3.*dr(3))
        enddo
		
		!jw interaction
		ig = iwhich_cgp(is)
		
		boxshiftx = x(j3+1) -lrf(ig)%cgp_cent(1)
		boxshifty = x(j3+2) -lrf(ig)%cgp_cent(2)
		boxshiftz = x(j3+3) -lrf(ig)%cgp_cent(3)

		boxshiftx = boxlength(1)*nint( boxshiftx*inv_boxl(1) )
		boxshifty = boxlength(2)*nint( boxshifty*inv_boxl(2) )
		boxshiftz = boxlength(3)*nint( boxshiftz*inv_boxl(3) )

        do ja = 1, 3
          j = nat_solute + 3*jw-3 + ja
          j3 = j*3-3
          !ig = ncgp + iw
          
          dr(1) = x(j3+1) - lrf(ig)%cgp_cent(1) - boxshiftx
          dr(2) = x(j3+2) - lrf(ig)%cgp_cent(2) - boxshifty
          dr(3) = x(j3+3) - lrf(ig)%cgp_cent(3) - boxshiftz

          r2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)

          field0=crg(j)/(r2*sqrt(r2))
          lrf(ig)%phi0=lrf(ig)%phi0+field0*r2
          lrf(ig)%phi1(1)=lrf(ig)%phi1(1)-field0*dr(1)
          lrf(ig)%phi1(2)=lrf(ig)%phi1(2)-field0*dr(2)
          lrf(ig)%phi1(3)=lrf(ig)%phi1(3)-field0*dr(3)
          field1=3.*field0/r2
          lrf(ig)%phi2(1)=lrf(ig)%phi2(1)+field1*dr(1)*dr(1)-field0
          lrf(ig)%phi2(2)=lrf(ig)%phi2(2)+field1*dr(1)*dr(2)
          lrf(ig)%phi2(3)=lrf(ig)%phi2(3)+field1*dr(1)*dr(3)
          lrf(ig)%phi2(4)=lrf(ig)%phi2(4)+field1*dr(2)*dr(1)
          lrf(ig)%phi2(5)=lrf(ig)%phi2(5)+field1*dr(2)*dr(2)-field0
          lrf(ig)%phi2(6)=lrf(ig)%phi2(6)+field1*dr(2)*dr(3)
          lrf(ig)%phi2(7)=lrf(ig)%phi2(7)+field1*dr(3)*dr(1)
          lrf(ig)%phi2(8)=lrf(ig)%phi2(8)+field1*dr(3)*dr(2)
          lrf(ig)%phi2(9)=lrf(ig)%phi2(9)+field1*dr(3)*dr(3)-field0
          field2=-field1/r2
          lrf(ig)%phi3(1 )=lrf(ig)%phi3(1 )+field2*(5.*dr(1)*dr(1)*dr(1)-r2*3.*dr(1))
          lrf(ig)%phi3(2 )=lrf(ig)%phi3(2 )+field2*(5.*dr(1)*dr(1)*dr(2)-r2*dr(2))
          lrf(ig)%phi3(3 )=lrf(ig)%phi3(3 )+field2*(5.*dr(1)*dr(1)*dr(3)-r2*dr(3))
          lrf(ig)%phi3(4 )=lrf(ig)%phi3(4 )+field2*(5.*dr(1)*dr(2)*dr(1)-r2*dr(2))
          lrf(ig)%phi3(5 )=lrf(ig)%phi3(5 )+field2*(5.*dr(1)*dr(2)*dr(2)-r2*dr(1))
          lrf(ig)%phi3(6 )=lrf(ig)%phi3(6 )+field2*(5.*dr(1)*dr(2)*dr(3))
          lrf(ig)%phi3(7 )=lrf(ig)%phi3(7 )+field2*(5.*dr(1)*dr(3)*dr(1)-r2*dr(3))
          lrf(ig)%phi3(8 )=lrf(ig)%phi3(8 )+field2*(5.*dr(1)*dr(3)*dr(2))
          lrf(ig)%phi3(9 )=lrf(ig)%phi3(9 )+field2*(5.*dr(1)*dr(3)*dr(3)-r2*dr(1))
          lrf(ig)%phi3(10)=lrf(ig)%phi3(10)+field2*(5.*dr(2)*dr(1)*dr(1)-r2*dr(2))
          lrf(ig)%phi3(11)=lrf(ig)%phi3(11)+field2*(5.*dr(2)*dr(1)*dr(2)-r2*dr(1))
          lrf(ig)%phi3(12)=lrf(ig)%phi3(12)+field2*(5.*dr(2)*dr(1)*dr(3))
          lrf(ig)%phi3(13)=lrf(ig)%phi3(13)+field2*(5.*dr(2)*dr(2)*dr(1)-r2*dr(1))
          lrf(ig)%phi3(14)=lrf(ig)%phi3(14)+field2*(5.*dr(2)*dr(2)*dr(2)-r2*3.*dr(2))
          lrf(ig)%phi3(15)=lrf(ig)%phi3(15)+field2*(5.*dr(2)*dr(2)*dr(3)-r2*dr(3))
          lrf(ig)%phi3(16)=lrf(ig)%phi3(16)+field2*(5.*dr(2)*dr(3)*dr(1))
          lrf(ig)%phi3(17)=lrf(ig)%phi3(17)+field2*(5.*dr(2)*dr(3)*dr(2)-r2*dr(3))
          lrf(ig)%phi3(18)=lrf(ig)%phi3(18)+field2*(5.*dr(2)*dr(3)*dr(3)-r2*dr(2))
          lrf(ig)%phi3(19)=lrf(ig)%phi3(19)+field2*(5.*dr(3)*dr(1)*dr(1)-r2*dr(3))
          lrf(ig)%phi3(20)=lrf(ig)%phi3(20)+field2*(5.*dr(3)*dr(1)*dr(2))
          lrf(ig)%phi3(21)=lrf(ig)%phi3(21)+field2*(5.*dr(3)*dr(1)*dr(3)-r2*dr(1))
          lrf(ig)%phi3(22)=lrf(ig)%phi3(22)+field2*(5.*dr(3)*dr(2)*dr(1))
          lrf(ig)%phi3(23)=lrf(ig)%phi3(23)+field2*(5.*dr(3)*dr(2)*dr(2)-r2*dr(3))
          lrf(ig)%phi3(24)=lrf(ig)%phi3(24)+field2*(5.*dr(3)*dr(2)*dr(3)-r2*dr(2))
          lrf(ig)%phi3(25)=lrf(ig)%phi3(25)+field2*(5.*dr(3)*dr(3)*dr(1)-r2*dr(1))
          lrf(ig)%phi3(26)=lrf(ig)%phi3(26)+field2*(5.*dr(3)*dr(3)*dr(2)-r2*dr(2))
          lrf(ig)%phi3(27)=lrf(ig)%phi3(27)+field2*(5.*dr(3)*dr(3)*dr(3)-r2*3.*dr(3))
        enddo
  end if

  ! if out of space then make more space
  if (nbww_pair .ge. calculation_assignment%ww%max) call reallocate_nonbondlist_ww

end do jwloop
 !now mark the end of the list of molecules interacting with molecule iw
 !by means of a zero element in the list
 nbww_pair = nbww_pair + 1
 nbww(nbww_pair) = 0	

end do iwloop
#if defined (PROFILING)
profile(2)%time = profile(2)%time + rtime() - start_loop_time
#endif

end subroutine nbwwlist_box_lrf
!---------------------------------------------------------------------------
subroutine nbmonitorlist
!set LJ code for the atom pairs in the selected atom groups to be monitored
! local variables
integer         :: i,j,ig,jg,ia,ja,i3,j3,nl,istate,LJ_code,maxingroups,par, atomnri
integer         :: grpi,grpj,atomi,atomj,qq_pair,aLJ,bLJ


if (monitor_group_pairs == 0) return

!check the size of the largest group
maxingroups=maxval(monitor_atom_group(:)%n)
allocate (special_LJcod(nstates,maxingroups,maxingroups,monitor_group_pairs))

do par=1,monitor_group_pairs
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
              do istate=1,nstates
                      !starting guess = use LJ_code matrix for topology atom types
                      LJ_code = ljcod(iac(atomi),iac(atomj))
                      if(qvdw_flag .and. LJ_code == 2) then
                         if((iqatom(atomi) /= 0 .and. iqatom(atomj) == 0) .or.  &
                            (iqatom(atomj)/= 0 .and. iqatom(atomi) == 0)) then
                            !can't use code 2 between q and non-q when Q-atom
                            !types are used. Q-Type 2 params are for exp. 
                            !repulsion, not LJ!
                            LJ_code = 1
                          end if
                      end if
                      !Are atoms of pair in 1-4 position?
                      if(iqatom(atomi) == 0 .and. iqatom(atomj) == 0) then !neither atom is q_atom
                        if (abs(atomj-atomi) .le. max_nbr_range ) then
                          if (atomi .gt. atomj ) then
                            if ( list14(atomi-atomj, atomi) ) LJ_code = 3 !3 means 1-4
                          else
                            if ( list14(atomj-atomi, atomj) ) LJ_code = 3
                          end if
                        else
                          do nl = 1, n14long
                            if ((list14long(1, nl) .eq. atomi &
                              .and. list14long(2, nl) .eq. atomj) .or. &
                               (list14long(1, nl) .eq. atomj &
                              .and. list14long(2, nl) .eq. atomi)) then
                                  LJ_code = 3
                            endif
                          end do
                        endif	!kolla 1-4 interaktioner
                      else	!at least one is Q-atom
                            !check Q-Q pairlist to find LJ-code
                        do qq_pair = 1, nbqq_pair(istate)
                          atomnri=iqseq(nbqq(qq_pair,istate)%iq)         !Find atom number from qatom number
                          if((atomnri == atomi .and. &     
                              nbqq(qq_pair,istate)%j == atomj ) .or. &
                             (nbqq(qq_pair,istate)%j == atomi .and. &
                              atomnri == atomj)) then   
                              LJ_code = nbqq(qq_pair,istate)%LJcod
                                 exit	
                          end if
                        end do
                           !if not found here then the first guess should be used
                      end if
                                special_LJcod(istate,i,j,par)=LJ_code
              end do !nstates
      end do   ! monitor_atom_group j
    end do  !monitor_atom_group i
end do  !par
end subroutine nbmonitorlist

!---------------------------------------------------------------------------------------------------------		   

subroutine nonbond_monitor
!monitor nonbonded energies between selected groups of atoms

real(8)  :: dx1,dx2,dx3,r,r2,r6
integer	 :: i,j,istate,LJ_code,par
integer	 :: grpi,grpj,atomi,atomj,qatomi,qatomj,qq_pair, iaci,iacj
real(8)  :: aLJi,bLJi,aLJj,bLJj,qi,qj, Vel,Vvdw,Vwel,Vwvdw,Vwsum
real(8)  :: r6_hc    !  softcore variables
integer  :: sc_1,sc_2     !  softcore variables, sc_1 is the first index in sc_lookup (the qatom)
logical  :: do_sc    !  softcore variables,   do_sc is a boolean to determine if softcore should be done
					! do_sc is true when atom i or j is a qatom  (and qvdw is true)


do par=1,monitor_group_pairs

        grpi=monitor_group_pair(par)%i      
        grpj=monitor_group_pair(par)%j 
        monitor_group_pair(par)%Vel(:)=0
        monitor_group_pair(par)%Vlj(:)=0
        monitor_group_pair(par)%Vwel = 0
        monitor_group_pair(par)%Vwlj = 0
        monitor_group_pair(par)%Vwsum= 0     

        do i=1,monitor_atom_group(grpi)%n
            
			atomi=monitor_atom_group(grpi)%atom(i)
            qatomi = iqatom(atomi)
			qi   = crg(atomi)
			iaci = iac(atomi)
            
			    do j=1,monitor_atom_group(grpj)%n
            
			            atomj=monitor_atom_group(grpj)%atom(j)  
                        qatomj = iqatom(atomj)
                        iacj = iac(atomj)
                        qj   = crg(atomj)
                        dx1  = x(3*atomi-2)-x(3*atomj-2)      ! calculate the distance
                        dx2  = x(3*atomi-1)-x(3*atomj-1)
                        dx3  = x(3*atomi)-x(3*atomj)
						
						
						if (use_PBC) then
							dx1 = dx1 - boxlength(1)*nint( dx1*inv_boxl(1) )
							dx2 = dx2 - boxlength(2)*nint( dx2*inv_boxl(2) )
							dx3 = dx3 - boxlength(3)*nint( dx3*inv_boxl(3) )
						end if

                        r2   = dx1*dx1 + dx2*dx2 + dx3*dx3
                        r    = SQRT(1/r2) 
                        r6   = r2*r2*r2
						r6_hc= r6   !needed for softcore
  						r6   = 1._8/r6
            
			            do istate=1,nstates			
								do_sc = .false.  !default is no softcore

								LJ_code  = special_LJcod(istate,i,j,par)
                                aLJi=iaclib(iaci)%avdw(LJ_code)
                                bLJi=iaclib(iaci)%bvdw(LJ_code)
                                aLJj=iaclib(iacj)%avdw(LJ_code)
                                bLJj=iaclib(iacj)%bvdw(LJ_code)   
                                if (qatomi /= 0) qi = qcrg(qatomi,istate) 
                                if (qatomj /= 0) qj = qcrg(qatomj,istate) 
                                if (qvdw_flag) then
                                        if (qatomi/=0) then
                                                iaci = qiac(qatomi,istate)
                                                aLJi = qavdw(iaci,LJ_code)
                                                bLJi = qbvdw(iaci,LJ_code)


												do_sc = .true.   ! atom i is a q-atom, softcore on
												sc_1 = qatomi    ! the first index in sc_lookup should be a qatom number

										else
												
												sc_2 = iaci      ! atom i was not a qatom, put atom code i in the second sc_lookup index

                                        endif	


                                        if (qatomj/=0) then					
                                                iacj  = qiac(qatomj,istate)			
                                                aLJj  = qavdw(iacj,LJ_code)			
                                                bLJj  = qbvdw(iacj,LJ_code)

												if (do_sc) then   ! do_sc is true if atom i is a qatom
													sc_2 = qatomj + natyps  ! qatomi is sc_1
												else
													do_sc = .true.  ! atom i was not a qatom but j is, softcore on
													sc_1 = qatomj	! qatom j should be the first index
												end if

										else
											sc_2 = iacj    ! atom j is not a qatom, should be index 2
										
												
                                        endif
										
										if (do_sc) then  ! calculate softcore r6
											r6 = r6_hc + sc_lookup (sc_1,sc_2,istate)
											r6 = 1._8/r6
										end if	 
                                endif										
                                Vel  = qi*qj*r	
                                if(ivdw_rule==1) then !geometric comb. rule
                                        Vvdw = aLJi*aLJj*r6*r6-bLJi*bLJj*r6
                                else !arithmetic
                                        Vvdw = bLJi * bLJj * (aLJi+aLJj)**6 * r6 * &
                                                ((aLJi+aLJj)**6 * r6 - 2.0)
                                endif
                                !add up for this pair of atom groups
                                monitor_group_pair(par)%Vel(istate)= monitor_group_pair(par)%Vel(istate)+Vel
                                monitor_group_pair(par)%Vlj(istate)= monitor_group_pair(par)%Vlj(istate)+Vvdw
                        end do	!nstates
                end do   ! monitor_atom_group j
        end do  !monitor_atom_group i
        !calc lambda-weighted sum
        monitor_group_pair(par)%Vwel=dot_product(monitor_group_pair(par)%Vel(1:nstates),EQ(1:nstates)%lambda)
        monitor_group_pair(par)%Vwlj=dot_product(monitor_group_pair(par)%Vlj(1:nstates),EQ(1:nstates)%lambda)
        monitor_group_pair(par)%Vwsum= monitor_group_pair(par)%Vwlj+monitor_group_pair(par)%Vwel
end do !par

end subroutine nonbond_monitor

!-------------------------------------------------------------------------

subroutine nonbon2_pp
! local variables
integer						:: ip
real(8)						:: aLJa,bLJa,dx1a,dx2a,dx3a,r2a,ra,r6a
real(8)						:: aLJb,bLJb,dx1b,dx2b,dx3b,r2b,rb,r6b
real(8)						:: Vela,V_aa,V_ba,dva
real(8)						:: Velb,V_ab,V_bb,dvb
type(NB_TYPE), pointer		:: pa
type(NB_TYPE), pointer		:: pb

! global variables used:
!  iaclib, x, crg, el14_scale, d, E

do ip = 1, nbpp_pair - 1, 2
! for every second pair:

! init pointers
pa => nbpp(ip)
pb => nbpp(ip+1)

! calculate aLJ and bLJ
aLJa  = iaclib(iac(pa%i))%avdw(pa%LJcod)+iaclib(iac(pa%j))%avdw(pa%LJcod)
aLJb  = iaclib(iac(pb%i))%avdw(pb%LJcod)+iaclib(iac(pb%j))%avdw(pb%LJcod)
bLJa  = iaclib(iac(pa%i))%bvdw(pa%LJcod)*iaclib(iac(pa%j))%bvdw(pa%LJcod)
bLJb  = iaclib(iac(pb%i))%bvdw(pb%LJcod)*iaclib(iac(pb%j))%bvdw(pb%LJcod)
aLJa  = aLJa*aLJa
aLJb  = aLJb*aLJb
aLJa  = aLJa*aLJa*aLJa
aLJb  = aLJb*aLJb*aLJb

! calculate dx, r and r2
dx1a  = x(pa%j*3-2) - x(pa%i*3-2)
dx1b  = x(pb%j*3-2) - x(pb%i*3-2)
dx2a  = x(pa%j*3-1) - x(pa%i*3-1)
dx2b  = x(pb%j*3-1) - x(pb%i*3-1)
dx3a  = x(pa%j*3-0) - x(pa%i*3-0)
dx3b  = x(pb%j*3-0) - x(pb%i*3-0)

r2a   = 1./(dx1a*dx1a + dx2a*dx2a + dx3a*dx3a)
r2b   = 1./(dx1b*dx1b + dx2b*dx2b + dx3b*dx3b)
ra = sqrt(r2a)
rb = sqrt(r2b) 

r6a   = r2a*r2a*r2a
r6b   = r2b*r2b*r2b

! calculate Vel and dv
Vela  = crg(pa%i)*crg(pa%j)*ra  
Velb  = crg(pb%i)*crg(pb%j)*rb  
if ( pa%LJcod .eq. 3 ) then
  Vela = Vela*el14_scale
end if
if ( pb%LJcod .eq. 3 ) then
  Velb = Velb*el14_scale
end if
V_aa  = bLJa*aLJa*aLJa*r6a*r6a 
V_ab  = bLJb*aLJb*aLJb*r6b*r6b 
V_ba  = 2.0*bLJa*aLJa*r6a
V_bb  = 2.0*bLJb*aLJb*r6b
dva   = r2a*( -Vela -12.*V_aa +6.*V_ba )
dvb   = r2b*( -Velb -12.*V_ab +6.*V_bb )

! update d
d(pa%i*3-2) = d(pa%i*3-2) - dva*dx1a
d(pb%i*3-2) = d(pb%i*3-2) - dvb*dx1b
d(pa%i*3-1) = d(pa%i*3-1) - dva*dx2a
d(pb%i*3-1) = d(pb%i*3-1) - dvb*dx2b
d(pa%i*3-0) = d(pa%i*3-0) - dva*dx3a
d(pb%i*3-0) = d(pb%i*3-0) - dvb*dx3b
d(pa%j*3-2) = d(pa%j*3-2) + dva*dx1a
d(pb%j*3-2) = d(pb%j*3-2) + dvb*dx1b
d(pa%j*3-1) = d(pa%j*3-1) + dva*dx2a
d(pb%j*3-1) = d(pb%j*3-1) + dvb*dx2b
d(pa%j*3-0) = d(pa%j*3-0) + dva*dx3a
d(pb%j*3-0) = d(pb%j*3-0) + dvb*dx3b

! update energies
E%pp%el  = E%pp%el + Vela + Velb
E%pp%vdw = E%pp%vdw + V_aa - V_ba + V_ab - V_bb
end do

if (ip .eq. nbpp_pair) then
! the last pair:

pa => nbpp(ip)
aLJa  = iaclib(iac(pa%i))%avdw(pa%LJcod)+iaclib(iac(pa%j))%avdw(pa%LJcod)
bLJa  = iaclib(iac(pa%i))%bvdw(pa%LJcod)*iaclib(iac(pa%j))%bvdw(pa%LJcod)
aLJa  = aLJa*aLJa
aLJa  = aLJa*aLJa*aLJa

dx1a  = x(pa%j*3-2) - x(pa%i*3-2)
dx2a  = x(pa%j*3-1) - x(pa%i*3-1)
dx3a  = x(pa%j*3-0) - x(pa%i*3-0)

r2a   = 1./(dx1a*dx1a + dx2a*dx2a + dx3a*dx3a)
ra = sqrt(r2a)

r6a   = r2a*r2a*r2a

Vela  = crg(pa%i)*crg(pa%j)*ra
if ( pa%LJcod .eq. 3 ) Vela = Vela*el14_scale

V_aa  = bLJa*aLJa*aLJa*r6a*r6a 
V_ba  = 2.0*bLJa*aLJa*r6a
dva   = r2a*( -Vela -12.*V_aa +6.*V_ba )

d(pa%i*3-2) = d(pa%i*3-2) - dva*dx1a
d(pa%i*3-1) = d(pa%i*3-1) - dva*dx2a
d(pa%i*3-0) = d(pa%i*3-0) - dva*dx3a
d(pa%j*3-2) = d(pa%j*3-2) + dva*dx1a
d(pa%j*3-1) = d(pa%j*3-1) + dva*dx2a
d(pa%j*3-0) = d(pa%j*3-0) + dva*dx3a

E%pp%el  = E%pp%el + Vela 
E%pp%vdw = E%pp%vdw + V_aa - V_ba 
end if

end subroutine nonbon2_pp

!------------------------------------------------------------------------
subroutine nonbon2_pp_box
  ! local variables
  integer						:: ip, ga, gb, group
  real(8)						:: aLJa,bLJa,dx1a,dx2a,dx3a,r2a,ra,r6a
  real(8)						:: aLJb,bLJb,dx1b,dx2b,dx3b,r2b,rb,r6b
  real(8)						:: Vela,V_aa,V_ba,dva
  real(8)						:: Velb,V_ab,V_bb,dvb
  type(NB_TYPE), pointer		:: pa
  type(NB_TYPE), pointer		:: pb

 ! global variables used:
  !  iaclib, x, crg, el14_scale, d, E


	do group = 1, nbpp_cgp_pair
		ga = nbpp_cgp(group)%i !atom index for the two switching atoms
		gb = nbpp_cgp(group)%j

		!the distance between the two switching atoms
		dx1a = x(3*gb-2) - x(3*ga-2)
		dx2a = x(3*gb-1) - x(3*ga-1)
		dx3a = x(3*gb  ) - x(3*ga  )

		nbpp_cgp(group)%x = boxlength(1)*nint( dx1a*inv_boxl(1) )
		nbpp_cgp(group)%y = boxlength(2)*nint( dx2a*inv_boxl(2) )	
		nbpp_cgp(group)%z = boxlength(3)*nint( dx3a*inv_boxl(3) )

	end do

  do ip = 1, nbpp_pair - 1, 2
    ! for every second pair:

	! init pointers
	pa => nbpp(ip)
	pb => nbpp(ip+1)
	ga = pa%cgp_pair
	gb = pb%cgp_pair

	! calculate aLJ and bLJ
    aLJa  = iaclib(iac(pa%i))%avdw(pa%LJcod)+iaclib(iac(pa%j))%avdw(pa%LJcod)
    aLJb  = iaclib(iac(pb%i))%avdw(pb%LJcod)+iaclib(iac(pb%j))%avdw(pb%LJcod)
    bLJa  = iaclib(iac(pa%i))%bvdw(pa%LJcod)*iaclib(iac(pa%j))%bvdw(pa%LJcod)
    bLJb  = iaclib(iac(pb%i))%bvdw(pb%LJcod)*iaclib(iac(pb%j))%bvdw(pb%LJcod)
    aLJa  = aLJa*aLJa
    aLJb  = aLJb*aLJb
    aLJa  = aLJa*aLJa*aLJa
    aLJb  = aLJb*aLJb*aLJb

	! calculate dx, r and r2
    dx1a  = x(pa%j*3-2) - x(pa%i*3-2)
    dx1b  = x(pb%j*3-2) - x(pb%i*3-2)
    dx2a  = x(pa%j*3-1) - x(pa%i*3-1)
    dx2b  = x(pb%j*3-1) - x(pb%i*3-1)
    dx3a  = x(pa%j*3-0) - x(pa%i*3-0)
    dx3b  = x(pb%j*3-0) - x(pb%i*3-0)
	dx1a = dx1a - nbpp_cgp(ga)%x         
	dx1b = dx1b - nbpp_cgp(gb)%x 
	dx2a = dx2a - nbpp_cgp(ga)%y    
	dx2b = dx2b - nbpp_cgp(gb)%y                
	dx3a = dx3a - nbpp_cgp(ga)%z  
	dx3b = dx3b - nbpp_cgp(gb)%z   
	


    r2a   = 1./(dx1a*dx1a + dx2a*dx2a + dx3a*dx3a)
    r2b   = 1./(dx1b*dx1b + dx2b*dx2b + dx3b*dx3b)
	ra = sqrt(r2a)
	rb = sqrt(r2b) 

    r6a   = r2a*r2a*r2a
    r6b   = r2b*r2b*r2b

	! calculate Vel and dv
    Vela  = crg(pa%i)*crg(pa%j)*ra  
    Velb  = crg(pb%i)*crg(pb%j)*rb  
    if ( pa%LJcod .eq. 3 ) then
	  Vela = Vela*el14_scale
	end if
    if ( pb%LJcod .eq. 3 ) then
	  Velb = Velb*el14_scale
	end if
    V_aa  = bLJa*aLJa*aLJa*r6a*r6a 
    V_ab  = bLJb*aLJb*aLJb*r6b*r6b 
    V_ba  = 2.0*bLJa*aLJa*r6a
    V_bb  = 2.0*bLJb*aLJb*r6b
    dva   = r2a*( -Vela -12.*V_aa +6.*V_ba )
    dvb   = r2b*( -Velb -12.*V_ab +6.*V_bb )

	! update d
    d(pa%i*3-2) = d(pa%i*3-2) - dva*dx1a
    d(pb%i*3-2) = d(pb%i*3-2) - dvb*dx1b
    d(pa%i*3-1) = d(pa%i*3-1) - dva*dx2a
    d(pb%i*3-1) = d(pb%i*3-1) - dvb*dx2b
    d(pa%i*3-0) = d(pa%i*3-0) - dva*dx3a
    d(pb%i*3-0) = d(pb%i*3-0) - dvb*dx3b
    d(pa%j*3-2) = d(pa%j*3-2) + dva*dx1a
    d(pb%j*3-2) = d(pb%j*3-2) + dvb*dx1b
    d(pa%j*3-1) = d(pa%j*3-1) + dva*dx2a
    d(pb%j*3-1) = d(pb%j*3-1) + dvb*dx2b
    d(pa%j*3-0) = d(pa%j*3-0) + dva*dx3a
    d(pb%j*3-0) = d(pb%j*3-0) + dvb*dx3b

	! update energies
    E%pp%el  = E%pp%el + Vela + Velb
    E%pp%vdw = E%pp%vdw + V_aa - V_ba + V_ab - V_bb
  end do

  if (ip .eq. nbpp_pair) then
    ! the last pair:

	pa => nbpp(ip)
	ga = pa%cgp_pair

    aLJa  = iaclib(iac(pa%i))%avdw(pa%LJcod)+iaclib(iac(pa%j))%avdw(pa%LJcod)
    bLJa  = iaclib(iac(pa%i))%bvdw(pa%LJcod)*iaclib(iac(pa%j))%bvdw(pa%LJcod)
    aLJa  = aLJa*aLJa
    aLJa  = aLJa*aLJa*aLJa

    dx1a  = x(pa%j*3-2) - x(pa%i*3-2)
    dx2a  = x(pa%j*3-1) - x(pa%i*3-1)
    dx3a  = x(pa%j*3-0) - x(pa%i*3-0)
	dx1a = dx1a - nbpp_cgp(ga)%x  
	dx2a = dx2a - nbpp_cgp(ga)%y  
	dx3a = dx3a - nbpp_cgp(ga)%z  


    r2a   = 1./(dx1a*dx1a + dx2a*dx2a + dx3a*dx3a)
	ra = sqrt(r2a)
    r6a   = r2a*r2a*r2a

   Vela  = crg(pa%i)*crg(pa%j)*ra
    if ( pa%LJcod .eq. 3 ) Vela = Vela*el14_scale
    V_aa  = bLJa*aLJa*aLJa*r6a*r6a 
    V_ba  = 2.0*bLJa*aLJa*r6a
    dva   = r2a*( -Vela -12.*V_aa +6.*V_ba )

    d(pa%i*3-2) = d(pa%i*3-2) - dva*dx1a
    d(pa%i*3-1) = d(pa%i*3-1) - dva*dx2a
    d(pa%i*3-0) = d(pa%i*3-0) - dva*dx3a
    d(pa%j*3-2) = d(pa%j*3-2) + dva*dx1a
    d(pa%j*3-1) = d(pa%j*3-1) + dva*dx2a
    d(pa%j*3-0) = d(pa%j*3-0) + dva*dx3a

    E%pp%el  = E%pp%el + Vela 
    E%pp%vdw = E%pp%vdw + V_aa - V_ba 
  end if

end subroutine nonbon2_pp_box
!----------------------------------------------------------------------

subroutine nonbon2_pw
! local variables
integer						:: ip,i,j,i3,j3,iaci,iacj,iLJ
real(8)						:: aLJ,bLJ,dx1,dx2,dx3,r2,r,r6,r12
real(8)						:: Vel,V_a,V_b,dv

! global variables used:
!  iac, crg, iaclib, x, d, E

do ip = 1, nbpw_pair
! for every assigned pair:
i    = nbpw(ip)%i
j    = nbpw(ip)%j
i3   = i*3-3
j3   = j*3-3
iaci = iac(i)
iacj = iac(j)
iLJ  = nbpw(ip)%LJcod
crg(i)   = crg(i)
crg(j)   = crg(j)
aLJ  = iaclib(iaci)%avdw(iLJ)+iaclib(iacj)%avdw(iLJ)
bLJ  = iaclib(iaci)%bvdw(iLJ)*iaclib(iacj)%bvdw(iLJ)
aLJ  = aLJ*aLJ
aLJ  = aLJ*aLJ*aLJ

! calculate dx and r
dx1  = x(j3+1) - x(i3+1)
dx2  = x(j3+2) - x(i3+2)
dx3  = x(j3+3) - x(i3+3)

r2   = dx1*dx1 + dx2*dx2 + dx3*dx3
r2   = 1./r2
r    = sqrt ( r2 ) 
r6   = r2*r2*r2
r12  = r6*r6

! calculate Vel and dv
Vel  = crg(i)*crg(j)*r
V_a  = bLJ*aLJ*aLJ*r12 
V_b  = 2.0*bLJ*aLJ*r6
dv   = r2*( -Vel -12.*V_a +6.*V_b )

! update forces
d(i3+1) = d(i3+1) - dv*dx1
d(i3+2) = d(i3+2) - dv*dx2
d(i3+3) = d(i3+3) - dv*dx3
d(j3+1) = d(j3+1) + dv*dx1
d(j3+2) = d(j3+2) + dv*dx2
d(j3+3) = d(j3+3) + dv*dx3

! update energies
E%pw%el  = E%pw%el + Vel       
E%pw%vdw = E%pw%vdw + V_a - V_b
end do

end subroutine nonbon2_pw

!-----------------------------------------------------------------------

subroutine nonbon2_pw_box
  ! local variables
  integer						:: ip,i,j,i3,j3,iaci,iacj,iLJ
  real(8)						:: aLJ,bLJ,dx1,dx2,dx3,r2,r,r6,r12
  real(8)						:: Vel,V_a,V_b,dv
  integer						:: group, ga, gb

  ! global variables used:
  !  iac, crg, iaclib, x, d, E


   !compute the peridocal shift for every charge group pair
	do group = 1, nbpw_cgp_pair
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


  do ip = 1, nbpw_pair
	! for every assigned pair:

	i    = nbpw(ip)%i !solute atom 
	j    = nbpw(ip)%j !solvent atom
	group = nbpw(ip)%cgp_pair
	i3   = i*3-3
	j3   = j*3-3
	iaci = iac(i)
	iacj = iac(j)
	iLJ  = nbpw(ip)%LJcod
	crg(i)   = crg(i)
	crg(j)   = crg(j)
	aLJ  = iaclib(iaci)%avdw(iLJ)+iaclib(iacj)%avdw(iLJ)
	bLJ  = iaclib(iaci)%bvdw(iLJ)*iaclib(iacj)%bvdw(iLJ)
	aLJ  = aLJ*aLJ
	aLJ  = aLJ*aLJ*aLJ

	! calculate dx and r
	dx1  = x(j3+1) - x(i3+1)
	dx2  = x(j3+2) - x(i3+2)
	dx3  = x(j3+3) - x(i3+3)
	dx1 = dx1 - nbpw_cgp(group)%x  
	dx2 = dx2 - nbpw_cgp(group)%y  
	dx3 = dx3 - nbpw_cgp(group)%z  
	

	r2   = dx1*dx1 + dx2*dx2 + dx3*dx3
	r2   = 1./r2
	r    = sqrt ( r2 ) 
	r6   = r2*r2*r2
	r12  = r6*r6

	! calculate Vel and dv
	Vel  = crg(i)*crg(j)*r
	V_a  = bLJ*aLJ*aLJ*r12 
	V_b  = 2.0*bLJ*aLJ*r6
	dv   = r2*( -Vel -12.*V_a +6.*V_b )

	! update forces
	d(i3+1) = d(i3+1) - dv*dx1
	d(i3+2) = d(i3+2) - dv*dx2
	d(i3+3) = d(i3+3) - dv*dx3
	d(j3+1) = d(j3+1) + dv*dx1
	d(j3+2) = d(j3+2) + dv*dx2
	d(j3+3) = d(j3+3) + dv*dx3

	! update energies
	E%pw%el  = E%pw%el + Vel       
	E%pw%vdw = E%pw%vdw + V_a - V_b
  end do

end subroutine nonbon2_pw_box

!-----------------------------------------------------------------------

subroutine nonbon2_qq
! local variables
integer						:: istate,jj
integer						:: ip,iq,jq,i,j,k,i3,j3,iaci,iacj,iLJ
real(8)						:: qi,qj,aLJ,bLJ,dx1,dx2,dx3,r2,r,r6,r12,r6_hc
real(8)						:: Vel,V_a,V_b,dv,el_scale

do istate = 1, nstates
! for every state:

do ip = 1, nbqq_pair(istate)
  ! for every pair:

  iq   = nbqq(ip,istate)%iq
  i    = iqseq(iq)
  j    = nbqq(ip,istate)%j
  jq   = nbqq(ip,istate)%jq
  i3   = i*3-3
  j3   = j*3-3
  iLJ  = nbqq(ip,istate)%LJcod
  qi   = qcrg(iq,istate)
  el_scale = nbqq(ip,istate)%el_scale


  if (.not. qvdw_flag) then
        iaci = iac(i)
        iacj = iac(j)
        aLJ  = iaclib(iaci)%avdw(iLJ)+iaclib(iacj)%avdw(iLJ)
        bLJ  = iaclib(iaci)%bvdw(iLJ)*iaclib(iacj)%bvdw(iLJ)
        if ( jq /= 0) then
          qj = qcrg(jq,istate)
        else
          qj = crg(j)
        end if
  else
        iaci = qiac(iq,istate)
        aLJ  = qavdw(iaci,iLJ)
        bLJ  = qbvdw(iaci,iLJ)
        if ( jq /= 0) then
          iacj = qiac(jq,istate)
          qj   = qcrg(jq,istate)
          if ( iLJ .eq. 2 ) then
                aLJ = aLJ*qavdw(iacj,iLJ)
          else
                aLJ = aLJ+qavdw(iacj,iLJ)
          end if
          bLJ = bLJ*qbvdw(iacj,iLJ)
        else
          if ( iLJ .eq. 2 ) then
                aLJ  = qavdw(iaci,1)
                bLJ  = qbvdw(iaci,1)
          end if
          iacj = iac(j)
          aLJ = aLJ+iaclib(iacj)%avdw(iLJ)
          bLJ = bLJ*iaclib(iacj)%bvdw(iLJ)
          qj = crg(j)
        end if
  end if

  ! calculate dx and r
  dx1  = x(j3+1) - x(i3+1)
  dx2  = x(j3+2) - x(i3+2)
  dx3  = x(j3+3) - x(i3+3)

  r2   = dx1*dx1 + dx2*dx2 + dx3*dx3
  r6_hc = r2*r2*r2  !for softcore
  r6   = r6_hc + sc_lookup(iq,jq+natyps,istate)  !softcore
  r6   = 1._8/r6
  r2   = 1./r2
  r    = sqrt ( r2 ) 
  r12  = r6*r6

  ! calculate Vel, V_a, V_b and dv
  Vel  = qi*qj*r*el_scale
  if ( iLJ .eq. 3 ) Vel = Vel*el14_scale
  if (qvdw_flag .and. jq /= 0 .and. iLJ .eq. 2 ) then
        V_a = aLJ*exp(-bLJ/r)
        V_b = 0.0
        dv  = r2*( -Vel -bLJ*V_a/r )*EQ(istate)%lambda
  else
        aLJ  = aLJ*aLJ
        aLJ  = aLJ*aLJ*aLJ
        V_a  = bLJ*aLJ*aLJ*r12 
        V_b  = 2.0*bLJ*aLJ*r6
        dv  = r2*( -Vel -(12.*V_a -6.*V_b)*r6*r6_hc )*EQ(istate)%lambda
  endif

  ! update forces
        d(i3+1) = d(i3+1) - dv*dx1
        d(i3+2) = d(i3+2) - dv*dx2
        d(i3+3) = d(i3+3) - dv*dx3
        d(j3+1) = d(j3+1) + dv*dx1
        d(j3+2) = d(j3+2) + dv*dx2
        d(j3+3) = d(j3+3) + dv*dx3

  ! update energies
  if ( jq /= 0 ) then
        EQ(istate)%qq%el  = EQ(istate)%qq%el + Vel
        EQ(istate)%qq%vdw = EQ(istate)%qq%vdw + V_a - V_b
	if (use_excluded_groups) then
                do jj = 1, ngroups_gc
                call set_gc_energies(i,j,Vel,(V_a-V_b),EQ_gc(jj,istate)%qq%el &
                ,EQ_gc(jj,istate)%qq%vdw,ST_gc(jj)%gcmask%mask)
                end do
	end if
  else
        EQ(istate)%qp%el  = EQ(istate)%qp%el + Vel
        EQ(istate)%qp%vdw = EQ(istate)%qp%vdw + V_a - V_b
	if (use_excluded_groups) then
                do jj = 1, ngroups_gc
                call set_gc_energies(i,j,Vel,(V_a-V_b),EQ_gc(jj,istate)%qp%el &
                ,EQ_gc(jj,istate)%qp%vdw,ST_gc(jj)%gcmask%mask)
                end do
	end if
  end if
end do ! ip

end do ! istate
end subroutine nonbon2_qq

!-----------------------------------------------------------------------
subroutine nonbon2_qq_lib_charges
! local variables
integer						:: istate,jj
integer						:: ip,iq,jq,i,j,k,i3,j3,iaci,iacj,iLJ
real(8)						:: qi,qj,aLJ,bLJ,dx1,dx2,dx3,r2,r,r6,r12,r6_hc
real(8)						:: Vel,V_a,V_b,dv,el_scale

do istate = 1, nstates
! for every state:

do ip = 1, nbqq_pair(istate)
  ! for every pair:

  iq   = nbqq(ip,istate)%iq
  i    = iqseq(iq)
  j    = nbqq(ip,istate)%j
  jq   = nbqq(ip,istate)%jq
  i3   = i*3-3
  j3   = j*3-3
  iLJ  = nbqq(ip,istate)%LJcod
  qi   = crg(i)
  el_scale = nbqq(ip,istate)%el_scale

  if (.not. qvdw_flag) then
        iaci = iac(i)
        iacj = iac(j)
        aLJ  = iaclib(iaci)%avdw(iLJ)+iaclib(iacj)%avdw(iLJ)
        bLJ  = iaclib(iaci)%bvdw(iLJ)*iaclib(iacj)%bvdw(iLJ)
        qj = crg(j)
  else
        iaci = qiac(iq,istate)
        aLJ  = qavdw(iaci,iLJ)
        bLJ  = qbvdw(iaci,iLJ)
        if ( jq /= 0) then
          iacj = qiac(jq,istate)
          if ( iLJ .eq. 2 ) then
                aLJ = aLJ*qavdw(iacj,iLJ)
          else
                aLJ = aLJ+qavdw(iacj,iLJ)
          end if
          bLJ = bLJ*qbvdw(iacj,iLJ)
        else
          if ( iLJ .eq. 2 ) then
                aLJ  = qavdw(iaci,1)
                bLJ  = qbvdw(iaci,1)
          end if
          iacj = iac(j)
          aLJ = aLJ+iaclib(iacj)%avdw(iLJ)
          bLJ = bLJ*iaclib(iacj)%bvdw(iLJ)
        end if
        qj = crg(j)
  end if

  ! calculate dx and r
  dx1  = x(j3+1) - x(i3+1)
  dx2  = x(j3+2) - x(i3+2)
  dx3  = x(j3+3) - x(i3+3)

  r2   = dx1*dx1 + dx2*dx2 + dx3*dx3
  r6_hc = r2*r2*r2   !needed for softcore
  r6   = r6_hc + sc_lookup(iq,jq+natyps,istate)  !softcore
  r6   = 1._8/r6
  r12  = r6*r6
  r2   = 1./r2
  r    = sqrt ( r2 ) 

  ! calculate Vel, V_a, V_b and dv
  Vel  = qi*qj*r*el_scale
  if ( iLJ .eq. 3 ) Vel = Vel*el14_scale
  if (qvdw_flag .and. jq /= 0 .and. iLJ .eq. 2 ) then
        V_a = aLJ*exp(-bLJ/r)
        V_b = qbvdw(iaci,1)*qbvdw(iacj,1)*r6
        dv  = r2*( -Vel -bLJ*V_a/r +6.*V_b )*EQ(istate)%lambda
        !
        ! ---       change here to exclude 1/r6 attraction
        !
        !	       V_b = 0.0
        !	       dv  = r2*( -Vel -bLJ*V_a/r )*EQ(istate)%lambda
  else
        aLJ  = aLJ*aLJ
        aLJ  = aLJ*aLJ*aLJ
        V_a  = bLJ*aLJ*aLJ*r12 
        V_b  = 2.0*bLJ*aLJ*r6
        dv  = r2*( -Vel -(12.*V_a -6.*V_b)*r6*r6_hc )*EQ(istate)%lambda
  endif

  ! update forces
  d(i3+1) = d(i3+1) - dv*dx1
  d(i3+2) = d(i3+2) - dv*dx2
  d(i3+3) = d(i3+3) - dv*dx3
  d(j3+1) = d(j3+1) + dv*dx1
  d(j3+2) = d(j3+2) + dv*dx2
  d(j3+3) = d(j3+3) + dv*dx3

  ! update energies
  if ( jq /= 0 ) then
        EQ(istate)%qq%el  = EQ(istate)%qq%el + Vel
        EQ(istate)%qq%vdw = EQ(istate)%qq%vdw + V_a - V_b
	if (use_excluded_groups) then
                do jj = 1, ngroups_gc
                call set_gc_energies(i,j,Vel,(V_a-V_b),EQ_gc(jj,istate)%qq%el &
                ,EQ_gc(jj,istate)%qq%vdw,ST_gc(jj)%gcmask%mask)
                end do
	end if
  else
        EQ(istate)%qp%el  = EQ(istate)%qp%el + Vel
        EQ(istate)%qp%vdw = EQ(istate)%qp%vdw + V_a - V_b
	if (use_excluded_groups) then
                do jj = 1, ngroups_gc
                call set_gc_energies(i,j,Vel,(V_a-V_b),EQ_gc(jj,istate)%qp%el &
                ,EQ_gc(jj,istate)%qp%vdw,ST_gc(jj)%gcmask%mask)
                end do
	end if
  end if
end do ! ip

end do ! istate
end subroutine nonbon2_qq_lib_charges

!-----------------------------------------------------------------------

subroutine nonbon2_qp
! local variables
integer						:: ip,iq,i,j,i3,j3,iaci,iacj,iLJ
integer						:: istate,jj
real(8)						:: aLJ,bLJ,dx1,dx2,dx3,r2,r,r6,r6_hc
real(8)						:: Vel,V_a,V_b,dv

! global variables used:
!  iqseq, iac, crg, x, nstates, qvdw_flag, iaclib, qiac, qavdw, qbvdw, qcrg, el14_scale, EQ, d, nat_solute

do ip = 1, nbqp_pair
! for every assigned q-s pair:

! init state-invariant variables:
iq   = nbqp(ip)%i
i    = iqseq(iq)
j    = nbqp(ip)%j
i3   = i*3-3
j3   = j*3-3
iacj = iac(j)
iLJ  = nbqp(ip)%LJcod


dx1  = x(j3+1) - x(i3+1)
dx2  = x(j3+2) - x(i3+2)
dx3  = x(j3+3) - x(i3+3)


r2   = dx1*dx1 + dx2*dx2 + dx3*dx3

r6   = r2*r2*r2
r6_hc = r6     !for softcore

r2   = 1._8/r2
r    = sqrt(r2)



do istate = 1, nstates
  ! for every state:

  ! calculate iaci, aLJ and bLJ
  if (.not. qvdw_flag) then
        iaci = iac(i)
        aLJ  = iaclib(iaci)%avdw(iLJ)
        bLJ  = iaclib(iaci)%bvdw(iLJ)
		r6 = 1._8/r6_hc
  else
        iaci = qiac(iq,istate)
		if ( iLJ .eq. 2 ) then
          aLJ  = qavdw(iaci,1)
          bLJ  = qbvdw(iaci,1)
        else
          aLJ  = qavdw(iaci,iLJ)
          bLJ  = qbvdw(iaci,iLJ)
		end if

		r6 = r6_hc + sc_lookup(iq,iacj,istate) !this is softcore
		r6 = 1._8/r6
  end if
  aLJ = aLJ+iaclib(iacj)%avdw(iLJ)
  bLJ = bLJ*iaclib(iacj)%bvdw(iLJ)
  aLJ  = aLJ*aLJ
  aLJ  = aLJ*aLJ*aLJ

  ! calculate qi, Vel, V_a, V_b and dv
  Vel  = qcrg(iq,istate)*crg(j)*r
  if ( iLJ .eq. 3 ) Vel = Vel*el14_scale
  V_a  = bLJ*aLJ*aLJ*r6*r6
  V_b  = 2.0*bLJ*aLJ*r6
  dv   = r2*( -Vel -(12.*V_a -6.*V_b)*r6*r6_hc )*EQ(istate)%lambda   !softcore r6*r6_hc is (r^6/(r^6+alpha))

  ! update forces
  d(i3+1) = d(i3+1) - dv*dx1
  d(i3+2) = d(i3+2) - dv*dx2
  d(i3+3) = d(i3+3) - dv*dx3
  d(j3+1) = d(j3+1) + dv*dx1
  d(j3+2) = d(j3+2) + dv*dx2
  d(j3+3) = d(j3+3) + dv*dx3

  ! update q-protein or q-water energies
        EQ(istate)%qp%el  = EQ(istate)%qp%el + Vel
        EQ(istate)%qp%vdw = EQ(istate)%qp%vdw + V_a - V_b
        if (use_excluded_groups) then
                do jj = 1, ngroups_gc
                call set_gc_energies(i,j,Vel,(V_a-V_b),EQ_gc(jj,istate)%qp%el &
                ,EQ_gc(jj,istate)%qp%vdw,ST_gc(jj)%gcmask%mask)
                end do 
        end if
end do ! istate

end do
end subroutine nonbon2_qp
!-----------------------------------------------------------------------

!******PWadded 2001-10-23
subroutine nonbon2_qp_box
  ! local variables
  integer						:: ip,iq,i,j,i3,j3,iaci,iacj,iLJ
  integer						:: istate,jj
  real(8)						:: aLJ,bLJ,dx1,dx2,dx3,r2,r,r6,r6_hc
  real(8)						:: Vel,V_a,V_b,dv
  integer						:: group, gr, ia

  ! global variables used:
  !  iqseq, iac, crg, x, nstates, qvdw_flag, iaclib, qiac, qavdw, qbvdw, qcrg, el14_scale, EQ, d, nat_solute


	!compute the peridocal shift for every charge group pair
	do gr = 1, nbqp_cgp_pair
		ia = nbqp_cgp(gr)%i !atom index for the atom
		
		!the distance between the two switching atoms
		dx1 = x(3*ia-2) - x(3*qswitch-2)
		dx2 = x(3*ia-1) - x(3*qswitch-1)
		dx3 = x(3*ia  ) - x(3*qswitch  )

		nbqp_cgp(gr)%x = boxlength(1)*nint( dx1*inv_boxl(1) )
		nbqp_cgp(gr)%y = boxlength(2)*nint( dx2*inv_boxl(2) )	
		nbqp_cgp(gr)%z = boxlength(3)*nint( dx3*inv_boxl(3) )

	end do

  do ip = 1, nbqp_pair
	! for every assigned q-s pair:

   	! init state-invariant variables:
	iq   = nbqp(ip)%i
	i    = iqseq(iq)
	j    = nbqp(ip)%j
	i3   = i*3-3
	j3   = j*3-3
	iacj = iac(j)
	iLJ  = nbqp(ip)%LJcod
	group = nbqp(ip)%cgp_pair

	dx1  = x(j3+1) - x(i3+1)
	dx2  = x(j3+2) - x(i3+2)
	dx3  = x(j3+3) - x(i3+3)
	dx1 = dx1 - nbqp_cgp(group)%x  
	dx2 = dx2 - nbqp_cgp(group)%y  
	dx3 = dx3 - nbqp_cgp(group)%z  
	

	r2   = dx1*dx1 + dx2*dx2 + dx3*dx3
	r6_hc = r2*r2*r2  !for softcore

	r2   = 1._8/r2
	r    = sqrt(r2)

	do istate = 1, nstates
	  ! for every state:

	  ! calculate iaci, aLJ and bLJ
	  if (.not. qvdw_flag) then
		iaci = iac(i)
		aLJ  = iaclib(iaci)%avdw(iLJ)
		bLJ  = iaclib(iaci)%bvdw(iLJ)
		r6 = 1._8/r6_hc
	  else
		iaci = qiac(iq,istate)
        if ( iLJ .eq. 2 ) then
		  aLJ  = qavdw(iaci,1)
		  bLJ  = qbvdw(iaci,1)
		else
		  aLJ  = qavdw(iaci,iLJ)
		  bLJ  = qbvdw(iaci,iLJ)
        end if
		r6 = r6_hc + sc_lookup(iq,iacj,istate)
		r6 = 1._8/r6
	  end if
	  aLJ = aLJ+iaclib(iacj)%avdw(iLJ)
	  bLJ = bLJ*iaclib(iacj)%bvdw(iLJ)
	  aLJ  = aLJ*aLJ
	  aLJ  = aLJ*aLJ*aLJ

	  ! calculate qi, Vel, V_a, V_b and dv
	  Vel  = qcrg(iq,istate)*crg(j)*r
	  if ( iLJ .eq. 3 ) Vel = Vel*el14_scale
	  V_a  = bLJ*aLJ*aLJ*r6*r6
	  V_b  = 2.0*bLJ*aLJ*r6
	  dv   = r2*( -Vel -(12.*V_a -6.*V_b)*r6*r6_hc )*EQ(istate)%lambda

	  ! update forces
	  d(i3+1) = d(i3+1) - dv*dx1
	  d(i3+2) = d(i3+2) - dv*dx2
	  d(i3+3) = d(i3+3) - dv*dx3
	  d(j3+1) = d(j3+1) + dv*dx1
	  d(j3+2) = d(j3+2) + dv*dx2
	  d(j3+3) = d(j3+3) + dv*dx3
	
	  ! update q-protein or q-water energies
		EQ(istate)%qp%el  = EQ(istate)%qp%el + Vel
		EQ(istate)%qp%vdw = EQ(istate)%qp%vdw + V_a - V_b
        if (use_excluded_groups) then
                do jj = 1, ngroups_gc
                call set_gc_energies(i,j,Vel,(V_a-V_b),EQ_gc(jj,istate)%qp%el &
                ,EQ_gc(jj,istate)%qp%vdw,ST_gc(jj)%gcmask%mask)
                end do
        end if
	end do ! istate

  end do
end subroutine nonbon2_qp_box

!-----------------------------------------------------------------------

subroutine nonbon2_qw

! local variables
integer						:: jw,iq,i,j,iLJO, iLJH, iaci
integer						:: istate
real(8)						::	aLJO, bLJO, aLJH, bLJH
real(8)						::	dxO, dyO, dzO, dxH1, dyH1, dzH1, dxH2, dyH2, dzH2
real(8)						::	rO, r2O, r6O, rH1, r2H1, r6H1, rH2, r2H2, r6H2,r6O_hc,r6H1_hc,r6H2_hc
real(8)						::	VelO, VelH1, VelH2, dvO, dvH1, dvH2
real(8)						:: V_ao, V_bo, V_ah1, V_bh1, V_ah2, V_bh2 
real(8), save				::	aO(2), bO(2), aH(2), bH(2)
integer, save				::	iac_ow, iac_hw
! global variables used:
!  iqseq, iac, crg, x, nstates, qvdw_flag, iaclib, qiac, qavdw, qbvdw, qcrg, el14_scale, EQ, d, nat_solute

if(iac_ow == 0) then !set first time
        iac_ow = iac(nat_solute + 1)
        iac_hw = iac(nat_solute + 2)
        aO(1:2) = iaclib(iac_ow)%avdw(1:2)
        bO(1:2) = iaclib(iac_ow)%bvdw(1:2)
        aH(1:2) = iaclib(iac_hw)%avdw(1:2)
        bH(1:2) = iaclib(iac_hw)%bvdw(1:2)
end if

!loop over listed waters
do jw = 1, nbqw_pair
        j = nbqw(jw) !top. # of O in water iw
        !loop over all q-atoms
        do iq = 1, nqat
                i = iqseq(iq)
                dxO  = x(3*j-2) - x(3*i-2)
                dyO  = x(3*j-1) - x(3*i-1)
                dzO  = x(3*j  ) - x(3*i  )
                dxH1 = x(3*j+1) - x(3*i-2)
                dyH1 = x(3*j+2) - x(3*i-1)
                dzH1 = x(3*j+3) - x(3*i  )
                dxH2 = x(3*j+4) - x(3*i-2)
                dyH2 = x(3*j+5) - x(3*i-1)
                dzH2 = x(3*j+6) - x(3*i  )
                r2O  = 1._8/(dxO*dxO + dyO*dyO + dzO*dzO)
                r2H1 = 1._8/(dxH1*dxH1 + dyH1*dyH1 + dzH1*dzH1)
                r2H2 = 1._8/(dxH2*dxH2 + dyH2*dyH2 + dzH2*dzH2)
                rO = sqrt(r2O)
                r6O = r2O*r2O*r2O
				r6O_hc = r6O	!softcore
                rH1 = sqrt(r2H1)
                r6H1 = r2H1*r2H1*r2H1
				r6H1_hc = r6H1	!softcore
                rH2 = sqrt(r2H2)
                r6H2 = r2H2*r2H2*r2H2
				r6H2_hc = r6H2	!softcore

                iaci = iac(i)

                !reset potential
                dvO = 0
                dvH1 = 0
                dvH2 = 0
                iLJO = LJcod(iac_ow, iaci)
                iLJH = LJcod(iac_hw, iaci)
                if(.not. qvdw_flag) then
                        !use same LJ params for all states
                        aLJO  = iaclib(iaci)%avdw(iLJO)+aO(iLJO)
                        bLJO  = iaclib(iaci)%bvdw(iLJO)*bO(iLJO)
                        aLJH  = iaclib(iaci)%avdw(iLJH)+aH(iLJH)
                        bLJH  = iaclib(iaci)%bvdw(iLJH)*bH(iLJH) 
                        aLJO = aLJO * aLJO
                        aLJO = aLJO * aLJO * aLJO
                        aLJH  = aLJH  * aLJH 
                        aLJH  = aLJH  * aLJH  * aLJH 
                        V_aO = bLJO*aLJO*aLJO*r6O*r6O 
                        V_bO = 2.0*bLJO*aLJO*r6O
                        V_aH1= bLJH*aLJH*aLJH*r6H1*r6H1
                        V_bH1= 2.0*bLJH*aLJH*r6H1
                        V_aH2= bLJH*aLJH*aLJH*r6H2*r6H2
                        V_bH2= 2.0*bLJH*aLJH*r6H2
                end if
                do istate = 1, nstates ! for every state:

                        ! set new LJ params if Q-atom types are used
                        if (qvdw_flag) then

								r6O = 1._8/r6O_hc
								r6O = r6O + sc_lookup(iq,iac_ow,istate)   !softcore
								r6O = 1._8/r6O
								r6H1 = 1._8/r6H1_hc
								r6H1 = r6H1 + sc_lookup(iq,iac_hw,istate)   !softcore
								r6H1 = 1._8/r6H1
								r6H2 = 1._8/r6H2_hc
								r6H2 = r6H2 + sc_lookup(iq,iac_hw,istate)   !softcore
								r6H2 = 1._8/r6H2
                                aLJO  = qavdw(qiac(iq,istate),1)+aO(iLJO)
                                bLJO  = qbvdw(qiac(iq,istate),1)*bO(iLJO)
                                aLJH  = qavdw(qiac(iq,istate),1)+aH(iLJH)
                                bLJH  = qbvdw(qiac(iq,istate),1)*bH(iLJH)
                                aLJO = aLJO * aLJO
                                aLJO = aLJO * aLJO * aLJO
                                aLJH  = aLJH  * aLJH 
                                aLJH  = aLJH  * aLJH  * aLJH 
                                V_aO = bLJO*aLJO*aLJO*r6O*r6O 
                                V_bO = 2.0*bLJO*aLJO*r6O
                                V_aH1= bLJH*aLJH*aLJH*r6H1*r6H1
                                V_bH1= 2.0*bLJH*aLJH*r6H1
                                V_aH2= bLJH*aLJH*aLJH*r6H2*r6H2
                                V_bH2= 2.0*bLJH*aLJH*r6H2
                        end if


                        ! calculate qi, Vel, V_a, V_b and dv
                        VelO = crg_ow*qcrg(iq,istate)*rO
                        VelH1 = crg_hw*qcrg(iq,istate)*rH1
                        VelH2 = crg_hw*qcrg(iq,istate)*rH2
                        dvO  = dvO  + r2O *( -VelO  -(12.*V_aO  -6.*V_bO )*r6O/r6O_hc)*EQ(istate)%lambda
                        dvH1 = dvH1 + r2H1*( -VelH1 -(12.*V_aH1 -6.*V_bH1)*r6H1/r6H1_hc)*EQ(istate)%lambda
                        dvH2 = dvH2 + r2H2*( -VelH2 -(12.*V_aH2 -6.*V_bH2)*r6H2/r6H2_hc)*EQ(istate)%lambda
                        ! update q-water energies
                        EQ(istate)%qw%el  = EQ(istate)%qw%el + VelO + VelH1 + VelH2
                        EQ(istate)%qw%vdw = EQ(istate)%qw%vdw + V_aO + V_aH1 + V_aH2 - V_bO - V_bH1 - V_bH2 
                end do !istate

                ! update forces on Q-atom
                d(3*i-2) = d(3*i-2) - dvO*dxO - dvH1*dxH1 - dvH2*dxH2
                d(3*i-1) = d(3*i-1) - dvO*dyO - dvH1*dyH1 - dvH2*dyH2
                d(3*i  ) = d(3*i  ) - dvO*dzO - dvH1*dzH1 - dvH2*dzH2

                ! update forces on water
                d(3*j-2) = d(3*j-2) + dvO*dxO
                d(3*j-1) = d(3*j-1) + dvO*dyO
                d(3*j  ) = d(3*j  ) + dvO*dzO
                d(3*j+1) = d(3*j+1) + dvH1*dxH1
                d(3*j+2) = d(3*j+2) + dvH1*dyH1
                d(3*j+3) = d(3*j+3) + dvH1*dzH1
                d(3*j+4) = d(3*j+4) + dvH2*dxH2
                d(3*j+5) = d(3*j+5) + dvH2*dyH2
                d(3*j+6) = d(3*j+6) + dvH2*dzH2

        end do !iq
end do !jw
end subroutine nonbon2_qw

!-----------------------------------------------------------------------
!******PWadded 2001-10-23
subroutine nonbon2_qw_box
	! local variables
	integer						:: jw,iq,i,j,iLJO, iLJH, iaci
	integer						:: istate
	real(8)						::	aLJO, bLJO, aLJH, bLJH
	real(8)						::	dxO, dyO, dzO, dxH1, dyH1, dzH1, dxH2, dyH2, dzH2
	real(8)						::	rO, r2O, r6O, rH1, r2H1, r6H1, rH2, r2H2, r6H2,r6O_hc,r6H1_hc,r6H2_hc
	real(8)						::	VelO, VelH1, VelH2, dvO, dvH1, dvH2
	real(8)						:: V_ao, V_bo, V_ah1, V_bh1, V_ah2, V_bh2
	real(8)						:: boxshiftx, boxshifty, boxshiftz, dx, dy, dz 
	real(8), save				::	aO(2), bO(2), aH(2), bH(2)
	integer, save				::	iac_ow, iac_hw
	! global variables used:
	!  iqseq, iac, crg, x, nstates, qvdw_flag, iaclib, qiac, qavdw, qbvdw, qcrg, el14_scale, EQ, d, nat_solute
  
	if(iac_ow == 0) then !set first time
		iac_ow = iac(nat_solute + 1)
		iac_hw = iac(nat_solute + 2)
		aO(1:2) = iaclib(iac_ow)%avdw(1:2)
		bO(1:2) = iaclib(iac_ow)%bvdw(1:2)
		aH(1:2) = iaclib(iac_hw)%avdw(1:2)
		bH(1:2) = iaclib(iac_hw)%bvdw(1:2)
	end if

	!loop over listed waters
	do jw = 1, nbqw_pair
		j = nbqw(jw) !top. # of O in water iw

		!compute the periodical shift
		dx = x(3*j-2) - x(3*qswitch-2)
		dy = x(3*j-1) - x(3*qswitch-1)
		dz = x(3*j  ) - x(3*qswitch  )
		boxshiftx = boxlength(1)*nint(dx*inv_boxl(1))
		boxshifty = boxlength(2)*nint(dy*inv_boxl(2))
		boxshiftz = boxlength(3)*nint(dz*inv_boxl(3))

		!loop over all q-atoms
		do iq = 1, nqat
			i = iqseq(iq)
			dxO  = x(3*j-2) - x(3*i-2)
			dyO  = x(3*j-1) - x(3*i-1)
			dzO  = x(3*j  ) - x(3*i  )
			dxH1 = x(3*j+1) - x(3*i-2)
			dyH1 = x(3*j+2) - x(3*i-1)
			dzH1 = x(3*j+3) - x(3*i  )
			dxH2 = x(3*j+4) - x(3*i-2)
			dyH2 = x(3*j+5) - x(3*i-1)
			dzH2 = x(3*j+6) - x(3*i  )
			dxO = dxO - boxshiftx
			dyO = dyO - boxshifty  
			dzO = dzO - boxshiftz  
			dxH1 = dxH1 - boxshiftx  
			dyH1 = dyH1 - boxshifty  
			dzH1 = dzH1 - boxshiftz  
			dxH2 = dxH2 - boxshiftx 
			dyH2 = dyH2 - boxshifty  
			dzH2 = dzH2 - boxshiftz  
			r2O  = 1._8/(dxO*dxO + dyO*dyO + dzO*dzO)
			r2H1 = 1._8/(dxH1*dxH1 + dyH1*dyH1 + dzH1*dzH1)
			r2H2 = 1._8/(dxH2*dxH2 + dyH2*dyH2 + dzH2*dzH2)
			rO = sqrt(r2O)
			r6O = r2O*r2O*r2O
			r6O_hc = r6O    !needed for softcore
			rH1 = sqrt(r2H1)
			r6H1 = r2H1*r2H1*r2H1
			r6H1_hc = r6H1    !needed for softcore
			rH2 = sqrt(r2H2)
			r6H2 = r2H2*r2H2*r2H2
			r6H2_hc = r6H2    !needed for softcore
			iaci = iac(i)
			
			!reset potential
			dvO = 0
			dvH1 = 0
			dvH2 = 0
			iLJO = LJcod(iac_ow, iaci)
			iLJH = LJcod(iac_hw, iaci)
			if(.not. qvdw_flag) then
				!use same LJ params for all states
				aLJO  = iaclib(iaci)%avdw(iLJO)+aO(iLJO)
				bLJO  = iaclib(iaci)%bvdw(iLJO)*bO(iLJO)
				aLJH  = iaclib(iaci)%avdw(iLJH)+aH(iLJH)
				bLJH  = iaclib(iaci)%bvdw(iLJH)*bH(iLJH) 
				aLJO = aLJO * aLJO
				aLJO = aLJO * aLJO * aLJO
				aLJH  = aLJH  * aLJH 
				aLJH  = aLJH  * aLJH  * aLJH 
				V_aO = bLJO*aLJO*aLJO*r6O*r6O 
				V_bO = 2.0*bLJO*aLJO*r6O
				V_aH1= bLJH*aLJH*aLJH*r6H1*r6H1
				V_bH1= 2.0*bLJH*aLJH*r6H1
				V_aH2= bLJH*aLJH*aLJH*r6H2*r6H2
				V_bH2= 2.0*bLJH*aLJH*r6H2
			end if
			do istate = 1, nstates ! for every state:
			
				! set new LJ params if Q-atom types are used
				if (qvdw_flag) then
					r6O = 1._8/r6O_hc
					r6O = r6O + sc_lookup(iq,iac_ow,istate)   !softcore
					r6O = 1._8/r6O
					r6H1 = 1._8/r6H1_hc
					r6H1 = r6H1 + sc_lookup(iq,iac_hw,istate)   !softcore
					r6H1 = 1._8/r6H1
					r6H2 = 1._8/r6H2_hc
					r6H2 = r6H2 + sc_lookup(iq,iac_hw,istate)   !softcore
					r6H2 = 1._8/r6H2
					aLJO  = qavdw(qiac(iq,istate),1)+aO(iLJO)
					bLJO  = qbvdw(qiac(iq,istate),1)*bO(iLJO)
					aLJH  = qavdw(qiac(iq,istate),1)+aH(iLJH)
					bLJH  = qbvdw(qiac(iq,istate),1)*bH(iLJH)
					aLJO = aLJO * aLJO
					aLJO = aLJO * aLJO * aLJO
					aLJH  = aLJH  * aLJH 
					aLJH  = aLJH  * aLJH  * aLJH 
					V_aO = bLJO*aLJO*aLJO*r6O*r6O 
					V_bO = 2.0*bLJO*aLJO*r6O
					V_aH1= bLJH*aLJH*aLJH*r6H1*r6H1
					V_bH1= 2.0*bLJH*aLJH*r6H1
					V_aH2= bLJH*aLJH*aLJH*r6H2*r6H2
					V_bH2= 2.0*bLJH*aLJH*r6H2
				end if


				! calculate qi, Vel, V_a, V_b and dv
				VelO = crg_ow*qcrg(iq,istate)*rO
				VelH1 = crg_hw*qcrg(iq,istate)*rH1
				VelH2 = crg_hw*qcrg(iq,istate)*rH2

				dvO  = dvO  + r2O *( -VelO  -(12.*V_aO  -6.*V_bO )*r6O/r6O_hc)*EQ(istate)%lambda    !r6O/r6O_hc softcore
				dvH1 = dvH1 + r2H1*( -VelH1 -(12.*V_aH1 -6.*V_bH1)*r6H1/r6H1_hc)*EQ(istate)%lambda  !r6H1/r6H1_hc softcore
				dvH2 = dvH2 + r2H2*( -VelH2 -(12.*V_aH2 -6.*V_bH2)*r6H2/r6H2_hc)*EQ(istate)%lambda  !r6H2/r6H2_hc softcore
				! update q-water energies
				EQ(istate)%qw%el  = EQ(istate)%qw%el + VelO + VelH1 + VelH2
				EQ(istate)%qw%vdw = EQ(istate)%qw%vdw + V_aO + V_aH1 + V_aH2 - V_bO - V_bH1 - V_bH2 
			end do !istate
			
			! update forces on Q-atom
			d(3*i-2) = d(3*i-2) - dvO*dxO - dvH1*dxH1 - dvH2*dxH2
			d(3*i-1) = d(3*i-1) - dvO*dyO - dvH1*dyH1 - dvH2*dyH2
			d(3*i  ) = d(3*i  ) - dvO*dzO - dvH1*dzH1 - dvH2*dzH2

			! update forces on water
			d(3*j-2) = d(3*j-2) + dvO*dxO
			d(3*j-1) = d(3*j-1) + dvO*dyO
			d(3*j  ) = d(3*j  ) + dvO*dzO
			d(3*j+1) = d(3*j+1) + dvH1*dxH1
			d(3*j+2) = d(3*j+2) + dvH1*dyH1
			d(3*j+3) = d(3*j+3) + dvH1*dzH1
			d(3*j+4) = d(3*j+4) + dvH2*dxH2
			d(3*j+5) = d(3*j+5) + dvH2*dyH2
			d(3*j+6) = d(3*j+6) + dvH2*dzH2

		end do !iq
	end do !jw
end subroutine nonbon2_qw_box

!----------------------------------------------------------------------------------------------------

subroutine nonbon2_ww
! local variables
integer						:: iw,ip,i,j,i3,j3,ia
integer						:: iaci,iacj,iLJ,ja
real(8)						:: aLJ,bLJ
integer						:: ipstart
real(8)						:: dx1,dx2,dx3,r2,r,r6,r12
real(8)						:: Vel,V_a,V_b,dv

! global variables used:
!  nat_solute, iac, crg, ljcod, iaclib, x, d, E

ipstart = 1

do iw = calculation_assignment%ww%start, calculation_assignment%ww%end
! for every assigned water molecule:

do ia = 1, 3
  ! for every atom of the current water molecule:
  i    = nat_solute+3*(iw-1)+ia
  i3   = i*3-3
  iaci = iac(i)
  crg(i)   = crg(i)

  ip = ipstart
  do while (nbww(ip) .ne. 0)
        ! loop over the interactions with other water molecules

        ! X-O
        j    = nbww(ip)
        j3   = j*3-3
        iacj = iac(j)
        iLJ  = ljcod(iac(i),iac(j))
        crg(j)   = crg(j)
        aLJ  = iaclib(iaci)%avdw(iLJ)+iaclib(iacj)%avdw(iLJ)
        bLJ  = iaclib(iaci)%bvdw(iLJ)*iaclib(iacj)%bvdw(iLJ)
        aLJ  = aLJ*aLJ
        aLJ  = aLJ*aLJ*aLJ
        dx1  = x(j3+1) - x(i3+1)
        dx2  = x(j3+2) - x(i3+2)
        dx3  = x(j3+3) - x(i3+3)
        r2   = dx1*dx1 + dx2*dx2 + dx3*dx3
        r2   = 1./r2
        r    = sqrt ( r2 ) 
        r6   = r2*r2*r2
        r12  = r6*r6
        Vel  = crg(i)*crg(j)*r
        V_a  = bLJ*aLJ*aLJ*r12 
        V_b  = 2.0*bLJ*aLJ*r6
        dv   = r2*( -Vel -12.*V_a +6.*V_b )
        d(i3+1) = d(i3+1) - dv*dx1
        d(i3+2) = d(i3+2) - dv*dx2
        d(i3+3) = d(i3+3) - dv*dx3
        d(j3+1) = d(j3+1) + dv*dx1
        d(j3+2) = d(j3+2) + dv*dx2
        d(j3+3) = d(j3+3) + dv*dx3
        E%ww%el  = E%ww%el + Vel       
        E%ww%vdw = E%ww%vdw + V_a - V_b

        ! X-H1
        j    = j + 1
        j3   = j3 + 3
        iacj = iac(j)
        iLJ  = ljcod(iac(i),iac(j))
        crg(j)   = crg(j)
        aLJ  = iaclib(iaci)%avdw(iLJ)+iaclib(iacj)%avdw(iLJ)
        bLJ  = iaclib(iaci)%bvdw(iLJ)*iaclib(iacj)%bvdw(iLJ)
        aLJ  = aLJ*aLJ
        aLJ  = aLJ*aLJ*aLJ
        dx1  = x(j3+1) - x(i3+1)
        dx2  = x(j3+2) - x(i3+2)
        dx3  = x(j3+3) - x(i3+3)
        r2   = dx1*dx1 + dx2*dx2 + dx3*dx3
        r2   = 1./r2
        r    = sqrt ( r2 ) 
        r6   = r2*r2*r2
        r12  = r6*r6
        Vel  = crg(i)*crg(j)*r
        V_a  = bLJ*aLJ*aLJ*r12 
        V_b  = 2.0*bLJ*aLJ*r6
        dv   = r2*( -Vel -12.*V_a +6.*V_b )
        d(i3+1) = d(i3+1) - dv*dx1
        d(i3+2) = d(i3+2) - dv*dx2
        d(i3+3) = d(i3+3) - dv*dx3
        d(j3+1) = d(j3+1) + dv*dx1
        d(j3+2) = d(j3+2) + dv*dx2
        d(j3+3) = d(j3+3) + dv*dx3
        E%ww%el  = E%ww%el + Vel       
        E%ww%vdw = E%ww%vdw + V_a - V_b

        ! X-H2
        j    = j + 1
        j3   = j3 + 3
        iacj = iac(j)
        iLJ  = ljcod(iac(i),iac(j))
        crg(j)   = crg(j)
        aLJ  = iaclib(iaci)%avdw(iLJ)+iaclib(iacj)%avdw(iLJ)
        bLJ  = iaclib(iaci)%bvdw(iLJ)*iaclib(iacj)%bvdw(iLJ)
        aLJ  = aLJ*aLJ
        aLJ  = aLJ*aLJ*aLJ
        dx1  = x(j3+1) - x(i3+1)
        dx2  = x(j3+2) - x(i3+2)
        dx3  = x(j3+3) - x(i3+3)
        r2   = dx1*dx1 + dx2*dx2 + dx3*dx3
        r2   = 1./r2
        r    = sqrt ( r2 ) 
        r6   = r2*r2*r2
        r12  = r6*r6
        Vel  = crg(i)*crg(j)*r
        V_a  = bLJ*aLJ*aLJ*r12 
        V_b  = 2.0*bLJ*aLJ*r6
        dv   = r2*( -Vel -12.*V_a +6.*V_b )
        d(i3+1) = d(i3+1) - dv*dx1
        d(i3+2) = d(i3+2) - dv*dx2
        d(i3+3) = d(i3+3) - dv*dx3
        d(j3+1) = d(j3+1) + dv*dx1
        d(j3+2) = d(j3+2) + dv*dx2
        d(j3+3) = d(j3+3) + dv*dx3
        E%ww%el  = E%ww%el + Vel       
        E%ww%vdw = E%ww%vdw + V_a - V_b

        ip = ip + 1
  end do
end do

ipstart = ip +1						! skip over the 0
end do

end subroutine nonbon2_ww

!----------------------------------------------------------------------------------------------------
subroutine nonbon2_ww_box
  ! local variables
  integer						:: iw,ip,i,j,i3,j3,ia
  integer						:: iaci,iacj,iLJ,ja
  real(8)						:: aLJ,bLJ
  integer						:: ipstart
  real(8)						:: dx1,dx2,dx3,r2,r,r6,r12
  real(8)						:: Vel,V_a,V_b,dv
  real(8)						:: ds1, ds2, ds3

  ! global variables used:
  !  nat_solute, iac, crg, ljcod, iaclib, x, d, E

  ipstart = 1

  do iw = calculation_assignment%ww%start, calculation_assignment%ww%end
	! for every assigned water molecule:

	do ia = 1, 3
	  ! for every atom of the current water molecule:
	  i    = nat_solute+3*(iw-1)+ia
	  i3   = i*3-3
	  iaci = iac(i)
	  crg(i)   = crg(i)

	  ip = ipstart
	  do while (nbww(ip) .ne. 0)
		! loop over the interactions with other water molecules

		! X-O
		j    = nbww(ip)
		j3   = j*3-3
		iacj = iac(j)
		iLJ  = ljcod(iac(i),iac(j))
		crg(j)   = crg(j)
		aLJ  = iaclib(iaci)%avdw(iLJ)+iaclib(iacj)%avdw(iLJ)
		bLJ  = iaclib(iaci)%bvdw(iLJ)*iaclib(iacj)%bvdw(iLJ)
		aLJ  = aLJ*aLJ
		aLJ  = aLJ*aLJ*aLJ

		!distance between this oxygen atom and the oxygen atom of the above watermolecule, iw
		ds1 = x(j3+1) - x( 3*(nat_solute+3*(iw-1)+1) - 2 )
		ds2 = x(j3+2) - x( 3*(nat_solute+3*(iw-1)+1) - 1 )
		ds3 = x(j3+3) - x( 3*(nat_solute+3*(iw-1)+1)     )
		!the peridic shift
		ds1 = boxlength(1)*nint( ds1*inv_boxl(1) )
		ds2 = boxlength(2)*nint( ds2*inv_boxl(2) )
		ds3 = boxlength(3)*nint( ds3*inv_boxl(3) )


		dx1  = x(j3+1) - x(i3+1)
		dx2  = x(j3+2) - x(i3+2)
		dx3  = x(j3+3) - x(i3+3)
		dx1 = dx1 - ds1
		dx2 = dx2 - ds2
		dx3 = dx3 - ds3

		r2   = dx1*dx1 + dx2*dx2 + dx3*dx3
		r2   = 1./r2
		r    = sqrt ( r2 ) 
		r6   = r2*r2*r2
		r12  = r6*r6
		Vel  = crg(i)*crg(j)*r
		V_a  = bLJ*aLJ*aLJ*r12 
		V_b  = 2.0*bLJ*aLJ*r6
		dv   = r2*( -Vel -12.*V_a +6.*V_b )
		d(i3+1) = d(i3+1) - dv*dx1
		d(i3+2) = d(i3+2) - dv*dx2
		d(i3+3) = d(i3+3) - dv*dx3
		d(j3+1) = d(j3+1) + dv*dx1
		d(j3+2) = d(j3+2) + dv*dx2
		d(j3+3) = d(j3+3) + dv*dx3
		E%ww%el  = E%ww%el + Vel       
		E%ww%vdw = E%ww%vdw + V_a - V_b

		! X-H1
		j    = j + 1
		j3   = j3 + 3
		iacj = iac(j)
		iLJ  = ljcod(iac(i),iac(j))
		crg(j)   = crg(j)
		aLJ  = iaclib(iaci)%avdw(iLJ)+iaclib(iacj)%avdw(iLJ)
		bLJ  = iaclib(iaci)%bvdw(iLJ)*iaclib(iacj)%bvdw(iLJ)
		aLJ  = aLJ*aLJ
		aLJ  = aLJ*aLJ*aLJ
		dx1  = x(j3+1) - x(i3+1)
		dx2  = x(j3+2) - x(i3+2)
		dx3  = x(j3+3) - x(i3+3)
		dx1 = dx1 - ds1
		dx2 = dx2 - ds2
		dx3 = dx3 - ds3
		
		r2   = dx1*dx1 + dx2*dx2 + dx3*dx3
		r2   = 1./r2
		r    = sqrt ( r2 ) 
		r6   = r2*r2*r2
		r12  = r6*r6
		Vel  = crg(i)*crg(j)*r
		V_a  = bLJ*aLJ*aLJ*r12 
		V_b  = 2.0*bLJ*aLJ*r6
		dv   = r2*( -Vel -12.*V_a +6.*V_b )
		d(i3+1) = d(i3+1) - dv*dx1
		d(i3+2) = d(i3+2) - dv*dx2
		d(i3+3) = d(i3+3) - dv*dx3
		d(j3+1) = d(j3+1) + dv*dx1
		d(j3+2) = d(j3+2) + dv*dx2
		d(j3+3) = d(j3+3) + dv*dx3
		E%ww%el  = E%ww%el + Vel       
		E%ww%vdw = E%ww%vdw + V_a - V_b
		
		! X-H2
		j    = j + 1
		j3   = j3 + 3
		iacj = iac(j)
		iLJ  = ljcod(iac(i),iac(j))
		crg(j)   = crg(j)
		aLJ  = iaclib(iaci)%avdw(iLJ)+iaclib(iacj)%avdw(iLJ)
		bLJ  = iaclib(iaci)%bvdw(iLJ)*iaclib(iacj)%bvdw(iLJ)
		aLJ  = aLJ*aLJ
		aLJ  = aLJ*aLJ*aLJ
		dx1  = x(j3+1) - x(i3+1)
		dx2  = x(j3+2) - x(i3+2)
		dx3  = x(j3+3) - x(i3+3)
		dx1 = dx1 - ds1
		dx2 = dx2 - ds2
		dx3 = dx3 - ds3
	
		r2   = dx1*dx1 + dx2*dx2 + dx3*dx3
		r2   = 1./r2
		r    = sqrt ( r2 ) 
		r6   = r2*r2*r2
		r12  = r6*r6
		Vel  = crg(i)*crg(j)*r
		V_a  = bLJ*aLJ*aLJ*r12 
		V_b  = 2.0*bLJ*aLJ*r6
		dv   = r2*( -Vel -12.*V_a +6.*V_b )
		d(i3+1) = d(i3+1) - dv*dx1
		d(i3+2) = d(i3+2) - dv*dx2
		d(i3+3) = d(i3+3) - dv*dx3
		d(j3+1) = d(j3+1) + dv*dx1
		d(j3+2) = d(j3+2) + dv*dx2
		d(j3+3) = d(j3+3) + dv*dx3
		E%ww%el  = E%ww%el + Vel       
		E%ww%vdw = E%ww%vdw + V_a - V_b

		ip = ip + 1
	  end do
	end do

	ipstart = ip +1						! skip over the 0
  end do

end subroutine nonbon2_ww_box

!-----------------------------------------------------------------------

subroutine nonbond_pp
! local variables
integer						:: ip
type(NB_TYPE), pointer		:: pa, pb
real(8)						:: dx1a,dx2a,dx3a,r2a,ra,r6a
real(8)						:: Vela,V_aa,V_ba,dva
real(8)						:: dx1b,dx2b,dx3b,r2b,rb,r6b
real(8)						:: Velb,V_ab,V_bb,dvb

! global variables used:

!  x, crg, el14_scale, iaclib, d, E, 


do ip = 1, nbpp_pair - 1, 2
! for every second pair (two parallel runs to improve performance):


! set up pointers
pa => nbpp(ip)
pb => nbpp(ip+1)

! calculate the distance r
dx1a  = x(pa%j*3-2) - x(pa%i*3-2)
dx1b  = x(pb%j*3-2) - x(pb%i*3-2)
dx2a  = x(pa%j*3-1) - x(pa%i*3-1)
dx2b  = x(pb%j*3-1) - x(pb%i*3-1)
dx3a  = x(pa%j*3-0) - x(pa%i*3-0)
dx3b  = x(pb%j*3-0) - x(pb%i*3-0)
r2a   = 1./(dx1a*dx1a + dx2a*dx2a + dx3a*dx3a)
r2b   = 1./(dx1b*dx1b + dx2b*dx2b + dx3b*dx3b)
ra    = sqrt ( r2a ) 
rb    = sqrt ( r2b ) 
r6a   = r2a*r2a*r2a
r6b   = r2b*r2b*r2b

! calculate Vel and dv
Vela  = crg(pa%i)*crg(pa%j)*ra
Velb  = crg(pb%i)*crg(pb%j)*rb
if ( pa%LJcod .eq. 3 ) then 
  Vela = Vela*el14_scale
end if
if ( pb%LJcod .eq. 3 ) then
 Velb = Velb*el14_scale
end if
V_aa = r6a*r6a*iaclib(iac(pa%i))%avdw(pa%LJcod) &
        *iaclib(iac(pa%j))%avdw(pa%LJcod)
V_ab = r6b*r6b*iaclib(iac(pb%i))%avdw(pb%LJcod) &
        *iaclib(iac(pb%j))%avdw(pb%LJcod)
V_ba = r6a*iaclib(iac(pa%i))%bvdw(pa%LJcod) &
        *iaclib(iac(pa%j))%bvdw(pa%LJcod)
V_bb = r6b*iaclib(iac(pb%i))%bvdw(pb%LJcod) &
        *iaclib(iac(pb%j))%bvdw(pb%LJcod)

dva   = r2a*( -Vela -12.*V_aa +6.*V_ba )
dvb   = r2b*( -Velb -12.*V_ab +6.*V_bb )

! update d
d(pa%i*3-2) = d(pa%i*3-2) - dva*dx1a
d(pb%i*3-2) = d(pb%i*3-2) - dvb*dx1b
d(pa%i*3-1) = d(pa%i*3-1) - dva*dx2a
d(pb%i*3-1) = d(pb%i*3-1) - dvb*dx2b
d(pa%i*3-0) = d(pa%i*3-0) - dva*dx3a
d(pb%i*3-0) = d(pb%i*3-0) - dvb*dx3b
d(pa%j*3-2) = d(pa%j*3-2) + dva*dx1a
d(pb%j*3-2) = d(pb%j*3-2) + dvb*dx1b
d(pa%j*3-1) = d(pa%j*3-1) + dva*dx2a
d(pb%j*3-1) = d(pb%j*3-1) + dvb*dx2b
d(pa%j*3-0) = d(pa%j*3-0) + dva*dx3a
d(pb%j*3-0) = d(pb%j*3-0) + dvb*dx3b

! update energies
E%pp%el  = E%pp%el + Vela + Velb
E%pp%vdw = E%pp%vdw + V_aa + V_ab - V_ba - V_bb
end do

if (ip .eq. nbpp_pair) then
! odd #pairs, handle the last pair
pa => nbpp(ip)

dx1a  = x(pa%j*3-2) - x(pa%i*3-2)
dx2a  = x(pa%j*3-1) - x(pa%i*3-1)
dx3a  = x(pa%j*3-0) - x(pa%i*3-0)
r2a   = 1./(dx1a*dx1a + dx2a*dx2a + dx3a*dx3a)
ra   = sqrt(r2a)
r6a   = r2a*r2a*r2a

Vela  = crg(pa%i)*crg(pa%j)*ra
if ( pa%LJcod .eq. 3 ) Vela = Vela*el14_scale
V_aa = r6a*r6a*iaclib(iac(pa%i))%avdw(pa%LJcod)*iaclib(iac(pa%j))%avdw(pa%LJcod)
V_ba = r6a*iaclib(iac(pa%i))%bvdw(pa%LJcod)*iaclib(iac(pa%j))%bvdw(pa%LJcod)
dva   = r2a*( -Vela -12.*V_aa +6.*V_ba )

d(pa%i*3-2) = d(pa%i*3-2) - dva*dx1a
d(pa%i*3-1) = d(pa%i*3-1) - dva*dx2a
d(pa%i*3-0) = d(pa%i*3-0) - dva*dx3a
d(pa%j*3-2) = d(pa%j*3-2) + dva*dx1a
d(pa%j*3-1) = d(pa%j*3-1) + dva*dx2a
d(pa%j*3-0) = d(pa%j*3-0) + dva*dx3a

E%pp%el  = E%pp%el + Vela 
E%pp%vdw = E%pp%vdw + V_aa - V_ba 
end if

end subroutine nonbond_pp

!-----------------------------------------------------------------------
!******PWadded 2001-10-23
subroutine nonbond_pp_box
  ! local variables
  integer						:: ip, group, ga, gb
  type(NB_TYPE), pointer		:: pa, pb
  real(8)						:: dx1a,dx2a,dx3a,r2a,ra,r6a
  real(8)						:: Vela,V_aa,V_ba,dva
  real(8)						:: dx1b,dx2b,dx3b,r2b,rb,r6b
  real(8)						:: Velb,V_ab,V_bb,dvb

  ! global variables used:
  !  x, crg, el14_scale, iaclib, d, E, 

  !compute the peridocal shift for every charge group pair
	do group = 1, nbpp_cgp_pair
		ga = nbpp_cgp(group)%i !atom index for the two switching atoms
		gb = nbpp_cgp(group)%j

		!the distance between the two switching atoms
		dx1a = x(3*gb-2) - x(3*ga-2)
		dx2a = x(3*gb-1) - x(3*ga-1)
		dx3a = x(3*gb  ) - x(3*ga  )

		nbpp_cgp(group)%x = boxlength(1)*nint( dx1a*inv_boxl(1) )
		nbpp_cgp(group)%y = boxlength(2)*nint( dx2a*inv_boxl(2) )	
		nbpp_cgp(group)%z = boxlength(3)*nint( dx3a*inv_boxl(3) )

	end do

   do ip = 1, nbpp_pair - 1, 2
    ! for every second pair (two parallel runs to improve performance):

	! set up pointers
	pa => nbpp(ip)
	pb => nbpp(ip+1)
	ga = pa%cgp_pair
	gb = pb%cgp_pair

	! calculate the distance r
    dx1a  = x(pa%j*3-2) - x(pa%i*3-2)
    dx1b  = x(pb%j*3-2) - x(pb%i*3-2)
    dx2a  = x(pa%j*3-1) - x(pa%i*3-1)
    dx2b  = x(pb%j*3-1) - x(pb%i*3-1)
    dx3a  = x(pa%j*3-0) - x(pa%i*3-0)
    dx3b  = x(pb%j*3-0) - x(pb%i*3-0)
	dx1a = dx1a - nbpp_cgp(ga)%x         
	dx1b = dx1b - nbpp_cgp(gb)%x 
	dx2a = dx2a - nbpp_cgp(ga)%y    
	dx2b = dx2b - nbpp_cgp(gb)%y                
	dx3a = dx3a - nbpp_cgp(ga)%z  
	dx3b = dx3b - nbpp_cgp(gb)%z   
	
    r2a   = 1./(dx1a*dx1a + dx2a*dx2a + dx3a*dx3a)
    r2b   = 1./(dx1b*dx1b + dx2b*dx2b + dx3b*dx3b)
    ra    = sqrt ( r2a ) 
    rb    = sqrt ( r2b ) 
    r6a   = r2a*r2a*r2a
    r6b   = r2b*r2b*r2b

	! calculate Vel and dv
    Vela  = crg(pa%i)*crg(pa%j)*ra
    Velb  = crg(pb%i)*crg(pb%j)*rb
    if ( pa%LJcod .eq. 3 ) then 
	  Vela = Vela*el14_scale
	end if
	if ( pb%LJcod .eq. 3 ) then
	  Velb = Velb*el14_scale
	end if
	V_aa = r6a*r6a*iaclib(iac(pa%i))%avdw(pa%LJcod) &
		*iaclib(iac(pa%j))%avdw(pa%LJcod)
	V_ab = r6b*r6b*iaclib(iac(pb%i))%avdw(pb%LJcod) &
		*iaclib(iac(pb%j))%avdw(pb%LJcod)
	V_ba = r6a*iaclib(iac(pa%i))%bvdw(pa%LJcod) &
		*iaclib(iac(pa%j))%bvdw(pa%LJcod)
	V_bb = r6b*iaclib(iac(pb%i))%bvdw(pb%LJcod) &
		*iaclib(iac(pb%j))%bvdw(pb%LJcod)
    dva   = r2a*( -Vela -12.*V_aa +6.*V_ba )
    dvb   = r2b*( -Velb -12.*V_ab +6.*V_bb )

	! update d
    d(pa%i*3-2) = d(pa%i*3-2) - dva*dx1a
    d(pb%i*3-2) = d(pb%i*3-2) - dvb*dx1b
    d(pa%i*3-1) = d(pa%i*3-1) - dva*dx2a
    d(pb%i*3-1) = d(pb%i*3-1) - dvb*dx2b
    d(pa%i*3-0) = d(pa%i*3-0) - dva*dx3a
    d(pb%i*3-0) = d(pb%i*3-0) - dvb*dx3b
    d(pa%j*3-2) = d(pa%j*3-2) + dva*dx1a
    d(pb%j*3-2) = d(pb%j*3-2) + dvb*dx1b
    d(pa%j*3-1) = d(pa%j*3-1) + dva*dx2a
    d(pb%j*3-1) = d(pb%j*3-1) + dvb*dx2b
    d(pa%j*3-0) = d(pa%j*3-0) + dva*dx3a
    d(pb%j*3-0) = d(pb%j*3-0) + dvb*dx3b

	! update energies
    E%pp%el  = E%pp%el + Vela + Velb
    E%pp%vdw = E%pp%vdw + V_aa + V_ab - V_ba - V_bb
  end do

  if (ip .eq. nbpp_pair) then
    ! odd #pairs, handle the last pair
	pa => nbpp(ip)
	ga = pa%cgp_pair

	dx1a  = x(pa%j*3-2) - x(pa%i*3-2)
	dx2a  = x(pa%j*3-1) - x(pa%i*3-1)
	dx3a  = x(pa%j*3-0) - x(pa%i*3-0)
	dx1a = dx1a - nbpp_cgp(ga)%x  
	dx2a = dx2a - nbpp_cgp(ga)%y  
	dx3a = dx3a - nbpp_cgp(ga)%z  
	
	r2a   = 1./(dx1a*dx1a + dx2a*dx2a + dx3a*dx3a)
	ra   = sqrt(r2a)


	r6a   = r2a*r2a*r2a

	Vela  = crg(pa%i)*crg(pa%j)*ra
	if ( pa%LJcod .eq. 3 ) Vela = Vela*el14_scale
	V_aa = r6a*r6a*iaclib(iac(pa%i))%avdw(pa%LJcod)*iaclib(iac(pa%j))%avdw(pa%LJcod)
	V_ba = r6a*iaclib(iac(pa%i))%bvdw(pa%LJcod)*iaclib(iac(pa%j))%bvdw(pa%LJcod)
	dva   = r2a*( -Vela -12.*V_aa +6.*V_ba )

	d(pa%i*3-2) = d(pa%i*3-2) - dva*dx1a
	d(pa%i*3-1) = d(pa%i*3-1) - dva*dx2a
	d(pa%i*3-0) = d(pa%i*3-0) - dva*dx3a
	d(pa%j*3-2) = d(pa%j*3-2) + dva*dx1a
	d(pa%j*3-1) = d(pa%j*3-1) + dva*dx2a
	d(pa%j*3-0) = d(pa%j*3-0) + dva*dx3a

	E%pp%el  = E%pp%el + Vela 
	E%pp%vdw = E%pp%vdw + V_aa - V_ba 
  end if

end subroutine nonbond_pp_box

!-----------------------------------------------------------------------

subroutine nonbond_pw
! local variables
integer						:: ip,i,j,i3,j3,iaci,iacj,iLJ
real(8)						:: aLJ,bLJ,dx1,dx2,dx3,r2,r,r6,r12
real(8)						:: Vel,V_a,V_b,dv

! global variables used:
!  iac, crg, iaclib, x, d, E

do ip = 1, nbpw_pair
! for every assigned pair:

i    = nbpw(ip)%i
j    = nbpw(ip)%j
i3   = i*3-3
j3   = j*3-3
iaci = iac(i)
iacj = iac(j)

iLJ  = nbpw(ip)%LJcod
aLJ  = iaclib(iaci)%avdw(iLJ)*iaclib(iacj)%avdw(iLJ)
bLJ  = iaclib(iaci)%bvdw(iLJ)*iaclib(iacj)%bvdw(iLJ)

! calculate dx and r
dx1  = x(j3+1) - x(i3+1)
dx2  = x(j3+2) - x(i3+2)
dx3  = x(j3+3) - x(i3+3)
r2   = dx1*dx1 + dx2*dx2 + dx3*dx3
r2   = 1./r2
r    = sqrt ( r2 )
r6   = r2*r2*r2
r12  = r6*r6

! calculate Vel and dv
Vel  = crg(i)*crg(j)*r
V_a  = aLJ*r12 
V_b  = bLJ*r6
dv   = r2*( -Vel -12.*V_a +6.*V_b )

! update forces
d(i3+1) = d(i3+1) - dv*dx1
d(i3+2) = d(i3+2) - dv*dx2
d(i3+3) = d(i3+3) - dv*dx3
d(j3+1) = d(j3+1) + dv*dx1
d(j3+2) = d(j3+2) + dv*dx2
d(j3+3) = d(j3+3) + dv*dx3

! update energies
E%pw%el  = E%pw%el + Vel       
E%pw%vdw = E%pw%vdw + V_a - V_b
end do
end subroutine nonbond_pw

!-----------------------------------------------------------------------
!******PWadded 2001-10-23
subroutine nonbond_pw_box
  ! local variables
  integer						:: ip,i,j,i3,j3,iaci,iacj,iLJ
  real(8)						:: aLJ,bLJ,dx1,dx2,dx3,r2,r,r6,r12
  real(8)						:: Vel,V_a,V_b,dv
  integer						:: group, ga, gb

  ! global variables used:
  !  iac, crg, iaclib, x, d, E

  !compute the peridocal shift for every charge group pair
	do group = 1, nbpw_cgp_pair
		ga = nbpw_cgp(group)%i !atom index for the solute switching atoms
		gb = nbpw_cgp(group)%j !atom index for the solvent switching atom 

		!the distance between the two switching atoms
		dx1 = x(3*gb-2) - x(3*ga-2)
		dx2 = x(3*gb-1) - x(3*ga-1)
		dx3 = x(3*gb  ) - x(3*ga  )

		nbpw_cgp(group)%x = boxlength(1)*nint( dx1*inv_boxl(1) )
		nbpw_cgp(group)%y = boxlength(2)*nint( dx2*inv_boxl(2) )	
		nbpw_cgp(group)%z = boxlength(3)*nint( dx3*inv_boxl(3) )

	end do


  do ip = 1, nbpw_pair
	! for every assigned pair:

	i    = nbpw(ip)%i  ! solute atom 
	j    = nbpw(ip)%j	! solvent atom 
	group = nbpw(ip)%cgp_pair 
	i3   = i*3-3
	j3   = j*3-3
	iaci = iac(i)
	iacj = iac(j)
	iLJ  = nbpw(ip)%LJcod
	aLJ  = iaclib(iaci)%avdw(iLJ)*iaclib(iacj)%avdw(iLJ)
	bLJ  = iaclib(iaci)%bvdw(iLJ)*iaclib(iacj)%bvdw(iLJ)

	! calculate dx and r
	dx1  = x(j3+1) - x(i3+1)
	dx2  = x(j3+2) - x(i3+2)
	dx3  = x(j3+3) - x(i3+3)
	dx1 = dx1 - nbpw_cgp(group)%x  
	dx2 = dx2 - nbpw_cgp(group)%y  
	dx3 = dx3 - nbpw_cgp(group)%z 
	
	r2   = dx1*dx1 + dx2*dx2 + dx3*dx3
	r2   = 1./r2
	r    = sqrt ( r2 )
	r6   = r2*r2*r2
	r12  = r6*r6

	! calculate Vel and dv
	Vel  = crg(i)*crg(j)*r
	V_a  = aLJ*r12 
	V_b  = bLJ*r6
	dv   = r2*( -Vel -12.*V_a +6.*V_b )

	! update forces
	d(i3+1) = d(i3+1) - dv*dx1
	d(i3+2) = d(i3+2) - dv*dx2
	d(i3+3) = d(i3+3) - dv*dx3
	d(j3+1) = d(j3+1) + dv*dx1
	d(j3+2) = d(j3+2) + dv*dx2
	d(j3+3) = d(j3+3) + dv*dx3

	! update energies
	E%pw%el  = E%pw%el + Vel       
	E%pw%vdw = E%pw%vdw + V_a - V_b
  end do
end subroutine nonbond_pw_box

!-----------------------------------------------------------------------

subroutine nonbond_qq
! local variables
integer					:: istate, is,jj
integer					:: ip,iq,jq,i,j,k,i3,j3,iaci,iacj,iLJ
real(8)					:: qi,qj,aLJ,bLJ,dx1,dx2,dx3,r2,r,r6,r12,r6_hc
real(8)					:: Vel,V_a,V_b,dv,el_scale

do istate = 1, nstates
! for every state:

do ip = 1, nbqq_pair(istate)
  ! for every pair:

  iq   = nbqq(ip,istate)%iq   !q-atom number
  i    = iqseq(iq)           !atom number
  j    = nbqq(ip,istate)%j   !atom number
  jq   = nbqq(ip,istate)%jq  !q-atom number (if any)
  i3   = i*3-3
  j3   = j*3-3
  iLJ  = nbqq(ip,istate)%LJcod
  qi   = qcrg(iq,istate)
  el_scale = nbqq(ip,istate)%el_scale

  if (.not. qvdw_flag) then
        iaci = iac(i)
        iacj = iac(j)
        aLJ  = iaclib(iaci)%avdw(iLJ)*iaclib(iacj)%avdw(iLJ)
        bLJ  = iaclib(iaci)%bvdw(iLJ)*iaclib(iacj)%bvdw(iLJ)

        if (jq /= 0) then
          qj = qcrg(jq,istate)
        else
          qj = crg(j)
        end if
  else
        iaci = qiac(iq,istate)
        aLJ  = qavdw(iaci,iLJ)
        bLJ  = qbvdw(iaci,iLJ)
        if (jq /= 0) then
          iacj = qiac(jq,istate)
          aLJ = aLJ*qavdw(iacj,iLJ)
          bLJ = bLJ*qbvdw(iacj,iLJ)
          qj = qcrg(jq,istate)
        else
          if ( iLJ .eq. 2 ) then
                aLJ  = qavdw(iaci,1)
                bLJ  = qbvdw(iaci,1)
          end if
          iacj = iac(j)
          aLJ = aLJ*iaclib(iacj)%avdw(iLJ)
          bLJ = bLJ*iaclib(iacj)%bvdw(iLJ)
          qj = crg(j)
        end if


  end if




  ! calculate dx and r
  dx1  = x(j3+1) - x(i3+1)
  dx2  = x(j3+2) - x(i3+2)
  dx3  = x(j3+3) - x(i3+3)
  r2   = dx1*dx1 + dx2*dx2 + dx3*dx3

  r6_hc   = r2*r2*r2   !hardcore
  r6   = r2*r2*r2+sc_lookup(iq,natyps+jq,istate)   !Use softcore instead. sc is 0 for hardcore MPA
  r6   = 1./r6
  r12  = r6*r6

  r2   = 1./r2
  r    = sqrt ( r2 ) 

  ! calculate Vel, V_a, V_b and dv
  Vel  = qi*qj*r*el_scale
  if ( iLJ .eq. 3 ) Vel = Vel*el14_scale
  if (qvdw_flag .and. jq /= 0 .and. iLJ .eq. 2 ) then
        V_a = aLJ*exp(-bLJ/r)
        V_b = 0.0
        dv  = r2*( -Vel -bLJ*V_a/r )*EQ(istate)%lambda
  else
        V_a = aLJ*r12 
        V_b = bLJ*r6
        dv  = r2*( -Vel - (12.*V_a - 6.*V_b)*r6_hc*r6 )*EQ(istate)%lambda
  endif

  ! update forces
  d(i3+1) = d(i3+1) - dv*dx1
  d(i3+2) = d(i3+2) - dv*dx2
  d(i3+3) = d(i3+3) - dv*dx3
  d(j3+1) = d(j3+1) + dv*dx1
  d(j3+2) = d(j3+2) + dv*dx2
  d(j3+3) = d(j3+3) + dv*dx3

  ! update energies
  if ( jq /= 0 ) then
        EQ(istate)%qq%el  = EQ(istate)%qq%el + Vel
        EQ(istate)%qq%vdw = EQ(istate)%qq%vdw + V_a - V_b
	if (use_excluded_groups) then
                do jj = 1, ngroups_gc
                call set_gc_energies(i,j,Vel,(V_a-V_b),EQ_gc(jj,istate)%qq%el &
                ,EQ_gc(jj,istate)%qq%vdw,ST_gc(jj)%gcmask%mask)
                end do
	end if
  else
        EQ(istate)%qp%el  = EQ(istate)%qp%el + Vel
        EQ(istate)%qp%vdw = EQ(istate)%qp%vdw + V_a - V_b
        if (use_excluded_groups) then
                do jj = 1, ngroups_gc
                call set_gc_energies(i,j,Vel,(V_a-V_b),EQ_gc(jj,istate)%qp%el &
                ,EQ_gc(jj,istate)%qp%vdw,ST_gc(jj)%gcmask%mask)
                end do
        end if
  end if
end do
end do
end subroutine nonbond_qq

!-----------------------------------------------------------------------

subroutine nonbond_qq_lib_charges
!special version that uses library charges - for transformation of
!solute-surrounding only

! local variables
integer					:: istate,jj
integer					:: ip,iq,jq,i,j,k,i3,j3,iaci,iacj,iLJ
real(8)					:: qi,qj,aLJ,bLJ,dx1,dx2,dx3,r2,r,r6,r12,r6_hc
real(8)					:: Vel,V_a,V_b,dv,el_scale

do istate = 1, nstates
! for every state:

do ip = 1, nbqq_pair(istate)
  ! for every pair:

  iq   = nbqq(ip,istate)%iq
  i    = iqseq(iq)
  j    = nbqq(ip,istate)%j
  jq   = nbqq(ip,istate)%jq
  i3   = i*3-3
  j3   = j*3-3
  iLJ  = nbqq(ip,istate)%LJcod

  if (jq /= 0) then
	  qi   = crg(i)   !Use library charges for i
  else
	  qi   = qcrg(iq,istate)   !Since j is not a qatom we need the FEP file charges for i
  end if

  el_scale = nbqq(ip,istate)%el_scale

  if (.not. qvdw_flag) then
        iaci = iac(i)
        iacj = iac(j)
        aLJ  = iaclib(iaci)%avdw(iLJ)*iaclib(iacj)%avdw(iLJ)
        bLJ  = iaclib(iaci)%bvdw(iLJ)*iaclib(iacj)%bvdw(iLJ)
        qj = crg(j)
  else
        iaci = qiac(iq,istate)
        aLJ  = qavdw(iaci,iLJ)
        bLJ  = qbvdw(iaci,iLJ)
        if (jq /= 0) then
          iacj = qiac(jq,istate)
          aLJ = aLJ*qavdw(iacj,iLJ)
          bLJ = bLJ*qbvdw(iacj,iLJ)
        else
          if ( iLJ .eq. 2 ) then
                aLJ  = qavdw(iaci,1)
                bLJ  = qbvdw(iaci,1)
          end if
          iacj = iac(j)
          aLJ = aLJ*iaclib(iacj)%avdw(iLJ)
          bLJ = bLJ*iaclib(iacj)%bvdw(iLJ)
        end if
    qj = crg(j)
  end if

  ! calculate dx and r
  dx1  = x(j3+1) - x(i3+1)
  dx2  = x(j3+2) - x(i3+2)
  dx3  = x(j3+3) - x(i3+3)
  r2   = dx1*dx1 + dx2*dx2 + dx3*dx3

  r6_hc = r2*r2*r2   !hardcore
  r6   = r2*r2*r2+sc_lookup(iq,natyps+jq,istate)   !sc_lookup is softcore fix MPA
  r6   = 1./r6
  r12  = r6*r6

  r2   = 1./r2
  r    = sqrt ( r2 ) 

  ! calculate Vel, V_a, V_b and dv
  Vel  = qi*qj*r*el_scale
  if ( iLJ .eq. 3 ) Vel = Vel*el14_scale
  if (qvdw_flag .and. jq /= 0 .and. iLJ .eq. 2 ) then
        V_a = aLJ*exp(-bLJ/r)
        V_b = 0.0
        dv  = r2*( -Vel -bLJ*V_a/r )*EQ(istate)%lambda
  else
        V_a = aLJ*r12 
        V_b = bLJ*r6
        dv  = r2*( -Vel -(12.*V_a -6.*V_b)*r6_hc*r6 )*EQ(istate)%lambda
  endif

  ! update forces
  d(i3+1) = d(i3+1) - dv*dx1
  d(i3+2) = d(i3+2) - dv*dx2
  d(i3+3) = d(i3+3) - dv*dx3
  d(j3+1) = d(j3+1) + dv*dx1
  d(j3+2) = d(j3+2) + dv*dx2
  d(j3+3) = d(j3+3) + dv*dx3

  ! update energies
  if ( jq /= 0 ) then
        EQ(istate)%qq%el  = EQ(istate)%qq%el + Vel
        EQ(istate)%qq%vdw = EQ(istate)%qq%vdw + V_a - V_b
	if (use_excluded_groups) then
                do jj = 1, ngroups_gc
                call set_gc_energies(i,j,Vel,(V_a-V_b),EQ_gc(jj,istate)%qq%el &
                ,EQ_gc(jj,istate)%qq%vdw,ST_gc(jj)%gcmask%mask)
                end do
	end if
  else  ! j is not a qatom
        EQ(istate)%qp%el  = EQ(istate)%qp%el + Vel
        EQ(istate)%qp%vdw = EQ(istate)%qp%vdw + V_a - V_b
        if (use_excluded_groups) then
                do jj = 1, ngroups_gc
                call set_gc_energies(i,j,Vel,(V_a-V_b),EQ_gc(jj,istate)%qp%el &
                ,EQ_gc(jj,istate)%qp%vdw,ST_gc(jj)%gcmask%mask)
                end do
        end if
  end if
end do
end do
end subroutine nonbond_qq_lib_charges


!-----------------------------------------------------------------------

subroutine nonbond_qp

!calculate non-bonded interactions between q-atom i and  non-q-atoms j
!using standard atom types

! local variables
integer						:: ip,iq,i,j,i3,j3,iaci,iacj,iLJ
integer						:: istate,jj
real(8)						:: aLJ,bLJ,dx1,dx2,dx3,r2,r,r6
real(8)						:: Vel,V_a,V_b,dv

! global variables used:
!  iqseq, iac, crg, x, nstates, qvdw_flag, iaclib, qiac, 
! qcrg, el14_scale, EQ, d, nat_solute

do ip = 1, nbqp_pair
! for every assigned q-s pair:

! init state-invariant variables:
iq   = nbqp(ip)%i
i    = iqseq(iq)
j    = nbqp(ip)%j
i3   = i*3-3
j3   = j*3-3
iacj = iac(j)
iLJ  = nbqp(ip)%LJcod
dx1  = x(j3+1) - x(i3+1)
dx2  = x(j3+2) - x(i3+2)
dx3  = x(j3+3) - x(i3+3)
r2   = dx1*dx1 + dx2*dx2 + dx3*dx3
!softcore not needed here since nonbond_qp_qvdw is called instead
r2   = 1./r2
r6   = r2*r2*r2
r    = sqrt(r2)
do istate = 1, nstates
  ! for every state:

  ! calculate iaci, aLJ and bLJ
        iaci = iac(i)
        aLJ  = iaclib(iaci)%avdw(iLJ)*iaclib(iacj)%avdw(iLJ)
        bLJ  = iaclib(iaci)%bvdw(iLJ)*iaclib(iacj)%bvdw(iLJ)

  ! calculate qi, Vel, V_a, V_b and dv
  Vel  = qcrg(iq,istate)*crg(j)*r
  if ( iLJ .eq. 3 ) Vel = Vel*el14_scale


  V_a  = aLJ*r6*r6 
  V_b  = bLJ*r6        
  dv   = r2*( -Vel -12.*V_a +6.*V_b )*EQ(istate)%lambda

  ! update forces
  d(i3+1) = d(i3+1) - dv*dx1
  d(i3+2) = d(i3+2) - dv*dx2
  d(i3+3) = d(i3+3) - dv*dx3
  d(j3+1) = d(j3+1) + dv*dx1
  d(j3+2) = d(j3+2) + dv*dx2
  d(j3+3) = d(j3+3) + dv*dx3

  ! update q-protein energies
        EQ(istate)%qp%el  = EQ(istate)%qp%el + Vel
        EQ(istate)%qp%vdw = EQ(istate)%qp%vdw + V_a - V_b
        if (use_excluded_groups) then
                do jj = 1, ngroups_gc
                call set_gc_energies(i,j,Vel,(V_a-V_b),EQ_gc(jj,istate)%qp%el &
                ,EQ_gc(jj,istate)%qp%vdw,ST_gc(jj)%gcmask%mask)
                end do
        end if
end do ! istate

end do
end subroutine nonbond_qp

!-----------------------------------------------------------------------
!******PWadded 2001-10-23
subroutine nonbond_qp_box

	!calculate non-bonded interactions between q-atom i and  non-q-atoms j
	!using standard atom types

  ! local variables
  integer						:: ip,iq,i,j,i3,j3,iaci,iacj,iLJ
  integer						:: istate,jj
  real(8)						:: aLJ,bLJ,dx1,dx2,dx3,r2,r,r6
  real(8)						:: Vel,V_a,V_b,dv
  integer						:: group, gr, ia
  ! global variables used:
  !  iqseq, iac, crg, x, nstates, qvdw_flag, iaclib, qiac, qavdw, qbvdw, 
  ! qcrg, el14_scale, EQ, d, nat_solute


	!compute the peridocal shift for every charge group pair
	do gr = 1, nbqp_cgp_pair
		ia = nbqp_cgp(gr)%i !atom index for the switching atom
		
		!the distance between the two switching atoms
		dx1 = x(3*ia-2) - x(3*qswitch-2)
		dx2 = x(3*ia-1) - x(3*qswitch-1)
		dx3 = x(3*ia  ) - x(3*qswitch  )

		nbqp_cgp(gr)%x = boxlength(1)*nint( dx1*inv_boxl(1) )
		nbqp_cgp(gr)%y = boxlength(2)*nint( dx2*inv_boxl(2) )	
		nbqp_cgp(gr)%z = boxlength(3)*nint( dx3*inv_boxl(3) )

	end do


  do ip = 1, nbqp_pair
	! for every assigned q-s pair:
	! init state-invariant variables:
	iq   = nbqp(ip)%i
	i    = iqseq(iq)
	j    = nbqp(ip)%j
	i3   = i*3-3
	j3   = j*3-3
	iacj = iac(j)
	iLJ  = nbqp(ip)%LJcod
	group = nbqp(ip)%cgp_pair
	dx1  = x(j3+1) - x(i3+1)
	dx2  = x(j3+2) - x(i3+2)
	dx3  = x(j3+3) - x(i3+3)
	dx1 = dx1 - nbqp_cgp(group)%x  
	dx2 = dx2 - nbqp_cgp(group)%y  
	dx3 = dx3 - nbqp_cgp(group)%z  
	
	r2   = dx1*dx1 + dx2*dx2 + dx3*dx3
	r2   = 1./r2
	r6   = r2*r2*r2   !softcore not needed here. taken care of in nonbond_qp_qvdw_box
	r    = sqrt(r2)

	do istate = 1, nstates
	  ! for every state:

	  ! calculate iaci, aLJ and bLJ
		iaci = iac(i)
		aLJ  = iaclib(iaci)%avdw(iLJ)*iaclib(iacj)%avdw(iLJ)
		bLJ  = iaclib(iaci)%bvdw(iLJ)*iaclib(iacj)%bvdw(iLJ)

	  ! calculate qi, Vel, V_a, V_b and dv
	  Vel  = qcrg(iq,istate)*crg(j)*r
	  if ( iLJ .eq. 3 ) Vel = Vel*el14_scale

	  V_a  = aLJ*r6*r6 
	  V_b  = bLJ*r6
	  dv   = r2*( -Vel -12.*V_a +6.*V_b )*EQ(istate)%lambda

	  ! update forces
	  d(i3+1) = d(i3+1) - dv*dx1
	  d(i3+2) = d(i3+2) - dv*dx2
	  d(i3+3) = d(i3+3) - dv*dx3
	  d(j3+1) = d(j3+1) + dv*dx1
	  d(j3+2) = d(j3+2) + dv*dx2
	  d(j3+3) = d(j3+3) + dv*dx3

	  ! update q-protein energies
		EQ(istate)%qp%el  = EQ(istate)%qp%el + Vel
		EQ(istate)%qp%vdw = EQ(istate)%qp%vdw + V_a - V_b
        if (use_excluded_groups) then
                do jj = 1, ngroups_gc
                call set_gc_energies(i,j,Vel,(V_a-V_b),EQ_gc(jj,istate)%qp%el &
                ,EQ_gc(jj,istate)%qp%vdw,ST_gc(jj)%gcmask%mask)
                end do
        end if
	end do ! istate

  end do
end subroutine nonbond_qp_box

!----------------------------------------------------------------------

subroutine nonbond_qp_qvdw
!calculate nonbonded interactions between q-atom i and non-q-atom j
!using q-atom types

! local variables
integer						:: ip,iq,i,j,i3,j3,iaci,iacj,iLJ,qLJ
integer						:: istate,jj
real(8)						:: aLJ,bLJ,dx1,dx2,dx3,r2,r,r6
real(8)						:: Vel,V_a,V_b,dv,r6_sc

! global variables used:
!  iqseq, iac, crg, x, nstates, qvdw_flag, iaclib, qiac, qavdw, qbvdw, qcrg, el14_scale, EQ, d, nat_solute


do ip = 1, nbqp_pair
! for every assigned q-s pair:

! init state-invariant variables:
iq   = nbqp(ip)%i
i    = iqseq(iq)
j    = nbqp(ip)%j
i3   = i*3-3
j3   = j*3-3
iacj = iac(j)
iLJ  = nbqp(ip)%LJcod
qLJ  = nbqp(ip)%qLJcod
dx1  = x(j3+1) - x(i3+1)
dx2  = x(j3+2) - x(i3+2)
dx3  = x(j3+3) - x(i3+3)
r2   = dx1*dx1 + dx2*dx2 + dx3*dx3

r6   = r2*r2*r2
r2   = 1./r2
r    = sqrt(r2)

do istate = 1, nstates
        ! for every state:

        ! calculate iaci, aLJ and bLJ
        iaci = qiac(iq,istate)
        aLJ  = qavdw(iaci,qLJ)*iaclib(iacj)%avdw(iLJ)
        bLJ  = qbvdw(iaci,qLJ)*iaclib(iacj)%bvdw(iLJ)

  ! calculate qi, Vel, V_a, V_b and dv
  Vel  = qcrg(iq,istate)*crg(j)*r
  if ( iLJ .eq. 3 ) Vel = Vel*el14_scale

  r6_sc = r6 + sc_lookup(iq,iacj,istate) !sc_lookup is softcore fix MPA
  V_a  = aLJ/(r6_sc*r6_sc)
  V_b  = bLJ/(r6_sc)
  dv   = r2*( -Vel -(12.*V_a -6.*V_b)*(r6/r6_sc) )*EQ(istate)%lambda  !r6 is r^6 not 1/r^6, r6_sc is r^6+sc not 1/(r^6+sc)

  ! update forces
  d(i3+1) = d(i3+1) - dv*dx1
  d(i3+2) = d(i3+2) - dv*dx2
  d(i3+3) = d(i3+3) - dv*dx3
  d(j3+1) = d(j3+1) + dv*dx1
  d(j3+2) = d(j3+2) + dv*dx2
  d(j3+3) = d(j3+3) + dv*dx3

  ! update q-protein energies
        EQ(istate)%qp%el  = EQ(istate)%qp%el + Vel
        EQ(istate)%qp%vdw = EQ(istate)%qp%vdw + V_a - V_b
	if (use_excluded_groups) then
                do jj = 1, ngroups_gc
                call set_gc_energies(i,j,Vel,(V_a-V_b),EQ_gc(jj,istate)%qp%el &
                ,EQ_gc(jj,istate)%qp%vdw,ST_gc(jj)%gcmask%mask)
                end do
	end if
end do ! istate

end do
end subroutine nonbond_qp_qvdw
!-----------------------------------------------------------------------
!******PWadded 2001-10-23
subroutine nonbond_qp_qvdw_box
	!calculate nonbonded interactions between q-atom i and non-q-atom j
	!using q-atom types

  ! local variables
  integer						:: ip,iq,i,j,i3,j3,iaci,iacj,iLJ,qLJ
  integer						:: istate,jj
  real(8)						:: aLJ,bLJ,dx1,dx2,dx3,r2,r,r6
  real(8)						:: Vel,V_a,V_b,dv,r6_sc
  integer						:: group, gr, ia


  ! global variables used:
  !  iqseq, iac, crg, x, nstates, qvdw_flag, iaclib, qiac, qavdw, qbvdw, qcrg, el14_scale, EQ, d, nat_solute


	!compute the peridocal shift for every charge group pair
	do gr = 1, nbqp_cgp_pair
		ia = nbqp_cgp(gr)%i !atom index for the atom
		
		!the distance between the two switching atoms
		dx1 = x(3*ia-2) - x(3*qswitch-2)
		dx2 = x(3*ia-1) - x(3*qswitch-1)
		dx3 = x(3*ia  ) - x(3*qswitch  )

		nbqp_cgp(gr)%x = boxlength(1)*nint( dx1*inv_boxl(1) )
		nbqp_cgp(gr)%y = boxlength(2)*nint( dx2*inv_boxl(2) )	
		nbqp_cgp(gr)%z = boxlength(3)*nint( dx3*inv_boxl(3) )

	end do


  do ip = 1, nbqp_pair
	! for every assigned q-s pair:

	! init state-invariant variables:
	iq   = nbqp(ip)%i
	i    = iqseq(iq)
	j    = nbqp(ip)%j
	i3   = i*3-3
	j3   = j*3-3
	iacj = iac(j)
	iLJ  = nbqp(ip)%LJcod
	qLJ  = nbqp(ip)%qLJcod
	group = nbqp(ip)%cgp_pair
	
	dx1  = x(j3+1) - x(i3+1)
	dx2  = x(j3+2) - x(i3+2)
	dx3  = x(j3+3) - x(i3+3)
	dx1 = dx1 - nbqp_cgp(group)%x  
	dx2 = dx2 - nbqp_cgp(group)%y  
	dx3 = dx3 - nbqp_cgp(group)%z  
	
	r2   = dx1*dx1 + dx2*dx2 + dx3*dx3

	r6   = r2*r2*r2
	r2   = 1./r2
	r    = sqrt(r2)

	do istate = 1, nstates
		! for every state:

		! calculate iaci, aLJ and bLJ
		iaci = qiac(iq,istate)
		aLJ  = qavdw(iaci,qLJ)*iaclib(iacj)%avdw(iLJ)
		bLJ  = qbvdw(iaci,qLJ)*iaclib(iacj)%bvdw(iLJ)

	  ! calculate qi, Vel, V_a, V_b and dv
	  Vel  = qcrg(iq,istate)*crg(j)*r
	  if ( iLJ .eq. 3 ) Vel = Vel*el14_scale


	  r6_sc = r6 + sc_lookup(iq,iacj,istate)        !sc_lookup is softcore fix MPA
	  V_a  = aLJ/(r6_sc*r6_sc)
	  V_b  = bLJ/r6_sc         !sc_lookup is softcore fix MPA
	  dv   = r2*( -Vel - ( (12.*V_a - 6.*V_b)*(r6/r6_sc) ) )*EQ(istate)%lambda  !r6 is r^6 not 1/r^6, r6_sc is r^6+sc not 1/(r^6+sc)

	  ! update forces
	  d(i3+1) = d(i3+1) - dv*dx1
	  d(i3+2) = d(i3+2) - dv*dx2
	  d(i3+3) = d(i3+3) - dv*dx3
	  d(j3+1) = d(j3+1) + dv*dx1
	  d(j3+2) = d(j3+2) + dv*dx2
	  d(j3+3) = d(j3+3) + dv*dx3

	  ! update q-protein energies
		EQ(istate)%qp%el  = EQ(istate)%qp%el + Vel
		EQ(istate)%qp%vdw = EQ(istate)%qp%vdw + V_a - V_b
        if (use_excluded_groups) then
                do jj = 1, ngroups_gc
                call set_gc_energies(i,j,Vel,(V_a-V_b),EQ_gc(jj,istate)%qp%el &
                ,EQ_gc(jj,istate)%qp%vdw,ST_gc(jj)%gcmask%mask)
                end do 
        end if
	end do ! istate

  end do
end subroutine nonbond_qp_qvdw_box

!-----------------------------------------------------------------------

subroutine nonbond_qw_spc
!calculate non-bonded interactions between Q-atoms and SPC water molecules
!(optimisations rely on LJ params = 0 for water H) using geometric comb. rule

! local variables
integer						:: jw,iq,i,j,iLJ
integer						:: istate
real(8)						::	aLJ, bLJ
real(8)						::	dxO, dyO, dzO, dxH1, dyH1, dzH1, dxH2, dyH2, dzH2
real(8)						::	rO, r2O, r6O, rH1, r2H1, r6H1, rH2, r2H2, r6H2
real(8)						::	VelO, VelH1, VelH2, dvO, dvH1, dvH2
real(8)						::  V_a, V_b, r6O_sc, r6O_hc
real(8), save				::	aO(2), bO(2)
integer, save				::	iac_ow = 0, iac_hw = 0

! global variables used:
!  iqseq, iac, crg, x, nstates, qvdw_flag, iaclib, qiac, qavdw, qbvdw, qcrg, el14_scale, EQ, d, nat_solute

if(iac_ow == 0) then !set first time
        iac_ow = iac(nat_solute + 1)
        iac_hw = iac(nat_solute + 2)
        aO(1:2) = iaclib(iac(nat_solute + 1))%avdw(1:2)
        bO(1:2) = iaclib(iac(nat_solute + 1))%bvdw(1:2)
end if

!loop over listed waters
do jw = 1, nbqw_pair
        j = nbqw(jw) !top. # of O in water iw
        !loop over all q-atoms
        do iq = 1, nqat
                i = iqseq(iq)
                dxO  = x(3*j-2) - x(3*i-2)
                dyO  = x(3*j-1) - x(3*i-1)
                dzO  = x(3*j  ) - x(3*i  )
                dxH1 = x(3*j+1) - x(3*i-2)
                dyH1 = x(3*j+2) - x(3*i-1)
                dzH1 = x(3*j+3) - x(3*i  )
                dxH2 = x(3*j+4) - x(3*i-2)
                dyH2 = x(3*j+5) - x(3*i-1)
                dzH2 = x(3*j+6) - x(3*i  )
                r2O  = dxO*dxO + dyO*dyO + dzO*dzO
                rH1  = sqrt(1._8/(dxH1*dxH1 + dyH1*dyH1 + dzH1*dzH1))
                rH2  = sqrt(1._8/(dxH2*dxH2 + dyH2*dyH2 + dzH2*dzH2))
                r6O_hc  = r2O*r2O*r2O   !will set r6O to 1/r6O later, needed for softcore
                r2O  = 1._8/r2O
                rO   = sqrt(r2O)
                r2H1 = rH1*rH1
                r6H1 = r2H1*r2H1*r2H1
                r2H2 = rH2*rH2
                r6H2 = r2H2*r2H2*r2H2

		        r6O_sc = r6O_hc   !default is hardcore (i.e. not softcore)  MPA

                !reset potential
                dvO = 0
                dvH1 = 0
                dvH2 = 0
                iLJ = LJcod(iac_ow, iac(i))
                if(.not. qvdw_flag) then
                        !use same LJ params for all states
						r6O  = 1./r6O_hc     !softcore hack, see comment 15 lines up
                        aLJ  = iaclib(iac(i))%avdw(iLJ)
                        bLJ  = iaclib(iac(i))%bvdw(iLJ)
                        V_a  = aLJ*aO(iLJ)*r6O*r6O 
                        V_b  = bLJ*bO(iLJ)*r6O
                end if
                do istate = 1, nstates
                        ! for every state:
                        ! calculate iaci, aLJ and bLJ
                        if (qvdw_flag) then
                                aLJ  = qavdw(qiac(iq,istate),1)
                                bLJ  = qbvdw(qiac(iq,istate),1)
				                r6O_sc = r6O_hc + sc_lookup(iq,iac_ow,istate)   !softcore  MPA

                                V_a  = aLJ*aO(iLJ)/(r6O_sc*r6O_sc)
                                V_b  = bLJ*bO(iLJ)/(r6O_sc)
                        end if
                        ! calculate qi, Vel, V_a, V_b and dv
                        VelO = crg_ow*qcrg(iq,istate)*rO
                        VelH1 = crg_hw*qcrg(iq,istate)*rH1
                        VelH2 = crg_hw*qcrg(iq,istate)*rH2
                        dvO  = dvO  + r2O*( -VelO -( (12.*V_a - 6.*V_b)*(r6O_hc/r6O_sc) ))*EQ(istate)%lambda
                        dvH1 = dvH1 - r2H1*VelH1*EQ(istate)%lambda
                        dvH2 = dvH2 - r2H2*VelH2*EQ(istate)%lambda
                        ! update q-water energies
                        EQ(istate)%qw%el  = EQ(istate)%qw%el + VelO + VelH1 + VelH2
                        EQ(istate)%qw%vdw = EQ(istate)%qw%vdw + V_a - V_b
                end do !istate
				
				! if qvdw_flag is true, then r6O is not the usual 1/rO^6, but rather rO^6. be careful!!! MPA
				
							
                ! update forces on Q-atom
                d(3*i-2) = d(3*i-2) - dvO*dxO - dvH1*dxH1 - dvH2*dxH2
                d(3*i-1) = d(3*i-1) - dvO*dyO - dvH1*dyH1 - dvH2*dyH2
                d(3*i  ) = d(3*i  ) - dvO*dzO - dvH1*dzH1 - dvH2*dzH2
                ! update forces on water
                d(3*j-2) = d(3*j-2) + dvO*dxO
                d(3*j-1) = d(3*j-1) + dvO*dyO
                d(3*j  ) = d(3*j  ) + dvO*dzO
                d(3*j+1) = d(3*j+1) + dvH1*dxH1
                d(3*j+2) = d(3*j+2) + dvH1*dyH1
                d(3*j+3) = d(3*j+3) + dvH1*dzH1
                d(3*j+4) = d(3*j+4) + dvH2*dxH2
                d(3*j+5) = d(3*j+5) + dvH2*dyH2
                d(3*j+6) = d(3*j+6) + dvH2*dzH2

        end do !iq
end do !jw
end subroutine nonbond_qw_spc

!-----------------------------------------------------------------------
!******PWadded 2001-10-23
subroutine nonbond_qw_spc_box
	!calculate non-bonded interactions between Q-atoms and SPC water molecules
	!(optimisations rely on LJ params = 0 for water H) using geometric comb. rule

	! local variables
	integer						:: jw,iq,i,j,iLJ
	integer						:: istate
	real(8)						::	aLJ, bLJ
	real(8)						::	dxO, dyO, dzO, dxH1, dyH1, dzH1, dxH2, dyH2, dzH2
	real(8)						::	rO, r2O, r6O, rH1, r2H1, r6H1, rH2, r2H2, r6H2
	real(8)						::	VelO, VelH1, VelH2, dvO, dvH1, dvH2
	real(8)						:: V_a, V_b, r6O_sc, r6O_hc
	real(8)						:: dx, dy, dz, boxshiftx, boxshifty, boxshiftz
	real(8), save				::	aO(2), bO(2)
	integer, save				::	iac_ow = 0, iac_hw = 0

	! global variables used:
	!  iqseq, iac, crg, x, nstates, qvdw_flag, iaclib, qiac, qavdw, qbvdw, qcrg, el14_scale, EQ, d, nat_solute

  	if(iac_ow == 0) then !set first time   
		iac_ow = iac(nat_solute + 1)
		iac_hw = iac(nat_solute + 2)    ! not used !?
		aO(1:2) = iaclib(iac(nat_solute + 1))%avdw(1:2)
		bO(1:2) = iaclib(iac(nat_solute + 1))%bvdw(1:2)
	end if

	!loop over listed waters
	do jw = 1, nbqw_pair
		j = nbqw(jw) !top. # of O in water iw
	
		!compute the periodical shift
		dx = x(3*j-2) - x(3*qswitch-2)
		dy = x(3*j-1) - x(3*qswitch-1)
		dz = x(3*j  ) - x(3*qswitch  )
		boxshiftx = boxlength(1)*nint(dx*inv_boxl(1))
		boxshifty = boxlength(2)*nint(dy*inv_boxl(2))
		boxshiftz = boxlength(3)*nint(dz*inv_boxl(3))	
	
		!loop over all q-atoms
		do iq = 1, nqat
			i = iqseq(iq)
			dxO  = x(3*j-2) - x(3*i-2)
			dyO  = x(3*j-1) - x(3*i-1)
			dzO  = x(3*j  ) - x(3*i  )
			dxH1 = x(3*j+1) - x(3*i-2)
			dyH1 = x(3*j+2) - x(3*i-1)
			dzH1 = x(3*j+3) - x(3*i  )
			dxH2 = x(3*j+4) - x(3*i-2)
			dyH2 = x(3*j+5) - x(3*i-1)
			dzH2 = x(3*j+6) - x(3*i  )
			dxO = dxO - boxshiftx  
			dyO = dyO - boxshifty  
			dzO = dzO - boxshiftz  
			dxH1 = dxH1 - boxshiftx  
			dyH1 = dyH1 - boxshifty  
			dzH1 = dzH1 - boxshiftz  
			dxH2 = dxH2 - boxshiftx  
			dyH2 = dyH2 - boxshifty  
			dzH2 = dzH2 - boxshiftz  

			r2O = dxO*dxO + dyO*dyO + dzO*dzO
			rH1 = sqrt(1._8/(dxH1*dxH1 + dyH1*dyH1 + dzH1*dzH1))
			rH2 = sqrt(1._8/(dxH2*dxH2 + dyH2*dyH2 + dzH2*dzH2))
			r6O_hc = r2O*r2O*r2O   !will set r6O = 1/r6O later, need for softcore MPA
			r2O = 1._8/r2O
			rO  = sqrt(r2O)
			r2H1 = rH1*rH1
			r6H1 = r2H1*r2H1*r2H1
			r2H2 = rH2*rH2
			r6H2 = r2H2*r2H2*r2H2

			r6O_sc = r6O_hc   !default is hardcore (i.e. not softcore)  MPA

			!reset potential
			dvO = 0
			dvH1 = 0
			dvH2 = 0
			iLJ = LJcod(iac_ow, iac(i))
			if(.not. qvdw_flag) then
				!use same LJ params for all states
				aLJ  = iaclib(iac(i))%avdw(iLJ)
				bLJ  = iaclib(iac(i))%bvdw(iLJ)
				r6O  = 1._8/r6O_hc   !softcore hack  MPA
				V_a  = aLJ*aO(iLJ)*r6O*r6O 
				V_b  = bLJ*bO(iLJ)*r6O
			end if
			do istate = 1, nstates
				! for every state:

				! calculate iaci, aLJ and bLJ
				if (qvdw_flag) then
					aLJ  = qavdw(qiac(iq,istate),1)
					bLJ  = qbvdw(qiac(iq,istate),1)
					r6O_sc = r6O_hc + sc_lookup(iq,iac_ow,istate)   !softcore  MPA

					V_a  = aLJ*aO(iLJ)/(r6O_sc*r6O_sc)
					V_b  = bLJ*bO(iLJ)/r6O_sc
				end if
				! calculate qi, Vel, V_a, V_b and dv
				VelO = crg_ow*qcrg(iq,istate)*rO
				VelH1 = crg_hw*qcrg(iq,istate)*rH1
				VelH2 = crg_hw*qcrg(iq,istate)*rH2
				dvO  = dvO  + r2O*( -VelO -( (12.*V_a - 6.*V_b)*(r6O_hc/r6O_sc) ))*EQ(istate)%lambda
				dvH1 = dvH1 - r2H1*VelH1*EQ(istate)%lambda
				dvH2 = dvH2 - r2H2*VelH2*EQ(istate)%lambda
				! update q-water energies
				EQ(istate)%qw%el  = EQ(istate)%qw%el + VelO + VelH1 + VelH2
				EQ(istate)%qw%vdw = EQ(istate)%qw%vdw + V_a - V_b
			end do !istate			
			! update forces on Q-atom
			d(3*i-2) = d(3*i-2) - dvO*dxO - dvH1*dxH1 - dvH2*dxH2
			d(3*i-1) = d(3*i-1) - dvO*dyO - dvH1*dyH1 - dvH2*dyH2
			d(3*i  ) = d(3*i  ) - dvO*dzO - dvH1*dzH1 - dvH2*dzH2
			! update forces on water
			d(3*j-2) = d(3*j-2) + dvO*dxO
			d(3*j-1) = d(3*j-1) + dvO*dyO
			d(3*j  ) = d(3*j  ) + dvO*dzO
			d(3*j+1) = d(3*j+1) + dvH1*dxH1
			d(3*j+2) = d(3*j+2) + dvH1*dyH1
			d(3*j+3) = d(3*j+3) + dvH1*dzH1
			d(3*j+4) = d(3*j+4) + dvH2*dxH2
			d(3*j+5) = d(3*j+5) + dvH2*dyH2
			d(3*j+6) = d(3*j+6) + dvH2*dzH2

		end do !iq
	end do !jw
end subroutine nonbond_qw_spc_box

!-----------------------------------------------------------------------

subroutine nonbond_qw_3atom
!calculate non-bonded interactions between Q-atoms and 3-atom solvent molecules
!using geometric comb. rule

! local variables
integer						:: jw,iq,i,j,iLJ1, iLJ2, iLJ3, iaci
integer						:: istate
real(8)						::	dx1, dy1, dz1, dx2, dy2, dz2, dx3, dy3, dz3
real(8)						::	r_1, r2_1, r6_1, r_2, r2_2, r6_2, r_3, r2_3, r6_3
real(8)						::	Vel1, Vel2, Vel3, dv1, dv2, dv3
real(8)						:: V_a1,V_b1, V_a2, V_b2, V_a3, V_b3,r6_1_sc,r6_2_sc,r6_3_sc,r6_1_hc,r6_2_hc,r6_3_hc
real(8), save				::	a1(2), b1(2), a2(2), b2(2), a3(2), b3(2)
integer, save				::	iac1, iac2, iac3
real, save					::	crg1, crg2, crg3
! global variables used:
!  iqseq, iac, crg, x, nstates, qvdw_flag, iaclib, qiac, qavdw, qbvdw, qcrg, el14_scale, EQ, d, nat_solute

if(a1(1) == 0.) then !set first time
        iac1 = iac(nat_solute+1)
        iac2 = iac(nat_solute+2)
        iac3 = iac(nat_solute+3)
        crg1 = crg(nat_solute+1)
        crg2 = crg(nat_solute+2)
        crg3 = crg(nat_solute+3)
        a1(1:2) = iaclib(iac1)%avdw(1:2)
        b1(1:2) = iaclib(iac1)%bvdw(1:2)
        a2(1:2) = iaclib(iac2)%avdw(1:2)
        b2(1:2) = iaclib(iac2)%bvdw(1:2)
        a3(1:2) = iaclib(iac3)%avdw(1:2)
        b3(1:2) = iaclib(iac3)%bvdw(1:2)
end if

!loop over listed waters
do jw = 1, nbqw_pair
        j = nbqw(jw) !top. # of O in water iw
        !loop over all q-atoms
        do iq = 1, nqat
                i = iqseq(iq)
                dx1 = x(3*j-2) - x(3*i-2)
                dy1 = x(3*j-1) - x(3*i-1)
                dz1 = x(3*j  ) - x(3*i  )
                dx2 = x(3*j+1) - x(3*i-2)
                dy2 = x(3*j+2) - x(3*i-1)
                dz2 = x(3*j+3) - x(3*i  )
                dx3 = x(3*j+4) - x(3*i-2)
                dy3 = x(3*j+5) - x(3*i-1)
                dz3 = x(3*j+6) - x(3*i  )
                r2_1 = dx1*dx1 + dy1*dy1 + dz1*dz1
                r2_2 = dx2*dx2 + dy2*dy2 + dz2*dz2
                r2_3 = dx3*dx3 + dy3*dy3 + dz3*dz3
                r6_1_hc = r2_1*r2_1*r2_1   !will set r6 to 1/r6 later, needed for softcore
                r2_1 = 1._8/r2_1
				r_1  = sqrt(r2_1)
                r6_2_hc = r2_2*r2_2*r2_2   !will set r6 to 1/r6 later, needed for softcore
                r2_2 = 1._8/r2_2
				r_2  = sqrt(r2_2)
                r6_3_hc = r2_3*r2_3*r2_3   !will set r6 to 1/r6 later, needed for softcore
                r2_3 = 1._8/r2_3
				r_3  = sqrt(r2_3)

                !reset potential
                dv1 = 0
                dv2 = 0
                dv3 = 0
                iaci = iac(i)
                iLJ1 = LJcod(iac1, iaci)
                iLJ2 = LJcod(iac2, iaci)
                iLJ3 = LJcod(iac3, iaci)
                if(.not. qvdw_flag) then
                        !use same LJ params for all states
						r6_1 = 1._8/r6_1_hc	!softcore hack
						r6_2 = 1._8/r6_2_hc	!softcore hack
						r6_3 = 1._8/r6_3_hc	!softcore hack

                        V_a1  = iaclib(iaci)%avdw(iLJ1)*a1(iLJ1)*r6_1*r6_1 
                        V_b1  = iaclib(iaci)%bvdw(iLJ1)*b1(iLJ1)*r6_1
                        V_a2 = iaclib(iaci)%avdw(iLJ2)*a2(iLJ2)*r6_2*r6_2
                        V_b2 = iaclib(iaci)%bvdw(iLJ2)*b2(iLJ2)*r6_2
                        V_a3 = iaclib(iaci)%avdw(iLJ3)*a3(iLJ3)*r6_3*r6_3
                        V_b3 = iaclib(iaci)%bvdw(iLJ3)*b3(iLJ3)*r6_3
                end if
                do istate = 1, nstates
                        ! for every state:

                        ! calculate V_a:s and V_b:s for each state
                        if (qvdw_flag) then
				r6_1_sc = r6_1_hc + sc_lookup(iq,iac1,istate)
				r6_2_sc = r6_2_hc + sc_lookup(iq,iac2,istate)
				r6_3_sc = r6_3_hc + sc_lookup(iq,iac3,istate)
				V_a1 = qavdw(qiac(iq,istate),1)*a1(iLJ1)/(r6_1_sc*r6_1_sc)
                                V_b1 = qbvdw(qiac(iq,istate),1)*b1(iLJ1)/(r6_1_sc)
                                V_a2 = qavdw(qiac(iq,istate),1)*a2(iLJ2)/(r6_2_sc*r6_2_sc)
                                V_b2 = qbvdw(qiac(iq,istate),1)*b2(iLJ2)/(r6_2_sc)
                                V_a3 = qavdw(qiac(iq,istate),1)*a3(iLJ3)/(r6_3_sc*r6_3_sc)
                                V_b3 = qbvdw(qiac(iq,istate),1)*b3(iLJ3)/(r6_3_sc)
                        end if
                        ! calculate  Vel, V_a, V_b and dv
                        Vel1 = crg1*qcrg(iq,istate)*r_1
                        Vel2 = crg2*qcrg(iq,istate)*r_2
                        Vel3 = crg3*qcrg(iq,istate)*r_3
                        dv1 = dv1 + r2_1*(-Vel1- ((12.*V_a1-6.*V_b1)*(r6_1_hc/r6_1_sc)) )*EQ(istate)%lambda
                        dv2 = dv2 + r2_2*(-Vel2- ((12.*V_a2-6.*V_b2)*(r6_2_hc/r6_2_sc)) )*EQ(istate)%lambda
                        dv3 = dv3 + r2_3*(-Vel3- ((12.*V_a3-6.*V_b3)*(r6_3_hc/r6_3_sc)) )*EQ(istate)%lambda
                        ! update q-water energies
                        EQ(istate)%qw%el  = EQ(istate)%qw%el + Vel1 + Vel2 + Vel3
                        EQ(istate)%qw%vdw = EQ(istate)%qw%vdw + V_a1 - V_b1 &
                                + V_a2 - V_b2 + V_a3 - V_b3
                end do !istate			
                ! update forces on Q-atom
                d(3*i-2) = d(3*i-2) - dv1*dx1 - dv2*dx2 - dv3*dx3
                d(3*i-1) = d(3*i-1) - dv1*dy1 - dv2*dy2 - dv3*dy3
                d(3*i  ) = d(3*i  ) - dv1*dz1 - dv2*dz2 - dv3*dz3
                ! update forces on water
                d(3*j-2) = d(3*j-2) + dv1*dx1
                d(3*j-1) = d(3*j-1) + dv1*dy1
                d(3*j  ) = d(3*j  ) + dv1*dz1
                d(3*j+1) = d(3*j+1) + dv2*dx2
                d(3*j+2) = d(3*j+2) + dv2*dy2
                d(3*j+3) = d(3*j+3) + dv2*dz2
                d(3*j+4) = d(3*j+4) + dv3*dx3
                d(3*j+5) = d(3*j+5) + dv3*dy3
                d(3*j+6) = d(3*j+6) + dv3*dz3

        end do !iq
end do !jw
end subroutine nonbond_qw_3atom

!-----------------------------------------------------------------------
!******PWadded 2001-10-23
subroutine nonbond_qw_3atom_box
	!calculate non-bonded interactions between Q-atoms and 3-atom solvent molecules
	!using geometric comb. rule
	! local variables
	integer						:: jw,iq,i,j,iLJ1, iLJ2, iLJ3, iaci
	integer						:: istate
	real(8)						::	dx1, dy1, dz1, dx2, dy2, dz2, dx3, dy3, dz3
	real(8)						::	r_1, r2_1, r6_1, r_2, r2_2, r6_2, r_3, r2_3, r6_3
	real(8)						::	Vel1, Vel2, Vel3, dv1, dv2, dv3
	real(8)						:: V_a1,V_b1, V_a2, V_b2, V_a3, V_b3,r6_1_sc,r6_2_sc,r6_3_sc,r6_1_hc,r6_2_hc,r6_3_hc
	real(8)						:: boxshiftx, boxshifty, boxshiftz, dx, dy, dz
	real(8), save				::	a1(2), b1(2), a2(2), b2(2), a3(2), b3(2)
	integer, save				::	iac1, iac2, iac3
	real, save					::	crg1, crg2, crg3
	! global variables used:
	!  iqseq, iac, crg, x, nstates, qvdw_flag, iaclib, qiac, qavdw, qbvdw, qcrg, el14_scale, EQ, d, nat_solute
  
	if(a1(1) == 0.) then !set first time
		iac1 = iac(nat_solute+1)
		iac2 = iac(nat_solute+2)
		iac3 = iac(nat_solute+3)
		crg1 = crg(nat_solute+1)
		crg2 = crg(nat_solute+2)
		crg3 = crg(nat_solute+3)
		a1(1:2) = iaclib(iac1)%avdw(1:2)
		b1(1:2) = iaclib(iac1)%bvdw(1:2)
		a2(1:2) = iaclib(iac2)%avdw(1:2)
		b2(1:2) = iaclib(iac2)%bvdw(1:2)
		a3(1:2) = iaclib(iac3)%avdw(1:2)
		b3(1:2) = iaclib(iac3)%bvdw(1:2)
	end if

	!loop over listed waters
	do jw = 1, nbqw_pair
		j = nbqw(jw) !top. # of O in water iw
		
		!compute the periodical shift
		dx = x(3*j-2) - x(3*qswitch-2)
		dy = x(3*j-1) - x(3*qswitch-1)
		dz = x(3*j  ) - x(3*qswitch  )
		boxshiftx = boxlength(1)*nint(dx*inv_boxl(1))
		boxshifty = boxlength(2)*nint(dy*inv_boxl(2))
		boxshiftz = boxlength(3)*nint(dz*inv_boxl(3))
				
		!loop over all q-atoms
		do iq = 1, nqat
			i = iqseq(iq)
			dx1 = x(3*j-2) - x(3*i-2)
			dy1 = x(3*j-1) - x(3*i-1)
			dz1 = x(3*j  ) - x(3*i  )
			dx2 = x(3*j+1) - x(3*i-2)
			dy2 = x(3*j+2) - x(3*i-1)
			dz2 = x(3*j+3) - x(3*i  )
			dx3 = x(3*j+4) - x(3*i-2)
			dy3 = x(3*j+5) - x(3*i-1)
			dz3 = x(3*j+6) - x(3*i  )
			dx1 = dx1 - boxshiftx  
			dy1 = dy1 - boxshifty  
			dz1 = dz1 - boxshiftz  
			dx2 = dx2 - boxshiftx 
			dy2 = dy2 - boxshifty 
			dz2 = dz2 - boxshiftz 		
			dx3 = dx3 - boxshiftx 
			dy3 = dy3 - boxshifty 
			dz3 = dz3 - boxshiftz  
			
            r2_1 = dx1*dx1 + dy1*dy1 + dz1*dz1
            r2_2 = dx2*dx2 + dy2*dy2 + dz2*dz2
            r2_3 = dx3*dx3 + dy3*dy3 + dz3*dz3
            r6_1_hc = r2_1*r2_1*r2_1   !will set r6 to 1/r6 later, needed for softcore
            r2_1 = 1._8/r2_1
			r_1  = sqrt(r2_1)
			r6_2_hc = r2_2*r2_2*r2_2   !will set r6 to 1/r6 later, needed for softcore
			r2_2 = 1._8/r2_2
			r_2  = sqrt(r2_2)
			r6_3_hc = r2_3*r2_3*r2_3   !will set r6 to 1/r6 later, needed for softcore
			r2_3 = 1._8/r2_3
			r_3  = sqrt(r2_3)

			!reset potential
			dv1 = 0
			dv2 = 0
			dv3 = 0
			iaci = iac(i)
			iLJ1 = LJcod(iac1, iaci)
			iLJ2 = LJcod(iac2, iaci)
			iLJ3 = LJcod(iac3, iaci)
			if(.not. qvdw_flag) then
				!use same LJ params for all states
				r6_1 = 1._8/r6_1_hc	!softcore hack
				r6_2 = 1._8/r6_2_hc	!softcore hack
				r6_3 = 1._8/r6_3_hc	!softcore hack

				V_a1  = iaclib(iaci)%avdw(iLJ1)*a1(iLJ1)*r6_1*r6_1 
				V_b1  = iaclib(iaci)%bvdw(iLJ1)*b1(iLJ1)*r6_1
				V_a2 = iaclib(iaci)%avdw(iLJ2)*a2(iLJ2)*r6_2*r6_2
				V_b2 = iaclib(iaci)%bvdw(iLJ2)*b2(iLJ2)*r6_2
				V_a3 = iaclib(iaci)%avdw(iLJ3)*a3(iLJ3)*r6_3*r6_3
				V_b3 = iaclib(iaci)%bvdw(iLJ3)*b3(iLJ3)*r6_3
			end if
			do istate = 1, nstates
				! for every state:

				! calculate V_a:s and V_b:s for each state
				if (qvdw_flag) then
					r6_1_sc = r6_1_hc + sc_lookup(iq,iac1,istate)
					r6_2_sc = r6_2_hc + sc_lookup(iq,iac2,istate)
					r6_3_sc = r6_3_hc + sc_lookup(iq,iac3,istate)
					V_a1 = qavdw(qiac(iq,istate),1)*a1(iLJ1)/(r6_1_sc*r6_1_sc)
					V_b1 = qbvdw(qiac(iq,istate),1)*b1(iLJ1)/(r6_1_sc)
					V_a2 = qavdw(qiac(iq,istate),1)*a2(iLJ2)/(r6_2_sc*r6_2_sc)
					V_b2 = qbvdw(qiac(iq,istate),1)*b2(iLJ2)/(r6_2_sc)
					V_a3 = qavdw(qiac(iq,istate),1)*a3(iLJ3)/(r6_3_sc*r6_3_sc)
					V_b3 = qbvdw(qiac(iq,istate),1)*b3(iLJ3)/(r6_3_sc)

				end if
				! calculate  Vel, V_a, V_b and dv
				Vel1 = crg1*qcrg(iq,istate)*r_1
				Vel2 = crg2*qcrg(iq,istate)*r_2
				Vel3 = crg3*qcrg(iq,istate)*r_3
				dv1 = dv1 + r2_1*(-Vel1- ((12.*V_a1-6.*V_b1)*(r6_1_hc/r6_1_sc)) )*EQ(istate)%lambda
				dv2 = dv2 + r2_2*(-Vel2- ((12.*V_a2-6.*V_b2)*(r6_2_hc/r6_2_sc)) )*EQ(istate)%lambda
				dv3 = dv3 + r2_3*(-Vel3- ((12.*V_a3-6.*V_b3)*(r6_3_hc/r6_3_sc)) )*EQ(istate)%lambda
				! update q-water energies
				EQ(istate)%qw%el  = EQ(istate)%qw%el + Vel1 + Vel2 + Vel3
				EQ(istate)%qw%vdw = EQ(istate)%qw%vdw + V_a1 - V_b1 &
					+ V_a2 - V_b2 + V_a3 - V_b3
			end do !istate			
			! update forces on Q-atom
			d(3*i-2) = d(3*i-2) - dv1*dx1 - dv2*dx2 - dv3*dx3
			d(3*i-1) = d(3*i-1) - dv1*dy1 - dv2*dy2 - dv3*dy3
			d(3*i  ) = d(3*i  ) - dv1*dz1 - dv2*dz2 - dv3*dz3
			! update forces on water
			d(3*j-2) = d(3*j-2) + dv1*dx1
			d(3*j-1) = d(3*j-1) + dv1*dy1
			d(3*j  ) = d(3*j  ) + dv1*dz1
			d(3*j+1) = d(3*j+1) + dv2*dx2
			d(3*j+2) = d(3*j+2) + dv2*dy2
			d(3*j+3) = d(3*j+3) + dv2*dz2
			d(3*j+4) = d(3*j+4) + dv3*dx3
			d(3*j+5) = d(3*j+5) + dv3*dy3
			d(3*j+6) = d(3*j+6) + dv3*dz3

		end do !iq
	end do !jw
end subroutine nonbond_qw_3atom_box

!-----------------------------------------------------------------------

subroutine nonbond_ww_spc
! local variables
integer						:: iw,ip,i,j,i3,j3, ia
integer						:: ipstart
real(8)						:: rOX, rH1X, rH2X, r2
real(8)						:: dxOX,  dyOX,  dzOX
real(8)						:: dxH1X, dyH1X, dzH1X
real(8)						:: dxH2X, dyH2X, dzH2X
real(8)						:: Vel,V_a,V_b,dv
real(8), save					:: A_OO, B_OO
integer						::	iac_ow, iac_hw

! global variables used:
!  iaclib, nat_solute, x, crg_ow, E, d, crg_hw

if(A_OO == 0) then !initialize 'em!
iac_ow = iac(nat_solute + 1)
iac_hw = iac(nat_solute + 2)
A_OO = iaclib(iac_ow)%avdw(ljcod(iac_ow, iac_ow)) &
        *iaclib(iac_ow)%avdw(ljcod(iac_ow, iac_ow))
B_OO = iaclib(iac_ow)%bvdw(ljcod(iac_ow, iac_ow)) &
        *iaclib(iac_ow)%bvdw(ljcod(iac_ow, iac_ow))
end if

ipstart = 1

do iw = calculation_assignment%ww%start, calculation_assignment%ww%end
ip = ipstart

do while (nbww(ip) .ne. 0)
  ! consider the pair (nat_solute+3*(iw-1)) - nbww(ip)

  ! --- O - (O,H1,H2) ---
  i3 = (nat_solute+3*(iw-1))*3+1 !point to x for O in mol. iw
  j3 = nbww(ip)*3-2 !point to O in interacting mol.
  ! O - O (X=O)
  dxOX = x((j3  ))-x((i3  ))
  dyOX = x((j3+1))-x((i3+1))
  dzOX = x((j3+2))-x((i3+2))
  rOX = dxOX*dxOX+dyOX*dyOX+dzOX*dzOX
  ! O-H1 (X=H1)
  j3 = j3 + 3

  dxH1X = x((j3  ))-x((i3  ))
  dyH1X = x((j3+1))-x((i3+1))
  dzH1X = x((j3+2))-x((i3+2))
  rH1X = dxH1X*dxH1X+dyH1X*dyH1X+dzH1X*dzH1X
  ! O-H2 (X=H2)
  j3 = j3 + 3
  dxH2X = x((j3  ))-x((i3  ))
  dyH2X = x((j3+1))-x((i3+1))
  dzH2X = x((j3+2))-x((i3+2))
  rH2X = dxH2X*dxH2X+dyH2X*dyH2X+dzH2X*dzH2X
  rOX=	sqrt(1._8/rOX  )
  rH1X=	sqrt(1._8/rH1X )
  rH2X=	sqrt(1._8/rH2X )
  ! O-O 
  ! LJ only for O-O
  r2 = rOX * rOX
  Vel = crg_ow*crg_ow*rOX
  V_a  = A_OO*(r2*r2*r2)*(r2*r2*r2)
  V_b  = B_OO*(r2*r2*r2)
  E%ww%vdw = E%ww%vdw + V_a-V_B
  E%ww%el = E%ww%el + Vel
  dv = r2*( -Vel -12.*V_a +6.*V_b )
  j3 = j3 - 6 !move pointer back to O in interacting mol.
  d((i3  )) = d((i3  )) -(dv*dxOX)
  d((j3  )) = d((j3  )) +(dv*dxOX)
  d((i3+1)) = d((i3+1)) -(dv*dyOX)
  d((j3+1)) = d((j3+1)) +(dv*dyOX)
  d((i3+2)) = d((i3+2)) -(dv*dzOX)
  d((j3+2)) = d((j3+2)) +(dv*dzOX)
  ! O-H1 
  r2 = rH1X * rH1X
  Vel = crg_ow*crg_hw*rH1X
  E%ww%el = E%ww%el + Vel
  dv   = r2*( -Vel)
  j3 = j3 +3 !point to H1 in j-molecule
  d((i3  )) = d((i3  )) -(dv*dxH1X)
  d((j3  )) = d((j3  )) +(dv*dxH1X)
  d((i3+1)) = d((i3+1)) -(dv*dyH1X)
  d((j3+1)) = d((j3+1)) +(dv*dyH1X)
  d((i3+2)) = d((i3+2)) -(dv*dzH1X)
  d((j3+2)) = d((j3+2)) +(dv*dzH1X)
  ! O-H2
  r2 = rH2X * rH2X
  Vel = crg_ow*crg_hw*rH2X
  E%ww%el = E%ww%el + Vel
  dv   = r2*( -Vel)
  j3 = j3 +3 !point to H2 in j-molecule
  d((i3  )) = d((i3  )) -(dv*dxH2X)
  d((j3  )) = d((j3  )) +(dv*dxH2X)
  d((i3+1)) = d((i3+1)) -(dv*dyH2X)
  d((j3+1)) = d((j3+1)) +(dv*dyH2X)
  d((i3+2)) = d((i3+2)) -(dv*dzH2X)
  d((j3+2)) = d((j3+2)) +(dv*dzH2X)

  ! --- H1 - (O,H1,H2) ---
  i3 = i3 + 3 !point to x for H1 in mol. iw
  j3 = nbww(ip)*3-2 !point to O in j-mol.
  ! H1 - O (X=O)
  dxOX = x((j3  ))-x((i3  ))
  dyOX = x((j3+1))-x((i3+1))
  dzOX = x((j3+2))-x((i3+2))
  rOX = dxOX*dxOX+dyOX*dyOX+dzOX*dzOX
  ! H1-H1 (X=H1)
  j3 = j3 + 3
  dxH1X = x((j3  ))-x((i3  ))
  dyH1X = x((j3+1))-x((i3+1))
  dzH1X = x((j3+2))-x((i3+2))
  rH1X = dxH1X*dxH1X+dyH1X*dyH1X+dzH1X*dzH1X
  ! H1-H2 (X=H2)
  j3 = j3 + 3
  dxH2X = x((j3  ))-x((i3  ))
  dyH2X = x((j3+1))-x((i3+1))
  dzH2X = x((j3+2))-x((i3+2))
  rH2X = dxH2X*dxH2X+dyH2X*dyH2X+dzH2X*dzH2X
  rOX=	sqrt(1._8/rOX  )
  rH1X=	sqrt(1._8/rH1X )
  rH2X=	sqrt(1._8/rH2X )
  ! H1-O 
  r2 = rOX * rOX
  Vel = crg_hw*crg_ow*rOX
  E%ww%el = E%ww%el + Vel
  dv   = r2*( -Vel)
  j3 = j3 - 6 !move pointer back to O in j-mol.
  d((i3  )) = d((i3  )) -(dv*dxOX)
  d((j3  )) = d((j3  )) +(dv*dxOX)
  d((i3+1)) = d((i3+1)) -(dv*dyOX)
  d((j3+1)) = d((j3+1)) +(dv*dyOX)
  d((i3+2)) = d((i3+2)) -(dv*dzOX)
  d((j3+2)) = d((j3+2)) +(dv*dzOX)
  ! H1-H1 
  r2 = rH1X * rH1X
  Vel = crg_hw*crg_hw*rH1X
  E%ww%el = E%ww%el + Vel
  dv   = r2*( -Vel)
  j3 = j3 +3 !point to H1 in j-molecule
  d((i3  )) = d((i3  )) -(dv*dxH1X)
  d((j3  )) = d((j3  )) +(dv*dxH1X)
  d((i3+1)) = d((i3+1)) -(dv*dyH1X)
  d((j3+1)) = d((j3+1)) +(dv*dyH1X)
  d((i3+2)) = d((i3+2)) -(dv*dzH1X)
  d((j3+2)) = d((j3+2)) +(dv*dzH1X)
  ! H1-H2
  r2 = rH2X * rH2X
  Vel = crg_hw*crg_hw*rH2X
  E%ww%el = E%ww%el + Vel
  dv   = r2*( -Vel)
  j3 = j3 +3 !point to H2 in j-molecule
  d((i3  )) = d((i3  )) -(dv*dxH2X)
  d((j3  )) = d((j3  )) +(dv*dxH2X)
  d((i3+1)) = d((i3+1)) -(dv*dyH2X)
  d((j3+1)) = d((j3+1)) +(dv*dyH2X)
  d((i3+2)) = d((i3+2)) -(dv*dzH2X)
  d((j3+2)) = d((j3+2)) +(dv*dzH2X)

  ! --- H2 - (O,H1,H2) ---
  i3 = i3 + 3 !point to x for H2 in mol. iw
  j3 = nbww(ip)*3-2 !point to O in j-mol.
  ! H2 - O (X=O)
  dxOX = x((j3  ))-x((i3  ))
  dyOX = x((j3+1))-x((i3+1))
  dzOX = x((j3+2))-x((i3+2))
  rOX = dxOX*dxOX+dyOX*dyOX+dzOX*dzOX
  ! H2-H1 (X=H1)
  j3 = j3 + 3
  dxH1X = x((j3  ))-x((i3  ))
  dyH1X = x((j3+1))-x((i3+1))
  dzH1X = x((j3+2))-x((i3+2))
  rH1X = dxH1X*dxH1X+dyH1X*dyH1X+dzH1X*dzH1X
  ! H2-H2 (X=H2)
  j3 = j3 + 3
  dxH2X = x((j3  ))-x((i3  ))
  dyH2X = x((j3+1))-x((i3+1))
  dzH2X = x((j3+2))-x((i3+2))
  rH2X = dxH2X*dxH2X+dyH2X*dyH2X+dzH2X*dzH2X
  rOX=	sqrt(1._8/rOX  )
  rH1X=	sqrt(1._8/rH1X )
  rH2X=	sqrt(1._8/rH2X )
  ! H2-O 
  r2 = rOX * rOX
  Vel = crg_hw*crg_ow*rOX
  E%ww%el = E%ww%el + Vel
  dv   = r2*( -Vel)
  j3 = j3 - 6 !move pointer back to O in j-mol.
  d((i3  )) = d((i3  )) -(dv*dxOX)
  d((j3  )) = d((j3  )) +(dv*dxOX)
  d((i3+1)) = d((i3+1)) -(dv*dyOX)
  d((j3+1)) = d((j3+1)) +(dv*dyOX)
  d((i3+2)) = d((i3+2)) -(dv*dzOX)
  d((j3+2)) = d((j3+2)) +(dv*dzOX)
  ! H2-H1 
  r2 = rH1X * rH1X
  Vel = crg_hw*crg_hw*rH1X
  E%ww%el = E%ww%el + Vel
  dv   = r2*( -Vel)
  j3 = j3 +3 !point to H1 in j-molecule
  d((i3  )) = d((i3  )) -(dv*dxH1X)
  d((j3  )) = d((j3  )) +(dv*dxH1X)
  d((i3+1)) = d((i3+1)) -(dv*dyH1X)
  d((j3+1)) = d((j3+1)) +(dv*dyH1X)
  d((i3+2)) = d((i3+2)) -(dv*dzH1X)
  d((j3+2)) = d((j3+2)) +(dv*dzH1X)
  ! H2-H2
  r2 = rH2X * rH2X
  Vel = crg_hw*crg_hw*rH2X
  E%ww%el = E%ww%el + Vel
  dv   = r2*( -Vel)
  j3 = j3 +3 !point to H2 in j-molecule
  d((i3  )) = d((i3  )) -(dv*dxH2X)
  d((j3  )) = d((j3  )) +(dv*dxH2X)
  d((i3+1)) = d((i3+1)) -(dv*dyH2X)
  d((j3+1)) = d((j3+1)) +(dv*dyH2X)
  d((i3+2)) = d((i3+2)) -(dv*dzH2X)
  d((j3+2)) = d((j3+2)) +(dv*dzH2X)

  ip = ip + 1
end do ! while ip

! skip the gap
ipstart = ip + 1			
end do ! iw

end subroutine nonbond_ww_spc

!-----------------------------------------------------------------------
!******PWadded 2001-10-23
subroutine nonbond_ww_spc_box
  ! local variables
  integer						:: iw,ip,i,j,i3,j3, ia
  integer						:: ipstart
  real(8)						:: rOX, rH1X, rH2X, r2
  real(8)						:: dxOX,  dyOX,  dzOX
  real(8)						:: dxH1X, dyH1X, dzH1X
  real(8)						:: dxH2X, dyH2X, dzH2X
  real(8)						:: Vel,V_a,V_b,dv
  real(8), save					:: A_OO, B_OO
  integer						::	iac_ow, iac_hw
  real(8)						:: boxshiftx, boxshifty, boxshiftz

  ! global variables used:
  !  iaclib, nat_solute, x, crg_ow, E, d, crg_hw

  if(A_OO == 0) then !initialize 'em!
	iac_ow = iac(nat_solute + 1)
	iac_hw = iac(nat_solute + 2)
	A_OO = iaclib(iac_ow)%avdw(ljcod(iac_ow, iac_ow)) &
		*iaclib(iac_ow)%avdw(ljcod(iac_ow, iac_ow))
	B_OO = iaclib(iac_ow)%bvdw(ljcod(iac_ow, iac_ow)) &
		*iaclib(iac_ow)%bvdw(ljcod(iac_ow, iac_ow))
  end if

  ipstart = 1

  do iw = calculation_assignment%ww%start, calculation_assignment%ww%end
	ip = ipstart

	do while (nbww(ip) .ne. 0)
	  ! consider the pair (nat_solute+3*(iw-1)) - nbww(ip)

	  ! --- O - (O,H1,H2) ---
	  i3 = (nat_solute+3*(iw-1))*3+1 !point to x for O in mol. iw
	  j3 = nbww(ip)*3-2 !point to O in interacting mol.
	  ! O - O (X=O)
	  dxOX = x((j3  ))-x((i3  ))
	  dyOX = x((j3+1))-x((i3+1))
	  dzOX = x((j3+2))-x((i3+2))
	  !the periodical shift
	  boxshiftx = boxlength(1)*nint( dxOX*inv_boxl(1) )
	  boxshifty = boxlength(2)*nint( dyOX*inv_boxl(2) )
	  boxshiftz = boxlength(3)*nint( dzOX*inv_boxl(3) )
	  ! O-H1 (X=H1)
	  j3 = j3 + 3
	  dxH1X = x((j3  ))-x((i3  ))
	  dyH1X = x((j3+1))-x((i3+1))
	  dzH1X = x((j3+2))-x((i3+2))
	  ! O-H2 (X=H2)
	  j3 = j3 + 3
	  dxH2X = x((j3  ))-x((i3  ))
	  dyH2X = x((j3+1))-x((i3+1))
	  dzH2X = x((j3+2))-x((i3+2))
	  dxOX = dxOX - boxshiftx
	  dyOX = dyOX - boxshifty
	  dzOX = dzOX - boxshiftz
	  dxH1X = dxH1X - boxshiftx
	  dyH1X = dyH1X - boxshifty
	  dzH1X = dzH1X - boxshiftz
	  dxH2X = dxH2X - boxshiftx
	  dyH2X = dyH2X - boxshifty
	  dzH2X = dzH2X - boxshiftz

	  rOX = dxOX*dxOX + dyOX*dyOX + dzOX*dzOX
	  rH1X = dxH1X*dxH1X + dyH1X*dyH1X + dzH1X*dzH1X
	  rH2X = dxH2X*dxH2X + dyH2X*dyH2X + dzH2X*dzH2X

	  rOX=	sqrt(1._8/rOX  )
	  rH1X=	sqrt(1._8/rH1X )
	  rH2X=	sqrt(1._8/rH2X )
	  ! O-O 
	  ! LJ only for O-O
	  r2 = rOX * rOX
	  Vel = crg_ow*crg_ow*rOX
	  V_a  = A_OO*(r2*r2*r2)*(r2*r2*r2)
	  V_b  = B_OO*(r2*r2*r2)
	  E%ww%vdw = E%ww%vdw + V_a-V_B
	  E%ww%el = E%ww%el + Vel
	  dv = r2*( -Vel -12.*V_a +6.*V_b )
	  j3 = j3 - 6 !move pointer back to O in interacting mol.
	  d((i3  )) = d((i3  )) -(dv*dxOX)
	  d((j3  )) = d((j3  )) +(dv*dxOX)
	  d((i3+1)) = d((i3+1)) -(dv*dyOX)
	  d((j3+1)) = d((j3+1)) +(dv*dyOX)
	  d((i3+2)) = d((i3+2)) -(dv*dzOX)
	  d((j3+2)) = d((j3+2)) +(dv*dzOX)
	  ! O-H1 
	  r2 = rH1X * rH1X
	  Vel = crg_ow*crg_hw*rH1X
	  E%ww%el = E%ww%el + Vel
	  dv   = r2*( -Vel)
	  j3 = j3 +3 !point to H1 in j-molecule
	  d((i3  )) = d((i3  )) -(dv*dxH1X)
	  d((j3  )) = d((j3  )) +(dv*dxH1X)
	  d((i3+1)) = d((i3+1)) -(dv*dyH1X)
	  d((j3+1)) = d((j3+1)) +(dv*dyH1X)
	  d((i3+2)) = d((i3+2)) -(dv*dzH1X)
	  d((j3+2)) = d((j3+2)) +(dv*dzH1X)
	  ! O-H2
	  r2 = rH2X * rH2X
	  Vel = crg_ow*crg_hw*rH2X
	  E%ww%el = E%ww%el + Vel
	  dv   = r2*( -Vel)
	  j3 = j3 +3 !point to H2 in j-molecule
	  d((i3  )) = d((i3  )) -(dv*dxH2X)
	  d((j3  )) = d((j3  )) +(dv*dxH2X)
	  d((i3+1)) = d((i3+1)) -(dv*dyH2X)
	  d((j3+1)) = d((j3+1)) +(dv*dyH2X)
	  d((i3+2)) = d((i3+2)) -(dv*dzH2X)
	  d((j3+2)) = d((j3+2)) +(dv*dzH2X)
			
	  ! --- H1 - (O,H1,H2) ---
	  i3 = i3 + 3 !point to x for H1 in mol. iw
	  j3 = nbww(ip)*3-2 !point to O in j-mol.
	  ! H1 - O (X=O)
	  dxOX = x((j3  ))-x((i3  ))
	  dyOX = x((j3+1))-x((i3+1))
	  dzOX = x((j3+2))-x((i3+2))
	  ! H1-H1 (X=H1)
	  j3 = j3 + 3
	  dxH1X = x((j3  ))-x((i3  ))
	  dyH1X = x((j3+1))-x((i3+1))
	  dzH1X = x((j3+2))-x((i3+2))
	  ! H1-H2 (X=H2)
	  j3 = j3 + 3
	  dxH2X = x((j3  ))-x((i3  ))
	  dyH2X = x((j3+1))-x((i3+1))
	  dzH2X = x((j3+2))-x((i3+2))
	  dxOX = dxOX - boxshiftx
	  dyOX = dyOX - boxshifty
	  dzOX = dzOX - boxshiftz
	  dxH1X = dxH1X - boxshiftx
	  dyH1X = dyH1X - boxshifty
	  dzH1X = dzH1X - boxshiftz
	  dxH2X = dxH2X - boxshiftx
	  dyH2X = dyH2X - boxshifty
	  dzH2X = dzH2X - boxshiftz

	  rOX = dxOX*dxOX + dyOX*dyOX + dzOX*dzOX
	  rH1X = dxH1X*dxH1X + dyH1X*dyH1X + dzH1X*dzH1X
	  rH2X = dxH2X*dxH2X + dyH2X*dyH2X + dzH2X*dzH2X
	  rOX=	sqrt(1._8/rOX  )
	  rH1X=	sqrt(1._8/rH1X )
	  rH2X=	sqrt(1._8/rH2X )
	  ! H1-O 
	  r2 = rOX * rOX
	  Vel = crg_hw*crg_ow*rOX
	  E%ww%el = E%ww%el + Vel
	  dv   = r2*( -Vel)
	  j3 = j3 - 6 !move pointer back to O in j-mol.
	  d((i3  )) = d((i3  )) -(dv*dxOX)
	  d((j3  )) = d((j3  )) +(dv*dxOX)
	  d((i3+1)) = d((i3+1)) -(dv*dyOX)
	  d((j3+1)) = d((j3+1)) +(dv*dyOX)
	  d((i3+2)) = d((i3+2)) -(dv*dzOX)
	  d((j3+2)) = d((j3+2)) +(dv*dzOX)
  ! H1-H1 
	  r2 = rH1X * rH1X
	  Vel = crg_hw*crg_hw*rH1X
	  E%ww%el = E%ww%el + Vel
	  dv   = r2*( -Vel)
	  j3 = j3 +3 !point to H1 in j-molecule
	  d((i3  )) = d((i3  )) -(dv*dxH1X)
	  d((j3  )) = d((j3  )) +(dv*dxH1X)
	  d((i3+1)) = d((i3+1)) -(dv*dyH1X)
	  d((j3+1)) = d((j3+1)) +(dv*dyH1X)
	  d((i3+2)) = d((i3+2)) -(dv*dzH1X)
	  d((j3+2)) = d((j3+2)) +(dv*dzH1X)
	  ! H1-H2
	  r2 = rH2X * rH2X
	  Vel = crg_hw*crg_hw*rH2X
	  E%ww%el = E%ww%el + Vel
	  dv   = r2*( -Vel)
	  j3 = j3 +3 !point to H2 in j-molecule
	  d((i3  )) = d((i3  )) -(dv*dxH2X)
	  d((j3  )) = d((j3  )) +(dv*dxH2X)
	  d((i3+1)) = d((i3+1)) -(dv*dyH2X)
	  d((j3+1)) = d((j3+1)) +(dv*dyH2X)
	  d((i3+2)) = d((i3+2)) -(dv*dzH2X)
	  d((j3+2)) = d((j3+2)) +(dv*dzH2X)
			
	  ! --- H2 - (O,H1,H2) ---
	  i3 = i3 + 3 !point to x for H2 in mol. iw
	  j3 = nbww(ip)*3-2 !point to O in j-mol.
	  ! H2 - O (X=O)
	  dxOX = x((j3  ))-x((i3  ))
	  dyOX = x((j3+1))-x((i3+1))
	  dzOX = x((j3+2))-x((i3+2))
	  ! H2-H1 (X=H1)
	  j3 = j3 + 3
	  dxH1X = x((j3  ))-x((i3  ))
	  dyH1X = x((j3+1))-x((i3+1))
	  dzH1X = x((j3+2))-x((i3+2))
	  ! H2-H2 (X=H2)
	  j3 = j3 + 3
	  dxH2X = x((j3  ))-x((i3  ))
	  dyH2X = x((j3+1))-x((i3+1))
	  dzH2X = x((j3+2))-x((i3+2))
	  dxOX = dxOX - boxshiftx
	  dyOX = dyOX - boxshifty
	  dzOX = dzOX - boxshiftz
	  dxH1X = dxH1X - boxshiftx
	  dyH1X = dyH1X - boxshifty
	  dzH1X = dzH1X - boxshiftz
	  dxH2X = dxH2X - boxshiftx
	  dyH2X = dyH2X - boxshifty
	  dzH2X = dzH2X - boxshiftz

	  rOX = dxOX*dxOX + dyOX*dyOX + dzOX*dzOX
	  rH1X = dxH1X*dxH1X + dyH1X*dyH1X + dzH1X*dzH1X
	  rH2X = dxH2X*dxH2X + dyH2X*dyH2X + dzH2X*dzH2X
	  rOX=	sqrt(1._8/rOX  )
	  rH1X=	sqrt(1._8/rH1X )
	  rH2X=	sqrt(1._8/rH2X )
	  ! H2-O 
	  r2 = rOX * rOX
	  Vel = crg_hw*crg_ow*rOX
	  E%ww%el = E%ww%el + Vel
	  dv   = r2*( -Vel)
	  j3 = j3 - 6 !move pointer back to O in j-mol.
	  d((i3  )) = d((i3  )) -(dv*dxOX)
	  d((j3  )) = d((j3  )) +(dv*dxOX)
	  d((i3+1)) = d((i3+1)) -(dv*dyOX)
	  d((j3+1)) = d((j3+1)) +(dv*dyOX)
	  d((i3+2)) = d((i3+2)) -(dv*dzOX)
	  d((j3+2)) = d((j3+2)) +(dv*dzOX)
	  ! H2-H1 
	  r2 = rH1X * rH1X
	  Vel = crg_hw*crg_hw*rH1X
	  E%ww%el = E%ww%el + Vel
	  dv   = r2*( -Vel)
	  j3 = j3 +3 !point to H1 in j-molecule
	  d((i3  )) = d((i3  )) -(dv*dxH1X)
	  d((j3  )) = d((j3  )) +(dv*dxH1X)
	  d((i3+1)) = d((i3+1)) -(dv*dyH1X)
	  d((j3+1)) = d((j3+1)) +(dv*dyH1X)
	  d((i3+2)) = d((i3+2)) -(dv*dzH1X)
	  d((j3+2)) = d((j3+2)) +(dv*dzH1X)
	  ! H2-H2
	  r2 = rH2X * rH2X
	  Vel = crg_hw*crg_hw*rH2X
	  E%ww%el = E%ww%el + Vel
	  dv   = r2*( -Vel)
	  j3 = j3 +3 !point to H2 in j-molecule
	  d((i3  )) = d((i3  )) -(dv*dxH2X)
	  d((j3  )) = d((j3  )) +(dv*dxH2X)
	  d((i3+1)) = d((i3+1)) -(dv*dyH2X)
	  d((j3+1)) = d((j3+1)) +(dv*dyH2X)
	  d((i3+2)) = d((i3+2)) -(dv*dzH2X)
	  d((j3+2)) = d((j3+2)) +(dv*dzH2X)
			
	  ip = ip + 1
	end do ! while ip

	! skip the gap
	ipstart = ip + 1			
  end do ! iw

end subroutine nonbond_ww_spc_box

!-----------------------------------------------------------------------

subroutine nonbond_3atomsolvent
! local variables
integer						:: iw,ip,i,j,i3,j3, ia
integer						:: ipstart
real(8)						::	r1X, r2X, r3X, r2
real(8)						::	dx1X,  dy1X,  dz1X
real(8)						::	dx2X, dy2X, dz2X
real(8)						::	dx3X, dy3X, dz3X
real(8)						::	Vel,V_a,V_b,dv
real(8), save					::	A_11, B_11, A_12, B_12, A_13, B_13
real(8), save					::	A_22, B_22, A_23, B_23, A_33, B_33
real, save					::	crg1, crg2, crg3
integer						::	iac1, iac2, iac3

! global variables used:
!  iaclib, nat_solute, x, E, d

if(A_11 == 0.) then !initialize static variables
iac1 = iac(nat_solute+1)
iac2 = iac(nat_solute+2)
iac3 = iac(nat_solute+3)
crg1 = crg(nat_solute+1)
crg2 = crg(nat_solute+2)
crg3 = crg(nat_solute+3)

A_11 = iaclib(iac1)%avdw(ljcod(iac1, iac1)) &
        *iaclib(iac1)%avdw(ljcod(iac1, iac1))
B_11 = iaclib(iac1)%bvdw(ljcod(iac1, iac1)) &
        *iaclib(iac1)%bvdw(ljcod(iac1, iac1))
A_12 = iaclib(iac1)%avdw(ljcod(iac1, iac2)) &
        *iaclib(iac2)%avdw(ljcod(iac1, iac2))
B_12 = iaclib(iac1)%bvdw(ljcod(iac1, iac2)) &
        *iaclib(iac2)%bvdw(ljcod(iac1, iac2))
A_13 = iaclib(iac1)%avdw(ljcod(iac1, iac3)) &
        *iaclib(iac3)%avdw(ljcod(iac1, iac3))
B_13 = iaclib(iac1)%bvdw(ljcod(iac1, iac3)) &
        *iaclib(iac3)%bvdw(ljcod(iac1, iac3))
A_22 = iaclib(iac2)%avdw(ljcod(iac2, iac2)) &
        *iaclib(iac2)%avdw(ljcod(iac2, iac2))
B_22 = iaclib(iac2)%bvdw(ljcod(iac2, iac2)) &
        *iaclib(iac2)%bvdw(ljcod(iac2, iac2))
A_23 = iaclib(iac2)%avdw(ljcod(iac2, iac3)) &
        *iaclib(iac3)%avdw(ljcod(iac2, iac3))
B_23 = iaclib(iac2)%bvdw(ljcod(iac2, iac3)) &
        *iaclib(iac3)%bvdw(ljcod(iac2, iac3))
A_33 = iaclib(iac3)%avdw(ljcod(iac3, iac3)) &
        *iaclib(iac3)%avdw(ljcod(iac3, iac3))
B_33 = iaclib(iac3)%bvdw(ljcod(iac3, iac3)) &
        *iaclib(iac3)%bvdw(ljcod(iac3, iac3))
end if

ipstart = 1

do iw = calculation_assignment%ww%start, calculation_assignment%ww%end
ip = ipstart

do while (nbww(ip) /= 0)
  ! consider the pair (nat_solute+3*(iw-1)) - nbww(ip)

  ! --- 1 - (1,2,3) ---
  i3 = (nat_solute+3*(iw-1))*3+1 !point to x for atom 1 in mol. iw
  j3 = nbww(ip)*3-2 !point to 1 in interacting mol.
  ! 1 - 1
  dx1X = x((j3  ))-x((i3  ))
  dy1X = x((j3+1))-x((i3+1))
  dz1X = x((j3+2))-x((i3+2))
  r1X = dx1X*dx1X+dy1X*dy1X+dz1X*dz1X
  ! 1-2
  j3 = j3 + 3
  dx2X = x((j3  ))-x((i3  ))
  dy2X = x((j3+1))-x((i3+1))
  dz2X = x((j3+2))-x((i3+2))
  r2X = dx2X*dx2X+dy2X*dy2X+dz2X*dz2X
  ! 1-3

  j3 = j3 + 3
  dx3X = x((j3  ))-x((i3  ))
  dy3X = x((j3+1))-x((i3+1))
  dz3X = x((j3+2))-x((i3+2))
  r3X = dx3X*dx3X+dy3X*dy3X+dz3X*dz3X
  r1X=	sqrt(1._8/r1X  )
  r2X=	sqrt(1._8/r2X )
  r3X=	sqrt(1._8/r3X )
  ! 1-1 
  r2 = r1X * r1X
  Vel = crg1*crg1*r1X
  V_a  = A_11*(r2*r2*r2)*(r2*r2*r2)
  V_b  = B_11*(r2*r2*r2)
  E%ww%vdw = E%ww%vdw + V_a-V_B
  E%ww%el = E%ww%el + Vel
  dv = r2*( -Vel -12.*V_a +6.*V_b )
  j3 = j3 - 6 !move pointer back to 1 in interacting mol.
  d((i3  )) = d((i3  )) -(dv*dx1X)
  d((j3  )) = d((j3  )) +(dv*dx1X)
  d((i3+1)) = d((i3+1)) -(dv*dy1X)
  d((j3+1)) = d((j3+1)) +(dv*dy1X)
  d((i3+2)) = d((i3+2)) -(dv*dz1X)
  d((j3+2)) = d((j3+2)) +(dv*dz1X)
  ! 1-2 
  r2 = r2X * r2X
  Vel = crg1*crg2*r2X
  V_a  = A_12*(r2*r2*r2)*(r2*r2*r2)
  V_b  = B_12*(r2*r2*r2)
  E%ww%vdw = E%ww%vdw + V_a-V_B
  E%ww%el = E%ww%el + Vel
  dv = r2*( -Vel -12.*V_a +6.*V_b )
  j3 = j3 +3 !point to 2 in j-molecule
  d((i3  )) = d((i3  )) -(dv*dx2X)
  d((j3  )) = d((j3  )) +(dv*dx2X)
  d((i3+1)) = d((i3+1)) -(dv*dy2X)
  d((j3+1)) = d((j3+1)) +(dv*dy2X)
  d((i3+2)) = d((i3+2)) -(dv*dz2X)
  d((j3+2)) = d((j3+2)) +(dv*dz2X)
  ! 1-3
  r2 = r3X * r3X
  Vel = crg1*crg3*r3X
  V_a  = A_13*(r2*r2*r2)*(r2*r2*r2)
  V_b  = B_13*(r2*r2*r2)
  E%ww%vdw = E%ww%vdw + V_a-V_B
  E%ww%el = E%ww%el + Vel
  dv = r2*( -Vel -12.*V_a +6.*V_b )
  j3 = j3 +3 !point to 3 in j-molecule
  d((i3  )) = d((i3  )) -(dv*dx3X)
  d((j3  )) = d((j3  )) +(dv*dx3X)
  d((i3+1)) = d((i3+1)) -(dv*dy3X)
  d((j3+1)) = d((j3+1)) +(dv*dy3X)
  d((i3+2)) = d((i3+2)) -(dv*dz3X)
  d((j3+2)) = d((j3+2)) +(dv*dz3X)

  ! --- 2 - (1,2,3) ---
  i3 = i3 + 3 !point to x for 2 in mol. iw
  j3 = nbww(ip)*3-2 !point to o in j-mol.
  ! 2 - 1 (X=1)
  dx1X = x((j3  ))-x((i3  ))
  dy1X = x((j3+1))-x((i3+1))
  dz1X = x((j3+2))-x((i3+2))
  r1X = dx1X*dx1X+dy1X*dy1X+dz1X*dz1X
  ! 2-2 (X=2)
  j3 = j3 + 3
  dx2X = x((j3  ))-x((i3  ))
  dy2X = x((j3+1))-x((i3+1))
  dz2X = x((j3+2))-x((i3+2))
  r2X = dx2X*dx2X+dy2X*dy2X+dz2X*dz2X
  ! 2-3 (X=3)
  j3 = j3 + 3
  dx3X = x((j3  ))-x((i3  ))
  dy3X = x((j3+1))-x((i3+1))
  dz3X = x((j3+2))-x((i3+2))
  r3X = dx3X*dx3X+dy3X*dy3X+dz3X*dz3X
  r1X=	sqrt(1._8/r1X  )
  r2X=	sqrt(1._8/r2X )
  r3X=	sqrt(1._8/r3X )
  ! 2-1 
  r2 = r1X * r1X
  Vel = crg2*crg1*r1X
  V_a  = A_12*(r2*r2*r2)*(r2*r2*r2)
  V_b  = B_12*(r2*r2*r2)
  E%ww%vdw = E%ww%vdw + V_a-V_B
  E%ww%el = E%ww%el + Vel
  dv = r2*( -Vel -12.*V_a +6.*V_b )
  j3 = j3 - 6 !move pointer back to o in j-mol.
  d((i3  )) = d((i3  )) -(dv*dx1X)
  d((j3  )) = d((j3  )) +(dv*dx1X)
  d((i3+1)) = d((i3+1)) -(dv*dy1X)
  d((j3+1)) = d((j3+1)) +(dv*dy1X)
  d((i3+2)) = d((i3+2)) -(dv*dz1X)
  d((j3+2)) = d((j3+2)) +(dv*dz1X)
  ! 2-2 
  r2 = r2X * r2X
  Vel = crg2*crg2*r2X
  V_a  = A_22*(r2*r2*r2)*(r2*r2*r2)
  V_b  = B_22*(r2*r2*r2)
  E%ww%vdw = E%ww%vdw + V_a-V_B
  E%ww%el = E%ww%el + Vel
  dv = r2*( -Vel -12.*V_a +6.*V_b )
  j3 = j3 +3 !point to 2 in j-molecule
  d((i3  )) = d((i3  )) -(dv*dx2X)
  d((j3  )) = d((j3  )) +(dv*dx2X)
  d((i3+1)) = d((i3+1)) -(dv*dy2X)
  d((j3+1)) = d((j3+1)) +(dv*dy2X)
  d((i3+2)) = d((i3+2)) -(dv*dz2X)
  d((j3+2)) = d((j3+2)) +(dv*dz2X)
  ! 2-3
  r2 = r3X * r3X
  Vel = crg2*crg3*r3X
  V_a  = A_23*(r2*r2*r2)*(r2*r2*r2)
  V_b  = B_23*(r2*r2*r2)
  E%ww%vdw = E%ww%vdw + V_a-V_B
  E%ww%el = E%ww%el + Vel
  dv = r2*( -Vel -12.*V_a +6.*V_b )
  j3 = j3 +3 !point to 3 in j-molecule
  d((i3  )) = d((i3  )) -(dv*dx3X)
  d((j3  )) = d((j3  )) +(dv*dx3X)
  d((i3+1)) = d((i3+1)) -(dv*dy3X)
  d((j3+1)) = d((j3+1)) +(dv*dy3X)
  d((i3+2)) = d((i3+2)) -(dv*dz3X)
  d((j3+2)) = d((j3+2)) +(dv*dz3X)

  ! --- 3 - (1,2,3) ---
  i3 = i3 + 3 !point to x for 3 in mol. iw
  j3 = nbww(ip)*3-2 !point to o in j-mol.
  ! 3 - 1 (X=1)
  dx1X = x((j3  ))-x((i3  ))
  dy1X = x((j3+1))-x((i3+1))
  dz1X = x((j3+2))-x((i3+2))
  r1X = dx1X*dx1X+dy1X*dy1X+dz1X*dz1X
  ! 3-2 (X=2)
  j3 = j3 + 3
  dx2X = x((j3  ))-x((i3  ))
  dy2X = x((j3+1))-x((i3+1))
  dz2X = x((j3+2))-x((i3+2))
  r2X = dx2X*dx2X+dy2X*dy2X+dz2X*dz2X
  ! 3-3 (X=3)
  j3 = j3 + 3
  dx3X = x((j3  ))-x((i3  ))
  dy3X = x((j3+1))-x((i3+1))
  dz3X = x((j3+2))-x((i3+2))
  r3X = dx3X*dx3X+dy3X*dy3X+dz3X*dz3X
  r1X=	sqrt(1._8/r1X  )
  r2X=	sqrt(1._8/r2X )
  r3X=	sqrt(1._8/r3X )
  ! 3-1 
  r2 = r1X * r1X
  Vel = crg3*crg1*r1X
  V_a  = A_13*(r2*r2*r2)*(r2*r2*r2)
  V_b  = B_13*(r2*r2*r2)
  E%ww%vdw = E%ww%vdw + V_a-V_B
  E%ww%el = E%ww%el + Vel
  dv = r2*( -Vel -12.*V_a +6.*V_b )
  j3 = j3 - 6 !move pointer back to o in j-mol.
  d((i3  )) = d((i3  )) -(dv*dx1X)
  d((j3  )) = d((j3  )) +(dv*dx1X)
  d((i3+1)) = d((i3+1)) -(dv*dy1X)
  d((j3+1)) = d((j3+1)) +(dv*dy1X)
  d((i3+2)) = d((i3+2)) -(dv*dz1X)
  d((j3+2)) = d((j3+2)) +(dv*dz1X)
  ! 3-2 
  r2 = r2X * r2X
  Vel = crg3*crg2*r2X
  V_a  = A_23*(r2*r2*r2)*(r2*r2*r2)
  V_b  = B_23*(r2*r2*r2)
  E%ww%vdw = E%ww%vdw + V_a-V_B
  E%ww%el = E%ww%el + Vel
  dv = r2*( -Vel -12.*V_a +6.*V_b )
  j3 = j3 +3 !point to 2 in j-molecule
  d((i3  )) = d((i3  )) -(dv*dx2X)
  d((j3  )) = d((j3  )) +(dv*dx2X)
  d((i3+1)) = d((i3+1)) -(dv*dy2X)
  d((j3+1)) = d((j3+1)) +(dv*dy2X)
  d((i3+2)) = d((i3+2)) -(dv*dz2X)
  d((j3+2)) = d((j3+2)) +(dv*dz2X)
  ! 3-3
  r2 = r3X * r3X
  Vel = crg3*crg3*r3X
  V_a  = A_33*(r2*r2*r2)*(r2*r2*r2)
  V_b  = B_33*(r2*r2*r2)
  E%ww%vdw = E%ww%vdw + V_a-V_B
  E%ww%el = E%ww%el + Vel
  dv = r2*( -Vel -12.*V_a +6.*V_b )
  j3 = j3 +3 !point to 3 in j-molecule
  d((i3  )) = d((i3  )) -(dv*dx3X)
  d((j3  )) = d((j3  )) +(dv*dx3X)
  d((i3+1)) = d((i3+1)) -(dv*dy3X)
  d((j3+1)) = d((j3+1)) +(dv*dy3X)
  d((i3+2)) = d((i3+2)) -(dv*dz3X)
  d((j3+2)) = d((j3+2)) +(dv*dz3X)

  ip = ip + 1
end do ! while ip

! skip the gap
ipstart = ip + 1			
end do ! iw

end subroutine nonbond_3atomsolvent

!-----------------------------------------------------------------------
!******PWadded 2001-10-23
subroutine nonbond_3atomsolvent_box
  ! local variables
  integer						:: iw,ip,i,j,i3,j3, ia
  integer						:: ipstart
  real(8)						::	r1X, r2X, r3X, r2
  real(8)						::	dx1X,  dy1X,  dz1X
  real(8)						::	dx2X, dy2X, dz2X
  real(8)						::	dx3X, dy3X, dz3X
  real(8)						::	Vel,V_a,V_b,dv
  real(8), save					::	A_11, B_11, A_12, B_12, A_13, B_13
  real(8), save					::	A_22, B_22, A_23, B_23, A_33, B_33
  real, save					::	crg1, crg2, crg3
  integer						::	iac1, iac2, iac3
  real(8)						::  boxshiftx, boxshifty, boxshiftz

  ! global variables used:
  !  iaclib, nat_solute, x, E, d

  if(A_11 == 0.) then !initialize static variables
	iac1 = iac(nat_solute+1)
	iac2 = iac(nat_solute+2)
	iac3 = iac(nat_solute+3)
	crg1 = crg(nat_solute+1)
	crg2 = crg(nat_solute+2)
	crg3 = crg(nat_solute+3)

	A_11 = iaclib(iac1)%avdw(ljcod(iac1, iac1)) &
		*iaclib(iac1)%avdw(ljcod(iac1, iac1))
	B_11 = iaclib(iac1)%bvdw(ljcod(iac1, iac1)) &
		*iaclib(iac1)%bvdw(ljcod(iac1, iac1))
	A_12 = iaclib(iac1)%avdw(ljcod(iac1, iac2)) &
		*iaclib(iac2)%avdw(ljcod(iac1, iac2))
	B_12 = iaclib(iac1)%bvdw(ljcod(iac1, iac2)) &
		*iaclib(iac2)%bvdw(ljcod(iac1, iac2))
	A_13 = iaclib(iac1)%avdw(ljcod(iac1, iac3)) &
		*iaclib(iac3)%avdw(ljcod(iac1, iac3))
	B_13 = iaclib(iac1)%bvdw(ljcod(iac1, iac3)) &
		*iaclib(iac3)%bvdw(ljcod(iac1, iac3))
	A_22 = iaclib(iac2)%avdw(ljcod(iac2, iac2)) &
		*iaclib(iac2)%avdw(ljcod(iac2, iac2))
	B_22 = iaclib(iac2)%bvdw(ljcod(iac2, iac2)) &
		*iaclib(iac2)%bvdw(ljcod(iac2, iac2))
	A_23 = iaclib(iac2)%avdw(ljcod(iac2, iac3)) &
		*iaclib(iac3)%avdw(ljcod(iac2, iac3))
	B_23 = iaclib(iac2)%bvdw(ljcod(iac2, iac3)) &
		*iaclib(iac3)%bvdw(ljcod(iac2, iac3))
	A_33 = iaclib(iac3)%avdw(ljcod(iac3, iac3)) &
		*iaclib(iac3)%avdw(ljcod(iac3, iac3))
	B_33 = iaclib(iac3)%bvdw(ljcod(iac3, iac3)) &
		*iaclib(iac3)%bvdw(ljcod(iac3, iac3))
  end if

  ipstart = 1

  do iw = calculation_assignment%ww%start, calculation_assignment%ww%end
	ip = ipstart

	do while (nbww(ip) /= 0)
	  ! consider the pair (nat_solute+3*(iw-1)) - nbww(ip)

	  ! --- 1 - (1,2,3) ---
	  i3 = (nat_solute+3*(iw-1))*3+1 !point to x for atom 1 in mol. iw
										!corresponds to oxygen in water
	  j3 = nbww(ip)*3-2 !point to 1 in interacting mol.
						!corresponds to oxygen in water	
	  ! 1 - 1
	  dx1X = x((j3  ))-x((i3  ))
	  dy1X = x((j3+1))-x((i3+1))
	  dz1X = x((j3+2))-x((i3+2))
	  !the periodical shift
	  boxshiftx = boxlength(1)*nint( dx1X*inv_boxl(1) )
	  boxshifty = boxlength(2)*nint( dy1X*inv_boxl(2) )
	  boxshiftz = boxlength(3)*nint( dz1X*inv_boxl(3) )
	  ! 1-2
	  j3 = j3 + 3
	  dx2X = x((j3  ))-x((i3  ))
	  dy2X = x((j3+1))-x((i3+1))
	  dz2X = x((j3+2))-x((i3+2))
	  ! 1-3
	  j3 = j3 + 3
	  dx3X = x((j3  ))-x((i3  ))
	  dy3X = x((j3+1))-x((i3+1))
	  dz3X = x((j3+2))-x((i3+2))
	  dx1X = dx1X - boxshiftx
	  dy1X = dy1X - boxshifty
	  dz1X = dz1X - boxshiftz
	  dx2X = dx2X - boxshiftx
	  dy2X = dy2X - boxshifty
	  dz2X = dz2X - boxshiftz
	  dx3X = dx3X - boxshiftx
	  dy3X = dy3X - boxshifty
	  dz3X = dz3X - boxshiftz

	  r1X = dx1X*dx1X + dy1X*dy1X + dz1X*dz1X
	  r2X = dx2X*dx2X + dy2X*dy2X + dz2X*dz2X
	  r3X = dx3X*dx3X + dy3X*dy3X + dz3X*dz3X

	  r1X=	sqrt(1._8/r1X  )
	  r2X=	sqrt(1._8/r2X )
	  r3X=	sqrt(1._8/r3X )
	  ! 1-1 
	  r2 = r1X * r1X
	  Vel = crg1*crg1*r1X
	  V_a  = A_11*(r2*r2*r2)*(r2*r2*r2)
	  V_b  = B_11*(r2*r2*r2)
	  E%ww%vdw = E%ww%vdw + V_a-V_B
	  E%ww%el = E%ww%el + Vel
	  dv = r2*( -Vel -12.*V_a +6.*V_b )
	  j3 = j3 - 6 !move pointer back to 1 in interacting mol.
	  d((i3  )) = d((i3  )) -(dv*dx1X)
	  d((j3  )) = d((j3  )) +(dv*dx1X)
	  d((i3+1)) = d((i3+1)) -(dv*dy1X)
	  d((j3+1)) = d((j3+1)) +(dv*dy1X)
	  d((i3+2)) = d((i3+2)) -(dv*dz1X)
	  d((j3+2)) = d((j3+2)) +(dv*dz1X)
	  ! 1-2 
	  r2 = r2X * r2X
	  Vel = crg1*crg2*r2X
	  V_a  = A_12*(r2*r2*r2)*(r2*r2*r2)
	  V_b  = B_12*(r2*r2*r2)
	  E%ww%vdw = E%ww%vdw + V_a-V_B
	  E%ww%el = E%ww%el + Vel
	  dv = r2*( -Vel -12.*V_a +6.*V_b )
	  j3 = j3 +3 !point to 2 in j-molecule
	  d((i3  )) = d((i3  )) -(dv*dx2X)
	  d((j3  )) = d((j3  )) +(dv*dx2X)
	  d((i3+1)) = d((i3+1)) -(dv*dy2X)
	  d((j3+1)) = d((j3+1)) +(dv*dy2X)
	  d((i3+2)) = d((i3+2)) -(dv*dz2X)
	  d((j3+2)) = d((j3+2)) +(dv*dz2X)
	  ! 1-3
	  r2 = r3X * r3X
	  Vel = crg1*crg3*r3X
	  V_a  = A_13*(r2*r2*r2)*(r2*r2*r2)
	  V_b  = B_13*(r2*r2*r2)
	  E%ww%vdw = E%ww%vdw + V_a-V_B
	  E%ww%el = E%ww%el + Vel
	  dv = r2*( -Vel -12.*V_a +6.*V_b )
	  j3 = j3 +3 !point to 3 in j-molecule
	  d((i3  )) = d((i3  )) -(dv*dx3X)
	  d((j3  )) = d((j3  )) +(dv*dx3X)
	  d((i3+1)) = d((i3+1)) -(dv*dy3X)
	  d((j3+1)) = d((j3+1)) +(dv*dy3X)
	  d((i3+2)) = d((i3+2)) -(dv*dz3X)
	  d((j3+2)) = d((j3+2)) +(dv*dz3X)
			
	  ! --- 2 - (1,2,3) ---
	  i3 = i3 + 3 !point to x for 2 in mol. iw
	  j3 = nbww(ip)*3-2 !point to o in j-mol.
	  ! 2 - 1 (X=1)
	  dx1X = x((j3  ))-x((i3  ))
	  dy1X = x((j3+1))-x((i3+1))
	  dz1X = x((j3+2))-x((i3+2))
	  ! 2-2 (X=2)
	  j3 = j3 + 3
	  dx2X = x((j3  ))-x((i3  ))
	  dy2X = x((j3+1))-x((i3+1))
	  dz2X = x((j3+2))-x((i3+2))
	  ! 2-3 (X=3)
	  j3 = j3 + 3
	  dx3X = x((j3  ))-x((i3  ))
	  dy3X = x((j3+1))-x((i3+1))
	  dz3X = x((j3+2))-x((i3+2))
	  dx1X = dx1X - boxshiftx
	  dy1X = dy1X - boxshifty
	  dz1X = dz1X - boxshiftz
	  dx2X = dx2X - boxshiftx
	  dy2X = dy2X - boxshifty
	  dz2X = dz2X - boxshiftz
	  dx3X = dx3X - boxshiftx
	  dy3X = dy3X - boxshifty
	  dz3X = dz3X - boxshiftz

	  r1X = dx1X*dx1X + dy1X*dy1X + dz1X*dz1X
	  r2X = dx2X*dx2X + dy2X*dy2X + dz2X*dz2X
	  r3X = dx3X*dx3X + dy3X*dy3X + dz3X*dz3X
	  r1X=	sqrt(1._8/r1X  )
	  r2X=	sqrt(1._8/r2X )
	  r3X=	sqrt(1._8/r3X )
	  ! 2-1 
	  r2 = r1X * r1X
	  Vel = crg2*crg1*r1X
	  V_a  = A_12*(r2*r2*r2)*(r2*r2*r2)
	  V_b  = B_12*(r2*r2*r2)
	  E%ww%vdw = E%ww%vdw + V_a-V_B
	  E%ww%el = E%ww%el + Vel
	  dv = r2*( -Vel -12.*V_a +6.*V_b )
	  j3 = j3 - 6 !move pointer back to o in j-mol.
	  d((i3  )) = d((i3  )) -(dv*dx1X)
	  d((j3  )) = d((j3  )) +(dv*dx1X)
	  d((i3+1)) = d((i3+1)) -(dv*dy1X)
	  d((j3+1)) = d((j3+1)) +(dv*dy1X)
	  d((i3+2)) = d((i3+2)) -(dv*dz1X)
	  d((j3+2)) = d((j3+2)) +(dv*dz1X)
	  ! 2-2 
	  r2 = r2X * r2X
	  Vel = crg2*crg2*r2X
	  V_a  = A_22*(r2*r2*r2)*(r2*r2*r2)
	  V_b  = B_22*(r2*r2*r2)
	  E%ww%vdw = E%ww%vdw + V_a-V_B
	  E%ww%el = E%ww%el + Vel
	  dv = r2*( -Vel -12.*V_a +6.*V_b )
	  j3 = j3 +3 !point to 2 in j-molecule
	  d((i3  )) = d((i3  )) -(dv*dx2X)
	  d((j3  )) = d((j3  )) +(dv*dx2X)
	  d((i3+1)) = d((i3+1)) -(dv*dy2X)
	  d((j3+1)) = d((j3+1)) +(dv*dy2X)
	  d((i3+2)) = d((i3+2)) -(dv*dz2X)
	  d((j3+2)) = d((j3+2)) +(dv*dz2X)
	  ! 2-3
	  r2 = r3X * r3X
	  Vel = crg2*crg3*r3X
	  V_a  = A_23*(r2*r2*r2)*(r2*r2*r2)
	  V_b  = B_23*(r2*r2*r2)
	  E%ww%vdw = E%ww%vdw + V_a-V_B
	  E%ww%el = E%ww%el + Vel
	  dv = r2*( -Vel -12.*V_a +6.*V_b )
	  j3 = j3 +3 !point to 3 in j-molecule
	  d((i3  )) = d((i3  )) -(dv*dx3X)
	  d((j3  )) = d((j3  )) +(dv*dx3X)
	  d((i3+1)) = d((i3+1)) -(dv*dy3X)
	  d((j3+1)) = d((j3+1)) +(dv*dy3X)
	  d((i3+2)) = d((i3+2)) -(dv*dz3X)
	  d((j3+2)) = d((j3+2)) +(dv*dz3X)
			
	  ! --- 3 - (1,2,3) ---
	  i3 = i3 + 3 !point to x for 3 in mol. iw
	  j3 = nbww(ip)*3-2 !point to o in j-mol.
	  ! 3 - 1 (X=1)
	  dx1X = x((j3  ))-x((i3  ))
	  dy1X = x((j3+1))-x((i3+1))
	  dz1X = x((j3+2))-x((i3+2))
	  ! 3-2 (X=2)
	  j3 = j3 + 3
	  dx2X = x((j3  ))-x((i3  ))
	  dy2X = x((j3+1))-x((i3+1))
	  dz2X = x((j3+2))-x((i3+2))
	  ! 3-3 (X=3)
	  j3 = j3 + 3
	  dx3X = x((j3  ))-x((i3  ))
	  dy3X = x((j3+1))-x((i3+1))
	  dz3X = x((j3+2))-x((i3+2))
	  dx1X = dx1X - boxshiftx
	  dy1X = dy1X - boxshifty
	  dz1X = dz1X - boxshiftz
	  dx2X = dx2X - boxshiftx
	  dy2X = dy2X - boxshifty
	  dz2X = dz2X - boxshiftz
	  dx3X = dx3X - boxshiftx
	  dy3X = dy3X - boxshifty
	  dz3X = dz3X - boxshiftz

	  r1X = dx1X*dx1X + dy1X*dy1X + dz1X*dz1X
	  r2X = dx2X*dx2X + dy2X*dy2X + dz2X*dz2X
	  r3X = dx3X*dx3X + dy3X*dy3X + dz3X*dz3X
	  r1X=	sqrt(1._8/r1X  )
	  r2X=	sqrt(1._8/r2X )
	  r3X=	sqrt(1._8/r3X )
	  ! 3-1 
	  r2 = r1X * r1X
	  Vel = crg3*crg1*r1X
	  V_a  = A_13*(r2*r2*r2)*(r2*r2*r2)
	  V_b  = B_13*(r2*r2*r2)
	  E%ww%vdw = E%ww%vdw + V_a-V_B
	  E%ww%el = E%ww%el + Vel
	  dv = r2*( -Vel -12.*V_a +6.*V_b )
	  j3 = j3 - 6 !move pointer back to o in j-mol.
	  d((i3  )) = d((i3  )) -(dv*dx1X)
	  d((j3  )) = d((j3  )) +(dv*dx1X)
	  d((i3+1)) = d((i3+1)) -(dv*dy1X)
	  d((j3+1)) = d((j3+1)) +(dv*dy1X)
	  d((i3+2)) = d((i3+2)) -(dv*dz1X)
	  d((j3+2)) = d((j3+2)) +(dv*dz1X)
	  ! 3-2 
	  r2 = r2X * r2X
	  Vel = crg3*crg2*r2X
	  V_a  = A_23*(r2*r2*r2)*(r2*r2*r2)
	  V_b  = B_23*(r2*r2*r2)
	  E%ww%vdw = E%ww%vdw + V_a-V_B
	  E%ww%el = E%ww%el + Vel
	  dv = r2*( -Vel -12.*V_a +6.*V_b )
	  j3 = j3 +3 !point to 2 in j-molecule
	  d((i3  )) = d((i3  )) -(dv*dx2X)
	  d((j3  )) = d((j3  )) +(dv*dx2X)
	  d((i3+1)) = d((i3+1)) -(dv*dy2X)
	  d((j3+1)) = d((j3+1)) +(dv*dy2X)
	  d((i3+2)) = d((i3+2)) -(dv*dz2X)
	  d((j3+2)) = d((j3+2)) +(dv*dz2X)
	  ! 3-3
	  r2 = r3X * r3X
	  Vel = crg3*crg3*r3X
	  V_a  = A_33*(r2*r2*r2)*(r2*r2*r2)
	  V_b  = B_33*(r2*r2*r2)
	  E%ww%vdw = E%ww%vdw + V_a-V_B
	  E%ww%el = E%ww%el + Vel
	  dv = r2*( -Vel -12.*V_a +6.*V_b )
	  j3 = j3 +3 !point to 3 in j-molecule
	  d((i3  )) = d((i3  )) -(dv*dx3X)
	  d((j3  )) = d((j3  )) +(dv*dx3X)
	  d((i3+1)) = d((i3+1)) -(dv*dy3X)
	  d((j3+1)) = d((j3+1)) +(dv*dy3X)
	  d((i3+2)) = d((i3+2)) -(dv*dz3X)
	  d((j3+2)) = d((j3+2)) +(dv*dz3X)

	  ip = ip + 1
	end do ! while ip

	! skip the gap
	ipstart = ip + 1			
  end do ! iw

end subroutine nonbond_3atomsolvent_box

!-----------------------------------------------------------------------

subroutine offdiag
! local variables
integer						:: io,i,j,k,l,k3,l3
real(8)						:: r

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
real(8) :: fk,r2,erst,Edum,x2,y2,z2,wgt,b,db,dv, totmass, theta,rij,r2ij,rjk,r2jk, scp,f1
real(8) :: dr(3), dr2(3), ctr(3), di(3), dk(3)
real(8)						::	fexp

! global variables used:
!  E, nstates, EQ, nrstr_seq, rstseq, heavy, x, xtop, d, nrstr_pos, rstpos, nrstr_dist, 
!  rstdis, nrstr_wall, rstang, nrstr_ang, rstwal, xwcent, excl, freeze

! sequence restraints (independent of Q-state)
do ir = 1, nrstr_seq
fk = rstseq(ir)%fk

if(rstseq(ir)%to_centre == 1) then     ! Put == 1, then equal to 2
  ! restrain to geometrical centre

  ! reset dr & atom counter
  dr(:) = 0.
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
        erst    = 0.5*fk*r2
        E%restraint%protein  = E%restraint%protein + erst

        ! apply same force to all atoms
        do i = rstseq(ir)%i, rstseq(ir)%j
        if ( heavy(i) .or. rstseq(ir)%ih .eq. 1 .and.(.not.(excl(i).and.freeze))) then
                d(i*3-2:i*3) = d(i*3-2:i*3) + fk*dr(:)*iaclib(iac(i))%mass/12.010
          end if
        end do
  end if
 
else if(rstseq(ir)%to_centre == 2) then     ! Put == 1, then equal to 2
  ! restrain to mass centre
  ! reset dr & variable to put masses
  dr(:) = 0.
  totmass = 0.
  
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
        erst    = 0.5*fk*r2
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

          erst    = 0.5*fk*r2
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
wgt = 1.0
end if

x2 = dr(1)**2
y2 = dr(2)**2
z2 = dr(3)**2

Edum = 0.5*rstpos(ir)%fk(1)*x2 + &
   0.5*rstpos(ir)%fk(2)*y2 + &
   0.5*rstpos(ir)%fk(3)*z2

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
wgt = 1.0
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

Edum   = 0.5*rstdis(ir)%fk*db**2
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
wgt = 1.0
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
if ( scp .gt.  1.0 ) scp =  1.0
if ( scp .lt. -1.0 ) scp = -1.0

theta = acos(scp)
db    = theta - (rstang(ir)%ang)*deg2rad

! dv is the force to be added in module
Edum   = 0.5*rstang(ir)%fk*db**2
dv     = wgt*rstang(ir)%fk*db

! calculate sin(theta) to use in forces
f1 = sin ( theta )
if ( abs(f1) .lt. 1.e-12 ) then
        ! avoid division by zero
        f1 = -1.e12
else
        f1 =  -1.0 / f1
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

                        if(db > 0.) then
                                erst =  0.5 * fk * db**2 - rstwal(ir)%Dmorse
                                dv = fk*db/b
                        else
                                fexp = exp(rstwal(ir)%aMorse*db)
                                erst = rstwal(ir)%dMorse*(fexp*fexp-2.*fexp)
                                dv=-2.*rstwal(ir)%dMorse*rstwal(ir)%aMorse*(fexp-fexp*fexp)/b
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

subroutine pot_energy
! local variables


integer					:: istate, i, nat3, numrequest
integer					:: is, j
#if defined (PROFILING)
real(8)					:: start_loop_time1
real(8)					:: start_loop_time2
real(8)					:: start_loop_time3
#endif

! --- reset all energies

E%potential = 0.0
!E%kinetic = 0.0	! no need to reset because it will be assigned its final value at once
E%LRF = 0.0
E%p%bond  = 0.0
E%p%angle = 0.0
E%p%torsion = 0.0
E%p%improper = 0.0
E%w%bond  = 0.0
E%w%angle = 0.0
E%w%torsion = 0.0
E%w%improper = 0.0
E%q%bond  = 0.0
E%q%angle = 0.0
E%q%torsion   = 0.0
E%q%improper   = 0.0
E%pp%el  = 0.0
E%pp%vdw = 0.0
E%pw%el  = 0.0
E%pw%vdw = 0.0
E%ww%el  = 0.0
E%ww%vdw = 0.0
E%qx%el    = 0.0
E%qx%vdw   = 0.0
!E%restraint%total = 0.0	! will be assigned its final value later
E%restraint%fix = 0.0
E%restraint%shell = 0.0
E%restraint%protein = 0.0
E%restraint%solvent_radial = 0.0
E%restraint%water_pol = 0.0
do istate = 1, nstates
!EQ(istate)%lambda set by initialize
!EQ(istate)%total assigned its final value later
EQ(istate)%q%bond = 0.0
EQ(istate)%q%angle = 0.0
EQ(istate)%q%torsion = 0.0
EQ(istate)%q%improper = 0.0
!EQ(istate)%qx%el = 0.0	! assigned its final value later
!EQ(istate)%qx%vdw = 0.0	! assigned its final value later
EQ(istate)%qq%el = 0.0
EQ(istate)%qq%vdw = 0.0
EQ(istate)%qp%el = 0.0
EQ(istate)%qp%vdw = 0.0
EQ(istate)%qw%el = 0.0
EQ(istate)%qw%vdw = 0.0
EQ(istate)%restraint = 0.0
if (use_excluded_groups) then
        do j=1,ngroups_gc
        EQ_gc(j,istate)%qq%el = 0.0
        EQ_gc(j,istate)%qq%vdw = 0.0
        EQ_gc(j,istate)%qp%el = 0.0
        EQ_gc(j,istate)%qp%vdw = 0.0
        end do
end if
end do

!reset derivatives ---
d(:) = 0.

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
if(.not. qq_use_library_charges) then
        if(ivdw_rule .eq. 1 ) then
                call nonbond_qq
        elseif ( ivdw_rule .eq. 2 ) then
                call nonbon2_qq
        end if
else
        if ( ivdw_rule .eq. 1 ) then
                call nonbond_qq_lib_charges
        else if ( ivdw_rule .eq. 2 ) then
                call nonbon2_qq_lib_charges
        end if
end if

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
if (use_excluded_groups) then
        numrequest = 4
else
        numrequest = 3
end if
do i = 1, numrequest
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
  EQ(1:nstates)%qp%el  = EQ(1:nstates)%qp%el  + EQ_recv(1:nstates,i)%qp%el
  EQ(1:nstates)%qp%vdw = EQ(1:nstates)%qp%vdw + EQ_recv(1:nstates,i)%qp%vdw
  EQ(1:nstates)%qw%el  = EQ(1:nstates)%qw%el  + EQ_recv(1:nstates,i)%qw%el
  EQ(1:nstates)%qw%vdw = EQ(1:nstates)%qw%vdw + EQ_recv(1:nstates,i)%qw%vdw
if (use_excluded_groups) then
!	do j=1,ngroups_gc
	EQ_gc(1:ngroups_gc,1:nstates)%qp%el = EQ_gc(1:ngroups_gc,1:nstates)%qp%el + &
						EQ_gc_recv(1:ngroups_gc,1:nstates,i)%qp%el
	EQ_gc(1:ngroups_gc,1:nstates)%qp%vdw = EQ_gc(1:ngroups_gc,1:nstates)%qp%vdw + &
						EQ_gc_recv(1:ngroups_gc,1:nstates,i)%qp%vdw
!	end do
end if
end do
#endif
if (use_excluded_groups) then
	do j=1,ngroups_gc
	EQ_gc(j,1:nstates)%qw%el = EQ(1:nstates)%qw%el
	EQ_gc(j,1:nstates)%qw%vdw = EQ(1:nstates)%qw%vdw
	EQ_gc(j,1:nstates)%q = EQ(1:nstates)%q
	EQ_gc(j,1:nstates)%restraint = EQ(1:nstates)%restraint
	do istate=1,nstates
	EQ_gc(j,istate)%qx%el = EQ_gc(j,istate)%qq%el + EQ_gc(j,istate)%qp%el + EQ_gc(j,istate)%qw%el
	EQ_gc(j,istate)%qx%vdw = EQ_gc(j,istate)%qq%vdw + EQ_gc(j,istate)%qp%vdw + EQ_gc(j,istate)%qw%vdw
	EQ_gc(j,istate)%total = EQ_gc(j,istate)%q%bond + EQ_gc(j,istate)%q%angle &
 + EQ_gc(j,istate)%q%torsion + EQ_gc(j,istate)%q%improper + EQ_gc(j,istate)%qx%el &
 + EQ_gc(j,istate)%qx%vdw + EQ_gc(j,istate)%restraint
	end do
	end do
end if

! q-atom energy summary
do istate = 1, nstates
! update EQ
EQ(istate)%qx%el  = EQ(istate)%qq%el +EQ(istate)%qp%el +EQ(istate)%qw%el
EQ(istate)%qx%vdw = EQ(istate)%qq%vdw+EQ(istate)%qp%vdw+EQ(istate)%qw%vdw

EQ(istate)%total =  EQ(istate)%q%bond + EQ(istate)%q%angle   &
  + EQ(istate)%q%torsion  + EQ(istate)%q%improper + EQ(istate)%qx%el &
  + EQ(istate)%qx%vdw  + EQ(istate)%restraint

! update E with an average of all states
E%q%bond  = E%q%bond  + EQ(istate)%q%bond *EQ(istate)%lambda
E%q%angle = E%q%angle + EQ(istate)%q%angle*EQ(istate)%lambda
E%q%torsion   = E%q%torsion   + EQ(istate)%q%torsion  *EQ(istate)%lambda
E%q%improper   = E%q%improper   + EQ(istate)%q%improper  *EQ(istate)%lambda
E%qx%el    = E%qx%el    + EQ(istate)%qx%el   *EQ(istate)%lambda
E%qx%vdw   = E%qx%vdw   + EQ(istate)%qx%vdw  *EQ(istate)%lambda

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

if( use_PBC ) then !periodic box
      
		select case(ivdw_rule)
        case(VDW_GEOMETRIC)
                call nonbond_pp_box
                call nonbond_pw_box
                if(qvdw_flag) then
                        call nonbond_qp_qvdw_box
                else
                        call nonbond_qp_box
                end if
                if(natom > nat_solute) then !if any solvent
                        if(solvent_type == SOLVENT_SPC) then
                                !use the optimised SPC routine when possible
                                call nonbond_ww_spc_box
                                call nonbond_qw_spc_box
                        elseif(solvent_type == SOLVENT_3ATOM) then !otherwise calc. LJ with all atoms
                                call nonbond_3atomsolvent_box
                                call nonbond_qw_3atom_box
                        end if
                end if
        case(VDW_ARITHMETIC)
                call nonbon2_pp_box
                call nonbon2_qp_box
                if(natom > nat_solute) then !if any solvent
                        call nonbon2_pw_box
                        call nonbon2_qw_box !no SPC-specific optimised routines here
                        call nonbon2_ww_box
                end if
        end select
        
		!LRF PBC
		if (use_LRF) then
                call lrf_taylor
        end if

else !simulation sphere

        select case(ivdw_rule)
        case(VDW_GEOMETRIC)
                call nonbond_pp   
                call nonbond_pw   
                if(qvdw_flag) then
                        call nonbond_qp_qvdw
                else
                        call nonbond_qp
                end if
                if(natom > nat_solute) then !if any solvent
                        if(solvent_type == SOLVENT_SPC) then
                                !use the optimised SPC routine when possible
                                call nonbond_ww_spc
                                call nonbond_qw_spc
                        elseif(solvent_type == SOLVENT_3ATOM) then !otherwise calc. LJ with all atoms
                                call nonbond_3atomsolvent
                                call nonbond_qw_3atom
                        end if
                end if
        case(VDW_ARITHMETIC)
                call nonbon2_pp
                call nonbon2_qp
                if(natom > nat_solute) then !if any solvent
                        call nonbon2_pw
                        call nonbon2_qw !no SPC-specific optimised routines here
                        call nonbon2_ww
                end if	
        end select

        ! on demand: taylor expansion of the electric field from charge groups beyond rcutoff
        if (use_LRF) then
                call lrf_taylor
        end if

end if

end subroutine pot_energy_nonbonds

!-----------------------------------------------------------------------

subroutine prep_coord


! local variables
integer(4)          :: i,nat3

! --- Refresh topology coords. if needed (external restraints file)
if ( implicit_rstr_from_file .eq. 1 ) then
write (*,'(/,a,/)') 'Refreshing topology coords for restraining...'
read (12) nat3,(xtop(i),i=1,nat_pro*3)
end if

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
        read (2) nat3
        rewind(2)
        if(nat3 /= 3*natom) then
                write(*,100) nat3/3, natom
100			format('>>>>> ERROR:',i5,' atoms in restart file not equal to',i5,&
                        ' in topology.')
                call die('wrong number of atoms in restart file')
        end if
        read (2,err=112,end=112) nat3, (x(i),i=1,nat3)
        read (2,err=112,end=112) nat3, (v(i),i=1,nat3)
        write (*,'(a30,i8)')   'Total number of atoms        =',natom
        write (*,'(a30,i8,/)') 'Number of waters encountered =',nwat

        if( use_PBC) then
                read(2,err=112,end=112) boxlength(:)
                read(2,err=112,end=112) boxcentre(:)
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
	real(8)						::	rin2,r2

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
	real(8)						::	rout2,rin2,r2
	real(8), allocatable		::	cgp_cent(:,:)
	nshellats = 0
	rin2  = rexcl_i**2

	shell(:) = .false.

	allocate(cgp_cent(3,ncgp+nwat))

	cgp_cent(:,:) = 0.

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
subroutine prep_sim
! local variables
integer						:: i, j, ig, istate

if (nodeid .eq. 0) then	
        write(*,*)
        call centered_heading('Initialising dynamics', '-')
end if

! Set parameters (bonds, angles, charges,...) & restraints for water   
if(nwat > 0) then
        select case (solvent_type)
        case (SOLVENT_SPC, SOLVENT_3ATOM)
                crg_ow = crg(nat_solute+1)
                crg_hw = -crg_ow/2.0
        case(SOLVENT_GENERAL)
                !add appropriate code here
                call die('Topology contains mixed or non-3-atomic solvent. This feature is not implemented yet.')
        end select

        if( .not. use_PBC ) then
                call wat_sphere
                if (wpol_restr) call wat_shells

        else !compute charges of the system for box case 
        !(done in subroutine wat_sphere for sphere case)
                !calc. total charge of non-Q-atoms
                crgtot = 0.0
                do i = 1, nat_solute
                        if ( iqatom(i)==0 ) crgtot = crgtot + crg(i)
                end do
                write (*,60) crgtot
60 format ('Total charge of non-Q atoms             = ',f10.2)


        !calc effective charge of whole system at this lambda
                crgQtot = 0.0
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
end if

!	Prepare an array of inverse masses
winv(:) = 1./iaclib(iac(:))%mass


if(use_PBC .and. control_box) then
        boxlength(:) = new_boxl(:)
	if ( put_solute_back_in_box .or. put_solvent_back_in_box ) then !only put back in box if either solute or solvent should be put back (qdyn input option)
        	call put_back_in_box
	end if
        write(*,'(a)') 'Boxsize changed. Equilibration may be needed'
end if

if( use_PBC ) then 
        !compute masses of all molecules
        allocate(mol_mass(nmol))
        mol_mass(:) = 0.0

        do i = 1,nmol-1 !all molecules but the last
                do j = istart_mol(i), istart_mol(i+1)-1 !all atoms of molecule
                        mol_mass(i) = mol_mass(i) + iaclib(iac(j))%mass
                end do
        end do

        do j = istart_mol(nmol), natom !last molecule
                mol_mass(nmol) = mol_mass(nmol) + iaclib(iac(j))%mass
        end do

        mol_mass(:) = 1./mol_mass(:)

        !prepare array of masses
        allocate(mass(natom))
!make this a duplicate of the iaclib array to avoid truncation errors
!        mass(:) = 1.0/winv(:)
        mass(:)=iaclib(iac(:))%mass

end if

!initialization for the Nosé-Hoover chain thermostat
if ( thermostat == NOSEHOOVER ) then
	do i=1,numchain
	   qnh(i)=nhq
	   vnh(i)=0
	   xnh(i)=0
	end do
	shake_constraints = sum(shake_mol(1:shake_molecules)%nconstraints)
	Ndegf=3*natom-shake_constraints	
	qnh(1)=Ndegf*nhq
end if

!scale charges by sqrt(coulomb_constant) 
crg(:) = crg(:) * sqrt(coulomb_constant)
crg_ow = crg_ow * sqrt(coulomb_constant)
crg_hw = crg_hw * sqrt(coulomb_constant)

if(nqat > 0) then
        qcrg(:,:) = qcrg(:,:) * sqrt(coulomb_constant)
end if
if (use_excluded_groups) then
        do i=1,ngroups_gc
                call mask_initialize(ST_gc(i)%gcmask)
                call gc_make_array(ST_gc(i))
                do j=1,ST_gc(i)%count
                stat = mask_add(ST_gc(i)%gcmask,ST_gc(i)%sendtomask(j))
                end do
                do j=1,nstates
                EQ_gc(i,j)%lambda=EQ(j)%lambda
                end do

        end do
end if
end subroutine prep_sim

!-----------------------------------------------------------------------

subroutine qangle (istate)
! arguments
integer						:: istate

! local variables
integer						:: ia,i,j,k,ic,i3,j3,k3,im,icoupl,ib
real(8)						:: bji,bjk,scp,ang,da,ae,dv,gamma
real(8)						:: rji(3),rjk(3),f1,di(3),dk(3)


do ia = 1, nqangle


ic = qang(ia)%cod(istate)
 !skip if angle not present (code 0)
if ( ic > 0 ) then

gamma = 1.0
icoupl = 0

do im = 1, nang_coupl
   if ( iang_coupl(1,im) .eq. ia ) then
      icoupl = im
      ib     = iang_coupl(2,im)
      gamma = EMorseD(ib)
	  !couple improper to bond breaking not making
	  if ( iang_coupl(3,im) .eq. 1) gamma = 1.0_8 - gamma
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
if ( scp .gt.  1.0 ) scp =  1.0
if ( scp .lt. -1.0 ) scp = -1.0
ang = acos(scp)
da = ang - qanglib(ic)%ang0
ae = 0.5*qanglib(ic)%fk*da**2
EQ(istate)%q%angle = EQ(istate)%q%angle + ae*gamma

dv = gamma*qanglib(ic)%fk*da*EQ(istate)%lambda
f1 = sin ( ang )
if ( abs(f1) .lt. 1.e-12 ) f1 = 1.e-12
f1 =  -1.0 / f1
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
real(8)						::	gamma
real(8)						::	rik(3), dik, du, ru, Eurey

do ia = 1, nqangle
        ic = qang(ia)%cod(istate)
        !skip if angle not present (code 0)
        if ( ic == 0  .or. qanglib(ic)%ureyfk == 0.) cycle

gamma = 1.0
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
real(8)						::	b,db,be,dv,fexp
real(8)						::	rij(3)

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
        be = qbondlib(ic)%Dmz*(fexp*fexp-2.*fexp) + 0.5*qbondlib(ic)%fk*db**2
        EMorseD(ib) = -(fexp*fexp-2.*fexp)
EQ(istate)%q%bond = EQ(istate)%q%bond + be
dv = (2.*qbondlib(ic)%Dmz*qbondlib(ic)%amz*(fexp-fexp*fexp) + qbondlib(ic)%fk*db)*EQ(istate)%lambda/b

d(i3+1) = d(i3+1) - dv*rij(1)
d(i3+2) = d(i3+2) - dv*rij(2)
d(i3+3) = d(i3+3) - dv*rij(3)
d(j3+1) = d(j3+1) + dv*rij(1)
d(j3+2) = d(j3+2) + dv*rij(2)
d(j3+3) = d(j3+3) + dv*rij(3)

!Force scaling factor to be 1 when distance is smaller than r0
 if ( db > 0 ) then
 	EMorseD(ib) = -(fexp*fexp-2.*fexp)
 	dMorse_i(:,ib) = +2.*qbondlib(ic)%amz*(fexp-fexp*fexp)*EQ(istate)%lambda/b*rij(:)
 	dMorse_j(:,ib) = -2.*qbondlib(ic)%amz*(fexp-fexp*fexp)*EQ(istate)%lambda/b*rij(:)
 else
 	EMorseD(ib) = 1
 	dMorse_i(:,ib) = 0
 	dMorse_j(:,ib) = 0
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
real(8)						::	bj,bk,scp,phi,sgn,pe,dv,arg,f1,gamma
real(8)						::	rji(3),rjk(3),rkl(3),rnj(3),rnk(3)
real(8)						::	rki(3),rlj(3),dp(12),di(3),dl(3)

do ip = 1,nqimp

ic = qimp(ip)%cod(istate)
!ic = qimpcod(ip,istate)

if ( ic > 0 ) then

gamma = 1.0
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
!i  = iqimp(ip)
!j  = jqimp(ip)
!k  = kqimp(ip)
!l  = lqimp(ip)

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
if ( scp .gt.  1.0 ) scp =  1.0
if ( scp .lt. -1.0 ) scp = -1.0
phi = acos ( scp )
sgn =  rjk(1)*(rnj(2)*rnk(3)-rnj(3)*rnk(2)) &
     +rjk(2)*(rnj(3)*rnk(1)-rnj(1)*rnk(3)) &
     +rjk(3)*(rnj(1)*rnk(2)-rnj(2)*rnk(1))
if ( sgn .lt. 0 ) phi = -phi

! ---       energy

arg = phi - qimplib(ic)%imp0
!arg = phi - qimp0(ic)
arg = arg - 2.*pi*nint(arg/(2.*pi))
dv = qimplib(ic)%fk*arg
!dv  = qfkimp(ic)*arg
pe  = 0.5*dv*arg
EQ(istate)%q%improper = EQ(istate)%q%improper + pe*gamma
dv = dv*gamma*EQ(istate)%lambda

! ---       forces

f1 = sin ( phi ) 
if ( abs(f1) .lt. 1.e-12 ) f1 = 1.e-12
f1 =  -1.0 / f1
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
real(8)						::	bj,bk,scp,phi,sgn,pe,dv,arg,f1,gamma
real(8)						::	rji(3),rjk(3),rkl(3),rnj(3),rnk(3)
real(8)						::	rki(3),rlj(3),dp(12),di(3),dl(3)

do ip = 1,nqtor

ic = qtor(ip)%cod(istate)
!ic = qtorcod(ip,istate)

if ( ic > 0 ) then

gamma = 1.0
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
!i  = iqtor(ip)
!j  = jqtor(ip)
!k  = kqtor(ip)
!l  = lqtor(ip)

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
if ( scp .gt.  1.0 ) scp =  1.0
if ( scp .lt. -1.0 ) scp = -1.0
phi = acos ( scp )
sgn =  rjk(1)*(rnj(2)*rnk(3)-rnj(3)*rnk(2)) &
     +rjk(2)*(rnj(3)*rnk(1)-rnj(1)*rnk(3)) &


    +rjk(3)*(rnj(1)*rnk(2)-rnj(2)*rnk(1))
if ( sgn .lt. 0 ) phi = -phi




! ---       energy

arg = qtorlib(ic)%rmult*phi-qtorlib(ic)%deltor
pe = qtorlib(ic)%fk*(1.0+cos(arg))
!arg = qrmult(ic)*phi-qdeltor(ic)
!pe  = qfktor(ic)*(1.0+cos(arg))
EQ(istate)%q%torsion = EQ(istate)%q%torsion + pe*gamma
dv = -qtorlib(ic)%fk*sin(arg)*gamma*EQ(istate)%lambda
!dv = -qrmult(ic)*qfktor(ic)*sin(arg)*gamma*EQ(istate)%lambda

! ---       forces

f1 = sin ( phi ) 
if ( abs(f1) .lt. 1.e-12 ) f1 = 1.e-12
f1 =  -1.0 / f1
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


real(8) function randm (ig)
! arguments
integer					::	ig

! local variables
integer, parameter		::	m = 100000000
integer, parameter		::  m1 = 10000
integer, parameter		::	mult=31415821
integer					::	irandh,irandl,multh,multl
real(8)					::	r
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
r = real(irand / 10) * 10 / real(m)
if ((r .le. 0.e0) .or. (r .gt. 1.e0)) r = 0.e0
randm = r
ig = irand

end function randm

!-----------------------------------------------------------------------

integer function shake(xx, x)
!arguments
real(8)						::	xx(:), x(:)
!	returns no. of iterations

! *** local variables
integer						::	i,j,i3,j3,mol,ic,nits
real(8)						::	xij2,diff,corr,scp,xxij2
real(8)						::	xij(3),xxij(3)
#if defined (PROFILING)
real(8)                                         :: start_loop_time
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
        do !iteration loop
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
                                corr = diff/(2.*scp*(winv(i)+winv(j)))

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
                                        ! repeat until done:

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
                        call die('shake failure')
                end if
100         format ('>>> Shake failed, i,j,d,d0 = ',2i6,2f10.5)
        end do

        ! update niter
        shake = shake+nits
end do

! set niter to the average number of iterations per molecule
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
real(8)						::	rmass,totmass
real(8)						::	vcm(3)

! calculate totmass and vcm
totmass = 0.0
do i=1,3
vcm(i) = 0.0
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
real(8)					::	box_min, vtemp, vtemp1

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

nwat = (natom - nat_solute) / 3
!add extra element to molecule start atom array to keep track of last atom
istart_mol(nmol+1) = nat_pro + 1

! abort if no atoms
if (natom .eq. 0) call die('zero particles to simulate')

! convert libraries from degrees to radians
anglib(1:nangcod)%ang0 = deg2rad*anglib(1:nangcod)%ang0
torlib(1:ntorcod)%paths = 1.0/torlib(1:ntorcod)%paths
torlib(1:ntorcod)%deltor = deg2rad*torlib(1:ntorcod)%deltor
implib(1:nimpcod)%imp0 = deg2rad*implib(1:nimpcod)%imp0

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
if( any(boxlength(:) == 0.0 ) ) then
        inv_boxl(:) = 0.0
else
        inv_boxl(:) = 1.0/boxlength(:)
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

real(8) function torsion(istart, iend)
!arguments
integer						::	istart, iend

! local variables
integer						::	ip
real(8)						::	scp,phi,dv,arg,f1
real(8)						::	bjinv, bkinv, bj2inv, bk2inv
real(8), save					::	rji(3),rjk(3),rkl(3),rnj(3),rnk(3)
real(8), save					::	rki(3),rlj(3),dp(12),di(3),dl(3)
type(TOR_TYPE), pointer		::  t
type(TORLIB_TYPE), pointer	::	lib

! global variables used:
!  tor, torlib, x, d

! calculate the total energy of all torsion angles
! updates d

torsion = 0.

do ip = iend, istart,-1
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

        bj2inv = 1./(rnj(1)**2 + rnj(2)**2 + rnj(3)**2 )
bk2inv = 1./(rnk(1)**2 + rnk(2)**2 + rnk(3)**2 )
        bjinv = sqrt(bj2inv)
        bkinv = sqrt(bk2inv)

        ! calculate scp and phi
scp = (rnj(1)*rnk(1)+rnj(2)*rnk(2)+rnj(3)*rnk(3))*(bjinv*bkinv)
if ( scp .gt.  1.0 ) then
                scp =  1.0
                phi = acos (1.0) ! const
        else if ( scp .lt. -1.0 ) then
                scp = -1.0
                phi = acos (-1.0) ! const
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
torsion = torsion + lib%fk*(1.0+cos(arg))*lib%paths   !lib%paths is previously inverted 
dv = -lib%rmult*lib%fk*sin(arg)*lib%paths

! ---       forces

f1 = sin ( phi ) 
if ( abs(f1) .lt. 1.e-12 ) f1 = 1.e-12
f1 =  -1.0 / f1
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

end function torsion

!-----------------------------------------------------------------------
subroutine restrain_solvent 
! local variables
integer						::	iw,i,i3
real(8)						::	b,db,erst,dv,fexp
real(8), save					::	dr(3)
real(8)						::	shift

! global variables used:
!  E, Boltz, Tfree, fk_wsphere, nwat, nat_pro, x, xwcent, rwat, Dwmz, awmz, d

if(fk_wsphere /= 0.) then
        shift = sqrt (Boltz*Tfree/fk_wsphere)
else
        shift = 0.
end if
do iw = ncgp_solute + 1, ncgp
        i  = cgp(iw)%iswitch
        if (excl(i)) cycle ! skip excluded topology waters
        i3 = 3*i-3

        dr(1) = x(i3+1) - xwcent(1)
        dr(2) = x(i3+2) - xwcent(2)
        dr(3) = x(i3+3) - xwcent(3)
        b = sqrt ( dr(1)**2 + dr(2)**2 + dr(3)**2 )
        db = b - (rwat - shift)

        ! calculate erst and dv
        if ( db > 0 ) then
                erst = 0.5 * fk_wsphere * db**2 - Dwmz
                dv = fk_wsphere*db/b
        else
                if (b > 0.0) then
                  fexp = exp ( awmz*db )
                  erst = Dwmz*(fexp*fexp-2.*fexp)
                  dv = -2.*Dwmz*awmz*(fexp-fexp*fexp)/b
                else
                  dv = 0
                  erst = 0
                end if
        end if

        ! update energy and forces
        E%restraint%solvent_radial = E%restraint%solvent_radial + erst
        d(i3+1) = d(i3+1) + dv*dr(1)
        d(i3+2) = d(i3+2) + dv*dr(2)
        d(i3+3) = d(i3+3) + dv*dr(3)
end do
end subroutine restrain_solvent

!-----------------------------------------------------------------------

subroutine wat_sphere
! local variables
integer					:: i,i3,kr,isort,int_wat,istate
real(8)					:: rc,rnwat
real(8), save				:: dr(3)
real(8)					:: crgexcl

!possibly override target sphere radius from topology
if(rwat_in > 0.) rwat = rwat_in


!calc. total charge of non-excluded non-Q-atoms and excluded atoms
crgtot = 0.0
crgexcl = 0.0
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
crgQtot = 0.0
do i = 1, nqat
	do istate = 1, nstates
		crgtot = crgtot + qcrg(i,istate)*EQ(istate)%lambda
		crgQtot = crgQtot + qcrg(i,istate)*EQ(istate)%lambda
	end do
end do
write (*,70) crgtot
70 format ('Total charge of system                  = ',f10.2)

if (.not. wpol_born) crgQtot = 0.0 !ignore total Q if Born corr. is off
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
        Dwmz = 0.26*exp(-0.19*(rwat-15.))+0.74
end if
if(awmz == -1) then !use magic for the reach of the Morse potential
        awmz = 0.2/(1.+exp(0.4*(rwat-25.)))+0.3
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
real(8)						::	rout, dr, ri, Vshell, rshell, drs
integer						::	is, n_insh


integer						::	nwpolr_shell_restart, filestat
integer						::	bndcodw, angcodw

!calc mu_w
!look up bond code for water
bndcodw = bnd(nbonds)%cod
angcodw = ang(nangles)%cod
!find charge of water O = charge of 1st solvent atom
crg_ow = crg(nat_solute + 1)
mu_w = -crg_ow*bondlib(bndcodw)%bnd0*cos(anglib(angcodw)%ang0/2)

! shell widths are drout, 2drout, 3drout
drs = wpolr_layer / drout !number of drouts 

! calc number of shells based on arithmetic series sum formula
nwpolr_shell = int(-0.5 + sqrt(2*drs + 0.25)) 
allocate(wshell(nwpolr_shell), stat=alloc_status)
call check_alloc('water polarisation shell array')

write(*, 100) nwpolr_shell
100	format(/,'Setting up ', i1, ' water shells for polarisation restraints.')

if(restart) then !try to load theta_corr from restart file
        read(2, iostat=filestat) nwpolr_shell_restart
        if(filestat /= 0 .or. nwpolr_shell_restart /= nwpolr_shell) then
                write(*,102) 
                wshell(:)%theta_corr = 0.
        else
                backspace(2)
                read(2) nwpolr_shell_restart, wshell(:)%theta_corr
                write(*,103)
        end if
else
        wshell(:)%theta_corr = 0.
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
rshell = (0.5*(rout**3+ri**3))**(1./3.)


        ! --- Note below: 0.98750 = (1-1/epsilon) for water
wshell(is)%cstb = crgQtot*0.98750/(rho_wat*mu_w*4.*pi*rshell**2)
        write(*, 110) is, rout, ri
        rout = rout - dr
end do

n_max_insh = n_max_insh * 1.5 !take largest and add some extra
call allocate_watpol_arrays

end subroutine wat_shells

!-----------------------------------------------------------------------

subroutine watpol
! local variables
integer						:: iw,is,i,i3,il,jl,jw,imin,jmin
real(8)						:: dr,rw,rshell,rm,rc,scp
real(8)						:: tmin,arg,avtdum,dv,f0
real(8), save					:: f1(9),f2(3)
real(8), save					:: rmu(3),rcu(3)

! global variables used:
!  E, wshell, bndw0, deg2rad, angw0, nwat, theta, theta0, nat_pro, x, xwcent,
!  tdum, nwpolr_shell, list_sh, pi, nsort, istep, itdis_update, fkwpol, d

! reset wshell%n_insh
wshell(:)%n_insh = 0

! calculate theta(:), tdum(:), wshell%n_insh
do iw = 1, nwat
theta(iw)  = 0.0
theta0(iw) = 0.0

i  = nat_solute + iw*3-2
if(excl(i)) cycle ! skip excluded topology waters
i3 = i*3-3

rmu(1) = x(i3+4) + x(i3+7) - 2.*x(i3+1)   !Water vector
rmu(2) = x(i3+5) + x(i3+8) - 2.*x(i3+2)
rmu(3) = x(i3+6) + x(i3+9) - 2.*x(i3+3)
rm = sqrt ( rmu(1)**2 + rmu(2)**2 + rmu(3)**2 )
rmu(1) = rmu(1)/rm
rmu(2) = rmu(2)/rm
rmu(3) = rmu(3)/rm

rcu(1) = x(i3+1) - xwcent(1)    !Radial vector to OW
rcu(2) = x(i3+2) - xwcent(2) 
rcu(3) = x(i3+3) - xwcent(3) 
rc = sqrt ( rcu(1)**2 + rcu(2)**2 + rcu(3)**2 )
rcu(1) = rcu(1)/rc
rcu(2) = rcu(2)/rc
rcu(3) = rcu(3)/rc

scp = rmu(1)*rcu(1)+rmu(2)*rcu(2)+rmu(3)*rcu(3)   !Calculate angle between water vector and radial vector
if ( scp .gt.  1.0 ) scp =  1.0
if ( scp .lt. -1.0 ) scp = -1.0
theta(iw) = acos( scp )
tdum(iw) = theta(iw)

if ( rc > wshell(nwpolr_shell)%rout-wshell(nwpolr_shell)%dr ) then
  do is = nwpolr_shell, 2, -1
        if(rc <= wshell(is)%rout) exit
  end do
  wshell(is)%n_insh = wshell(is)%n_insh + 1
  list_sh(wshell(is)%n_insh,is) = iw
end if
end do

! sort the waters according to theta
do is = 1, nwpolr_shell
imin = 0
do il = 1, wshell(is)%n_insh
tmin = 2.*pi
do jl = 1, wshell(is)%n_insh
jw = list_sh(jl,is)
if ( tdum(jw) .lt. tmin ) then
  jmin = jw
  tmin = theta(jw)
end if
end do
imin = imin+1
nsort(imin,is) = jmin
  tdum(jmin) = 99999.

end do

end do

! calculate energy and force
if ( istep .ne. 0 .and. mod(istep,itdis_update) .eq. 0) then
call centered_heading('Water polarisation restraint data', '-')
write(*,'(a)') 'shell    <n>    <theta>    theta_0 theta_corr'
do is = 1, nwpolr_shell
wshell(is)%avtheta = wshell(is)%avtheta / real (itdis_update)
wshell(is)%avn_insh = wshell(is)%avn_insh / real (itdis_update)
wshell(is)%theta_corr = wshell(is)%theta_corr + wshell(is)%avtheta-acos(wshell(is)%cstb)
write (*,10) is,wshell(is)%avn_insh,wshell(is)%avtheta/deg2rad, &
     acos(wshell(is)%cstb)/deg2rad,wshell(is)%theta_corr/deg2rad
10	  format(i5,1x,f6.1,3x,f8.3,3x,f8.3,3x,f8.3)
wshell(is)%avtheta = 0.0
wshell(is)%avn_insh = 0.0
end do
end if

do is = 1, nwpolr_shell
if(wshell(is)%n_insh == 0) cycle !skip empty shell
avtdum = 0.0
do il = 1, wshell(is)%n_insh
iw = nsort(il,is)
arg = 1. + (1. - 2.*real(il))/real(wshell(is)%n_insh)
theta0(il) = acos ( arg )
theta0(il) = theta0(il)-3.*sin(theta0(il))*wshell(is)%cstb/2.
if ( theta0(il) .lt. 0.0 ) theta0(il) = 0.0
if ( theta0(il) .gt. pi)   theta0(il) = pi

avtdum = avtdum + theta(iw)

E%restraint%water_pol = E%restraint%water_pol + 0.5*fkwpol* &
     (theta(iw)-theta0(il)+wshell(is)%theta_corr)**2

dv = fkwpol*(theta(iw)-theta0(il)+wshell(is)%theta_corr)

i  = nat_solute + iw*3-2
i3 = i*3-3

rmu(1) = x(i3+4) + x(i3+7) - 2.*x(i3+1)
rmu(2) = x(i3+5) + x(i3+8) - 2.*x(i3+2)
rmu(3) = x(i3+6) + x(i3+9) - 2.*x(i3+3)
rm = sqrt ( rmu(1)**2 + rmu(2)**2 + rmu(3)**2 )
rmu(1) = rmu(1)/rm
rmu(2) = rmu(2)/rm
rmu(3) = rmu(3)/rm

rcu(1) = x(i3+1) - xwcent(1)
rcu(2) = x(i3+2) - xwcent(2)
rcu(3) = x(i3+3) - xwcent(3)
rc = sqrt ( rcu(1)**2 + rcu(2)**2 + rcu(3)**2 )
rcu(1) = rcu(1)/rc
rcu(2) = rcu(2)/rc
rcu(3) = rcu(3)/rc



scp = rmu(1)*rcu(1)+rmu(2)*rcu(2)+rmu(3)*rcu(3)
if ( scp .gt.  1.0 ) scp =  1.0
if ( scp .lt. -1.0 ) scp = -1.0
f0 = sin ( acos(scp) )
if ( abs(f0) .lt. 1.e-12 ) f0 = 1.e-12
f0 = -1.0 / f0
f0 = dv*f0

f1(1) = -2.*(rcu(1)-rmu(1)*scp)/rm
f1(2) = -2.*(rcu(2)-rmu(2)*scp)/rm
f1(3) = -2.*(rcu(3)-rmu(3)*scp)/rm
f1(4) =     (rcu(1)-rmu(1)*scp)/rm
f1(5) =     (rcu(2)-rmu(2)*scp)/rm
f1(6) =     (rcu(3)-rmu(3)*scp)/rm
f1(7) =     (rcu(1)-rmu(1)*scp)/rm
f1(8) =     (rcu(2)-rmu(2)*scp)/rm
f1(9) =     (rcu(3)-rmu(3)*scp)/rm

f2(1) = ( rmu(1)-rcu(1)*scp)/rc
f2(2) = ( rmu(2)-rcu(2)*scp)/rc
f2(3) = ( rmu(3)-rcu(3)*scp)/rc

d(i3+1) = d(i3+1) + f0 * ( f1(1) + f2(1) )
d(i3+2) = d(i3+2) + f0 * ( f1(2) + f2(2) )
d(i3+3) = d(i3+3) + f0 * ( f1(3) + f2(3) )
d(i3+4) = d(i3+4) + f0 * ( f1(4) )
d(i3+5) = d(i3+5) + f0 * ( f1(5) )
d(i3+6) = d(i3+6) + f0 * ( f1(6) )
d(i3+7) = d(i3+7) + f0 * ( f1(7) )
d(i3+8) = d(i3+8) + f0 * ( f1(8) )
d(i3+9) = d(i3+9) + f0 * ( f1(9) )
end do

wshell(is)%avtheta = wshell(is)%avtheta + avtdum/real(wshell(is)%n_insh)
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
6 format(A,T17, 6F10.2)

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
  write (*,32) 'Q-Q', istate, EQ(istate)%lambda, EQ(istate)%qq%el, EQ(istate)%qq%vdw
end do
write(*,*)
if(nat_solute > nqat) then !only if there is something else than Q-atoms in topology
  do istate =1, nstates
        write (*,32) 'Q-prot', istate,EQ(istate)%lambda, EQ(istate)%qp%el, EQ(istate)%qp%vdw
  end do
  write(*,*)
end if

if(nwat > 0) then
  do istate =1, nstates
        write (*,32) 'Q-wat', istate, EQ(istate)%lambda, EQ(istate)%qw%el, EQ(istate)%qw%vdw
  end do
  write(*,*)
end if

do istate =1, nstates
  write (*,32) 'Q-surr.',istate, EQ(istate)%lambda, &
                EQ(istate)%qp%el + EQ(istate)%qw%el, EQ(istate)%qp%vdw &
                + EQ(istate)%qw%vdw
end do
write(*,*)

do istate = 1, nstates
  write (*,36) 'Q-any', istate, EQ(istate)%lambda, EQ(istate)%qx%el,&
                EQ(istate)%qx%vdw, EQ(istate)%q%bond, EQ(istate)%q%angle,&
                EQ(istate)%q%torsion, EQ(istate)%q%improper
end do
write(*,*)

write(*,22) 'total', 'restraint'
do istate = 1, nstates
  write (*,32) 'Q-SUM', istate, EQ(istate)%lambda,&
                EQ(istate)%total, EQ(istate)%restraint
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
        write(*,47) 'Attempts', volume_try, volume_acc, real(volume_acc)/volume_try
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

nat3 = natom*3

rewind (3)
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

real(8)				::	boxc(1:3)
integer				::	i, j, starten, slutet
!the borders of the periodic box
real(8)				::	x_max, x_min, y_max, y_min, z_max, z_min
real(8)				:: cm(1:3)
integer				::  mvd_mol(1:nmol) !moved molecule 1=yes, 0=no
integer				::  k, ig
integer				::  pbib_start, pbib_stop

if( .not. rigid_box_centre ) then !if the box is allowed to float around, center on solute if present, otherwise center on solvent
        if( nat_solute > 0) then
                slutet = nat_solute
                !starten = ncgp_solute + 1
                
        else !if no solute present, centre box around solvent
                slutet = natom
                !starten = 1
        end if
        !find center
        boxc(:) = 0.0
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

if (nodeid .eq. 0) then
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
	pbib_start = nmol - (natom - nat_solute)/(istart_mol(nmol) - istart_mol(nmol-1)) + 1    !(number of mol - number of solvent molecules + 1)
end if
if ( .not. put_solvent_back_in_box ) then !we're not putting solvent back in box
        pbib_stop  = nmol - (natom - nat_solute)/(istart_mol(nmol) - istart_mol(nmol-1))        !(number of mol - number of solvent molecules)
end if          


do i=pbib_start,pbib_stop
        cm(:) =0.0
        do j = istart_mol(i),istart_mol(i+1)-1 !loop over all atoms in molecule i
                cm(1) = cm(1) + x(j*3-2)*mass(j)
                cm(2) = cm(2) + x(j*3-1)*mass(j)
                cm(3) = cm(3) + x(j*3  )*mass(j)
        end do
        cm(:) = cm(:) * mol_mass(i) !centre of mass of molecule i
		
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

!LRF: if molecule moved update all cgp_centers from first charge group to the last one
if (use_LRF) then

!Broadcast mvd_mol(:) & x(:)
#if defined(USE_MPI)
call MPI_Bcast(mvd_mol, nmol, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast mvd_mol(k)')
call MPI_Bcast(x, nat3, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast x')
#endif

do k=pbib_start,pbib_stop
		if (mvd_mol(k) == 1) then  
						do ig=iwhich_cgp(istart_mol(k)),iwhich_cgp(istart_mol(k+1)-1)
								lrf(ig)%cgp_cent(:) = 0
								do i  = cgp(ig)%first, cgp(ig)%last
										lrf(ig)%cgp_cent(:) = lrf(ig)%cgp_cent(:) + x(cgpatom(i)*3-2:cgpatom(i)*3)
								end do

						lrf(ig)%cgp_cent(:) = lrf(ig)%cgp_cent(:)/real(cgp(ig)%last - cgp(ig)%first +1)
						end do
		end if 
end do
end if

end subroutine put_back_in_box
!----------------------------------------------------------------------------
subroutine MC_volume()

real(8)									:: old_x(1:3*nat_pro), old_xx(1:3*nat_pro), x_move(1:3*nat_pro)
type(ENERGIES)							:: old_E, previous_E
type(Q_ENERGIES), dimension(1:nstates)	:: old_EQ
real(8)									:: old_boxl(1:3), old_inv(1:3)
real(8)									:: old_V, new_V, deltaLength
real(8)									:: deltaV, deltaE, deltaW
real(8)									:: box_min
integer									:: starten, slutet, i, j, sw_no !indeces
real(8)									:: randomno !random number
real(8)									:: new_x, new_y, new_z
real(8)									:: move_x, move_y, move_z
logical									:: acc
integer									:: longest , niter
real(8)									:: cubr_vol_ratio
real(8)									:: cm(1:3)
real(8)									::	old_EMorseD(max_qat)
real(8)									::	old_dMorse_i(3,max_qat)
real(8)									::	old_dMorse_j(3,max_qat)
type(NB_TYPE), allocatable						::old_nbpp(:)
type(NB_TYPE), allocatable						::old_nbpw(:)
integer(AI), allocatable						::old_nbww(:)
type(NBQP_TYPE), allocatable						::old_nbqp(:)
type(CGP_PAIR_TYPE), allocatable					:: old_nbpp_cgp(:)
type(CGP_PAIR_TYPE), allocatable					:: old_nbpw_cgp(:)
type(CGP_PAIR_TYPE), allocatable					:: old_nbww_cgp(:)
type(CGP_PAIR_TYPE), allocatable					:: old_nbqp_cgp(:)
integer									:: old_nbww_pair, old_nbww_true_pair , old_nbpw_pair , old_nbpw_cgp_pair 
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
old_EMorseD = EMorseD
old_dMorse_i = dMorse_i
old_dMorse_j = dMorse_j

!shake stuff
old_xx(:) = xx(:)



if (use_LRF) then
	old_nbww_pair = nbww_pair
	old_nbpw_pair = nbpw_pair
	old_nbpp_pair = nbpp_pair
	old_nbqp_pair = nbqp_pair
	old_nbww_true_pair = nbww_true_pair
	old_nbpw_cgp_pair = nbpw_cgp_pair
	old_nbpp_cgp_pair = nbpp_cgp_pair
	old_nbqp_cgp_pair = nbqp_cgp_pair
	allocate(old_nbww(calculation_assignment%ww%max))
	allocate(old_nbpw(calculation_assignment%pw%max))
	allocate(old_nbpp(calculation_assignment%pp%max))
	allocate(old_nbqp(calculation_assignment%qp%max))
	allocate(old_nbpw_cgp(size(nbpw_cgp, 1)))
	allocate(old_nbpp_cgp(size(nbpp_cgp, 1)))
	allocate(old_nbqp_cgp(size(nbqp_cgp, 1)))
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
	old_EQ = EQ(1:nstates)
	old_V = old_boxl(1) * old_boxl(2) * old_boxl(3)


	!new volume randomized
	randomno = randm(pressure_seed) ! 0<=randomno<=1
	randomno = randomno*2 - 1    !-1 <= randomno <= 1
	deltaV = randomno * max_vol_displ
	new_V = deltaV + old_V
	cubr_vol_ratio = (new_V/old_V)**(1./3.)
	write(*,4) 'old', 'new', 'delta'
	write(*,6) 'Volume', old_V, new_V, deltaV
	write(*,*)

	!compute new boxlenth and inv_boxl
	boxlength(1) = boxlength(1)*cubr_vol_ratio
	boxlength(2) = boxlength(2)*cubr_vol_ratio
	boxlength(3) = boxlength(3)*cubr_vol_ratio
	inv_boxl(:) = 1.0/boxlength(:)
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
		do i=1,nmol-1 !looping over all molecules but the last one
			cm(:) =0.0
			do j = istart_mol(i),istart_mol(i+1)-1 !loop over all atoms in molecule i
				cm(1) = cm(1) + x(j*3-2)*mass(j)
				cm(2) = cm(2) + x(j*3-1)*mass(j)
				cm(3) = cm(3) + x(j*3  )*mass(j)
			end do
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

		! the last molecule
		cm(:) = 0.0
		do j = istart_mol(nmol),natom
			cm(1) = cm(1) + x(j*3-2)*mass(j)
			cm(2) = cm(2) + x(j*3-1)*mass(j)
			cm(3) = cm(3) + x(j*3  )*mass(j)
		end do

		cm(:) = cm(:) * mol_mass(nmol)

		move_x = ( ( cm(1)-boxcentre(1) )*boxlength(1)/old_boxl(1) + boxcentre(1) ) - cm(1)
		move_y = ( ( cm(2)-boxcentre(2) )*boxlength(2)/old_boxl(2) + boxcentre(2) ) - cm(2)
		move_z = ( ( cm(3)-boxcentre(3) )*boxlength(3)/old_boxl(3) + boxcentre(3) ) - cm(3)

		do j=istart_mol(nmol) , natom
			x(j*3-2) = x(j*3-2) + move_x
			x(j*3-1) = x(j*3-1) + move_y
			x(j*3  ) = x(j*3  ) + move_z		
		end do

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
end if !use_LRF



!compute the new potential, in parallel if possible
call new_potential( previous_E ) !do we need a broadcast before this if using LRF (i.e. if the nonbond lists have been updated)

if (nodeid .eq. 0) then
	!Jamfor nya med gamla
	deltaE = E%potential - old_E%potential
	deltaW = deltaE + pressure * deltaV - nmol*Boltz*Temp0*log(new_V/old_V)
	write(*,4) 'old', 'new', 'delta'
	write(*,6) 'Potential', old_E%potential, E%potential, deltaE
	write(*,*)

	!accept or reject
	if( deltaW<=0.0 ) then
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
        EMorseD = old_EMorseD
	dMorse_i = old_dMorse_i
	dMorse_j = old_dMorse_j

	xx(:) = old_xx(:)

	if (use_LRF) then
		nbww_pair = old_nbww_pair
		nbpw_pair = old_nbpw_pair
		nbpp_pair = old_nbpp_pair
		nbqp_pair = old_nbqp_pair
		nbww_true_pair = old_nbww_true_pair
		nbpw_cgp_pair = old_nbpw_cgp_pair
		nbpp_cgp_pair = old_nbpp_cgp_pair
		nbqp_cgp_pair = old_nbqp_cgp_pair
		deallocate(nbww)
		deallocate(nbpw)
		deallocate(nbpp)
		deallocate(nbqp)
		deallocate(nbpw_cgp)
		deallocate(nbpp_cgp)
		deallocate(nbqp_cgp)
		calculation_assignment%ww%max = old_ww_max
		calculation_assignment%pw%max = old_pw_max
		calculation_assignment%pp%max = old_pp_max
		calculation_assignment%qp%max = old_qp_max
		allocate(nbww(calculation_assignment%ww%max))
		allocate(nbpw(calculation_assignment%pw%max))
		allocate(nbpp(calculation_assignment%pp%max))
		allocate(nbqp(calculation_assignment%qp%max))
		allocate(nbpw_cgp(size(old_nbpw_cgp, 1)))
		allocate(nbpp_cgp(size(old_nbpp_cgp, 1)))
		allocate(nbqp_cgp(size(old_nbqp_cgp, 1)))
		nbww(1:old_ww_max) = old_nbww(1:old_ww_max)
		nbpw(1:old_pw_max) = old_nbpw(1:old_pw_max)
		nbpp(1:old_pp_max) = old_nbpp(1:old_pp_max)
		nbqp(1:old_qp_max) = old_nbqp(1:old_qp_max)
		nbpw_cgp(:) = old_nbpw_cgp(:)
		nbpp_cgp(:) = old_nbpp_cgp(:)
		nbqp_cgp(:) = old_nbqp_cgp(:)
		lrf(:) = old_lrf(:)
		deallocate(old_nbww)
		deallocate(old_nbpw)
		deallocate(old_nbpp)
		deallocate(old_nbqp)
		deallocate(old_nbpw_cgp)
		deallocate(old_nbpp_cgp)
		deallocate(old_nbqp_cgp)
	end if
end if

	write(*,11) boxlength(1)*boxlength(2)*boxlength(3)
	write(*,12) boxlength
	write(*,13) sum(mass(:))/(boxlength(1)*boxlength(2)*boxlength(3)*1E-24*6.02E23)
	11 format('Final volume: ', f10.3)
	12 format('Final boxlength: ', 3f10.3)
	13 format('Final density (g/cm^3): ', f10.3)
	write(*,*)
	write(*,4) 'total', 'accepted', 'ratio'
	write(*,7) 'Attempts', volume_try, volume_acc, real(volume_acc)/volume_try
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
integer							:: istate,i, numrequest

!zero all energies
E%potential = 0.0
E%pp%el  = 0.0
E%pp%vdw = 0.0
E%pw%el  = 0.0
E%pw%vdw = 0.0
E%ww%el  = 0.0
E%ww%vdw = 0.0
E%qx%el    = 0.0
E%qx%vdw   = 0.0
E%restraint%protein = 0.0
E%LRF      = 0.0 
E%p%bond =  0.0
E%w%bond =  0.0
E%p%angle =  0.0
E%w%angle =  0.0
E%p%angle =  0.0
E%w%angle =  0.0
E%p%torsion =  0.0
E%w%torsion =  0.0
E%p%improper =  0.0
E%w%improper =  0.0

do istate = 1, nstates
	EQ(istate)%qq%el = 0.0
	EQ(istate)%qq%vdw = 0.0
	EQ(istate)%qp%el = 0.0
	EQ(istate)%qp%vdw = 0.0
	EQ(istate)%qw%el = 0.0
	EQ(istate)%qw%vdw = 0.0
	EQ(istate)%restraint = 0.0
end do

!reset derivatives ---
d(:) = 0.0

#if defined (USE_MPI)
!First post recieves for gathering data from slaves
if (nodeid .eq. 0) then
	call gather_nonbond
end if
#endif

!compute the new potential
select case(ivdw_rule)
        case(VDW_GEOMETRIC)
                call nonbond_pp_box
                call nonbond_pw_box
                if(qvdw_flag) then
                        call nonbond_qp_qvdw_box
                else
                        call nonbond_qp_box
                end if
                if(natom > nat_solute) then !if any solvent
                if(solvent_type == SOLVENT_SPC) then
                        !use the optimised SPC routine when possible
                        call nonbond_ww_spc_box
                        call nonbond_qw_spc_box
                elseif(solvent_type == SOLVENT_3ATOM) then !otherwise calc. LJ with all atoms
                        call nonbond_3atomsolvent_box
                call nonbond_qw_3atom_box
                end if
                end if
        case(VDW_ARITHMETIC)
                call nonbon2_pp_box
                call nonbon2_qp_box
                if(natom > nat_solute) then !if any solvent
                        call nonbon2_pw_box
                        call nonbon2_qw_box !no SPC-specific optimised routines here
                        call nonbon2_ww_box
                end if
end select


 if (use_LRF) then
     call lrf_taylor
 end if

if (nodeid .eq. 0) then
 call pot_energy_bonds   !added by Almlof. Will not lead to a significant performance decrease, but will reduce bugs and perhaps increase compiler optimizations? (esepecially once atom-based-scaling is introduced)
 call p_restrain

 if(.not. qq_use_library_charges) then
        if(ivdw_rule .eq. 1 ) then
                call nonbond_qq
        elseif ( ivdw_rule .eq. 2 ) then
                call nonbon2_qq
        end if
 else
        if ( ivdw_rule .eq. 1 ) then
                call nonbond_qq_lib_charges
        else if ( ivdw_rule .eq. 2 ) then
                call nonbon2_qq_lib_charges
        end if
 end if
end if

#if defined(USE_MPI)
if (nodeid .ne. 0) then  !Slave nodes
	call gather_nonbond
end if
#endif

if (nodeid .eq. 0) then 
#if (USE_MPI)
if (use_excluded_groups) then
        numrequest = 4
else
        numrequest = 3
end if
do i = 1, numrequest
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
	  EQ(1:nstates)%qp%el  = EQ(1:nstates)%qp%el  + EQ_recv(1:nstates,i)%qp%el
	  EQ(1:nstates)%qp%vdw = EQ(1:nstates)%qp%vdw + EQ_recv(1:nstates,i)%qp%vdw
	  EQ(1:nstates)%qw%el  = EQ(1:nstates)%qw%el  + EQ_recv(1:nstates,i)%qw%el
	  EQ(1:nstates)%qw%vdw = EQ(1:nstates)%qw%vdw + EQ_recv(1:nstates,i)%qw%vdw
	end do
#endif

	!summation of energies
	do istate = 1, nstates
			! update EQ
			EQ(istate)%qx%el  = EQ(istate)%qq%el +EQ(istate)%qp%el +EQ(istate)%qw%el
			EQ(istate)%qx%vdw = EQ(istate)%qq%vdw+EQ(istate)%qp%vdw+EQ(istate)%qw%vdw
			E%qx%el    = E%qx%el    + EQ(istate)%qx%el   *EQ(istate)%lambda
			E%qx%vdw   = E%qx%vdw   + EQ(istate)%qx%vdw  *EQ(istate)%lambda

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
! and now also the EQ_gc variables from the group contributions
! Allocate  - status


subroutine gather_nonbond()

integer,parameter                       :: vars=3
integer,dimension(4,numnodes-1)         :: tag
integer,dimension(vars)	                :: blockcnt,ftype 
integer(kind=MPI_ADDRESS_KIND), dimension(vars)	:: fdisp, base
integer                                 :: mpitype_package,mpitype_send
integer                                 :: i,j,istate

do i=1,numnodes-1
tag(1,i)=numnodes*100+i
tag(2,i)=numnodes*200+i
tag(3,i)=numnodes*300+i
tag(4,i)=numnodes*400+i
end do

if (nodeid .eq. 0) then        !master


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
  call MPI_IRecv(EQ_recv(1,i), nstates*2*2, MPI_REAL8, i, tag(3,i), MPI_COMM_WORLD, &
       request_recv(i,3),ierr)
  if (use_excluded_groups) then
	call MPI_IRecv(EQ_gc_recv(1,1,i), nstates*2*2*ngroups_gc, MPI_REAL8, i, tag(4,i), MPI_COMM_WORLD, &
        request_recv(i,4),ierr)
  end if
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
EQ_send(1:nstates)%qp%el  = EQ(1:nstates)%qp%el
EQ_send(1:nstates)%qp%vdw = EQ(1:nstates)%qp%vdw
EQ_send(1:nstates)%qw%el  = EQ(1:nstates)%qw%el
EQ_send(1:nstates)%qw%vdw = EQ(1:nstates)%qw%vdw
if (use_excluded_groups) then
        EQ_gc_send(1:ngroups_gc,1:nstates)%qp%el  = EQ_gc(1:ngroups_gc,1:nstates)%qp%el
        EQ_gc_send(1:ngroups_gc,1:nstates)%qp%vdw  = EQ_gc(1:ngroups_gc,1:nstates)%qp%vdw
end if 

! See comments above on the IRecv part
call MPI_Send(d, natom*3, MPI_REAL8, 0, tag(1,nodeid), MPI_COMM_WORLD,ierr) 
if (ierr .ne. 0) call die('gather_nonbond/Send d')
call MPI_Send(E_send, 3*2+1, MPI_REAL8, 0, tag(2,nodeid), MPI_COMM_WORLD,ierr) 
if (ierr .ne. 0) call die('gather_nonbond/Send E_send')
call MPI_Send(EQ_send, nstates*2*2, MPI_REAL8, 0, tag(3,nodeid), MPI_COMM_WORLD,ierr) 
if (use_excluded_groups) then
	call MPI_Send(EQ_gc_send(1,1),nstates*2*2*ngroups_gc, MPI_REAL8, 0, tag(4,nodeid), MPI_COMM_WORLD,ierr)
end if
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




