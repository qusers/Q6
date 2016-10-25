! (C) 2014 Uppsala Molekylmekaniska HB, Uppsala, Sweden
! qalloc.f90
! based on
! md.f90
! by Johan Åqvist, John Marelius, Anders Kaplan, Isabella Feierberg, Martin Nervall & Martin Almlöf
! memory management

module QALLOC

! used modules
!use PROFILING
use SIZES
use MPIGLOB
use GLOBALS
use QATOM
use QMATH

contains


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

allocate(x(natom), &
                xx(natom), &
                v(natom), &
                d(natom), &
                winv(natom), &
                iqatom(natom), &
                stat=alloc_status)
call check_alloc('atom data arrays')
!set arrays to zero after allocation
!we can not trust the compiler to do this
        x(:)  = x(:)  * zero
        xx(:) = xx(:) * zero
        v(:)  = v(:)  * zero
        d(:)  = d(:)  * zero
	winv(:)= zero
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

end module QALLOC




