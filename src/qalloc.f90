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
use TRJ
implicit none
#ifdef USE_MPI
include "mpif.h"
#endif
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
        x(:)     = zero
        xx(:)    = zero
        v(:)     = zero
        d(:)     = zero
	winv(:)  = zero
	iqatom(:)= 0
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
         d_recv(natom,numnodes-1), &
         E_recv(numnodes-1), &
         EQ_recv(nstates,numnodes-1), &
         stat=alloc_status)
 call check_alloc('MPI data arrays head node')
else
 allocate(E_send(1), &
          EQ_send(nstates), &
          stat=alloc_status)
 call check_alloc('MPI data arrays slave nodes')
end if
#ifdef _OPENMP
allocate(qomp_elec(nstates),qomp_vdw(nstates),stat=alloc_status)
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
	integer			:: ii
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
if(allocated(const_mol)) then
	do ii=1,nmol
                if (associated(const_mol(ii)%bond))     nullify(const_mol(ii)%bond)
                if (associated(const_mol(ii)%linc))     nullify(const_mol(ii)%linc)
                if (allocated(const_bonds(ii)%bond))    deallocate(const_bonds(ii)%bond)
                if (allocated(const_lincs)) then
                        if (allocated(const_lincs(ii)%con))     deallocate(const_lincs(ii)%con)
                        if (allocated(const_lincs(ii)%length))  deallocate(const_lincs(ii)%length)
                        if (allocated(const_lincs(ii)%length2)) deallocate(const_lincs(ii)%length2)
                        if (allocated(const_lincs(ii)%S))       deallocate(const_lincs(ii)%S)
                        if (allocated(const_lincs(ii)%A))       deallocate(const_lincs(ii)%A)
                        if (allocated(const_lincs(ii)%B))       deallocate(const_lincs(ii)%B)
                        if (allocated(const_lincs(ii)%coef))    deallocate(const_lincs(ii)%coef)
                        if (allocated(const_lincs(ii)%rhs))     deallocate(const_lincs(ii)%rhs)
                        if (allocated(const_lincs(ii)%sol))     deallocate(const_lincs(ii)%sol)
                end if
	end do
	deallocate(const_mol)
        if(allocated(const_lincs)) deallocate(const_lincs)
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

! new temperature stuff
if(allocated(tscale)) deallocate(tscale)

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
write(*,'(a)') 'ABNORMAL TERMINATION of Q'
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
stop 'Q terminated abnormally'
#endif

end subroutine die

!-----------------------------------------------------------------------
! need writing out routines here because they are linked directly in die
! TODO
! needs to split apart and made much easier to handle

subroutine write_out
! arguments
! local variables
integer					::	i,istate,ngrms

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
        ! sum up all forces, can't do dot product as before
        grms = zero
        do ngrms = 1,natom
        grms = grms + q_dotprod(d(ngrms),d(ngrms)/real(3*natom,kind=prec))
        end do
        grms = q_sqrt(grms)
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
        write(*,45) boxlength%x*boxlength%y*boxlength%z
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
#ifdef HAVEQUAD
elseif (prec .eq. quadprecision) then
canary = -13337
#endif
else
stop 'No such precision'
end if
rewind (3)
!new canary on top of file
write (3) canary
write (3) nat3, (x(i),i=1,natom)
write (3) nat3, (v(i),i=1,natom)
!save dynamic polarisation restraint data
	if(wpol_restr .and. allocated(wshell)) then
        write (3) nwpolr_shell, wshell(:)%theta_corr
end if

if( use_PBC )then
        write(3) boxlength
        write(3) boxcentre
end if
end subroutine write_xfin

!----------------------------------------------------------------------------

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

!----------------------------------------------------------------------

end module QALLOC




