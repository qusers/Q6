! (C) 2014 Uppsala Molekylmekaniska HB, Uppsala, Sweden
! potene.f90 
! based on md.f90
! by Johan Åqvist, John Marelius, Anders Kaplan, Isabella Feierberg, Martin Nervall & Martin Almlöf
! wrapper for all potential energy calculation
! by Paul Bauer

module POTENE

! used modules
use NONBONDENE
use BONDENE
use SIZES
use GLOBALS
!$ use omp_lib
implicit none



contains

subroutine pot_energy(E_loc,EQ_loc)
! arguments
TYPE(ENERGIES)                          :: E_loc
TYPE(OQ_ENERGIES)                       :: EQ_loc(:)
! local variables


integer					:: istate, i, nat3
integer					:: is, j
#if defined (PROFILING)
real(kind=prec)					:: start_loop_time1
real(kind=prec)					:: start_loop_time2
real(kind=prec)					:: start_loop_time3
#endif

! --- reset all energies

E_loc%potential = zero
!E%kinetic = zero	! no need to reset because it will be assigned its final value at once
E_loc%LRF = zero
E_loc%p%bond  = zero
E_loc%p%angle = zero
E_loc%p%torsion = zero
E_loc%p%improper = zero
E_loc%w%bond  = zero
E_loc%w%angle = zero
E_loc%w%torsion = zero
E_loc%w%improper = zero
E_loc%q%bond  = zero
E_loc%q%angle = zero
E_loc%q%torsion   = zero
E_loc%q%improper   = zero
E_loc%pp%el  = zero
E_loc%pp%vdw = zero
E_loc%pw%el  = zero
E_loc%pw%vdw = zero
E_loc%ww%el  = zero
E_loc%ww%vdw = zero
E_loc%qx%el    = zero
E_loc%qx%vdw   = zero
!E%restraint%total = zero	! will be assigned its final value later
E_loc%restraint%fix = zero
E_loc%restraint%shell = zero
E_loc%restraint%protein = zero
E_loc%restraint%solvent_radial = zero
E_loc%restraint%water_pol = zero
do istate = 1, nstates
!EQ(istate)%lambda set by initialize
!EQ(istate)%total assigned its final value later
EQ_loc(istate)%q%bond = zero
EQ_loc(istate)%q%angle = zero
EQ_loc(istate)%q%torsion = zero
EQ_loc(istate)%q%improper = zero
!EQ(istate)%qx%el = zero	! assigned its final value later
!EQ(istate)%qx%vdw = zero	! assigned its final value later
EQ_loc(istate)%qq%el = zero
EQ_loc(istate)%qq%vdw = zero
EQ_loc(istate)%qp%el = zero
EQ_loc(istate)%qp%vdw = zero
EQ_loc(istate)%qw%el = zero
EQ_loc(istate)%qw%vdw = zero
EQ_loc(istate)%restraint = zero
end do

!reset derivatives ---
d = d * zero

! --- calculate the potential energy and derivatives ---
! *** nonbonds distribueras

#if defined (USE_MPI)
if (nodeid .eq. 0) then
!First post recieves for gathering data from slaves
call gather_nonbond(E_loc,EQ_loc(:))
end if
#endif

#if defined (PROFILING)
start_loop_time2 = rtime()
#endif

! classical nonbonds
call pot_energy_nonbonds(E_loc,EQ_loc(:))
#if defined (PROFILING)
profile(10)%time = profile(10)%time + rtime() - start_loop_time2
#endif

if (nodeid .eq. 0) then


! classical bond interactions (master node only)
#if defined (PROFILING)
start_loop_time1 = rtime()
#endif
call pot_energy_bonds(E_loc)
#if defined (PROFILING)
profile(8)%time = profile(8)%time + rtime() - start_loop_time1
#endif

! various restraints
if( .not. use_PBC ) then
   call fix_shell     !Restrain all excluded atoms plus heavy solute atoms in the inner shell.
end if

call p_restrain(E_loc%restraint%protein,EQ_loc(:)%restraint,EQ(:)%lambda)
!Seq. restraints, dist. restaints, etc

if( .not. use_PBC ) then
        if(nwat > 0) then
          call restrain_solvent(E_loc%restraint)
        if (wpol_restr) call watpol(E_loc%restraint)
        end if
end if

! q-q nonbonded interactions
call nonbond_qq(EQ_loc(:)%qq,EQ_loc(:)%lambda)
call nonbond_qqp(EQ_loc(:)%qp,EQ_loc(:)%lambda)

! q-atom bonded interactions: loop over q-atom states
do istate = 1, nstates
  ! bonds, angles, torsions and impropers
  call qbond (EQ_loc(istate)%q%bond,EQ(istate)%lambda,istate)
  call qangle (EQ_loc(istate)%q%angle,EQ(istate)%lambda,istate)
  if(ff_type == FF_CHARMM) call qurey_bradley(EQ_loc(istate)%q%angle,EQ_loc(istate)%lambda,istate)
  call qtorsion (EQ_loc(istate)%q%torsion,EQ_loc(istate)%lambda,istate)
  call qimproper (EQ_loc(istate)%q%improper,EQ_loc(istate)%lambda,istate)
end do
#if defined (PROFILING)
profile(9)%time = profile(9)%time + rtime() - start_loop_time1 - profile(8)%time
#endif
#if defined(USE_MPI)
else  !Slave nodes
call gather_nonbond(E_loc,EQ_loc(:))
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
  E_loc%pp%el   = E_loc%pp%el  + E_recv(i)%pp%el
  E_loc%pp%vdw  = E_loc%pp%vdw + E_recv(i)%pp%vdw
  E_loc%pw%el   = E_loc%pw%el  + E_recv(i)%pw%el
  E_loc%pw%vdw  = E_loc%pw%vdw + E_recv(i)%pw%vdw
  E_loc%ww%el   = E_loc%ww%el  + E_recv(i)%ww%el
  E_loc%ww%vdw  = E_loc%ww%vdw + E_recv(i)%ww%vdw
  E_loc%lrf     = E_loc%lrf    + E_recv(i)%lrf
	do istate=1,nstates
  EQ_loc(istate)%qp%el  = EQ_loc(istate)%qp%el  + EQ_recv(istate,i)%qp%el
  EQ_loc(istate)%qp%vdw = EQ_loc(istate)%qp%vdw + EQ_recv(istate,i)%qp%vdw
  EQ_loc(istate)%qw%el  = EQ_loc(istate)%qw%el  + EQ_recv(istate,i)%qw%el
  EQ_loc(istate)%qw%vdw = EQ_loc(istate)%qw%vdw + EQ_recv(istate,i)%qw%vdw
	end do
end do
#endif
! q-atom energy summary
do istate = 1, nstates
! update EQ
EQ_loc(istate)%qx%el  = EQ_loc(istate)%qq%el +EQ_loc(istate)%qp%el +EQ_loc(istate)%qw%el
EQ_loc(istate)%qx%vdw = EQ_loc(istate)%qq%vdw+EQ_loc(istate)%qp%vdw+EQ_loc(istate)%qw%vdw

EQ_loc(istate)%total =  EQ_loc(istate)%q%bond + EQ_loc(istate)%q%angle   &
  + EQ_loc(istate)%q%torsion  + EQ_loc(istate)%q%improper + EQ_loc(istate)%qx%el &
  + EQ_loc(istate)%qx%vdw  + EQ_loc(istate)%restraint
end do

! update E with an average of all states
E_loc%q%bond     = E_loc%q%bond     + EQ_loc(istate)%q%bond     * EQ_loc(istate)%lambda
E_loc%q%angle    = E_loc%q%angle    + EQ_loc(istate)%q%angle    * EQ_loc(istate)%lambda
E_loc%q%torsion  = E_loc%q%torsion  + EQ_loc(istate)%q%torsion  * EQ_loc(istate)%lambda
E_loc%q%improper = E_loc%q%improper + EQ_loc(istate)%q%improper * EQ_loc(istate)%lambda
E_loc%qx%el    = E_loc%qx%el    + EQ_loc(istate)%qx%el   *EQ_loc(istate)%lambda
E_loc%qx%vdw   = E_loc%qx%vdw   + EQ_loc(istate)%qx%vdw  *EQ_loc(istate)%lambda

! update E%restraint%protein with an average of all states
E_loc%restraint%protein = E_loc%restraint%protein + EQ_loc(istate)%restraint*EQ_loc(istate)%lambda

! total energy summary
E_loc%restraint%total = E_loc%restraint%fix + E_loc%restraint%shell + &
E_loc%restraint%protein + E_loc%restraint%solvent_radial + E_loc%restraint%water_pol

E_loc%potential = E_loc%p%bond + E_loc%w%bond + E_loc%p%angle + E_loc%w%angle + E_loc%p%torsion + &
E_loc%p%improper + E_loc%pp%el + E_loc%pp%vdw + E_loc%pw%el + E_loc%pw%vdw + E_loc%ww%el + &
E_loc%ww%vdw + E_loc%q%bond + E_loc%q%angle + E_loc%q%torsion + &
E_loc%q%improper + E_loc%qx%el + E_loc%qx%vdw + E_loc%restraint%total + E_loc%LRF
end if

end subroutine pot_energy

!-----------------------------------------------------------------------

subroutine pot_energy_bonds(E_loc)
! arguments
TYPE(ENERGIES)                  :: E_loc
! bond, angle, torsion and improper potential energy
select case(ff_type)
case(FF_GROMOS)
        E_loc%p%bond = bond(1, nbonds_solute,bnd,bondlib)
        E_loc%w%bond = bond(nbonds_solute+1, nbonds,bnd,bondlib)
        E_loc%p%angle = angle(1, nangles_solute,ang,anglib)
        E_loc%w%angle = angle(nangles_solute+1, nangles,ang,anglib)
        E_loc%p%torsion = torsion(1, ntors_solute,tor,torlib)
        E_loc%w%torsion = torsion(ntors_solute+1, ntors,tor,torlib)
        E_loc%p%improper = improper(1, nimps_solute,imp,implib)
        E_loc%w%improper = improper(nimps_solute+1, nimps,imp,implib)
case(FF_AMBER)
        E_loc%p%bond = bond(1, nbonds_solute,bnd,bondlib)
        E_loc%w%bond = bond(nbonds_solute+1, nbonds,bnd,bondlib)
        E_loc%p%angle = angle(1, nangles_solute,ang,anglib)
        E_loc%w%angle = angle(nangles_solute+1, nangles,ang,anglib)
        E_loc%p%torsion = torsion(1, ntors_solute,tor,torlib)
        E_loc%w%torsion = torsion(ntors_solute+1, ntors,tor,torlib)
        E_loc%p%improper = improper2(1, nimps_solute,imp,implib)
        E_loc%w%improper = improper2(nimps_solute+1, nimps,imp,implib)
case(FF_CHARMM)
        E_loc%p%bond = bond(1, nbonds_solute,bnd,bondlib)
        E_loc%w%bond = bond(nbonds_solute+1, nbonds,bnd,bondlib)
        E_loc%p%angle = angle(1, nangles_solute,ang,anglib)
        E_loc%w%angle = angle(nangles_solute+1, nangles,ang,anglib)
        E_loc%p%angle = E%p%angle + urey_bradley(1, nangles_solute,ang,anglib)
        E_loc%w%angle = E%w%angle + urey_bradley(nangles_solute+1, nangles,ang,anglib)
        E_loc%p%torsion = torsion(1, ntors_solute,tor,torlib)
        E_loc%w%torsion = torsion(ntors_solute+1, ntors,tor,torlib)
        E_loc%p%improper = improper(1, nimps_solute,imp,implib)
        E_loc%w%improper = improper(nimps_solute+1, nimps,imp,implib)
end select
end subroutine pot_energy_bonds

!-----------------------------------------------------------------------
subroutine pot_energy_nonbonds(E_loc,EQ_loc)
! arguments
TYPE(ENERGIES)                  :: E_loc
TYPE(OQ_ENERGIES)               :: EQ_loc(:)
!nonbonded interactions
!$omp parallel default(none) shared(use_PBC,ivdw_rule,natom,nat_solute,solvent_type,use_LRF,ntors,ntors_solute) &
!$omp reduction(+:d)
if( use_PBC ) then !periodic box
        if(natom > nat_solute) then
                if((ivdw_rule.eq.VDW_GEOMETRIC).and. &
                        (solvent_type == SOLVENT_SPC)) then
                        call nonbond_qw_spc_box(EQ_loc(:)%qw,EQ_loc(:)%lambda)
                        call nonbond_ww_spc_box(E_loc%ww)
                else
                        call nonbond_qw_box(EQ_loc(:)%qw,EQ_loc(:)%lambda)
                        call nonbond_ww_box(E_loc%ww)
                end if
                if(ntors>ntors_solute) then
                        call nonbond_solvent_internal_box(E_loc%ww)
                end if
                call nonbond_pw_box(E_loc%pw)
        end if
        call nonbond_pp_box(E_loc%pp)
        call nonbond_qp_box(EQ_loc(:)%qp,EQ_loc(:)%lambda)
else
        if(natom > nat_solute) then
                if((ivdw_rule.eq.VDW_GEOMETRIC).and. &
                        (solvent_type == SOLVENT_SPC)) then
                        call nonbond_qw_spc(EQ_loc(:)%qw,EQ_loc(:)%lambda)
                        call nonbond_ww_spc(E_loc%ww)
                else
                        call nonbond_qw(EQ_loc(:)%qw,EQ_loc(:)%lambda)
                        call nonbond_ww(E_loc%ww)
                end if
! now we get started on the funky stuff with solvent torsions
! only run this if we have torsions in the solvent
! because it means we have at least 1-4 interactions
                if(ntors>ntors_solute) then
                        call nonbond_solvent_internal(E_loc%ww)
                end if
                call nonbond_pw(E_loc%pw)
        end if
        call nonbond_pp(E_loc%pp)
        call nonbond_qp(EQ_loc(:)%qp,EQ_loc(:)%lambda)
end if

!LRF
!$omp single
if (use_LRF) then
        call lrf_taylor(E_loc%lrf)
end if
!$omp end single
!$omp end parallel
end subroutine pot_energy_nonbonds

!-----------------------------------------------------------------------


!Restrain all excluded atoms plus heavy solute atoms in the inner shell.
subroutine fix_shell
! local variables
integer						::	i,i3
real(kind=prec)					::	fk,r2,erst
TYPE(qr_dist)                                   :: dist
! global variables used:
!  E, nat_pro, excl, shell, heavy, fk_fix, fk_pshell, x, xtop, d

do i = 1, nat_pro
if (excl(i) .or. shell(i)) then
! decide which fk to use
if ( excl(i) ) then 
fk = fk_fix
if ( freeze ) then
v(i) = v(i) * zero
x(i) = xtop(i)
fk = 0
endif
else
fk = fk_pshell
end if

! calculate drift from topology
dist = q_dist(x(i),xtop(i))
r2   = dist%r2
erst = 0.5_prec*fk*r2

! update restraint energies
if ( excl(i) ) E%restraint%fix   = E%restraint%fix + erst
if ( shell(i) ) E%restraint%shell = E%restraint%shell + erst 

! update forces
d(i) = d(i) + dist%vec * fk
end if
end do
end subroutine fix_shell

!----------------------------------------------------------------------
end module POTENE




