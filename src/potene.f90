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
!$ use omp_lib
implicit none



contains

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
d = d * zero

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
        E%p%bond = bond(1, nbonds_solute,bnd,bondlib)
        E%w%bond = bond(nbonds_solute+1, nbonds,bnd,bondlib)
        E%p%angle = angle(1, nangles_solute,ang,anglib)
        E%w%angle = angle(nangles_solute+1, nangles,ang,anglib)
        E%p%torsion = torsion(1, ntors_solute,tor,torlib)
        E%w%torsion = torsion(ntors_solute+1, ntors,tor,torlib)
        E%p%improper = improper(1, nimps_solute,imp,implib)
        E%w%improper = improper(nimps_solute+1, nimps,imp,implib)
case(FF_AMBER)
        E%p%bond = bond(1, nbonds_solute,bnd,bondlib)
        E%w%bond = bond(nbonds_solute+1, nbonds,bnd,bondlib)
        E%p%angle = angle(1, nangles_solute,ang,anglib)
        E%w%angle = angle(nangles_solute+1, nangles,ang,anglib)
        E%p%torsion = torsion(1, ntors_solute,tor,torlib)
        E%w%torsion = torsion(ntors_solute+1, ntors,tor,torlib)
        E%p%improper = improper2(1, nimps_solute,imp,implib)
        E%w%improper = improper2(nimps_solute+1, nimps,imp,implib)
case(FF_CHARMM)
        E%p%bond = bond(1, nbonds_solute,bnd,bondlib)
        E%w%bond = bond(nbonds_solute+1, nbonds,bnd,bondlib)
        E%p%angle = angle(1, nangles_solute,ang,anglib)
        E%w%angle = angle(nangles_solute+1, nangles,ang,anglib)
        E%p%angle = E%p%angle + urey_bradley(1, nangles_solute,ang,anglib)
        E%w%angle = E%w%angle + urey_bradley(nangles_solute+1, nangles,ang,anglib)
        E%p%torsion = torsion(1, ntors_solute,tor,torlib)
        E%w%torsion = torsion(ntors_solute+1, ntors,tor,torlib)
        E%p%improper = improper(1, nimps_solute,imp,implib)
        E%w%improper = improper(nimps_solute+1, nimps,imp,implib)
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

end module POTENE




