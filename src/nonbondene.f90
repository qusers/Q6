! (C) 2014 Uppsala Molekylmekaniska HB, Uppsala, Sweden
! nonbondene.f90 
! based on md.f90
! by Johan Åqvist, John Marelius, Anders Kaplan, Isabella Feierberg, Martin Nervall & Martin Almlöf
! calculation of nonbonded energies and list generation
! by Paul Bauer

module NONBONDENE

! used modules
use SIZES
use TRJ
use QATOM
use EXC
use QALLOC
use QSIMPREP
!$ use omp_lib
implicit none



contains

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

subroutine nbpp_count(npp, nppcgp, Rcut)
! arguments
integer						:: npp
integer						:: nppcgp(:)
real(kind=prec)					:: Rcut
! local variables
integer						:: i,j,ig,jg,ia,ja,nl
TYPE(qr_vec)					:: shift
real(kind=prec)						:: rcut2,r2
integer						:: LJ_code
! This routine counts non-bonded solute-solute atom pairs 
! excluding any Q-atoms.

! uses the global variables:
!  Rcpp, ncgp, cgp, excl, x, cgpatom, iqatom, ljcod, crg, 
!  iaclib, max_nbr_range, listex, nexlong, listexlong


npp = 0
rcut2 = Rcut*Rcut

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
  r2 = q_dist4(x(ia),x(ja))
  else
	shift = qvec_sub(x(ia),x(ja))
	r2 = q_dist4(shift,q_vecscale(boxlength,nint(q_vecscale(shift,inv_boxl))))
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


subroutine nbpplis2(Rcut)
! args
real(kind=prec)					:: Rcut
! local variables
integer						:: i,j,ig,jg,ia,ja,nl,inside,jgr
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
rcut2 = Rcut*Rcut

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
                        ja = cgp(jgr)%first
                        do while ((ja .le. cgp(jgr)%last) .and. (inside .eq. 0))
                                j = cgpatom(ja)
                                r2 = q_dist4(x(i),x(j))
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

subroutine nbpplis2_box(Rcut)
! args
real(kind=prec)					:: Rcut
! local variables
integer						:: i,j,ig,jg,ia,ja,nl,inside,jgr
real(kind=prec)					:: rcut2,r2
TYPE(qr_vec)					:: shift

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
  rcut2 = Rcut*Rcut
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
					shift = qvec_sub(x(i),x(j))
					r2 = q_dist4(shift,q_vecscale(boxlength,nint(q_vecscale(shift,inv_boxl))))
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
subroutine nbpplis2_box_lrf(Rcut,RLRF)
! args
real(kind=prec)						:: Rcut,RLRF
! local variables
integer						:: i,j,ig,jg,ia,ja,nl,inside,jgr
TYPE(qr_vec)					:: shift
real(kind=prec)						:: rcut2,r2

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
  rcut2 = Rcut*Rcut
  RcLRF2 = RLRF*RLRF
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
                                shift = qvec_sub(x(i),x(j))
                                r2 = q_dist4(shift,q_vecscale(boxlength,nint(q_vecscale(shift,inv_boxl))))
                                if ( r2 .le. rcut2 ) then
! one atom pair is within cutoff: set inside
                                        inside = 1
!$omp critical
                                        if (nbpp_cgp_pair .eq. size(nbpp_cgp, 1) )  call reallocate_nbpp_cgp
                                        nbpp_cgp_pair = nbpp_cgp_pair + 1
                                        nbpp_cgp(nbpp_cgp_pair)%i = i
                                        nbpp_cgp(nbpp_cgp_pair)%j = j
!$omp end critical
                                elseif ((r2 .le. RcLRF2) .or. (RLRF .eq. -one)) then
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
subroutine nbpplis2_lrf(Rcut,RLRF)
! args
real(kind=prec)					:: Rcut,RLRF
! local variables
integer						:: i,j,ig,jg,ia,ja,nl,is,jgr
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
rcut2 = Rcut*Rcut
RcLRF2 = RLRF*RLRF

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
			ja = cgp(jgr)%first
                        do while (( ja.le.cgp(jgr)%last) .and. (inside.eqv..false.))
                                j = cgpatom(ja)
				r2 = q_dist4(x(i),x(j))

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

subroutine nbpplist(Rcut)
! args
real(kind=prec)					:: Rcut
! local variables
integer						:: i,j,ig,jg,ia,ja,nl,jgr
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
		r2 = q_dist4(x(ia),x(ja))

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

subroutine nbpplist_box(Rcut)
! args
real(kind=prec)						:: Rcut
! local variables
integer						:: i,j,ig,jg,ia,ja,nl,ig_sw, jg_sw,jgr
TYPE(qr_vec)					:: shift
real(kind=prec)						:: rcut2,r2
integer						:: LJ_code
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
	  shift = qvec_sub(x(ig_sw),x(jg_sw))
	  r2 = q_dist4(shift,q_vecscale(boxlength,nint(q_vecscale(shift,inv_boxl))))

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

subroutine nbpplist_lrf(Rcut,RLRF)
! args
real(kind=prec)					:: Rcut,RLRF
! local variables
integer						:: i,j,ig,jg,ia,ja,nl,is,is3,jgr
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
		r2 = q_dist4(x(is),x(ja))

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
subroutine nbpplist_box_lrf(Rcut,RLRF)
! args
real(kind=prec)						:: Rcut,RLRF
! local variables
integer						:: i,j,ig,jg,ia,ja,nl,ig_sw, jg_sw, is3, jgr
real(kind=prec)						:: rcut2,r2
TYPE(qr_vec)					:: shift
integer						:: LJ_code

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
  rcut2 = Rcut*Rcut
  RcLRF2 = RLRF*RLRF
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
	 	shift = qvec_sub(x(ig_sw),x(jg_sw))
	 	r2 = q_dist4(shift,q_vecscale(boxlength,nint(q_vecscale(shift,inv_boxl))))
	
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
                elseif ((r2 .le. RcLRF2) .or. (RLRF .eq. -one)) then
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
subroutine nbpw_count(npw, npwcgp,Rcut)
! arguments
integer						:: npw
integer						:: npwcgp(:)
real(kind=prec)					:: Rcut
! local variables
integer						:: i,ig,jg,ia,ja
real(kind=prec)						:: rcut2,r2
integer						:: LJ_code
TYPE(qr_vec)					:: shift

!	This routine makes a list of non-bonded solute-solvent atom pairs
!	excluding Q-atoms.

! uses the global variables:
!  Rcpw, ncgp, cgp, excl, nwat, nat_solute, x, cgpatom, iqatom, ljcod, crg, iaclib

npw = 0
rcut2 = Rcut*Rcut

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
	r2 = q_dist3(x(ia),x(ja))
  else
	shift = qvec_sub(x(ia),x(ja))
	r2 = q_dist4(shift,q_vecscale(boxlength,nint(q_vecscale(shift,inv_boxl))))
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

subroutine nbpwlis2(Rcut)
! args
real(kind=prec)					:: Rcut
! local variables
integer						:: i,ig,jg,ia,ja, inside,jgr,j
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
			r2 = q_dist4(x(ia),x(ja))

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
subroutine nbpwlis2_box(Rcut)
! args
real(kind=prec)						:: Rcut
! local variables
integer						:: i,ig,jg,ia,ja,ig_atom, inside, jgr, LJ_code,j
TYPE(qr_vec)					:: shift
real(kind=prec)						:: rcut2,r2
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
rcut2 = Rcut*Rcut

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
!	      --- outside cutoff ? ---
		inside = 0
		ig_atom = cgp(ig)%first
		do while ((ig_atom .le. cgp(ig)%last) .and. (inside .eq. 0))
                   	i = cgpatom(ig_atom)
			shift = qvec_sub(x(i),x(ja))
			r2 = q_dist4(shift,q_vecscale(boxlength,nint(q_vecscale(shift,inv_boxl))))

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
subroutine nbpwlis2_box_lrf(Rcut,RLRF)
! args
real(kind=prec)					:: Rcut,RLRF
! local variables
integer						:: i,ig,jg,ia,ja,ig_atom, inside, jgr, LJ_code, j
real(kind=prec)					:: rcut2,r2
TYPE(qr_vec)					:: shift
!LRF
real(kind=prec)					:: RcLRF2
integer						:: jg_cgp, inside_LRF, is3

integer                                         ::iagrid, igrid, jgrid, kgrid, gridnum

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
  rcut2 = Rcut*Rcut
  RcLRF2 = RLRF*RLRF

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
!	      --- outside cutoff ? ---
                inside = 0
                inside_LRF = 0
                ig_atom = cgp(ig)%first
                do while ((ig_atom .le. cgp(ig)%last) .and. (inside .eq. 0))
                        i = cgpatom(ig_atom)
                        shift = qvec_sub(x(i),x(ja))
                        r2 = q_dist4(shift,q_vecscale(boxlength,nint(q_vecscale(shift,inv_boxl))))
                        if ( r2 .le. rcut2 ) then
! inside cutoff, raise the flag
                                inside = 1
!$omp critical
                                if (nbpw_cgp_pair .eq. size(nbpw_cgp, 1)) call reallocate_nbpw_cgp
                                nbpw_cgp_pair = nbpw_cgp_pair + 1
                                nbpw_cgp(nbpw_cgp_pair)%i = i
                                nbpw_cgp(nbpw_cgp_pair)%j = ja
!$omp end critical
                        elseif ((r2 .le. RcLRF2) .or. (RLRF .eq. -one)) then
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
subroutine nbpwlis2_lrf(Rcut,RLRF)
! args
real(kind=prec)					:: Rcut, RLRF
! local variables
integer						:: i,j,ig,jg,jg_cgp,ia,ja,inside,is,jgr
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
rcut2 = Rcut*Rcut
RcLRF2 = RLRF*RLRF
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
			r2 = q_dist4(x(is),x(ja))
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

subroutine nbpwlist(Rcut)
! args
real(kind=prec)					:: Rcut
! local variables
integer						:: i,ig,jg,ia,ja,jgr,j
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
rcut2 = Rcut*Rcut
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
		r2 = q_dist4(x(ia),(x(ja))

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

subroutine nbpwlist_box(Rcut)
! args
real(kind=prec)					:: Rcut
! local variables
integer						:: i,ig,jg,ia,ja,ig_sw,jg_sw,jgr,j
real(kind=prec)						:: rcut2,r2
TYPE(qr_vec)					:: shift
integer						:: LJ_code
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
rcut2 = Rcut*Rcut
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
		shift = qvec_sub(x(ig_sw),x(jg_sw))
		r2 = q_dist4(shift,q_vecscale(boxlength,nint(q_vecscale(shift,inv_boxl))))

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

subroutine nbpwlist_lrf(Rcut,RLRF)
! args
real(kind=prec)					:: Rcut,RLRF
! local variables
integer						:: i,j,ig,jg,jg_cgp,ia,ja,is,jgr
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
rcut2 = Rcut*Rcut
RcLRF2 = RLRF*RLRF
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
		r2 = q_dist4(x(is),(ja))
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
subroutine nbpwlist_box_lrf(Rcut,RLRF)
! args
real(kind=prec)					:: Rcut,RLRF
! local variables
integer						:: i,ig,jg,ia,ja,ig_sw,jg_sw,jgr
real(kind=prec)						:: rcut2,r2
integer						:: LJ_code
TYPE(qr_vec)					:: shift
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
  rcut2 = Rcut*Rcut
  RcLRF2 = RLRF*RLRF
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
		shift = qvec_sub(x(ig_sw),x(jg_sw))
		r2 = q_dist4(shift,q_vecscale(boxlength,nint(q_vecscale(shift,inv_boxl))))

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
subroutine nbqp_count(Rq,nqp, nqpcgp)
! arguments
real(kind=prec)					:: Rq
integer						:: nqp
integer						:: nqpcgp(:)

! local variables
integer						:: ig,ia,i,j,iq,i3
real(kind=prec)					:: rcut2,r2
TYPE(qr_vec)					:: shift

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
rcut2 = Rq*Rq



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
	r2 = q_dist4(x(ia),xpcent)
else
	if (Rq .gt. zero ) then
		shift = qvec_sub(x(i),x(qswitch))
		r2 = q_dist4(shift,q_vecscale(boxlength,nint(q_vecscale(shift,inv_boxl))))
	end if
end if

! skip if outside cutoff
if ( ( ( use_PBC ) .and. (r2 .gt. rcut2) .and. (Rq .gt. zero ) ) .or. & 
	( ( r2 .gt. rcut2 ) ) ) cycle igloop

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
subroutine nbqw_count(Rq,nqw, nqwmol)
! arguments
real(kind=prec)					:: Rq
integer						:: nqw
integer						:: nqwmol(:)

! local variables
integer						:: ig,ia,i,j,iq,i3
real(kind=prec)					:: rcut2,r2
TYPE(qr_vec)					:: shift

!	This routine counts water molecules that interact with q-atoms
!	Note that in PBC the default is no cutoff (rcq < 0), so we just count all of them
nqw = 0
rcut2 = Rq*Rq


if(nqat .eq. 0) return

! --- solvent - Q-atoms

iwloop: do ig = 1, nwat
nqwmol(ig) = 0
ia = nat_solute + solv_atom*ig-(solv_atom-1)
if(.not. use_PBC .and. excl(ia)) cycle iwloop ! skip excluded waters
i3 = 3*ia-3

!******PWadded if-statement 2001-10-01
if( .not. use_PBC ) then
	r2 = q_dist4(x(ia),xpcent)
else
	if (Rq .gt. zero) then
		shift = qvec_sub(x(ia),x(qswitch))
		r2 = q_dist4(shift,q_vecscale(boxlength,nint(q_vecscale(shift,inv_boxl))))
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

subroutine nbqplis2(Rq)
! args
real(kind=prec)					:: Rq
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
if(list_done .and. (Rq > rexcl_o)) return
!$omp single
nbqp_pair = 0
!$omp end single
!$omp barrier
rcut2 = Rq*Rq


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

		r2 = q_dist4(x(i),xpcent)
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

subroutine nbqplis2_box(Rq)
! args
real(kind=prec)					:: Rq
!  ! local variables
integer						:: ig,ia,i,j,iq,i3,nl,inside,ig_atom,is
real(kind=prec)					:: rcut2,r2
integer						:: xspec
logical,save					:: list_done = .false.
TYPE(qr_vec)					:: shift

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
rcut2 = Rq*Rq

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
	if (Rq.lt.zero) then
	        inside = 1
	else
	        inside = 0
	end if
	ig_atom = cgp(ig)%first
	do while ((ig_atom .le. cgp(ig)%last) .and. (inside .eq. 0))
                i = cgpatom(ig_atom)
                shift = qvec_sub(x(i),x(qswitch))
                r2 = q_dist4(shift,q_vecscale(boxlength,nint(q_vecscale(shift,inv_boxl))))
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

subroutine nbqplist(Rq)
! args
real(kind=prec)					:: Rq
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
if(list_done .and. (Rq .gt. rexcl_o)) return
if(nqat .eq. 0) return
!$omp single
nbqp_pair = 0
!$omp end single
!$omp barrier
rcut2 = Rq*Rq


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
	r2 = q_dist4(x(ia),xpcent)

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
subroutine nbqplist_box(Rq)
! args
real(kind=prec)					:: Rq
! local variables
integer						:: ig,ia,i,j,iq,i3,nl,ig_atom,is
real(kind=prec)					:: rcut2,r2
integer						::	xspec
integer						:: ga, gb, inside
logical,save					:: list_done = .false.
real(kind=prec)					:: shift
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
rcut2 = Rq*Rq
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
        if (Rcq.lt.zero) then
                inside = 1
        else
                inside = 0
        end if
        ig_atom = cgp(ig)%first
        do while ((ig_atom .le. cgp(ig)%last) .and. (inside .eq. 0))
                i = cgpatom(ig_atom)
                shift = qvec_sub(x(i),x(qswitch))
                r2 = q_dist4(shift,q_vecscale(boxlength,nint(q_vecscale(shift,inv_boxl))))
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

subroutine nbqwlist(Rq)
! args
real(kind=prec)					:: Rq
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
if( list_done .and.  (Rq .gt. rexcl_o)) return

if(nqat .eq. 0) return
!$omp single
nbqw_pair = 0
!$omp end single
!$omp barrier
rcut2 = Rq*Rq



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
	r2 = q_dist4(x(ia),xpcent)

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
subroutine nbqwlist_box(Rq)
! args
real(kind=prec)					:: Rq
! local variables
integer						:: ig,ia,i,i3,is,iq
real(kind=prec)					:: rcut2,r2
logical,save					:: list_done = .false.
TYPE(qr_vec)					:: shift
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
rcut2 = Rq*Rq

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
	if (Rq.gt.zero) then
	shift = qvec_sub(x(ia),x(qswitch))
	r2 = q_dist4(shift,q_vecscale(boxlength,nint(q_vecscale(shift,inv_boxl))))
	end if
	! store if inside cutoff
	if ( ( r2 .le. rcut2 ).or. (Rq.lt.zero)) then
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
subroutine nbww_count(Rcut,nww, nwwmol)
! arguments
real(kind=prec)					:: Rcut
integer						:: nww
integer						:: nwwmol(:)

! local variables
integer						:: iw,jw,ia,ja,i3,j3
real(kind=prec)					:: rcut2,r2
TYPE(qr_vec)					:: shift

! This routine counts non-bonded solvent-solvent atom pairs.

nww = 0
rcut2 = Rcut*Rcut

iwloop: do iw = 1, nwat
nwwmol(iw) = 0

ia = nat_solute + solv_atom*iw-(solv_atom-1)
if(.not. use_PBC .and. excl(ia)) cycle iwloop ! skip excluded waters


jwloop: do jw = 1, nwat
  ja = nat_solute + solv_atom*jw-(solv_atom-1)
  if(.not. use_PBC .and. excl(ja)) cycle jwloop ! skip excluded waters


  ! count each w-w pair once only
  if ( ((iw .gt. jw) .and. (mod(iw+jw,2) .eq. 0)) .or. &
           ((iw .lt. jw) .and. (mod(iw+jw,2) .eq. 1)) .or. &
           (iw .eq. jw)) &
           cycle jwloop

  if( use_PBC ) then
  	shift = qvec_sub(x(ia),x(ja))
  	r2 = q_dist4(shift,q_vecscale(boxlength,nint(q_vecscale(shift,inv_boxl))))
  else
        r2 = q_dist4(x(ia),x(ja))
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
subroutine nbwwlist(Rcut)
! args
real(kind=prec)					:: Rcut
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
rcut2 = Rcut*Rcut
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

		r2 = q_dist4(x(ia),x(ja))

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
subroutine nbwwlist_box(Rcut)
! args
real(kind=prec)					:: Rcut
! local variables
integer						:: iw,jw,ia,ja,i3,j3,la,ka,jwr
real(kind=prec)					:: rcut2,r2
integer                                         :: iagrid, igrid, jgrid, kgrid, gridnum
TYPE(qr_vec)					:: shift
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
rcut2 = Rcut*Rcut

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
! count each w-w pair once only
                if ( ((iw .gt. jwr) .and. (mod(iw+jwr,2) .eq. 0)) .or. &
                        ((iw .lt. jwr) .and. (mod(iw+jwr,2) .eq. 1)) .or. &
                        (iw .eq. jwr)) &
                        cycle jwloop

                shift = qvec_sub(x(ia),x(ja))
                r2 = q_dist4(shift,q_vecscale(boxlength,nint(q_vecscale(shift,inv_boxl))))
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

subroutine nbwwlist_lrf(Rcut,RLRF)
! args
real(kind=prec)					:: Rcut,RLRF
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
rcut2 = Rcut*Rcut
RcLRF2 = RLRF*RLRF

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
! count each w-w pair once only
                if ( ((iw .gt. jwr) .and. (mod(iw+jwr,2) .eq. 0)) .or. &
                        ((iw .lt. jwr) .and. (mod(iw+jwr,2) .eq. 1)) .or. &
                        (iw .eq. jwr) ) &
                        cycle jwloop

		r2 = q_dist4(x(is),x(ja))
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
subroutine nbwwlist_box_lrf(Rcut,RLRF)
! args
real(kind=prec)					:: Rcut,RLRF
! local variables
integer						:: i,j,ig,jg,iw,jw,ia,ja,i3,j3,is,is3,la,ka,jwr
real(kind=prec)					:: rcut2,r2
real(kind=prec)					::	RcLRF2
integer                                         :: iagrid, igrid, jgrid, kgrid, gridnum
TYPE(qr_vec)					:: shift
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
rcut2 = Rcut*Rcut
RcLRF2 = RLRF*RLRF

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

                shift = qvec_sub(x(is),x(ja))
                r2 = q_dist4(shift,q_vecscale(boxlength,nint(q_vecscale(shift,inv_boxl))))
                jg = iwhich_cgp(ja)
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

real(kind=prec)  :: r,r2,r6,r12,r6_hc
integer	 :: i,j,istate,par
integer	 :: atomi,atomj,cstart,cend,calcgroup
real(kind=prec)  :: V_a,V_b,Vel,Vvdw,Vwel,Vwvdw,Vwsum
TYPE(qr_vec) :: shift
TYPE(qr_dist3) :: distance

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

			if(use_PBC) then
				shift = qvec_sub(x(i),x(j))
				distance = q_dist3(shift,q_vescale(boxlength,nint(q_vecscale(shift,inv_boxl))))
			else
				distance = q_dist3(x(i),x(j))
			end if



                        r2    = distance%r2
                        r6_hc = one/distance%r6  !for softcore
                        r6    = r6_hc + monitor_group_int(calcgroup,istate)%score !softcore, default = zero
                        r6    = one/r6
                        r     = distance%r
                        r12   = r6*r6
        
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
! because it is its own module, they get passed the list they need


subroutine nonbond_pp
! local variables
integer						:: ip,i,j
real(kind=prec)						:: r2,r,r6,r12
real(kind=prec)						:: Vel,V_a,V_b,dv
TYPE(qr_dist2)					:: distance
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

distance = q_dist2(x(i),x(j))
! q_dist2 already gives distances as one/r*
! so no need for further conversion
r2    = distance%r2
r     = distance%r
r6    = distance%r6
r12   = distance%r12

! calculate Vel and dv
Vel  = nbpp(ip)%elec * r
V_a  = nbpp(ip)%vdWA *r12 
V_b  = nbpp(ip)%vdWB *r6
dv   = r2*( -Vel -12.0_prec*V_a +6.0_prec*V_b )
! update d
d(i) = d(i) - dv*distance%vec
d(j) = d(j) + dv*distance%vec

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
real(kind=prec)						:: r2,r,r6,r12
real(kind=prec)						:: Vel,V_a,V_b,dv
TYPE(qr_dist2)					:: distance
TYPE(qr_vec)					:: shift
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
shift = qvec_sub(x(ga),x(gb))

nbpp_cgp(group) = q_vecscale(boxlength,nint(q_vecscale(shift,inv_boxl)))

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
shift = qvec_sub(x(i),x(j))
distance = q_dist2(shift,nbpp_cgp(ga))


r2   = distance%r2
r    = distance%r
r6   = distance%r6
r12  = distance%r12

Vel  = nbpp(ip)%elec * r
V_a  = nbpp(ip)%vdWA *r12 
V_b  = nbpp(ip)%vdWB *r6
dv   = r2*( -Vel -12.0_prec*V_a +6.0_prec*V_b )

! update d
d(i) = d(i) - dv*distance%vec
d(j) = d(j) + dv*distance%vec

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
integer						:: ip,i,j
real(kind=prec)					:: r2,r,r6,r12
real(kind=prec)					:: Vel,V_a,V_b,dv
TYPE(qr_dist2)					:: distance

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

! calculate dx and r
distance = q_dist2(x(i),x(j))

r2   = distance%r2
r    = distance%r 
r6   = distance%r6
r12  = distance%r12

! calculate Vel and dv
Vel  = nbpw(ip)%elec *r
V_a  = nbpw(ip)%vdWA *r12 
V_b  = nbpw(ip)%vdWB *r6
dv   = r2*( -Vel -12.0_prec*V_a +6.0_prec*V_b )

! update forces
d(i) = d(i) - dv*distance%vec
d(j) = d(j) + dv*distance%vec

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
integer						:: ip,i,j
real(kind=prec)					:: r2,r,r6,r12
real(kind=prec)					:: Vel,V_a,V_b,dv
integer						:: group, ga, gb
TYPE(qr_vec)					:: shift
TYPE(qr_dist2)					:: distance
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
shift = qvec_sub(x(ga),x(gb))
nbpp_cgp(group) = q_vecscale(boxlength,nint(q_vecscale(shift,inv_boxl)))

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

! calculate dx and r
shift = qvec_sub(x(i),x(j))
distance = q_dist2(shift,nbpp_cgp(group))

r2   = distance%r2
r    = distance%r
r6   = distance%r6
r12  = distance%r12

! calculate Vel and dv
Vel  = nbpw(ip)%elec *r
V_a  = nbpw(ip)%vdWA *r12 
V_b  = nbpw(ip)%vdWB *r6
dv   = r2*( -Vel -12.0_prec*V_a +6.0_prec*V_b )
! update forces
d(i) = d(i) - dv*distance%vec
d(j) = d(j) + dv*distance%vec

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
integer						:: ip,iq,jq,i,j,k
real(kind=prec)					:: qi,qj,r2,r,r6,r12,r6_hc
real(kind=prec)					:: Vel,V_a,V_b,dv
TYPE(qr_dist3)					:: distance
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

! calculate dx and r
distance = q_dist3(x(i),x(j))

r2    = distance%r2
r6_hc = one/distance%r6
r     = distance%r
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
d(i) = d(i) - dv*distance%vec
d(j) = d(j) + dv*distance%vec
! update energies
#ifdef _OPENMP
qomp_elec(istate,1)  = qomp_elec(istate,1) + Vel
qomp_vdw(istate,1)   = qomp_vdw(istate,1) + V_a - V_b
#else
EQ(istate)%qq(1)%el  = EQ(istate)%qq(1)%el + Vel
EQ(istate)%qq(1)%vdw = EQ(istate)%qq(1)%vdw + V_a - V_b
#endif
end do ! ip

end do ! istate
#ifdef _OPENMP
do istate=1, nstates
!$omp atomic update
EQ(istate)%qq(1)%el  = EQ(istate)%qq(1)%el + qomp_elec(istate,1)
!$omp atomic update
EQ(istate)%qq(1)%vdw = EQ(istate)%qq(1)%vdw + qomp_vdw(istate,1)
end do
#endif
end subroutine nonbond_qq

!-----------------------------------------------------------------------
subroutine nonbond_qqp
! local variables
integer                                         :: istate,jj
integer                                         :: ip,iq,jq,i,j,k
real(kind=prec)                                 :: qi,qj,r2,r,r6,r12,r6_hc
real(kind=prec)                                 :: Vel,V_a,V_b,dv
TYPE(qr_dist3)					:: distance
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
! calculate dx and r
distance = q_dist3(x(i),x(j))

r2    = distance%r2
r6_hc = one/distance%r6
r     = distance%r
r6    = r6_hc + nbqq(ip,istate)%score !softcore
r6    = one/r6
r12   = r6*r6


! calculate Vel, V_a, V_b and dv
Vel  = nbqqp(ip,istate)%elec *r
V_a  = nbqqp(ip,istate)%vdWA *r12
V_b  = nbqqp(ip,istate)%vdWB *r6
dv  = r2*( -Vel -(12.0_prec*V_a -6.0_prec*V_b)*r6*r6_hc )*EQ(istate)%lambda

! update forces
d(i) = d(i) - dv*distance%vec
d(j) = d(j) + dv*distance%vec
! update energies
#ifdef _OPENMP
qomp_elec(istate,1) = qomp_elec(istate,1) + Vel
qomp_vdw(istate,1) = qomp_vdw(istate,1) + V_a - V_b
#else
EQ(istate)%qp(1)%el  = EQ(istate)%qp(1)%el + Vel
EQ(istate)%qp(1)%vdw = EQ(istate)%qp(1)%vdw + V_a - V_b
#endif
end do ! ip

end do ! istate
#ifdef _OPENMP
do istate=1, nstates
!$omp atomic update
EQ(istate)%qp(1)%el  = EQ(istate)%qp(1)%el + qomp_elec(istate,1)
!$omp atomic update
EQ(istate)%qp(1)%vdw = EQ(istate)%qp(1)%vdw + qomp_vdw(istate,1)
end do
#endif
end subroutine nonbond_qqp

!-----------------------------------------------------------------------

subroutine nonbond_qp
! local variables
integer						:: ip,iq,i,j
integer						:: istate,jj
real(kind=prec)					:: r2,r,r6,r6_hc,r12
real(kind=prec)					:: Vel,V_a,V_b,dv
TYPE(qr_dist3)					:: distance
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

distance = q_dist3(x(i),x(j))

r2    = distance%r2
r6_hc = one/distance%r6
r     = distance%r

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
d(i) = d(i) - dv*distance%vec
d(j) = d(j) + dv*distance%vec

  ! update q-protein or q-water energies
#ifdef _OPENMP
qomp_elec(istate,1) = qomp_elec(istate,1) + Vel
qomp_vdw(istate,1)  = qomp_vdw(istate,1) + V_a - V_b
#else
EQ(istate)%qp(1)%el  = EQ(istate)%qp(1)%el + Vel
EQ(istate)%qp(1)%vdw = EQ(istate)%qp(1)%vdw + V_a - V_b
#endif
end do ! istate

end do

#ifdef _OPENMP
do istate=1, nstates
!$omp atomic update
EQ(istate)%qp(1)%el  = EQ(istate)%qp(1)%el + qomp_vdw(istate,1)
!$omp atomic update
EQ(istate)%qp(1)%vdw = EQ(istate)%qp(1)%vdw + qomp_elec(istate,1)
end do
#endif
end subroutine nonbond_qp
!-----------------------------------------------------------------------

!******PWadded 2001-10-23
subroutine nonbond_qp_box
! local variables
integer						:: ip,iq,i,j
integer						:: istate,jj
real(kind=prec)					:: r2,r,r6,r6_hc,r12
real(kind=prec)					:: Vel,V_a,V_b,dv
integer						:: group, gr, ia
TYPE(qr_dist3)					:: distance
TYPE(qr_vec)					:: shift
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
shift = qvec_sub(x(ga),x(qswitch))
nbqp_cgp(gr) = q_vecscale(boxlength,nint(q_vecscale(shift,inv_boxl)))

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
group = nbqp(ip,1)%cgp_pair

shift = qvec_sub(x(i),x(j))
distance = q_dist3(shift,nbqp_cgp(group))

r2    = distance%r2
r6_hc = one/distance%r6
r     = distance%r

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
d(i) = d(i) - dv*distance%vec
d(j) = d(j) + dv*distance%vec

! update q-protein or q-water energies
#ifdef _OPENMP
qomp_elec(istate,1) = qomp_elec(istate,1) + Vel
qomp_vdw(istate,1)  = qomp_vdw(istate,1) + V_a - V_b
#else
EQ(istate)%qp(1)%el  = EQ(istate)%qp(1)%el + Vel
EQ(istate)%qp(1)%vdw = EQ(istate)%qp(1)%vdw + V_a - V_b
#endif
end do ! istate

end do
#ifdef _OPENMP
do istate=1, nstates
!$omp atomic update
EQ(istate)%qp(1)%el  = EQ(istate)%qp(1)%el + qomp_elec(istate,1)
!$omp atomic update
EQ(istate)%qp(1)%vdw = EQ(istate)%qp(1)%vdw + qomp_vdw(istate,1)
end do
#endif
end subroutine nonbond_qp_box

!-----------------------------------------------------------------------

subroutine nonbond_qw

! local variables
integer						:: jw,iq,i,j,jj,iw
integer						:: istate
real(kind=prec)					:: r2,r,r6,r6_hc,dv,Vel,V_a,V_b,r12
TYPE(qr_dist2)					:: distance
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
! init state-invariant variables:
iq   = nbqw(iw,1)%i
i = iqseq(iq)


do jw = 0, solv_atom-1
! for every assigned q-w pair:

! init state-invariant variables:
j    = nbqw(iw+jw,1)%j

distance = q_dist3(x(i),x(j))

r2    = distance%r2
r6_hc = one/distance%r6
r     = distance%r



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
d(i) = d(i) - dv*distance%vec
d(j) = d(j) + dv*distance%vec

#ifdef _OPENMP
qomp_elec(istate,1) = qomp_elec(istate,1) + Vel
qomp_vdw(istate,1)  = qomp_vdw(istate,1) + V_a - V_b
#else
EQ(istate)%qw(1)%el  = EQ(istate)%qw(1)%el + Vel
EQ(istate)%qw(1)%vdw = EQ(istate)%qw(1)%vdw + V_a - V_b 
#endif

end do ! nstates
end do !jw

end do ! nbqw
#ifdef _OPENMP
do istate=1, nstates
!$omp atomic update
EQ(istate)%qp(1)%el  = EQ(istate)%qp(1)%el + qomp_elec(istate,1)
!$omp atomic update
EQ(istate)%qp(1)%vdw = EQ(istate)%qp(1)%vdw + qomp_vdw(istate,1)
end do
#endif
end subroutine nonbond_qw

!-----------------------------------------------------------------------
!******PWadded 2001-10-23
subroutine nonbond_qw_box
	! local variables
integer                                         :: jw,iq,i,j,jj,iw
integer						:: istate
real(kind=prec)                                 :: r2,r,r6,r6_hc,dv,Vel,V_a,V_b,r12
TYPE(qr_dist3)					:: distance
TYPE(qr_vec)					:: shift
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
do jw = 0, solv_atom-1

! init state-invariant variables:
j    = nbqw(iw+jw,1)%j

!compute the periodical shift
shift = qvec_sub(x(qswitch),x(j))
shift = q_vecscale(boxlength,nint(q_vecscale(shift,inv_boxl)))

distance = q_dist3(qvec_sub(x(i),x(j)),shift)

r2    = distance%r2
r6_hc = one/distance%r6
r     = distance%r



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
d(i) = d(i) - dv*distance%vec
d(j) = d(j) + dv*distance%vec

#ifdef _OPENMP
qomp_elec(istate,1) = qomp_elec(istate,1) + Vel
qomp_vdw(istate,1)  = qomp_vdw(istate,1) + V_a - V_b
#else
EQ(istate)%qw(1)%el  = EQ(istate)%qw(1)%el + Vel
EQ(istate)%qw(1)%vdw = EQ(istate)%qw(1)%vdw + V_a - V_b
#endif

end do ! nstates
end do ! jw

end do ! nbqw
#ifdef _OPENMP
do istate=1, nstates
!$omp atomic update
EQ(istate)%qp(1)%el  = EQ(istate)%qp(1)%el + qomp_elec(istate,1)
!$omp atomic update
EQ(istate)%qp(1)%vdw = EQ(istate)%qp(1)%vdw + qomp_vdw(istate,1)
end do
#endif
end subroutine nonbond_qw_box

!----------------------------------------------------------------------------------------------------

subroutine nonbond_ww
! local variables
integer						:: iw,ip,i,j,ia,ichg
integer						:: ja
integer						:: ipstart
real(kind=prec)					:: r2,r,r6,r12
real(kind=prec)					:: Vel,V_a,V_b,dv
TYPE(qr_dist2)					:: distance
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
do ip = 0, solv_atom - 1
j    = nbww(iw+ip)%j

distance = q_dist2(x(i),x(j))


r2   = distance%r2
r    = distance%r 
r6   = distance%r6
r12  = distance%r12
Vel  = ww_precomp(ichg,ip+1)%elec *r
V_a  = ww_precomp(ichg,ip+1)%vdWA *r12 
V_b  = ww_precomp(ichg,ip+1)%vdWB *r6
dv   = r2*( -Vel -12.0_prec*V_a +6.0_prec*V_b )
d(i) = d(i) - dv*distance%vec
d(j) = d(j) + dv*distance%vec
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
integer						:: iw,ip,i,j,ia,ichg
integer						:: ja
integer						:: ipstart
real(kind=prec)					:: r2,r,r6,r12
real(kind=prec)					:: Vel,V_a,V_b,dv
TYPE(qr_dist2)					:: distance
TYPE(qr_vec)					:: shift
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
j    = nbww(iw)%j

!distance between this oxygen atom and the oxygen atom of the above watermolecule, iw
!this is always the first interaction :P
!only need to calculate this for ichg = 1
if (ichg .eq.1) then
! periodic shift
shift = qvec_sub(x(i),x(j))
shift = q_vecscale(boxlength,nint(q_vecscale(shift,inv_boxl)))
end if

do ip = 0, solv_atom - 1
j    = nbww(iw+ip)%j

distance = q_dist2(qvec_sub(x(i),x(j)),shift)

r2   = distance%r2
r    = distance%r 
r6   = distance%r6
r12  = distance%r12

Vel  = ww_precomp(ichg,ip+1)%elec *r
V_a  = ww_precomp(ichg,ip+1)%vdWA *r12
V_b  = ww_precomp(ichg,ip+1)%vdWB *r6
dv   = r2*( -Vel -12.0_prec*V_a +6.0_prec*V_b )

d(i) = d(i) - dv*distance%vec
d(j) = d(j) + dv*distance%vec

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
integer						:: iw,jw,iq,i,j,jj,ip,ja
integer						:: istate
real(kind=prec)					:: r, r2, r6, r12
real(kind=prec)					:: Vel, dv
real(kind=prec)					:: V_a, V_b, r6_hc
TYPE(qr_dist3)					:: distance
TYPE(qr_dist)					:: eldist
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
j    = nbqw(iw,1)%j

! only q - O distance first, only one with vdW
distance = q_dist3(x(i),x(j))

r2    = distance%r2
r6_hc = one/distance%r6
r     = distance%r


dv = zero
do istate = 1, nstates
r6    = r6_hc + nbqw(iw,istate)%score !softcore
r6    = one/r6
r12   = r6*r6
! calculate qi, Vel, V_a, V_b and dv
V_a = nbqw(iw,istate)%vdWA *r12
V_b = nbqw(iw,istate)%vdWB *r6
Vel = nbqw(iw,istate)%elec *r
dv  = dv  + r2*( -Vel -( (12.0_prec*V_a - 6.0_prec*V_b)*r6*r6_hc ))*EQ(istate)%lambda
!softcore r6*r6_hc is (r^6/(r^6+alpha))
! update q-water energies
#ifdef _OPENMP
qomp_elec(istate,1) = qomp_elec(istate,1) + Vel
qomp_vdw(istate,1)  = qomp_vdw(istate,1) + V_a - V_b
#else
EQ(istate)%qw(1)%el  = EQ(istate)%qw(1)%el + Vel
EQ(istate)%qw(1)%vdw = EQ(istate)%qw(1)%vdw + V_a - V_b
#endif
end do !istate
! update forces on q atom
d(i) = d(i) - dv*distance%vec
! update forces on water
d(j) = d(j) + dv*distance%vec

!now calculate only charge-charge interaction for hydrogens
do ip = 1, solv_atom-1

ja   = j + ip

eldist = q_dist(x(i),x(ja))

r2   = eldist%r2
r    = eldist%r

dv = 0
do istate = 1, nstates
Vel = nbqw(iw+ip,istate)%elec *r
#ifdef _OPENMP
qomp_elec(istate,1) = qomp_elec(istate,1) + Vel
#else
EQ(istate)%qw(1)%el  = EQ(istate)%qw(1)%el + Vel
#endif
dv = dv - r2*Vel*EQ(istate)%lambda
end do ! istate
! update forces on q atom
d(i) = d(i) - dv*eldist%vec
! update forces on water
d(ja) = d(ja) + dv*eldist%vec
end do ! ip

end do ! nbqw
#ifdef _OPENMP
do istate = 1, nstates
!$omp atomic update
EQ(istate)%qw(1)%vdw = EQ(istate)%qw(1)%vdw + qomp_vdw(istate,1)
!$omp atomic update
EQ(istate)%qw(1)%el  = EQ(istate)%qw(1)%el + qomp_elec(istate,1)
end do ! nstates
#endif
end subroutine nonbond_qw_spc

!-----------------------------------------------------------------------
!******PWadded 2001-10-23
subroutine nonbond_qw_spc_box
	!calculate non-bonded interactions between Q-atoms and SPC water molecules
	!(optimisations rely on LJ params = 0 for water H) using geometric comb. rule

	! local variables
integer                                         :: iw,iq,i,j,iLJ,jj,ip,ja
integer						:: istate
real(kind=prec)                                 :: r, r2, r6, r12
real(kind=prec)                                 :: Vel, dv
real(kind=prec)                                 :: V_a, V_b, r6_hc
TYPE(qr_dist3)					:: distance
TYPE(qr_dist)					:: eldist
TYPE(qr_vec)					:: shift
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
j    = nbqw(iw,1)%j

shift = qvec_sub(x(qswitch),x(j))
shift = q_vecscale(boxlength,nint(q_vecscale(shift,inv_boxl)))
distance = q_dist3(qvec_sub(x(i),x(j)),shift)

r2    = distance%r2
r6_hc = one/distance%r6
r     = distance%r

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
qomp_elec(istate,1) = qomp_elec(istate,1) + Vel
qomp_vdw(istate,1)  = qomp_vdw(istate,1) + V_a - V_b
#else
EQ(istate)%qw(1)%el  = EQ(istate)%qw(1)%el + Vel
EQ(istate)%qw(1)%vdw = EQ(istate)%qw(1)%vdw + V_a - V_b
#endif
end do !istate			
! update forces on q atom
d(i) = d(i) - dv*distance%vec
! update forces on water
d(j) = d(j) + dv*distance%vec

!now calculate only charge-charge interaction for hydrogens
do ip = 1, solv_atom-1

ja   = j + ip

eldist = q_dist2(qvec_sub(x(i),x(ja)),shift)

r2   = eldist%r2
r    = eldist%r

dv = 0
do istate = 1, nstates
Vel = nbqw(iw+ip,istate)%elec *r
#ifdef _OPENMP
qomp_elec(istate,1) = qomp_elec(istate,1) + Vel
#else
EQ(istate)%qw(1)%el  = EQ(istate)%qw(1)%el + Vel
#endif
dv = dv - r2*Vel*EQ(istate)%lambda
end do ! istate
! update forces on q atom
d(i) = d(i) - dv*eldist%vec
! update forces on water
d(ja) = d(ja) + dv*eldist%vec
end do ! ip

end do ! nbqw
#ifdef _OPENMP
do istate = 1, nstates
!$omp atomic update
EQ(istate)%qw(1)%vdw = EQ(istate)%qw(1)%vdw + qomp_vdw(istate,1)
!$omp atomic update
EQ(istate)%qw(1)%el  = EQ(istate)%qw(1)%el + qomp_elec(istate,1)
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
integer                                         :: iw,ip,i,j,ia,ichg,ja,stat,ic
integer						:: ipstart
real(kind=prec)                                 :: r2,r,r6,r12
real(kind=prec)					:: Vel,V_a,V_b,dv
TYPE(qr_dist2)					:: distance
TYPE(qr_dist)					:: eldist
#ifdef _OPENMP
integer :: quotient, remainder
real(kind=prec) :: Vel_omp, Vvdw_omp
#endif

! the loop is now organised in the following way (as before for the arithmetic case)
! we always jump solv_atom positions, to calculate atomAn-atomB(1-solv_atom)
! and iterate a separate counter to get the current position of the atomA(1-solv_atom) position
! vdW is only calculated if atomA1 - atomB1

#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (nbww_pair/solv_atom**2)/threads_num
remainder = MOD(nbww_pair/solv_atom**2, threads_num)
mp_start = thread_id * quotient + 1 + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
mp_start = (mp_start*(solv_atom**2)) - (solv_atom**2) + 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif
mp_end = (mp_end*(solv_atom**2)) 

Vel_omp = zero
Vvdw_omp = zero
do iw = mp_start,mp_end,solv_atom**2
#else
do iw = 1, nbww_pair, solv_atom**2
#endif
ichg=1
i    = nbww(iw)%i
j    = nbww(iw)%j

! very first one is the only calc with vdW
distance = q_dist2(x(i),x(j))
r2  = distance%r2
r   = distance%r
r6  = distance%r6
r12 = distance%r12
V_a = ww_precomp(ichg,1)%vdWA * r12
V_b = ww_precomp(ichg,1)%vdWB * r6
Vel = ww_precomp(ichg,1)%elec * r

#ifdef _OPENMP
Vvdw_omp = Vvdw_omp + V_a - V_b
Vel_omp = Vel_omp + Vel
#else
E%ww%vdw = E%ww%vdw + V_a - V_b
E%ww%el = E%ww%el + Vel
#endif
dv = r2*(-12.0_prec*V_a +6.0_prec*V_b ) + r2*( -Vel)
d(i) = d(i) - dv*distance%vec
d(j) = d(j) + dv*distance%vec


! now calc the rest of the interactions
! remember that we need to switch the atoms every solv_atom steps
! so ip needs to be reset with stat
! ic counts where we are globally to set the array index right
stat = 1
ic   = 0
do while (ichg .le. solv_atom )
do ip = stat, solv_atom - 1
i = nbww(iw+ic+ip)%i
j = nbww(iw+ic+ip)%j

eldist = q_dist(x(i),x(j))
r2 = eldist%r2
r  = eldist%r
Vel  = ww_precomp(ichg,ip+1)%elec *r
#ifdef _OPENMP
Vel_omp = Vel_omp + Vel
#else
E%ww%el = E%ww%el + Vel
#endif

dv = r2*( -Vel)
d(i) = d(i) - dv*eldist%vec
d(j) = d(j) + dv*eldist%vec

ic = ic + 1
end do
stat = 0
ichg = ichg + 1
end do
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
integer                                         :: iw,ip,i,j,ia,ichg,ja,stat,ic
integer						:: ipstart
real(kind=prec)					:: r2,r,r6,r12
real(kind=prec)					:: Vel,V_a,V_b,dv
TYPE(qr_dist2)					:: distance
TYPE(qr_dist)					:: eldist
TYPE(qr_vec)					:: shift
#ifdef _OPENMP
integer :: quotient, remainder
real(kind=prec) :: Vel_omp, Vvdw_omp
#endif

! the loop is now organised in the following way (as before for the arithmetic case)
! we always jump solv_atom positions, to calculate atomAn-atomB(1-solv_atom)
! and iterate a separate counter to get the current position of the atomA(1-solv_atom) position
! vdW is only calculated if atomA1 - atomB1

#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (nbww_pair/solv_atom**2)/threads_num
remainder = MOD(nbww_pair/solv_atom**2, threads_num)
mp_start = thread_id * quotient + 1 + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
mp_start = (mp_start*(solv_atom**2)) - (solv_atom**2) + 1
if (remainder .gt. thread_id) then
mp_end = mp_end + 1
endif
mp_end = (mp_end*(solv_atom**2)) 

Vel_omp = zero
Vvdw_omp = zero

do iw = mp_start,mp_end,solv_atom
#else
do iw = 1, nbww_pair, solv_atom
#endif
ichg=1
i    = nbww(iw)%i
j    = nbww(iw)%j


! very first one is the only calc with vdW
! and the one defining the box shift
shift = qvec_sub(x(i),x(j))
shift = q_vecscale(boxlength,nint(q_vecscale(shift,inv_boxl)))

distance = q_dist2(qvec_sub(x(i),x(j)),shift)
r2  = distance%r2
r   = distance%r
r6  = distance%r6
r12 = distance%r12
V_a = ww_precomp(ichg,1)%vdWA * r12
V_b = ww_precomp(ichg,1)%vdWB * r6
Vel = ww_precomp(ichg,1)%elec * r

#ifdef _OPENMP
Vvdw_omp = Vvdw_omp + V_a - V_b
Vel_omp = Vel_omp + Vel
#else
E%ww%vdw = E%ww%vdw + V_a - V_b
E%ww%el = E%ww%el + Vel
#endif
dv = r2*(-12.0_prec*V_a +6.0_prec*V_b ) + r2*( -Vel)
d(i) = d(i) - dv*distance%vec
d(j) = d(j) + dv*distance%vec


! now calc the rest of the interactions
! remember that we need to switch the atoms every solv_atom steps
! so ip needs to be reset with stat
! ic counts where we are globally to set the array index right
stat = 1
ic   = 0
do while (ichg .le. solv_atom )
do ip = stat, solv_atom - 1
i = nbww(iw+ic+ip)%i
j = nbww(iw+ic+ip)%j

eldist = q_dist(qvec_sub(x(i),x(j)),shift)
r2 = eldist%r2
r  = eldist%r
Vel  = ww_precomp(ichg,ip+1)%elec *r
#ifdef _OPENMP
Vel_omp = Vel_omp + Vel
#else
E%ww%el = E%ww%el + Vel
dv = r2*( -Vel)
d(i) = d(i) - dv*eldist%vec
d(j) = d(j) + dv*eldist%vec

ic = ic + 1
end do
stat = 0
ichg = ichg + 1
end do
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
integer                                         :: iw,is,ip,ia,ib
real(kind=prec)					:: r2,r6,r,r12
real(kind=prec)					:: Vel,V_a,V_b,dv
TYPE(qr_dist2)					:: distance
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

distance = q_dist2(x(ia),x(ib))

r2  = distance%r2 
r   = distance%r 
r6  = distance%r6 
r12 = distance%r12 

Vel  = nonbnd_solv_int(ip)%elec * r
V_a  = nonbnd_solv_int(ip)%vdWA * r12
V_b  = nonbnd_solv_int(ip)%vdWB * r6

dv   = r2*( -Vel) + r2*(-12.0_prec*V_a +6.0_prec*V_b )
#ifdef _OPENMP
Vel_omp  = Vel_omp + Vel
Vvdw_omp = Vvdw_omp + V_a - V_b
#else
E%ww%el  = E%ww%el + Vel
E%ww%vdw = E%ww%vdw + V_a - V_b
#endif
d(ia) = d(ia) - dv*distance%vec
d(ib) = d(ib) + dv*distance%vec
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
integer                                         :: iw,is,ip,ia,ib
real(kind=prec)                                 :: r2,r6,r,r12
real(kind=prec)                                 :: Vel,V_a,V_b,dv
TYPE(qr_dist2)					:: distance
TYPE(qr_vec)					:: shift
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

shift = qvec_sub(x(ia),x(ib))
shift = q_vecscale(boxlength,nint(q_vecscale(shift,inv_boxl)))
distance = q_dist2(qvec_sub(x(ia),x(ib)),shift)

r2  = distance%r2
r   = distance%r
r6  = distance%r6
r12 = distance%r12

Vel  = nonbnd_solv_int(ip)%elec * r
V_a  = nonbnd_solv_int(ip)%vdWA * r12
V_b  = nonbnd_solv_int(ip)%vdWB * r6

dv   = r2*( -Vel) + r2*(-12.0_prec*V_a +6.0_prec*V_b )
#ifdef _OPENMP
Vel_omp  = Vel_omp + Vel
Vvdw_omp = Vvdw_omp + V_a - V_b
#else
E%ww%el  = E%ww%el + Vel
E%ww%vdw = E%ww%vdw + V_a - V_b
#endif
d(ia) = d(ia) - dv*distance%vec
d(ib) = d(ib) + dv*distance%vec
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
integer						:: io,i,j,k,l
real(kind=prec)					:: r
TYPE(qr_dist)					:: distance

! global variables used:
!  offd, noffd, iqseq, x, Hij, offd2

do io = 1, noffd
! for every offd:

i  = offd(io)%i
j  = offd(io)%j
k  = iqseq(offd2(io)%k)
l  = iqseq(offd2(io)%l)
distance = q_dist(x(k),x(l))

r = distance%r

Hij(i,j) = offd2(io)%A * exp(-offd2(io)%mu*r)
offd(io)%Hij = Hij(i,j)	! store for save
offd(io)%rkl = r
end do
end subroutine offdiag

!-----------------------------------------------------------------------
subroutine p_restrain
! *** Local variables
integer :: ir,i,j,k,istate, n_ctr
real(kind=prec) :: fk,r2,erst,Edum,wgt,b,db,dv, totmass 
real(kind=prec)						::	fexp
TYPE(qr_vec)	:: dr,shift
TYPE(bond_val)	:: distres
TYPE(angl_val)	:: anglres
! global variables used:
!  E, nstates, EQ, nrstr_seq, rstseq, heavy, x, xtop, d, nrstr_pos, rstpos, nrstr_dist, 
!  rstdis, nrstr_wall, rstang, nrstr_ang, rstwal, xwcent, excl, freeze

! sequence restraints (independent of Q-state)
do ir = 1, nrstr_seq
fk = rstseq(ir)%fk

if(rstseq(ir)%to_centre == 1) then     ! Put == 1, then equal to 2
  ! restrain to geometrical centre

  ! reset dr & atom counter
  dr    = zero
  n_ctr = 0

  ! calculate deviation from center
  do i = rstseq(ir)%i, rstseq(ir)%j
        if ( heavy(i) .or. rstseq(ir)%ih .eq. 1 .and.(.not.(excl(i).and.freeze))) then
          n_ctr = n_ctr + 1
          dr = qvec_add(dr,qvec_sub(xtop(i),x(i)))
        end if
  end do

  if(n_ctr > 0) then 
    ! only if atoms were found:

        ! form average
        dr = dr / n_ctr 
        r2      = qvec_square(dr)
        erst    = 0.5_prec*fk*r2
        E%restraint%protein  = E%restraint%protein + erst

        ! apply same force to all atoms
        do i = rstseq(ir)%i, rstseq(ir)%j
        if ( heavy(i) .or. rstseq(ir)%ih .eq. 1 .and.(.not.(excl(i).and.freeze))) then
                d(i) = qvec_add(d(i),fk*dr*iaclib(iac(i))%mass/12.010_prec)
          end if
        end do
  end if
 
else if(rstseq(ir)%to_centre == 2) then     ! Put == 1, then equal to 2
  ! restrain to mass centre
  ! reset dr & variable to put masses
  dr = zero
  totmass = zero
  
! calculate deviation from mass center
  do i = rstseq(ir)%i, rstseq(ir)%j
        if ( heavy(i) .or. rstseq(ir)%ih .eq. 1 .and.(.not.(excl(i).and.freeze))) then
          totmass = totmass + iaclib(iac(i))%mass                              ! Add masses
          dr = qvec_add(dr,iaclib(iac(i))%mass*qvec_sub(xtop(i),x(i)))
		end if
  end do

 if(totmass > zero) then 
    ! only if atoms were found: (i.e has a total mass)

        ! form average
        dr = dr/totmass                                  ! divide by total mass
        r2      = qvec_square(dr)
        erst    = 0.5_prec*fk*r2
        E%restraint%protein  = E%restraint%protein + erst

        ! apply same force to all atoms
        do i = rstseq(ir)%i, rstseq(ir)%j
        if ( heavy(i) .or. rstseq(ir)%ih .eq. 1 .and.(.not.(excl(i).and.freeze))) then
                d(i) = qvec_add(d(i),fk*dr)
          end if
        end do
  end if

else 
  ! restrain each atom to its topology co-ordinate
  do i = rstseq(ir)%i, rstseq(ir)%j
        if ( heavy(i) .or. rstseq(ir)%ih .eq. 1 .and.(.not.(excl(i).and.freeze))) then



          dr = qvec_sub(xtop(i),x(i))

       !use the periodically minimal distance:
       if( use_PBC ) then
               dr = qvec_sub(q_vecscale(boxlength,nint(q_vecscale(dr,inv_boxl))),dr)
       end if
          r2      = qvec_square(dr)

          erst    = 0.5_prec*fk*r2
          E%restraint%protein  = E%restraint%protein + erst
          d(i) = qvec_add(d(i),fk*dr)
        end if
  end do
end if
end do

! extra positional restraints (Q-state dependent)
do ir = 1, nrstr_pos
istate = rstpos(ir)%ipsi
i      = rstpos(ir)%i

dr = qvec_sub(rstpos(ir)%x,x(i))

if ( istate .ne. 0 ) then
wgt = EQ(istate)%lambda
else
wgt = one
end if

Edum = 0.5_prec * q_dotprod(rstpos(ir)%fk,q_vecscale(dr,dr))

d(i) = qvec_add(d(i),wgt*q_vecscale(rstpos(ir)%fk,dr))


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

! if PBC then adjust lengths according to periodicity - MA
if( use_PBC ) then
	dr = qvec_sub(x(i),x(j))
	distres = bond_calc(dr,q_vecscale(boxlength,nint(q_vecscale(dr,inv_boxl))))
else
	distres = bond_calc(x(i),x(j))
end if

if ( istate .ne. 0 ) then
wgt = EQ(istate)%lambda
else
wgt = one
end if



if(distres%dist < rstdis(ir)%d1) then !shorter than d1
        db     = distres%dist - rstdis(ir)%d1
elseif(b > rstdis(ir)%d2) then !longer than d2
        db     = distres%dist - rstdis(ir)%d2
else
        db = 0
        cycle !skip zero force calculation
endif

Edum   = 0.5_prec*rstdis(ir)%fk*db**2
dv     = wgt*rstdis(ir)%fk*db/distres%dist

d(i) = qvec_add(d(i),dv*distres%avec)
d(j) = qvec_add(d(i),dv*distres%bvec)

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

! if PBC then adjust lengths according to periodicity - MA
if( use_PBC ) then
	anglres = box_angle_calc((x(i),x(j),x(k),boxlength,inv_boxl)
else
	anglres = angle_calc(x(i),x(j),x(k))

end if

if ( istate .ne. 0 ) then
wgt = EQ(istate)%lambda
else
wgt = one
end if

db    = anglres%angl - (rstang(ir)%ang)*deg2rad

! dv is the force to be added in module
Edum   = 0.5_prec*rstang(ir)%fk*db**2
dv     = wgt*rstang(ir)%fk*db


        ! update d
d(i) = qvec_add(d(i),dv*anglres%a_vec)
d(j) = qvec_add(d(j),dv*anglres%b_vec)
d(k) = qvec_add(d(k),dv*anglres%c_vec)

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

			distres = bond_calc(x(i),xwcent)
                        db = distres%dist - rstwal(ir)%d

                        if(db > zero) then
                                erst =  0.5_prec * fk * db**2 - rstwal(ir)%Dmorse
                                dv = fk*db/distres%dist                        else
                                fexp = exp(rstwal(ir)%aMorse*db)
                                erst = rstwal(ir)%dMorse*(fexp*fexp-2.0_prec*fexp)
                                dv=-2.0_prec*rstwal(ir)%dMorse*rstwal(ir)%aMorse*(fexp-fexp*fexp)/distres%dist
                        end if
                        E%restraint%protein = E%restraint%protein + erst
			d(i) = qvec_add(d(i),dv*distres%a_vec)
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
integer(4)			:: i,nat3
logical				:: old_restart = .false.
TYPE(qr_vec),allocatable	:: x_old(:),v_old(:)
TYPE(qr_vec)			:: old_boxlength, old_boxcentre
integer				:: headercheck,myprec
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
         read (12,err=112,end=112) nat3, (x_single(i),i=1,nat_pro)
         xtop(1:nat_pro) = x_single(1:nat_pro)
         deallocate(x_single)
       else if (headercheck .eq. -1337) then
         allocate(x_double(natom))
         read (12,err=112,end=112) nat3, (x_double(i),i=1,nat_pro)
         xtop(1:nat_pro) = x_double(1:nat_pro)
         deallocate(x_double)
       else if (headercheck .eq. -13337) then
#ifndef PGI
         allocate(x_quad(natom))
         read (12,err=112,end=112) nat3, (x_quad(i),i=1,nat_pro)
         xtop(1:nat_pro) = x_quad(1:nat_pro)
         deallocate(x_quad)
#else
         call die('Quadruple precision not supported in PGI')
#endif
       end if
     else
       read (12,err=112,end=112) nat3, (xtop(i),i=1,nat_pro)
     end if
  else
     allocate(x_old(natom))
     read (12,err=112,end=112) nat3, (x_old(i),i=1,nat_pro)
     xtop(1:nat_pro) = x_old(1:nat_pro)
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
                 read(2,err=112,end=112) boxl_single(:)
                 read(2,err=112,end=112) boxc_single(:)
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
                 read(2,err=112,end=112) boxl_double(:)
                 read(2,err=112,end=112) boxc_double(:)
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
                 read(2,err=112,end=112) boxl_quad(:)
                 read(2,err=112,end=112) boxc_quad(:)
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
        x(1:nat_pro) = xtop(1:nat_pro)
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

subroutine restrain_solvent 
! local variables
integer						::	iw,i,isolv,jsolv
real(kind=prec)						::	b,db,erst,dv,fexp
TYPE(qr_vec)						::	dr,dcent
real(kind=prec)						::	shift

! global variables used:
!  E, Boltz, Tfree, fk_wsphere, nwat, nat_pro, x, xwcent, rwat, Dwmz, awmz, d

if(fk_wsphere .ne. zero) then
        shift = sqrt (Boltz*Tfree/fk_wsphere)
else
        shift = zero
end if

! make decision here to use old or new code, aviod branch points during loop
if(solv_atom.eq.3) then
do iw = ncgp_solute + 1, ncgp
        i  = cgp(iw)%iswitch
        if (excl(i)) cycle ! skip excluded topology waters
        dr = qvec_sub(x(i),xwcent)
        b = qvec_square(dr)
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
        d(i) = qvec_add(d(i),dv*dr)
        E%restraint%solvent_radial = E%restraint%solvent_radial + erst
else
do iw = ncgp_solute + 1, ncgp
        i  = cgp(iw)%iswitch
        if (excl(i)) cycle ! skip excluded topology waters
        dcent=zero
        jsolv = iw - ncgp_solute
        i = nat_solute + solv_atom*jsolv - (solv_atom-1)
        do isolv = 0 , solv_atom - 1
                dcent = qvec_add(dcent,x(i+isolv))
        end do
        dcent = dcent/solv_atom
        dr = qvec_sub(dcent,xwcent)
        b = qvec_square(dr)
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
        dv = dv / solv_atom

        ! update energy and forces
        E%restraint%solvent_radial = E%restraint%solvent_radial + erst
        jsolv = iw - ncgp_solute
        i = nat_solute + solv_atom*jsolv - (solv_atom-1)
        do isolv = 0 , solv_atom - 1
               d(i+isolv) = qvec_add(d(i+isolv),dv*dr) 
        end do
end do
end if
end subroutine restrain_solvent

!-----------------------------------------------------------------------

subroutine wat_sphere
! local variables
integer					:: i,i3,kr,isort,int_wat,istate
real(kind=prec)					:: rc,rnwat
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




