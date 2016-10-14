! (C) 2014 Uppsala Molekylmekaniska HB, Uppsala, Sweden
! bonded.f90
! based on md.f90
! by Johan Åqvist, John Marelius, Anders Kaplan, Isabella Feierberg, Martin Nervall & Martin Almlöf
! calculation of all bonded interactions

module BONDENE

! used modules
!use PROFILING
use SIZES
use TRJ
use QBOND
!$ use omp_lib
implicit none

real(kind=prec) function angle(angles,istart, iend)
! *** arguments
integer						::	istart, iend
TYPE(ANG_TYPE),pointer				::	angles(:)

! *** local variables
TYPE(angl_val)					::	calc
integer						::	i,j,k,ia,ic
real(kind=prec)					::	da,dv

#ifdef _OPENMP
integer :: quotient, remainder
#endif
! global variables used:
! ang, x, anglib, d

! calculate the total energy of all protein or water angles, depending
! updates d

! reset Eangle
angle = zero			!zero = 0.0_prec

!$omp parallel default(none) shared(threads_num, istart, iend, angles, x, anglib, d, angle) private(remainder, quotient,i,j,k,ia,ic,da,dv)
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
        i  = angles(ia)%i
        j  = angles(ia)%j
        k  = angles(ia)%k
        ic = angles(ia)%cod
! use new function to calculate angle
! and get all variables
	calc = angle_calc(x(i),x(j),x(k))
! get da and dv
	da = calc%angl - anglib(ic)%ang0
#ifdef _OPENMP
    mp_real_tmp = mp_real_tmp + 0.5_prec*anglib(ic)%fk*da**2
#else
    angle = angle + 0.5_prec*anglib(ic)%fk*da**2
#endif
	dv = anglib(ic)%fk*da
! calculate forces from function vectors
! and update force array
d(i) = d(i) + dv*calc%a_vec
d(j) = d(j) + dv*calc%b_vec
d(k) = d(k) + dv*calc%c_vec

end do
#ifdef _OPENMP
!$omp atomic update
angle = angle + mp_real_tmp
#endif
!$omp end parallel

end function angle

!-----------------------------------------------------------------------
real(kind=prec) function urey_bradley(angles,istart, iend)
! *** arguments
integer						::	istart, iend
TYPE(ANG_TYPE),pointer                          ::      angles(:)
! *** local variables
TYPE(bond_val)					::	calc
integer						::	i,k,ia,ic
real(kind=prec)					::	ru,du 

#ifdef _OPENMP
integer :: quotient, remainder
#endif
! global variables used:
! ang, x, anglib, d

! calculate the total energy of all protein or water angles, depending
! updates d
! reset energy
urey_bradley = zero                    !zero = 0.0_prec

!$omp parallel default(none) shared(threads_num, istart, iend, angles, x, anglib, d, urey_bradley) private(remainder, quotient,i,k,ia,ic,calc,ru,du)
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
! for each angle in range and not zero:
	ic = angles(ia)%cod
	if(anglib(ic)%ureyfk .gt. zero) then
		i  = angles(ia)%i
		k  = angles(ia)%k
		calc = bonc_calc(x(i),x(k))
		ru = calc%dist - anglib(ic)%ureyr0
#ifdef _OPENMP
		mp_real_tmp = mp_real_tmp + anglib(ic)%ureyfk*ru**2
#else 
                urey_bradley = urey_bradley + anglib(ic)%ureyfk*ru**2
#endif
                du = 2*anglib(ic)%ureyfk*ru/calc%dist
		d(k) = d(k) + du*calc%a_vec
		d(i) = d(i) + du*calc%b_vec
        end if
end do
#ifdef _OPENMP
!$omp atomic update
urey_bradley = urey_bradley + mp_real_tmp
#endif
!$omp end parallel

end function urey_bradley

!-----------------------------------------------------------------------

real(kind=prec) function bond(bonds,istart, iend)
! *** arguments
integer						::	istart, iend
TYPE(BND_LYB),pointer				::	bonds(:)
! *** local variables
integer						::	i,j,ib,ic
TYPE(bond_val)					::	calc
real(kind=prec)					::	db,dv

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

        i  = bonds(ib)%i
        j  = bonds(ib)%j
        ic = bonds(ib)%cod
	calc = bond_calc(x(i),x(j))
        db = calc%dist - bondlib(ic)%bnd0
#ifdef _OPENMP
    mp_real_tmp = mp_real_tmp + 0.5_prec*bondlib(ic)%fk*db**2
#else  
        bond = bond + 0.5_prec*bondlib(ic)%fk*db**2
#endif

        ! calculate dv and update d
        dv = bondlib(ic)%fk*db/calc%dist
	d(i) = d(i) + calc%a_vec*dv
	d(j) = d(j) + calc%b_vec*dv
end do
#ifdef _OPENMP
!$omp atomic update
bond = bond + mp_real_tmp
#endif
!$omp end parallel

end function bond
!-----------------------------------------------------------------------

real(kind=prec) function improper(impropers,istart, iend)
!arguments
integer						::	istart, iend
TYPE(TOR_TYPE),pointer				::	impropers
! evaluate harmonic impropers
! local variables
integer						::	ip,i,j,k,l,ic
TYPE(tors_val)					::	calc
real(kind=prec)					::	dv,arg
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
	i  = impropers(ip)%i
	j  = impropers(ip)%j
	k  = impropers(ip)%k
	l  = impropers(ip)%l
	ic = impropers(ip)%cod
	calc = improper_calc(x(i),x(j),x(k),x(l)) 

! ---       energy

	arg = phi - implib(ic)%imp0
	arg = arg - 2.0_prec*pi*nint(arg/(2.0_prec*pi))
	dv  = implib(ic)%fk*arg
	#ifdef _OPENMP
	    mp_real_tmp = mp_real_tmp + 0.5_prec*dv*arg
	#else 
	    improper = improper + 0.5_prec*dv*arg
	#endif

! ---       forces
	d(i) = d(i) + dv*calc%a_vec
	d(j) = d(j) + dv*calc%b_vec
	d(k) = d(k) + dv*calc%c_vec
	d(l) = d(l) + dv*calc%d_vec

end do
#ifdef _OPENMP
!$omp atomic update
improper = improper + mp_real_tmp
#endif
!$omp end parallel
end function improper

!-----------------------------------------------------------------------

real(kind=prec) function improper2(impropers,istart, iend)
!evaluate periodic impropers
!arguments
integer						::	istart, iend
TYPE(TOR_LIB),pointer				::	impropers
! local variables
integer						::	ip,i,j,k,l,ic
real(kind=prec)					::	dv,arg
TYPE(tors_cval)					::	calc
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
	i  = impropers(ip)%i
	j  = impropers(ip)%j
	k  = impropers(ip)%k
	l  = impropers(ip)%l
	ic = impropers(ip)%cod
	calc = improper_calc(x(i),x(j),x(k),x(l))

	arg = 2*phi - implib(ic)%imp0
#ifdef _OPENMP
	mp_real_tmp = mp_real_tmp + lib%fk * (1 + cos(arg))
#else
	improper2 = improper2 + lib%fk * (1 + cos(arg))
#endif
	dv  = -2*implib(ic)%fk * sin(arg)

! ---       forces
	d(i) = d(i) + dv*calc%a_vec
	d(j) = d(j) + dv*calc%b_vec
	d(k) = d(k) + dv*calc%c_vec
	d(l) = d(l) + dv*calc%d_vec
end do
#ifdef _OPENMP
!$omp atomic update
improper2 = improper2 + mp_real_tmp
#endif
!$omp end parallel
end function improper2

!-----------------------------------------------------------------------

subroutine qangle (istate)
! arguments
integer						:: istate

! local variables
integer						:: ia,i,j,k,ic,im,icoupl,ib
TYPE(angl_val)					:: calc
real(kind=prec)					:: da,ae,dv


do ia = 1, nqangle


ic = qang(ia)%cod(istate)
 !skip if angle not present (code 0)
if ( ic .gt. 0 ) then

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

calc = angle_calc(x(i),x(j),x(k))

da = calc%angl - qanglib(ic)%ang0
ae = 0.5_prec*qanglib(ic)%fk*da**2
EQ(istate)%q%angle = EQ(istate)%q%angle + ae*gamma

dv = gamma*qanglib(ic)%fk*da*EQ(istate)%lambda

d(i) = d(i) + calc%a_vec*dv 
d(j) = d(j) + calc%b_vec*dv
d(k) = d(k) + calc%c_vec*dv

if ( icoupl .ne. 0 ) then

   i  = qbnd(ib)%i
   j  = qbnd(ib)%j

   d(i) = d(i) + dMorse_i(ib)*ae
   d(j) = d(j) + dMorse_j(ib)*ae

end if
end if
end do
end subroutine qangle

!-----------------------------------------------------------------------

subroutine qurey_bradley (istate)
! arguments
integer						:: istate



! local variables
integer						::	ia,i,j,k,ic,im,icoupl,ib
TYPE(bond_val)					::	calc
real(kind=prec)					::	gamma
real(kind=prec)					::	du, ru, Eurey

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
	k  = qang(ia)%k
	calc = bond_calc(x(i),x(k))
        ru = calc%dist - qanglib(ic)%ureyr0
        Eurey = qanglib(ic)%ureyfk*ru**2
        EQ(istate)%q%angle = EQ(istate)%q%angle + Eurey*gamma
        du = gamma*2*(qanglib(ic)%ureyfk*ru/calc%dist)*EQ(istate)%lambda

	d(k) = d(k) + du*calc%a_vec
	d(i) = d(i) + du*calc%b_vec

if ( icoupl .ne. 0 ) then

   i  = qbnd(ib)%i
   j  = qbnd(ib)%j

   d(i) = d(i) + dMorse_i(ib)*Eurey
   d(j) = d(j) + dMorse_j(ib)*Eurey

end if
end do
end subroutine qurey_bradley

!-----------------------------------------------------------------------

subroutine qbond (istate)
! arguments
integer						::	istate

! local variables
integer						::	ib,i,j,ic
TYPE(bond_val)					::	calc
real(kind=prec)					::	b,db,be,dv,fexp

do ib = 1, nqbond

        ic = qbnd(ib)%cod(istate)
        !code 0 means bond not present
	if ( ic > 0 ) then

	i  = qbnd(ib)%i
	j  = qbnd(ib)%j

	calc = bond_calc(x(i),x(j))
	db = calc%dist - qbondlib(ic)%r0

        fexp = exp ( -qbondlib(ic)%amz*db ) 
        be = qbondlib(ic)%Dmz*(fexp*fexp-2.0_prec*fexp) + 0.5_prec*qbondlib(ic)%fk*db**2
        EMorseD(ib) = -(fexp*fexp-2.0_prec*fexp)
	EQ(istate)%q%bond = EQ(istate)%q%bond + be
	dv = (2.0_prec*qbondlib(ic)%Dmz*qbondlib(ic)%amz*(fexp-fexp*fexp) &
		+ qbondlib(ic)%fk*db)*EQ(istate)%lambda/calc%dist

	d(i) = d(i) + dv*calc%a_vec
	d(j) = d(j) + dv*calc%b_vec

!Force scaling factor to be 1 when distance is smaller than r0
 if ( db > 0 ) then
 	EMorseD(ib) = -(fexp*fexp-2.0_prec*fexp)
 	dMorse_i(ib) = 2.0_prec*qbondlib(ic)%amz*(fexp-fexp*fexp)*EQ(istate)%lambda/(calc%dist*calc%a_vec)
 	dMorse_j(ib) = 2.0_prec*qbondlib(ic)%amz*(fexp-fexp*fexp)*EQ(istate)%lambda/(calc%dist*calc%b_vec)
 else
 	EMorseD(ib) = one
 	dMorse_i(ib) = zero
 	dMorse_j(ib) = zero
 end if

end if
end do
end subroutine qbond

!-----------------------------------------------------------------------

subroutine qimproper (istate)
! arguments
integer						::	istate

! local variables
integer						::	i,j,k,l,ip,ic
integer						::	icoupl,im,ib
TYPE(tors_val)					::	calc
real(kind=prec)					::	pe,dv,arg,gamma

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

calc = improper_calc(x(i),x(j),x(k),x(l))

! ---       energy
arg = calc%angl - qimplib(ic)%imp0
arg = arg - 2.0_prec*pi*nint(arg/(2.0_prec*pi))
dv = qimplib(ic)%fk*arg
pe  = 0.5_prec*dv*arg
EQ(istate)%q%improper = EQ(istate)%q%improper + pe*gamma
dv = dv*gamma*EQ(istate)%lambda

! ---       forces

d(i) = d(i) + dv*calc%a_vec
d(j) = d(j) + dv*calc%b_vec
d(k) = d(k) + dv*calc%c_vec
d(l) = d(l) + dv*calc%d_vec

if ( icoupl .ne. 0 ) then

   i  = qbnd(ib)%i
   j  = qbnd(ib)%j

   d(i) = d(i) + dMorse_i(ib)*pe
   d(j) = d(j) + dMorse_j(ib)*pe

end if
end if
end do
end subroutine qimproper

!-----------------------------------------------------------------------

subroutine qtorsion (istate)
! arguments
integer						::	istate

! local variables
integer						::	i,j,k,l,ip,ic
integer						::	icoupl,im,ib
TYPE(tors_val)					::	calc
real(kind=prec)					::	pe,dv,arg,gamma

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

calc = torsion_calc(x(i),x(j),x(k),x(l))
! ---       energy

arg = qtorlib(ic)%rmult*calc%angl-qtorlib(ic)%deltor
pe = qtorlib(ic)%fk*(one+cos(arg))
EQ(istate)%q%torsion = EQ(istate)%q%torsion + pe*gamma
dv = -qtorlib(ic)%rmult*qtorlib(ic)%fk*sin(arg)*gamma*EQ(istate)%lambda

! ---       forces

d(i) = d(i) + dv*calc%a_vec
d(j) = d(j) + dv*calc%b_vec
d(k) = d(k) + dv*calc%c_vec
d(l) = d(l) + dv*calc%d_vec

if ( icoupl .ne. 0 ) then

   i  = qbnd(ib)%i
   j  = qbnd(ib)%j

   d(i) = d(i) + dMorse_i(ib)*pe
   d(j) = d(j) + dMorse_j(ib)*pe

end if
end if
end do
end subroutine qtorsion

!-----------------------------------------------------------------------

real(kind=prec) function torsion(istart, iend)
!arguments
integer						::	istart, iend

! local variables
integer						::	ip,ic,i,j,k,l,ic
real(kind=prec)					::	scp,phi,dv,arg,f1
TYPE(tors_val)					::	calc
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

i  = tor(ip)%i
j  = tor(ip)%j
k  = tor(ip)%k
l  = tor(ip)%l
ic = tor(ip)%cod

calc = torsion_calc(x(i),x(j),x(k),x(l))


        t => tor(ip)
        lib => torlib(t%cod)

! ---       energy

arg = torlib(ic)%rmult*calc%angl-torlib(ic)%deltor
#ifdef _OPENMP
    mp_real_tmp = mp_real_tmp + torlib(ic)%fk*(one+cos(arg))*torlib(ic)%paths   !lib%paths is previously inverted 
#else    
    torsion = torsion + torlib(ic)%fk*(one+cos(arg))*torlib(ic)%paths   !lib%paths is previously inverted 
#endif
dv = -torlib(ic)%rmult*torlib(ic)%fk*sin(arg)*torlib(ic)%paths

! ---       forces
d(i) = d(i) + dv*calc%a_vec
d(j) = d(j) + dv*calc%b_vec
d(k) = d(k) + dv*calc%c_vec
d(l) = d(l) + dv*calc%d_vec
end do
#ifdef _OPENMP
!$omp atomic update
torsion = torsion + mp_real_tmp
#endif
!$omp end parallel

end function torsion

!-----------------------------------------------------------------------

integer function shake(xx, x)
!arguments
TYPE(qr_vec)	::	xx(:), x(:)
!       returns no. of iterations

! *** local variables
integer				:: i,j,mol,ic,nits
real(kind=prec)			:: diff,corr,scp,dist
TYPE(qr_dist5)			:: xij
TYPE(qr_vec)			:: xxij
logical				:: dead_shake =.false.

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
                                xij     = q_dist5(x(i),x(j))
                                diff    = shake_mol(mol)%bond(ic)%dist2 - xij%r2
                                if(abs(diff) < SHAKE_TOL*shake_mol(mol)%bond(ic)%dist2) then
                                        shake_mol(mol)%bond(ic)%ready = .true. ! in range
                                end if
                                xxij = qvec_sub(xx(i),xx(j))
                                scp  = q_dotprod(xij%vec,xxij)
                                corr = diff/(2.0_prec*scp*(winv(i)+winv(j)))
                                x(i) = qvec_add(x(i),corr*winv(i)*xxij)
                                x(j) = qvec_add(x(i),-corr*winv(i)*xxij)
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
                                        dist = q_dist4(xx(i),xx(j))
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


end module BONDENE




