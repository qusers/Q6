! (C) 2014 Uppsala Molekylmekaniska HB, Uppsala, Sweden
! bonded.f90
! based on md.f90
! by Johan qvist, John Marelius, Anders Kaplan, Isabella Feierberg, Martin Nervall & Martin Almlf
! calculation of all bonded interactions

module BONDENE

! used modules
!use PROFILING
use SIZES
use QMATH
use BONDED
use GLOBALS
use TOPO
use QATOM
use QALLOC
!$ use omp_lib
implicit none

contains

real(kind=prec) function angle(istart, iend,angles,ang_lib)
! *** arguments
integer						::	istart, iend
TYPE(ANG_TYPE)                                  :: angles(:)
TYPE(ANGLIB_TYPE)                               :: ang_lib(:)

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
	da = calc%angl - ang_lib(ic)%ang0
#ifdef _OPENMP
    mp_real_tmp = mp_real_tmp + 0.5_prec*ang_lib(ic)%fk*da**2
#else
    angle = angle + 0.5_prec*ang_lib(ic)%fk*da**2
#endif
	dv = ang_lib(ic)%fk*da
! calculate forces from function vectors
! and update force array
d(i) = d(i) + calc%a_vec*dv
d(j) = d(j) - calc%b_vec*dv
d(k) = d(k) + calc%c_vec*dv

end do
#ifdef _OPENMP
!$omp atomic update
angle = angle + mp_real_tmp
#endif
!$omp end parallel

end function angle

!-----------------------------------------------------------------------
real(kind=prec) function urey_bradley(istart, iend,angles,ang_lib)
! *** arguments
integer						::	istart, iend
TYPE(ANG_TYPE)                                  :: angles(:)
TYPE(ANGLIB_TYPE)                               :: ang_lib(:)
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
	if(ang_lib(ic)%ureyfk .gt. zero) then
		i  = angles(ia)%i
		k  = angles(ia)%k
		calc = bond_calc(x(i),x(k))
		ru = calc%dist - ang_lib(ic)%ureyr0
#ifdef _OPENMP
		mp_real_tmp = mp_real_tmp + ang_lib(ic)%ureyfk*ru**2
#else 
                urey_bradley = urey_bradley + ang_lib(ic)%ureyfk*ru**2
#endif
                du = 2*ang_lib(ic)%ureyfk*ru/calc%dist
		d(k) = d(k) - calc%a_vec*du
		d(i) = d(i) + calc%a_vec*du
        end if
end do
#ifdef _OPENMP
!$omp atomic update
urey_bradley = urey_bradley + mp_real_tmp
#endif
!$omp end parallel

end function urey_bradley

!-----------------------------------------------------------------------

real(kind=prec) function bond(istart, iend, bonds, bnd_lib)
! *** arguments
integer						::	istart, iend
TYPE(BOND_TYPE)                                 :: bonds(:)
TYPE(BONDLIB_TYPE)                              :: bnd_lib(:)
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
        db = calc%dist - bnd_lib(ic)%bnd0
#ifdef _OPENMP
    mp_real_tmp = mp_real_tmp + 0.5_prec*bnd_lib(ic)%fk*db**2
#else  
        bond = bond + 0.5_prec*bnd_lib(ic)%fk*db**2
#endif

        ! calculate dv and update d
        dv = bnd_lib(ic)%fk*db/calc%dist
	d(i) = d(i) - calc%a_vec*dv
	d(j) = d(j) + calc%a_vec*dv
end do
#ifdef _OPENMP
!$omp atomic update
bond = bond + mp_real_tmp
#endif
!$omp end parallel

end function bond
!-----------------------------------------------------------------------

real(kind=prec) function improper(istart, iend, impropers, imp_lib)
!arguments
integer						::	istart, iend
TYPE(TOR_TYPE)                                  :: impropers(:)
TYPE(IMPLIB_TYPE)                               :: imp_lib(:)
! evaluate harmonic impropers
! local variables
integer						::	ip,i,j,k,l,ic
TYPE(tors_val)					::	calc
real(kind=prec)					::	dv,arg,phi
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

	arg = calc%angl - imp_lib(ic)%imp0
	arg = arg - 2.0_prec*pi*nint(arg/(2.0_prec*pi))
	dv  = imp_lib(ic)%fk*arg
#ifdef _OPENMP
	    mp_real_tmp = mp_real_tmp + 0.5_prec*dv*arg
#else 
	    improper = improper + 0.5_prec*dv*arg
#endif

! ---       forces
	d(i) = d(i) + calc%a_vec*dv
	d(j) = d(j) + calc%b_vec*dv
	d(k) = d(k) + calc%c_vec*dv
	d(l) = d(l) + calc%d_vec*dv

end do
#ifdef _OPENMP
!$omp atomic update
improper = improper + mp_real_tmp
#endif
!$omp end parallel
end function improper

!-----------------------------------------------------------------------

real(kind=prec) function improper2(istart, iend, impropers, imp_lib)
!evaluate periodic impropers
!arguments
integer						::	istart, iend
TYPE(TOR_TYPE)                                  :: impropers(:)
TYPE(IMPLIB_TYPE)                               :: imp_lib(:)
! local variables
integer						::	ip,i,j,k,l,ic
real(kind=prec)					::	dv,arg,phi
TYPE(tors_val)					::	calc
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

	arg = 2*calc%angl - imp_lib(ic)%imp0
#ifdef _OPENMP
	mp_real_tmp = mp_real_tmp + imp_lib(ic)%fk * (1 + q_cos(arg))
#else
	improper2 = improper2 + imp_lib(ic)%fk * (1 + q_cos(arg))
#endif
	dv  = -2*imp_lib(ic)%fk * q_sin(arg)

! ---       forces
	d(i) = d(i) + calc%a_vec*dv
	d(j) = d(j) + calc%b_vec*dv
	d(k) = d(k) + calc%c_vec*dv
	d(l) = d(l) + calc%d_vec*dv
end do
#ifdef _OPENMP
!$omp atomic update
improper2 = improper2 + mp_real_tmp
#endif
!$omp end parallel
end function improper2

!-----------------------------------------------------------------------

subroutine qangle (E_loc,lambda,istate)
! arguments
real(kind=prec)                                 :: E_loc,lambda
integer						:: istate

! local variables
integer						:: ia,i,j,k,ic,im,icoupl,ib
TYPE(angl_val)					:: calc
real(kind=prec)					:: da,ae,dv,gamma


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
E_loc = E_loc + ae*gamma

dv = gamma*qanglib(ic)%fk*da*lambda

d(i) = d(i) + calc%a_vec*dv 
d(j) = d(j) - calc%b_vec*dv
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

subroutine qurey_bradley (E_loc,lambda,istate)
! arguments
real(kind=prec)                         :: E_loc,lambda
integer					:: istate



! local variables
integer					:: ia,i,j,k,ic,im,icoupl,ib
TYPE(bond_val)				:: calc
real(kind=prec)				:: gamma
real(kind=prec)				:: du, ru, Eurey

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
        E_loc = E_loc + Eurey*gamma
        du = gamma*2*(qanglib(ic)%ureyfk*ru/calc%dist)*lambda

	d(i) = d(i) - calc%a_vec*du
	d(k) = d(k) + calc%a_vec*du

if ( icoupl .ne. 0 ) then

   i  = qbnd(ib)%i
   j  = qbnd(ib)%j

   d(i) = d(i) + dMorse_i(ib)*Eurey
   d(j) = d(j) + dMorse_j(ib)*Eurey

end if
end do
end subroutine qurey_bradley

!-----------------------------------------------------------------------

subroutine qbond (E_loc,lambda,istate)
! arguments
real(kind=prec)                                 :: E_loc,lambda
integer						::	istate

! local variables
integer						::	ib,i,j,ic
TYPE(bond_val)					::	calc
real(kind=prec)					::	b,db,be,dv,fexp,tmp

do ib = 1, nqbond

        ic = qbnd(ib)%cod(istate)
        !code 0 means bond not present
	if ( ic > 0 ) then

	i  = qbnd(ib)%i
	j  = qbnd(ib)%j

	calc = bond_calc(x(i),x(j))
	db = calc%dist - qbondlib(ic)%r0

        fexp = q_exp ( -qbondlib(ic)%amz*db ) 
        be = qbondlib(ic)%Dmz*(fexp*fexp-2.0_prec*fexp) + 0.5_prec*qbondlib(ic)%fk*db**2
        EMorseD(ib) = -(fexp*fexp-2.0_prec*fexp)
	E_loc = E_loc + be
	dv = (2.0_prec*qbondlib(ic)%Dmz*qbondlib(ic)%amz*(fexp-fexp*fexp) &
		+ qbondlib(ic)%fk*db)*lambda/calc%dist

	d(i) = d(i) - calc%a_vec*dv
	d(j) = d(j) + calc%a_vec*dv

!Force scaling factor to be 1 when distance is smaller than r0
 if ( db > zero ) then
 	EMorseD(ib) = -(fexp*fexp-2.0_prec*fexp)
        dMorse_i(ib) = 2.0_prec*qbondlib(ic)%amz*(fexp-fexp*fexp)*EQ(istate)%lambda/(calc%a_vec*calc%dist)
 	dMorse_j(ib) = -2.0_prec*qbondlib(ic)%amz*(fexp-fexp*fexp)*EQ(istate)%lambda/(calc%a_vec*calc%dist)
 else
 	EMorseD(ib) = one
 	dMorse_i(ib) = dMorse_i(ib) * zero
 	dMorse_j(ib) = dMorse_j(ib) * zero
 end if

end if
end do
end subroutine qbond

!-----------------------------------------------------------------------

subroutine qimproper (E_loc,lambda,istate)
! arguments
real(kind=prec)                                 :: E_loc,lambda
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
E_loc = E_loc + pe*gamma
dv = dv*gamma*lambda

! ---       forces      

d(i) = d(i) + calc%a_vec*dv
d(j) = d(j) + calc%b_vec*dv
d(k) = d(k) + calc%c_vec*dv
d(l) = d(l) + calc%d_vec*dv

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

subroutine qtorsion (E_loc,lambda,istate)
! arguments
real(kind=prec)                                 :: E_loc,lambda
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
pe = qtorlib(ic)%fk*(one+q_cos(arg))
E_loc = E_loc + pe*gamma
dv = -qtorlib(ic)%rmult*qtorlib(ic)%fk*q_sin(arg)*gamma*lambda

! ---       forces

d(i) = d(i) + calc%a_vec*dv
d(j) = d(j) + calc%b_vec*dv
d(k) = d(k) + calc%c_vec*dv
d(l) = d(l) + calc%d_vec*dv

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

real(kind=prec) function torsion(istart, iend, torsions, tor_lib)
!arguments
integer						::	istart, iend
TYPE(TOR_TYPE)                                  :: torsions(:)
TYPE(TORLIB_TYPE)                               :: tor_lib(:)
! local variables
integer						::	ip,ic,i,j,k,l
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

i  = torsions(ip)%i
j  = torsions(ip)%j
k  = torsions(ip)%k
l  = torsions(ip)%l
ic = torsions(ip)%cod

calc = torsion_calc(x(i),x(j),x(k),x(l))

! ---       energy

arg = tor_lib(ic)%rmult*calc%angl-tor_lib(ic)%deltor
#ifdef _OPENMP
    mp_real_tmp = mp_real_tmp + tor_lib(ic)%fk*(one+q_cos(arg))*tor_lib(ic)%paths   !lib%paths is previously inverted 
#else    
    torsion = torsion + tor_lib(ic)%fk*(one+q_cos(arg))*tor_lib(ic)%paths   !lib%paths is previously inverted 
#endif
dv = -tor_lib(ic)%rmult*tor_lib(ic)%fk*q_sin(arg)*tor_lib(ic)%paths

! ---       forces
d(i) = d(i) + calc%a_vec*dv
d(j) = d(j) + calc%b_vec*dv
d(k) = d(k) + calc%c_vec*dv
d(l) = d(l) + calc%d_vec*dv
end do
#ifdef _OPENMP
!$omp atomic update
torsion = torsion + mp_real_tmp
#endif
!$omp end parallel

end function torsion

!-----------------------------------------------------------------------

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

!-----------------------------------------------------------------------

subroutine p_restrain(E_loc,EQ_rest,lambda)
! arguments
real(kind=prec)                         :: E_loc, EQ_rest(:),lambda(:)
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
          dr = dr + x(i) - xtop(i)
        end if
  end do

  if(n_ctr > 0) then 
    ! only if atoms were found:

        ! form average
        dr = dr / real(n_ctr,kind=prec)
        r2      = qvec_square(dr)
        erst    = 0.5_prec*fk*r2
        E_loc  = E_loc + erst

        ! apply same force to all atoms
        do i = rstseq(ir)%i, rstseq(ir)%j
        if ( heavy(i) .or. rstseq(ir)%ih .eq. 1 .and.(.not.(excl(i).and.freeze))) then
                d(i) = d(i) + dr*(fk*iaclib(iac(i))%mass/12.010_prec)
          end if
        end do
  end if
#ifdef DEBUG
call q_vecsum(d,nat_pro,'p_restrain 1')
#endif

 
else if(rstseq(ir)%to_centre == 2) then     ! Put == 1, then equal to 2
  ! restrain to mass centre
  ! reset dr & variable to put masses
  dr = zero
  totmass = zero
  
! calculate deviation from mass center
  do i = rstseq(ir)%i, rstseq(ir)%j
        if ( heavy(i) .or. rstseq(ir)%ih .eq. 1 .and.(.not.(excl(i).and.freeze))) then
          totmass = totmass + iaclib(iac(i))%mass                              ! Add masses
          dr = dr + ((x(i)-xtop(i))*iaclib(iac(i))%mass)
		end if
  end do

 if(totmass > zero) then 
    ! only if atoms were found: (i.e has a total mass)

        ! form average
        dr = dr/totmass                                  ! divide by total mass
        r2      = qvec_square(dr)
        erst    = 0.5_prec*fk*r2
        E_loc  = E_loc + erst

        ! apply same force to all atoms
        do i = rstseq(ir)%i, rstseq(ir)%j
        if ( heavy(i) .or. rstseq(ir)%ih .eq. 1 .and.(.not.(excl(i).and.freeze))) then
                d(i) = d(i) + dr*fk
          end if
        end do
  end if
#ifdef DEBUG
call q_vecsum(d,nat_pro,'p_restrain 2')
#endif

else 
  ! restrain each atom to its topology co-ordinate
  do i = rstseq(ir)%i, rstseq(ir)%j
        if ( heavy(i) .or. rstseq(ir)%ih .eq. 1 .and.(.not.(excl(i).and.freeze))) then



          dr = x(i)-xtop(i)

       !use the periodically minimal distance:
       if( use_PBC ) then
               dr = (boxlength*q_nint(dr*inv_boxl))-dr
       end if
          r2      = qvec_square(dr)

          erst    = 0.5_prec*fk*r2
          E_loc  = E_loc + erst
          d(i) = d(i) + dr*fk
        end if
  end do
#ifdef DEBUG
call q_vecsum(d,nat_pro,'p_restrain 3')
#endif

end if
end do

! extra positional restraints (Q-state dependent)
do ir = 1, nrstr_pos
istate = rstpos(ir)%ipsi
i      = rstpos(ir)%i

dr = x(i) - rstpos(ir)%x

if ( istate .ne. 0 ) then
wgt = lambda(istate)
else
wgt = one
end if

Edum = 0.5_prec * q_dotprod(rstpos(ir)%fk,(dr*dr))

d(i) = d(i) + (rstpos(ir)%fk*dr)*wgt


if ( istate .eq. 0 ) then
do k = 1, nstates
EQ_rest(k) = EQ_rest(k) + Edum
end do
if ( nstates .eq. 0 ) E_loc = E_loc + Edum
else
EQ_rest(istate) = EQ_rest(istate) + Edum
end if
end do

! atom-atom distance restraints (Q-state dependent)
do ir = 1, nrstr_dist
istate = rstdis(ir)%ipsi
i      = rstdis(ir)%i
j      = rstdis(ir)%j

! if PBC then adjust lengths according to periodicity - MA
if( use_PBC ) then
	dr = x(j) - x(i)
	distres = bond_calc(dr,boxlength*q_nint(dr*inv_boxl))
else
	distres = bond_calc(x(i),x(j))
end if

if ( istate .ne. 0 ) then
wgt = lambda(istate)
else
wgt = one
end if



if(distres%dist < rstdis(ir)%d1) then !shorter than d1
        db     = distres%dist - rstdis(ir)%d1
elseif(distres%dist > rstdis(ir)%d2) then !longer than d2
        db     = distres%dist - rstdis(ir)%d2
else
        db = 0
        cycle !skip zero force calculation
endif

Edum   = 0.5_prec*rstdis(ir)%fk*db**2
dv     = wgt*rstdis(ir)%fk*db/distres%dist

d(i) = d(i) - distres%a_vec*dv
d(j) = d(j) + distres%a_vec*dv

if ( istate .eq. 0 ) then
do k = 1, nstates
EQ_rest(k) = EQ_rest(k) + Edum
end do
if ( nstates .eq. 0 ) E_loc = E_loc + Edum
else
EQ_rest(istate) = EQ_rest(istate) + Edum
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
	anglres = box_angle_calc(x(i),x(j),x(k),boxlength,inv_boxl)
else
	anglres = angle_calc(x(i),x(j),x(k))

end if

if ( istate .ne. 0 ) then
wgt = lambda(istate)
else
wgt = one
end if

db    = anglres%angl - (rstang(ir)%ang)*deg2rad

! dv is the force to be added in module
Edum   = 0.5_prec*rstang(ir)%fk*db**2
dv     = wgt*rstang(ir)%fk*db


        ! update d
d(i) = d(i) + anglres%a_vec*dv
d(j) = d(j) - anglres%b_vec*dv
d(k) = d(k) + anglres%c_vec*dv

if ( istate .eq. 0 ) then
do k = 1, nstates
EQ_rest(k) = EQ_rest(k) + Edum
end do
if ( nstates .eq. 0 ) E_loc = E_loc + Edum
else
EQ_rest(istate) = EQ_rest(istate) + Edum
end if
end do
if( .not. use_PBC ) then
! extra half-harmonic wall restraints
do ir = 1, nrstr_wall
        fk = rstwal(ir)%fk
        do i = rstwal(ir)%i, rstwal(ir)%j
                if ( heavy(i) .or. rstwal(ir)%ih .eq. 1 ) then

			distres = bond_calc(xwcent,x(i))
                        db = distres%dist - rstwal(ir)%d

                        if(db > zero) then
                                erst =  0.5_prec * fk * db**2 - rstwal(ir)%Dmorse
                                dv = fk*db/distres%dist
                        else
                                fexp = q_exp(rstwal(ir)%aMorse*db)
                                erst = rstwal(ir)%dMorse*(fexp*fexp-2.0_prec*fexp)
                                dv=-2.0_prec*rstwal(ir)%dMorse*rstwal(ir)%aMorse*(fexp-fexp*fexp)/distres%dist
                        end if
                        E_loc = E_loc + erst
			d(i) = d(i) + distres%a_vec*dv
                end if
        end do
end do

end if

end subroutine p_restrain
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

xx(:)=x(:)-v(:)*dt
niter=shake(x, xx)	
write(*,100) 'v', niter

v(:)=(x(:)-xx(:))/dt


end subroutine initial_shaking

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
                                xij     = q_dist5(x(j),x(i))
                                diff    = shake_mol(mol)%bond(ic)%dist2 - xij%r2
                                if(q_abs(diff) < SHAKE_TOL*shake_mol(mol)%bond(ic)%dist2) then
                                        shake_mol(mol)%bond(ic)%ready = .true. ! in range
                                end if
                                xxij = xx(i) - xx(j)
                                scp  = q_dotprod(xij%vec,xxij)
                                corr = diff/(2.0_prec*scp*(winv(i)+winv(j)))
                                x(i) = x(i) + ( xxij*corr*winv(i))
                                x(j) = x(j) + (-xxij*corr*winv(j))
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
                                        write (*,100) i,j,dist,&
                                                q_sqrt(shake_mol(mol)%bond(ic)%dist2)
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




