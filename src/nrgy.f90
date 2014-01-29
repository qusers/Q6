!	(C) 2000 Uppsala Molekylmekaniska HB, Uppsala, Sweden

!	nrgy.f90
!	by John Marelius

!	energy data and energy file I/O

module NRGY
use SIZES

implicit none

	character(*), parameter	::	NRGY_VERSION = '5.01'
	character(*), parameter	::	NRGY_DATE = '2003-06-02'

	type BONDED_ENERGIES
		sequence
		real(8)					::	bond, angle, torsion, improper
	end type BONDED_ENERGIES

	type NB_ENERGIES
		sequence
		real(8)					::	el, vdw
	end type NB_ENERGIES

	type RESTRAINT_ENERGIES
		sequence
		real(8)					::	total, fix, shell, protein
		real(8)					::	solvent_radial, water_pol
	end type RESTRAINT_ENERGIES

	type ENERGIES
		sequence
		real(8)					::	potential, kinetic, LRF
		type(BONDED_ENERGIES)	::	p, w, q
		type(NB_ENERGIES)		::	pp, pw, ww, qx
		type(RESTRAINT_ENERGIES)::	restraint
	end type ENERGIES

	type Q_ENERGIES
		sequence
		real(8)					::	lambda
		real(8)					::	total
		type(BONDED_ENERGIES)	::	q
		type(NB_ENERGIES)		::	qx, qq, qp, qw
		real(8)					::	restraint
	end type Q_ENERGIES

	type OFFDIAG_SAVE
		sequence	
		!integer(4) avoids unaligned access & is compatible with qdyn v2
		integer(4)				::	i,j			
		real(8)					::	Hij, rkl
	end type OFFDIAG_SAVE

	type OFFDIAG_AUX
		integer(4)				::	k,l
		real(8)					::	A, mu
	end type OFFDIAG_AUX

	interface operator(+)
		module procedure add_ene
	end interface

	interface operator(*)
		module procedure scale_ene
	end interface
contains

!----------------------------------------------------------------------

subroutine nrgy_startup

end subroutine nrgy_startup

!----------------------------------------------------------------------

subroutine put_ene(unit, e2, OFFD)
!arguments
	integer						::	unit
	type(Q_ENERGIES), intent(in)::	e2(:)
	type(OFFDIAG_SAVE), intent(in)::OFFD(:)

!local variables
	integer						::	i, bound(1), first, last

	bound = lbound(e2)
	first = bound(1)
	bound = ubound(e2)
	last = bound(1)
	
	do i=first, last
		write (unit) i, e2(i)
	end do

	bound = lbound(OFFD)
	first = bound(1)
	bound = ubound(OFFD)
	last = bound(1)
	write(unit) OFFD(first:last)
	
!	write (unit) (OFFD(i)%i, OFFD(i)%j, OFFD(i)%Hij, OFFD(i)%rkl, i=first, last)

end subroutine put_ene

!----------------------------------------------------------------------

integer function get_ene(unit, e2, OFFD, nstates, noffd)
!arguments
	integer					::	unit
	type(Q_ENERGIES), intent(out)::	e2(:)
	type(OFFDIAG_SAVE), intent(out)::	OFFD(:)
	integer, optional				:: nstates, noffd

!local variables
	integer						::	i, bound(1), first, last, dummy
	if(present(nstates)) then
		first = 1
		last = nstates
	else
		bound = lbound(e2)
		first = bound(1)
		bound = ubound(e2)
		last = bound(1)
	end if
		
	do i=first, last
		read (unit, end=10) dummy, e2(i)
	end do
	
	if(present(noffd)) then	
		first = 1
		last = noffd
	else
		bound = lbound(OFFD)
		first = bound(1)
		bound = ubound(OFFD)
		last = bound(1)
	end if
!	read(unit, end=20) (OFFD(i)%i, OFFD(i)%j, OFFD(i)%Hij, OFFD(i)%rkl, i=first, last)
	read(unit, end=20) OFFD(first:last)

	get_ene = 0 !it's OK
	return

10	get_ene = i !failed energies
	return	
20	get_ene = -1 !failed offd
	return	

end function get_ene

!----------------------------------------------------------------------

function add_ene (e1, e2)
   type(Q_ENERGIES), INTENT (IN) :: e1 (:), e2 (SIZE (e1))
   type(Q_ENERGIES) :: add_ene (SIZE (e1))

			add_ene(:)%total=e1(:)%total+e2(:)%total
			add_ene(:)%q%bond=e1(:)%q%bond+e2(:)%q%bond
			add_ene(:)%q%angle=e1(:)%q%angle+e2(:)%q%angle
			add_ene(:)%q%torsion=e1(:)%q%torsion+e2(:)%q%torsion
			add_ene(:)%q%improper=e1(:)%q%improper+e2(:)%q%improper
			add_ene(:)%qx%el=e1(:)%qx%el+e2(:)%qx%el
			add_ene(:)%qx%vdw=e1(:)%qx%vdw+e2(:)%qx%vdw
			add_ene(:)%qq%el=e1(:)%qq%el+e2(:)%qq%el
			add_ene(:)%qq%vdw=e1(:)%qq%vdw+e2(:)%qq%vdw
			add_ene(:)%qp%el=e1(:)%qp%el+e2(:)%qp%el
			add_ene(:)%qp%vdw=e1(:)%qp%vdw+e2(:)%qp%vdw
			add_ene(:)%qw%el=e1(:)%qw%el+e2(:)%qw%el
			add_ene(:)%qw%vdw=e1(:)%qw%vdw+e2(:)%qw%vdw
			add_ene(:)%restraint=e1(:)%restraint+e2(:)%restraint

end function add_ene

!----------------------------------------------------------------------

function scale_ene (e1, k)
   type(Q_ENERGIES), INTENT (IN):: e1 (:)
   real, intent(in)				:: k
   type(Q_ENERGIES)				:: scale_ene (SIZE (e1))

			scale_ene(:)%total=e1(:)%total*k
			scale_ene(:)%q%bond=e1(:)%q%bond*k
			scale_ene(:)%q%angle=e1(:)%q%angle*k
			scale_ene(:)%q%torsion=e1(:)%q%torsion*k
			scale_ene(:)%q%improper=e1(:)%q%improper*k
			scale_ene(:)%qx%el=e1(:)%qx%el*k
			scale_ene(:)%qx%vdw=e1(:)%qx%vdw*k
			scale_ene(:)%qq%el=e1(:)%qq%el*k
			scale_ene(:)%qq%vdw=e1(:)%qq%vdw*k
			scale_ene(:)%qp%el=e1(:)%qp%el*k
			scale_ene(:)%qp%vdw=e1(:)%qp%vdw*k
			scale_ene(:)%qw%el=e1(:)%qw%el*k
			scale_ene(:)%qw%vdw=e1(:)%qw%vdw*k
			scale_ene(:)%restraint=e1(:)%restraint*k

end function scale_ene

!----------------------------------------------------------------------

real(8) function sum_bonded(r, eb)
	real(8), intent(in)			::	r
	type(BONDED_ENERGIES), intent(in)::	eb

	sum_bonded = r + eb%bond + eb%angle + eb%torsion + eb%improper
end function sum_bonded

!----------------------------------------------------------------------

real(8) function sum_non_bonded(r, enb)
	real(8), intent(in)			::	r
	type(NB_ENERGIES), intent(in)::	enb

	sum_non_bonded = r + enb%el + enb%vdw
end function sum_non_bonded

!----------------------------------------------------------------------

end module NRGY 
