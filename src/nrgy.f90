! (C) 2014 Uppsala Molekylmekaniska HB, Uppsala, Sweden
! nrgy.f90
! by John Marelius
! energy data and energy file I/O
!TODO: remove default real statmend - in accordance with best practice

module NRGY
use SIZES

implicit none

	character(*), parameter	::	NRGY_VERSION = '5.06'
	character(*), parameter	::	NRGY_DATE    = '2014-04-21'

	type BONDED_ENERGIES
		sequence
		real(kind=prec)					::	bond, angle, torsion, improper
	end type BONDED_ENERGIES

	type NB_ENERGIES
		sequence
		real(kind=prec)					::	el, vdw
	end type NB_ENERGIES

	type RESTRAINT_ENERGIES
		sequence
		real(kind=prec)					::	total, fix, shell, protein
		real(kind=prec)					::	solvent_radial, water_pol
	end type RESTRAINT_ENERGIES

	type ENERGIES
		sequence
		real(kind=prec)					::	potential, kinetic, LRF
		type(BONDED_ENERGIES)	::	p, w, q
		type(NB_ENERGIES)		::	pp, pw, ww, qx
		type(RESTRAINT_ENERGIES)::	restraint
	end type ENERGIES

	type Q_ENERGIES
		sequence
		real(kind=prec)				::	lambda
		real(kind=prec),allocatable		::	total(:)
		type(BONDED_ENERGIES),allocatable	::	q(:)
		type(NB_ENERGIES),allocatable		::	qx(:), qq(:), qp(:), qw(:)
		real(kind=prec),allocatable		::	restraint(:)
	end type Q_ENERGIES
	type OQ_ENERGIES
                sequence
                real(kind=prec)                                 ::      lambda
                real(kind=prec)                                 ::      total
                type(BONDED_ENERGIES)   ::      q
                type(NB_ENERGIES)               ::      qx, qq, qp, qw
                real(kind=prec)                                 ::      restraint
        end type OQ_ENERGIES

	type OFFDIAG_SAVE
		sequence	
		!integer(4) avoids unaligned access & is compatible with qdyn v2
		integer(4)				::	i,j			
		real(kind=prec)					::	Hij, rkl
	end type OFFDIAG_SAVE

	type OFFDIAG_AUX
		integer(4)				::	k,l
		real(kind=prec)					::	A, mu
	end type OFFDIAG_AUX
!Header type for energy file data structure
	type Q_ENE_HEAD
		sequence
		integer(4)					::	arrays
		integer(4)					::	totresid
		integer(4),allocatable			::	types(:)
		integer(4),allocatable			:: 	numres(:)
		integer(4),allocatable			::	resid(:)
		integer(4),allocatable			::	gcnum(:)
		character(80)				:: version
	end type Q_ENE_HEAD

	interface operator(+)
		module procedure add_ene,add_classical_ene
	end interface

	interface operator(*)
		module procedure scale_ene,scale_classical_ene
	end interface
        interface operator(/)
                module procedure div_ene, div_classical_ene
        end interface
        interface operator(-)
                module procedure sub_ene,sub_classical_ene
        end interface

public :: assignment(=)

interface assignment(=)
        module procedure ene_set, ene_set_classical
end interface


contains

!----------------------------------------------------------------------

subroutine nrgy_startup

end subroutine nrgy_startup

!----------------------------------------------------------------------

subroutine put_ene(unit, e2, OFFD,arrays,nstates)
!arguments
	integer						::	unit
	type(Q_ENERGIES), intent(in)::	e2(:)
	type(OFFDIAG_SAVE), intent(in)::OFFD(:)

!local variables
	integer						::	i, bound(1), first, last,arrays,nstates

	bound = lbound(e2)
	first = bound(1)
	bound = ubound(e2)
	last = bound(1)
	
!	do i=first, last
	do i=1,nstates
		write (unit) i, e2(i)%lambda,e2(i)%total(1:arrays),e2(i)%q(1:arrays),&
			e2(i)%qx(1:arrays),e2(i)%qq(1:arrays),e2(i)%qp(1:arrays),&
			e2(i)%qw(1:arrays),e2(i)%restraint(1:arrays)
	end do

	bound = lbound(OFFD)
	first = bound(1)
	bound = ubound(OFFD)
	last = bound(1)
	write(unit) OFFD(first:last)
	
!	write (unit) (OFFD(i)%i, OFFD(i)%j, OFFD(i)%Hij, OFFD(i)%rkl, i=first, last)

end subroutine put_ene

!----------------------------------------------------------------------

integer function get_ene(unit, e2, OFFD, nstates, noffd,arrays)
!arguments
	integer					::	unit
	type(Q_ENERGIES)		:: 	e2(:)
	type(OFFDIAG_SAVE), intent(out)::	OFFD(:)
	integer, optional				:: nstates, noffd

!local variables
	integer						::	i, bound(1), first, last, dummy,arrays
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
		read (unit, end=10) dummy,e2(i)%lambda,e2(i)%total(1:arrays),e2(i)%q,&
                        e2(i)%qx(1:arrays),e2(i)%qq(1:arrays),e2(i)%qp(1:arrays),&
                        e2(i)%qw(1:arrays),e2(i)%restraint
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
   type(OQ_ENERGIES), INTENT (IN) :: e1 (:), e2 (SIZE (e1))
   type(OQ_ENERGIES) :: add_ene (SIZE (e1))
        add_ene(:)%lambda       =e1(:)%lambda
	add_ene(:)%total	=e1(:)%total		+e2(:)%total
	add_ene(:)%q%bond	=e1(:)%q%bond		+e2(:)%q%bond
	add_ene(:)%q%angle	=e1(:)%q%angle		+e2(:)%q%angle
	add_ene(:)%q%torsion	=e1(:)%q%torsion	+e2(:)%q%torsion
	add_ene(:)%q%improper	=e1(:)%q%improper	+e2(:)%q%improper
	add_ene(:)%qx%el	=e1(:)%qx%el		+e2(:)%qx%el
	add_ene(:)%qx%vdw	=e1(:)%qx%vdw		+e2(:)%qx%vdw
	add_ene(:)%qq%el	=e1(:)%qq%el		+e2(:)%qq%el
	add_ene(:)%qq%vdw	=e1(:)%qq%vdw		+e2(:)%qq%vdw
	add_ene(:)%qp%el	=e1(:)%qp%el		+e2(:)%qp%el
	add_ene(:)%qp%vdw	=e1(:)%qp%vdw		+e2(:)%qp%vdw
	add_ene(:)%qw%el	=e1(:)%qw%el		+e2(:)%qw%el
	add_ene(:)%qw%vdw	=e1(:)%qw%vdw		+e2(:)%qw%vdw
	add_ene(:)%restraint	=e1(:)%restraint	+e2(:)%restraint

end function add_ene

!----------------------------------------------------------------------

function sub_ene (e1, e2)
   type(OQ_ENERGIES), INTENT (IN) :: e1 (:), e2 (SIZE (e1))
   type(OQ_ENERGIES) :: sub_ene (SIZE (e1))

        sub_ene(:)%lambda       =e1(:)%lambda
	sub_ene(:)%total	=e1(:)%total		-e2(:)%total
	sub_ene(:)%q%bond	=e1(:)%q%bond		-e2(:)%q%bond
	sub_ene(:)%q%angle	=e1(:)%q%angle		-e2(:)%q%angle
	sub_ene(:)%q%torsion	=e1(:)%q%torsion	-e2(:)%q%torsion
	sub_ene(:)%q%improper	=e1(:)%q%improper	-e2(:)%q%improper
	sub_ene(:)%qx%el	=e1(:)%qx%el		-e2(:)%qx%el
	sub_ene(:)%qx%vdw	=e1(:)%qx%vdw		-e2(:)%qx%vdw
	sub_ene(:)%qq%el	=e1(:)%qq%el		-e2(:)%qq%el
	sub_ene(:)%qq%vdw	=e1(:)%qq%vdw		-e2(:)%qq%vdw
	sub_ene(:)%qp%el	=e1(:)%qp%el		-e2(:)%qp%el
	sub_ene(:)%qp%vdw	=e1(:)%qp%vdw		-e2(:)%qp%vdw
	sub_ene(:)%qw%el	=e1(:)%qw%el		-e2(:)%qw%el
	sub_ene(:)%qw%vdw	=e1(:)%qw%vdw		-e2(:)%qw%vdw
	sub_ene(:)%restraint	=e1(:)%restraint	-e2(:)%restraint

end function sub_ene

!----------------------------------------------------------------------

function scale_ene (e1, k)
   type(OQ_ENERGIES), INTENT (IN):: e1 (:)
   real(kind=prec), intent(in)				:: k
   type(OQ_ENERGIES)				:: scale_ene (SIZE (e1))

        scale_ene(:)%lambda     =e1(:)%lambda
	scale_ene(:)%total	=e1(:)%total*k
	scale_ene(:)%q%bond	=e1(:)%q%bond*k
	scale_ene(:)%q%angle	=e1(:)%q%angle*k
	scale_ene(:)%q%torsion	=e1(:)%q%torsion*k
	scale_ene(:)%q%improper	=e1(:)%q%improper*k
	scale_ene(:)%qx%el	=e1(:)%qx%el*k
	scale_ene(:)%qx%vdw	=e1(:)%qx%vdw*k
	scale_ene(:)%qq%el	=e1(:)%qq%el*k
	scale_ene(:)%qq%vdw	=e1(:)%qq%vdw*k
	scale_ene(:)%qp%el	=e1(:)%qp%el*k
	scale_ene(:)%qp%vdw	=e1(:)%qp%vdw*k
	scale_ene(:)%qw%el	=e1(:)%qw%el*k
	scale_ene(:)%qw%vdw	=e1(:)%qw%vdw*k
	scale_ene(:)%restraint	=e1(:)%restraint*k

end function scale_ene
!----------------------------------------------------------------------
function div_ene(e1,k)
        type(OQ_ENERGIES),INTENT(IN)    :: e1(:)
        real(kind=prec),intent(in)      :: k
        type(OQ_ENERGIES)               :: div_ene (SIZE(e1))

        div_ene(:)%lambda       =e1(:)%lambda
	div_ene(:)%total	=e1(:)%total            /k
	div_ene(:)%q%bond	=e1(:)%q%bond           /k
	div_ene(:)%q%angle	=e1(:)%q%angle          /k
	div_ene(:)%q%torsion	=e1(:)%q%torsion        /k
	div_ene(:)%q%improper	=e1(:)%q%improper       /k
	div_ene(:)%qx%el	=e1(:)%qx%el            /k
	div_ene(:)%qx%vdw	=e1(:)%qx%vdw           /k
	div_ene(:)%qq%el	=e1(:)%qq%el            /k
	div_ene(:)%qq%vdw	=e1(:)%qq%vdw           /k
	div_ene(:)%qp%el	=e1(:)%qp%el            /k
	div_ene(:)%qp%vdw	=e1(:)%qp%vdw           /k
	div_ene(:)%qw%el	=e1(:)%qw%el            /k
	div_ene(:)%qw%vdw	=e1(:)%qw%vdw           /k
	div_ene(:)%restraint	=e1(:)%restraint        /k

end function div_ene

!----------------------------------------------------------------------

subroutine ene_set (e1, k)
   type(OQ_ENERGIES), INTENT (INOUT):: e1 (:)
   real(kind=prec), intent(in)				:: k

        e1(:)%lambda       =  k 
	e1(:)%total	=  k
	e1(:)%q%bond	=  k
	e1(:)%q%angle	=  k
	e1(:)%q%torsion	=  k
	e1(:)%q%improper	=  k
	e1(:)%qx%el	=  k
	e1(:)%qx%vdw	=  k
	e1(:)%qq%el	=  k
	e1(:)%qq%vdw	=  k
	e1(:)%qp%el	=  k
	e1(:)%qp%vdw	=  k
	e1(:)%qw%el	=  k
	e1(:)%qw%vdw	=  k
	e1(:)%restraint	=  k

end subroutine ene_set

!----------------------------------------------------------------------
function scale_classical_ene(e1,k)
type(ENERGIES),INTENT(IN)       :: e1
real(kind=prec),intent(IN)      :: k
type(ENERGIES)                  :: scale_classical_ene

scale_classical_ene%potential  = e1%potential  * k
scale_classical_ene%kinetic    = e1%kinetic    * k
scale_classical_ene%LRF        = e1%LRF        * k
scale_classical_ene%p%bond     = e1%p%bond     * k
scale_classical_ene%p%angle    = e1%p%angle    * k
scale_classical_ene%p%torsion  = e1%p%torsion  * k
scale_classical_ene%p%improper = e1%p%improper * k
scale_classical_ene%w%bond     = e1%w%bond     * k
scale_classical_ene%w%angle    = e1%w%angle    * k
scale_classical_ene%w%torsion  = e1%w%torsion  * k
scale_classical_ene%w%improper = e1%w%improper * k
scale_classical_ene%q%bond     = e1%q%bond     * k
scale_classical_ene%q%angle    = e1%q%angle    * k
scale_classical_ene%q%torsion  = e1%q%torsion  * k
scale_classical_ene%q%improper = e1%q%improper * k
scale_classical_ene%pp%el      = e1%pp%el      * k
scale_classical_ene%pp%vdw     = e1%pp%vdw     * k
scale_classical_ene%pw%el      = e1%pw%el      * k
scale_classical_ene%pw%vdw     = e1%pw%vdw     * k
scale_classical_ene%ww%el      = e1%ww%el      * k
scale_classical_ene%ww%vdw     = e1%ww%vdw     * k
scale_classical_ene%qx%el      = e1%qx%el      * k
scale_classical_ene%qx%vdw     = e1%qx%vdw     * k

scale_classical_ene%restraint%total          = e1%restraint%total          * k
scale_classical_ene%restraint%fix            = e1%restraint%fix            * k
scale_classical_ene%restraint%shell          = e1%restraint%shell          * k
scale_classical_ene%restraint%protein        = e1%restraint%protein        * k
scale_classical_ene%restraint%solvent_radial = e1%restraint%solvent_radial * k
scale_classical_ene%restraint%water_pol      = e1%restraint%water_pol      * k


end function scale_classical_ene
!----------------------------------------------------------------------
function add_classical_ene(e1,e2)
type(ENERGIES),INTENT(IN)       :: e1,e2
type(ENERGIES)                  :: add_classical_ene

add_classical_ene%potential      = e1%potential       + e2%potential      
add_classical_ene%kinetic        = e1%kinetic         + e2%kinetic        
add_classical_ene%LRF            = e1%LRF             + e2%LRF            
add_classical_ene%p%bond         = e1%p%bond          + e2%p%bond         
add_classical_ene%p%angle        = e1%p%angle         + e2%p%angle        
add_classical_ene%p%torsion      = e1%p%torsion       + e2%p%torsion      
add_classical_ene%p%improper     = e1%p%improper      + e2%p%improper     
add_classical_ene%w%bond         = e1%w%bond          + e2%w%bond         
add_classical_ene%w%angle        = e1%w%angle         + e2%w%angle        
add_classical_ene%w%torsion      = e1%w%torsion       + e2%w%torsion      
add_classical_ene%w%improper     = e1%w%improper      + e2%w%improper     
add_classical_ene%q%bond         = e1%q%bond          + e2%q%bond         
add_classical_ene%q%angle        = e1%q%angle         + e2%q%angle        
add_classical_ene%q%torsion      = e1%q%torsion       + e2%q%torsion      
add_classical_ene%q%improper     = e1%q%improper      + e2%q%improper     
add_classical_ene%pp%el          = e1%pp%el           + e2%pp%el          
add_classical_ene%pp%vdw         = e1%pp%vdw          + e2%pp%vdw         
add_classical_ene%pw%el          = e1%pw%el           + e2%pw%el          
add_classical_ene%pw%vdw         = e1%pw%vdw          + e2%pw%vdw         
add_classical_ene%ww%el          = e1%ww%el           + e2%ww%el          
add_classical_ene%ww%vdw         = e1%ww%vdw          + e2%ww%vdw         
add_classical_ene%qx%el          = e1%qx%el           + e2%qx%el          
add_classical_ene%qx%vdw         = e1%qx%vdw          + e2%qx%vdw         
add_classical_ene%restraint%total          = e1%restraint%total           + e2%restraint%total          
add_classical_ene%restraint%fix            = e1%restraint%fix             + e2%restraint%fix            
add_classical_ene%restraint%shell          = e1%restraint%shell           + e2%restraint%shell          
add_classical_ene%restraint%protein        = e1%restraint%protein         + e2%restraint%protein        
add_classical_ene%restraint%solvent_radial = e1%restraint%solvent_radial  + e2%restraint%solvent_radial 
add_classical_ene%restraint%water_pol      = e1%restraint%water_pol       + e2%restraint%water_pol      

end function add_classical_ene
!----------------------------------------------------------------------
function div_classical_ene(e1,k)
TYPE(ENERGIES),INTENT(IN)       :: e1
real(kind=prec),INTENT(IN)      :: k
TYPE(ENERGIES)                  :: div_classical_ene

div_classical_ene%potential      = e1%potential        / k 
div_classical_ene%kinetic        = e1%kinetic          / k 
div_classical_ene%LRF            = e1%LRF              / k 
div_classical_ene%p%bond         = e1%p%bond           / k 
div_classical_ene%p%angle        = e1%p%angle          / k 
div_classical_ene%p%torsion      = e1%p%torsion        / k 
div_classical_ene%p%improper     = e1%p%improper       / k 
div_classical_ene%w%bond         = e1%w%bond           / k 
div_classical_ene%w%angle        = e1%w%angle          / k 
div_classical_ene%w%torsion      = e1%w%torsion        / k 
div_classical_ene%w%improper     = e1%w%improper       / k 
div_classical_ene%q%bond         = e1%q%bond           / k 
div_classical_ene%q%angle        = e1%q%angle          / k 
div_classical_ene%q%torsion      = e1%q%torsion        / k 
div_classical_ene%q%improper     = e1%q%improper       / k 
div_classical_ene%pp%el          = e1%pp%el            / k 
div_classical_ene%pp%vdw         = e1%pp%vdw           / k 
div_classical_ene%pw%el          = e1%pw%el            / k 
div_classical_ene%pw%vdw         = e1%pw%vdw           / k 
div_classical_ene%ww%el          = e1%ww%el            / k 
div_classical_ene%ww%vdw         = e1%ww%vdw           / k 
div_classical_ene%qx%el          = e1%qx%el            / k 
div_classical_ene%qx%vdw         = e1%qx%vdw           / k 
div_classical_ene%restraint%total          = e1%restraint%total            / k 
div_classical_ene%restraint%fix            = e1%restraint%fix              / k 
div_classical_ene%restraint%shell          = e1%restraint%shell            / k 
div_classical_ene%restraint%protein        = e1%restraint%protein          / k 
div_classical_ene%restraint%solvent_radial = e1%restraint%solvent_radial   / k 
div_classical_ene%restraint%water_pol      = e1%restraint%water_pol        / k 

end function div_classical_ene

!----------------------------------------------------------------------

function sub_classical_ene(e1,e2)
type(ENERGIES),INTENT(IN)       :: e1,e2
type(ENERGIES)                  :: sub_classical_ene

sub_classical_ene%potential      = e1%potential       - e2%potential      
sub_classical_ene%kinetic        = e1%kinetic         - e2%kinetic        
sub_classical_ene%LRF            = e1%LRF             - e2%LRF            
sub_classical_ene%p%bond         = e1%p%bond          - e2%p%bond         
sub_classical_ene%p%angle        = e1%p%angle         - e2%p%angle        
sub_classical_ene%p%torsion      = e1%p%torsion       - e2%p%torsion      
sub_classical_ene%p%improper     = e1%p%improper      - e2%p%improper     
sub_classical_ene%w%bond         = e1%w%bond          - e2%w%bond         
sub_classical_ene%w%angle        = e1%w%angle         - e2%w%angle        
sub_classical_ene%w%torsion      = e1%w%torsion       - e2%w%torsion      
sub_classical_ene%w%improper     = e1%w%improper      - e2%w%improper     
sub_classical_ene%q%bond         = e1%q%bond          - e2%q%bond         
sub_classical_ene%q%angle        = e1%q%angle         - e2%q%angle        
sub_classical_ene%q%torsion      = e1%q%torsion       - e2%q%torsion      
sub_classical_ene%q%improper     = e1%q%improper      - e2%q%improper     
sub_classical_ene%pp%el          = e1%pp%el           - e2%pp%el          
sub_classical_ene%pp%vdw         = e1%pp%vdw          - e2%pp%vdw         
sub_classical_ene%pw%el          = e1%pw%el           - e2%pw%el          
sub_classical_ene%pw%vdw         = e1%pw%vdw          - e2%pw%vdw         
sub_classical_ene%ww%el          = e1%ww%el           - e2%ww%el          
sub_classical_ene%ww%vdw         = e1%ww%vdw          - e2%ww%vdw         
sub_classical_ene%qx%el          = e1%qx%el           - e2%qx%el          
sub_classical_ene%qx%vdw         = e1%qx%vdw          - e2%qx%vdw         
sub_classical_ene%restraint%total          = e1%restraint%total           - e2%restraint%total          
sub_classical_ene%restraint%fix            = e1%restraint%fix             - e2%restraint%fix            
sub_classical_ene%restraint%shell          = e1%restraint%shell           - e2%restraint%shell          
sub_classical_ene%restraint%protein        = e1%restraint%protein         - e2%restraint%protein        
sub_classical_ene%restraint%solvent_radial = e1%restraint%solvent_radial  - e2%restraint%solvent_radial 
sub_classical_ene%restraint%water_pol      = e1%restraint%water_pol       - e2%restraint%water_pol      

end function sub_classical_ene
!-------------------------------------------------------------------------
subroutine ene_set_classical(e1,k)
type(ENERGIES),INTENT(INOUT)    :: e1
real(kind=prec),intent(IN)      :: k

e1%potential  = k
e1%kinetic    = k
e1%LRF        = k
e1%p%bond     = k
e1%p%angle    = k
e1%p%torsion  = k
e1%p%improper = k
e1%w%bond     = k
e1%w%angle    = k
e1%w%torsion  = k
e1%w%improper = k
e1%q%bond     = k
e1%q%angle    = k
e1%q%torsion  = k
e1%q%improper = k
e1%pp%el      = k
e1%pp%vdw     = k
e1%pw%el      = k
e1%pw%vdw     = k
e1%ww%el      = k
e1%ww%vdw     = k
e1%qx%el      = k
e1%qx%vdw     = k

e1%restraint%total          = k
e1%restraint%fix            = k
e1%restraint%shell          = k
e1%restraint%protein        = k
e1%restraint%solvent_radial = k
e1%restraint%water_pol      = k


end subroutine ene_set_classical
!----------------------------------------------------------------------












!----------------------------------------------------------------------

real(kind=prec) function sum_bonded(r, eb)
	real(kind=prec), intent(in)			::	r
	type(BONDED_ENERGIES), intent(in)::	eb

	sum_bonded = r + eb%bond + eb%angle + eb%torsion + eb%improper
end function sum_bonded

!----------------------------------------------------------------------

real(kind=prec) function sum_non_bonded(r, enb)
	real(kind=prec), intent(in)			::	r
	type(NB_ENERGIES), intent(in)::	enb

	sum_non_bonded = r + enb%el + enb%vdw
end function sum_non_bonded

!----------------------------------------------------------------------

end module NRGY 
