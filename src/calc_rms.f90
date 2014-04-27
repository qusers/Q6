! (C) 2014 Uppsala Molekylmekaniska HB, Uppsala, Sweden
! calc_rms.f90
! by John Marelius
! root mean square coordinate deviation

module CALC_RMS
	use CALC_BASE
	use MASKMANIP
	implicit none

!constants
	integer, parameter			::	MAX_MASKS = 10
!module variables
	type(MASK_TYPE), private, target	::	masks(MAX_MASKS)
	integer, private			::	Nmasks = 0
	type RMS_COORD_TYPE
		real, pointer		::	x(:), x0(:)
	end type RMS_COORD_TYPE
	type(RMS_COORD_TYPE), private	::	coords(MAX_MASKS)
contains

subroutine RMS_initialize
end subroutine RMS_initialize

subroutine RMS_finalize(i)
	integer						::	i

	call mask_finalize(masks(i))
end subroutine RMS_finalize

integer function RMS_add(desc)
	!arguments
	character(*)				::	desc
	integer						::	ats
	if(Nmasks == MAX_MASKS) then
		write(*,10) MAX_MASKS
		return
	end if
10	format('Sorry, the maximum number of RMS calculations is ',i2)
	!add a new RMS mask
	Nmasks = Nmasks + 1
	call mask_initialize(masks(Nmasks))
	ats =  maskmanip_make(masks(Nmasks))
	!discard if no atoms in mask
	if(ats == 0) then
		call mask_finalize(masks(Nmasks))
		Nmasks = Nmasks - 1
		RMS_add = 0
		return
	end if

	allocate(coords(Nmasks)%x(3*ats), coords(Nmasks)%x0(3*ats))

	call RMS_make_ref(Nmasks)
	RMS_add = Nmasks
	write(desc, 20) masks(Nmasks)%included
20	format('RMS coordinate deviation for ',i6,' atoms')
 end function RMS_add


subroutine RMS_calc(i)
	!arguments
	integer, intent(in)			::	i

	!locals
	real(8)						::	r

	if(i < 1 .or. i > Nmasks) return

	call mask_get(masks(i), xin, coords(i)%x)                            
	r = sqrt(  sum((coords(i)%x-coords(i)%x0)**2)/(masks(i)%included) )  ! removed 3 in front of mask(i)   
	write(100,*) r
 	write(*,100, advance='no') r
100	format(f10.3)
end subroutine RMS_calc


subroutine RMS_make_ref(i)
	integer						::	i,at

	if(i < 1 .or. i > Nmasks) return
	call mask_get(masks(i), xtop, coords(i)%x0)

end subroutine RMS_make_ref

subroutine RMS_heading(i)
	integer						::	i

	write(*,'(a)', advance='no') 'RMSd(A)'
end subroutine RMS_heading
end module CALC_RMS
