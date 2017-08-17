! Q6: A comprehensive simulation package for molecular dynamics simulations and 
! free energy calculations, including empirical valence bond simulations, 
! linear interaction energy calculations, and free energy perturbation.
! 
! Copyright © 2017 Johan Åqvist, John Marelius, Shina Caroline Lynn Kamerlin and Paul Bauer
! 
! This program is free software; you can redistribute it and/or modify it under the 
! terms of the GNU General Public License as published by the Free 
! Software Foundation; either version 2 of the License, or any later version.
! 
! This program is distributed in the hope that it will be useful, 
! but WITHOUT ANY WARRANTY; without even the implied warranty of 
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
! See the GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License along with 
! this program; if not, write to the Free Software Foundation, Inc., 51 Franklin 
! Street, Fifth Floor, Boston, MA  02110-1301, USA. Also add information on 
! how to contact you by electronic and paper mail.
! calc_rmsf.f90
! by Martin Almlöf
! root mean square coordinate deviation

module CALC_RMSF
	use CALC_BASE
	use MASKMANIP
	implicit none

!constants
	integer, parameter			::	MAX_MASKS = 10
!module variables
	integer, private			:: frames(MAX_MASKS), apa
	real(kind=prec), allocatable		:: rmsf(:)
	type(MASK_TYPE), private, target	::	masks(MAX_MASKS)
	integer, private			::	Nmasks = 0
	type RMSF_COORD_TYPE
		real(kind=prec), pointer		::	x(:), x0(:), x2(:)
	end type RMSF_COORD_TYPE
	type(RMSF_COORD_TYPE), private	::	coords(MAX_MASKS)
contains

subroutine RMSF_initialize
end subroutine RMSF_initialize

subroutine RMSF_finalize(i)
	integer						::	i

	call mask_finalize(masks(i))
end subroutine RMSF_finalize

integer function RMSF_add(desc)
	!arguments
	character(*)				::	desc
	integer						::	ats
	if(Nmasks == MAX_MASKS) then
		write(*,10) MAX_MASKS
		return
	end if
10	format('Sorry, the maximum number of RMSF calculations is ',i2)
	!add a new RMSF mask
	Nmasks = Nmasks + 1
	call mask_initialize(masks(Nmasks))
	ats =  maskmanip_make(masks(Nmasks))
	!discard if no atoms in mask
	if(ats == 0) then
		call mask_finalize(masks(Nmasks))
		Nmasks = Nmasks - 1
		RMSF_add = 0
		return
	end if

	allocate(coords(Nmasks)%x(3*ats), coords(Nmasks)%x0(3*ats),coords(Nmasks)%x2(3*ats), rmsf(masks(Nmasks)%included))
	coords(Nmasks)%x2(:) = 0
	coords(Nmasks)%x0(:) = 0
	frames(Nmasks) = 0
	call RMSF_make_ref(Nmasks)
	RMSF_add = Nmasks
	write(desc, 20) masks(Nmasks)%included
20	format('RMSF coordinate deviation for ',i6,' atoms')
 end function RMSF_add


subroutine RMSF_calc(i)
	!arguments
	integer, intent(in)			::	i

	!locals
	real(kind=prec)						::	r
        integer :: ats

	if(i < 1 .or. i > Nmasks) return

	frames(i)=frames(i) + 1
        ats = masks(i)%included
        allocate(tmp_coords(ats))
        tmp_coords = translate_coords(coords(i)%x)
	call mask_get(masks(i), translate_coords(xin), tmp_coords)
        coords(i)%x = btranslate_coords(tmp_coords)
        deallocate(tmp_coords)
	!calculate the sum of the squared coordinates up to this frame
	coords(i)%x2 = coords(i)%x2 + (coords(i)%x)**2

	!calculate the average coordinates up to this frame
	coords(i)%x0 = coords(i)%x0*(1._8 - 1._8/frames(i)) + coords(i)%x/frames(i)

	do apa = 1, masks(i)%included
		rmsf(apa) = sqrt(1._8/frames(i) * sum(coords(i)%x2(3*apa-2:3*apa)) - sum((coords(i)%x0(3*apa-2:3*apa))**2))
	end do
	
	
	
	r = 1._8/masks(i)%included * sum(rmsf)

	
 	write(*,100, advance='no') r
100	format(f10.3)
end subroutine RMSF_calc


subroutine RMSF_make_ref(i)
	integer						::	i,at,ats

	if(i < 1 .or. i > Nmasks) return
        ats = masks(i)%included
        allocate(tmp_coords(ats))
        tmp_coords = translate_coords(coords(i)%x0)
	call mask_get(masks(i), xtop, tmp_coords)
        coords(i)%x0 = btranslate_coords(tmp_coords)
        deallocate(tmp_coords)

end subroutine RMSF_make_ref

subroutine RMSF_heading(i)
	integer						::	i

	write(*,'(a)', advance='no') 'RMSFd(A)'
end subroutine RMSF_heading
end module CALC_RMSF
