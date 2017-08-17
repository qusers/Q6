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
! calc_fit.f90
! by John Marelius
! parts (C) 1979 W. F. van Gunsteren
! pair-wise superposition of coordidinate sets

module CALC_FIT
	use CALC_BASE
	use MASKMANIP
        use QEIGEN
	implicit none

!constants
	integer, parameter			::	MAX_MASKS = 10
!module variables
	type(MASK_TYPE), private, target	::	masks(MAX_MASKS)
	integer, private			::	Nmasks = 0

	type FIT_COORD_TYPE
		real(kind=prec), pointer		::	xr(:)
		real(kind=prec)					::	xrcm(3)
	end type FIT_COORD_TYPE
	type(FIT_COORD_TYPE), private	::	coords(MAX_MASKS)

contains
	
subroutine fit_initialize
        calc_xtop(1:3*nat_pro-2:3) = xtop(1:nat_pro)%x
        calc_xtop(2:3*nat_pro-1:3) = xtop(1:nat_pro)%y
        calc_xtop(3:3*nat_pro  :3) = xtop(1:nat_pro)%z
end subroutine fit_initialize

subroutine fit_finalize(i)
	integer						::	i

	call mask_finalize(masks(i))
end subroutine fit_finalize

integer function fit_add(desc)
	!arguments
	character(*)				::	desc
	integer						::	ats
	if(Nmasks == MAX_MASKS) then
		write(*,10) MAX_MASKS
		return
	end if
10	format('Sorry, the maximum number of LSQ fits is ',i2)

	!add a new RMS mask
	Nmasks = Nmasks + 1
	call mask_initialize(masks(Nmasks))
	ats =  maskmanip_make(masks(Nmasks))

	!discard if no atoms in mask
	if(ats == 0) then
		call mask_finalize(masks(Nmasks))
		Nmasks = Nmasks - 1
		fit_add = 0
		return
	end if

	allocate(coords(Nmasks)%xr(3*nat_pro))

	call fit_make_ref(Nmasks)

	fit_add = Nmasks
	write(desc, 20) masks(Nmasks)%included
20	format('Least squares fit of ',i6,' atoms')

end function fit_add

subroutine fit_calc(i)
	!arguments
	integer, intent(in)			::	i

	!locals
	integer						::	at
	integer						::	lsqstatus
	real(kind=prec)						::	xcm(3), totmass, U(3,3), E

    if(i < 1 .or. i > Nmasks) return
	
	!calc. mass centre of xin
	xcm(:) = zero
	totmass = zero
	do at = 1, nat_pro
		if(masks(i)%mask(at)) then
	      xcm(:) = xcm(:) + xin(3*at-2:3*at)*iaclib(iac(at))%mass 
		  totmass = totmass + iaclib(iac(at))%mass 
		end if
	end do
	
  	! Calculate center of mass (/masks(i)%included for geometric)
	  xcm(:) = xcm(:)/totmass   
    
    ! shift xin to origin and massweight
	do at = 1,nat_pro
	  xin(3*at-2) = (xin(3*at-2) - xcm(1))*sqrt(iaclib(iac(at))%mass)  !sqrt(iaclib(iac(at))%mass)
	  xin(3*at-1) = (xin(3*at-1) - xcm(2))*sqrt(iaclib(iac(at))%mass)  !sqrt(iaclib(iac(at))%mass)
	  xin(3*at  ) = (xin(3*at  ) - xcm(3))*sqrt(iaclib(iac(at))%mass)  !sqrt(iaclib(iac(at))%mass)
	end do
	
    !rotate xin to fit with coords(i)%xr
	lsqstatus = LSQSTR(nat_pro, masks(i)%mask(:), coords(i)%xr(:), xin(:),E, U)

   	! Remove massweighting and 
	! Add the translation vector
    if( lsqstatus == 1 ) then
      print*, "Failed!"
	endif

	do at = 1,nat_pro
	  xin(3*at-2) = xin(3*at-2)/sqrt(iaclib(iac(at))%mass) + coords(i)%xrcm(1) 
	  xin(3*at-1) = xin(3*at-1)/sqrt(iaclib(iac(at))%mass) + coords(i)%xrcm(2) 
	  xin(3*at  ) = xin(3*at  )/sqrt(iaclib(iac(at))%mass) + coords(i)%xrcm(3)  
	enddo

 

end subroutine fit_calc

subroutine fit_make_ref(i)
	!arguments
	integer						::	i
	!locals
	integer						::	at
    real(kind=prec)                     :: totmass

	if(i < 1 .or. i > Nmasks) return
	!calc centre vector
	coords(i)%xrcm(:) = zero
	totmass = zero
	do at = 1, nat_pro
		if(masks(i)%mask(at)) then
		   coords(i)%xrcm(:) = coords(i)%xrcm(:) + calc_xtop(3*at-2:3*at)*iaclib(iac(at))%mass 
		   totmass = totmass + iaclib(iac(at))%mass 
		end if
	end do
	
	coords(i)%xrcm(:) = coords(i)%xrcm(:)/totmass !/for geometric masks(i)%included

	! shift to origin and massweight coordinates
	do at = 1,nat_pro
	 	coords(i)%xr(3*at-2) = (calc_xtop(3*at-2) - coords(i)%xrcm(1))*sqrt(iaclib(iac(at))%mass)
		coords(i)%xr(3*at-1) = (calc_xtop(3*at-1) - coords(i)%xrcm(2))*sqrt(iaclib(iac(at))%mass)
		coords(i)%xr(3*at  ) = (calc_xtop(3*at  ) - coords(i)%xrcm(3))*sqrt(iaclib(iac(at))%mass)
	end do

end subroutine fit_make_ref

end module CALC_FIT
