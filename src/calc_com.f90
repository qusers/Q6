! (C) 2014 Uppsala Molekylmekaniska HB, Uppsala, Sweden
! calc_com.f90
! by Martin Almlöf
! calculates the center of mass coordinates of whatever atom mask is specified
!TODO: precision not fixed

module CALC_COM
	use CALC_BASE
	use MASKMANIP
	implicit none

!constants
	integer, parameter			::	MAX_MASKS = 10
!module variables
	type(MASK_TYPE), private, target	::	masks(MAX_MASKS)
	integer, private			::	Nmasks = 0
	type COM_COORD_TYPE
		real(kind=prec), pointer		::	x(:), y(:), z(:), mass(:)
	end type COM_COORD_TYPE
	type(COM_COORD_TYPE), private	::	coords_mass(MAX_MASKS)
	
	type COORD_TYPE
		real(kind=prec), pointer		::	xyz(:)
	end type COORD_TYPE
	type(COORD_TYPE), private	::	coords(MAX_MASKS)


	type MASS_AVE_TYPE
		real(kind=prec)		::	x,y,z
	end type MASS_AVE_TYPE
	type(MASS_AVE_TYPE), private	::	mass_ave(MAX_MASKS)
	
	real(kind=prec),private			:: tot_mass(MAX_MASKS)
	
contains

subroutine COM_initialize
end subroutine COM_initialize

subroutine COM_finalize(i)
	integer						::	i

	call mask_finalize(masks(i))
end subroutine COM_finalize

integer function COM_add(desc)
	!arguments
	character(*)				::	desc
	character(len=80)			::	line
	integer						::	readstat
	integer						::	ats,j
	if(Nmasks == MAX_MASKS) then
		write(*,10) MAX_MASKS
		return
	end if
10	format('Sorry, the maximum number of COM calculations is ',i2)



	!add a new COM mask
	Nmasks = Nmasks + 1
	call mask_initialize(masks(Nmasks))
	ats =  maskmanip_make(masks(Nmasks))
	!discard if no atoms in mask
	if(ats == 0) then
		call mask_finalize(masks(Nmasks))
		Nmasks = Nmasks - 1
		COM_add = 0
		return
	end if

	allocate(coords(Nmasks)%xyz(3*ats))
	allocate(coords_mass(Nmasks)%x(ats), coords_mass(Nmasks)%y(ats), coords_mass(Nmasks)%z(ats), coords_mass(Nmasks)%mass(ats))

	coords_mass(Nmasks)%x(:) = zero
	coords_mass(Nmasks)%y(:) = zero
	coords_mass(Nmasks)%z(:) = zero
	coords_mass(Nmasks)%mass(:) = zero
	coords(Nmasks)%xyz(:) = zero
	
	call COM_put_mass(Nmasks)
	COM_add = Nmasks
	write(desc, 20) masks(Nmasks)%included
20	format('Center of mass position ',i6,' atoms')
 end function COM_add


subroutine COM_calc(i)
	!arguments
	integer, intent(in)			::	i
        integer                 :: ats
	!locals

	if(i < 1 .or. i > Nmasks) return

        ats = masks(i)%included
        allocate(tmp_coords(ats))
        tmp_coords = translate_coords(coords(i)%xyz)
	call mask_get(masks(i), translate_coords(xin), tmp_coords)
        coords(i)%xyz = btranslate_coords(tmp_coords)
        deallocate(tmp_coords)

	!split coords into x, y, and z coords
	coords_mass(i)%x = coords(i)%xyz(1::3)
	coords_mass(i)%y = coords(i)%xyz(2::3)
	coords_mass(i)%z = coords(i)%xyz(3::3)
	

	!calculate center of mass
	mass_ave(i)%x = dot_product(coords_mass(i)%x(:),coords_mass(i)%mass)/tot_mass(i)
	mass_ave(i)%y = dot_product(coords_mass(i)%y(:),coords_mass(i)%mass)/tot_mass(i)
	mass_ave(i)%z = dot_product(coords_mass(i)%z(:),coords_mass(i)%mass)/tot_mass(i)
	
	

	
 	write(*,100, advance='no') mass_ave(i)%x,mass_ave(i)%y,mass_ave(i)%z
100	format(3f9.4)
end subroutine COM_calc


subroutine COM_put_mass(i)
	integer						::	k,j,i,at
	real(kind=prec)						::	mass

	if(i < 1 .or. i > Nmasks) return

	tot_mass(i) = zero
	!put in masses into coords_mass
	k=1
	do j = 1, nat_pro
		if (masks(i)%MASK(j)) then
			mass = iaclib(iac(j))%mass
			coords_mass(i)%mass(k) = mass
			tot_mass(i) = tot_mass(i) + mass
			k = k+1
		end if
	end do

 	write(*,168) "Total mass: ",tot_mass(i)
168	format(a,f10.3)

end subroutine COM_put_mass

subroutine COM_heading(i)
	integer						::	i

	write(*,'(a)', advance='no') '    X        Y        Z    '
end subroutine COM_heading


end module CALC_COM
