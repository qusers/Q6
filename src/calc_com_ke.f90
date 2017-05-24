! (C) 2014 Uppsala Molekylmekaniska HB, Uppsala, Sweden
! calc_kineticenergy.f90
! by Martin Almlöf
! calculates the kinetic energy of the center of mass of whatever atom mask is specified
!TODO: precision not fixed

module CALC_COM_KE
	use CALC_BASE
	use MASKMANIP
        use QEIGEN
!	use LAPACKMINI
	implicit none

!constants
	integer, parameter			::	MAX_MASKS = 10
!module variables
	real(kind=prec), parameter, private :: conversion_factor = 2390.0574   ! gram/mol*Å^2/fs^2  -->  kcal/mol
	integer, private			:: frames(MAX_MASKS), apa
	real(kind=prec), allocatable		:: kineticenergy(:)
	type(MASK_TYPE), private, target	::	masks(MAX_MASKS)
	integer, private			::	Nmasks = 0
	type COM_KE_COORD_TYPE
		real(kind=prec), pointer		::	x(:), y(:), z(:), mass(:)
	end type COM_KE_COORD_TYPE
	type(COM_KE_COORD_TYPE), private	::	coords_mass(MAX_MASKS), prev_coords_mass(MAX_MASKS)

	type COM_KE_VELOCITY_TYPE
		real(kind=prec), pointer		::	x(:), y(:), z(:)
	end type COM_KE_VELOCITY_TYPE
	type(COM_KE_VELOCITY_TYPE), private	::	velocity(MAX_MASKS), rel_coords(MAX_MASKS), prev_rel_coords(MAX_MASKS),rad_vec(MAX_MASKS,3) !rel_coords is not velocities
	
	type COORD_TYPE
		real(kind=prec), pointer		::	xyz(:)
	end type COORD_TYPE
	type(COORD_TYPE), private	::	coords(MAX_MASKS), prev_coords(MAX_MASKS)

	type DP_TYPE
		real(kind=prec), pointer		::	dp(:)
	end type DP_TYPE
	type(DP_TYPE), private	::	dp_vect(MAX_MASKS)



	type MASS_AVE_TYPE
		real(kind=prec)		::	x,y,z
	end type MASS_AVE_TYPE
	type(MASS_AVE_TYPE), private	::	mass_ave(MAX_MASKS), prev_mass_ave(MAX_MASKS) , ang_momentum(MAX_MASKS,3)
	
	type EIGEN_STUFF_TYPE
		real(kind=prec)		::	evalue(3),evector(3,3)
	end type EIGEN_STUFF_TYPE
	type(EIGEN_STUFF_TYPE), private	::	eigen_stuff(MAX_MASKS)




	real(kind=prec),private			:: tot_mass(MAX_MASKS), KE_rot(MAX_MASKS,3)
	logical,private			:: first_frame(MAX_MASKS) = .true.
!	real(kind=prec),private			:: previous_mass_center(3,MAX_MASKS)
	real(kind=prec), private			:: frame_length = 0
	
contains

subroutine COM_KE_initialize
end subroutine COM_KE_initialize

subroutine COM_KE_finalize(i)
	integer						::	i

	call mask_finalize(masks(i))
end subroutine COM_KE_finalize

integer function COM_KE_add(desc)
	!arguments
	character(*)				::	desc
	character(len=80)			::	line
	integer						::	readstat
	integer						::	ats,j
	if(Nmasks == MAX_MASKS) then
		write(*,10) MAX_MASKS
		return
	end if
10	format('Sorry, the maximum number of COM_KE calculations is ',i2)


	!get frame time length
	if (frame_length == 0) then
		write(*,'(a)', advance='no') 'Enter time span between each trajectory frame (fs): '
		read(*,'(a)', iostat=readstat) line
		read(line, *, iostat=readstat) frame_length
		if(readstat /= 0 .or. frame_length <= 0) then
			write(*, 900) 
	900		format('>>>>> ERROR: Invalid time span.')
			COM_KE_add = 0
			return
		end if
	end if

	!add a new COM_KE mask
	Nmasks = Nmasks + 1
	call mask_initialize(masks(Nmasks))
	ats =  maskmanip_make(masks(Nmasks))
	!discard if no atoms in mask
	if(ats == 0) then
		call mask_finalize(masks(Nmasks))
		Nmasks = Nmasks - 1
		COM_KE_add = 0
		return
	end if

	allocate(coords(Nmasks)%xyz(3*ats), prev_coords(Nmasks)%xyz(3*ats))
	allocate(coords_mass(Nmasks)%x(ats), coords_mass(Nmasks)%y(ats), coords_mass(Nmasks)%z(ats), coords_mass(Nmasks)%mass(ats))
	allocate(prev_coords_mass(Nmasks)%x(ats), prev_coords_mass(Nmasks)%y(ats))
	allocate(prev_coords_mass(Nmasks)%z(ats), prev_coords_mass(Nmasks)%mass(ats))
	allocate(velocity(Nmasks)%x(ats), velocity(Nmasks)%y(ats), velocity(Nmasks)%z(ats))
	allocate(rel_coords(Nmasks)%x(ats), rel_coords(Nmasks)%y(ats), rel_coords(Nmasks)%z(ats))
	allocate(dp_vect(Nmasks)%dp(ats) ,prev_rel_coords(Nmasks)%x(ats))
	allocate(prev_rel_coords(Nmasks)%y(ats), prev_rel_coords(Nmasks)%z(ats))

	do j=1,3
		allocate(rad_vec(Nmasks,j)%x(ats),rad_vec(Nmasks,j)%y(ats),rad_vec(Nmasks,j)%z(ats))
	end do



	coords_mass(Nmasks)%x(:) = zero
	coords_mass(Nmasks)%y(:) = zero
	coords_mass(Nmasks)%z(:) = zero
	coords_mass(Nmasks)%mass(:) = zero
	coords(Nmasks)%xyz(:) = zero
	
	frames(Nmasks) = 0
	call COM_KE_put_mass(Nmasks)
	COM_KE_add = Nmasks
	write(desc, 20) masks(Nmasks)%included
20	format('Center of mass kinetic energy for ',i6,' atoms')
 end function COM_KE_add


subroutine COM_KE_calc(i)
	!arguments
	integer, intent(in)			::	i
	integer			:: info, IPIV(3,3),j
        integer                 :: ats
	real(kind=prec)	:: A(3,3), B(3), W(30), K(6), C(6)

	!locals
	real(kind=prec)	::	KE, IXX, IXY, IXZ, IYY, IYZ, IZZ, tot_KE_rot

	if(i < 1 .or. i > Nmasks) return

        ats = masks(i)%included
	frames(i)=frames(i) + 1
        allocate(tmp_coords(ats))
        tmp_coords = translate_coords(coords(i)%xyz)
	if (first_frame(i)) then
		call mask_get(masks(i), translate_coords(xin), tmp_coords)
                coords(i)%xyz = btranslate_coords(tmp_coords)
		prev_coords(i)%xyz = coords(i)%xyz
		first_frame(i) = .false.
	else
		prev_coords(i)%xyz = coords(i)%xyz
		call mask_get(masks(i), translate_coords(xin), tmp_coords)
                coords(i)%xyz = btranslate_coords(tmp_coords)
	end if		
        deallocate(tmp_coords)

	!split coords into x, y, and z coords
	coords_mass(i)%x = coords(i)%xyz(1::3)
	coords_mass(i)%y = coords(i)%xyz(2::3)
	coords_mass(i)%z = coords(i)%xyz(3::3)
	prev_coords_mass(i)%x = prev_coords(i)%xyz(1::3)
	prev_coords_mass(i)%y = prev_coords(i)%xyz(2::3)
	prev_coords_mass(i)%z = prev_coords(i)%xyz(3::3)
	

	!calculate center of mass
	mass_ave(i)%x = dot_product(coords_mass(i)%x(:),coords_mass(i)%mass)/tot_mass(i)
	mass_ave(i)%y = dot_product(coords_mass(i)%y(:),coords_mass(i)%mass)/tot_mass(i)
	mass_ave(i)%z = dot_product(coords_mass(i)%z(:),coords_mass(i)%mass)/tot_mass(i)
	prev_mass_ave(i)%x = dot_product(prev_coords_mass(i)%x(:),coords_mass(i)%mass)/tot_mass(i)
	prev_mass_ave(i)%y = dot_product(prev_coords_mass(i)%y(:),coords_mass(i)%mass)/tot_mass(i)
	prev_mass_ave(i)%z = dot_product(prev_coords_mass(i)%z(:),coords_mass(i)%mass)/tot_mass(i)


	!create coordinate set relative mass center
	rel_coords(i)%x = coords_mass(i)%x - mass_ave(i)%x
	rel_coords(i)%y = coords_mass(i)%y - mass_ave(i)%y
	rel_coords(i)%z = coords_mass(i)%z - mass_ave(i)%z
	prev_rel_coords(i)%x = prev_coords_mass(i)%x - prev_mass_ave(i)%x
	prev_rel_coords(i)%y = prev_coords_mass(i)%y - prev_mass_ave(i)%y
	prev_rel_coords(i)%z = prev_coords_mass(i)%z - prev_mass_ave(i)%z


	!calculate moment of inertia tensor
	IXX = dot_product( (rel_coords(i)%y)**2 + (rel_coords(i)%z)**2, coords_mass(i)%mass)
	IYY = dot_product( (rel_coords(i)%x)**2 + (rel_coords(i)%z)**2, coords_mass(i)%mass)
	IZZ = dot_product( (rel_coords(i)%y)**2 + (rel_coords(i)%x)**2, coords_mass(i)%mass)
	IXY = -one * sum( (rel_coords(i)%y) * (rel_coords(i)%x) * coords_mass(i)%mass )
	IXZ = -one * sum( (rel_coords(i)%x) * (rel_coords(i)%z) * coords_mass(i)%mass )
	IYZ = -one * sum( (rel_coords(i)%y) * (rel_coords(i)%z) * coords_mass(i)%mass )
	

	!calculate individual atom (relative?) velocities
	velocity(i)%x = (rel_coords(i)%x - prev_rel_coords(i)%x) / frame_length
	velocity(i)%y = (rel_coords(i)%y - prev_rel_coords(i)%y) / frame_length
	velocity(i)%z = (rel_coords(i)%z - prev_rel_coords(i)%z) / frame_length



!	write (*,*) ' '

	
	K(:) = (/IXX,IXY,IYY,IXZ,IYZ,IZZ/)

!	write(*,'(6f18.3)'), K(:)
	call EIGEN(K,A,3,0)
	
	do j = 1, 3
!		write (*,*) A(1,j),A(2,j),A(3,j),K(j*(j+1)/2)
!		write (*,*) info
		eigen_stuff(i)%evalue(j) = K(j*(j+1)/2)
		eigen_stuff(i)%evector(:,j) = A(:,j)
!		write (*,'(4f23.8)') eigen_stuff(i)%evalue(j), eigen_stuff(i)%evector(:,j)

	end do
	
	!eigen_stuff(i)%evector are the NORMALIZED principal axes of rotation

	!vector of dot_product(principal axis,atom coordinate)

	do j=1,3
		dp_vect(i)%dp = eigen_stuff(i)%evector(1,j) * &
		    rel_coords(i)%x + eigen_stuff(i)%evector(2,j) * &
		    rel_coords(i)%y + eigen_stuff(i)%evector(3,j) * &
		    rel_coords(i)%z

		rad_vec(i,j)%x = rel_coords(i)%x - dp_vect(i)%dp*eigen_stuff(i)%evector(1,j)
		rad_vec(i,j)%y = rel_coords(i)%y - dp_vect(i)%dp*eigen_stuff(i)%evector(2,j)
		rad_vec(i,j)%z = rel_coords(i)%z - dp_vect(i)%dp*eigen_stuff(i)%evector(3,j)

	end do
	!now rad_vec(i,j)%x,y,z contains all the vectors of each atom to the principal axis j
	
	
	!time to take the cross-product of each atoms rad_vec with it's linear momentum vector

!	do j=1,3	
!		write (*,'(3f20.12)') sum(rad_vec(i,j)%x),sum(rad_vec(i,j)%y),sum(rad_vec(i,j)%z)
!	end do
	



	do j= 1,3
		ang_momentum(i,j)%x = sum(rad_vec(i,j)%y*velocity(i)%z*coords_mass(i)%mass - &
		rad_vec(i,j)%z*velocity(i)%y*coords_mass(i)%mass)
		ang_momentum(i,j)%y = sum(rad_vec(i,j)%z*velocity(i)%x*coords_mass(i)%mass - &
		rad_vec(i,j)%x*velocity(i)%z*coords_mass(i)%mass)
		ang_momentum(i,j)%z = sum(rad_vec(i,j)%x*velocity(i)%y*coords_mass(i)%mass - &
		rad_vec(i,j)%y*velocity(i)%x*coords_mass(i)%mass)
!		write (*,'(3f23.15)') ang_momentum(i,j)%x,ang_momentum(i,j)%y,ang_momentum(i,j)%z
	end do

	tot_KE_rot = 0
	!Rotational kinetic energy = Lj^2/(2Ij)
	do j=1,3
		KE_rot(i,j) = ((ang_momentum(i,j)%x)**2+(ang_momentum(i,j)%y)**2+ &
		(ang_momentum(i,j)%z)**2)/(2*eigen_stuff(i)%evalue(j))*conversion_factor
		tot_KE_rot = tot_KE_rot + KE_rot(i,j)
	end do
	
!	write (*,'(3f23.15)') KE_rot(i,1),KE_rot(i,2),KE_rot(i,3)
	
	
!	calculate the kinetic energy
	KE = 0.5 * tot_mass(i) * ((prev_mass_ave(i)%x-mass_ave(i)%x)**2+ &
	    (prev_mass_ave(i)%y-mass_ave(i)%y)**2+(prev_mass_ave(i)%z- &
	    mass_ave(i)%z)**2) / frame_length**2 * conversion_factor


	
	

	
 	write(*,100, advance='no') KE, tot_KE_rot
100	format(2f12.3)
end subroutine COM_KE_calc


subroutine COM_KE_put_mass(i)
	integer						::	k,j,i,at
	real						::	mass

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

end subroutine COM_KE_put_mass

subroutine COM_KE_heading(i)
	integer						::	i

	write(*,'(a)', advance='no') 'Trans  Rot (kcal/mol)'
end subroutine COM_KE_heading
                                         

end module CALC_COM_KE
