! (C) 2014 Uppsala Molekylmekaniska HB, Uppsala, Sweden
! calc_rdf.f90
! by Martin Andér and Martin Almlöf (multiple atoms in first mask and PBC)
!       
! Calculates RDF for a center atom with respect to all atoms in
! target_mask, using Nbins bins from r=0 to r=rdf_radius
!TODO: precision not fixed

module CALC_RDF
	use CALC_BASE
	use MASKMANIP
	use PARSE
	implicit none

! constants
	integer, parameter							::	MAX_MASKS = 10
! module variables

	type(MASK_TYPE), private, target			::	temp_masks(2)	! masks for making pair lists
	
	integer, private							::	ats1,ats2,Nmasks = 0				! # target masks

	type RDF_CALC_TYPE
		integer									::	atoms_in_mask1					! center atom idx
		integer									::	frames = 0				! # frames used
		integer									::	Nbins					! Divide rdf_radius into Nbins bins
		real(kind=prec)									::	rdf_radius				! Calculate RDF within this radius
		real(kind=prec),pointer							::	bins(:)					! RDF bins
	end type RDF_CALC_TYPE

	type(RDF_CALC_TYPE), private, save				::  rdf_calcs(MAX_MASKS)	! stores data for indivdual RDF calcs
	type RDF_LIST_TYPE
		integer							::  number_of_pairs
		integer, pointer				::	atom1(:), atom2(:)
	end type RDF_LIST_TYPE
	type(RDF_LIST_TYPE), private, allocatable::  rdf_listan(:)

contains

subroutine RDF_initialize
	allocate(rdf_listan(MAX_MASKS))
end subroutine RDF_initialize

subroutine RDF_finalize(i)
	! arguments
	integer				::	i		! calculation #

	! locals			
	integer				::	j		! atom index
	real(kind=prec)				::	volume	! shell volume

	! write output header 1
	write (*,3) i
3	format('Results from RDF calulation' , i3 , ', normalized with respect to shell volume and frame count')

	! write output header 2
	write (*,4)
4	format('    Bin  r_max    RDF')
	
	! normalize the contents of each bin with respect to shell
	! volume and # processed frames and print the results
	do j = 1, rdf_calcs(i)%Nbins
		volume = 4.0_prec * pi / 3.0_prec * ((j * rdf_calcs(i)%rdf_radius / &
		rdf_calcs(i)%Nbins)**3 - ((j-1) * rdf_calcs(i)%rdf_radius / rdf_calcs(i)%Nbins)**3)
		rdf_calcs(i)%bins(j) = rdf_calcs(i)%bins(j) / volume / &
		rdf_calcs(i)%frames / rdf_calcs(i)%atoms_in_mask1
		write (*,5) j , j * rdf_calcs(i)%rdf_radius / rdf_calcs(i)%Nbins , rdf_calcs(i)%bins(j)
	end do
5	format(i7 , f10.6 , f10.6)

	! finalize and deallocate
	deallocate(rdf_calcs(i)%bins)			! deallocate bins
	deallocate(rdf_listan(i)%atom1)
	deallocate(rdf_listan(i)%atom2)
end subroutine RDF_finalize


!*********************************************************
! Qcalc calls this function to set up one RDF calculation,
! which is then added to the list of calculations in Qcalc
!*********************************************************
integer function RDF_add(desc)
	! arguments
	character(*)				::	desc

	! locals
	integer						::	ats

	! check the # of calculations
	if(Nmasks == MAX_MASKS) then
		write(*,10) MAX_MASKS
		return
	end if
10	format('Sorry, the maximum number of RDF calculations is ',i2)

	! start adding new calculation
	Nmasks = Nmasks + 1													! set current calculation
	write(*, 1011)
1011	format('Enter mask for first atom set in RDF calculation')
	call mask_initialize(temp_masks(1))
	ats1 =  maskmanip_make(temp_masks(1))
	if(ats1 == 0 ) then
		call mask_finalize(temp_masks(1))
		RDF_add = 0
		return
	end if
	
	write(*, 1012)
1012	format('Enter mask for second atom set in RDF calculation')
	call mask_initialize(temp_masks(2))
	ats2 =  maskmanip_make(temp_masks(2))

	if(ats2 == 0) then
		call mask_finalize(temp_masks(1))
		call mask_finalize(temp_masks(2))
		RDF_add = 0
		return
	end if


	!create rdf_listan
	allocate(rdf_listan(Nmasks)%atom1(ats1*ats2))
	allocate(rdf_listan(Nmasks)%atom2(ats1*ats2))
	!rdf_make_list returns number of pairs	
	rdf_listan(Nmasks)%number_of_pairs = rdf_make_list(temp_masks(1),temp_masks(2),rdf_listan(Nmasks))

	call mask_finalize(temp_masks(1))
	call mask_finalize(temp_masks(2))

	! set the parameters for this calculation
	rdf_calcs(Nmasks)%atoms_in_mask1 = ats1
	rdf_calcs(Nmasks)%rdf_radius = get_real_arg('RDF radius      :')	! get RDF radius
	rdf_calcs(Nmasks)%Nbins      =  get_int_arg('Number of bins  :')	! get # bins
	allocate(rdf_calcs(Nmasks)%bins(rdf_calcs(Nmasks)%Nbins))			! allocate bins
	rdf_calcs(Nmasks)%bins(:) = 0										! set bins(:) = 0

	RDF_add = Nmasks
	write(desc, 25) ats1
25	format('RDF calculation for ',i6,' atoms')
end function RDF_add


!*****************************************************************
! Calculates RDF data from the current trajectory frame
!*****************************************************************
subroutine RDF_calc(i)
	!arguments
	integer, intent(in)	::	i

	!locals
	integer				::	j		! atom index
	integer				::	binidx	! bin index
	real(kind=prec)			::	dist	! distance

!	do j = 1, nat_pro
!		if(target_masks(i)%mask(j))	then					! if atom j in this target_mask
!			dist = r(rdf_calcs(i)%center, j)				! calculate distance between center atom and atom j
!			if(dist < rdf_calcs(i)%rdf_radius) then			! if atom j is within rdf_radius of center atom
!				binidx = dist / rdf_calcs(i)%rdf_radius * rdf_calcs(i)%Nbins + 1
!				rdf_calcs(i)%bins(binidx) = rdf_calcs(i)%bins(binidx) + 1
!			end if
!		end if
!	end do


	do j = 1, rdf_listan(i)%number_of_pairs
		dist = r(rdf_listan(i)%atom1(j),rdf_listan(i)%atom2(j))
		if(dist < rdf_calcs(i)%rdf_radius) then	
			binidx = dist / rdf_calcs(i)%rdf_radius * rdf_calcs(i)%Nbins + 1
			rdf_calcs(i)%bins(binidx) = rdf_calcs(i)%bins(binidx) + 1
		end if
	end do






	rdf_calcs(i)%frames = rdf_calcs(i)%frames + 1			! update # processed frames

 	write(*,100, advance='no') 	rdf_calcs(i)%atoms_in_mask1
100	format(i10)
end subroutine RDF_calc


subroutine RDF_heading(i)
	integer						::	i

	write(*,'(a)', advance='no') 'RDF ats'
end subroutine RDF_heading

!*****************************************************
! r(a,b) calculates the distance between atoms a and b
!*****************************************************
real(kind=prec) function r(a,b)
	! arguments
	integer				::	a,b		!atom indices

	! locals
	real,dimension(3)	::	delta
	
	delta = xin(3*a-2:3*a)-xin(3*b-2:3*b)

	! if PBC then adjust lengths according to periodicity - MA
	if( use_PBC ) then
		delta(1) = delta(1) - boxlength%x*nint( delta(1)/boxlength%x )
		delta(2) = delta(2) - boxlength%y*nint( delta(2)/boxlength%y )
		delta(3) = delta(3) - boxlength%z*nint( delta(3)/boxlength%z )
	end if

	r = sqrt(dot_product(delta,delta))

end function r




!*****************************************************
! rdf_make_list makes a list of all atom atom pairs in the masks and return the number of pairs
!*****************************************************
integer function rdf_make_list(mask1, mask2, rdf_list)
	integer			:: i,j,k,l,m,nat1,nat2
	integer, allocatable	::	group1(:), group2(:)
	type(MASK_TYPE), intent(in) :: mask1, mask2
	type(RDF_LIST_TYPE)			:: rdf_list

	l = 1
	m = 1


	!make the list of atoms in first mask
	nat1 = mask1%included 
	allocate(group1(nat1))
	do j = 1, nat_pro
		if (mask1%MASK(j)) then
			group1(l) = j
			l = l+1
		end if
	end do

	!make the list of atoms in second mask
	nat2 = mask2%included 
	allocate(group2(nat2))
	do j = 1, nat_pro
		if (mask2%MASK(j)) then
			group2(m) = j
			m = m+1
		end if
	end do

	i=1
	do j = 1, nat1
		do k = 1, nat2
			if (group1(j) .ne. group2(k)) then
				rdf_list%atom1(i) = group1(j)
				rdf_list%atom2(i) = group2(k)
				i = i+1
			end if
		end do
	end do

	rdf_make_list = i-1    !set the return value

end function rdf_make_list






end module CALC_RDF
