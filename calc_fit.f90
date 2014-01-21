!	(C) 2000 Uppsala Molekylmekaniska HB, Uppsala, Sweden

!	calc_fit.f90
!	by John Marelius
!	parts (C) 1979 W. F. van Gunsteren
!	pair-wise superposition of coordidinate sets

module CALC_FIT
	use CALC_BASE
	use MASKMANIP
	implicit none

!constants
	integer, parameter			::	MAX_MASKS = 10
!module variables
	type(MASK_TYPE), private, target	::	masks(MAX_MASKS)
	integer, private			::	Nmasks = 0

	type FIT_COORD_TYPE
		real(8), pointer		::	xr(:)
		real(8)					::	xrcm(3)
	end type FIT_COORD_TYPE
	type(FIT_COORD_TYPE), private	::	coords(MAX_MASKS)

contains
	
subroutine fit_initialize
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
	real(8)						::	xcm(3), totmass, U(3,3), E

    if(i < 1 .or. i > Nmasks) return
	
	!calc. mass centre of xin
	xcm(:) = 0.
	totmass = 0.0
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
    real(8)                     :: totmass

	if(i < 1 .or. i > Nmasks) return
	!calc centre vector
	coords(i)%xrcm(:) = 0.
	totmass = 0.0
	do at = 1, nat_pro
		if(masks(i)%mask(at)) then
		   coords(i)%xrcm(:) = coords(i)%xrcm(:) + xtop(3*at-2:3*at)*iaclib(iac(at))%mass 
		   totmass = totmass + iaclib(iac(at))%mass 
		end if
	end do
	
	coords(i)%xrcm(:) = coords(i)%xrcm(:)/totmass !/for geometric masks(i)%included

	! shift to origin and massweight coordinates
	do at = 1,nat_pro
	 	coords(i)%xr(3*at-2) = (xtop(3*at-2) - coords(i)%xrcm(1))*sqrt(iaclib(iac(at))%mass)
		coords(i)%xr(3*at-1) = (xtop(3*at-1) - coords(i)%xrcm(2))*sqrt(iaclib(iac(at))%mass)
		coords(i)%xr(3*at  ) = (xtop(3*at  ) - coords(i)%xrcm(3))*sqrt(iaclib(iac(at))%mass)
	end do

end subroutine fit_make_ref


integer function LSQSTR (NR,W,XP,X,E,U)


!CCCCC W.F. VAN GUNSTEREN, CAMBRIDGE, JULY 1979 CCCCCCCCCCCCCCCCCCCCCCCC
!																	   C
!	  SUBROUTINE LSQSTR (NR,W,XP,X,E,IROT)							   C
!																	   C
!OMMENT   LSQSTR ROTATES THE ATOMS WITH COORDINATES X ABOUT THE ORIGIN C
!	  SUCH THAT THE FUNCTION E = 0.5 * SUM OVER ALL NR ATOMS OF 	   C
!	  W*(X-XP)**2 HAS A MINIMUM. HERE W DENOTES THE WEIGHT FACTORS AND C
!	  XP ARE THE REFERENCE COORDINATES. FOR A DESCRIPTION OF THE	   C
!	  METHOD SEE A.D. MCLACHLAN, J. MOL. BIOL. 128 (1979) 49.		   C
!	  IF THE SUBROUTINE FAILS, IT IS RETURNED WITH A MESSAGE AND	   C
!	  IROT=0.														   C
!																	   C
!	  NR = NUMBER OF ATOMS											   C
!	  W(1..NR) = ATOMIC WEIGHT FACTORS								   C
!	  XP(1..3*NR) = REFERENCE ATOM COORDINATES						   C
!	  X(1..3*NR) = ATOM COORDINATES; DELIVERED WITH THE ROTATED ONES   C
!	  If present,
!	  E is DELIVERED WITH THE MINIMUM VALUE OF THE FUNCTION E          C
!																	   C
!	  LSQSTR USES SUBR. EIGEN.										   C
!																	   C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!	converted to F90 function by John Marelius, July 2000
!	returns 0 on success, 1 if Det(U) < 1e-9, 2 if Det(U) < 0 and Omega has
!	degenerate eigenvalues
!	Weight vector W has been made logical

	implicit none
!arguments
	integer, intent(in)			::	NR
	logical(1), intent(in)		::	W(:)
	real(8), intent(in)			::	XP(:)
	real(8), intent(inout)		::	X(:)
	real(8), optional, intent(out)::	E
   
!locals
	real(8), optional           ::	U(3,3)  
	real(8)						::	COM(21), OM(6,6)
	real(8)						::	VH(3,3), VK(3,3)
	real(8)						::	DU
	real(8)						::	SIG
	integer						::	I, J, M, M1, M2


!
!*****CALCULATE THE MATRIX U AND ITS DETERMINANT
	U(:,:) = 0.
	DO J=1,NR
		if(W(J)) then
			U(:,:) = U(:,:) + matmul(reshape(X(3*J-2:3*J),(/3,1/)), &
				reshape(XP(3*J-2:3*J),(/1,3/)))
		end if
	END DO
!
	DU= U(1,1)*U(2,2)*U(3,3)+U(1,3)*U(2,1)*U(3,2) &
		+U(1,2)*U(2,3)*U(3,1)-U(3,1)*U(2,2)*U(1,3) &
		-U(3,3)*U(2,1)*U(1,2)-U(3,2)*U(2,3)*U(1,1)
	IF ( ABS(DU) < 1.E-9) then
		LSQSTR = 1
		return
	end if

!*****CONSTRUCT OMEGA, DIAGONALIZE IT AND DETERMINE H AND K
	M=0
	DO M1=1,6
		DO M2=1,M1
			M=M+1
			IF(M1 > 3 .AND. M2 <= 3) then
				COM(M)=U(M2,M1-3)
			else
				COM(M)=0.
			end if
		END DO
	END DO
!
	CALL EIGEN (COM,OM,6,0)
	IF (DU < 0. .AND. ABS(COM(3)-COM(6)) < 1.E-5) then
		LSQSTR = 2
		return
	end if
!
	VH(:,:) = sqrt(2.) * OM(1:3, 1:3)
	VK(:,:) = sqrt(2.) * OM(4:6, 1:3)

	SIG= (VH(2,1)*VH(3,2)-VH(3,1)*VH(2,2))*VH(1,3) &
		+(VH(3,1)*VH(1,2)-VH(1,1)*VH(3,2))*VH(2,3) &
		+(VH(1,1)*VH(2,2)-VH(2,1)*VH(1,2))*VH(3,3)
	IF(SIG <= 0.) then
		VH(:,3) = -VH(:,3)
		VK(:,3) = -VK(:,3)
	end if
!
!*****DETERMINE R AND ROTATE X
	DO M2=1,3
		DO M1=1,3
			U(M1,M2)=VK(M1,1)*VH(M2,1)+VK(M1,2)*VH(M2,2) &
				+sign(1._8,DU)*VK(M1,3)*VH(M2,3)
		end do
	END DO
!
	DO J=1,NR
		X(3*J-2:3*J) = matmul(U,X(3*J-2:3*J))
	END DO
	LSQSTR = 0 !OK
!
!*****CALCULATE E, WHEN REQUIRED
	if(present(E)) then
		E=0.
		DO J=1,NR
			if(W(J)) E=E+sum((X(3*J-2:3*J)-XP(3*J-2:3*J))**2)
		end do
		E=E/2.
	end if

END function LSQSTR


end module CALC_FIT
