! (C) 2014 Uppsala Molekylmekaniska HB, Uppsala, Sweden
! calc_entropy.f90
! By Jens Carlsson, PhD 
! Absolute entropies according to Quasiharmonic analysis, Schlitter's formula and RMS fluctuations
!TODO: precision not fixed

module CALC_ENTROPY
    
	use CALC_BASE
	use MASKMANIP
    use TRJ
	use CALC_FIT
	use CALC_GEOM
implicit none

! Constants

real(8)										:: pi_val, hs, R, ref_vol, e, u, kb, cm   ! pi, plancs constant, R, reference_volume

!module variables


	type(MASK_TYPE), private, target		:: masks(MAX_MASKS)
	integer, private						:: Nmasks = 0

    type COVARIANCE_MATRIX                                                  ! Numbers stored for covariance matrix calculation
	    real(8)								:: Exy, Ex, Ey					! C[X,Y] = E[XY] - E[X] -E[Y]
	end type COVARIANCE_MATRIX


	type ENTROPY_COORD_TYPE
		real, pointer						:: x(:), xref(:)				! , Reference structure for rotational fit 				
        real(8), pointer					:: xr(:)
		real(8)								:: xrcm(3) 
		integer								:: interval, Calc_type			! Interval of calculations, 
		character(len=9)					:: Calc_select
		real(8)								:: Temperature					! Temperature
	    integer								:: tot_no_frames				! Total number of frames
	    
	end type ENTROPY_COORD_TYPE
	
	type(ENTROPY_COORD_TYPE), private		:: coords(MAX_MASKS)
	integer :: entropy_frame, start

    type trajectory
	    type(COVARIANCE_MATRIX), pointer	 :: COV_DATA(:,:)													! E[XY], E[X], E[Y] - for calculation of Covariance matrix
		real(8), pointer					 :: C(:,:),  VIZ(:), massvector(:), Xcm(:,:), Xrot(:) , X(:,:)
	    integer								 :: startframe, endframe, no_atoms, no_frames, store_cov, Temperature, Calc_type 
		! Under development:
		real(8), dimension(3,3)			     :: ROT    
		real(8)								 :: masstot					! Total mass
     end type trajectory
   
contains

!---------------
!- Entropy Add -
!---------------

integer function entropy_add(desc)
	!arguments
	character(*)				::	desc
	integer						::	ats
	
	if(Nmasks == 1) then
		write(*,10) 1
		entropy_add = 0
		return
	end if
 
10	format('Sorry, the maximum number of entropy calculations is ',i2)
	

	! Open output file
	open(100, file='Entropy.out',err=999)		! Entropy output file
999 write(*,'(a)') 'WARNING: Unable to initialize entropy output file. If you are not doing entropy calculations and &
	&multiple instances of Q are running in the same directory, you may ignore this warning.'

	!Add a new entropy mask
	
	Nmasks = Nmasks + 1
	
	call mask_initialize(masks(Nmasks))
	ats =  maskmanip_make(masks(Nmasks))
    
	WRITE(*, '(a)') &
	'', &
	'===================== Entropy calculation settings 1 ==============================', &
	'   description                                                command   ', &
	' 1 Schlitters formula                                            1      ', &
	' 2 Quasiharmonic analysis                                        2      ', & 
	' 3 Schlitters formula for masscenter of molecule                 3      ', &
	' 4 Principal variances (translation)                             4      ', &
	' 5 RMS deviations in Euler angles (Rotation)                     5      '
	WRITE(*, '(a)', advance = 'no') &
	'Qcalc> '
	read(*,*) coords(Nmasks)%Calc_type
	
	WRITE(*, '(a)') &
	'', &
	'===================== Entropy calculation settings 2 ==============================', &
	'   description                                                command   ', &
	' 1 S[vibration,rotation,translation,configuration]            All       ', &
	' 2 S[vibration,rotation, conformation (vibration)]            No_transl ', &
	' 3 S[conformation (vibration)]                                No_rot    '
	call getlin(coords(Nmasks)%Calc_select, 'Qcalc> ')
	
	WRITE(*, '(a)') &
	'Number of frames:'
    
	WRITE(*, '(a)', advance = 'no') &
	'Qcalc> '
	read(*,*) coords(Nmasks)%interval
    
	! Temperature set to default 300 
	coords(Nmasks)%Temperature = 300

	WRITE(*, '(a)') &
    'Set no frames:'
	WRITE(*, '(a)', advance = 'no') &
	'Qcalc> '
	read(*,*) coords(Nmasks)%tot_no_frames
	
	if(ats == 0) then
		call mask_finalize(masks(Nmasks))
		Nmasks = Nmasks - 1
		entropy_add = 0
		return
	end if

    allocate(coords(Nmasks)%x(3*ats), coords(Nmasks)%xref(3*ats))
    
	entropy_add = Nmasks
	
	write(desc, 20) masks(Nmasks)%included
20	format('Entropy calculation for ',i6,' atoms')
    
end function entropy_add

!------------------------
!- Construct trajectory -
!------------------------

subroutine construct_trajectory(t, Temp, startframe, endframe, no_atoms, Calc_type)
 
    type(trajectory)		 :: t
	integer					 :: startframe, endframe, no_frames, no_atoms, i,j, Calc_type
    real(8)					 :: Temp
	character(len=20)		         :: dcd, top
    
	t%startframe = startframe						! Starting frame
    t%endframe = endframe							! Ending frame	
    t%no_frames = endframe-startframe+1				! Calculate total number of frames
    t%no_atoms = no_atoms							! Number of atoms in mask
		t%Temperature = Temp							! Temperature
    t%Calc_type = Calc_type							! Store choice of calculation
    
 end subroutine construct_trajectory

!-------------------------------------------
!- Subroutine: Entropy initialize/finalize -
!-------------------------------------------

subroutine ENTROPY_initialize

! Set constants
 pi_val = 4 * atan(1.)								 ! pi
 e = exp(1.0)										 ! Eulers number
 R = 1.9855026										 ! In cal/mol
 kb=1.38066											 ! Boltzmann konstant  *10**(-23) 
 hs=1.0545727										 ! Plancs constant     *10**(-34) 
 u=1.6605402										 ! Atomic mass constant *10**(-27) from u to Kg
 ref_vol = 1660.0									 ! Corresponds to a 1M concentration
 cm = 3.3356410/(2*pi_val)						     ! Constant for vibrational frequencies calculation 

 ! Other variables
 entropy_frame = 1
 start = 1

 return

end subroutine ENTROPY_initialize

subroutine ENTROPY_finalize(i)
 integer						::	i
	call mask_finalize(masks(i))
    deallocate(coords(i)%xr, coords(i)%x, coords(i)%xref)
end subroutine ENTROPY_finalize

!###################
!- Entropy heading -
!###################

subroutine entropy_heading(i)                      
 integer						::	i
	write(*,'(a)', advance='no') 'Entropy(cal/mol K)'
end subroutine entropy_heading

!#######################################
!- Main entropy calculation subroutine -
!#######################################

subroutine entropy_calc(i)

 integer, intent(in)		:: i
 integer					:: j,k,l,n,m,o
 real(8)					:: det, gauss, uniform !, ref_vol
 type(trajectory)			:: t                                                                     

 real(8), dimension(3)	    :: uin,vin        
 ! Go on with entropy calculation
 
  if( start == 1) then
    
	allocate(coords(i)%xr(3*nat_pro))
    allocate( t%Xrot(3) )
	 
    select case(coords(i)%Calc_type)
	 case(1:2)  
	  allocate(t%COV_DATA(1:masks(i)%included*3,masks(i)%included*3)) 
	 case(3:5) 
	  allocate(t%COV_DATA(3,3)) 
	end select 
   
    call fit_make_reference(i)			! Reference for rotational fit      
                  
    
    start=0
  end if                            
 
  call construct_trajectory(t, coords(i)%Temperature, 1, entropy_frame, masks(i)%included,coords(i)%Calc_type ) 
    
   select case(t%Calc_type)
      case(1:2)
	    if( coords(i)%Calc_select == 'No_rot' .OR. coords(i)%Calc_select == 'No_transl') then 
         call RotationAndTranslation(t,i)
        end if
	  case(3:5)
	   call RotationAndTranslation(t,i)
   end select
  
  ! Locate which atoms to make the calculation on
   call mask_get(masks(i), xin, coords(i)%x)    
   
  ! Uppdate COV_DATA for Schlitters formula and Quasiharmonics   
  if(  t%Calc_type == 1 .OR. t%Calc_type == 2 ) then 
  
      allocate(t%massvector(3*t%no_atoms))   
      n=1
     do l = 1,nat_pro                                  ! nat_pro is the total number of atoms in topology
       if(masks(Nmasks)%mask(l)) then
        t%massvector(3*n-2) = iaclib(iac(l))%mass                             
        t%massvector(3*n-1) = iaclib(iac(l))%mass  
	    t%massvector(3*n)   = iaclib(iac(l))%mass     
	    n=n+1	  
       endif
     end do

     j=1
     k=1
      
	 do j=1,t%no_atoms*3
      do k =1,j
   	
	  t%COV_DATA(j,k)%Exy = (t%COV_DATA(j,k)%Exy*(t%no_frames-1) + coords(i)%x(j)*coords(i)%x(k))/t%no_frames
  	  t%COV_DATA(j,k)%Ex  = (t%COV_DATA(j,k)%Ex*(t%no_frames-1)  + coords(i)%x(j))/t%no_frames
	  t%COV_DATA(j,k)%Ey  = (t%COV_DATA(j,k)%Ey*(t%no_frames-1)  + coords(i)%x(k))/t%no_frames
   	
      t%COV_DATA(k,j)%Exy = t%COV_DATA(j,k)%Exy
 	  t%COV_DATA(k,j)%Ex  = t%COV_DATA(j,k)%Ex
	  t%COV_DATA(k,j)%Ey  = t%COV_DATA(j,k)%Ey
    
	  enddo
     enddo
     
  endif  ! Schlitter or quasi?
 
  if(mod(entropy_frame,coords(i)%interval) == 0)  then          ! Check if entropy is to be calculated
   
   select case( t%Calc_type )									! Check which type of calculation is to be performed
   
    !---------------------
    ! Schlitter's formula 
    !---------------------
   
    case(1)                                              
      call CovarianceMatrix(t)
      call SchlittersFormula(t,det)  
	  write(*,*) 0.5*R*det
      write(100,*) entropy_frame, 0.5*1.9855026*det 
   
    !------------------------
    ! Quasiharmonic analysis 
    !------------------------
   	     
    case(2)                                                       
	  call CovarianceMatrix(t) 
      call Quasiharmonic_analysis(i,t,det)   
      write(*,*) R*det                                       ! 0.5* Used before to be compatible to schlitters formula
      write(100,*) R*det                                     ! 0.5*
   
    !------------------------------------------
    ! Schlitter's formula using center of mass 
    !------------------------------------------

    case(3)                                           
     
	 t%no_atoms = 1
	 call CovarianceMatrix(t)
	 call SchlittersFormula(t, det)   
     
	 write(*,*) 0.5*R*det
     write(100,*) 0.5*R*det 
   
    !-----------------------------------------------
    ! Principal RMS fluctuations for center of mass 
    !-----------------------------------------------  

   case(4)                                    
     
   ! ref_vol = 1660.0   !34.123**3        !1660.0         !133.565**3

    call PCA(t,det)
	gauss = R*log( sqrt(det)*(2*pi_val*exp(1.0))**(1.5)/ref_vol )
	uniform = R*log( sqrt(det)*(12)**(1.5)/ref_vol )
    write(*,*) 'Gaussian approximation:  dS(transl) =', gauss
    write(*,*) 'Uniform approximation:      dS(transl) =', uniform
    write(100,*) gauss, uniform

    !----------------------------------
    ! RMS fluctuations in Euler angles 
    !----------------------------------

   case(5) 
     
	 call PCA(t,det)
     gauss = R*log( sqrt(det)*sin(t%COV_DATA(1,1)%Ex)*(2*pi_val*exp(1.0))**(1.5)/(8*pi_val**2) )
	 uniform = R*log( sqrt(det)*sin(t%COV_DATA(1,1)%Ex)*(12)**(1.5)/(8*pi_val**2) )
	 write(*,*) 'Gaussian approximation:  dS(rot) =', gauss
     write(*,*) 'Uniform approximation:   dS(rot) =', uniform
	 write(100,*) entropy_frame, gauss, uniform
   
   end select 
  
 endif  


  ! Deallocate vector containing all coordinates
  if( entropy_frame == coords(i)%tot_no_frames ) then
    deallocate(t%massvector)
  end if

  entropy_frame = entropy_frame + 1

end subroutine entropy_calc

!#################################
!- Covariance matrix calculation -
!#################################
 
subroutine CovarianceMatrix(t)							 
 
  type(trajectory)	    :: t
  integer				:: i,j 
  
   ! Allocate covariance matrix
  allocate( t%C(1:t%no_atoms*3, 1:t%no_atoms*3)  )
  
  do i=1,t%no_atoms*3
    do j=1,i
   	t%C(i,j) = t%COV_DATA(i,j)%Exy  - t%COV_DATA(i,j)%Ex*t%COV_DATA(i,j)%Ey   	
	t%C(j,i) = t%C(i,j)
	enddo
  enddo

end subroutine CovarianceMatrix

!#####################################################
!- Calculation of Schlitter's formula type entropies -
!#####################################################

! Evaluates the determinant of a coleskyfactorized matrix

subroutine SchlittersFormula(t, det)

 type(trajectory)			:: t
 integer					:: i, k
 real(8)					:: det
   
 call MultiplyMatrices(t)            !  MCM + 1
 call CholeskyFactorization(t)       !  Coleskyfactorized matrix, L (C = L*L ), is stored in C 
 
 ! Calculate entropy
 det = 0.0
 do i=1,t%no_atoms*3
  det = det + 2*log(t%C(i,i))        !  det L*L =  2 * det L 
 end do                                               

 deallocate(t%C) ! Deallocate covariance matrix

end subroutine SchlittersFormula

!################################################
! Multiplies Matrice with Schlitter's constants -
!################################################

subroutine MultiplyMatrices(t)
 
   type(trajectory)				:: t
   integer						:: i,j,k,l,n
   real(8)						:: det, c !, kb, hs, e, c, u

  ! Set constants

 ! kb=1.38066       ! Boltzmann konstant  *10**(-23) 
 ! hs=1.0545727    ! Plancs constant     *10**(-34) 
 ! e=2.718281828       ! Euler number
 ! u=1.6605402      ! Atomic mass constant *10**(-27) from u to Kg
   
  c=(e**2)*kb*(t%Temperature)*u*0.01/(hs**2)        ! 0.01 is Å to m, u is u to kg (hs*hs) to hs**2
  
  do i=1,t%no_atoms*3
     do j=1,i                     !t%no_atoms*3	   
	       
       ! Multiply with constant
          t%C(i,j) = c*t%C(i,j)*sqrt(t%massvector(i)*t%massvector(j))                                  
     	  t%C(j,i) = t%C(i,j)     ! Symmetrix matrix
	   
       ! Add diagonal matrix to massweigted matrix        	   
  	   if( i == j ) then
         t%C(i,i) = t%C(i,i) + 1.0          
       endif

    enddo
   enddo

end subroutine MultiplyMatrices    

!#############################################
!- Choleskyfactorization of symmetric matrix -
!#############################################

! Coleskyfactorizes the covariance matrix and
! and stores the result (L) in C, C = L*L
! Heath, M., Scientific Computing - An Introduction Survey

subroutine CholeskyFactorization(t)           
											 
 type(trajectory)		:: t
 integer				:: i,j,k,l

 do j = 1,t%no_atoms*3
   do k = 1,j-1
     do i = j,t%no_atoms*3
       t%C(i,j) = t%C(i,j) - t%C(i,k)*t%C(j,k)
     enddo
   enddo
        if( t%C(j,j) < 0.0 ) then
          print*, "Matrix is not positive definite!"         ! Should hold.
        endif
        t%C(j,j) = sqrt(t%C(j,j)) 
   do k=j+1,t%no_atoms*3
        t%C(k,j) = t%C(k,j)/t%C(j,j)
   end do
 end do

end subroutine CholeskyFactorization

!#######################################################
!- Removal of Translational and/or rotational movement -
!#######################################################

subroutine RotationAndTranslation(t,i)
 
 integer, intent(in)   :: i
 integer               :: at, j, k
 integer               :: lsqstatus
 real(8)               :: totmass, error
 real(8)               :: xcm(3)	
 type(trajectory)      :: t
 
	!calc. mass centre of xin
	xcm(:) = 0.
	totmass = 0
	do at = 1, nat_pro
		if(masks(i)%mask(at)) then                          ! Take only the ones specified in mask
			xcm(:) = xcm(:) + xin(3*at-2:3*at)*iaclib(iac(at))%mass                 
		    totmass = totmass + iaclib(iac(at))%mass 
		end if
    end do

	xcm(:) = xcm(:)/totmass                       
   
  allocate(t%massvector(3)) 
	t%massvector(1)	= totmass  
    t%massvector(2) = totmass  
    t%massvector(3) = totmass  

  if( t%Calc_type == 3 .OR. t%Calc_type == 4) then    ! If we want to calculate the center of mass motion
    
	! Uppdate covariance matrix    
    do j=1,3
      do k =1,j
    	t%COV_DATA(j,k)%Exy = (t%COV_DATA(j,k)%Exy*(t%no_frames-1) + xcm(j)*xcm(k))/t%no_frames
  		t%COV_DATA(j,k)%Ex  = (t%COV_DATA(j,k)%Ex*(t%no_frames-1)  + xcm(j))/t%no_frames
		t%COV_DATA(j,k)%Ey  = (t%COV_DATA(j,k)%Ey*(t%no_frames-1)  + xcm(k))/t%no_frames
	
		t%COV_DATA(k,j)%Exy = t%COV_DATA(j,k)%Exy
		t%COV_DATA(k,j)%Ex  = t%COV_DATA(j,k)%Ex
		t%COV_DATA(k,j)%Ey  = t%COV_DATA(j,k)%Ey 
       enddo
    enddo

   endif

	if( coords(i)%Calc_select == 'No_transl' ) then               ! Translation is removed

     ! Remove translational motion, shift xin to origin 
	  do at = 1, nat_pro
		xin(3*at-2) = xin(3*at-2) - xcm(1)      
		xin(3*at-1) = xin(3*at-1) - xcm(2)
		xin(3*at  ) = xin(3*at  ) - xcm(3)
	  enddo
	
	elseif( coords(i)%Calc_select == 'No_rot' ) then              ! Translation and Rotation is removed
	
     ! Remove translational motion and massweight, shift xin to origin 
	  do at = 1, nat_pro
		xin(3*at-2) = (xin(3*at-2) - xcm(1))*sqrt(iaclib(iac(at))%mass)       
		xin(3*at-1) = (xin(3*at-1) - xcm(2))*sqrt(iaclib(iac(at))%mass)
		xin(3*at  ) = (xin(3*at  ) - xcm(3))*sqrt(iaclib(iac(at))%mass)
	  end do
	
      ! rotate xin to fit with coords(i)%xr
      lsqstatus = LSQSTR(nat_pro, masks(i)%mask(:), coords(i)%xr(:), xin(:), error, t%ROT )
      
	  if( lsqstatus == 1 ) then
 	    print*, 'Error: Least sqares fit failed'
      end if
	  
	  if( t%Calc_type == 5) then 
       call Euler_angles(t)  
      endif
	  
       ! Divide the output coordinates with the mass of the particles.
	   ! xrcm is not weighted in fit make reference and is not needed for entropy calculation (it keeps things at zero)
          
	  do at = 1, nat_pro
	    xin(3*at-2) = xin(3*at-2)/sqrt(iaclib(iac(at))%mass)  
	    xin(3*at-1) = xin(3*at-1)/sqrt(iaclib(iac(at))%mass) 
	    xin(3*at  ) = xin(3*at  )/sqrt(iaclib(iac(at))%mass)  
	  end do
 
    endif  

end subroutine RotationAndTranslation

!################################
!- Makes reference for RMS fit  -
!################################

subroutine fit_make_reference(i)
	!arguments
	integer						::	i
	!locals
	integer						::	at
    real(8)                     ::  totmass

	if(i < 1 .or. i > Nmasks) return
	!calc centre vector of mass
	 coords(i)%xrcm(:) = 0.

   	totmass = 0.0

	do at = 1, nat_pro
		if(masks(i)%mask(at)) then    
			coords(i)%xrcm(:) = coords(i)%xrcm(:) + xin(3*at-2:3*at)*iaclib(iac(at))%mass
		    totmass = totmass + iaclib(iac(at))%mass
		end if
	end do
	
	coords(i)%xrcm(:) = coords(i)%xrcm(:)/totmass    ! Center of mass
	
	! shift to origin and massweight coordinates in order to get massweighted RMS fit
	do at = 1, nat_pro
	  coords(i)%xr(3*at-2) = (xin(3*at-2) - coords(i)%xrcm(1))*sqrt(iaclib(iac(at))%mass) 
	  coords(i)%xr(3*at-1) = (xin(3*at-1) - coords(i)%xrcm(2))*sqrt(iaclib(iac(at))%mass) 
	  coords(i)%xr(3*at  ) = (xin(3*at  ) - coords(i)%xrcm(3))*sqrt(iaclib(iac(at))%mass) 
	end do

end subroutine fit_make_reference

!########################################################
!- Calculation of Quasiharmonic analysis type entropies -
!########################################################

! Performs Quasiharmonic analysis according to: 
! Andricioaei, I. and Karplus, M. (2001) Journal of Chemical Physics 115, 6289-6292.

subroutine Quasiharmonic_analysis(i, t, det)

type(trajectory)					:: t
real(8), dimension(:,:), pointer	:: R     ! can I remove this
real(8)								:: det
integer								:: i

call massweight(t)                      ! Massweighting of the Covariancematrix sqrt(M)*C*sqrt(M)
call VIZ(t)                             ! Transformation to storage for symmetrix matrix
call eigen(t%VIZ,R,t%no_atoms*3,1)      ! Calculation of eigenvalues (vibrational frequencies)
call entropy_ho(i,t,det)                ! Calculation of entropy from Harmonic oscillator

end subroutine Quasiharmonic_analysis

!#########################
!- Subroutine massweight -
!#########################

subroutine massweight(t)             ! Massweighting of coordinates
 
 type(trajectory)    :: t
 integer             :: i,j
 
 do i=1,t%no_atoms*3
   do j = 1,i               ! Symmetrix matrix
    t%C(i,j) = sqrt(t%massvector(i)*t%massvector(j))*t%C(i,j)
    t%C(j,i) = t%C(i,j)

!	if( i == j ) then
!     t%C(i,j) = t%C(i,j)          + 1.0
!	endif

   end do
 end do

 end subroutine massweight

!############################################
!- VIZ - smart storage for symmetric matrix -
!############################################

subroutine VIZ(t)

type(trajectory) :: t
integer :: i,j,k 
 
 allocate( t%VIZ(3*t%no_atoms*(3*t%no_atoms+1)/2) )
 
 do i=1,t%no_atoms*3
    do j=1,t%no_atoms*3      
     k = i*(i-1)/2+j
     t%VIZ(k) = t%C(i,j)   
    end do
 end do
 
 deallocate(t%C)          

end subroutine VIZ

!####################################
!- Entropy from harmonic oscillator -
!####################################

subroutine entropy_ho(m, t, det)

type(trajectory) :: t
integer :: i,k,no_freq,m
real(8) :: det, hswkT, w, cm, scaling ! , hs, kb, u, ,pi

! kb = 1.380658				! 10**(-23)
! hs = 1.05457266			! 10**(-34)
! u  = 1.6605402		    ! 10**(-27), u to Kg
! pi = acos(-1.0)			! Pi
 cm = 3.3356410/(2*pi)      ! 10**(-11), s**(-1) to cm**(-1)
 
! Different number of eigenvalues are used depending if
! any degrees of freedom have been removed

 if( coords(m)%Calc_select == 'No_rot') then
   no_freq = 3*t%no_atoms-6
 elseif( coords(m)%Calc_select == 'No_transl') then
   no_freq = 3*t%no_atoms-3
 else
   no_freq = 3*t%no_atoms
 endif
  
  WRITE(*, '(a)') ':'
  WRITE(*, '(a)') &
  '------------------------------------------------------------------------',&
  '         Frequency (cm**(-1))		  S (cal/mol Kelvin)'                                    

   det = 0.0
   do i= 1,no_freq       
     k = i*(i+1)/2                                          ! Transformation to VIZ format
     w = sqrt(kb*t%Temperature/(u*t%VIZ(k)))    ! Removed -1.0 since i don't have to shift  ! 10**(12) = sqrt(10**(-23, kb)/(10**(-27, u)*10**(-20, Å**2 to m**2))
     hswkT = 10*hs*w/(kb*t%Temperature)                     ! 10**(1) = 10**(-34, h)*10**(12, w)/10**(-23, kb)
     
    write(*,'(i10,a,3f8.3)', advance='no') i,'. ',10*cm*w   ! 10**(1) = 10**(12,w)*10**(-11,cm)
    write(*,'(a)', advance='no') '			'
	write(*,'(3f8.3)') 1.9855026*(hswkT/(exp(hswkT)-1) - log(1-exp(-hswkT)))    
	 
    det = det + hswkT/(exp(hswkT)-1) - log(1-exp(-hswkT))   
     
  enddo
  
  WRITE(*, '(a)') &
  '-------------------------------------------------------------------------'
  WRITE(*, '(a)', advance = 'no') &
  'S(Total) = ' 
  
  !det=det*2.0        ! Technical thing to make it compatible to calculations using Schlitter's formula - not needed anymore

  
end subroutine entropy_ho


!##################################
!- Principal components of motion -
!##################################

! Finds principal variances in translational motion

subroutine PCA(t,det)

type(trajectory) :: t
integer :: i,j,k, rotation
real(8), dimension(:,:), pointer :: R
real(8) :: det

! Can only be done if center off mass has been removed

if( t%Calc_type == 5 ) then
 !deallocate(t%X)
 !allocate(t%X(1:3,1:t%no_frames))
 !t%X = t%Xrot
 t%no_atoms = 1
 print*, '<theta> =', t%COV_DATA(1,1)%Ex
endif

! Remove masses from covariance calculation

 do i=1,3
  t%massvector(i) = 1.0
 enddo

 call CovarianceMatrix(t)
 call VIZ(t)
 
 if( t%Calc_type == 4 ) then
  call eigen(t%VIZ,R,3,1) 
 endif

 det=1
 print*, ''
 do j=1,3
  k = j*(j+1)/2
  print*, 'Variance in direction', j ,t%VIZ(k)
  det = det*t%VIZ(k)                       ! multiplying principal RMS deviations
 enddo 

    
end subroutine PCA

!###############################
!- Calculation of Euler angles -
!###############################

! Euler type angle from Classical dynamics of particles and systems, Marion Thornton 

subroutine Euler_angles(t)

type(trajectory) :: t
real(8) :: phi, psi, theta, delta
integer :: i,j,k
	
theta =  acos(t%ROT(3,3))											! Calculate Euler angles                   
phi   =  atan2(t%ROT(3,1)/sin(theta),-t%ROT(3,2)/sin(theta)) 
psi   =  atan2(t%ROT(1,3)/sin(theta),t%ROT(2,3)/sin(theta))   

t%Xrot(1) = theta
t%Xrot(2) = psi
t%Xrot(3) = phi  

! Uppdate covariance matrix

   do j=1,3
      do k =1,j
    	t%COV_DATA(j,k)%Exy = (t%COV_DATA(j,k)%Exy*(t%no_frames-1) + t%Xrot(j)*t%Xrot(k))/t%no_frames
  		t%COV_DATA(j,k)%Ex  = (t%COV_DATA(j,k)%Ex*(t%no_frames-1)  + t%Xrot(j))/t%no_frames
		t%COV_DATA(j,k)%Ey  = (t%COV_DATA(j,k)%Ey*(t%no_frames-1)  + t%Xrot(k))/t%no_frames
	
		t%COV_DATA(k,j)%Exy = t%COV_DATA(j,k)%Exy
		t%COV_DATA(k,j)%Ex  = t%COV_DATA(j,k)%Ex
		t%COV_DATA(k,j)%Ey  = t%COV_DATA(j,k)%Ey 
      enddo
    enddo
      
end subroutine Euler_angles

!#########################################
!- Subroutine to update Covariancematrix -
!#########################################

subroutine UpdateCovariancematrix(t,i)  ! NOTE: Not in use yet!

type(trajectory) :: t
integer :: i,j,k
	
	do j=1,t%no_atoms*3
      do k =1,j
   	
	  t%COV_DATA(j,k)%Exy = (t%COV_DATA(j,k)%Exy*(t%no_frames-1) + coords(i)%x(j)*coords(i)%x(k))/t%no_frames
  	  t%COV_DATA(j,k)%Ex  = (t%COV_DATA(j,k)%Ex*(t%no_frames-1)  + coords(i)%x(j))/t%no_frames
	  t%COV_DATA(j,k)%Ey  = (t%COV_DATA(j,k)%Ey*(t%no_frames-1)  + coords(i)%x(k))/t%no_frames
   	
      t%COV_DATA(k,j)%Exy = t%COV_DATA(j,k)%Exy
 	  t%COV_DATA(k,j)%Ex  = t%COV_DATA(j,k)%Ex
	  t%COV_DATA(k,j)%Ey  = t%COV_DATA(j,k)%Ey
    
	  enddo
     enddo

end subroutine UpdateCovariancematrix

    


end module CALC_ENTROPY
