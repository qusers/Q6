module QEIGEN
	use QMATH
        implicit none

contains

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
!arguments
	integer, intent(in)			::	NR
	logical(1), intent(in)		::	W(:)
	real(kind=prec), intent(in)			::	XP(:)
	real(kind=prec), intent(inout)		::	X(:)
	real(kind=prec), optional, intent(out)::	E
   
!locals
	real(kind=prec), optional           ::	U(3,3)  
	real(kind=prec)						::	COM(21), OM(6,6)
	real(kind=prec)						::	VH(3,3), VK(3,3)
	real(kind=prec)						::	DU
	real(kind=prec)						::	SIG
	integer						::	I, J, M, M1, M2


!
!*****CALCULATE THE MATRIX U AND ITS DETERMINANT
	U(:,:) = zero
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
	IF ( ABS(DU) < QREAL_EPS) then
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
				COM(M)=zero
			end if
		END DO
	END DO
!
	CALL EIGEN (COM,OM,6,0)
	IF (DU < zero .AND. ABS(COM(3)-COM(6)) < QREAL_EPS) then
		LSQSTR = 2
		return
	end if
!
	VH(:,:) = sqrt(2.0_prec) * OM(1:3, 1:3)
	VK(:,:) = sqrt(2.0_prec) * OM(4:6, 1:3)

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
				+sign(1.0_prec,DU)*VK(M1,3)*VH(M2,3)
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
		E=zero
		DO J=1,NR
			if(W(J)) E=E+sum((X(3*J-2:3*J)-XP(3*J-2:3*J))**2)
		end do
		E=E/2.0_prec
	end if

END function LSQSTR

      SUBROUTINE EIGEN (A,R,N,MV) 
!                                                                       
!CCCCC W.F. VAN GUNSTEREN, CAMBRIDGE, JUNE 1979 CCCCCCCCCCCCCCCCCCCCCCCC
!                                                                      C
!     SUBROUTINE EIGEN (A,R,N,MV)                                      C
!                                                                      C
!COMMENT:   EIGEN COMPUTES EIGENVALUES AND EIGENVECTORS OF THE REAL    C
!     SYMMETRIC N*N MATRIX A, USING THE DIAGONALIZATION METHOD         C
!     DESCRIBED IN "MATHEMATICAL METHODS FOR DIGITAL COMPUTERS", EDS.  C
!     A.RALSTON AND H.S.WILF, WILEY, NEW YORK, 1962, CHAPTER 7.        C
!     IT HAS BEEN COPIED FROM THE IBM SCIENTIFIC SUBROUTINE PACKAGE.   C
!                                                                      C
!     A(1..N*(N+1)/2) = MATRIX TO BE DIAGONALIZED, STORED IN SYMMETRIC C
!                       STORAGE MODE, VIZ. THE I,J-TH ELEMENT (I.GE.J) C
!                       IS STORED AT THE LOCATION K=I*(I-1)/2+J IN A;  C
!                       THE EIGENVALUES ARE DELIVERED IN DESCENDING    C
!                       ORDER ON THE DIAGONAL, VIZ. AT THE LOCATIONS   C
!                       K=I*(I+1)/2                                    C
!     R(1..N,1..N) = DELIVERED WITH THE CORRESPONDING EIGENVECTORS     C
!                    STORED COLUMNWISE                                 C
!     N = ORDER OF MATRICES A AND R                                    C
!     MV = 0 : EIGENVALUES AND EIGENVECTORS ARE COMPUTED               C
!        = 1 : ONLY EIGENVALUES ARE COMPUTED                           C
!                                                                      C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                                       
        !implicit double precision (A-H, O-Z)
        real(kind=prec) :: A,B,C,D,E,F,G,H,O,P,Q,R,S,T,U,V,W,X,Y,Z
	real(kind=prec) :: anorm, anrmx, thr, sinx, sinx2, cosx, cosx2, sincs


        !implicit integer (I-N)
        integer         :: I,J,K,L,M,N, range, mv, iq, ij, ia, ind, mq, lq, lm, ll, mm, ilq, imq, im, il, ilr, imr, jq


    DIMENSION A(1),R(1) 

!                                                                       
!*****GENERATE IDENTITY MATRIX                                          
    5 RANGE=1.E-12 
      IF (MV-1) 10,25,10 
   10 IQ=-N 
      DO 20 J=1,N 
      IQ=IQ+N 
      DO 20 I=1,N 
      IJ=IQ+I 
      R(IJ)=0.E0 
      IF (I-J) 20,15,20 
   15 R(IJ)=1.E0 
   20 CONTINUE 
!                                                                       
!*****COMPUTE INITIAL AND FINAL NORMS (ANORM AND ANRMX)                 
   25 ANORM=0.E0 
      DO 35 I=1,N 
      DO 35 J=I,N 
      IF (I-J) 30,35,30 
   30 IA=I+(J*J-J)/2 
      ANORM=ANORM+A(IA)*A(IA) 
   35 CONTINUE 
      IF (ANORM) 165,165,40 
   40 ANORM=1.414E0* SQRT(ANORM) 
      ANRMX=ANORM*RANGE/FLOAT(N) 
!                                                                       
!*****INITIALIZE INDICATORS AND COMPUTE THRESHOLD, THR                  
      IND=0 
      THR=ANORM 
   45 THR=THR/FLOAT(N) 
   50 L=1 
   55 M=L+1 
!                                                                       
!*****COMPUT SIN AND COS                                                
   60 MQ=(M*M-M)/2 
      LQ=(L*L-L)/2 
      LM=L+MQ 
   62 IF ( ABS(A(LM))-THR) 130,65,65 
   65 IND=1 
      LL=L+LQ 
      MM=M+MQ 
      X=0.5E0*(A(LL)-A(MM)) 
   68 Y=-A(LM)/ SQRT(A(LM)*A(LM)+X*X) 
      IF (X) 70,75,75 
   70 Y=-Y 
   75 SINX=Y/ SQRT(2.E0*(1.E0+( SQRT(1.E0-Y*Y)))) 
      SINX2=SINX*SINX 
   78 COSX= SQRT(1.E0-SINX2) 
      COSX2=COSX*COSX 
      SINCS=SINX*COSX 
!                                                                       
!*****ROTATE L AND M COLUMNS                                            
      ILQ=N*(L-1) 
      IMQ=N*(M-1) 
      DO 125 I=1,N 
      IQ=(I*I-I)/2 
      IF (I-L) 80,115,80 
   80 IF (I-M) 85,115,90 
   85 IM=I+MQ 
      GOTO 95 
   90 IM=M+IQ 
   95 IF (I-L) 100,105,105 
  100 IL=I+LQ 
      GOTO 110 
  105 IL=L+IQ 
  110 X=A(IL)*COSX-A(IM)*SINX 
      A(IM)=A(IL)*SINX+A(IM)*COSX 
      A(IL)=X 
  115 IF (MV-1) 120,125,120 
  120 ILR=ILQ+I 
      IMR=IMQ+I 
      X=R(ILR)*COSX-R(IMR)*SINX 
      R(IMR)=R(ILR)*SINX+R(IMR)*COSX 
      R(ILR)=X 
  125 END DO 
      X=2.E0*A(LM)*SINCS 
      Y=A(LL)*COSX2+A(MM)*SINX2-X 
      X=A(LL)*SINX2+A(MM)*COSX2+X 
      A(LM)=(A(LL)-A(MM))*SINCS+A(LM)*(COSX2-SINX2) 
      A(LL)=Y 
      A(MM)=X 
!                                                                       
!*****TESTS FOR COMPLETION                                              
!                                                                       
!*****TEST FOR M = LAST COLUMN                                          
  130 IF (M-N) 135,140,135 
  135 M=M+1 
      GOTO 60 
!                                                                       
!*****TEST FOR L = SECOND FROM LAST COLUMN                              
  140 IF (L-(N-1)) 145,150,145 
  145 L=L+1 
      GOTO 55 
  150 IF (IND-1) 160,155,160 
  155 IND=0 
      GOTO 50 
!                                                                       
!*****COMPARE THRESHOLD WITH FINAL NORM                                 
  160 IF (THR-ANRMX) 165,165,45 
!                                                                       
!*****SORT EIGENVALUES AND EIGENVECTORS                                 
  165 IQ=-N 
      DO 185 I=1,N 
      IQ=IQ+N 
      LL=I+(I*I-I)/2 
      JQ=N*(I-2) 
      DO 185 J=I,N 
      JQ=JQ+N 
      MM=J+(J*J-J)/2 
      IF (A(LL)-A(MM)) 170,185,185 
  170 X=A(LL) 
      A(LL)=A(MM) 
      A(MM)=X 
      IF (MV-1) 175,185,175 
  175 DO 180 K=1,N 
      ILR=IQ+K 
      IMR=JQ+K 
      X=R(ILR) 
      R(ILR)=R(IMR) 
  180 R(IMR)=X 
  185 CONTINUE 
!                                                                       
      RETURN 
      END subroutine eigen                                           
                                          

end module QEIGEN

