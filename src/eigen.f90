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
      SUBROUTINE EIGEN (A,R,N,MV) 

      implicit double precision (A-H, O-Z)
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
      END                                           
