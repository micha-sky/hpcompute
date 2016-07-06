C======================================================================C 
C 
	SUBROUTINE MATRIX  
C
C----------------------------------------------------------------------C
C
C     MATRIX: calculate MATRIX of linear equations system.
C              
C----------------------------------------------------------------------C 
      USE DRLINK 
      USE DRNODE 
      USE DRWFLD  
      USE FLDINX  
      USE LANDSF
      USE LINEQN 
      USE POINTS
C----------------------------------------------------------------------C
C     Implicit Double Precision.
C----------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C----------------------------------------------------------------------C 
C     Type Declarations. 
C----------------------------------------------------------------------C 
      INTEGER*4, POINTER :: KIN(:) => null()
      INTEGER*4, POINTER :: KOUT(:) => null() ! KIN(LMAXL),KOUT(LMAXL)    
C----------------------------------------------------------------------C
C     Initialize the Jacobian matrix.
C----------------------------------------------------------------------C
      DO 40 I = 1,NNODE
        DO 30 J = 1,NNODE
          ALU(I,J) = 0.D+0
   30   CONTINUE
        BLU(I) = 0.D+0
   40 CONTINUE   
C----------------------------------------------------------------------C 
      if( .not.associated(KIN) ) then  
        allocate(KIN(MaxLinks),KOUT(MaxLinks)) 
        KIN = 0.d0; KOUT = 0.d0;
      endif    
C----------------------------------------------------------------------C
C     Compute elements for the Jacobian matrix.
C----------------------------------------------------------------------C
      DO 600 I=1,NNODE   
C     Compute number of outflow and inflow links.  
        KKK = KK(I) 
        NOUT = 0  
        NIN = 0 
        DO K=1,KKK 
          IF( LK(I,K) .LT. 0 ) THEN 
            NIN = NIN +1 
            KIN(NIN) = -LK(I,K) 
          ELSE 
            NOUT = NOUT +1 
            KOUT(NOUT) = LK(I,K)  
          ENDIF   
        ENDDO                      
C      
!        IF( NOUT .NE. 0 ) THEN  
C          IF( NIN .NE. 0 ) THEN  
          SUM_IN = 0.D+0   
          S_IN = 0.D+0 
          DO K=1,NIN    
            KI = KIN(K) 
            KO = MU(KI) 
            NP = ND(KI,1)  
            ALU(I,KO) = CRH(NP,4) !inflow links 
            SUM_IN = SUM_IN +CRH(NP,6) 
            S_IN = S_IN +CRH(NP,5) +QH(ND(KI,IL(KI)))   
C            S_IN = S_IN +CRH(NP,5) +Q(2,KI,IL(KI))   
          ENDDO        
C          ENDIF  
C           
          SUM_OUT = 0.D+0  
          S_OUT = 0.D+0 
          DO K=1,NOUT     
            KO = KOUT(K) 
            KI = MD(KO) 
            NP = ND(KO,1) 
            ALU(I,KI) = -CRH(NP,3) !outflow links 
            SUM_OUT = SUM_OUT +CRH(NP,1) 
            S_OUT = S_OUT +CRH(NP,2) +QH(NP)  
C            S_OUT = S_OUT +CRH(NP,2) +Q(2,KO,1)  
          ENDDO 
C        
          ALU(I,I) = SUM_IN -SUM_OUT 
          BLU(I) = -S_IN +S_OUT               
C          BLU(I) = -QGRNO(i) -S_IN +S_OUT               
!        ELSE  
!          ALU(I,I) = 1.D+0 
!          KI = KIN(1) 
!          BLU(I) = -HR(ND(KI,IL(KI))) -TSF(ND(KI,IL(KI)))     
C          BLU(I) = YGRNO(i) -HR(ND(KI,IL(KI))) -TSF(ND(KI,IL(KI)))     
C          BLU(I) = YGRNO(i) -Y(2,KI,IL(KI))   
!        ENDIF             
  600 CONTINUE       
C----------------------------------------------------------------------C 
C     End of MATRIX group. 
C----------------------------------------------------------------------C 
C 
      RETURN 
      END 
 
C======================================================================C 
          