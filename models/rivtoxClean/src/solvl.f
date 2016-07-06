C======================================================================C 
C 
	SUBROUTINE SOLVL   
C
C----------------------------------------------------------------------C
C
C     SOLVL: SOLVer for system of Linear equations.
C           
C----------------------------------------------------------------------C  
      USE FLDINX
      USE DRNODE 
      USE LNDAMS 
      USE LINEQN
C----------------------------------------------------------------------C
C     Implicit Double Precision.
C----------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C----------------------------------------------------------------------C 
C     Type Declarations. 
C----------------------------------------------------------------------C 
C----------------------------------------------------------------------C 
      if( .not.associated(DY) ) then  
        allocate(DY(NNODE)); DY = 0.d0; 
      endif    
C----------------------------------------------------------------------C 
      NDAM = NGROUP 
      NY = NGSIZE(1)
      DO I=1,NY-1
        IF( ALU(I+1,I) .NE. 0 ) THEN
          alfa = ALU(I+1,I)/ALU(I,I)
          DO J=I,NY
            ALU(I+1,J) = ALU(I+1,J) -alfa*ALU(I,J)
          ENDDO    
        ENDIF
      ENDDO   
C----------------------------------------------------------------------C
      DY(NY) = BLU(NY)/ALU(NY,NY)
      DO I=1,NY-1
        J = NY -I
        SUM = 0.D+0
        DO K=J+1,NY
          SUM = SUM +ALU(J,K)*DY(K)
        ENDDO      
        DY(J) = (BLU(J) -SUM)/ALU(J,J)
      ENDDO 
C----------------------------------------------------------------------C
      IF( NDAM .NE. 0 ) THEN
        DO knode=1,NDAM
          NY1 = NGSIZE(knode)
          NY2 = NGSIZE(knode+1)
          DO I=NY1+1,NY2-1
            IF( ALU(I+1,I) .NE. 0 ) THEN
              alfa = ALU(I+1,I)/ALU(I,I)
              DO J=I,NY2
                ALU(I+1,J) = ALU(I+1,J) -alfa*ALU(I,J)
              ENDDO ! po j
            ENDIF
		ENDDO ! po i
        ENDDO  ! po knode
      ENDIF
C----------------------------------------------------------------------C
      IF( NDAM .NE. 0 ) THEN
        DO knode=1,NDAM
          NY1 = NGSIZE(knode)
          NY2 = NGSIZE(knode+1)
          DY(NY2) = BLU(NY2)/ALU(NY2,NY2)
          DO I=NY1+1,NY2-1
            J = NY2 +NY1 -I
            SUM = 0.D+0
            DO K=J+1,NY2
              SUM = SUM +ALU(J,K)*DY(K)
            ENDDO  ! po k
            DY(J) = (BLU(J) -SUM)/ALU(J,J)
          ENDDO ! po i
        ENDDO ! po knode
      ENDIF
C----------------------------------------------------------------------C 
C     End of SOLVL group. 
C----------------------------------------------------------------------C 
C 
      RETURN 
      END 
 
C======================================================================C 
