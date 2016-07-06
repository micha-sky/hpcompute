C======================================================================C 
C 
	SUBROUTINE CRSSEC
C
C----------------------------------------------------------------------C
C
C     CRSSEC: Calculate CRoSs-SECtional characteristics (area,width,
C             perimeter.
C           
C----------------------------------------------------------------------C 
      USE CRSECT 
      USE DRLINK 
      USE DRNAME 
      USE DRWFLD 
      USE FILINX 
      USE FLDINX 
      USE POINTS 
      USE NUMBRS
C----------------------------------------------------------------------C
C     Implicit Double Precision.
C----------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C----------------------------------------------------------------------C 
C     Type Declarations. 
C----------------------------------------------------------------------C 
C----------------------------------------------------------------------C
C     Calculate cross-sectional area, width, wetting perimeter.
C----------------------------------------------------------------------C
      DO 400 L=1,NLINK 
        DO 300 I=1,IL(L)  
          NP = ND(L,I)  
                     
C----------------------------------------------------------------------C
C     Define numbers of first and last records at cross-sectional table
C     for selected computational point. 
C----------------------------------------------------------------------C
          NRB = NSCT(1,MSEC(NP)) 
          NRE = NSCT(2,MSEC(NP)) 
C     Surface water depth exceeds of cross-sectional depth. 
          IF( SHG(NRE) .LT. HR(NP) ) THEN  
            WRITE(*,*) 'Error:',
     &    ' Water depth exceeds cross-sectional depth (See OUTPUT file)'
            WRITE(IWR,*) 'Error:',
     &    ' Water depth exceeds cross-sectional depth (See OUTPUT file)'
            WRITE (IWR,9008) I,LinkName(L) 
            STOP 
          ENDIF       
C             
          AW(NP) = Small; DAW(NP) = Small;
          PER(NP) = Small; DPER(NP) = Small; 
            
          Do j=NRB+1,NRE   
            IF( HR(NP) >= SHG(j-1) .and. HR(NP) <= SHG(j) ) THEN 
              DHr   = SHG(j) -SHG(j-1)   
              DelH = HR(NP) -SHG(j-1) 
              
              BW(NP) = SWD(j)*DelH/DHr +SWD(j-1)*(1.d0 -DelH/DHr)
              AW(NP) = DelH*(SWD(j-1) +0.5d0*(SWD(j)-SWD(j-1))*DelH/DHr)  
              AW(NP) = SAr(J-1) +AW(NP) 
              PER(NP) = DelH*(SPER(j) -SPER(j-1))/DHr +SPER(j-1)
              DAW(NP) = SWD(j-1) +(SWD(j) -SWD(j-1))*DelH/DHr 
              DPER(NP) = (SPER(j) -SPER(j-1))/DHr     
              GOTO 300   
            ENDIF   
          EndDo

  300   CONTINUE 
  400 CONTINUE         
C----------------------------------------------------------------------C 
C     Format Statements. 
C----------------------------------------------------------------------C 
 9008 FORMAT(/' Input or Compilation Error: Water depth exceeds ', 
     &        'cross-sectional depth for point No ',I6,
     &        ' on the link ',A, '.') 
C----------------------------------------------------------------------C 
C     End of CRSSEC group. 
C----------------------------------------------------------------------C 
C 
      RETURN 
      END 
 
C======================================================================C
C
      SUBROUTINE LININT( X,Y,DYDX,XTABL,YTABL,NTABL )
C
C----------------------------------------------------------------------C
C
C     LININT: LINear INTerpolation.
C
C     Subroutine computes Y and DYDX from a table of XTABL-YTABL values
C     with linear interpolation.  The XTABL-YTABL may be in ascending or
C     descending order.
C
C     MSTS (Multiphase Subsurface Transport Simulator) code. 
C     M.D. White, and W.E. Nichols. Pacific Northwest National Laboratory, 
C     Richland, WA, 1993.
C----------------------------------------------------------------------C 
      USE FILINX 
C----------------------------------------------------------------------C
C     Implicit Double Precision.
C----------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C----------------------------------------------------------------------C
C     Type Declarations.
C----------------------------------------------------------------------C
      REAL*8 XTABL(*),YTABL(*)
C----------------------------------------------------------------------C
C     Determine table order.
C----------------------------------------------------------------------C
      IF( XTABL(1) .LT. XTABL(NTABL) ) THEN
        ISTART = 2
        ISTOP  = NTABL
        INCR   = 1
C----------------------------------------------------------------------C
C     Limit Y to the table range.
C----------------------------------------------------------------------C
        XDUM = MIN( X,XTABL(NTABL) )
        XDUM = MAX( XDUM,XTABL(1) )
      ELSE
        ISTART = NTABL -1
        ISTOP  = 1
        INCR   = -1
C----------------------------------------------------------------------C
C     Limit Y to the table range.
C----------------------------------------------------------------------C
        XDUM = MIN( X,XTABL(1) )
        XDUM = MAX( XDUM,XTABL(NTABL) )
      ENDIF
C----------------------------------------------------------------------C
C     Loop through table to find bounding paramete.
C----------------------------------------------------------------------C 
      DO 200 I = ISTART,ISTOP,INCR
        IF( XDUM .LE. XTABL(I) .AND. XDUM .GE. XTABL(I-INCR) ) THEN 
          DYDX = (YTABL(I)-YTABL(I-INCR))/(XTABL(I)-XTABL(I-INCR))
          Y = (XDUM-XTABL(I-INCR))*DYDX +YTABL(I-INCR)
          GOTO 400
        ENDIF
  200 CONTINUE
      WRITE (IWR,'(2A)') 'NTABL = ',NTABL 
      WRITE (IWR,'(2A)') 'XTABL(1) = ',XTABL(1)
      WRITE (IWR,'(2A)') 'XTABL(NTABL)',XTABL(NTABL) 
      WRITE (IWR,'(2A)') 'YTABL(1) = ',YTABL(1)
      WRITE (IWR,'(2A)') 'YTABL(NTABL)',YTABL(NTABL) 
      WRITE (IWR,'(2A)') 'Table Error: A non-monotonically increasing ', 
     &  'or decreasing table has been encountered.'
      WRITE (*,'(/A)') 'FATAL TABLE ERROR: Check output file.' 
      STOP
  400 CONTINUE
C----------------------------------------------------------------------C  
C     End of LININT group.
C----------------------------------------------------------------------C
C
      RETURN
      END

C======================================================================C
C
      SUBROUTINE LQPREV
C
C----------------------------------------------------------------------C
C
C     LQPREV: save LiQuid phase variables at PREVious time step..
C
C----------------------------------------------------------------------C 
      USE DRWFLD  
      USE DRSPFL
      USE FLDINX  
      USE LANDSF
      USE POINTS 
      USE SOLVAR
C----------------------------------------------------------------------C
C     Implicit Double Precision.
C----------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C----------------------------------------------------------------------C
C     Assign the primary variables at previous time step.
C----------------------------------------------------------------------C
      DO 100 N  = 1,NPOINT
        HRO(N)  = HR(N)
        QHO(N)  = QH(N)
        AWO(N)  = AW(N) 
	  PERO(N) = PER(N)    
  100 CONTINUE  
      
      IF( ISOLVE(2) == 0 ) RETURN 
      DO N=1, NPOINT
        SHO(N)  = SH(N) 
        TSFO(N) = TSF(N) 
        ZSDO(N) = ZSD(N)   
      ENDDO  
C----------------------------------------------------------------------C
C     End of LQPREV group.
C----------------------------------------------------------------------C
C
      RETURN
      END

C======================================================================C
          