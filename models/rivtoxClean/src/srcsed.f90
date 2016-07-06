      SUBROUTINE SRCSED
!
!------------------------------------------------------------------------------C
!
!     SRCSED: SouRCe of SEDimnet for transport equation.
!
!     HR - channel water depth,
!     QH - channel discharge,
!     AR - cross-sectional area of channel flow.
!------------------------------------------------------------------------------C
!     Used Modules.                                                
!------------------------------------------------------------------------------C
!  USE CONSTS 
!  USE CRSECT   
  USE DRCATC 
  USE DREVAP
  USE DRLINK 
  USE DRNODE
  USE DRRAIN 
  USE DRSPFL
  USE DRWFLD   
  USE FLDINX
  USE LINEQN  
  USE NUMBRS  
  USE PLZONE 
  USE POINTS
!  USE RESIDL  
  USE SOLVAR  
!------------------------------------------------------------------------------C
!     Implicit Double Precision.
!------------------------------------------------------------------------------C
  IMPLICIT REAL*8 (A-H,O-Z)
  IMPLICIT INTEGER*4 (I-N)
!------------------------------------------------------------------------------C
!     Type Declarations.
!------------------------------------------------------------------------------C
!------------------------------------------------------------------------------C 
!------------------------------------------------------------------------------C
!     Compute elements for the species equation.
!------------------------------------------------------------------------------C
  DO 500 L=1,NLINK 
    DO 400 I=1,IL(L)  
      N  = ND(L,I) 
      IPLANT = IZPL(N)
      IZRN = IZ(N)

      IF( I == 1 ) THEN 
        NE  = ND(L,I+1) 
        VOL = 0.5d0*(X(N+1) -X(N)) 
      ELSEIF( I == IL(L) ) THEN             
        NW  = ND(L,I-1) 
        VOL = 0.5d0*(X(N) -X(N-1)) 
      ELSE 
        NW  = ND(L,I-1) 
        NE  = ND(L,I+1) 
        VOL = 0.5d0*(X(N+1) -X(N-1)) 
      ENDIF  

      HP = SIGMA*HR(N) + (1.D+0 -SIGMA)*HRO(N)
      QP = SIGMA*QH(N) + (1.D+0 -SIGMA)*QHO(N)
      AP = SIGMA*AW(N) + (1.D+0 -SIGMA)*AWO(N)
!----------------------------------------------------------------------!
!     Source term.
!----------------------------------------------------------------------!
     if( NSRAIN +NSEVAP > 0 ) then
       FF = VOL*(RAIN(N) -EVAP(N))
     else
       FF=0.d0
     endif
     FF = max(-AWO(N)/DT,FF)
     if( HR(N) <= HSMALL ) FF = max(0.d0,FF)
     NEL = iEq(N)
     ALC(NEL) = ALC(NEL) + FF

     SRC = 0.d0
!      Sediment transport.
     if( ILOAD > 0 ) then !--Skip if no Transport Capacity Equation is used
       SP  = SIGMA*SH(N)  + (1.d0 -SIGMA)*SHO(N)
       if( ILOAD /= 9 ) then
         if( SP >= SEQ(N) ) then
           SRC = VOL *VSD
         else
           tErod = EROD(IZPLN); 
           SRC = VOL *VSD *tErod
         endif
         ALC(NEL) = ALC(NEL) +SRC
         SRC = SEQ(N) *SRC
       else
!        Cohesive sediments
!        Tay_e = 2.5 N/m^2 (surficial material) and 6-7 N/m^2 (granular,
!        consolidated deposits) (for erosion)
!        Tay_d = 0.1 N/m^2 (for deposition)
!        Tay_e = 0.79*Tay_y**0.94 according to Mitura(1989) Tay_y - yield stress
!           TAY_E = 2.0D+0
!           WCM  = 1.7D-5
!           TAY_E = 2.5D+0 !--Now defined in input
!           TAY_D = 0.1D+0 !--Now defined in input
         UP = sqrt(2.d0)*HP*QP/sqrt(HP**4 +max(HP**4,small))
         TAY = UP*UP*RHOLQ
         if( TAY < TAY_D ) ALC(NEL) = ALC(NEL) +VOL*VSD*(1.d0 -TAY/TAY_D) !--Deposition
         if( TAY > TAY_E ) SRC = VOL *WCM *(TAY/TAY_E -1.d0) !--Erosion
       endif
     endif
!----------------------------------------------------------------------!
!     Solution vector.
!----------------------------------------------------------------------!
     BLC(N) = BLC(N) +SRC
400  CONTINUE 
500 CONTINUE    
!----------------------------------------------------------------------C 
!     End of SRCSED group. 
!----------------------------------------------------------------------C 
! 
    RETURN
    END
!======================================================================C 
 
!======================================================================C 
! 
      SUBROUTINE EXNER  
! 
!----------------------------------------------------------------------C 
! 
!     EXNER: land surface erosion/deposition equation (EXNER equation). 
! 
!----------------------------------------------------------------------C  
      USE CONSTS 
      USE DRLINK  
      USE DRSEDT 
      USE DRWFLD 
      USE DRSPFL
      USE FLDINX
      USE LANDSF
      USE NUMBRS 
      USE PLZONE  
      USE PROPER 
      USE SOLVAR 
!----------------------------------------------------------------------C 
!     Implicit Double Precision. 
!----------------------------------------------------------------------C 
      IMPLICIT REAL*8 (A-H,O-Z) 
      IMPLICIT INTEGER*4 (I-N) 
!----------------------------------------------------------------------C 
!     Compute thickness of upper bottom deposit layer . 
!----------------------------------------------------------------------C 
      DO 200 L = 1,NLINK 
        DO 100 I = 1,IL(L) 
          N = ND(L,I) 
          IZN = IZ(N) 
!----------------------------------------------------------------------C 
!     Thickness of upper bottom deposit. 
!----------------------------------------------------------------------C 
          QS = SDOWN(N)  
          QB = RSUP(N)  
!          ttws = TSF(N)+HR(N)
          dz = DT*(QS - QB)/(1.D+0 -POR(IZN))/RHOS(IZN)
          ZSD(N) = MAX(ZSDO(N) + dz,ZERO)
          TSF(N) = TSFO(N) + dz

!          IF( ISOLVE(1) == 1 ) THEN
!            HR(N)=max(ttws-TSF(N),ZERO)
!            if (HR(N)>=SMALL) then
!             UH(N)=QX(N)/HR(N)
!            else
!             QX(N)=0.0d+0
!             QY(N)=0.0d+0
!            endif 
!          ENDIF 

  100   CONTINUE 
  200 CONTINUE 
!----------------------------------------------------------------------C 
!     End of EXNER group. 
!----------------------------------------------------------------------! 
! 
      RETURN 
      END 
 
!======================================================================C 
!
      FUNCTION SDOWN ( N )
!
!----------------------------------------------------------------------C
!
!     SDOWN: Sedimentation rate in DOWN direction.
!
!----------------------------------------------------------------------C
      USE CONSTS 
      USE DRLINK  
      USE DRSEDT 
      USE DRWFLD 
      USE DRSPFL
      USE FLDINX
      USE LANDSF
      USE NUMBRS 
      USE PLZONE  
      USE PROPER 
      USE SOLVAR 
!----------------------------------------------------------------------C
!     Implicit Double Precision.
!----------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
!----------------------------------------------------------------------C 
!     The overland flow depth is equal to zero? 
!----------------------------------------------------------------------C
      HP  = SIGMA*HR(N)  + (1.D+0 -SIGMA)*HRO(N) 
      SDOWN = 0.D+0 
      IF( HP <= HSMALL ) THEN
        RETURN 
      ENDIF      
!----------------------------------------------------------------------C
      IZN = IZ(N) 
      SP  = SIGMA*SH(N)  + (1.D+0 -SIGMA)*SHO(N) 
      IF( ILOAD .NE. 9 ) THEN 
        SDOWN = MAX( 0.d+0,VSD*(SP-SEQ(N)) ) 
      ELSE 
!     Cohesive sediments
!     Tay_d = 0.1 N/m^2 (for deposition)  
!        TAY_D = 0.1D+0            
        VISKN = VISLQ/RHOLQ 
        QP = SIGMA*QH(N)  + (1.D+0 -SIGMA)*QHO(N) 
        UP = sqrt(2.d0)*HP*QP/sqrt(HP**4 +max(HP**4,small))
!        US = abs(UP) 
!        WG  = DSE*DSE*GRAV/18.D+0*(RHOSE-RHOLQ)/VISLQ 
!!        UCR = 0.06D+0*GRAV*(RHOSE/RHOLQ -1.D+0)*SQRT(DSE*VISKN) 
!        UCR = SQRT(TAY_D/RHOLQ)  
!        UR  = US/MAX(UCR,SMALL) 
!        SDOWN = MAX( 0.D+0,WG*SP*(1.D+0 - UR*UR) )  
        
        TAY = UP*UP *RHOLQ
        if( TAY < TAY_D ) SDOWN = VSD *SP *(1.d0 -TAY/TAY_D)
      ENDIF      
!----------------------------------------------------------------------C
!     End of SDOWN group.
!----------------------------------------------------------------------C
!
      RETURN
      END

!======================================================================C
!
      FUNCTION RSUP ( N )
!
!----------------------------------------------------------------------C
!
!     RSUP: ReSuspension rate in UP direction.
!
!----------------------------------------------------------------------C
      USE CONSTS  
      USE DRCATC
      USE DRLINK  
      USE DRSEDT 
      USE DRWFLD 
      USE DRSPFL
      USE FLDINX
      USE LANDSF
      USE NUMBRS 
      USE PLZONE  
      USE PROPER 
      USE SOLVAR 
!----------------------------------------------------------------------C
!     Implicit Double Precision.
!----------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
!----------------------------------------------------------------------C 
!     The overland flow depth is equal to zero? 
!----------------------------------------------------------------------C
      HP  = SIGMA*HR(N)  + (1.D+0 -SIGMA)*HRO(N) 
      RSUP = 0.D+0 
      IF( HP .LE. HSMALL ) THEN
        RETURN 
      ENDIF      
!----------------------------------------------------------------------C
      IZN = IZ(N)
      IZPLN = IZPL(N) 
      SP  = SIGMA*SH(N)  + (1.D+0 -SIGMA)*SHO(N)
      IF( ILOAD .NE. 9 ) THEN 
        RSUP = MAX( 0.d+0,EROD(IZPLN)*VSD*(SEQ(N)-SP) ) 
      ELSE 
        VISKN = VISLQ/RHOLQ 
        QP = SIGMA*QH(N)  + (1.d+0 -SIGMA)*QHO(N) 
        UP = sqrt(2.d0)*HP*QP/sqrt(HP**4 +max(HP**4,small))
!        US = abs(UP) 
!!        UCR = 0.06D+0*GRAV*(RHOSE/RHOLQ -1.D+0)*SQRT(DSE*VISKN) 
!        UCR = SQRT(TAY_E/RHOLQ)  
!        UR  = US/MAX(UCR,SMALL) 
!        RSUP = MAX( 0.D+0,WCM*(UR*UR - 1.d+0) )  
        
        TAY = UP*UP *RHOLQ
        if( TAY > TAY_E ) RSUP = WCM*(TAY/TAY_E -1.d0)
      ENDIF      
!----------------------------------------------------------------------C
!     End of RSUP group.
!----------------------------------------------------------------------C
!
      RETURN
      END

!======================================================================C
    
    
   
   


