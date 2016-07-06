!==============================================================================C 
!
      SUBROUTINE SPEC (C,CO)
!
!------------------------------------------------------------------------------C
!
!     SPEC: SPECies transport equation..
!
!     HR - channel water depth,
!     QH - channel discharge,
!     AR - cross-sectional area of channel flow.
!------------------------------------------------------------------------------C
!     Used Modules.                                                
!------------------------------------------------------------------------------C
!      USE CONSTS 
!      USE CRSECT   
      USE DRCATC
      USE DRLINK 
      USE DRNODE
      USE DRWFLD   
      USE FLDINX
      USE LINEQN  
      USE NUMBRS  
      USE PLZONE 
      USE POINTS
!      USE RESIDL  
      USE SOLVAR  
      USE DREVAP  
      USE DRRAIN
!------------------------------------------------------------------------------C
!     Implicit Double Precision.
!------------------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
!------------------------------------------------------------------------------C
!     Type Declarations.
!------------------------------------------------------------------------------C
      REAL*8 C(*),CO(*)
!------------------------------------------------------------------------------C 
      iDf = 3
      IERR = 0  
      NEL = 0 
      ALC = 0.d0; BLC = 0.d0; IA = 0; JA = 0;
!------------------------------------------------------------------------------C
!     Compute elements for the species equation.
!------------------------------------------------------------------------------C
      DO 500 L=1,NLINK 
        DO 400 I=1,IL(L)  
          N  = ND(L,I) 
!          NX = IX(L) +I -1 
          IZN = IZ(N)

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

!          UP = SIGMA*UH(N) + (1.D+0 -SIGMA)*UHO(N)
          QP = SIGMA*QH(N) + (1.D+0 -SIGMA)*QHO(N)
          AP = SIGMA*AW(N) + (1.D+0 -SIGMA)*AWO(N)
!------------------------------------------------------------------------------C
!     Storage terms.
!------------------------------------------------------------------------------C  
          SC = VOL/DT/SIGMA   
          NEL = iEq(N)
          ALC(NEL) = SC*AWO(N)  
          IA(NEL) = N  
          JA(NEL) = N
!------------------------------------------------------------------------------C
!     West side diffusion terms.
!------------------------------------------------------------------------------C 
          IF( I /= 1 ) THEN 
            IZNW = IZ(NW) 
            QW = SIGMA*QH(NW) + (1.D+0 -SIGMA)*QHO(NW)
!            UW = SIGMA*UH(NW) + (1.D+0 -SIGMA)*UHO(NW)
            QB = 0.5D+0 * (QW + QP)
            AWs = SIGMA*AW(NW) + (1.D+0 -SIGMA)*AWO(NW)
            AB = 0.5D+0 * (AWs + AP)

!     Diffusion terms 
!            DISW = APRXM( DFSL(IZNW)*AWs,DFSL(IZN)*AP,X(N)-X(N-1),            &
!                                 X(N+1)-X(N),IDMEAN(iDf) ) 
            DISW = APRXM( DFSL(IZNW)*AWs,DFSL(IZN)*AP,1.d0,1.d0,IDMEAN(iDf) ) 
            DLX = DISW/(X(N) -X(N-1)) 
! 
            DFP  = (DLX - MIN(QB,ZERO)) 
            DFW  = (DLX + MAX(QB,ZERO))   
!----------------------------------------------------------------------C 
!     Modify Matrix of Linear Equations. 
!----------------------------------------------------------------------C 
            NEL = iEq(N)
            ALC(NEL) = ALC(NEL) + DFW 
            NEL = iEq(N) +1
            ALC(NEL) = ALC(NEL) - DFW  
            IA(NEL) = N 
            JA(NEL) = NW  
          ELSE   
            inode = MU(L) 
            if( QH(n) > 0.d0 .and. KK(inode) > 1 ) then   
              KKK = KK(inode) 
              Qsum = 0.d0; Qout = 0.d0; 
              do j=1,KKK  
                Lin = -LK(inode,j) 
                if( Lin > 0 ) then
                  nIn = nd(Lin,IL(Lin))  
                else  
                  nIn = nd(-Lin,1)  
                endif   
!                
                if( Lin*QH(nIn) > 0 ) then  
                  Qsum = Qsum +abs(QH(nIn))   
                else  
                  Qout = Qout +abs(QH(nIn))  
                endif    
              enddo  
!              
              NEL = iEq(N)
              ALC(NEL) = ALC(NEL) + QH(n) 
              iew = 3;
              do j=1,KKK 
                Lin = -LK(inode,j) 
                if( Lin > 0 ) then
                  nIn = nd(Lin,IL(Lin))  
                else  
                  nIn = nd(-Lin,1)  
                endif  
!                
                if( Lin*QH(nIn) > 0 ) then  
                  DFW = abs(QH(nIn))/Qout *QH(n)   
!     Modify matrix of linear equations.                  
                  NEL = iEq(N) +iew
                  ALC(NEL) = ALC(NEL) - DFW  
                  IA(NEL) = N 
                  JA(NEL) = nIn 
                  iew = iew +1
                endif  
              enddo  
            endif
          ENDIF  
!------------------------------------------------------------------------------C
!     East side diffusion terms.
!------------------------------------------------------------------------------C 
          IF( I /= IL(L) ) THEN 
            IZNE = IZ(NE) 
!            
            QE = SIGMA*QH(NE) + (1.D+0 -SIGMA)*QHO(NE)
!            UE = SIGMA*UH(NE) + (1.D+0 -SIGMA)*UHO(NE)
            QB = 0.5D+0 * (QE + QP)
            AE = SIGMA*AW(NE) + (1.D+0 -SIGMA)*AWO(NE)
            AB = 0.5D+0 * (AE + AP)
!             
!     Diffusion terms 
!            DISE = APRXM( DFSL(IZN)*AP,DFSL(IZNE)*AE,X(N)-X(N-1),             &
!                           X(N+1)-X(N),IDMEAN(iDf) ) 
            DISE = APRXM( DFSL(IZN)*AP,DFSL(IZNE)*AE,1.d0,1.d0,IDMEAN(iDf) ) 
            DLX = DISE/(X(N+1) -X(N)) 
! 
            DFP  = (DLX + MAX(QB,ZERO)) 
            DFE  = (DLX - MIN(QB,ZERO)) 
!----------------------------------------------------------------------C 
!     Modify Matrix of Linear Equations. 
!----------------------------------------------------------------------C 
            NEL = iEq(N)
            ALC(NEL) = ALC(NEL) + DFE 
            NEL = iEq(N) +2
            ALC(NEL) = ALC(NEL) - DFE  
            IA(NEL) = N  
            JA(NEL) = NE 
          ELSE  
            inode = MD(L) 
            if( QH(n) < 0.d0 .and. KK(inode) > 1 ) then   
              KKK = KK(inode) 
              Qsum = 0.d0; Qout = 0.d0; 
              do j=1,KKK  
                Lin = LK(inode,j) 
                if( Lin < 0 ) then
                  nIn = nd(-Lin,IL(-Lin))  
                else  
                  nIn = nd(Lin,1)  
                endif  
!                
                if( Lin*QH(nIn) < 0 ) then  
                  Qsum = Qsum +abs(QH(nIn))   
                else   
                  Qout = Qout +abs(QH(nIn))  
                endif    
              enddo  
!              
              NEL = iEq(N)
              ALC(NEL) = ALC(NEL) -QH(n) 
              iew = 3;
              do j=1,KKK 
                Lin = LK(inode,j) 
                if( Lin < 0 ) then
                  nIn = nd(-Lin,IL(-Lin))  
                else  
                  nIn = nd(Lin,1)  
                endif  
!                
                if( Lin*QH(nIn) < 0 ) then  
                  DFW = abs(QH(nIn))/Qout *QH(n)    
!     Modify matrix of linear equations.                  
                  NEL = iEq(N) +iew
                  ALC(NEL) = ALC(NEL) +DFW  
                  IA(NEL) = N 
                  JA(NEL) = nIn 
                  iew = iew +1
                endif  
              enddo    
            endif
          ENDIF 
!----------------------------------------------------------------------!
!     Source term.
!----------------------------------------------------------------------!
     if( NSPREC +NSEVAP > 0) then
       FF = VOL*(RAIN(N) -EVAP(N))
     else
       FF=0.d0
     endif
     FF = max(-AWO(N)/DT,FF)
     if (HR(N) <= 0.) FF = max(0.d0,FF)
     NEL = iEq(N) 
     ALC(NEL) = ALC(NEL) + FF
!------------------------------------------------------------------------------C
!     Water equation functions.
!------------------------------------------------------------------------------C
          NEL = iEq(N) 
!          ALC(NEL)      = ALC(NEL) +FF(3)
!----------------------------------------------------------------------C 
!     Solution vector. 
!----------------------------------------------------------------------C 
          BLC(N) = BLC(N) +CO(N)*SC*AWO(N)   
  400   CONTINUE 
  500 CONTINUE 
!----------------------------------------------------------------------C 
!     End of SPEC group. 
!----------------------------------------------------------------------C 
! 
      RETURN 
      END 
 
!==============================================================================C 
!
      SUBROUTINE BCSEDM
!
!------------------------------------------------------------------------------C
!
!     BCSEDM: Boundary Condition for the SEDiMent transport equation.
!
!------------------------------------------------------------------------------C
!     Used Modules.                                                
!------------------------------------------------------------------------------C
      USE CONSTS 
!      USE CRSECT  
      USE DRLINK 
      USE DRBCSE  
      USE DRNODE
      USE DRWFLD  
      USE FILINX
      USE FLDINX
      USE NUMBRS  
      USE LINEQN  
!      USE RESIDL  
      USE POINTS
      USE SOLVAR  
!------------------------------------------------------------------------------C
!     Implicit Double Precision.
!------------------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
!------------------------------------------------------------------------------C
!     Type Declarations.
!------------------------------------------------------------------------------C
      REAL*8, POINTER :: VBC(:) => null()
!------------------------------------------------------------------------------C
!     Loop over boundary value tables.
!------------------------------------------------------------------------------C  
      IF( .NOT.ASSOCIATED(VBC) ) THEN  
        ALLOCATE(VBC(NBTS));   
      ENDIF  
      VBC = 0.d0;
      DO 100 NS = 1,NBTS
        NSB = NBRNS(1,NS)
        NSE = NBRNS(2,NS)
!------------------------------------------------------------------------------C
!     Constant source.
!------------------------------------------------------------------------------C
        IF( NSB .EQ. NSE ) THEN
          VBC(NS) = BCVS(NSB)
!------------------------------------------------------------------------------C
!     Tabular source.
!------------------------------------------------------------------------------C
        ELSE
          VBC(NS) = 0.D+0
          IF( TIME .GT. BCTS(NSB) .AND. TIME-DT .LT. BCTS(NSE)) THEN
            DO 210 I = NSB+1,NSE
              IF( TIME .GT. BCTS(I-1) .AND. TIME-DT .LT. BCTS(I)) THEN
!
                TMIN = MAX( TIME-DT, BCTS(I-1) )
                TMAX = MIN( TIME, BCTS(I) )
                DTBND = TMAX -TMIN
                TMID = 0.5*(TMIN +TMAX)
                NLEN = 2
!                CALL LININT( TIME,VALUE,ZZ,BCT(I-1),BCV(I-1),NLEN)
!                CALL LININT( TMID,VALUE,ZZ,BCTS(I-1),BCVS(I-1),NLEN) 
                VALUE = BCVS(I-1) 
                VBC(NS) = VBC(NS) + VALUE
!                VBC(NS) = VBC(NS) + VALUE*DTBND/DT
              ENDIF  
  210       CONTINUE    
          ENDIF
        ENDIF
100   CONTINUE   
!------------------------------------------------------------------------------C
!     Sediment Boundary Conditions.
!------------------------------------------------------------------------------C
      DO 300 NBC = 1,NBCS
        N   = IBCNS(NBC)      ! Node number.
        IBC = IBCTS(NBC)  
        BV  = VBC(MBCS(NBC))  ! Boundary field value. 
!        IZN = IZR(N)
!        IZPLN = IZPL(N)
   
        KKK = KK(N) 
        if( KKK > 1 ) then   
          write(*,'(a,i5,a)') ' Node ',N,' is not a boundary node'  
          write(IWR,'(a,i5,a)') ' Node ',N,' is not a boundary node'  
          stop
        endif    
!------------------------------------------------------------------------------C
!     Local grid node indices.
!------------------------------------------------------------------------------C
        L = LK(N,1) 
        if( L > 0 ) then
          i = ND(L,1); ! Number of computational node (first in link). 
        else  
          i = ND(-L,IL(-L)); ! Number of computational node (last in link).
        endif   
          
        QP = SIGMA*QH(i) + (1.D+0 -SIGMA)*QHO(i)
!        AP = SIGMA*AW(i) + (1.D+0 -SIGMA)*AWO(i)  
!------------------------------------------------------------------------------C
!     Boundary direction controller.
!------------------------------------------------------------------------------C 
        AB = 0.d0; AP = 0.d0;
        IF( L > 0 ) GOTO 240
        IF( L < 0 ) GOTO 220  
        GOTO 230
!------------------------------------------------------------------------------C
!    Left Link Boundary.
!------------------------------------------------------------------------------C
240       CONTINUE
!     Diffusion terms
!   DLP   = tDFSL*DXIW*HP 
        DLP = 0.d0
        if (IBC == 3) then
           DLP = 0.d0
        endif 

        AB = DLP + max(QP,ZERO)
        AP = DLP - min(QP,ZERO)
        GOTO 230
!------------------------------------------------------------------------------C
!     Right Link Boundary.
!------------------------------------------------------------------------------C
220     CONTINUE
!     Diffusion terms
!   DLP   = tDFSL*DXIE*HP 
        DLP = 0.d0
        if (IBC == 3) then
          DLP = 0.d0
        endif

        AB = DLP - min(QP,ZERO)
        AP = DLP + max(QP,ZERO)
230     CONTINUE
!------------------------------------------------------------------------------C
!     Modify the Equation Matrix and solution vector according to the
!     boundary conditions.
!------------------------------------------------------------------------------C 
        NEL = iEq(i) 
        if ((IBC == 1) .or. (IBC == 3)) then
          ALC(NEL) = ALC(NEL) + AB
          BLC(i)   = BLC(i) + AB*BV
        elseif (IBC == 2) then
          BLC(i) = BLC(i) - BV
        endif
  300 CONTINUE
!------------------------------------------------------------------------------C
!     End of BCSEDM group.
!------------------------------------------------------------------------------C
!
      RETURN
      END

!==============================================================================C 
!==============================================================================C 
!
      SUBROUTINE NOZERO( NVAR,NTOT,NELT )
!
!
!     NOZERO: store in matrix only NOn-ZERO elements.
!
!------------------------------------------------------------------------------C
!     Used Modules.                                                
!------------------------------------------------------------------------------C
      USE FILINX  
      USE LINEQN  
!------------------------------------------------------------------------------C
!     Implicit Double Precision.
!------------------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
!------------------------------------------------------------------------------C
!     Type Declarations.
!------------------------------------------------------------------------------C  
      CHARACTER SS*10 
!------------------------------------------------------------------------------C
!     Store only non-zero elements.
!------------------------------------------------------------------------------C
      NELT = 1  
      DO 100 N=1,NTOT  
        IF( ALC(N) .EQ. 0.D+0 ) GOTO 100 
        IF( N .NE. NELT ) THEN  
          ALC(NELT) = ALC(N) 
          IA(NELT)  = IA(N)   
          JA(NELT)  = JA(N)  
          IF( IA(N) .LT. 1 .OR. IA(N) .GT. NVAR ) THEN 
            WRITE(*,*)' ERROR: Incorrect index of the matrix:',     &
                 ' IA: ',IA(N),N
            WRITE(IWR,*) ' ERROR: Incorrect index of the matrix:',  &
                 ' IA: ',IA(N),N
            STOP    
          ENDIF    
          IF( JA(N) .LT. 1 .OR. JA(N) .GT. NVAR ) THEN  
            WRITE(*,*) ' ERROR: Incorrect index of the matrix',     &
                 ' JA: ',JA(N),N
            WRITE(IWR,*) ' ERROR: Incorrect index of the matrix',   &
                 ' JA: ',JA(N),N
            STOP    
          ENDIF    
        ENDIF      
        NELT = NELT +1     
  100 CONTINUE 
      NELT = NELT -1       
!------------------------------------------------------------------------------C
!     End of NOZERO group.
!------------------------------------------------------------------------------C
!
      RETURN
      END

!==============================================================================C 
!==============================================================================C 
! 
      FUNCTION APRXM ( FDL,FDH,FL,FH,IDMN )
! 
!------------------------------------------------------------------------------C 
!  
!      APRXM: APpRoXiMation stencil for equation coefficients. 
!  
!------------------------------------------------------------------------------C  
      USE NUMBRS
!------------------------------------------------------------------------------C 
!      Implicit Double Precision. 
!------------------------------------------------------------------------------C 
      IMPLICIT REAL*8 (A-H,O-Z) 
      IMPLICIT INTEGER*4 (I-N) 
!------------------------------------------------------------------------------C 
!      Harmonic mean: default mode. 
!------------------------------------------------------------------------------C 
      IF( IDMN .EQ. 1 ) THEN 
        APRXM = (2.*FDL*FDH)/(FDL+FDH+SMALL) 
!------------------------------------------------------------------------------C 
!      Geometric mean. 
!------------------------------------------------------------------------------C 
      ELSEIF( IDMN .EQ. 2 ) THEN 
        APRXM = SQRT( FDL*FDH ) 
!------------------------------------------------------------------------------C 
!      Arithmetic mean. 
!------------------------------------------------------------------------------C 
      ELSEIF( IDMN .EQ. 3 ) THEN 
        APRXM = 0.5*( FDH+FDL ) 
!------------------------------------------------------------------------------C 
!      Upwind mean. 
!------------------------------------------------------------------------------C 
      ELSEIF( IDMN .EQ. 4 ) THEN 
        IF( FL .GE. FH ) THEN 
          APRXM = FDL 
        ELSE 
          APRXM = FDH 
        ENDIF 
      ENDIF 
!------------------------------------------------------------------------------C 
!      End of APRXM group. 
!------------------------------------------------------------------------------C 
!  
      RETURN 
      END 
 
!==============================================================================C 
!  
      FUNCTION DAPRXM ( DDL,FDL,DDH,FDH,FL,FH,IDMN ) 
!  
!------------------------------------------------------------------------------C 
!  
!      DAPRXM: Derivative of APpRoXiMation stencil for equation  
!                            coefficients. 
!  
!------------------------------------------------------------------------------C
!     Used Modules.                                                
!------------------------------------------------------------------------------C
      USE NUMBRS  
!------------------------------------------------------------------------------C 
!      Implicit Double Precision. 
!------------------------------------------------------------------------------C 
      IMPLICIT REAL*8 (A-H,O-Z) 
      IMPLICIT INTEGER*4 (I-N) 
!------------------------------------------------------------------------------C 
!      Harmonic mean: default mode. 
!------------------------------------------------------------------------------C 
      IF( IDMN .EQ. 1 ) THEN 
        DAPRXM = 2.*(DDL*FDH*FDH+DDH*FDL*FDL)/((FDL+FDH)**2+SMALL) 
!------------------------------------------------------------------------------C 
!      Geometric mean. 
!------------------------------------------------------------------------------C 
      ELSEIF( IDMN .EQ. 2 ) THEN 
        DAPRXM = 0.5*(DDL*FDH+DDH*FDL)/(SQRT( FDL*FDH )+SMALL) 
!------------------------------------------------------------------------------C 
!      Arithmetic mean. 
!------------------------------------------------------------------------------C 
      ELSEIF( IDMN .EQ. 3 ) THEN 
        DAPRXM = 0.5*( DDH+DDL ) 
!------------------------------------------------------------------------------C 
!      Upwind mean. 
!------------------------------------------------------------------------------C 
      ELSEIF( IDMN .EQ. 4 ) THEN 
        IF( FL .GE. FH ) THEN 
          DAPRXM = DDL 
        ELSE 
          DAPRXM = DDH 
        ENDIF 
      ENDIF 
!------------------------------------------------------------------------------C 
!      End of DAPRXM group. 
!------------------------------------------------------------------------------C 
!  
      RETURN 
      END 
 
!==============================================================================C 
!==============================================================================C 
!==============================================================================C 


