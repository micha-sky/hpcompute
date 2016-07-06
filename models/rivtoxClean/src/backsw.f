C======================================================================C 
C 
	SUBROUTINE BACKSW  
C
C----------------------------------------------------------------------C
C
C     BACKSW: proceed BACKward SWeep.
C             compute dY(1,L,i),dQ(1,L,i) via DNY(node).
C           
C----------------------------------------------------------------------C 
      USE DRLINK
      USE DRWFLD
      USE FLDINX 
      USE LINEQN 
      USE SOLVAR
C----------------------------------------------------------------------C
C     Implicit Double Precision.
C----------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C----------------------------------------------------------------------C
C     Link Backward Sweep.
C----------------------------------------------------------------------C
      DO 200 L=1,NLINK
        NFIRST = MU(L)       ! upward node 
        NLAST  = MD(L)       ! downward node   
        N1 = ND(L,1) 
        NL = ND(L,IL(L))  
C     Define surface water level at the end points of link.
        DH(N1) = DY(NFIRST)       
        DH(NL) = DY(NLAST)   
C     Define water discharge at the end points of link.
        DQ(N1) = CRH(N1,1)*DH(N1) +CRH(N1,2) +CRH(N1,3)*DH(NL)
        DQ(NL) = CRH(N1,4)*DH(N1) +CRH(N1,5) +CRH(N1,6)*DH(NL)
C
        DO 100 I=1,IL(L)-2 
          NP = ND(L,I) 
          A  = CRA(NP,1)
          B  = CRA(NP,2)
          C  = CRA(NP,3)
          D  = CRA(NP,4)
          G  = CRA(NP,5)
          A1 = CRA(NP,6)
          B1 = CRA(NP,7)
          C1 = CRA(NP,8)
          D1 = CRA(NP,9)
          G1 = CRA(NP,10)
          DEN = A*B1 -A1*B
          RL = (C*B1 -C1*B)/DEN
          RMM = (D*B1 -D1*B)/DEN
          RN = (G*B1 -G1*B)/DEN 
          NE = ND(L,I+1) 
          DH(NE) = DH(NP)*RL +DQ(NP)*RMM +RN
          DQ(NE) = DH(NE)*CRH(NE,1) +CRH(NE,2) +DH(NL)*CRH(NE,3)
  100   ENDDO 
  200 ENDDO   
C----------------------------------------------------------------------C 
C     End of BACKSW group. 
C----------------------------------------------------------------------C 
C 
      RETURN 
      END 
 
C======================================================================C 
C 
      SUBROUTINE CONVIT( NITER,NDTR,IRCODE,IERR ) 
C 
C----------------------------------------------------------------------C 
C 
C     CONVIT: testing CONVergence of Newton's ITeration for an 
C             nonlinear equation solution. 
C 
C----------------------------------------------------------------------C 
      USE CHARAC  
      USE FILINX 
      USE FLDINX
      USE DRLINK
      USE DRWFLD
      USE FLDINX  
      USE LANDSF
      USE NUMBRS
      USE POINTS 
      USE RESIDL
      USE SOLVAR
C----------------------------------------------------------------------C 
C     Implicit Double Precision. 
C----------------------------------------------------------------------C 
      IMPLICIT REAL*8 (A-H,O-Z) 
      IMPLICIT INTEGER*4 (I-N) 
C----------------------------------------------------------------------C 
C     Type Declarations. 
C----------------------------------------------------------------------C 
      CHARACTER*6 UNITS(6)  
      character eqtype(2)*35
!      character ffff*20
C----------------------------------------------------------------------C 
C     Data Statements. 
C----------------------------------------------------------------------C
      SAVE UNITS
      DATA UNITS/'sec','min','hour','day','week','year'/
!----------------------------------------------------------------------C  
      eqtype(1) = 'momentum eqn.'
      eqtype(2) = 'continuity eqn.'
C----------------------------------------------------------------------C 
C     Define time units prescribed by user. 
C----------------------------------------------------------------------C
      IUNIT = 1
      IF( TSCL .EQ. 6.D+1 ) THEN
        IUNIT = 2
      ELSEIF( TSCL .EQ. 3.6D+3 ) THEN
        IUNIT = 3 
      ELSEIF( TSCL .EQ. 8.64D+4 ) THEN
        IUNIT = 4 
      ELSEIF( TSCL .EQ. 6.048D+5 ) THEN
        IUNIT = 5 
      ELSEIF( TSCL .EQ. 3.15576D+7 ) THEN
        IUNIT = 6
      ENDIF 
C----------------------------------------------------------------------C 
C     Residual in terms of departure from previous step. 
C----------------------------------------------------------------------C 
      IF( IERR .EQ. 4 ) THEN
        RSD(1) = NPOINT
        NDRSD(1) = 1
        NITER = NRSD
        GOTO 600
      ENDIF
 
      IRCODE = 0 
      NRES = 0 
 
      RSU = 0.D+0 
      RSH = 0.D+0 
      RMXU = 0.D+0 
      RMXH = 0.D+0    
      DO 10 N=1,2 
        NDRSD(N) = 0 
   10 CONTINUE   
      rxQ = 0.d0 
      rxH = 0.d0
!=========================================================   
!      ffff = 'step  .txt'  
!      write(ffff(5:6),'(i2)') niter  
!      open(223,file=trim(ffff))   
!              write(223,'(2e15.6)') 0.d0,Fprev 
!              a1 = 0.d0; da = .01d0; a2 = 0d0;
!              do i=1,100  
!                a2 = a2 +da; 
!                call HQMODIFY (a1,a2)  
!                call CRSSEC   
!                call HQRESDL (F1)   
!                write(223,'(2e15.6)') a2,f1 
!                a1 = a2
!              enddo
!              close(223)
!=========================================================      
C 
      DO 150 L = 1,NLINK 
        DO 100 I=1,IL(L)  
          N = ND(L,I)  
          RESU = 0.D+0 
          RESH = 0.D+0 
C     Shallow Water Depth Equation.
          HR(N) = HR(N) +DH(N) 
          RESH = ABS(DH(N)/MAX(ABS(HR(N)),1.D-0)) 
          HR(N) = MAX(HR(N),ZERO)    
          IF( RESH .GT. RMXH ) THEN 
            NDRSD(2) = N 
            RMXH = RESH 
          ENDIF  
          IF( RESH .GT. 1.d+1 .OR. DH(N) .GT. 1.d+1 ) THEN   
            NITER = NRSD 
          ENDIF    
          
          trH = abs(rh(n))/max(abs(hr(n)),1.d0) 
          if( trH > rxH ) then 
            rxH = trH  
          endif  
C     Shallow X-Velocity Equation. 
          QH(N) = QH(N) + DQ(N)
          RESU = ABS(DQ(N)/MAX(abs(QH(N)),1.D+0)) 
          IF( RESU .GT. RMXU ) THEN 
            NDRSD(1) = N 
            RMXU = RESU  
          ENDIF
          IF( RESU .GT. 1.D+2 .or. abs(DQ(n)) > 1.d6 ) THEN 
C          IF( RESU.GT.1.D+20 .OR. ABS(UH(N)).GT.1.D+2 ) THEN 
            NITER = NRSD  
          ENDIF   

          trQ = abs(rq(n))/max(abs(qh(n)),1.d0)  
          if( trQ > rxQ ) then    
            rxQ = trQ  
          endif    
C 
          IF( IRSD .EQ. 1 ) THEN   
            RSU  = RSU + RESU*RESU 
            RSH  = RSH + RESH*RESH 
          ELSEif( irsd == 2 ) then  
            RSU  = MAX( RSU, RESU ) 
            RSH  = MAX( RSH, RESH )   
          else  
            RSU  = MAX( RSU, RxQ ) 
            RSH  = MAX( RSH, RxH )   
          ENDIF 
  100   CONTINUE 
  150 CONTINUE 

      IF( IRSD .EQ. 1 ) THEN 
C----------------------------------------------------------------------C 
C     Least square criterion. 
C----------------------------------------------------------------------C 
         RSD(1) = SQRT(RSU)/REAL(NPOINT) 
         RSD(2) = SQRT(RSH)/REAL(NPOINT) 
      ELSE 
C----------------------------------------------------------------------C 
C     Maximum residual criterion. 
C----------------------------------------------------------------------C 
         RSD(1) = RSU 
         RSD(2) = RSH 
      ENDIF
  600 CONTINUE
C----------------------------------------------------------------------C 
C     Test for convergence. 
C----------------------------------------------------------------------C  
      IF( RSD(1) .GT. RSDMX(1) .OR. RSD(2) .GT. RSDMX(2) ) THEN 
C----------------------------------------------------------------------C 
C     Test for convergence; do at most NRSD steps. 
C----------------------------------------------------------------------C 
        IF( NITER .LT. NRSD ) THEN 
          IRCODE = 1  
          
          call CRSSEC   
          call HQRESDL (Fnew)   
          if( Fnew > Fprev ) then   
            Alf_old = 1.d0; 
 700        Alf_new = 0.5d0*Alf_old;  
            if( Alf_new < 1.d-4 ) return
            
            call HQMODIFY (Alf_old,Alf_new)  
            call CRSSEC   
            call HQRESDL (Fcost)   
            x_l = Alf_new; x_c = Alf_old -Alf_new; x_r = Alf_old;
            det = x_c*Fprev -x_r*Fcost +x_l*Fnew 
            if( det > 0 ) then
              Alf_c = 0.5d0 *((x_r*x_r -x_l*x_l)*Fprev -x_r*x_r*Fcost
     &              +x_l*x_l*Fnew)/(x_c*Fprev -x_r*Fcost +x_l*Fnew)   
              
              call HQMODIFY (Alf_new,Alf_c)  
              call CRSSEC   
              call HQRESDL (Fcor)   
              if( Fcost < Fcor ) then 
                call HQMODIFY (Alf_c,Alf_new)  
                Alf_c = Alf_new   
              else  
                Fcost = Fcor  
              endif  
            endif    
            if( Fcost > Fprev ) then     
              Alf_old = Alf_c  
              Fnew = Fcost
              goto 700 
            elseif( Fcost < 1.d-6 ) then  
              IRCODE = 0  
            endif    
          endif    
C----------------------------------------------------------------------C 
C     Reduce the time step. 
C----------------------------------------------------------------------C 
        ELSEIF( NDTR.LT.MXDTR ) THEN
          SCLH = 5.D+0*SCLH 
          TIME  = TIME - DT 
          IF( IBREAK .EQ. 0 ) DTOLD = DT 
          IBREAK = IBREAK +1  
          DO 305 N=1,2 
            IF( RSD(N) .GT. RSDMX(N) ) THEN 
              NDRSD(N) = MAX(1,NDRSD(N))  
              WRITE (*,200) eqtype(n),RSD(N), 
     &                      ILK(NDRSD(N)),JPT(NDRSD(N)),HR(NDRSD(N)),
     &                      QH(NDRSD(N))
              WRITE (IWR,200) eqtype(n),RSD(N), 
     &                      ILK(NDRSD(N)),JPT(NDRSD(N)),HR(NDRSD(N)),
     &                      QH(NDRSD(N)) 
 
            ENDIF 
  305     CONTINUE 
          WRITE (*,300) DT/TSCL,DTR*DT/TSCL,UNITS(IUNIT)
          WRITE (IWR,300) DT/TSCL,DTR*DT/TSCL,UNITS(IUNIT) 
          DT    = DTR*DT 
          TIME  = TIME + DT 
          NITER = 0 
          DO N = 1,NPOINT 
            HR(N)  = HRO(N) 
            QH(N)  = QHO(N) 
            AW(N)  = AWO(N)
            PER(N) = PERO(N) 
C            TSF(N) = TSFO(N) 
          ENDDO 
          NDTR = NDTR +1 
          IRCODE = 2
C----------------------------------------------------------------------C 
C     Solution unconverged after the limit on time step reductions. 
C----------------------------------------------------------------------C 
        ELSE 
          DO 310 N=1,2 
            IF( RSD(N) .GT. RSDMX(N) ) THEN 
              WRITE (*,200) eqtype(n),RSD(N), 
     &                      ILK(NDRSD(N)),JPT(NDRSD(N)) 
              WRITE (IWR,200) eqtype(n),RSD(N), 
     &                      ILK(NDRSD(N)),JPT(NDRSD(N)) 
            ENDIF 
  310     CONTINUE 
          WRITE (*,400) 
          WRITE (IWR,400) 
          IRCODE = 3 
        ENDIF 
      ENDIF   
C----------------------------------------------------------------------C 
C     Format Statements. 
C----------------------------------------------------------------------C 
  200 FORMAT(' ** Convergence Failure: ',A16,' Max. Residual = ', 
     &       1PG11.4,' Node(I,J) ',I5,',',I5,3E16.6) 
  300 FORMAT(' ** Time Step Reduction from ',1PG11.4, 
     &       ' to ',1PG11.4,1X,A) 
  400 FORMAT(' *** Unconverged Solution ') 
C----------------------------------------------------------------------C 
C     End of CONVIT group. 
C----------------------------------------------------------------------C 
C 
      RETURN 
      END 
 
C======================================================================C 
C======================================================================C 
C 
	SUBROUTINE HQRESDL (Fcost)
C
C----------------------------------------------------------------------C
C
C     HQRESDL: calculate residuals of the continuity equation and   
!              momentum conservation equation.      
C           
C----------------------------------------------------------------------C 
      USE CONSTS 
      USE CRSECT 
      USE DRLINK 
      USE DRNAME 
      USE DRSOLV 
      USE DRWFLD 
      USE FLDINX  
      USE LANDSF
      USE LINEQN 
      USE POINTS 
      USE SOLVAR
C----------------------------------------------------------------------C
C     Implicit Double Precision.
C----------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C----------------------------------------------------------------------C 
C     Type Declarations. 
C----------------------------------------------------------------------C 
C----------------------------------------------------------------------C
C     Link Forward Sweep.
C----------------------------------------------------------------------C
      Fcost = 0.d0;
      DO 500 L=1,NLINK 
        DO 400 I=1,IL(L)-1  
          NP = ND(L,I) 
          NE = ND(L,I+1) 
          DX = X(NE) -X(NP)  
C----------------------------------------------------------------------C
C     Calculate conveyances for the points (L,I) and (L,I+1).
C----------------------------------------------------------------------C 
          FRP  = (AW(NP)/PER(NP))**(2.D+0/3.D+0)/ROU(NP)  
          FRE  = (AW(NE)/PER(NE))**(2.D+0/3.D+0)/ROU(NE)  
          DFRP  = 5.D+0/3.D+0*DAW(NP)*FRP -
     &    (2.D+0/3.D+0)*(AW(NP)/PER(NP))**(5.D+0/3.D+0)*DPER(NP)/ROU(NP)
          DFRE  = 5.D+0/3.D+0*DAW(NE)*FRE - 
     &    (2.D+0/3.D+0)*(AW(NE)/PER(NE))**(5.D+0/3.D+0)*DPER(NE)/ROU(NE)
          FRP  = AW(NP)*FRP  
          FRE  = AW(NE)*FRE  
          FRPO = AWO(NP)*(AWO(NP)/PERO(NP))**(2.D+0/3.D+0)/ROU(NP) 
          FREO = AWO(NE)*(AWO(NE)/PERO(NE))**(2.D+0/3.D+0)/ROU(NE)  
C----------------------------------------------------------------------C
C     Saint Venant Water-Mass Conservation Equation.
C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C     Storage terms.
C----------------------------------------------------------------------C
          RH(np) = PHI/DT*(AW(NE) -AWO(NE))               ! G(i)
     &        +(1.D+0 -PHI)/DT*(AW(NP) -AWO(NP))             
C----------------------------------------------------------------------C
C     Upstream (West) advection terms.
C----------------------------------------------------------------------C
          RH(np) = RH(np) + ( THETA*(QH(NE) -QH(NP)) 
     &        +(1.D+0 -THETA)*(QHO(NE) -QHO(NP)) )/DX
C----------------------------------------------------------------------C
C     Source term.
C----------------------------------------------------------------------C 
          SRC = 0.d0; 
          if( I == 1 ) then 
            SRC = SRC +FLQ(np)  
          else 
            SRC = SRC +0.5d0*FLQ(np)  
          endif  
          if( I+1 == IL(L) ) then  
            SRC = SRC +FLQ(np) 
          else 
            SRC = SRC +0.5d0*FLQ(np)  
          endif    
          RH(np) = RH(np) -SRC/DX  
C            
C----------------------------------------------------------------------C
C     Saint Venant Momentum Conservation Equation.
C----------------------------------------------------------------------C
          YP = TSF(NP) +HR(NP) 
          YE = TSF(NE) +HR(NE) 
          YPO = TSF(NP) +HRO(NP) 
          YEO = TSF(NE) +HRO(NE)    
C----------------------------------------------------------------------C
C     Storage terms.
C----------------------------------------------------------------------C
          RQ(np) = PHI/DT*(QH(NE) -QHO(NE))
     &             + (1.D+0 -PHI)/DT*(QH(NP) -QHO(NP))          
C----------------------------------------------------------------------C
C     Advection terms.
C----------------------------------------------------------------------C
          RQ(np) = RQ(np) + ALPHA/DX*                   !! G1(i)  
     &               ((THETA*(QH(NP)/AW(NP) +QH(NE)/AW(NE))
     &        + (1.D+0 -THETA)*(QHO(NP)/AWO(NP)+ QHO(NE)/AWO(NE)))*
     &          (THETA*(QH(NE)-QH(NP))+(1.D+0 -THETA)*(QHO(NE)-QHO(NP)))
     &        - 0.25D+0*(THETA*(QH(NP)/AW(NP)+QH(NE)/AW(NE))**2 
     &        + (1.D+0 -THETA)*(QHO(NP)/AWO(NP)+QHO(NE)/AWO(NE))**2)*
     &         (THETA*(AW(NE)-AW(NP))+(1.D+0 -THETA)*(AWO(NE)-AWO(NP))))
     &        + 0.5D+0*GRAV/DX*
     &         (THETA*(AW(NP)+AW(NE))+(1.D+0 -THETA)*(AWO(NP)+AWO(NE)))*
     &          (THETA*(YE -YP) +(1.D+0 -THETA)*(YEO -YPO))

          RQ(np) = RQ(np) + 0.5D+0*GRAV* 
     &         (THETA*(AW(NP)+AW(NE))+(1.D+0 -THETA)*(AWO(NP)+AWO(NE)))*
     &         (THETA*(BETA*QH(NP)*ABS(QH(NP))/(FRP**2)
     &        + (1.D+0 -BETA)*QH(NE)*ABS(QH(NE))/(FRE**2))
     &        + (1.D+0 -THETA)*(BETA*QHO(NP)*ABS(QHO(NP))/(FRPO**2)
     &        + (1.D+0 -BETA)*QHO(NE)*ABS(QHO(NE))/(FREO**2)))     

  400   CONTINUE  
C----------------------------------------------------------------------C
C     Calculate coefficients E(i),F(i),H(i), and etc.
C----------------------------------------------------------------------C  
        DO 450 K=1,IL(L)-1
          I = IL(L) -K  
          NP = ND(L,I) 
          Fcost = Fcost +RH(np)*RH(np) +RQ(np)*RQ(np)
  450   CONTINUE
  500 CONTINUE     
C----------------------------------------------------------------------C 
C     End of HQRESDL group. 
C----------------------------------------------------------------------C 
C 
      RETURN 
      END 
 
C======================================================================C 
C======================================================================C 
C 
	SUBROUTINE HQMODIFY (Alf_old,Alf_new)
C
C----------------------------------------------------------------------C
C
C     HQMODIFY: modify HR and QH according to new alpha for the   
!               damped Newton's method.      
C           
C----------------------------------------------------------------------C 
      USE CONSTS 
      USE CRSECT 
      USE DRLINK 
      USE DRNAME 
      USE DRSOLV 
      USE DRWFLD 
      USE FLDINX  
      USE LANDSF
      USE LINEQN 
      USE POINTS 
      USE SOLVAR
C----------------------------------------------------------------------C
C     Implicit Double Precision.
C----------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C----------------------------------------------------------------------C 
C     Type Declarations. 
C----------------------------------------------------------------------C 
C----------------------------------------------------------------------C
C     Link Forward Sweep.
C----------------------------------------------------------------------C
      DO 500 L=1,NLINK 
        DO 400 I=1,IL(L)  
          np = ND(L,I) 
C----------------------------------------------------------------------C
C     Update HR and QH.
C----------------------------------------------------------------------C 
          HR(np) = HR(np) -(Alf_old -Alf_new)*DH(np)  
          QH(np) = QH(np) -(Alf_old -Alf_new)*DQ(np)   
  400   CONTINUE  
  500 CONTINUE   
!      Alf_old = Alf_new
C----------------------------------------------------------------------C 
C     End of HQMODIFY group. 
C----------------------------------------------------------------------C 
C 
      RETURN 
      END 
 
C======================================================================C 
 
