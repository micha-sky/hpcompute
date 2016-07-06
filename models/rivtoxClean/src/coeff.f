C======================================================================C 
C 
	SUBROUTINE COEFF 
C
C----------------------------------------------------------------------C
C
C     COEFF: calculate COEFFicients A(i),B(i),...,G1(i).
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
!      REAL*8, POINTER :: CF(:,:), CE(:,:)  !   CF(10,LMAXP),CE(6,LMAXP)    
      real*8 CF(10,MaxLPoints),CE(6,MaxLPoints)
C----------------------------------------------------------------------C
C     Allocate memory.
C----------------------------------------------------------------------C 
!      allocate(CF(10,MaxLPoints),CE(6,MaxLPoints)) 
      CF = 0.d0; CE = 0.d0;
C----------------------------------------------------------------------C
C     Link Forward Sweep.
C----------------------------------------------------------------------C 
      Fprev = 0.d0;
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
C     CF(5,i) - initial finite-difference scheme for mass conservation
C               equation with respect to Y and Q 
C     CF(1,i) - derivative of CF(5,i) with respect to Y(i+1) 
C     CF(2,i) - derivative of CF(5,i) with respect to Q(i+1) 
C     CF(3,i) - derivative of CF(5,i) with respect to Y(i) 
C     CF(4,i) - derivative of CF(5,i) with respect to Q(i) 
C----------------------------------------------------------------------C
C     Storage terms.
C----------------------------------------------------------------------C
          CF(1,I) = PHI*DAW(NE)/DT                         ! A(i) 
          CF(3,I) = (1.D+0 -PHI)*DAW(NP)/DT                ! C(i)
          CF(5,I) = PHI/DT*(AW(NE) -AWO(NE))               ! G(i)
     &        +(1.D+0 -PHI)/DT*(AW(NP) -AWO(NP))             
C----------------------------------------------------------------------C
C     Upstream (West) advection terms.
C----------------------------------------------------------------------C
          CF(2,I) =  THETA/DX
          CF(4,I) = -THETA/DX
          CF(5,I) = CF(5,I) + ( THETA*(QH(NE) -QH(NP)) 
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
          CF(5,i) = CF(5,i) -SRC/DX  
          RH(np) = CF(5,i)
C            
C----------------------------------------------------------------------C
C     Saint Venant Momentum Conservation Equation.
C----------------------------------------------------------------------C
C     CF(10,i) - initial finite-difference scheme for momentum 
C                conservation equation with respect to Y and Q 
C     CF(6,i)  - derivative of CF(10,i) with respect to Y(i+1) 
C     CF(7,i)  - derivative of CF(10,i) with respect to Q(i+1) 
C     CF(8,i)  - derivative of CF(10,i) with respect to Y(i) 
C     CF(9,i)  - derivative of CF(10,i) with respect to Q(i) 
          YP = TSF(NP) +HR(NP) 
          YE = TSF(NE) +HR(NE) 
          YPO = TSF(NP) +HRO(NP) 
          YEO = TSF(NE) +HRO(NE)    
C----------------------------------------------------------------------C
C     Storage terms.
C----------------------------------------------------------------------C
          CF(7,I)  = PHI/DT           
          CF(9,I)  = (1.D+0 -PHI)/DT          
          CF(10,I) = PHI/DT*(QH(NE) -QHO(NE))
     &             + (1.D+0 -PHI)/DT*(QH(NP) -QHO(NP))          
C----------------------------------------------------------------------C
C     Advection terms.
C----------------------------------------------------------------------C
          CF(6,I) = -ALPHA*THETA*DAW(NE)/DX *((THETA*(QH(NE) -QH(NP))              
     &    + (1.D+0 -THETA)*(QHO(NE) -QHO(NP))) /(AW(NE)**2)*QH(NE)
     &    + (THETA*(QH(NP)/AW(NP) +QH(NE)/AW(NP))**2            !! A1(i)
     &    + (1.D+0 -THETA)*(QHO(NP)/AWO(NP) +QHO(NE)/AWO(NE))**2)/4.D+0)                     

          CF(6,I) = CF(6,I) + ALPHA*THETA*DAW(NE)/2.D+0
     &              /(AW(NE)**2)*QH(NE)*(QH(NP)/AW(NP)
     &       + QH(NE)/AW(NE))/DX*(THETA*(AW(NE)-AW(NP))
     &                         + (1.D+0 -THETA)*(AWO(NE)-AWO(NP))) 
      
          CF(6,I) = CF(6,I) +0.5D+0*THETA*GRAV/DX* 
     &        ((THETA*(AW(NP)+AW(NE)) +(1.D+0 -THETA)*(AWO(NP)+AWO(NE)))
     &      + DAW(NE)*(THETA*(YE -YP) +(1.D+0 -THETA)*(YEO -YPO)))  
     &      - THETA*GRAV*(1.D+0 -BETA)*DFRE/(FRE**3)*QH(NE)
     &       *ABS(QH(NE))*(THETA*(AW(NP)+AW(NE))
     &         +(1.D+0 -THETA)*(AWO(NP)+AWO(NE)))

          CF(6,I) = CF(6,I) +THETA*GRAV*DAW(NE)/2.D+0      !! A1(i)
     &       *(THETA*(BETA*QH(NP)*ABS(QH(NP))/(FRP**2)
     &       + (1.D+0 -BETA)*QH(NE)*ABS(QH(NE))/(FRE**2))
     &       + (1D+0 -THETA)*(BETA*QHO(NP)*ABS(QHO(NP))/(FRPO**2)
     &       + (1D+0 -BETA)*QHO(NE)*ABS(QHO(NE))/(FREO**2))) 
C
          CF(7,I) = CF(7,I) +ALPHA*THETA/DX*                 !! B1(i)
     &            ((THETA*(QH(NP)/AW(NP) +QH(NE)/AW(NE)) 
     &          + (1.D+0 -THETA)*(QHO(NP)/AWO(NP) +QHO(NE)/AWO(NE)))
     &      + 1.D+0/AW(NE)*
     &        ((THETA*(QH(NE)-QH(NP)) +(1.D+0 -THETA)*(QHO(NE)-QHO(NP)))
     &      - 0.5D+0*(QH(NP)/AW(NP)+QH(NE)/AW(NE))*
     &       (THETA*(AW(NE)-AW(NP)) +(1.D+0 -THETA)*(AWO(NE)-AWO(NP)))))
     &      + GRAV*THETA*(1.D+0 -BETA)*ABS(QH(NE))/(FRE**2)*
     &        (THETA*(AW(NP)+AW(NE)) +(1.D+0 -THETA)*(AWO(NP)+AWO(NE))) 
C
          CF(8,I) = -ALPHA*THETA/DX*(DAW(NP)*(QH(NP)/(AW(NP)**2)*
     &         (THETA*(QH(NE)-QH(NP))+(1.D+0 -THETA)*(QHO(NE)-QHO(NP)))
     &      - 0.25D+0*(THETA*(QH(NP)/AW(NP) +QH(NE)/AW(NE))**2 
     &      +  (1.D+0 -THETA)*(QHO(NP)/AWO(NP) +QHO(NE)/AWO(NE))**2)) 
     &      - 0.5D+0*(DAW(NP)*QH(NP)**2/AW(NP)**3 
     &      + DAW(NP)*QH(NP)*QH(NE)/(AW(NE)*AW(NP)**2))*  
     &        (THETA*(AW(NE)-AW(NP)) +(1.D+0 -THETA)*(AWO(NE)-AWO(NP))))  

          CF(8,I) = CF(8,I) - THETA*GRAV/DX/2.D+0*          !! C1(i)
     &       ((THETA*(AW(NE)+AW(NP)) +(1.D+0 -THETA)*(AWO(NE)+AWO(NP)))  
     &      - DAW(NP)* 
     &        (THETA*(YE -YP) +(1.D+0 -THETA)*(YEO -YPO)))
     &      - GRAV*THETA*BETA*DFRP*QH(NP)*ABS(QH(NP))/(FRP**3)*
     &        (THETA*(AW(NP)+AW(NE)) +(1.D+0 -THETA)*(AWO(NP)+AWO(NE)))    

          CF(8,I) = CF(8,I) + THETA*GRAV*DAW(NP)/2.D+0*
     &            (THETA*(BETA*QH(NP)*ABS(QH(NP))/(FRP**2)  
     &          + (1.D+0 -BETA)*QH(NE)*ABS(QH(NE))/(FRE**2))
     &          + (1.D+0 -THETA)*(BETA*QHO(NP)*ABS(QHO(NP))/(FRPO**2)
     &          + (1.D+0 -BETA)*QHO(NE)*ABS(QHO(NE))/(FREO**2))) 
C
          CF(9,I) = CF(9,I) -ALPHA*THETA/DX*                 !! D1(i)
     &            ((THETA*(QH(NP)/AW(NP) +QH(NE)/AW(NE))
     &          + (1.D+0 -THETA)*(QHO(NP)/AWO(NP) +QHO(NE)/AWO(NE))
     &          - 1.D+0/AW(NP)*
     &      ((THETA*(QH(NE)-QH(NP)) +(1.D+0 -THETA)*(QHO(NE)-QHO(NP)))
     &          - 0.5D+0*(QH(NP)/AW(NP)+QH(NE)/AW(NE))*
     &      (THETA*(AW(NE)-AW(NP)) +(1.D+0 -THETA)*(AWO(NE)-AWO(NP))))))
     &          + GRAV*THETA*BETA*ABS(QH(NP))/(FRP**2)*
     &      (THETA*(AW(NP)+AW(NE)) +(1.D+0 -THETA)*(AWO(NP)+AWO(NE))) 
C
          CF(10,I) = CF(10,I) + ALPHA/DX*                   !! G1(i)  
     &               ((THETA*(QH(NP)/AW(NP) +QH(NE)/AW(NE))
     &        + (1.D+0 -THETA)*(QHO(NP)/AWO(NP)+ QHO(NE)/AWO(NE)))*
     &          (THETA*(QH(NE)-QH(NP))+(1.D+0 -THETA)*(QHO(NE)-QHO(NP)))
     &        - 0.25D+0*(THETA*(QH(NP)/AW(NP)+QH(NE)/AW(NE))**2 
     &        + (1.D+0 -THETA)*(QHO(NP)/AWO(NP)+QHO(NE)/AWO(NE))**2)*
     &         (THETA*(AW(NE)-AW(NP))+(1.D+0 -THETA)*(AWO(NE)-AWO(NP))))
     &        + 0.5D+0*GRAV/DX*
     &         (THETA*(AW(NP)+AW(NE))+(1.D+0 -THETA)*(AWO(NP)+AWO(NE)))*
     &          (THETA*(YE -YP) +(1.D+0 -THETA)*(YEO -YPO))

          CF(10,I) = CF(10,I) + 0.5D+0*GRAV* 
     &         (THETA*(AW(NP)+AW(NE))+(1.D+0 -THETA)*(AWO(NP)+AWO(NE)))*
     &         (THETA*(BETA*QH(NP)*ABS(QH(NP))/(FRP**2)
     &        + (1.D+0 -BETA)*QH(NE)*ABS(QH(NE))/(FRE**2))
     &        + (1.D+0 -THETA)*(BETA*QHO(NP)*ABS(QHO(NP))/(FRPO**2)
     &        + (1.D+0 -BETA)*QHO(NE)*ABS(QHO(NE))/(FREO**2)))     
          RQ(np) = CF(10,i)

C----------------------------------------------------------------------C
C     Source term.
C----------------------------------------------------------------------C 
C
          DO 300 J=3,5  
            CF(J,I) = -CF(J,I) 
            CF(J+5,I) = -CF(J+5,I)  
  300     CONTINUE  
C----------------------------------------------------------------------C
C     Store calculated coefficients.
C----------------------------------------------------------------------C 
          DO 350 J=1,10  
		  CRA(NP,J) = CF(J,I) 
  350     CONTINUE		  
  400   CONTINUE  
C----------------------------------------------------------------------C
C     Calculate coefficients E(i),F(i),H(i), and etc.
C----------------------------------------------------------------------C  
        DO 450 K=1,IL(L)-1
          I = IL(L) -K  
          NP = ND(L,I) 
          IF( K .GT. 1 ) GOTO 420
          DEN     = -CF(2,i)*CF(9,i) +CF(7,i)*CF(4,i)
          CE(1,i) = (CF(8,i)*CF(2,i) -CF(3,i)*CF(7,i))/DEN  ! E(IL-1)
          CE(2,i) = (CF(10,i)*CF(2,i) -CF(5,i)*CF(7,i))/DEN ! F(IL-1)
          CE(3,i) = (CF(1,i)*CF(7,i) -CF(6,i)*CF(2,i))/DEN  ! H(IL-1) 

          DEN     = -CF(4,i)*CF(7,i) +CF(9,i)*CF(2,i) 
          CE(4,i) = (CF(3,i)*CF(9,i) -CF(8,i)*CF(4,i))/DEN  ! E1(IL-1)
          CE(5,i) = (CF(5,i)*CF(9,i) -CF(10,i)*CF(4,i))/DEN ! F1(IL-1)
          CE(6,i) = (CF(6,i)*CF(4,i) -CF(1,i)*CF(9,i))/DEN  ! H1(IL-1)
          GOTO 430   

  420     DEN   =  CF(1,i)*CF(7,i) -CF(6,i)*CF(2,i)
          RL    = (CF(3,i)*CF(7,i) -CF(8,i)*CF(2,i))/DEN
          RMM   = (CF(4,i)*CF(7,i) -CF(9,i)*CF(2,i))/DEN
          RN    = (CF(5,i)*CF(7,i) -CF(10,i)*CF(2,i))/DEN
C 
          DEN = (CF(1,i) +CF(2,i)*CE(1,i+1))*RMM -CF(4,i)
          CE(1,i) = (CF(3,i) -RL*(CF(1,i) +CF(2,i)*CE(1,i+1)))/DEN
          CE(2,i) = (CF(5,i) -(CF(1,i) +CF(2,i)*CE(1,i+1))*RN
     &            - CF(2,i)*CE(2,i+1))/DEN
          CE(3,i) = -CF(2,i)*CE(3,i+1)/DEN
          CE(4,i) = CE(4,i+1)*(RL +RMM*CE(1,i))
          CE(5,i) = CE(5,i+1) +CE(4,i+1)*(RMM*CE(2,i) +RN)
          CE(6,i) = CE(6,i+1) +CE(3,i)*CE(4,i+1)*RMM
  430     CONTINUE
C----------------------------------------------------------------------C
C     Store calculated coefficients.
C----------------------------------------------------------------------C 
          DO 440 J=1,6  
            CRH(NP,J) = CE(J,I) 
  440     CONTINUE	  
          Fprev = Fprev +RH(np)*RH(np) +RQ(np)*RQ(np)
  450   CONTINUE
  500 CONTINUE     
C----------------------------------------------------------------------C
C     Free memory.
C----------------------------------------------------------------------C 
!      deallocate(CF,CE) 
C----------------------------------------------------------------------C 
C     End of COEFF group. 
C----------------------------------------------------------------------C 
C 
      RETURN 
      END 
 
C======================================================================C 
          