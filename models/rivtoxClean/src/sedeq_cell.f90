!======================================================================C 
!----------------------------------------------------------------------C 
!
!     MODULE CONTAINING GLOBAL VARIABLES FOR SEDEQ_CELL FILE.
!
!----------------------------------------------------------------------C 
MODULE SEDEQ_CELL
 REAL*8 ZERO,SMALL,PI,GRAV,VISLQ,RHOLQ,RHOSE,DSE,D90,VISKN,SN,DSTAR,VSD
 REAL*8 DSTAR0_3,UCR,theta_cri0,HSMALL
END MODULE
MODULE SEDEQ_CELL_VAR
 REAL*8 Cb,SLP,TLP,HP,QP
END MODULE
!======================================================================C  
    
!======================================================================C 
SUBROUTINE SEDEQR
!----------------------------------------------------------------------C 
!
!     Prepare Gradients and call SEDEQRCell procedure.
!
!----------------------------------------------------------------------C  
USE CONSTS  
USE FLDINX
USE DRLINK
USE DRWFLD 
USE NUMBRS 
USE PLZONE  
USE POINTS
USE PROPER  
USE DRCATC 
USE DRSEDT 
USE DRSPFL 
USE LANDSF
USE SEDEQ_CELL, &
         sdqZERO=>ZERO, sdqSMALL=>SMALL, sdqPI=>PI, sdqGRAV=>GRAV, sdqVISLQ=>VISLQ, &
         sdqRHOLQ=>RHOLQ, sdqRHOSE=>RHOSE, sdqDSE=>DSE, sdqD90=>D90, sdqVISKN=>VISKN, &
         sdqSN=>SN, sdqDSTAR=>DSTAR, sdqVSD=>VSD, &
         sdqDSTAR0_3=>DSTAR0_3, sdqUCR=>UCR, sdqtheta_cri0=>theta_cri0
USE SEDEQ_CELL_VAR, &
         sdqCb=>Cb, sdqSLP=>SLP, sdqTLP=>TLP, sdqHP=>HP, sdqQP=>QP
IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)

INTEGER*4 N

!----------------------------------------------------------------------C 
!REAL*8 ls1,ls2,ls3
!----------------------------------------------------------------------C 
!----------------------------------------------------------------------!
!     SET CONSTANT sedeq_cell VARS.
!----------------------------------------------------------------------!
!----------------------------------------------------------------------!
!     Kinematic viscosity coefficient.
!----------------------------------------------------------------------!
 VISKN = VISLQ/RHOLQ
!----------------------------------------------------------------------!
!     Specific density.
!----------------------------------------------------------------------!
 SN = max(0.d0,GRAV*(RHOSE/RHOLQ - 1.d0))
 DSTAR = DSE*(SN/VISKN/VISKN)**(1.d0/3.d0)
 sdqDSE   = DSE
 sdqVISKN = VISKN
 sdqDSTAR = DSTAR
!----------------------------------------------------------------------!
!     Compute fall velocity according KINEROS model.
!----------------------------------------------------------------------!
 if (ILOAD/=9) then !--Non-cohesive
   !VSD = VfallCheng()   ! fall velocity !change tag
    VSD = VfallCheng(VSD)  
 else !--Cohesive
   VSD = DSE*DSE*GRAV/18.d0*(RHOSE-RHOLQ)/VISLQ
 endif

 sdqZERO  = 0.d0
 sdqSMALL = 1.d-20
 sdqPI    = acos(-1.d0)
 sdqGRAV  = GRAV
 sdqVISLQ = VISLQ
 sdqRHOLQ = RHOLQ
 sdqRHOSE = RHOSE
 sdqD90   = D90
 sdqSN    = SN
 sdqVSD   = VSD
!----------------------------------------------------------------------!
!     Precalculate constants.
!----------------------------------------------------------------------!
 sdqDSTAR0_3 = DSTAR**(-0.3d0)
 sdqUCR = SHIEDR(DSTAR) !--Critical shear velocity (Rijn)
 sdqtheta_cri0 = 0.3d0/(1.d0 +1.2d0*DSTAR) +0.055d0*(1.d0 -exp(-0.02d0*DSTAR)) !--Critical shear stress (Soulsby)

!----------------------------------------------------------------------!
!     No equilibrium concentrations and loads for cohesive.
!----------------------------------------------------------------------!
 if (ILOAD==9) RETURN
!----------------------------------------------------------------------!
!     SET SPACE VARIABLE sedeq_cell VARS.
!----------------------------------------------------------------------!
!----------------------------------------------------------------------!
!     Loop through all active cells.
!----------------------------------------------------------------------!
 do l=1,NLIN
   do i=1,IL(l)
     N = ND(l,i)  
     
     tH = HR(N)
     if( tH < SMALL ) CYCLE
!----------------------------------------------------------------------C 
!     Calc bottom, WSEL gradient.
!----------------------------------------------------------------------C 
   if( i > 1 ) then
     NW=ND(l,i-1)
     tLSW=0.5*(TSF(NW)+TSF(N))
     tWSW=0.5*(TSF(NW)+HR(NW)+TSF(N)+HR(N)) 
     XW = 0.5d0 *(X(N) +X(NW))
   else
     tLSW=TSF(N)
     tWSW=TSF(N)+HR(N) 
     XW = X(N)
   endif
   if( i < IL(l) ) then
     NE = ND(l,i+1)
     tLSE=0.5*(TSF(N)+TSF(NE))
     tWSE=0.5*(TSF(N)+HR(N)+TSF(NE)+HR(NE))  
     XE = 0.5d0 *(X(N) +X(NE))
   else
     tLSE=TSF(N)
     tWSE=TSF(N)+HR(N) 
     XE = X(N)
   endif
   tLSx = (tLSE-tLSW)/(XE -XW)
   tWSx = (tWSE-tWSW)/(XE -XW)

     IPLANT = IZPL(N)

     rFr=FRIC(IPLANT); 
     tCb=rFr; if (IFRIC(IPLANT)==3) tCb=tCb/(tH**r1d3)
     tCb=tCb*(1.d0+exp(-10.d0*tH))

     sdqCb  = tCb
     sdqSLP = tWSx
     sdqTLP = tLSx
     sdqHP  = HR(N)
     sdqQP  = QH(N)
     
     if(.not. Associated(BLX)) then 
         Allocate (BLX(SIZE(SEQ)))
         Allocate (SLX(SIZE(SEQ)))
     end if
     
     call SEDEQRCell(ILOAD,SEQ(N),BLX(N),SLX(N))
   enddo
 enddo
!----------------------------------------------------------------------C 
!     End of SEDEQR group. 
!----------------------------------------------------------------------C 
! 
RETURN
END
!======================================================================C 

!======================================================================C 
SUBROUTINE SEDEQRCell(ILOAD,SEQ,QBT,QST)
!----------------------------------------------------------------------C 
!
!     SEDEQRCell: compute SEDiment concentration at EQuilibRium transport
!                         capacity in CELL.
!
!     Transport Capacity Equations
!       ILOAD =  1 : Yalin bedload equation
!       ILOAD =  2 : modified Yalin bedload equation
!       ILOAD =  3 : Engelund-Hansen total load equation
!       ILOAD =  4 : Einstein-Brown bedload equation
!       ILOAD =  5 : Bagnold total load equation
!       ILOAD =  6 : Ackers-White total load equation
!       ILOAD =  7 : Rijn total load equation (1993)
!       ILOAD =  8 : simplified Rijn total-load equation
!
!       ILOAD = 10 : Bijker total load equation (W&C)
!       ILOAD = 11 : Van Rijn total load equation (W&C)
!       ILOAD = 12 : Dibajnia-Watanabe total load equation (W&C)
!       ILOAD = 13 : Meyer-Peter & Muller (1947) bedload equation (W&C)
!       ILOAD = 14 : Nielson bedload equation (W&C)
!       ILOAD = 15 : Soulsby bedload equation (W&C)
!       ILOAD = 16 : Camenen & Larson bedload equation (W&C)
!       ILOAD = 17 : Yang's transport equation (W&C)
!       ILOAD = 18 : Yang's transport equation for river
!       ILOAD = 19 : Van Rijn total load equation similar mike21  
!       ILOAD = 20 : Van Rijn total load equation (2004)
!
!----------------------------------------------------------------------C 
USE SEDEQ_CELL
USE SEDEQ_CELL_VAR 
IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
!----------------------------------------------------------------------C 
INTEGER*4 ILOAD
REAL*8 SEQ
!----------------------------------------------------------------------C 
REAL*8 QBT,QST
REAL*8 MSHIED
!----------------------------------------------------------------------C 
!INCLUDE 'MathStatLib' 
 sqr(r) = r*r 
 pwr3(r) = r*r*r
 pwr3_2(r) = r*sqrt(r)
!----------------------------------------------------------------------C 
 SEQ = 0.d0
 QBT = 0.d0 !--bed-load transport rate
 QST = 0.d0 !--suspended-load transport rate
 PVS = 0.d0 !--total transport rate
!----------------------------------------------------------------------!
!     CHECK.
!----------------------------------------------------------------------!
 if (DSE < SMALL) RETURN
 D90=max(D90,DSE)
 if (HP < 2.d0*D90) RETURN
!----------------------------------------------------------------------!
!     Compute current speed.
!----------------------------------------------------------------------! 
 UP = sqrt(2.d0)*HP*QP/sqrt(HP**4 +max(HP**4,small))
 UFL  = abs(UP)
!----------------------------------------------------------------------!
!     Check for zero current velocity.
!     Otherwise infinite equilibrium concentrations can arise
!     because bed load is simulated by suspended sediments.
!----------------------------------------------------------------------!
 if (UFL < SMALL) RETURN

!----------------------------------------------------------------------!
!     Compute shear velocity.
!----------------------------------------------------------------------!
 if (ILOAD>=1 .and. ILOAD<=6) then 
   USH  = sqrt(GRAV*HP*SLP) ! according energy slope
 
   xapa = 0.41d0          ! according current flow
   slm = 1.d3*DSE 
   tnu = 0.074d0*slm**1.19 
   skb = 27.7d0*tnu**2/slm   
   z0 = skb/30.d0 
!   USH = xapa*UFL/log(15.d0*HP/skb)

!   Cb1 = 0.02d0**2 *grav/hp**(1./3.)
   USH = sqrt(Cb)*UFL  ! Manning shear stress 
 endif    

!----------------------------------------------------------------------!
!     Yalin bedload equation.
!----------------------------------------------------------------------!
 if (ILOAD==1) then
   RNUM = DSE*USH*RHOLQ/VISLQ
   FNUM = USH*USH/(SN*DSE)
   YCR = SHIELD(RNUM)
   SPR = max(ZERO,FNUM/YCR-1.d0)
   APR = 2.45d0*sqrt(YCR)*(RHOLQ/RHOSE)**0.4
   QBT = 0.635d0*SPR*(1.d0 - log(1.d0 + APR*SPR)/max(APR*SPR,SMALL))
   QBT = DSE*USH*QBT
!----------------------------------------------------------------------!
!     Modified Yalin bedload equation.
!----------------------------------------------------------------------!
 elseif (ILOAD==2) then
   YPR = (DSE/VISLQ)*sqrt(RHOLQ*RHOLQ*SN*DSE)
   YCR = MSHIED(YPR)
   TAUCR = YCR*RHOLQ*SN*DSE
   FC = 0.006d0
   TAUB = 0.5d0*FC*RHOLQ*sqr(UFL)
   SPR = max(ZERO,TAUB/TAUCR-1.d0)
   APR = 2.45d0*sqrt(YCR)*(RHOLQ/RHOSE)**0.4
   QBT = 0.635d0*SPR*(1.d0 - log(1.d0 + APR*SPR)/max(APR*SPR,SMALL))
   QBT = DSE*USH*QBT
!----------------------------------------------------------------------!
!     Engelund-Hansen total load equation.
!----------------------------------------------------------------------!
 elseif (ILOAD==3) then
   tt1=USH/SN
   QST = 0.05d0*sqr(UFL)*sqr(tt1)*USH/DSE

!   YPR = (DSE/VISKN)*SQRT(SN*DSE)
!   YCR1 = MSHIED( YPR )  
!   YCR = SHIELD1( YPR )  
!   QST = 0.05D+0*UFL*UFL*SQRT(DSE/SN*YCR**3)   
   
   snp = rhose/rholq -1.d0
   QST = 0.05d0 *UFL*UFL *(SLP*HP)**1.5 /sqrt(grav)/DSE/snp**2
   QST = 1.488d0 *QST  ! convert in SI

!  according to Rijn   
   skb = 3.0d0*D90
   skb = 5.0d0*DSE
   z0 = skb/30.d0
   Cd = (0.4d0/(log(HP/z0) -1.d0))**2 
!   QST = 0.05d0 *Cd**1.5 *UFL**5/DSE/SN**2
   
 ! according to Cheng
  taub = Cb *UFL**2 
  teta = taub/SN/DSE 
!  QST = 0.05d0*UFL**2 *Teta**2.5 *sqrt(SN*DSE) *DSE/grav/SLP/HP
!  QST = QST*rholq/rhose
!----------------------------------------------------------------------!
!     Einstein-Brown bedload equation.
!----------------------------------------------------------------------!
 elseif (ILOAD==4) then
   tt1=USH*USH/(DSE*SN)
   QBT = 40.d0*VSD*DSE*pwr3(tt1)
!----------------------------------------------------------------------!
!     Bagnold total load equation.
!----------------------------------------------------------------------!
 elseif (ILOAD==5) then
   YPR = (DSE/VISLQ)*RHOLQ*sqrt(SN*DSE)
   YCR = MSHIED(YPR)
   TCR = YCR*RHOLQ*SN*DSE
   TCW = RHOLQ*USH*USH
   if (TCW>=TCR) then
     PW = min(600.D+0,0.7D+0*(TCW-TCR)/TCR)
     PK  = 0.005D+0*exp(PW)
   else
     PK = 0.d0
   endif
   QST = PK*TCW*UFL/(RHOLQ*SN)
!----------------------------------------------------------------------!
!     Ackers-White total load equation.
!----------------------------------------------------------------------!
 elseif (ILOAD==6) then
   DGR = DSE*(RHOLQ*RHOLQ*SN/sqr(VISLQ))**(1.d0/3.d0)
   if (DGR<1.d0) DGR = 1.d0
   if (DGR>60.d0) then
     AGR = 0.17d0
     PN  = 0.d0
 !    PM  = 1.5d0     ! edition 1973
     PM  = 1.78d0    ! edition 1990
     PC  = 0.025d0
   elseif (DGR>=1.d0) then
     AGR = 0.23d0/sqrt(DGR) + 0.14d0
     DGL = log10(DGR)
!     DGL = log(DGR)
!     PN  = 1.d0 - 0.56d0*DGL
     PN  = 1.d0 - 0.243d0*DGL  ! according Rijn
 !    PM  = 9.66d0/DGR + 1.34d0      ! edition 1973 
 !    PCL = 2.86d0*DGL - DGL*DGL - 3.53d0    ! edition 1973
     PM  = 6.83d0/DGR + 1.67d0      ! edition 1990 
!     PCL = 2.79d0*DGL - 0.98d0*DGL*DGL - 3.46d0  ! edition 1990
     PCL = 2.79d0*DGL - 0.426d0*DGL*DGL - 7.97d0  ! edition 1990 (Rijn)
     PC  = 10.d0**PCL
!     PC  = exp(PCL)
   endif 
   Fc = 0.24d0  !  pure-wave friction
   Fc = 0.006d0 ! for current only 
   !  USH = max(SMALL, AGR*sqrt(2.d0*DSE*SN/0.006d0))   ! bad approximation
!   USH2 = max(SMALL, 5.75d0*AGR*sqrt(DSE*SN)*log10(4.d0*HP/D90))  ! bad approximation

   xapa = 0.41d0 
   slm = 1.d3*DSE 
   tnu = 0.074d0*slm**1.19 
   skb = 27.7d0*tnu**2/slm   
!   skb1 = 2.5d0*DSE
   z0 = skb/30.d0 
   Fc = max(small,(0.41d0/log(0.5d0*HP/z0))**2) 
   USH1 = xapa*UFL/log(15.d0*HP/skb)

!  calculate skb according Rijn  
   sxi = UFL*UFL/SN/DSE 
   if( sxi <= 50.d0 ) then 
     sks = 15.d1 *DSE  
   elseif( sxi <= 250.d0 ) then 
     sks = DSE *(182.5d0 -0.65d0*sxi)   
   else 
     sks = 2.d1 *DSE  
   endif
   Fc = 18.d0 *log10(12.d0*HP/sks)
   USH6 = xapa*UFL/log(15.d0*HP/sks)  
   
   USH = sqrt(GRAV*HP*SLP) 

   Cb1 = 0.04d0**2 *grav/hp**(1./3.)
   USH8 = sqrt(Cb)*UFL  ! Manning shear stress 


   if (HP<=HSMALL) then
     FGR = 0.d0
   else
     tt1=log10(10.d0*HP/DSE)
!     tt1=log(10.d0*HP/DSE)
     if (tt1>0.d0) then
!       FGR = USH**PN/sqrt(DSE*SN)*(UFL/(5.75d0*tt1))**(1.d0 -PN)
       FGR = USH**PN/sqrt(DSE*SN)*(UFL/(2.46d0*tt1))**(1.d0 -PN)    ! (Rijn)
       if (FGR>AGR) then
         GGR = PC*(FGR/AGR -1.d0)**PM
       else
         GGR = 0.d0
       endif
     else
       GGR = 0.d0
     endif
   endif
   QST = UFL*DSE*GGR*(UFL/USH)**PN
!   QST = QST/RHOSE
!----------------------------------------------------------------------!
!     Rijn total load equation (1993).
!----------------------------------------------------------------------!
 elseif (ILOAD==7) then
!   DSIG = 0.5d0*(D16 +D84)/DSE
   DSIG = 2.5d0
   BETR = 1.d0
   CKAR = 0.4d0
   CR0  = 0.65d0

   CMN = 18.d0*log10(4.d0*HP/D90)
   USH = sqrt(GRAV)*UFL/CMN
   tt1=USH/UCR
   TRN = max(sqr(tt1)-1.d0,ZERO)
   QBT = 0.053d0*sqrt(SN*DSE)*DSE*TRN**2.1d0*DSTAR0_3

   ARF = min(max(ROUGH, 0.01d0*HP),0.5d0*HP) !--set in input--
   if (ARF<DSE) ARF = DSE
   CA = 0.015d0*DSE/ARF*pwr3_2(TRN)*DSTAR0_3
   DS = DSE*(1.d0 +0.011d0*(DSIG -1.d0)*(TRN -25.d0))
   VSD = SETVLR(DS)
   USTAR = max(SMALL, sqrt(GRAV*HP*TLP))
   UWR = VSD/USTAR
!   if (UWR>1.d-1 .and. UWR<1.d0) then
   RBET = 1.d0 + 2.d0*UWR*UWR
!   else
!     RBET = 1.d0
!   endif
!   if (UWR>=1.d-2 .and. UWR<=1.d0) then
   RPHI = 2.5d0*UWR**0.8*(CA/CR0)**0.4
!   else
!     RPHI = 1.d0
!   endif
   ZR  = UWR/CKAR/RBET
   ZRS = ZR + RPHI
   if ((UWR>0.d0) .and. (UWR<=1.d0)) then
     AD = ARF/HP
     if (abs(1.2d0 -ZRS)<SMALL) then
       FQS = 0.d0
     else
       FQS = (AD**ZRS - AD**1.2)/(1.d0 -AD)**ZRS/(1.2d0 -ZRS)
     endif
   else
     FQS = 0.d0
   endif
   QST = FQS*UFL*HP*CA
   350 CONTINUE
!----------------------------------------------------------------------!
!     Unified Simplified Rijn total-load equation.
!----------------------------------------------------------------------!
 elseif (ILOAD==8) then  
   RHY = log10(4.d0*HP/D90)
   if ( DSE <= 0.5d-3 ) then
     UCRc = 0.19d0 *(DSE)**0.1d0 *RHY 
   else
     UCRc = 8.5d0 *(DSE)**0.6d0 *RHY
   endif  
   
   if (HP <= DSE) then
     RLD = 0.d0
   else
     RLD = DSE/HP
   endif
   BRC = (max(0.d0,UFL-UCRc)/sqrt(SN*DSE))
   QBT = 0.015d0*UFL*HP *RLD**1.2d0 *BRC**1.5d0  
   QST = 0.012d0*UFL*DSE *BRC**2.4d0 *DSTAR**(-0.6d0) 
   360 CONTINUE
!----------------------------------------------------------------------!
!     Rijn total load equation (according mike 21).
!----------------------------------------------------------------------!
 elseif (ILOAD==19) then
   DSIG = 2.5d0
   BETR = 1.d0
   CKAR = 0.4d0
   CR0  = 0.65d0

   CMN = 18.d0*log10(4.d0*HP/D90)
   USH = sqrt(GRAV)*UFL/CMN
   tt1=USH/UCR
   TRN = max(sqr(tt1)-1.d0,ZERO)
   QBT = 0.053d0*sqrt(SN*DSE)*DSE*TRN**2.1d0*DSTAR0_3

   Arf = max(2.d0*DSE, 0.01d0*HP)  
   Uf = sqrt(Cb)*UFL 
   
   if( (DSTAR<=1.d1 .and. Uf>4.d0*VSD/DSTAR) .or.                     &
       (DSTAR>1.d1 .and. Uf>0.4d0*VSD) ) then  
     Ca = 0.015d0*DSE/Arf*pwr3_2(TRN)*DSTAR0_3
   else  
     QST = 0.d0 
     return
   endif  
       
   UWR = VSD/Uf
   RBET = 1.d0 + 2.d0*UWR*UWR
   RPHI = 2.5d0*UWR**0.8*(Ca/CR0)**0.4

   ZR  = UWR/CKAR/RBET
   ZRS = ZR + RPHI
!   if ((UWR>0.d0) .and. (UWR<=1.d0)) then
     AD = Arf/HP
     if (abs(1.2d0 -ZRS)<SMALL) then
       FQS = 0.d0
     else
       FQS = (AD**ZRS - AD**1.2)/(1.d0 -AD)**ZRS/(1.2d0 -ZRS)
     endif
!   else
!     FQS = 0.d0
!   endif
   QST = FQS*UFL*HP*CA
!----------------------------------------------------------------------!
!     Rijn total load equation (2004) (for fine sand and silt).
!----------------------------------------------------------------------!
 elseif (ILOAD==20) then  !DG Rijkswaterstaat. Description of TRANSPOR2004 and
                          !Implementation in Delft3D-ONLINE
     !Predictiono f SuspendedB ed MaterialT ransporti n Flows Over Silt
     !and Very Fine Sand (JAN H. VAN DEN BERG AND ANDR]� VAN GELDER)
   DSIG = 2.5d0
   BETR = 1.d0
   CKAR = 0.4d0
   CR0  = 0.65d0 
 
   ! my modification
   xapa = 0.41d0 
   slm = 1.d3*DSE 
   tnu = 0.074d0*slm**1.19 
   skb = 27.7d0*tnu**2/slm   
   z0 = skb/30.d0 
   Fc = max(small,(0.41d0/log(0.5d0*HP/z0))**2) 
   USH1 = xapa*UFL/log(15.d0*HP/skb)

!  calculate skb according Rijn  
   sxi = UFL*UFL/SN/DSE  
!  current-related roughness due to small-scale ripples   
   if( sxi <= 50.d0 ) then 
     sks = 15.d1 *DSE  
   elseif( sxi <= 250.d0 ) then 
     sks = DSE *(182.5d0 -0.65d0*sxi)   
   else 
     sks = 2.d1 *DSE  
   endif   
   sks = max(0.002d0, min(0.075,sks)) 
   
!  current-related roughness due to mega ripples   
   skm = 0.d0 
   if( HP > 1.d0 .and. UFL > 0.3d0 ) then
     if( sxi <= 50.d0 ) then 
       skm = 0.01d0 *HP  
     elseif( sxi <= 550.d0 ) then 
       skm = HP *(0.011d0 -0.00002d0*sxi)   
     else 
       skm = 0.d0     
     endif
   endif   
   
!  current-related roughness due to dunes in rivers   
   skd = 0.d0 
   if( HP > 1.d0 .and. UFL > 0.3d0 ) then
     if( sxi <= 100.d0 ) then 
       skd = 0.0004d0 *HP *sxi  
     elseif( sxi <= 600.d0 ) then 
       skd = HP *(0.048d0 -0.0008d0*sxi)   
     else 
       skd = 0.d0     
     endif
   endif   
   skd = max(0.02d0,skd); skd = min(1.d0,skd);    
   
   Arf = min(0.2d0*HP, max(0.5d0*sks,0.01d0));
   sks = sqrt(sks*sks +skm*skm +skd*skd);
   
   Fc = 18.d0 *log10(12.d0*HP/sks)
!   USH6 = xapa*UFL/log(15.d0*HP/sks)  
   
!   USH8 = sqrt(Cb)*UFL  ! Manning shear stress 
   ! my modification
   
!   sks = 12.d0*HP /10**min(20.d0,UFL/(18.d0*sqrt(HP*SLP))) 
!   sks = min(0.03d0,sks)
!   Arf = max(sks, 0.01d0*HP)  
!   Arf = min(Arf, 0.05d0*HP) 
!   Arf = min(Arf, 0.015d0)

   Csh = 18.d0*log10(12.d0*HP/sks)
   Cmn = 18.d0*log10(4.d0*HP/D90)
   USH = sqrt(GRAV)*UFL/Csh
!   USH2 = sqrt(GRAV)*UFL/Cmn

   teta = grav/Cmn**2 *UFL**2/SN/DSTAR  
!   teta = UFL**2/SN/DSTAR   ! my modification not corrected

!   if( DSTAR <= 6.d0 ) then
!     teta_cr = 0.109d0/sqrt(DSTAR)
   if( DSTAR <= 4.d0 ) then
     teta_cr = 0.115d0/sqrt(DSTAR)
   elseif( DSTAR <= 10.d0 ) then
     teta_cr = 0.14d0/DSTAR**0.64
   elseif( DSTAR <= 20.d0 ) then
     teta_cr = 0.04d0/DSTAR**0.1
   elseif( DSTAR <= 150.d0 ) then
     teta_cr = 0.013d0*DSTAR**(0.29)
   else
     teta_cr = 0.055
   endif 
   
   TRN = max(teta/teta_cr -1.d0,ZERO)
   QBT = 0.5d0 *DSE *USH *TRN *DSTAR0_3

   Ca = 0.015d0 *DSE/Arf *TRN**1.5 *DSTAR0_3  
   Ca = min(Ca, 0.05d0 *RHOSE)
   
   Az = 20.d0*Arf/HP
   Uf = VSD*(1.d0 -Az*Ca)**4 
   
   UWR = Uf/USH
   RBET = 1.d0 + 2.d0*UWR*UWR 
   RBET = min(2.d0,RBET)
   
   Za  = UWR/CKAR/RBET
   RPHI = 3.5d0*Za**0.8*(Az*Ca/CR0)**0.4

   Zrs = Za + RPHI
   AD = Arf/HP
   FQS = (AD**Zrs - AD**1.2)/(1.d0 -AD)**Zrs/(1.2d0 -Zrs)
   QST = FQS*UFL*HP*Ca
!----------------------------------------------------------------------!
!     Gravel (W&C).
!----------------------------------------------------------------------!
 elseif (ILOAD>=13 .and. ILOAD<=18) then
   if (HP<=0.001d0) GOTO 390
   xapa = 0.41d0 !--von Karman constant

   z0s = DSE/12.d0 !--total bed roughness height (FST2DH)

   !-----calculate bed shear stress under current-----!
   tt1=xapa/(1.d0 +log(z0s/HP))
   fc = 2.d0*tt1*tt1
   tau_c = 0.5d0*fc*sqr(UFL)

   !-----calculate bed shear stress under wave-----!
   tau_w = 0.d0
!   if (UBW <= SMALL) GOTO 381

!   fws = 1.d0; Rw = WAW/z0s;
!   if (Rw>=20 .and. Rw<200) then
!     fws = 18.d0/Rw
!   elseif (Rw>=200 .and. Rw<11000) then
!     fws = 1.39d0/Rw**0.52
!   elseif ( Rw>=11000) then
!     fws = 0.112d0/Rw**0.25
!   endif
!   tau_ws = 0.5d0*fws*UBW2
!   tau_ws = tau_ws/(SN*DSE)
!   z0t = 5.67d0*DSE*sqrt(max(0.d0,tau_ws -0.05d0))
!
!   xsi = UBW2/(SN*DSE) !--this is for sand ripples (FST2DH)
!   derw = 0.d0
!   dlyw = 0.d0
!   if (xsi<=10.d0) then
!!     derw = WAW*(0.275d0 - 0.022d0*sqrt(xsi))
!!     dlyw = derw/0.14d0 !--Grasmeijer & Kleinhans
!
!     derw = 0.22d0*WAW     ! van Rijn (1989)
!     dlyw = 1.25d0*WAW
!   elseif (xsi>10.d0 .and. xsi<=25.d1) then
!!     derw = WAW*2.d0/xsi     !--Grasmeijer & Kleinhans
!!     dlyw = derw/(-0.078d0 +0.355d0/xsi**0.221)
!
!     derw = 2.8d-13*(25.d1 -xsi)**5  *WAW     ! van Rijn (1989)
!     dlyw = 1.4d-6*(25.1 -xsi)**2.5 *WAW
!   endif
!!   if (tau_ws<0.831d0) then !--(FST2DH)
!!     dlyw = derw/(0.182d0 -0.022d0*pwr3_2(tau_ws))
!!   endif
!   if (dlyw<=SMALL) then
!     z0fw = 0.d0
!   else
!!    z0fw = 0.267d0*derw*derw/dlyw       !--Grasmeijer & Kleinhans
!     z0fw = 2.d0/3.d0*derw*derw/dlyw            ! van Rijn (1989)
!   endif
!
!   dlyc = 10d3*DSE
!   derc = dlyc/7.d0
!   z0fc = 0.267d0*derc*derc/dlyc
!
!!   z0 = max(z0s+z0fw+z0t,5.d-5)          !temp - account for sand ripples - Camenen & Larson gives strange results (d50 increases - SEQ also increases, and very strangely, see seq_d50_3m.emf)
!   z0 = max(z0s,5.d-5)                         !temp - don't account for sand ripples - Camenen & Larson SEQ slowly decrease while d50 increase
!   zz = z0*3.d1       !   it is ks

!   Rw = WAW/zz
!   if (Rw <= 0.5d0) then
!     fw = 1.43967d0*sqrt(Rw)
!   else
!     fw = exp(5.5d0/Rw**0.2 -6.3d0)       ! Nielson
!   endif
!   tau_w = 0.5d0*fw*UBW2

   381 CONTINUE
   !-----calculate bed shear stress under combined wave and current-----!
   tau_m = tau_c*(1.d0 +1.2d0*(tau_w/max(tau_w+tau_c,SMALL))**3.2)

   !-----Shield's parameter-----!
   theta = tau_m/(SN*DSE)

   !-----skin friction-----!
!   theta_skin = 0.06d0 +0.4d0*theta*theta

   !-----calculate bed slope-----!
   d_tlp=0.d0; if (UFL>1.d-10) d_tlp=(TLPx*UP +TLPy*VP)/UFL;

   !-----calculate critical Shields parameter-----!
!   YPR = (DSE/VISLQ)*sqrt(RHOLQ*RHOLQ*SN*DSE)
!   theta_cri0 = MSHIED(YPR)          ! Yalin's method

   !-----take into account bed slope-----!
   theta_cri = theta_cri0*(1.d0 +d_tlp)
   coef1 = DSE*sqrt(SN*DSE)

!   Ufric = UFL/(6.d0 +2.5d0*log(HP/z0s) )     ! friction velocity
!   theta = Ufric*Ufric/(SN*DSE)                    ! skin-friction
   !----------------------------------------------------------------------!
   !     Meyer_Peter & Muler bedload equation (W&C).
   !----------------------------------------------------------------------!
   theta_skin = theta
   if (ILOAD==13) then
     tt1=max(0.d0,theta_skin -theta_cri)
     QBT = 8.d0*coef1*pwr3_2(tt1)
   !----------------------------------------------------------------------!
   !     Nielson bedload equation (W&C).
   !----------------------------------------------------------------------!
   elseif (ILOAD==14) then
     QBT = 12.d0*coef1*max(0.d0,theta -theta_cri)*sqrt(theta)
   !----------------------------------------------------------------------!
   !     Soulsby bedload equation (W&C).
   !----------------------------------------------------------------------!
   elseif (ILOAD==15) then
     tt1=max(0.d0,theta_skin -theta_cri)
     QBT = 5.1d0*coef1*pwr3_2(tt1)
   !----------------------------------------------------------------------!
   !     Camenen & Larson bedload equation (W&C).
   !----------------------------------------------------------------------!
   elseif (ILOAD==16) then
     theta_c = tau_c/(SN*DSE)
     theta_w = tau_w/(SN*DSE)
     theta_wm = 0.5d0*theta_w
     theta_cw   = sqrt(sqr(theta_c) +sqr(theta_w) +2.d0*theta_w*theta_c*cosWCANG)
     theta_cwm = sqrt(sqr(theta_c) +sqr(theta_wm) +2.d0*theta_wm*theta_c*cosWCANG)
     QBT = 12.d0*coef1*sqrt(theta_c)*theta_cwm*exp(-4.5d0*theta_cri/theta_cw)
   !----------------------------------------------------------------------!
   !     Yang's transport equation (W&C).
   !----------------------------------------------------------------------!
   elseif (ILOAD==17) then
     U_skin = sqrt(theta_skin)
     U_energy = sqrt(theta)

     par_1 = U_skin*DSE/VISKN
     par_2 = VSD*DSE/VISKN
     par_3 = U_skin/VSD

     U_cri = 3.6d3*VSD
     if (par_1>1.2d0 .and. par_1<7.d1) then
       U_cri = (0.66d0 +2.5d0/(log10(par_1) -0.6d-1))*VSD
     elseif (par_1>=70.d0) then
       U_cri = 2.05d0*VSD
     endif

     if (DSE<2.d-3) then !--sand equation
       pow_m = 5.435d0 -0.286d0*log10(par_2) -0.457d0*log10(par_3)
       pow_n = 1.799d0 -0.409d0*log10(par_2) -0.314d0*log10(par_3)
     else !--gravel equation
       pow_m = 6.681d0 -0.633d0*log10(par_2) -4.816d0*log10(par_3)
       pow_n = 2.784d0 -0.305d0*log10(par_2) -0.282d0*log10(par_3)
     endif

!     QBT = HP*UFL*10.d0**(-6.d0 +pow_m +pow_n*log10(U_energy*(UFL-U_cri)/VSD))
     if (UFL>U_cri) QBT = HP*UFL*10.d0**(-9.d0 +pow_m +pow_n*log10(SLP*(UFL-U_cri)/VSD))
     QBT = RHOLQ/(RHOSE-RHOLQ)*QBT
!----------------------------------------------------------------------C  
!        Yang's transport equation for river (1979).
!----------------------------------------------------------------------C  
   elseif (ILOAD==18) then
     U_skin = sqrt(Grav*HP*SLP)
     U_energy = UFL*SLP

     par_1 = U_skin*DSE/VISKN
     par_2 = VSD*DSE/VISKN
     par_3 = U_skin/VSD

!     U_cri = 3.6d3*VSD
     U_cri = 1.31d2*VSD  ! my approximation
     if (par_1>=1.2d0 .and. par_1<7.d1) then
       U_cri = (0.66d0 +2.5d0/(log10(par_1) -0.6d-1))*VSD
     elseif (par_1>=70.d0) then
       U_cri = 2.05d0*VSD
     endif

     QBT = 0.d0
     if (DSE < 2.d-3) then !--sand equation
       pow_m = 5.435d0 -0.286d0*log10(par_2) -0.457d0*log10(par_3)
       pow_n = 1.799d0 -0.409d0*log10(par_2) -0.314d0*log10(par_3)
       if (UFL > U_cri) QBT = 10.d0**(pow_m +pow_n*log10(SLP*(UFL-U_cri)/VSD))
       if( QBT >= 1.d2 ) then 
         pow_m = 5.165d0 -0.153d0*log10(par_2) -0.297d0*log10(par_3)
         pow_n = 1.780d0 -0.360d0*log10(par_2) -0.480d0*log10(par_3) 
         QBT = 10.d0**(pow_m +pow_n*log10(SLP*UFL/VSD)) 
       endif    
     else !--gravel equation
       pow_m = 6.681d0 -0.633d0*log10(par_2) -4.816d0*log10(par_3)
       pow_n = 2.784d0 -0.305d0*log10(par_2) -0.282d0*log10(par_3)
       if (UFL > U_cri) QBT = 10.d0**(pow_m +pow_n*log10(SLP*(UFL-U_cri)/VSD))
     endif

!     QBT = HP*UFL*10.d0**(-6.d0 +pow_m +pow_n*log10(U_energy*(UFL-U_cri)/VSD))
     QBT = HP*UFL*QBT*1.d-6
!     QBT = RHOLQ/(RHOSE-RHOLQ)*QBT  
     QBT = RHOLQ/RHOSE*QBT  
   endif

   390 CONTINUE
!----------------------------------------------------------------------!
 endif
!----------------------------------------------------------------------!
 PVS = QBT + QST

 if ((ILOAD==7) .or. (ILOAD==8)) VSD = SETVLR(DSE)

! PTR = (RHOSE-RHOLQ)*PVS
 PTR = RHOSE*PVS  ! PVS in m^3 of sediment / m^3 of water
 SEQ = PTR/(HP*UFL)

!----------------------------------------------------------------------C 
!     End of SEDEQRCell group. 
!----------------------------------------------------------------------C 
! 
RETURN
END
!======================================================================C 

!======================================================================C 
SUBROUTINE FAL(D,G,AN,PL,A)
!----------------------------------------------------------------------C 
IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
!----------------------------------------------------------------------C 
!INCLUDE 'MathStatLib'
pwr3_2(r) = r*sqrt(r)
!----------------------------------------------------------------------C 
 A1=sqrt(G*D)
 AR=sqrt(G)*pwr3_2(D)/AN
 A=A1*(sqrt(2.d0*PL/3.d0+36.d0/AR/AR)-6.d0/AR)
!----------------------------------------------------------------------C 
!     End of FAL group. 
!----------------------------------------------------------------------C 
! 
RETURN
END
!======================================================================C 


!======================================================================C 
SUBROUTINE FAL_1(D,A)
!----------------------------------------------------------------------C 
!  Wychislenie gidrawlicheskoj krupnosti a   (cm/sek)
!  po
!  diametru chasticy D  (cm),uskoreniyu sw.pad.G(cm/sek)
!  kinemat.wqzkosti aN (cm/sek/sek) i otnositelxnoj
!  plotnosti chastic PL
!----------------------------------------------------------------------C 
IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
!----------------------------------------------------------------------C 
 dsm=d*1000.d0
 if (dsm<=0.12d0) then
   a=64.1217d0*dsm*dsm+0.3948d0*dsm-0.0097d0
 else
   a=11.3916d0*dsm-0.3374d0
 endif
 a=0.01d0*a
!----------------------------------------------------------------------C 
!     End of FAL_1 group. 
!----------------------------------------------------------------------C 
! 
RETURN
END
!======================================================================C 

!======================================================================C 
SUBROUTINE SETVL
!----------------------------------------------------------------------C 
!
!     SETVL: compute SETtling VeLocity.
!
!----------------------------------------------------------------------C 
USE SEDEQ_CELL
IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
!----------------------------------------------------------------------C 
!INCLUDE 'MathStatLib'
 pwr3(r) = r*r*r
 sqr(r) = r*r
!----------------------------------------------------------------------C 
!----------------------------------------------------------------------!
!     Vanoni(1977), Raudkivi(1976) and Graf(1977) formulation for the
!     drag coefficient.
!----------------------------------------------------------------------! 
!goto 200
!----------------------------------------------------------------------!
!     Initial guess for the Reynolds number.
!----------------------------------------------------------------------!
 REYN = 1.d0
!----------------------------------------------------------------------!
!     Iteration loop for calculation of the Reynolds number.
!----------------------------------------------------------------------!
 100 CONTINUE
  tt1=REYN**0.687
  RS = 18.d0*sqr(VISKN)*REYN*(1.d0 + 0.15d0*tt1) - SN*pwr3(DSE)
  DERS = 18.d0*sqr(VISKN)*(1.d0 + 0.25305d0*tt1)
  DREY = -RS/max(DERS,SMALL)
  REYN = REYN + DREY
 if (abs(DREY/max(REYN,SMALL))>1.d-6) GOTO 100
!----------------------------------------------------------------------!
!     Compute settling velocity.
!----------------------------------------------------------------------!
 VSD = REYN*VISKN/DSE   
 RETURN
200 continue 
!----------------------------------------------------------------------!
!     Gibbs, Mathews, and Link(1971) formula for settling velocity.
!----------------------------------------------------------------------!
 VSD = (-3.d0*VISLQ + sqrt(9.d0*sqr(VISLQ) + &
       0.25d0*sqr(DSE*RHOLQ)*SN* &
      (0.01555d0 +0.0992*DSE)))/(RHOLQ*(0.0116d0 +0.0744d0*DSE))
!----------------------------------------------------------------------C 
!     End of SETVL group. 
!----------------------------------------------------------------------C 
! 
RETURN
END
!======================================================================C 


!======================================================================C 
REAL*8 FUNCTION SETVLR(DS)
!----------------------------------------------------------------------C 
!
!     SETVLR: SETtling VeLocity for the Rijn formula.
!
!----------------------------------------------------------------------C 
USE SEDEQ_CELL
IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
!----------------------------------------------------------------------C 
REAL*8 DS
!----------------------------------------------------------------------C 
!----------------------------------------------------------------------!
!     Specific density.
!----------------------------------------------------------------------!
 RS = SN*DS
!----------------------------------------------------------------------!
!     Stokes-range.
!----------------------------------------------------------------------!
 if (DS<1.d-4) then
   SETVLR = RS*DS/18.d0/VISKN
 elseif (DS<=1.d-3) then
!----------------------------------------------------------------------!
!     Zanke equation.
!----------------------------------------------------------------------!
   DSR = DS/VISKN
   SETVLR = 10.d0*VISKN/DS*(sqrt(1.d0 + 1.d-2*RS*DSR*DSR) -1.d0)
 else
   SETVLR = 1.1d0*sqrt(RS)
 endif
!----------------------------------------------------------------------C 
!     End of SETVLR group. 
!----------------------------------------------------------------------C 
! 
RETURN
END
!======================================================================C 

!======================================================================C 
REAL*8 FUNCTION SETVLRB(DS)
!----------------------------------------------------------------------C 
!
!     SETVLR: SETtling VeLocity for the Rubey's (1933) formula.
!
!----------------------------------------------------------------------C 
USE SEDEQ_CELL
IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
!----------------------------------------------------------------------C 
REAL*8 DS
!----------------------------------------------------------------------C   
  Temp = 20.d0; 
  VISKN1 = 1.0d0 +Temp*(0.0337d0 +0.000221*Temp) 
  VISKN1 = 1.792d-6/VISKN1;
!----------------------------------------------------------------------!
!     Specific coefficient.
!----------------------------------------------------------------------! 
 if( DS <= 1.d-3 ) then
   Fs =  36.d0*(VISKN1/DS)**2/SN/DS 
   Fs = sqrt(2.d0/3.d0 +Fs) -sqrt(Fs)
 else
   Fs = 0.79d0 
 endif  
 SETVLRB = Fs*sqrt(DS*SN)
!----------------------------------------------------------------------C 
!     End of SETVLRB group. 
!----------------------------------------------------------------------C 
! 
RETURN
END
!======================================================================C 

!======================================================================C 
REAL*8 FUNCTION VfallCheng(Res)
!----------------------------------------------------------------------C 
!
!     VfallCheng: Fall VeLocity for the Cheng's (1997) formula.
!
!----------------------------------------------------------------------C  
!USE DRSEDT
USE SEDEQ_CELL
IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
!----------------------------------------------------------------------C 
!----------------------------------------------------------------------C   
  Re = (sqrt(25.d0 +1.2d0*DSTAR**2) -5.d0)**1.5 
  !VfallCheng = Re*VISKN/DSE
  Res = Re*VISKN/DSE
!----------------------------------------------------------------------C 
!     End of VfallCheng group. 
!----------------------------------------------------------------------C 
! 
RETURN
END
!======================================================================C 

!======================================================================C 
REAL*8 FUNCTION SHIELD(REYN)
!----------------------------------------------------------------------C 
!
!     SHIELD: SHIELDs's diagram for particle bedload reactive force.
!
!----------------------------------------------------------------------C 
USE SEDEQ_CELL
IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
!----------------------------------------------------------------------C 
REAL*8 REYN
!----------------------------------------------------------------------C 
!----------------------------------------------------------------------!
!     Compute the critical tractive force from Shields diagram.
!----------------------------------------------------------------------!
 if (REYN<0.1d0) then
   SHIELD = 0.118d0*(0.1d0**(-0.973d0))
 elseif (REYN<=2.d0) then
   SHIELD = 0.118d0*(REYN**(-0.973d0))
 elseif (REYN<=4.d0) then
   SHIELD = 0.090d0*(REYN**(-0.585d0))
 elseif (REYN<=10.d0) then
   SHIELD = 0.0434d0*(REYN**(-0.119d0))
 elseif (REYN<=30.d0) then
   SHIELD = 0.0275d0*(REYN**(0.0792d0))
 elseif (REYN<=500.d0) then
   SHIELD = 0.0194d0*(REYN**(0.181d0))
 else
   SHIELD = 0.0194d0*(5.d2**(0.181d0))
 endif
!----------------------------------------------------------------------C 
!     End of SHIELD group. 
!----------------------------------------------------------------------C 
! 
RETURN
END
!======================================================================C 


!======================================================================C 
REAL*8 FUNCTION MSHIED(REYN)
!----------------------------------------------------------------------C 
!
!     MSHIED: Modified SHIElDs's diagram for particle bedload
!             reactive force.
!
!----------------------------------------------------------------------C 
USE SEDEQ_CELL
IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
!----------------------------------------------------------------------C 
REAL*8 REYN
!----------------------------------------------------------------------C 
!----------------------------------------------------------------------!
!     Compute the critical tractive force from Shields diagram.
!----------------------------------------------------------------------!
 if (REYN<=1.d2) then
   RVAL = log10(REYN)
   RFN = 0.041d0*RVAL*RVAL - 0.356d0*RVAL - 0.977d0
   MSHIED = 10.d0**RFN
 elseif (REYN<=3.d3) then
   RVAL = log10(REYN)
   RFN = 0.132d0*RVAL - 1.804d0
   MSHIED = 10.d0**RFN
 else
   MSHIED = 0.045d0
 endif
!----------------------------------------------------------------------C 
!     End of MSHIED group. 
!----------------------------------------------------------------------C 
! 
RETURN
END
!======================================================================C 

!======================================================================C 
REAL*8 FUNCTION SHIELD1(REYN)
!----------------------------------------------------------------------C 
!
!     SHIELD: SHIELDs's diagram for particle bedload reactive force.
!
!----------------------------------------------------------------------C 
USE SEDEQ_CELL
IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
!----------------------------------------------------------------------C 
REAL*8 REYN
!----------------------------------------------------------------------C 
!----------------------------------------------------------------------!
!     Compute the critical tractive force from Shields diagram.
!----------------------------------------------------------------------!
 if (REYN<6.61d0) then
   SHIELD1 = 0.1414d0*(REYN**(-0.2306d0))
 elseif (REYN<=282.84d0) then
   SHIELD1 = (1.d0 +(0.0223d0*REYN)**(2.8358d0))**0.3542/(3.0946*REYN**0.6769)
 else
   SHIELD1 = 0.045d0
 endif
!----------------------------------------------------------------------C 
!     End of SHIELD1 group. 
!----------------------------------------------------------------------C 
! 
RETURN
END
!======================================================================C 


!======================================================================C 
REAL*8 FUNCTION SHIEDR(aDSTAR)
!----------------------------------------------------------------------C 
!
!     SHIEDR: SHIElDs's diagram for critical bed-shear velocity in
!             Rijn equation.
!
!----------------------------------------------------------------------C 
USE SEDEQ_CELL
IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
!----------------------------------------------------------------------C 
REAL*8 aDSTAR
!----------------------------------------------------------------------C 
!----------------------------------------------------------------------!
!     Specific density.
!----------------------------------------------------------------------!
 RS = SN*DSE
!----------------------------------------------------------------------!
!     Compute the critical mobility parameter from Shields diagram.
!----------------------------------------------------------------------!
 if (aDSTAR<=4.d0) then
!   TETA = 0.24d0/max(aDSTAR,SMALL)
   TETA = 0.115d0/sqrt(max(aDSTAR,SMALL))
 elseif (aDSTAR<=10.d0) then
   TETA = 0.14d0/aDSTAR**0.64d0
 elseif (aDSTAR<=20.d0) then
   TETA = 0.04d0/aDSTAR**0.1d0
 elseif (aDSTAR<=150.d0) then
   TETA = 0.013d0*aDSTAR**0.29d0
 else
   TETA = 0.055d0
 endif
!----------------------------------------------------------------------!
!     Compute critical bed-shear velocity according to Shields diagram.
!----------------------------------------------------------------------!
! SHIEDR = sqrt(RS*TETA)
 SHIEDR = TETA
!----------------------------------------------------------------------C 
!     End of SHIEDR group. 
!----------------------------------------------------------------------C 
! 
RETURN
END
!======================================================================C 
