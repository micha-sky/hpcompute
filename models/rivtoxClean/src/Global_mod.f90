!==============================================================================!
!------------------------------------------------------------------------------!
!
!     MODULES CONTAINING GLOBAL VARIABLES.
!
!-------------------------INTERFACEMOD--MODULES--------------------------------!
MODULE CHARAC
 CHARACTER NAMEQN(10)*16
END MODULE

MODULE CONSTS
 REAL*8 RHOLQ,VISLQ,GRAV,PATM,PCR,PMN,PI,VKARMAN
END MODULE   

MODULE CRSECT  
 INTEGER*4 NSEC 
 INTEGER*4, POINTER :: NSCT(:,:)
 REAL*8, POINTER :: SWD(:),SHG(:),SAR(:),SPer(:)  
 CHARACTER*10, POINTER :: TSEC(:)
END MODULE   

MODULE DRLINK 
 INTEGER*4 MaxLPoints, mPoints, GRange(4), icomma
 INTEGER*4, POINTER :: IL(:), MU(:), MD(:)
END MODULE

MODULE DRNAME   
 CHARACTER*10, POINTER  :: LinkName(:)
 CHARACTER*5, POINTER   :: NodeName(:), NMU(:), NMD(:)
 CHARACTER*6, POINTER :: NLK(:,:)
END MODULE   

MODULE DRNODE 
 INTEGER*4 NGroup, MaxLinks
 INTEGER*4, POINTER :: KK(:), NGSize(:), NGRP(:), NPos(:) 
 INTEGER*4, POINTER :: LK(:,:)
END MODULE   

MODULE DRSOLV   
 INTEGER*4 iIFRIC
 REAL*8 ALPHA, BETA, PHI, THETA, rFRIC  
END MODULE   

!MODULE DRWFLD   
! REAL*8, POINTER :: HD(:), HDO(:), QH(:), QHO(:) 
! REAL*8, POINTER :: DH(:), DQ(:)
!END MODULE   

MODULE FILINX
 INTEGER*4 IRD,IWR,IPL,IRS,ISF,ITS,IGF
END MODULE

MODULE HEADNG
 CHARACTER TITLE*150,SNOTE*1500,FileName(20)*150,NDATE*10,NTIME*10,USER*80,COMP*40,DATVAL*40,TIMVAL*40
 CHARACTER DirName*80
! CHARACTER*70 BigHeader(5)
! DATA BigHeader / &
!  'RRRR    III    AAA    SSSS  TTTTT   OOO   X   X     W       W', &
!  'R   R    I    A   A  S        T    O   O   X X      W       W', &
!  'RRRR     I    AAAAA   SSS     T    O   O    X   --  W   W   W', &
!  'R  R     I    A   A      S    T    O   O   X X       W W W W ', &
!  'R   R   III   A   A  SSSS     T     OOO   X   X       W   W  '/
! CHARACTER*(*) Head1,Head2
! PARAMETER(Head1='RIVTOX-M')
! PARAMETER(Head2='R I V T O X - M')
END MODULE

MODULE LOGCLS
 INTEGER*4 NPLOTEP !--Number of output variables for each set
 LOGICAL*4 SCREEN,INITL,IRNV(50),IPLOT(50),IPLOTEP(50),RESTART,ZSTART,STEADY
 LOGICAL*4, POINTER :: IPLOTB(:,:) 
 LOGICAL LQSOUR,SDSOUR,CLSOUR,CPSOUR
 LOGICAL LQSRND,SDSRND,CLSRND,CPSRND
END MODULE

MODULE NUMBRS
 INTEGER*4 ISMALL,IBIG
 REAL*8 SMALL,BIG,ZERO,HSMALL,EPSMIN,SIGMA
END MODULE

MODULE REFERN
 INTEGER*4 IFQEP,NREF(4),IREF(4),JREF(4)
 INTEGER*4 IFQREF,IFQMV,IFQRST,ICREF(4),ICREFL(4),ICREFG(4),ICREFO(4),IHREF,NMREF
 INTEGER*4 iRSTType !--RESTART file format {ascii,binary}
 INTEGER*4 IFQEPType !--Extraction points output frequency type {1-timestep based, 2-time based}
 REAL*8 DTOUT_EP !--Extraction points output time step (used when IFQEPType=2)
END MODULE

MODULE SOLVAR
 INTEGER*4 MXDTR,ISOLVE(10),IDMEAN(10),NSTEP,IBREAK,NHDS!,NITER
 REAL*8 TIME,TIMAX,DT,DTMIN,DTMAX,DTFAC,DTR,DTOLD,TSCL,TMVL
END MODULE

MODULE RESIDL
 INTEGER*4 IRSD,NRSD,NFRSD,NDRSD(10)
 REAL*8 RSDMX(10),RSD(10)
END MODULE

MODULE INPVAR
 INTEGER*4 NROCK,NPLANT
 CHARACTER*80, DIMENSION(:), POINTER :: ROCTYP,PLNTYP
END MODULE

MODULE PLZONE
 INTEGER*4, DIMENSION(:), POINTER :: IZPL,IZ
END MODULE

MODULE PROPER
 REAL*8, DIMENSION(:), POINTER :: RHOS,POR 
END MODULE

MODULE DRCATC
 REAL*8 HMIN
 INTEGER*4, POINTER :: IFRIC(:), IDIFF(:)
 REAL*8, POINTER :: EROD(:),EVC(:),EVW(:),DFSL(:),ROUGH(:),FRIC(:)
 REAL*8, POINTER :: PCSD(:),PCSB(:),EXSD(:),EXSB(:)
END MODULE

MODULE DRSEDT
 INTEGER*4 ILOAD,ISHMODE
 REAL*8 DSE,D90,RHOSE,VSD
 REAL*8 TAY_D,TAY_E,WCM !--For cohesive
END MODULE

MODULE SPECIE
 REAL*8, POINTER :: SMDB(:),SMDL(:),TPKD(:),PCKD(:),HFLF(:)
END MODULE

MODULE DRBCLQ
 INTEGER*4 NBCL,NBTU,NBTL
 INTEGER*4, POINTER :: IBCND(:),IBCTY(:),MBCND(:) !
 INTEGER*4, POINTER :: NBCT(:,:) !
 REAL*8, POINTER :: BCV(:),BCT(:) !
END MODULE

MODULE DRBCSE
 INTEGER*4 NBCS,NBTS,NBTLS
 INTEGER*4, POINTER :: IBCNS(:),IBCTS(:),MBCS(:),NBRNS(:,:) !--(NBC)S; (2,NBTS); NBTS, NBCS = NFB
 REAL*8, POINTER :: BCVS(:),BCTS(:) !--NBTLS <= NumOfBC input lines + NumOfBC Table input lines + 1
END MODULE
  
MODULE DRBCCL
 INTEGER*4 NBCT,NBTT,NBTLT
 INTEGER*4, POINTER :: IBCNT(:),IBCDT(:),IBCTT(:),MBCT(:),NBTST(:,:) !--(NBC)S; (2,NBTS); NBTS, NBCS = NFB
 REAL*8, POINTER :: BCVT(:),BCTT(:) !--NBTLS <= NumOfBC input lines + NumOfBC Table input lines + 1
END MODULE
!      COMMON /RNBCSP/ NBCT,BCVT(LBTSE),BCTT(LBTSE),IBCNT(LBCSE), 
!     &                IBCDT(LBCSE),IBCTT(LBCSE),MBST(LBCSE),NBTT, 
!     &                NBTST(2,LBCSE) 
  
MODULE DRBCCP
 INTEGER*4 NBCR,NBTR,NBTLR
 INTEGER*4, POINTER :: IBCNR(:),IBCDR(:),IBCTR(:),MBCR(:),NBTSR(:,:) !--(NBC)S; (2,NBTS); NBTS, NBCS = NFB
 REAL*8, POINTER :: BCVR(:),BCTR(:) !--NBTLS <= NumOfBC input lines + NumOfBC Table input lines + 1
END MODULE
!      COMMON /RNBCSS/ NBCR,BCVR(LBTSE),BCTR(LBTSE),IBCNR(LBCSE), 
!     &                IBCDR(LBCSE),IBCTR(LBCSE),MBSR(LBCSE),NBTR, 
!     &                NBTSR(2,LBCSE) 
  
    
MODULE DRRAIN
 INTEGER*4 NSPREC,NRAIN !--Num of Rains, Num of Rain Table Lines
 INTEGER*4, POINTER :: IPRCV(:),KSPREC(:,:),ISPREC(:,:) !
 REAL*8, POINTER :: SPREC(:),SPREC_TM(:) !
 REAL*8, POINTER :: RAIN(:)
END MODULE

MODULE DREVAP
 INTEGER*4 NSEVAP,NEVAP !--Num of Evaps, Num of Evap Table Lines
 INTEGER*4, POINTER :: IEVPV(:),KSEVAP(:,:),ISEVAP(:,:) !--(NSEVAP),(2,NSEVAP):  NSEVAP = NumOfEvap input lines
 REAL*8, POINTER :: SEVAP(:),SEVAP_TM(:) !--(NEVAP,2),(NEVAP),(4,NSEVAP): NEVAP = 2*NumOfEvap input lines + NumOfEvap Table input lines
 REAL*8, POINTER :: EVAP(:)
END MODULE

MODULE WIND
 LOGICAL*4 fWind
 INTEGER*4 iWDCType !--Wind formula
 REAL*8 RHAIR, ANEMHEIGHT, PreW10, WindPreCoef
 INTEGER*4 NSWIND,NWIND !--Num of Winds, Num of Wind Table Lines
 INTEGER*4, POINTER :: IWINDV(:),KSWIND(:,:),WINDR(:,:) !--(NSWIND),(2,NSWIND):  NSWIND = NumOfWind input lines
 REAL*8, POINTER :: SWIND(:,:),SWIND_TM(:) !--(NWIND,2),(NWIND),(4,NSWIND): NWIND = 2*NumOfWind input lines + NumOfWind Table input lines
 REAL*8, POINTER :: WINDX(:),WINDY(:) !--(NC)
END MODULE

MODULE OUTPLT
 INTEGER*4 NPRTM,iOutInaMode,iOutDryMode,NFOut,IParOut,IParInp
 REAL*8 OutInaWSL,OutDryWSL 
 !INTEGER*4 IOutType,IOutMode
 INTEGER*4, POINTER :: IOutType(:),IOutMode(:) !--Binary output type and mode
 LOGICAL*4, POINTER :: fNodeOutputSolBin(:) !--Binary output at nodes or at cells
 CHARACTER*150, POINTER :: IOutFile(:) !--Binary output file name
 REAL*8, POINTER :: TMIN_OUT(:),TMAX_OUT(:),DT_OUT(:), TOuts(:)
 INTEGER*4, POINTER :: NSolFrames(:,:)  
 LOGICAL*4, POINTER :: fFirstSolWrite(:), FOuts(:)
 REAL*8, POINTER :: SolStartT(:,:),SolEndT(:,:)
 REAL*8, POINTER :: PRTM(:)
END MODULE

MODULE PVSPNL
 PARAMETER(MNVARS=50) !--Max Number of output variables
 INTEGER*4 NVARS !--Number of output variables
 CHARACTER*10 FMTS(3)
 CHARACTER*3 PNS(MNVARS)
 CHARACTER*3 PNS2(2,MNVARS)
 CHARACTER*40 PNL(MNVARS)
 CHARACTER*10 PNU(MNVARS)
 CHARACTER*50 PNLU(MNVARS)
 INTEGER*4 PNF(MNVARS) !--Flag indicationg that variable is calculatable (and its value can be contained in or derived from RESTART): 0-not calc, 1-calc.
 INTEGER*4 PNDim(MNVARS) !--Dimensionality of variable
 DATA FMTS /'bouss2d','sms','immsp'/
END MODULE

MODULE EP_MOD
 INTEGER*4 EPN,iFUEP
 CHARACTER*150 EPFName
! INTEGER*4 IFQEP
! INTEGER*4 IFQEPType !--Extraction points output frequency type {1-timestep based, 2-time based}
! REAL*8 DTOUT_EP !--Extraction points output time step (used when IFQEPType=2)
 INTEGER*4, POINTER :: EPIJ(:,:)
 REAL*8, POINTER :: EPXY(:,:)
 CHARACTER*30, POINTER :: EPName(:)
END MODULE

!--------------------------------GLOBAL--MODULES-------------------------------!
MODULE FLDINX
 INTEGER*4 NNode,NLink,NPoint,NXP,nBound
 INTEGER*4, POINTER :: ILK(:),JPT(:),ND(:,:),NSX(:),NSY(:),IXP(:)
END MODULE

!MODULE GEOMTR
! REAL*8 XX0,YY0,XX1,YY1
! REAL*8, DIMENSION(:), POINTER :: X,Y, DX,DY, DXI,DYI, VOL, AFX,AFY, XBD,YBD
! REAL*8, DIMENSION(:), POINTER :: CX,CY !--X,Y for all cells
!END MODULE

MODULE DRWFLD   
 REAL*8 FNew, FPrev, FMid
 REAL*8, POINTER :: HR(:), HRO(:), QH(:), QHO(:), YH(:)  
 REAL*8, POINTER :: UH(:), UHO(:)
 REAL*8, POINTER :: DH(:), DQ(:) 
 REAL*8, POINTER :: RH(:), RQ(:)  ! residuals of the continuity equation and momentum conservation equation 
END MODULE   

MODULE DRSPFL
 REAL*8, DIMENSION(:), POINTER :: SH,SHO,SEQ
 REAL*8, DIMENSION(:), POINTER :: CH,CHO,CSH,CSHO,CB,CBO
 REAL*8, DIMENSION(:), POINTER :: BLX,BLY,SLX,SLY
END MODULE

MODULE LANDSF
 REAL*8, DIMENSION(:), POINTER :: TSF,TSFO,ZSD,ZSDO
END MODULE

MODULE DIFF
 LOGICAL*4 fDIFF
 REAL*8, POINTER :: DiffC(:)
END MODULE

MODULE LINEQN
 INTEGER*4 NDIG,NCOUP,NEL,nALC  
 REAL*8, POINTER :: CRA(:,:), CRH(:,:), DY(:) 
 REAL*8, POINTER :: ALU(:,:), BLU(:) !--(5*LEQN,LCOUP),(LCOUP)
 INTEGER*4, POINTER :: IA(:),JA(:) 
 REAL*8, POINTER :: ALC(:), BLC(:)  
 INTEGER*4, POINTER :: iEq(:)
END MODULE

MODULE POINTS
 INTEGER*4, POINTER :: MSEC(:)
 REAL*8, POINTER :: X(:), ROU(:), FLQ(:) 
 REAL*8, POINTER :: AW(:), AWO(:), PER(:), PERO(:), DPER(:), DAW(:), BW(:)
END MODULE

MODULE LNDAMS
 INTEGER*4 NDAM
 INTEGER*4, POINTER :: IUPDAM(:), IDNDAM(:)
END MODULE

MODULE WORK
 INTEGER*4 NRWORK,NIWORK
 REAL*8, POINTER :: RWORK(:) !--(2*(4+LDIG)*LCOUP)
 INTEGER*4, POINTER :: IWORK(:) !--(2+2*(2+LDIG)*LCOUP)
 REAL*8, DIMENSION(:), POINTER :: VV1,VV2,VV3
 REAL*8, DIMENSION(:), POINTER :: PA1,PA2,PA3,PA4 !--Shared pointer arrays
END MODULE

MODULE RVDISH
 REAL*8, DIMENSION(:), POINTER :: UC,VC
END MODULE

MODULE LQSORC
 INTEGER*4 NSRCLQ,NCLQ 
 INTEGER*4, POINTER :: ISLQV(:),KSRCLQ(:,:),ISRCLQ(:,:)  
 REAL*8, POINTER :: SlkLQ(:)
 REAL*8, POINTER :: SRCLQ(:),SRCLQ_TM(:) 
END MODULE

MODULE SPSOUR  
 INTEGER*4 NSRCSP(3),NCSP(3) 
 INTEGER*4, POINTER :: ISSPV(:,:),KSRCSP(:,:,:),ISRCSP(:,:,:),ISPTP(:,:)  
 INTEGER*4, POINTER :: ISVCL(:),KSRCCL(:,:),ISRCCL(:,:),ISTYCL(:) 
 INTEGER*4, POINTER :: ISVCS(:),KSRCCS(:,:),ISRCCS(:,:),ISTYCS(:) 
 INTEGER*4, POINTER :: ISVSH(:),KSRCSH(:,:),ISRCSH(:,:),ISTYSH(:) 
 REAL*8, POINTER :: SRCSP1(:),SRCSP_TM1(:) 
 REAL*8, POINTER :: SRCSP2(:),SRCSP_TM2(:) 
 REAL*8, POINTER :: SRCSP3(:),SRCSP_TM3(:)  
 REAL*8, POINTER :: SlkCL(:),SlkCS(:),SlkSH(:)
END MODULE

MODULE LQNODE
 INTEGER*4 NDSRLQ,NDLQ 
 INTEGER*4, POINTER :: INDLQV(:),KNDSLQ(:,:),INDSLQ(:,:)  
 REAL*8, POINTER :: SNDLQ(:),SNDLQ_TM(:),TEMP_AR(:) 
 REAL*8, POINTER :: SRCndLQ(:)
END MODULE

MODULE SPNODE
 INTEGER*4 NDSRSP(3),NDSP(3) 
 INTEGER*4, POINTER :: INDSPV(:,:),KNDSSP(:,:,:),INDSSP(:,:,:),INDTP(:,:)   
 INTEGER*4, POINTER :: INDVCL(:),KNDSCL(:,:),INDSCL(:,:), INDTYCL(:) 
 INTEGER*4, POINTER :: INDVCS(:),KNDSCS(:,:),INDSCS(:,:), INDTYCS(:) 
 INTEGER*4, POINTER :: INDVSH(:),KNDSSH(:,:),INDSSH(:,:), INDTYSH(:)   
 REAL*8, POINTER :: SndCL(:),SndCS(:),SndSH(:)
 REAL*8, POINTER :: SNDSP1(:),SNDSP_TM1(:) 
 REAL*8, POINTER :: SNDSP2(:),SNDSP_TM2(:) 
 REAL*8, POINTER :: SNDSP3(:),SNDSP_TM3(:) 
 REAL*8, POINTER :: SRCndSP(:),SPndConc(:) 
    END MODULE

!------------------------------------------------------------------------------!
!======================================================================C   
!     Add to solution files    
!======================================================================C   
MODULE ArrayPointers
 TYPE T1DR4
  REAL*4, POINTER :: V(:)=>NULL()
 END TYPE

 TYPE T2DR4
  REAL*4, POINTER :: V(:,:)=>NULL()
 END TYPE

 TYPE T3DR4
  REAL*4, POINTER :: V(:,:,:)=>NULL()
 END TYPE

 TYPE T4DR4
  REAL*4, POINTER :: V(:,:,:,:)=>NULL()
 END TYPE

 TYPE T5DR4
  REAL*4, POINTER :: V(:,:,:,:,:)=>NULL()
 END TYPE

 TYPE T1DR8
  REAL*8, POINTER :: V(:)=>NULL()
 END TYPE

 TYPE T2DR8
  REAL*8, POINTER :: V(:,:)=>NULL()
 END TYPE

 TYPE T3DR8
  REAL*8, POINTER :: V(:,:,:)=>NULL()
 END TYPE

 TYPE T4DR8
  REAL*8, POINTER :: V(:,:,:,:)=>NULL()
 END TYPE

 TYPE T5DR8
  REAL*8, POINTER :: V(:,:,:,:,:)=>NULL()
 END TYPE

 TYPE T1DI4
  INTEGER*4, POINTER :: V(:)=>NULL()
 END TYPE

 TYPE T2DI4
  INTEGER*4, POINTER :: V(:,:)=>NULL()
 END TYPE

 TYPE T3DI4
  INTEGER*4, POINTER :: V(:,:,:)=>NULL()
 END TYPE

 TYPE T4DI4
  INTEGER*4, POINTER :: V(:,:,:,:)=>NULL()
 END TYPE

 TYPE T5DI4
  INTEGER*4, POINTER :: V(:,:,:,:,:)=>NULL()
 END TYPE
 END MODULE
    
 MODULE WTSOL
 USE ArrayPointers
 INTEGER*4 NWTSol !--Num of Solution files
 INTEGER*4 iFUWTSol !--Starting index of solution file units
 CHARACTER*150, POINTER :: WTSolFN(:) !--Grid names for solutions
 INTEGER*4, POINTER :: IWTSolIFU(:),IWTSolVR(:),IWTSolFT(:)
 REAL*8, POINTER :: WTSolCU(:) !--Solution unit values
 LOGICAL*4, POINTER :: IWTSolUpd(:) !--Mark solutions that were updated
 INTEGER*4, POINTER :: IWTSolExtrMode0(:),IWTSolExtrMode1(:) !--Time-Extrapolation mode before the start and after the end of solution
 LOGICAL*4, POINTER :: fWTSolExtr1(:) !--Mark solutions that were time-extrapolated after the end of solution
 INTEGER*4, POINTER :: WTSolNDim(:),WTSolNN(:)
 INTEGER*4, POINTER :: iWTSolFr(:),WTSolNFr(:) !--Solution current frame indices;Solution nums of frames;
 INTEGER*4, POINTER :: WTSoliPar(:,:) !--iPar values for solution files

 TYPE(T1DR8), POINTER :: WTSolTimes(:) !--Frame Time arrays for all solutions

 INTEGER*4, POINTER :: WTSolVDim(:) !--1D\2D mode for all solutions (used by WTSolV0)
 TYPE(T1DR8), POINTER :: WTSolV0(:,:),WTSolV1(:,:) !--Frame arrays for all solutions


 INTEGER*4, POINTER :: IWTSolGridInd(:) !--Solution indices to corresponding grid (see below)

 INTEGER*4 NWTSolGrids !--Num of different Solution file grids
 CHARACTER*150, POINTER :: WTSolGridFN(:) !--Grid names for solutions
 INTEGER*4, POINTER :: IWTSolGridRI(:),IWTSolGridFmt(:) !--Grid RI types {0-model grid, 1-irreg, 2-rect}; file types {1 - ADCIRC 14 , 2 - xyz};
 INTEGER*4, POINTER :: WTSolGridNN(:),WTSolGridNX(:),WTSolGridNY(:),WTSolGridNTR(:) !--Grid NN,NX,NY, NTR

 TYPE(T1DR8), POINTER :: WTSolGridX(:),WTSolGridY(:) !--X,Y arrays for all grids
 TYPE(T2DI4), POINTER :: WTSolGridTri(:) !--Tri arrays for all grids
 TYPE(T1DI4), POINTER :: WTSolGridIWK(:) !--IWK arrays for all grids
 TYPE(T1DR8), POINTER :: WTSolGridRWK(:) !--RWK arrays for all grids
END MODULE



MODULE RW_MOD
 !-----MUST BE INITIALIZED BY CALLING PROGRAM BEFORE USING THE SUBROUTINES OF THIS MODULE!-----!
 INTEGER*4 :: RW_IFUTemp=1001 !--FileUnit for temporary file operations
 INTEGER*4 :: RW_IFUEcho=-1   !--FileUnit to echo information to  
END MODULE 
!======================================================================C   

MODULE INTERFACEMOD
 USE NUMBRS
 USE CONSTS
 USE FILINX
 USE HEADNG
 USE LOGCLS
 USE CHARAC
 USE REFERN
 USE SOLVAR
 USE DRSOLV
 USE RESIDL
 USE INPVAR
 USE PROPER
 USE DRCATC 
 USE DRSEDT 
 USE SPECIE
 USE DRBCLQ
 USE DRBCSE
 USE DRRAIN
 USE DREVAP
! USE WIND
 USE OUTPLT
END MODULE
!------------------------------------------------------------------------------!
MODULE mGLOBAL
 USE INTERFACEMOD
 USE FLDINX
! USE GEOMTR
 USE DRWFLD
! USE RVSPFL
 USE LANDSF
! USE DIFF 
 USE POINTS
 USE LINEQN
 USE WORK
 USE RVDISH
END MODULE
!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
!==============================================================================!
