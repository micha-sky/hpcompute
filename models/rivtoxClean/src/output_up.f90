
    
    
    
    
    
    
    
    
    
    !======================================================================C 
!
SUBROUTINE WriteSolutionFiles(kFOut)
!----------------------------------------------------------------------C 
!
!     Writes requested variables to field time-series files.
!     Use VV1N,VV2N as a temporary storage.
!
!----------------------------------------------------------------------C 
!USE GLOBAL  
USE DRLINK
USE LOGCLS
USE OUTPLT
USE PVSPNL  
USE WORK 
USE SOLVAR
IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
!----------------------------------------------------------------------C 
INTEGER*4 kFOut !--index of solution files set--
!----------------------------------------------------------------------C 
CHARACTER FName*50, DimStr*3, FOutSuffStr*2
!------------------------------------------------------------------------------C
!----------------------------------------------------------------------!
!     Create a filename numbering suffix.
!----------------------------------------------------------------------!
 write(FOutSuffStr,'(I2.2)') kFOut

!----------------------------------------------------------------------!
!     Use VV1N,VV2N as a temporary storage !!!
!----------------------------------------------------------------------! 
 if( .not.associated(VV1) ) then  
   allocate(VV1(mPoints)); VV1 = 0.d0;  
 endif    
 if( .not.associated(VV2) ) then  
   allocate(VV2(mPoints)); VV2 = 0.d0;  
 endif    
 do iv=1,NVARS
   if (.not. IPLOTB(iv,kFOut)) CYCLE
!----------------------------------------------------------------------!
!     Create file dimension.
!----------------------------------------------------------------------!
   DimStr='dat' !--default extension
   select case (IOutType(kFOut))
   case (2) !--bouss2d--
     select case (iv)
     case (1:24)
       DimStr='eta'
     case (25)
       DimStr='uv'
     end select
   case (3,4) !--sms,immsp--
     DimStr='dat'
   end select
!----------------------------------------------------------------------!
!     Create file name.
!----------------------------------------------------------------------!
   FName=trim(FMTS(IOutType(kFOut)-1))//'_'//trim(PNS(iv))//'_'//FOutSuffStr//'.'//trim(DimStr)//CHAR(0)
!----------------------------------------------------------------------!
!     Fill arrays.
!----------------------------------------------------------------------!
   call FillVarValuesAtNodes(iv)
!----------------------------------------------------------------------!
!     Write.
!----------------------------------------------------------------------!
   call WriteSolutionFilesProc(kFOut,iv,FName, VV1,VV2, PNDim(iv), SolStartT(iv,kFOut),SolEndT(iv,kFOut),NSolFrames(iv,kFOut))
 enddo

!----------------------------------------------------------------------!
!     Update params after all variables are written.
!     Do this if frame was really written.
!     Inside write routines they are initialized only once at the 
!     very start. And then left unchanged because otherwise
!     it creates problems with multiple variables output.
!----------------------------------------------------------------------!
 do iv=1,NVARS
   if (.not. IPLOTB(iv,kFOut)) CYCLE
   if (TIME>SolEndT(iv,kFOut)) then
     SolEndT(iv,kFOut)=TIME
     NSolFrames(iv,kFOut)=NSolFrames(iv,kFOut)+1
   endif
 enddo

 fFirstSolWrite(kFOut)=.false.
!----------------------------------------------------------------------C 
RETURN
END

!======================================================================C 
!======================================================================C 
!  
SUBROUTINE WriteSolutionFilesProc(kFOut,iv,FName, V1,V2, NDim, StartT,EndT,NFrames)
!----------------------------------------------------------------------C 
!USE GLOBAL  
USE LOGCLS
USE OUTPLT 
USE DRLINK 
USE SOLVAR  
USE PVSPNL
IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
!----------------------------------------------------------------------C  
Interface    
!  subroutine wSMSGenericBinary(FName,fRewrite,fFirstWrite,ParNameL,frTime,NDim,NN,V1,V2,IObT,NC,ITU,StartT,EndT,NFrames)  
!    CHARACTER*(*) FName
!    LOGICAL*4 fRewrite,fFirstWrite
!    CHARACTER*(*) ParNameL
!    REAL*8 frTime
!    INTEGER*4 NDim, NN
!    REAL*8 V1(NN),V2(NN)
!    INTEGER*4 IObT,NC,ITU
!    REAL*8 StartT,EndT
!    INTEGER*4 NFrames
!  end subroutine wSMSGenericBinary
end Interface
!----------------------------------------------------------------------C 
INTEGER*4 kFOut,iv !--index of solution files set and of variable
CHARACTER*(*) FName
REAL*8 V1(*),V2(*)
INTEGER*4 NDim
REAL*8 StartT,EndT
INTEGER*4 NFrames
LOGICAL*4 fRewrite,fFirstWrite
REAL*8 DelX,DelY,Grid_Orientation
!----------------------------------------------------------------------C 
 fRewrite=.not.RESTART
 fFirstWrite=fFirstSolWrite(kFOut)
 select case (IOutType(kFOut))
!   DelX=(X(IFLD+1) - X(2))/(IFLD-1)
!   DelY=(Y(JFLD+1) - Y(2))/(JFLD-1)
!   Grid_Orientation=0.d0
!   call wBouss2DBinary(FName, fRewrite,fFirstWrite, PNLU(iv), TIME, NDim, &
!    IFLD,JFLD, DelX, DelY, Grid_Orientation, X(2), Y(2), DT_OUT(kFOut), &
!    V1,V2, StartT,EndT,NFrames)
 case(3)
   call wSMSGenericBinary(FName, fRewrite,fFirstWrite, PNLU(iv), TIME, NDim, mPoints, V1,V2, 5,mPoints, 2,StartT,EndT,NFrames)
 case(4)
   call wIMMSPBinary(FName, fRewrite,fFirstWrite, PNLU(iv), TIME, NDim, mPoints, V1,V2, StartT,EndT,NFrames)
 end select
!----------------------------------------------------------------------C 
    RETURN
    END

!======================================================================C 
!======================================================================C 
!======================================================================C 
!
    SUBROUTINE FillVarArray(icn,IVAR,ICF,VARF,CHRF)
!----------------------------------------------------------------------C 
!
!     FillVarArray: Fill Variables Array.
!
!     Subroutine fills array with variables values that have to be outputted, at cells or nodes.
!
!----------------------------------------------------------------------C  
USE DRWFLD  
USE DRSPFL 
USE PVSPNL 
USE LANDSF 
USE POINTS
!USE GLOBAL
IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
!----------------------------------------------------------------------C 
INTEGER*4 icn,ICF
LOGICAL IVAR(MNVARS)
REAL*8 VARF(MNVARS)
CHARACTER*3 CHRF(MNVARS)
!----------------------------------------------------------------------C 
! IZRC = IZR(icn)
 ICF = 0
 do iv=1,NVARS
   if (IVAR(iv)) then
     ICF = ICF +1
     select case (iv)
     case(1)
       VARF(ICF) = AW(icn)
     case(2)
       VARF(ICF) = ROU(icn)

     case(3)
       VARF(ICF) = TSF(icn)
     case(4)
       VARF(ICF) = TSF(icn) + HR(icn)
     case(5)
       VARF(ICF) = HR(icn)
     case(6)
       VARF(ICF) = QH(icn)
     case(7)
       VARF(ICF) = ZSD(icn)

     case(8)
       VARF(ICF) = SH(icn)
     case(9)
       VARF(ICF) = SEQ(icn)
     case(10)
       VARF(ICF) = CH(icn)
     case(11)
       VARF(ICF) = CSH(icn)
     case(12)
       VARF(ICF) = CB(icn)

     case(13)
       VARF(ICF) = SDOWN(icn)
     case(14)
       VARF(ICF) = RSUP(icn)

     case(15)
       VARF(ICF) = ZSD(icn)
     case(16)
!       SDOWN(icn)
     case(17)
!       RSUP(icn)
     end select

   endif
 enddo
!----------------------------------------------------------------------C 
    RETURN
    END  
    
!======================================================================C   
!    
SUBROUTINE FillVarValuesAtNodes(iv)
!----------------------------------------------------------------------C 
!
!     FillVarValuesAtNodes: Fill Variable Values At all Nodes.
!
!     Subroutine fills VV1N,VV2N arrays with output variable values.
!
!----------------------------------------------------------------------C 
!USE GLOBAL 
USE DRLINK
USE WORK
IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
!----------------------------------------------------------------------C 
INTEGER*4 iv
!----------------------------------------------------------------------C 
 call PointToVarArray(iv)

!----------------------------------------------------------------------!
!     1D.
!----------------------------------------------------------------------!
 select case (iv)
 case (1:3, 5:12) !--No adjustments
   VV1(1:mPoints) = PA1(1:mPoints)
   RETURN  
 case (4)  
   VV1(1:mPoints) = PA1(1:mPoints) !!+PA2(1:mPoints) 
   return  
 case (13)  
   do n=1,mPoints  
     VV1(n) = SDOWN(n)  
   enddo  
   return
 case (14)  
   do n=1,mPoints  
     VV1(n) = RSUP(n)  
   enddo  
   return
 end select

!----------------------------------------------------------------------!
!     2D.
!----------------------------------------------------------------------!
 select case (iv)
 case (25, 26, 28, 35, 39) !--No adjustments
   VV1(1:mPoints) = PA1(1:mPoints)
   VV2(1:mPoints) = PA2(1:mPoints)
   RETURN
 case (32) !-SHR
   VV1(1:mPoints) = PA1(1:mPoints)*PA2(1:mPoints)
   RETURN
 end select

!----------------------------------------------------------------------!
!     Variables defined in cells, more >2D.
!----------------------------------------------------------------------!
 select case (iv)
 case (30) !--Total load
   VV1(1:mPoints) = PA1(1:mPoints)+PA2(1:mPoints)
   VV2(1:mPoints) = PA3(1:mPoints)+PA4(1:mPoints)
   RETURN
 case (31) !--Total load Magnitude
   VV1(1:mPoints) = PA1(1:mPoints)+PA2(1:mPoints)
   VV2(1:mPoints) = PA3(1:mPoints)+PA4(1:mPoints)
   VV1(1:mPoints) = sqrt(VV1(1:mPoints)*VV1(1:mPoints) + VV2(1:mPoints)*VV2(1:mPoints))
   RETURN
 end select
!----------------------------------------------------------------------C 
    RETURN
    END 
    
!======================================================================C   
!    
    
LOGICAL*4 FUNCTION MAKEDIRQQ(path)
use iso_c_binding

  interface
    function mkdir(path,mode) bind(c,name="mkdir")
      use iso_c_binding
      integer(c_int) :: mkdir
      character(kind=c_char,len=1) :: path(*)
      integer(c_int16_t), value :: mode
    end function mkdir
  end interface
  
  integer i
  character(*) path
  
  
  i = mkdir(path, int(o'772',c_int16_t))
  if (i == 0) then
     MAKEDIRQQ = .true.
  else
     MAKEDIRQQ = .false.
  endif
  
  
END FUNCTION    

! =======================================
    
SUBROUTINE InitExtractionPoints
!----------------------------------------------------------------------C 
!USE GLOBAL 
USE PVSPNL 
USE LOGCLS 
USE EP_MOD
IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
!----------------------------------------------------------------------C 
CHARACTER*500 str,strH,FName
CHARACTER*20 FStatus
LOGICAL*4 fFExist,fFlag,fRes,fDupl, MAKEDIRQQ

CHARACTER*2 frmt

!----------------------------------------------------------------------C 
!-----------------------------------------------------------------------! 
!     Calc number of extraction point variables.
!-----------------------------------------------------------------------! 
 NPLOTEP = 0
 do iv=1,NVARS
   if (IPLOTEP(iv)) NPLOTEP = NPLOTEP +1
 enddo
!----------------------------------------------------------------------! 
!     Create Header string.
!----------------------------------------------------------------------! 
 strH='       Time'
 lnH=len_trim(strH)
 do iv=1,NVARS
   if (.not. IPLOTEP(iv)) CYCLE
   
   if (PNDim(iv)==2) then !--2D
     do j=1,2
       ln=len_trim(PNS2(j,iv))
       ICNT=14-ln
       
       frmt = '  '
       write(frmt,'(I0)') ICNT
       write(strH(lnH+1:),'('//TRIM(frmt)//'(X),A)') PNS2(j,iv)(1:ln) !! <ICNT>(X),A
       
       lnH=len_trim(strH)
     enddo
   else !--Scalar
     ln=len_trim(PNS(iv))
     ICNT=14-ln

     frmt = '  '
     write(frmt,'(I0)') ICNT
     write(strH(lnH+1:),'('//TRIM(frmt)//'(X),A)') PNS(iv)(1:ln) !! <ICNT>(X),A
   endif

   lnH=len_trim(strH)
 enddo

!----------------------------------------------------------------------!
!     Output.
!----------------------------------------------------------------------!
 do kl=1,EPN

!----------------------------------------------------------------------!
!     If outside the domain issue a warning.
!----------------------------------------------------------------------!
   if ((EPIJ(kl,1)==0).or.(EPIJ(kl,2)==0)) then
     if (EPName(kl)/='') then
       write(IWR,'(A,I0,2A)') & 
       'INPUT WARNING! - Coordinates of extraction point No.',kl,' ('//trim(EPName(kl)),') located outside of the domain!'
       write(*,'(A,I0,2A)') & 
       'INPUT WARNING! - Coordinates of extraction point No.',kl,' ('//trim(EPName(kl)),') located outside of the domain!'
     else
       write(IWR,'(A,I0,A)') 'INPUT WARNING! - Coordinates of extraction point No.',kl,' located outside of the domain!'
       write(*,'(A,I0,A)') 'INPUT WARNING! - Coordinates of extraction point No.',kl,' located outside of the domain!'
     endif
   endif

!----------------------------------------------------------------------!
!     If gauge name exists use it to create outloc filename.
!----------------------------------------------------------------------!
   if (EPName(kl)/='') then
     fDupl=.false.
     ln=len_trim(EPName(kl))
!----------------------------------------------------------------------!
!     Check for duplicate names.
!     Rename outloc file by adding ID if needed.
!     Issue a warning if needed.
!----------------------------------------------------------------------!
     do k=1,EPN
       if (k==kl) CYCLE
       if (EPName(k)==EPName(kl)) then
         FName=EPName(kl)(1:ln)//'_ID0000.txt'//CHAR(0)
         write(FName(ln+4:ln+7),'(I4.4)') kl
         write(IWR,'(/3A,I0,A,I0)') & 
         'INPUT WARNING! - Duplicate Gauge Name "',trim(EPName(kl)),'" for extraction points No.',kl,' and No.',k
         write(IWR,'(A)') 'HINT! - Use quotation marks "" to specify gauge names with delimiters!'
         write(IWR,'(3A)') 'Outloc file renamed to ',trim(FName),'!'
         write(*,'(/3A,I0,A,I0)') & 
         'INPUT WARNING! - Duplicate Gauge Name "',trim(EPName(kl)),'" for extraction points No.',kl,' and No.',k
         write(*,'(A)') 'HINT! - Use quotation marks "" to specify gauge names with delimiters!'
         write(*,'(3A)') 'Outloc file renamed to ',trim(FName),'!'
         fDupl=.true.
         EXIT
       endif
     enddo
     if (.not. fDupl) FName=trim(EPName(kl))//'.txt'//CHAR(0)
!----------------------------------------------------------------------!
!     Otherwise use its ID.
!----------------------------------------------------------------------!
   else
     FName='outloc.0000'//CHAR(0)
     write(FName(8:),'(I4.4)') kl
   endif

   fFExist=.false.
   INQUIRE(FILE='outlocs/'//FName, EXIST = fFExist)
   fFlag=.not.(RESTART .and. fFExist)
   if (fFlag) then
     FStatus='UNKNOWN'
     fRes = MAKEDIRQQ('outlocs'//CHAR(0))
   else
     FStatus='UNKNOWN'
   endif
   OPEN(UNIT=iFUEP+kl-1, FILE='outlocs/'//FName, STATUS=FStatus, &
           ACCESS='STREAM', FORM='FORMATTED') !OPENMARK

   if (fFlag) then
     write(iFUEP+kl-1,'(A)') '##########################################################'
!     write(iFUEP+kl-1,'(A,2(1PG16.9,A))') '# (x,y) = (',EPXY(kl,1),',',EPXY(kl,2),')'
     if ((EPIJ(kl,1)/=0).and.(EPIJ(kl,2)/=0)) then
       write(iFUEP+kl-1,'(A,2(I0,A))') '# (iLink,jPoint) = (',EPIJ(kl,1),',',EPIJ(kl,2),')'
       write(iFUEP+kl-1,'(A)') '##########################################################'
       write(iFUEP+kl-1,'(A)') trim(strH)
     else
       write(iFUEP+kl-1,'(A)') '# WARNING! - Coordinates outside of the domain!'
       write(iFUEP+kl-1,'(A)') '##########################################################'
     endif
   endif
 enddo
 write(*,'(A)') 'Done'
!----------------------------------------------------------------------!
!     Format Statements.
!----------------------------------------------------------------------!
!!!  1111 FORMAT(*(X),A)
!----------------------------------------------------------------------C 
    RETURN
    END

!======================================================================C 
!
INTEGER*4 FUNCTION len_trimm(str)
!----------------------------------------------------------------------C 
CHARACTER*(*) str
INTEGER*4 n
!----------------------------------------------------------------------C 
 n = len(str)
 do while(n >= 1)
   if (str(n:n) == ' ' .or. ichar(str(n:n)) == 0) then
     n=n-1
   else
     EXIT
   endif
 enddo

 len_trimm=n
!----------------------------------------------------------------------C 
    RETURN
    END FUNCTION 
    
!======================================================================C   
!    
SUBROUTINE FindVarName(str,iv)
!----------------------------------------------------------------------C 
!
!     FindVarName: Set Flag Corresponding to Variable Name to True.
!
!----------------------------------------------------------------------C 
IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
!----------------------------------------------------------------------C 
CHARACTER*(*) str
INTEGER*4 iv
!----------------------------------------------------------------------C 
 iv=0
 ln=len_trim(str)
 if (str(1:13) == 'cross-section'  )  iv=1
 if (str(1:12) == 'manning coef'   )  iv=2

 if (str(1:12) == 'thalweg elev'   )  iv=3
 if (str(1:13) == 'water surface'  )  iv=4
 if (str(1:11) == 'water depth'    )  iv=5
 if (str(1:14) == 'flow discharge' )  iv=6
 if (str(1:15) == 'top-layer depth')  iv=7

 if (str(1:14) == 'sediment conce' )  iv=8
 if (str(1:14) == 'sediment equil' )  iv=9
 if (str(1:14) == 'sol species co' )  iv=10
 if (str(1:15) == 'part species co')  iv=11
 if (str(1:14) == 'bot species co' )  iv=12

 if (str(1:15) == 'deposition rate')  iv=13
 if (str(1:12) == 'erosion rate'   )  iv=14
 !if (str(1:15) == 'top-layer depth') iv=15
!----------------------------------------------------------------------C 
    RETURN
    END

!======================================================================C   
!    
SUBROUTINE PointToVarArray(iv)
!----------------------------------------------------------------------C 
!
!     PointToVarArray: Assign Pointers to Variable Array or Arrays.
!
!----------------------------------------------------------------------C   
USE WORK
USE DRWFLD  
USE DRSPFL 
USE PVSPNL 
USE LANDSF 
USE POINTS
!USE GLOBAL
IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
!----------------------------------------------------------------------C 
INTEGER*4 iv
!----------------------------------------------------------------------C 
 PA1=>NULL()
 PA2=>NULL()
 PA3=>NULL()
 PA4=>NULL()
 select case (iv)
 case (1)
   PA1=>AW
 case (2)
   PA1=>ROU
 case (3)
   PA1=>TSF
 case (4)
!   PA2=>TSF 
   PA1=>YH
 case (5)
   PA1=>HR
 case (6)
   PA1=>QH
 case (7)
   PA1=>ZSD
 case (8)
   PA1=>SH
 case (9)
   PA1=>SEQ
 case (10)
   PA1=>CH
 case (11)
   PA1=>CSH
 case (12)
   PA1=>CB
 end select
!----------------------------------------------------------------------C 
    RETURN
    END 
     
!======================================================================C   
!    
SUBROUTINE WriteExtractionPoints
!----------------------------------------------------------------------C 
!
!     Writes requested variables at extraction points.
!
!----------------------------------------------------------------------C 
!USE GLOBAL   
USE FLDINX
USE EP_MOD
USE LOGCLS
USE PVSPNL 
USE SOLVAR
IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
!----------------------------------------------------------------------C 
REAL*8 F(MNVARS)
CHARACTER*3 CHRF(MNVARS)

CHARACTER*10 frmt
!----------------------------------------------------------------------C 
 if (NPLOTEP==0) RETURN
!----------------------------------------------------------------------!
!     Print field variables by node number to the files.
!----------------------------------------------------------------------!
 do kl=1,EPN
   if ((EPIJ(kl,1)==0).or.(EPIJ(kl,2)==0)) CYCLE !--Don't output if outside of domain
   N=ND(EPIJ(kl,1),EPIJ(kl,2))

   call FillVarArray(N,IPLOTEP,ICNT,F,CHRF)

   do k=1,ICNT
     if (abs(F(k))<=1.d-99) F(k)=0.d0
   enddo

   frmt(1:10) = ' '
   write(frmt,'(I0)') ICNT
   write(iFUEP+kl-1,'(1PG17.10,'//TRIM(frmt)//'(1PG14.6))') TIME/TSCL, F(1:ICNT)
 enddo
!----------------------------------------------------------------------!
!     Format Statements.
!----------------------------------------------------------------------!
!!! 1111 FORMAT(1PG17.10,<ICNT>(1PG14.6))
!----------------------------------------------------------------------C 
    RETURN
    END   
    
!======================================================================C   
!    
SUBROUTINE LoadWaterSolutionsFrames
!----------------------------------------------------------------------C 
!
! LoadWaterSolutionsFrames: Load Frames for each Water Solution.
!
!----------------------------------------------------------------------C  
USE WTSOL  
USE SOLVAR
!USE GLOBAL
IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
!----------------------------------------------------------------------C 
!----------------------------------------------------------------------!
!     Simply read all needed frames except 2 last.
!----------------------------------------------------------------------!
 do i=1,NWTSol
   if (iWTSolFr(i)<=WTSolNFr(i)-1) then
     call FindTimeInterval(TIME,WTSolNFr(i),WTSolTimes(i)%V,iWTSolFr(i),iWTSolFr1)
     do k=1,(iWTSolFr1-iWTSolFr(i))-2
       iWTSolFr(i)=iWTSolFr(i)+1
       call LoadWaterSolutionFrame(i,.false.)
     enddo
   endif
 enddo
!----------------------------------------------------------------------!
!     For 2 last frames do also: a) convert to model grid b) make specific
!     adjustments.
!----------------------------------------------------------------------!
 do kk=1,2
   do i=1,NWTSol
     IWTSolUpd(i)=.false.
     if (iWTSolFr(i)<=WTSolNFr(i)-1) then
       call FindTimeInterval(TIME,WTSolNFr(i),WTSolTimes(i)%V,iWTSolFr(i),iWTSolFr1)
       do k=1,min(1,iWTSolFr1-iWTSolFr(i))
         iWTSolFr(i) = iWTSolFr(i)+1
         call LoadWaterSolutionFrame(i,.true.)
       enddo
     endif
   enddo
   call AdjustSpecificWaterSolutions
 enddo
!----------------------------------------------------------------------C 
    RETURN
    END   
    
!======================================================================C  
!    
SUBROUTINE CalcWaterSolutionsAtTime
!----------------------------------------------------------------------C 
!
! CalcWaterSolutionsAtTime: Calculates each Water Solution At the Time
!                           by interpolating between 2 fields.
!
!----------------------------------------------------------------------C  
USE DRWFLD 
USE DRLINK
USE WORK  
USE WTSOL  
USE SOLVAR
!USE GLOBAL
IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
!----------------------------------------------------------------------C 
LOGICAl*4 fVel,fDisch
!----------------------------------------------------------------------C 
 fVel=.false.
 fDisch=.false.
 do iNS=1,NWTSol

   if (iWTSolFr(iNS)==0) then !--Time-before-extrapolation
     if (IWTSolExtrMode0(iNS)==1) then !--Leave unchanged
       CYCLE
     elseif (IWTSolExtrMode0(iNS)==2) then !--Constant
       r0=0.d0
       r1=1.d0
     endif
   elseif (iWTSolFr(iNS)>=WTSolNFr(iNS)) then !--Time-after-extrapolation
     if (fWTSolExtr1(iNS)) CYCLE !--Don't make the same extrapolation twice
     if (IWTSolExtrMode1(iNS)==1) then !--Zero
       r0=0.d0
       r1=0.d0
     elseif (IWTSolExtrMode1(iNS)==2) then !--Constant
       r0=1.d0
       r1=0.d0
     endif
   else !--Time-interpolation inside
     Time0=WTSolTimes(iNS)%V(iWTSolFr(iNS))
     Time1=WTSolTimes(iNS)%V(iWTSolFr(iNS)+1)
     r0=(Time1-TIME)/(Time1-Time0)
     r1=1.d0-r0
   endif

   call PointToVarArray(IWTSolVR(iNS))

   if (IWTSolVR(iNS)==4) then !--Water Surface 
!     YH(1:mPoints) = PA1(1:mPoints) +TSF(1:mPoints)  
!     PA1=>YH
   endif

   !----------------------------------------------------------------------!
   !     Define specific mode for Big Arrays: 1D or 2D.
   !----------------------------------------------------------------------!
   select case (IWTSolVR(iNS))
!   case (40) !--Use 2D array instead of default 1D, cause 2 derivatives have been calculated
!     PA1=>PaDX
!     PA2=>PaDY
   end select

   !----------------------------------------------------------------------!
   !     Interpolate Big Arrays to current time.
   !----------------------------------------------------------------------!
   if (WTSolVDim(iNS)==2) then !--2D arrays
     do ic=1,mPoints
!       if (IXP(ic)<=0) CYCLE
       PA1(ic)=r0*WTSolV0(1,iNS)%V(ic)+r1*WTSolV1(1,iNS)%V(ic)
       PA2(ic)=r0*WTSolV0(2,iNS)%V(ic)+r1*WTSolV1(2,iNS)%V(ic)
     enddo
   else !--1D arrays
     do ic=1,mPoints
!       if (IXP(ic)<=0) CYCLE
       PA1(ic)=r0*WTSolV0(1,iNS)%V(ic)+r1*WTSolV1(1,iNS)%V(ic)
     enddo
   endif

   !----------------------------------------------------------------------!
   !     Do variable specific operations.
   !     Raise flags for more global variable specific operations that can use 
   !     several Solutions (2D vectors for example) to do later.
   !----------------------------------------------------------------------!
   select case (IWTSolVR(iNS))
   case (4) !--Water Surface Elevation
!     if (iSHM==2) call WSToHAll(WS,HR) !--godunov
   case (5) !--Water Depth
!     if (iSHM==2) call HToWSAll(HR,WS) !--godunov
   case (25) !--Flow Velocity
     fVel=.true.
   case (23:24) !--Flow Discharge
     fDisch=.true.
   end select

   !----------------------------------------------------------------------!
   !     Mark the solution as been time-after-extrapolated (for next times).
   !----------------------------------------------------------------------!
   if (iWTSolFr(iNS)>=WTSolNFr(iNS)) fWTSolExtr1(iNS)=.true.
 enddo

!----------------------------------------------------------------------!
!     If Velocity was mentioned, calc discharges.
!----------------------------------------------------------------------!
 if (fVel) then
   do ic=1,mPoints
!     if (IXP(ic)<=0) CYCLE
     QH(ic)=UH(ic)*HR(ic)
   enddo
 endif
!----------------------------------------------------------------------!
!     If Discharge was mentioned, calc velocities.
!----------------------------------------------------------------------!
 if (fDisch) then
   do ic=1,mPoints
     if (HR(ic)>=SMALL) then
       UH(ic)=QH(ic)/HR(ic)
     endif
   enddo
 endif
!----------------------------------------------------------------------C 
    RETURN
    END  
    
!======================================================================C  
    
SUBROUTINE LoadWaterSolutionFrame(iWTSol,fConv)
!----------------------------------------------------------------------C 
!
! LoadWaterSolutionFrame: Load one iWTSol Water Solution frame from solution file,
!  reinterpolate if needed.
!
!----------------------------------------------------------------------C 
!USE GLOBAL 
USE WTSOL 
USE WORK
IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
!----------------------------------------------------------------------C 
INTEGER*4 iWTSol
LOGICAL*4 fConv
!----------------------------------------------------------------------C 
REAL*8, DIMENSION(:), POINTER :: tV1,tV2 !--Shared pointer arrays
!----------------------------------------------------------------------C 
!----------------------------------------------------------------------!
!     Fill Before-Values Array.
!----------------------------------------------------------------------!
 if (fConv) then !--Skip if not converting
   WTSolV0(1,iWTSol)%V=WTSolV1(1,iWTSol)%V
   if (WTSolVDim(iWTSol)==2) WTSolV0(2,iWTSol)%V=WTSolV1(2,iWTSol)%V
 else !--Skip frame
   call SkipSolutionFrameProc(IWTSolIFU(iWTSol),IWTSolFT(iWTSol),  WTSolNDim(iWTSol),WTSolNN(iWTSol), WTSoliPar(:,iWTSol))
   RETURN
 endif
 if (iWTSolFr(iWTSol)>=WTSolNFr(iWTSol)) RETURN !--Use constant extrapolation after solution times

 iGrid=IWTSolGridInd(iWTSol)
!----------------------------------------------------------------------!
!     Set work arrays for reading solution frame.
!----------------------------------------------------------------------!
 if (iGrid>0) then !--Non model grid
   tV1=>RWORK
   tV2=>RWORK(WTSolNN(iWTSol)+1:2*WTSolNN(iWTSol))
 else !--Model grid
   tV1=>WTSolV1(1,iWTSol)%V
   tV2=>WTSolV1(2,iWTSol)%V
 endif
!----------------------------------------------------------------------!
!     Load new solution frame.
!----------------------------------------------------------------------!
 call ReadSolutionFrameProc(IWTSolIFU(iWTSol),IWTSolFT(iWTSol), WTSolNDim(iWTSol),WTSolNN(iWTSol), frTime,tV1,tV2, &
 WTSoliPar(:,iWTSol))

 if (.not. fConv) RETURN !--Skip if not converting
!----------------------------------------------------------------------!
!     Reinterpolate for non-model grid,
!     Use MD=2 cause we can reinterp to cells or nodes.
!----------------------------------------------------------------------!
 if (iGrid>0) then !--Non model grid
!   call Reinterpolate(3,IWTSolGridRI(iGrid), WTSolNDim(iWTSol), &
!    WTSolGridNTR(iGrid), WTSolGridTri(iGrid)%V, &
!    WTSolGridNN(iGrid), WTSolGridNX(iGrid), WTSolGridNY(iGrid), &
!    WTSolGridX(iGrid)%V,WTSolGridY(iGrid)%V, tV1,tV2, &
!    NC, CX,CY, WTSolV1(1,iWTSol)%V,WTSolV1(2,iWTSol)%V, &
!    WTSolGridIWK(iGrid)%V,WTSolGridRWK(iGrid)%V, -1,0.d0)
 endif
!----------------------------------------------------------------------!
!     Adjust by var unit.
!----------------------------------------------------------------------!
 CU=WTSolCU(iWTSol)
 WTSolV1(1,iWTSol)%V=CU*WTSolV1(1,iWTSol)%V
 if (WTSolVDim(iWTSol)==2) WTSolV1(2,iWTSol)%V=CU*WTSolV1(2,iWTSol)%V
!----------------------------------------------------------------------!
!     Mark this Water Solution as been updated.
!----------------------------------------------------------------------!
 IWTSolUpd(iWTSol)=.true.
!----------------------------------------------------------------------C 
    RETURN
    END 
    
!======================================================================C   
    
SUBROUTINE RDWaterSolutions
!----------------------------------------------------------------------C 
!
! RDWaterSolutions: ReaD Water Solutions input group
!                 and get solution files information.
!
!----------------------------------------------------------------------C 
!USE GLOBAL 
USE DRLINK
USE FILINX
USE PVSPNL
USE SOLVAR
USE WTSOL  
USE WORK
use intersub
IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
!----------------------------------------------------------------------C 
CHARACTER*500 str,str1
CHARACTER*80 Units
CHARACTER sID*150
INTEGER*4 sNFr,sNC,sNN,sNDim,sIObT,iPar(10)
REAL*8 Times(1000000)

INTEGER*4 NTRg,NNg,NXg,NYg
REAL*8, POINTER :: xg(:),yg(:),hg(:)=>NULL()
INTEGER*4, POINTER :: Trig(:,:)
INTEGER*4 NIWK,NWK
INTEGER*4, POINTER :: iwk(:)
REAL*8, POINTER :: wk(:)
!----------------------------------------------------------------------C 
!----------------------------------------------------------------------!
!     Write header to output file.
!----------------------------------------------------------------------!
 write(IWR,'(/A)') ' Water Solutions'
 write(IWR,'(A)')  ' ---------------'

!----------------------------------------------------------------------!
!     Set initial values.
!----------------------------------------------------------------------!
 NWTSolGrids=0
 maxNN=0
!----------------------------------------------------------------------!
!     Found 'Water Solutions' data.
!----------------------------------------------------------------------!
 read(IRD,*) NWTSol
 if (NWTSol==0) RETURN

!----------------------------------------------------------------------!
!     ALLOCATE ARRAYS.
!----------------------------------------------------------------------!
 ALLOCATE(WTSolFN(NWTSol)); WTSolFN='';
 ALLOCATE(IWTSolIFU(NWTSol), IWTSolVR(NWTSol), IWTSolFT(NWTSol)); IWTSolIFU=0; IWTSolVR=0; IWTSolFT=0;
 ALLOCATE(WTSolCU(NWTSol)); WTSolCU=0.d0;
 ALLOCATE(IWTSolUpd(NWTSol)); IWTSolUpd=.false.;
 ALLOCATE(IWTSolExtrMode0(NWTSol),IWTSolExtrMode1(NWTSol)); IWTSolExtrMode0=1; IWTSolExtrMode1=1;
 ALLOCATE(fWTSolExtr1(NWTSol)); fWTSolExtr1=.false.;
 ALLOCATE(WTSolNDim(NWTSol),WTSolNN(NWTSol)); WTSolNDim=0; WTSolNN=0;
 ALLOCATE(iWTSolFr(NWTSol), WTSolNFr(NWTSol)); iWTSolFr=-1; WTSolNFr=0;
 ALLOCATE(WTSoliPar(10,NWTSol)); WTSoliPar=0;
 ALLOCATE(WTSolTimes(NWTSol))
 ALLOCATE(WTSolVDim(NWTSol)); WTSolVDim=0;
 ALLOCATE(WTSolV0(2,NWTSol), WTSolV1(2,NWTSol))

 !-----Max sizes-----!
 ALLOCATE(IWTSolGridInd(NWTSol)); IWTSolGridInd=0;
 ALLOCATE(WTSolGridFN(NWTSol)); WTSolGridFN='';
 ALLOCATE(IWTSolGridRI(NWTSol), IWTSolGridFmt(NWTSol)); IWTSolGridRI=0; IWTSolGridFmt=0;

!----------------------------------------------------------------------!
!     READ INPUT LINES.
!----------------------------------------------------------------------!
 do iLNS = 1,NWTSol
   if (iLNS/=1) write(IWR,*)
   WTSolFN(iLNS)=''
   Units=''
   !----------------------------------------------------------------------!
   !     Read record into character buffer and convert to lower case.
   !----------------------------------------------------------------------!
   read(IRD,'(A)') str
   istart=1
   !----------------------------------------------------------------------!
   !     Read water solution variable name.
   !----------------------------------------------------------------------!
   ires = RDStr(',',istart,str,str1)

   call FindVarName(str1,iv)
   if (iv>0) then
     write(IWR,'(1X,A)') trim(PNLU(iv))
     IWTSolVR(iLNS)=iv
   else
     call Msg(0,IWR,'INPUT ERROR! - Incorrect variable name: '//str1)
     STOP
   endif

!----------------------------------------------------------------------!
!     Define default mode for Big Arrays: 1D or 2D.
!----------------------------------------------------------------------!
   WTSolVDim(iLNS)=PNDim(IWTSolVR(iLNS))
!----------------------------------------------------------------------!
!     Make some variable specific actions and skip unsupported variables.
!----------------------------------------------------------------------!
   select case (iv)
   case (4:6) !--Turn off HD Solver if HD water solution is used
     ISOLVE(1)=0
   case (33:35) !--Allocate arrays if wind water solution is used
     fWind=1.0 !! .true. !Logical to real
   case default
     CYCLE
   end select

!----------------------------------------------------------------------!
!     Read water solution file type.
!----------------------------------------------------------------------!
   ires = RDStr(',',istart,str,str1)
   call SolutionFileType(str1,IWTSolFT(iLNS))
   select case (IWTSolFT(iLNS))
   case (2)
     write(IWR,'(A)')' Solution File Type: Bouss2DBinary'
   case (3)
     write(IWR,'(A)')' Solution File Type: SMSGenericBinary'
   case (4)
     write(IWR,'(A)')' Solution File Type: IMMSPBinary'
   case (5)
     write(IWR,'(A)')' Solution File Type: IMMSPAscii'
   end select

!----------------------------------------------------------------------!
!     Read water solution file name.
!----------------------------------------------------------------------!
   ires = RDStr(',',istart,str,WTSolFN(iLNS),1)
   call Msg(-1,IWR,' Solution File Name: '//trim(WTSolFN(iLNS)))
!----------------------------------------------------------------------!
!     Read water solution variable units.
!----------------------------------------------------------------------!
   ires = RDStr(',',istart,str,Units)
   WTSolCU(iLNS)=1.d0
   call ConvToSI(WTSolCU(iLNS),Units)
!----------------------------------------------------------------------!
!     Read water solution time units and base value.
!----------------------------------------------------------------------!
   ires = RDR8(',',istart,str,WTSolBaseTime)
   ires = RDStr(',',istart,str,Units)
   WTSolTScale=1.d0
   call ConvToSI(WTSolTScale,Units)
   call Msg(-1,IWR,'  Solution File Time Units: '//trim(Units))
   write(IWR,'(A,1PG14.7,A)') '  Model Zero-Time in Solution File Time-System: ',WTSolBaseTime,','//trim(Units)

!----------------------------------------------------------------------!
!     Get information from new solution file.
!----------------------------------------------------------------------!
   write(*,'(A$)') 'Calculating number of frames in file:  '//trim(WTSolFN(iLNS))//' .....'
   IWTSolIFU(iLNS)=iFUWTSol+iLNS-1
   call OpenSolutionFileProc(IWTSolIFU(iLNS), WTSolFN(iLNS), IWTSolFT(iLNS))
   call ReadSolutionHeaderTimesProc(IWTSolIFU(iLNS), IWTSolFT(iLNS),sID,sIObT,sNC,sNFr,sNDim,sNN,Times,iPar)
   write(*,'(A,I0,A)') 'Done (',sNFr,' frames)'

   WTSolNFr(iLNS)=sNFr
   WTSolNDim(iLNS)=sNDim
   WTSolNN(iLNS)=sNN
   WTSoliPar(1:10,iLNS)=iPar(1:10)
!----------------------------------------------------------------------!
!     If not coincident dimensions invoke error message.
!----------------------------------------------------------------------!
   if (sNDim/=PNDim(IWTSolVR(iLNS))) then
     call Msg(0,IWR,'INPUT ERROR! - Number of dimensions of water solution file and corresponding variable aren''t equal.')
     STOP
   endif
!----------------------------------------------------------------------!
!     If no frames invoke error message.
!----------------------------------------------------------------------!
   if (sNFr==0) then
     call Msg(0,IWR,'INPUT ERROR! - No frames in Water Solution file.')
     STOP
   endif

!----------------------------------------------------------------------!
!     Output Time info.
!----------------------------------------------------------------------!
   write(IWR,'(A,I0)') '  Number of frames:  ', sNFr
   write(IWR,'(A,1PG14.7,A)') '  Start Time:', Times(1),','//trim(Units)
   write(IWR,'(A,1PG14.7,A)') '  End Time:', Times(sNFr),','//trim(Units)
   if (sNFr>1) then
     ADT=(Times(sNFr)-Times(1))/(sNFr-1)
     write(IWR,'(A,1PG14.7,A)') '  Average DT:', ADT,','//trim(Units)
   endif
!----------------------------------------------------------------------!
!     Convert water solution times to seconds with shift.
!----------------------------------------------------------------------!
   do i=1,sNFr
     Times(i)=(Times(i)-WTSolBaseTime)*WTSolTScale
   enddo

!----------------------------------------------------------------------!
!     Allocate Frame Arrays.
!----------------------------------------------------------------------!
   if (WTSolVDim(iLNS)==2) then !--Reserve 2D arrays
     ALLOCATE(WTSolV0(1,iLNS)%V(mPoints), WTSolV1(1,iLNS)%V(mPoints));  WTSolV0(1,iLNS)%V=0.d0; WTSolV1(1,iLNS)%V=0.d0;
     ALLOCATE(WTSolV0(2,iLNS)%V(mPoints), WTSolV1(2,iLNS)%V(mPoints));  WTSolV0(2,iLNS)%V=0.d0; WTSolV1(2,iLNS)%V=0.d0;
   else !--Reserve 1D arrays
     ALLOCATE(WTSolV0(1,iLNS)%V(mPoints), WTSolV1(1,iLNS)%V(mPoints));  WTSolV0(1,iLNS)%V=0.d0; WTSolV1(1,iLNS)%V=0.d0;
   endif

!----------------------------------------------------------------------!
!     Allocate and fill Time Array.
!----------------------------------------------------------------------!
   ALLOCATE(WTSolTimes(iLNS)%V(sNFr)); WTSolTimes(iLNS)%V=0.d0;
   WTSolTimes(iLNS)%V=Times

!----------------------------------------------------------------------!
!     Read time-before-extrapolation mode.
!----------------------------------------------------------------------!
   ires = RDStr(',',istart,str,str1)
   if (str1(1:5)=='const') then !--Constant
     IWTSolExtrMode0(iLNS)=2
     write(IWR,'(A)') ' Solution time-before-extrapolation mode: constant'
   else !--Leave unchanged
     IWTSolExtrMode0(iLNS)=1
     write(IWR,'(A)') ' Solution time-before-extrapolation mode: leave unchanged'
   endif

!----------------------------------------------------------------------!
!     Read time-after-extrapolation mode.
!----------------------------------------------------------------------!
   ires = RDStr(',',istart,str,str1)
   if (str1(1:5)=='const') then !--Constant
     IWTSolExtrMode1(iLNS)=2
     write(IWR,'(A)') ' Solution time-after-extrapolation mode: constant'
   else !--Zero
     IWTSolExtrMode1(iLNS)=1
     write(IWR,'(A)') ' Solution time-after-extrapolation mode: zero'
   endif

!----------------------------------------------------------------------!
!     Read water solution grid info.
!----------------------------------------------------------------------!
   ires = RDStr(',',istart,str,str1)
   if (str1(1:5)=='irreg') then !--Non model grid
     IGridRI=1
   elseif (str1(1:4)=='rect') then !--Non model grid
     IGridRI=2
   else !--Model grid
     IGridRI=0
     IWTSolGridInd(iLNS)=0
     write(IWR,'(A)') ' Solution Grid: model grid used'
   endif


   !----------------------------------------------------------------------!
   !     If not model grid then read grid filename and find if already mentioned.
   !----------------------------------------------------------------------!
   if (IGridRI>=1) then
     maxNN=max(maxNN,sNDim*sNN)
     ires = RDStr(',',istart,str,str1,1)

     !-----Look for a name-----!
     do k=1,NWTSolGrids
       if (str1==WTSolGridFN(k)) then
         IWTSolGridInd(iLNS)=k
         EXIT
       endif
     enddo
     call Msg(-1,IWR,' Solution Grid FileName: '//trim(str1))
     !-----Add new grid-----!
     if (IWTSolGridInd(iLNS)==0) then
       NWTSolGrids=NWTSolGrids+1
       WTSolGridFN(NWTSolGrids)=str1
       IWTSolGridRI(NWTSolGrids)=IGridRI
       IWTSolGridInd(iLNS)=NWTSolGrids
     else
       if (IGridRI/=IWTSolGridRI(IWTSolGridInd(iLNS))) then
         call Msg(0,IWR,'INPUT ERROR! - Grid type mismatches the previously mentioned type.')
         STOP
       endif
     endif
   endif

 enddo

!----------------------------------------------------------------------!
!     Reserve (Ndim*maxNN) elements in RWORK array for Solution
!     non-model grid frames.
!----------------------------------------------------------------------!
 NRWORK=max(NRWORK,2*maxNN)



 if (NWTSolGrids==0) RETURN
!----------------------------------------------------------------------!
!     ALLOCATING ARRAYS.
!----------------------------------------------------------------------!
 ALLOCATE(WTSolGridNN(NWTSolGrids), WTSolGridNX(NWTSolGrids), WTSolGridNY(NWTSolGrids), WTSolGridNTR(NWTSolGrids));
 WTSolGridNN=0; WTSolGridNX=0; WTSolGridNY=0; WTSolGridNTR=0;
 ALLOCATE(WTSolGridX(NWTSolGrids),WTSolGridY(NWTSolGrids),WTSolGridTri(NWTSolGrids),WTSolGridIWK(NWTSolGrids), &
 WTSolGridRWK(NWTSolGrids))
!----------------------------------------------------------------------!
!     LOADING GRIDS.
!----------------------------------------------------------------------!
 write(IWR,'(/A/A)') ' Non model Water Solution Grids:',' ------------------------------'

!----------------------------------------------------------------------!
!     Loop through grids.
!----------------------------------------------------------------------!
 do i=1,NWTSolGrids

   if (i/=1) write(IWR,*)
   call Msg(-1,IWR,'  Grid : '//WTSolGridFN(i))
   call GetExtension(WTSolGridFN(i),str1)
   if (str1(1:2)=='14') then !--ADCIRC 14 file
     IWTSolGridFmt(i)=1
   elseif (str1(1:3)=='xyz') then !--XYZ file
     IWTSolGridFmt(i)=2
   elseif (str1(1:2)=='xy') then !--XY file
     IWTSolGridFmt(i)=3
   else
     IWTSolGridFmt(i)=3
   endif
   !----------------------------------------------------------------------!
   !     Read grid and allocate grid arrays.
   !----------------------------------------------------------------------!
!!   call ReadAllocGrid(WTSolGridFN(i),IWTSolGridFmt(i),IWTSolGridRI(i), NTRg,Trig,NNg,NXg,NYg,xg,yg,hg)
   !----------------------------------------------------------------------!
   !     Prepare for reinterpolation and allocate work arrays enough for
   !     reinterpolation to nodes and cells.
   !----------------------------------------------------------------------!
!!   call PrepareReinterpolator(IWTSolGridFmt(i),IWTSolGridRI(i), NTRg,Trig,NNg,NXg,NYg,xg,yg,hg, NC, CX,CY,VV1, NIWK,NWK, iwk,wk, -1,0.d0)
   !----------------------------------------------------------------------!
   !     Put evrth to big arrays, deallocate not needed ones.
   !----------------------------------------------------------------------!
   if (ASSOCIATED(hg)) DEALLOCATE(hg) !--Not needed

   if (IWTSolGridFmt(i)==1) then
     write(IWR,'(A)') '   Solution Grid Format: ADCIRC 14 (already triangulated, no extrapolation possible)'
   elseif (IWTSolGridFmt(i)==2) then
     write(IWR,'(A)') '   Solution Grid Format: XYZ (not triangulated yet, extrapolation is possible)'
   elseif (IWTSolGridFmt(i)==3) then
     write(IWR,'(A)') &
     '   Solution Grid Format: XY (not triangulated yet, extrapolation is possible, no Z values (zero assumed if needed))'
   endif

   if (IWTSolGridRI(i)==1) then
     write(IWR,'(A)') '   Solution Grid Type  : irregular'
     write(IWR,'(A,I0)') '  Number of Nodes:  ',NNg
   elseif (IWTSolGridRI(i)==2) then
     write(IWR,'(A)') '   Solution Grid Type  : rectilinear'
     write(IWR,'(A,I0)') '  Number of Nodes in X-dir.:  ',NXg
     write(IWR,'(A,I0)') '  Number of Nodes in Y-dir.:  ',NYg
   endif
   write(IWR,'(A,I0)') '  Number of Cells:  ',NTRg

   WTSolGridNN(i)=NNg
   WTSolGridNX(i)=NXg
   WTSolGridNY(i)=NYg
   WTSolGridNTR(i)=NTRg

   WTSolGridX(i)%V=>xg !--Share array
   WTSolGridY(i)%V=>yg !--Share array
   WTSolGridTri(i)%V=>Trig !--Share array, (Trig can contain more than NTRg elements!) 
   WTSolGridIWK(i)%V=>iwk !--Share array
   WTSolGridRWK(i)%V=>wk !--Share array

 enddo

!----------------------------------------------------------------------!
!     Check for matching of dimensions.
!----------------------------------------------------------------------!
 do iWTSol=1,NWTSol
   iGrid=IWTSolGridInd(iWTSol)
   if (iGrid>0) then
     if (WTSolNN(iWTSol)/=WTSolGridNN(iGrid)) then
       call Msg(0,IWR,'INPUT ERROR! - Dimensions of Water Solution file and corresponding grid aren''t equal.')

       call Msg(0,IWR,'Solution File Name: '//trim(WTSolFN(iWTSol)),1)
       write(IWR,'(A,I0)') ', number of nodes =  ',WTSolNN(iWTSol)
       write(*,'(A,I0)') ', number of nodes =  ',WTSolNN(iWTSol)

       call Msg(0,IWR,'Grid File Name: '//trim(WTSolGridFN(iGrid)),1)
       write(IWR,'(A,I0)') ', number of nodes =  ',WTSolGridNN(iGrid)
       write(*,'(A,I0)') ', number of nodes =  ',WTSolGridNN(iGrid)
       STOP
     endif
   endif
 enddo
!----------------------------------------------------------------------C 
    RETURN
    END 
    
!======================================================================C   
    
SUBROUTINE AdjustSpecificWaterSolutions
!------------------------------------------------------------------------------C
!
! AdjustSpecificWaterSolutions: Make variable specific preparations.
!
!------------------------------------------------------------------------------C
!USE GLOBALALL   
USE DRLINK
USE WTSOL 
USE WORK
IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
!------------------------------------------------------------------------------C
!----------------------------------------------------------------------!
!     Find variables that need specific preparations.
!----------------------------------------------------------------------!
 do iWTSol=1,NWTSol
   if (.not.IWTSolUpd(iWTSol)) CYCLE !--Skip if solution hasn't been updated.
   select case (IWTSolVR(iWTSol))
   !----------------------------------------------------------------------!
   !     Convert Air Pressure to Air Pressure Derivatives.
   !     Use zero pressure as NoData flag to calc derivatives only in
   !     that part of model grid covered by solution grid.
   !----------------------------------------------------------------------!
   case (40)
     VV1(1:mPoints)=WTSolV1(1,isPa)%V(1:mPoints)
!     call CalcDerivativeXNoBC(VV1, WTSolV1(1,iWTSol)%V)
!     call CalcDerivativeYNoBC(VV1, WTSolV1(2,iWTSol)%V)
   end select
 enddo
!------------------------------------------------------------------------------C
    RETURN
    END  
    
!======================================================================C   
    
SUBROUTINE SkipSolutionFrameProc(IFU,IFType, NDim,NN, iPar)
!------------------------------------------------------------------------------C
IMPLICIT NONE
!------------------------------------------------------------------------------C
INTEGER*4 IFU,IFType,NDim,NN, iPar(10)
!------------------------------------------------------------------------------C
 select case (IFType)
 case (3)
   call SkipSolutionFrame_SMSGenericBinary(IFU,NDim,NN,iPar)
 case (4)
   call SkipSolutionFrame_IMMSPBinary(IFU,NDim,NN)
 case (5)
   call SkipSolutionFrame_IMMSPAscii(IFU)
 end select
!------------------------------------------------------------------------------C
    RETURN
    END  
    
!======================================================================C   
    
SUBROUTINE ReadSolutionFrameProc(IFU, IFType, NDim,NN, frTime,V1,V2, iPar)
!------------------------------------------------------------------------------C
IMPLICIT NONE
!------------------------------------------------------------------------------C
INTEGER*4 IFU,IFType,NDim,NN, iPar(10)
REAL*8 frTime,V1(NN),V2(NN)
!------------------------------------------------------------------------------C
 select case (IFType)
 case (3)
   call ReadSolutionFrame_SMSGenericBinary(IFU,NDim,NN,iPar, frTime,V1,V2)
 case (4)
   call ReadSolutionFrame_IMMSPBinary(IFU,NDim,NN, frTime,V1,V2)
 case (5)
   call ReadSolutionFrame_IMMSPAscii(IFU,NDim,NN, frTime,V1,V2)
 end select
!------------------------------------------------------------------------------C
    RETURN
    END  
    
!======================================================================C   
    
SUBROUTINE FindTimeInterval(Time,NT,Times,it0,it)
!------------------------------------------------------------------------------C
! FindTimeInterval: Find it>=it0 such that Times(it)<Time<=Times(it+1),
!                            if Time<=Times(1) => it=0,
!                            if Time>Times(NT) => it=NT.    
! my change:                 if Time < Times(1) => it=0,  
!------------------------------------------------------------------------------C
INTEGER*4 NT,it
REAL*8 Time,Times(NT)
!------------------------------------------------------------------------------C
 if( it0 == 0 .and. Time == Times(1) ) then 
   it = 1; return;  
 endif  
 do it=max(it0,0),NT-1   
   if (Time<=Times(it+1)) RETURN
 enddo
 it=NT
!------------------------------------------------------------------------------C
    RETURN
    END SUBROUTINE  
    
!======================================================================C  
!    
SUBROUTINE SolutionFileType(FTypeStr,IFType)
!------------------------------------------------------------------------------C
IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
!------------------------------------------------------------------------------C
CHARACTER*(*) FTypeStr
INTEGER*4 IFType
!------------------------------------------------------------------------------C
 IFType=3
 if (FTypeStr(1:6)=='bouss2d') then
!   IFType=2
 elseif (FTypeStr(1:6)=='smsgen') then
   IFType=3
 elseif (FTypeStr(1:8)=='immspbin') then
   IFType=4
 elseif (FTypeStr(1:8)=='immspasc') then
   IFType=5
 endif
!------------------------------------------------------------------------------C
    RETURN
    END 
    
!======================================================================C 
!  
SUBROUTINE GetExtension(FName,FExt)
!------------------------------------------------------------------------------C
CHARACTER*(*) FName,FExt
INTEGER*4 ln,i
!------------------------------------------------------------------------------C
 ln=len_trim(FName)
 i=index(FName,'.',BACK = .true.)
 FExt=''
 if ((i>0).and.(i<ln)) write(FExt,'(A)') FName(i+1:ln)
!------------------------------------------------------------------------------C
RETURN
END SUBROUTINE
    
!======================================================================C 
!  
SUBROUTINE OpenSolutionFileProc(IFU, FName, IFType)
!------------------------------------------------------------------------------C
IMPLICIT NONE
!------------------------------------------------------------------------------C
INTEGER*4 IFU,IFType
CHARACTER*(*) FName
!------------------------------------------------------------------------------C
 select case (IFType)
 case (-1,5) !--text files
   OPEN(IFU,FILE=FName, ACCESS = 'STREAM') !OPENMARK
 case (2,3,4) !--binary files
   OPEN(IFU,FILE=FName, ACCESS = 'STREAM', FORM='UNFORMATTED') !OPENMARK
 end select
!------------------------------------------------------------------------------C
RETURN
END
    
!======================================================================C 
!  
SUBROUTINE ReadSolutionHeaderProc(IFU, IFType, ID,IObT,NC,NFr,NDim,NN,StartT,EndT,iPar)
!------------------------------------------------------------------------------C
IMPLICIT NONE
!------------------------------------------------------------------------------C
INTEGER*4 IFU,IFType
CHARACTER*(*) ID
INTEGER*4 IObT,NC,NFr,NDim,NN,iPar(10)
REAL*8 StartT,EndT
!------------------------------------------------------------------------------C
INTEGER*4 ITU
!------------------------------------------------------------------------------C
 select case (IFType)
 case (3)
   call ReadSolutionHeader_SMSGenericBinary(IFU,ID,IObT,NC,NFr,NDim,NN,ITU,StartT,EndT,iPar)
 case (4)
   call ReadSolutionHeader_IMMSPBinary(IFU, ID,NFr,NDim,NN,StartT,EndT)
 case (5)
   call ReadSolutionHeader_IMMSPAscii(IFU, ID,NFr,NDim,NN,StartT,EndT)
 end select
!------------------------------------------------------------------------------C
RETURN
END
    
!======================================================================C 
!  
SUBROUTINE ReadSolutionHeaderTimesProc(IFU, IFType, ID,IObT,NC,NFr,NDim,NN,Times,iPar)
!------------------------------------------------------------------------------C
IMPLICIT NONE
!------------------------------------------------------------------------------C
INTEGER*4 IFU,IFType
CHARACTER*(*) ID
INTEGER*4 IObT,NC,NFr,NDim,NN,iPar(10)
REAL*8 Times(*)
!------------------------------------------------------------------------------C
INTEGER*4 ITU
!------------------------------------------------------------------------------C
 select case (IFType)
 case (3)
   call ReadSolutionHeaderTimes_SMSGenericBinary(IFU, ID,IObT,NC,NFr,NDim,NN,ITU,Times,iPar)
 case (4)
   call ReadSolutionHeaderTimes_IMMSPBinary(IFU, ID,NFr,NDim,NN,Times)
 case (5)
   call ReadSolutionHeaderTimes_IMMSPAscii(IFU, ID,NFr,NDim,NN,Times)
 end select
!------------------------------------------------------------------------------C
RETURN
END
    
!======================================================================C 
!======================================================================C 
!  
SUBROUTINE InitVariables
!------------------------------------------------------------------------------C
!
!     InitVariables: Initialize variable names and attributes.
!
!------------------------------------------------------------------------------C
USE POINTS
USE SOLVAR  
USE PVSPNL
!USE GLOBAL
IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
!------------------------------------------------------------------------------C
 i=1;  PNS(i)='A';   PNS2(1,i)='';    PNS2(2,i)='';    PNF(i)=0; PNL(i)='Cross-Sectional Area of Flow'; PNU(i)='m^2'
 i=2;  PNS(i)='MN';  PNS2(1,i)='';    PNS2(2,i)='';    PNF(i)=0; PNL(i)='Manning Roughness Coefficient'; PNU(i)=''

 i=3;  PNS(i)='tw';  PNS2(1,i)='';    PNS2(2,i)='';    PNF(i)=1; PNL(i)='Thalweg Elevation'; PNU(i)='m'
 i=4;  PNS(i)='ws';  PNS2(1,i)='';    PNS2(2,i)='';    PNF(i)=1; PNL(i)='Water Surface Elevation'; PNU(i)='m'
 i=5;  PNS(i)='h';   PNS2(1,i)='';    PNS2(2,i)='';    PNF(i)=1; PNL(i)='Water Depth'; PNU(i)='m'
 i=6;  PNS(i)='q';   PNS2(1,i)='';    PNS2(2,i)='';    PNF(i)=1; PNL(i)='Flow Discharge'; PNU(i)='m^3/s'
 i=7;  PNS(i)='zd';  PNS2(1,i)='';    PNS2(2,i)='';    PNF(i)=1; PNL(i)='Depth of Bottom Top-Layer'; PNU(i)='m'

 i=8;  PNS(i)='sh';  PNS2(1,i)='';    PNS2(2,i)='';    PNF(i)=1; PNL(i)='Sediment Concentration'; PNU(i)='kg/m^3'
 i=9;  PNS(i)='seq'; PNS2(1,i)='';    PNS2(2,i)='';    PNF(i)=0; PNL(i)='Sediment Equilibrium Concentration'; PNU(i)='kg/m^3'
 i=10; PNS(i)='ch';  PNS2(1,i)='';    PNS2(2,i)='';    PNF(i)=1; PNL(i)='Soluble Species Concentration'; PNU(i)='1/m^3'
 i=11; PNS(i)='csh'; PNS2(1,i)='';    PNS2(2,i)='';    PNF(i)=1; PNL(i)='Particulate Species Concentration'; PNU(i)='1/m^3'
 i=12; PNS(i)='cb';  PNS2(1,i)='';    PNS2(2,i)='';    PNF(i)=1; PNL(i)='Bottom Species Concentration'; PNU(i)='1/m^3'

 i=13; PNS(i)='qs';  PNS2(1,i)='';    PNS2(2,i)='';    PNF(i)=0; PNL(i)='Deposition Rate'; PNU(i)='kg/m^2/s'
 i=14; PNS(i)='qb';  PNS2(1,i)='';    PNS2(2,i)='';    PNF(i)=0; PNL(i)='Erosion Rate'; PNU(i)='kg/m^2/s'

! i=18; PNS(i)='vm';  PNS2(1,i)='';    PNS2(2,i)='';    PNF(i)=0; PNL(i)='Current Magnitude'; PNU(i)='m/s'
! i=19; PNS(i)='va';  PNS2(1,i)='';    PNS2(2,i)='';    PNF(i)=0; PNL(i)='Current Direction'; PNU(i)='deg'

 NVARS=i !--Total number of variables
!----------------------------------------------------------------------!
!     Post-processing.
!----------------------------------------------------------------------!
 do i=1,NVARS
   PNDim(i)=1; if (PNS2(1,i)/='') PNDim(i)=2; !--Dimensionality
   PNLU(i)= trim(PNL(i))//', '//trim(PNU(i)) !--Long name with units
 enddo
!------------------------------------------------------------------------------C
RETURN
END
    
!======================================================================C 
!======================================================================C 
!  
SUBROUTINE ReadConvTable(FName,delim,NHdrLines,NLines,iUnMode,NV,V1,V2)
!------------------------------------------------------------------------------C
! Read Table delimited with spaces and converts to SI system with units in the file header.
! File contains header of NHdrLines lines than 1 line with one unit for all vars if iUnMode=1, or
! units for each var if iUnMode=2. Then NLines lines with NV vars in each.
!------------------------------------------------------------------------------C 
USE RW_MOD  
use intersub
IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
!------------------------------------------------------------------------------C
CHARACTER*(*) FName
CHARACTER delim
INTEGER*4 NHdrLines,NLines,iUnMode,NV
REAL*8 V1(NLines),V2(NLines)
!------------------------------------------------------------------------------C
INTEGER*4 NU
CHARACTER*500 str
CHARACTER*80 Units(2)
REAL*8 ru(2)
REAL*8 Value(2)
!------------------------------------------------------------------------------C
 ru=1.d0

 OPEN(UNIT=RW_IFUTemp, FILE=FName)

 !----------------------------------------------------------------------!
 !     Skip header lines.
 !----------------------------------------------------------------------!
 do i=1,NHdrLines
   read(RW_IFUTemp,*)
 enddo

 !----------------------------------------------------------------------!
 !     Read unit line.
 !----------------------------------------------------------------------!
 if (iUnMode==1) then
   NU=1
 else
   NU=2
 endif

 read(RW_IFUTemp,'(A)') str
! istart=1
! if (delim==' ') call ElimConsecSpaces(str)


! do j = 1, NU
!   ires = RDStr(delim,istart,str,Units(j))
!   call ConvToSI(ru(j), Units(j))
! enddo
! do j = NU+1, 2
!   ru(j)=ru(1)
! enddo

 !----------------------------------------------------------------------!
 !     Read table.
 !----------------------------------------------------------------------!
 do i=1,NLines
   read(RW_IFUTemp,'(A)') str
   istart=1
   if (delim==' ') call ElimConsecSpaces(str)
   do j = 1, NV
     ires = RDR8(delim,istart,str, Value(j))
!     Value(j)=Value(j)*ru(j)
   enddo
   V1(i)=Value(1)
   if (NV>=2) V2(i) = Value(2)
 enddo

 CLOSE(RW_IFUTemp)
!------------------------------------------------------------------------------C
RETURN
END SUBROUTINE
    
!======================================================================C 
!======================================================================C 
!  
SUBROUTINE ElimConsecSpaces(str)
!------------------------------------------------------------------------------C
!     Eliminates leading and consecutive spaces.
!------------------------------------------------------------------------------C
CHARACTER*(*) str
CHARACTER ch,chprev
INTEGER*4 i,j,ln
!------------------------------------------------------------------------------C
 ln=len_trim(str)
 j=0
 chprev=' '
 do i=1,ln
   ch=str(i:i)
   if ((ch/=' ').or.(chprev/=' ')) then
     j=j+1
     str(j:j)=ch
     chprev=ch
   endif
 enddo
 str(j+1:ln)=' '
!------------------------------------------------------------------------------C
RETURN
END SUBROUTINE
    
!======================================================================C 
!======================================================================C 
!  
    
!======================================================================C 
!======================================================================C 
!  
    
!======================================================================C 
