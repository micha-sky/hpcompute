!======================================================================C 
! 
      SUBROUTINE REFNOD( NFCNV,NITER ) 
! 
!----------------------------------------------------------------------C  
! 
!     REFNOD: REFerence NODe ouput. 
! 
!     Subroutine prints diagnostics information for a single node to 
!     the screen and/or output file. 
 
!     MSTS (Multiphase Subsurface Transport Simulator) code.  
!     M.D. White, and W.E. Nichols. Pacific Northwest National Laboratory,   
!     Richland, WA, 1993. 
!----------------------------------------------------------------------C 
      USE CONSTS
      USE DRLINK  
      USE DRSPFL
      USE DRWFLD
      USE FILINX
      USE FLDINX  
      USE LANDSF
      USE LOGCLS
      USE NUMBRS
      USE POINTS 
      USE REFERN
      USE SOLVAR
!----------------------------------------------------------------------C 
!     Implicit Double Precision. 
!----------------------------------------------------------------------C 
      IMPLICIT REAL*8 (A-H,O-Z)  
      IMPLICIT INTEGER*4 (I-N) 
!----------------------------------------------------------------------C 
!     Type Declarations. 
!----------------------------------------------------------------------C 
      CHARACTER*3 CHRF(20),CHRS(20),CHRG(50) 
      REAL*8 VARF(20),VARS(20),VARG(50) 
!----------------------------------------------------------------------C 
!     Data Statements. 
!----------------------------------------------------------------------C 
      SAVE CHRG 
      DATA CHRG/                                                       &
         ' A ','MN ','TW ',' Y ',' H ',' Q ','ZD ','SH ','SEQ','CH ',  &
         'CSH','CB ','QS ','QB ','QX ','QY ','ZD ','|U|','U-A','HWV',  &
         'DWV','TWV','TX ','TY ','???','???','???','???','???','???',  &
         '???','???','???','???','???','???','???','???','???','???',  &
         '???','???','???','???','???','???','???','???','???','???'/ 
!----------------------------------------------------------------------C 
!     Compute properties for node NREF(J).      
!----------------------------------------------------------------------C 
!      NFCNV = 1 
      DO 400 J = 1,NMREF  
        N = NREF(J) 
        L = ILK(N) 
        Ip = JPT(N)  
        do i=1,50  
          if( IRNV(i) ) then  
            select case (i)   
            case(1)    
              VARG( 1) = AW(N) 
            case(2)  
              VARG( 2) = ROU(N)  
            case(3)
              VARG( 3) = TSF(N)   
            case(4) 
              VARG( 4) = TSF(N) + HR(N)   
            case(5)
              VARG( 5) = HR(N)   
            case(6)
              VARG( 6) = QH(N)  
            case(7)  
              VARG( 7) = ZSD(N)  
            case(8)
              VARG( 8) = SH(N)  
            case(9)  
              VARG( 9) = SEQ(N)  
            case(10)
              VARG(10) = CH(N)  
            case(11)
              VARG(11) = CSH(N)  
            case(12)  
              VARG(12) = CB(N)  
            case(13)
              VARG(13) = SDOWN(N)   
            case(14)
              VARG(14) = RSUP(N) 
!        VARG(13) = QHX(NPX) 
!        VARG(14) = QHY(NPY) 
!        VARG(15) = ZSD(N)    
            end select  
          endif  
        enddo
        ICF = 0  
        DO 100 I = 1,50 
          IF( IRNV(I) .AND. ICF .LT. 9 ) THEN   
            ICF = ICF +1 
            VARF(ICF) = VARG(I) 
            CHRF(ICF) = CHRG(I) 
          ENDIF 
  100   CONTINUE 
        ICS = 0 
        DO 200 I = 1,50 
          IF( IRNV(I) .AND. ICS .LT. 9 ) THEN   
            ICS = ICS +1 
            VARS(ICS) = VARG(I) 
            CHRS(ICS) = CHRG(I) 
          ENDIF 
  200   CONTINUE 
!----------------------------------------------------------------------C 
!     Print output header. 
!----------------------------------------------------------------------C 
      IF( IHREF .EQ. 25 ) THEN 
        IHREF = 24 
        WRITE (IWR,9030) 
      ENDIF 
! 
      IF( IHREF .EQ. 24 ) THEN 
        IHREF = 0 
        WRITE (IWR,9032) 
! 
        DO 300 K = 1,NMREF 
          WRITE (IWR,9034) IREF(K),JREF(K) 
  300   CONTINUE 
! 
        WRITE (IWR,'(/)') 
        WRITE (IWR,9036) (CHRF(I),I=1,ICF) 
! 
        IF( SCREEN ) THEN 
          WRITE (*,9032) 
          DO 310 K = 1,NMREF 
            WRITE (*,9034) IREF(K),JREF(K) 
  310     CONTINUE 
          WRITE (*,'(/)') 
          WRITE (*,9040) (CHRS(I),I=1,ICS) 
        ENDIF 
! 
      ENDIF 
!----------------------------------------------------------------------C 
!     Convert pressures to absolute pressures for printing. 
!----------------------------------------------------------------------C 
      IHREF = IHREF +1 
      WRITE(IWR,9130) NSTEP,TIME/TSCL,NITER,NFCNV,DT,(VARF(I),I=1,ICF) 
      IF( SCREEN ) WRITE (*,9140) NSTEP,TIME/TSCL,NITER,NFCNV,         &
        (VARS(I),I=1,ICS) 
  400 CONTINUE 
!----------------------------------------------------------------------C 
!     Format Statements. 
!----------------------------------------------------------------------C 
 9030 FORMAT(/' REFERENCE NODE OUTPUT RECORD',                         &
            /' ----------------------------') 
 9032 FORMAT(/' REFERENCE NODE(S)',$) 
 9034 FORMAT(' (',I3,',',I3,')',$) 
!----------------------------------------------------------------------C 
 9036 FORMAT(4X,'STEP',3X,'TIME',4X,'ITER',2X,'EQUA',3X,'TIMESTEP',7X, &
            9(A2,9X)) 
 9040 FORMAT(3X,'STEP',3X,'TIME',4X,'ITER',2X,'EQUA',5X,7(A3,7X)) 
 9130 FORMAT(1X,I6,2X,1PE10.3,1X,I3,4X,I2,3X,1PE10.3,2X,9(1PE10.3,1X)) 
 9140 FORMAT(I6,1X,1PE9.2,2X,I3,4X,I2,3X,7(1PE9.2,1X)) 
!----------------------------------------------------------------------C 
!     End of REFNOD group. 
!----------------------------------------------------------------------C 
! 
      RETURN 
      END 
 
!======================================================================C 
! 
      SUBROUTINE WRPLOT(FIPLOT)  
! 
!----------------------------------------------------------------------C  
! 
!     WRPLOT: WRite PLOT files. 
! 
!----------------------------------------------------------------------C  
      USE CONSTS
      USE DRLINK
      USE DRSPFL
      USE DRWFLD
      USE FILINX
      USE FLDINX
      USE HEADNG 
      USE LANDSF
      USE LOGCLS
      USE NUMBRS
      USE POINTS
      USE SOLVAR  
      USE PVSPNL 
      USE DRNAME
!----------------------------------------------------------------------C 
!     Implicit Double Precision. 
!----------------------------------------------------------------------C  
      IMPLICIT REAL*8 (A-H,O-Z) 
      IMPLICIT INTEGER*4 (I-N) 
!----------------------------------------------------------------------C 
!     Type Declarations. 
!----------------------------------------------------------------------C 
      CHARACTER*27 FORM1 
      CHARACTER*28 FORM2 
!      CHARACTER*32 FileName 
      CHARACTER*32 ADUM, ADUM0
!      CHARACTER USER*80,TITLE*80,COMP*40 
!      CHARACTER*40 DATVAL,TIMVAL 
!      REAL*8 F(28) 
      LOGICAL FIPLOT(MNVARS)
      CHARACTER*32 str, str0
      REAL*8 F(MNVARS) 
!!*************************************************************
      LOGICAL*4 fExist
!*************************************************************
!----------------------------------------------------------------------C 
!     Data Statements. 
!----------------------------------------------------------------------C 
      !! todo  SAVE FORM1, FORM2 
      DATA FORM1 /'(I6,2X,A,2I6, (1X,1PE11.4))' / 
      DATA FORM2 /'(I6,2X,A,2I6,  (1X,1PE14.7))'/ 
!      DATA FORM1 /'(3I6, (1X,1PE18.9))' / 
!      DATA FORM2 /'(3I6,  (1X,1PE18.9))'/ 
!----------------------------------------------------------------------C 
!     Count the number of plot file variables requested. 
!----------------------------------------------------------------------C 
      NPLOT = 0 
      DO 80 N = 1,NVARS 
        IF( fIPLOT(N) ) NPLOT = NPLOT +1 
   80 CONTINUE 
      IF( NPLOT .EQ. 0 ) GOTO 200 

      goto 81
!----------------------------------------------------------------------C 
!     Open a new plot file with NSTEP as an extention. 
!----------------------------------------------------------------------C 
      ADUM = FileName(IPL) 
      IF( NSTEP .LE. 9 ) THEN 
        WRITE (ADUM(6:),'(I1)') NSTEP 
      ELSEIF( NSTEP .LE. 99 ) THEN 
        WRITE (ADUM(6:),'(I2)') NSTEP 
      ELSEIF( NSTEP .LE. 999 ) THEN 
        WRITE (ADUM(6:),'(I3)') NSTEP 
      ELSEIF( NSTEP .LE. 9999 ) THEN 
        WRITE (ADUM(6:),'(I4)') NSTEP 
      ELSEIF( NSTEP .LE. 99999 ) THEN 
        WRITE (ADUM(6:),'(I5)') NSTEP 
      ELSEIF( NSTEP .LE. 999999 ) THEN 
        WRITE (ADUM(6:),'(I6)') NSTEP 
      ELSEIF( NSTEP .LE. 9999999 ) THEN 
        WRITE (ADUM(6:),'(I7)') NSTEP 
      ELSE 
        WRITE (ADUM(5:),'(I7)') NSTEP 
      ENDIF 
      ADUM(5:5) = '.' 

  81  continue
!----------------------------------------------------------------------C 
!     Open a new plot file with TIME in seconds as an extention. 
!----------------------------------------------------------------------C 
      ADUM = FileName(IPL) 
      WRITE (ADUM(6:),'(A9)') '000000000' 
      ITIME=anint(TIME)
      IF( ITIME <= 9 ) THEN 
        WRITE (ADUM(14:),'(I1)') ITIME 
      ELSEIF( ITIME <= 99 ) THEN 
        WRITE (ADUM(13:),'(I2)') ITIME 
      ELSEIF( ITIME <= 999 ) THEN 
        WRITE (ADUM(12:),'(I3)') ITIME 
      ELSEIF( ITIME <= 9999 ) THEN 
        WRITE (ADUM(11:),'(I4)') ITIME 
      ELSEIF( ITIME <= 99999 ) THEN 
        WRITE (ADUM(10:),'(I5)') ITIME 
      ELSEIF( ITIME <= 999999 ) THEN 
        WRITE (ADUM( 9:),'(I6)') ITIME 
      ELSEIF( ITIME <= 9999999 ) THEN 
        WRITE (ADUM( 8:),'(I7)') ITIME 
      ELSEIF( ITIME <= 99999999 ) THEN 
        WRITE (ADUM( 7:),'(I8)') ITIME 
      ELSEIF( ITIME <= 999999999 ) THEN 
        WRITE (ADUM( 6:),'(I9)') ITIME 
      ELSE 
        WRITE (ADUM( 6:),'(I9)') ITIME 
      ENDIF 
      ADUM(5:5) = '.' 
      
      
      !---If file already exists, add suffix and try again--
      ADUM0 = ADUM
      nLen=LEN_TRIM(ADUM0)
	ni=0
  85  continue
      fExist=.false.
      INQUIRE (FILE=ADUM, EXIST = fExist)
	if (fExist) then
	 ni=ni+1
       ADUM=ADUM0

       if (ni <= 9) then
        write (ADUM(nLen+1:),'(A2,I1)') '_0',ni
       else
        write (ADUM(nLen+1:),'(A1,I2)') '_',ni
       endif

       goto 85
      endif
      !--------------------------------------------------

      OPEN(UNIT=IPL, FILE=ADUM,STATUS='UNKNOWN',FORM='FORMATTED') 
      CLOSE(UNIT=IPL,STATUS='DELETE') 
      OPEN(UNIT=IPL, FILE=ADUM,STATUS='NEW',FORM='FORMATTED') 
!----------------------------------------------------------------------C 
!     Write credits and timing information to the 'plot' file. 
!----------------------------------------------------------------------C 
      WRITE (IPL,9110)  
      WRITE (IPL,9160) NDATE, NTIME 
!----------------------------------------------------------------------C 
!     Print header information to the plot file. 
!----------------------------------------------------------------------C 
      WRITE (IPL,'(/)') 
      WRITE (IPL,'(A)') 'Grid Geometry- Number of Links and Points;' 
      WRITE (IPL,'(2(1X,I5))') NLINK, NPOINT 
      WRITE (IPL,'(/)') 
      WRITE (IPL,'(A)') 'Field and Flux Variable Results @' 
      WRITE (IPL,'(I11,A)') NSTEP,' steps' 
      rMTIME = TIME/60. 
      HTIME = TIME/3600. 
      DTIME = HTIME/24. 
      WTIME = DTIME/7. 
      YTIME = DTIME/365.25 
      WRITE (IPL,'(6(1PE11.4,A,/))') TIME,' s',rMTIME,' min',HTIME,' hr'&
        ,DTIME,' day',WTIME,' wk',YTIME,' yr' 
!----------------------------------------------------------------------C 
!     Print field variable order information to the plot file. 
!----------------------------------------------------------------------C 
      WRITE (IPL,'(A)') 'Field Variables @ Node Centers: Plot Order' 
      WRITE (IPL,'(A)') 'Point Number' 
      WRITE (IPL,'(A)') 'Link Name' 
      WRITE (IPL,'(A)') 'Link Number' 
      WRITE (IPL,'(A)') 'Point Number on the Link' 
      WRITE (IPL,'(A)') 'X-Direction Node Position, m' 
      IF( fIPLOT(1) ) WRITE (IPL,'(A)')                                 &
                    'Cross-Sectional Area of Flow, m^2' 
      IF( fIPLOT(2) ) WRITE (IPL,'(A)') 'Manning Roughness Coefficient' 
      IF( fIPLOT(3) ) WRITE (IPL,'(A)') 'Thalweg Elevation, m' 
      IF( fIPLOT(4) ) WRITE (IPL,'(A)') 'Water Surface Elevation, m'
      IF( fIPLOT(5) ) WRITE (IPL,'(A)') 'Flow Depth, m' 
      IF( fIPLOT(6) ) WRITE (IPL,'(A)') 'Flow Discharge, m^3/s' 
      IF( fIPLOT(7) ) WRITE (IPL,'(A)') 'Depth of Bottom Top-Layer, m' 
      IF( fIPLOT(8) ) WRITE (IPL,'(A)') 'Sediment Concentration, kg/m^3' 
      IF( fIPLOT(9) ) WRITE (IPL,'(A)') 'Sediment Equilibrium Concentration, kg/m^3' 
      IF( fIPLOT(10) ) WRITE (IPL,'(A)')'Soluble Species Concentration, 1/m^3' 
      IF( fIPLOT(11) ) WRITE (IPL,'(A)')'Particulate Species Concentration, 1/m^3' 
      IF( fIPLOT(12) ) WRITE (IPL,'(A)')'Bottom Species Concentration, 1/m^3' 
      IF( fIPLOT(13) ) WRITE (IPL,'(A)') 'Deposition Rate, kg/m^2 s' 
      IF( fIPLOT(14) ) WRITE (IPL,'(A)') 'Erosion Rate, kg/m^2 s' 

!      IF( IPLOT(13) ) WRITE (IPL,'(A)')  
!     &    'X-Direction Water Discharge, m^3/s' 
!      IF( IPLOT(14) ) WRITE (IPL,'(A)')  
!     &    'Y-Direction Water Discharge, m^3/s' 

!----------------------------------------------------------------------C 
!     Print field variables by node number to the plot file. 
!----------------------------------------------------------------------C 
      DO 110 L = 1,NLINK 
        DO 100 I = 1,IL(L) 
          N = ND(L,I)  
!
          ICNT = 0 
          IF( IPLOT(1) ) THEN 
            ICNT = ICNT +1 
            F(ICNT) = AW(N) 
          ENDIF 
          IF( IPLOT(2) ) THEN 
            ICNT = ICNT +1 
            F(ICNT) = ROU(N) 
          ENDIF 
          IF( IPLOT(3) ) THEN 
            ICNT = ICNT +1 
            F(ICNT) = TSF(N) 
          ENDIF 
          IF( IPLOT(4) ) THEN 
            ICNT = ICNT +1 
            F(ICNT) = TSF(N) + HR(N) 
          ENDIF 
          IF( IPLOT(5) ) THEN 
            ICNT = ICNT +1 
            F(ICNT) = HR(N) 
          ENDIF 
          IF( IPLOT(6) ) THEN 
            ICNT = ICNT +1 
            F(ICNT) = QH(N) 
          ENDIF 
          IF( IPLOT(7) ) THEN 
            ICNT = ICNT +1 
            F(ICNT) = ZSD(N) 
          ENDIF 
          IF( IPLOT(8) ) THEN 
            ICNT = ICNT +1 
            F(ICNT) = SH(N) 
          ENDIF 
          IF( IPLOT(9) ) THEN 
            ICNT = ICNT +1 
            F(ICNT) = SEQ(N) 
          ENDIF 
          IF( IPLOT(10) ) THEN 
            ICNT = ICNT +1 
            F(ICNT) = CH(N) 
          ENDIF 
          IF( IPLOT(11) ) THEN 
            ICNT = ICNT +1 
            F(ICNT) = CSH(N) 
          ENDIF 
          IF( IPLOT(12) ) THEN 
            ICNT = ICNT +1 
            F(ICNT) = CB(N) 
          ENDIF 
          IF( IPLOT(13) ) THEN 
            ICNT = ICNT +1 
            F(ICNT) = SDOWN(N) 
          ENDIF 
          IF( IPLOT(14) ) THEN 
            ICNT = ICNT +1 
            F(ICNT) = RSUP(N) 
          ENDIF 

!
!           
          IFORM = ICNT + 4 
          DO 90 M = 1,ICNT 
            IF( ABS(F(M)) .LE. 1.0D-99 ) F(M) = 0.0D+0 
   90     CONTINUE 
          IF( IFORM .LE. 9 ) THEN 
            WRITE (FORM1(14:14),'(I1)') IFORM 
            WRITE (IPL,FORM1) N,LinkName(L),L,I,X(N),(F(M),M=1,ICNT) 
          ELSE 
            WRITE ( FORM2(14:15),'(I2)') IFORM 
            WRITE (IPL,FORM2) N,LinkName(L),L,I,X(N),(F(M),M=1,ICNT) 
          ENDIF 
  100   CONTINUE 
  110 CONTINUE 
!----------------------------------------------------------------------C 
!     Close plot file. 
!----------------------------------------------------------------------C 
      CLOSE(UNIT=IPL) 
  200 CONTINUE 
!----------------------------------------------------------------------C 
!     Format Statements. 
!----------------------------------------------------------------------C 
 9110 FORMAT(                                                          &
       10X,'                  RIVTOX-DR Plot File                  ',  &
      /10X,'------------------ R I V T O X - D R ------------------') 
 9160 FORMAT(16X,'Execution Date: ',A10,                               &
          /16X,'Execution Time: ',A10,/) 
!----------------------------------------------------------------------C 
!     End of WRPLOT group. 
!----------------------------------------------------------------------C 
! 
      RETURN 
      END 
 
!======================================================================C 
! 
      SUBROUTINE RDRST 
! 
!----------------------------------------------------------------------C 
! 
!     RDRST: ReaD ReSTart files. 
! 
!----------------------------------------------------------------------C   
      USE DRSPFL
      USE DRWFLD
      USE FILINX 
      USE FLDINX 
      USE HEADNG
      USE LOGCLS  
      USE LANDSF
      USE SOLVAR
!----------------------------------------------------------------------C 
!     Implicit Double Precision. 
!----------------------------------------------------------------------C 
      IMPLICIT REAL*8 (A-H,O-Z) 
      IMPLICIT INTEGER*4 (I-N) 
!----------------------------------------------------------------------C 
!     Type Declarations. 
!----------------------------------------------------------------------C 
!      LOGICAL SCREEN,INITL,RESTART,ZSTART,IPLOT,IRNV,STEADY
      CHARACTER*240 CHDUM 
!----------------------------------------------------------------------C 
!     Write time and field variables to the restart file. 
!----------------------------------------------------------------------C  
      DO 100 N = 1,5 
        READ (IRS,'(A)') CHDUM 
  100 CONTINUE 
      IF( ZSTART ) THEN 
        READ (IRS,*) YDUM,ZDUM 
      ELSE 
        READ (IRS,'(3(1PE22.15))') TIME,DT,TMVL 
      ENDIF  
      
!      DO 200 N = 1,NPOINT 
!        READ (IRS,*) QH(N),HR(N),CH(N),SH(N),CSH(N),CB(N),TSF(N),ZSD(N) 
200   CONTINUE  
      
 if( (ISOLVE(2) == 1) .and. (ISOLVE(3)+ISOLVE(4) > 0) ) then
   read(IRS,*) (QH(ic),HR(ic),TSF(ic),SH(ic),ZSD(ic),CH(ic),CSH(ic),CB(ic), ic=1,NPOINT)
 elseif( ISOLVE(3) +ISOLVE(4) > 0 ) then
   read(IRS,*) (QH(ic),HR(ic),TSF(ic),tt,ZSD(ic),CH(ic),CSH(ic),CB(ic), ic=1,NPOINT)
 elseif( ISOLVE(2) == 1 ) then
   read(IRS,*) (QH(ic),HR(ic),TSF(ic),SH(ic),ZSD(ic),tt,tt,tt, ic=1,NPOINT)
 else
   read(IRS,*) (QH(ic),HR(ic),TSF(ic),tt,tt,tt,tt,tt, ic=1,NPOINT)
 endif

  150 continue   
  800 CONTINUE  
!----------------------------------------------------------------------C 
!     End of RDRST group. 
!----------------------------------------------------------------------C 
! 
      RETURN 
      END 
 
!======================================================================C 
! 
      SUBROUTINE WRRST 
! 
!----------------------------------------------------------------------C 
! 
!     WRRST: WRite ReSTart files. 
! 
!----------------------------------------------------------------------C  
      USE DRSPFL
      USE DRWFLD
      USE FILINX  
      USE FLDINX 
      USE HEADNG  
      USE LANDSF
      USE SOLVAR
!----------------------------------------------------------------------C 
!     Implicit Double Precision. 
!----------------------------------------------------------------------C 
      IMPLICIT REAL*8 (A-H,O-Z) 
      IMPLICIT INTEGER*4 (I-N) 
!----------------------------------------------------------------------C 
!     Type Declarations. 
!----------------------------------------------------------------------C 
      LOGICAL*4 fFExist
!----------------------------------------------------------------------!
!     Rename old restart file to restart_prev if exist.
!----------------------------------------------------------------------!
      fFExist=.false.
      INQUIRE(FILE=FileName(IRS), EXIST = fFExist)
      if( fFExist ) then  
        close(IRS)  
        !!! call System('move '//FileName(IRS)//' '//trim(FileName(IRS))//'_prev >nul') 
        CALL RENAME(TRIM(FileName(IRS)), TRIM(FileName(IRS))//'_prev')
      endif  
!   if (fFExist) call System('mv '//trim(DirName)//FileName(IRS)//' '//trim(DirName)//trim(FileName(IRS))//'_prev')
!----------------------------------------------------------------------!
!     Create restart file.
!----------------------------------------------------------------------!
      OPEN(UNIT=IRS, FILE=FileName(IRS), STATUS='REPLACE')
!----------------------------------------------------------------------C 
!     Delete old restart file and open a new restart file. 
!----------------------------------------------------------------------C 
!      CLOSE (UNIT=IRS,STATUS='DELETE') 
!      OPEN (UNIT=IRS,FILE=FileName(IRS),STATUS='NEW',FORM='FORMATTED' ) 
!----------------------------------------------------------------------C 
!     Write credits and timing information to the restart file. 
!----------------------------------------------------------------------C 
      WRITE (IRS,9210) 
      WRITE (IRS,9220) NDATE, NTIME 
!----------------------------------------------------------------------C 
!     Write time and field variables to the restart file. 
!----------------------------------------------------------------------C 
      WRITE (IRS,'(3(1PE22.15))') TIME,DT,TMVL
!      DO 200 N = 1,NPOINT 
!        WRITE (IRS,'(24(1PE22.15,1X))') QH(N),HR(N),                 &
!                CH(N),SH(N),CSH(N),CB(N),TSF(N),ZSD(N)   
200   CONTINUE   
      
 if( (ISOLVE(2) == 1) .and. (ISOLVE(3) +ISOLVE(4) > 0) ) then
   write(IRS,'(8(1PG22.15,1X))') (QH(ic),HR(ic),TSF(ic),SH(ic),ZSD(ic),CH(ic),CSH(ic),CB(ic), ic=1,NPOINT)
 elseif( ISOLVE(3) +ISOLVE(4) >0 ) then
   write(IRS,'(8(1PG22.15,1X))') (QH(ic),HR(ic),TSF(ic),0.d0,ZSD(ic),CH(ic),CSH(ic),CB(ic), ic=1,NPOINT)
 elseif( ISOLVE(2) == 1 ) then
   write(IRS,'(8(1PG22.15,1X))') (QH(ic),HR(ic),TSF(ic),SH(ic),ZSD(ic),0.d0,0.d0,0.d0, ic=1,NPOINT)
 else
   write(IRS,'(8(1PG22.15,1X))') (QH(ic),HR(ic),TSF(ic),0.d0,0.d0,0.d0,0.d0,0.d0, ic=1,NPOINT)
 endif
      
!----------------------------------------------------------------------C 
!     Format Statements. 
!----------------------------------------------------------------------C 
 9000 FORMAT(5(1PE22.15,1X)) 
 9210 FORMAT(                                                         &
      10X,'                  RIVTOX-DR Restart File                ', &
     /10X,'------------------- R I V T O X - D R --------------------') 
 9220 FORMAT(16X,'Execution Date: ',A10,                              &
           /16X,'Execution Time: ',A10,/) 
!----------------------------------------------------------------------C 
!     End of WRRST group. 
!----------------------------------------------------------------------C 
! 
      RETURN 
      END 
 
!======================================================================C 
! 
      SUBROUTINE RDOUTP 
! 
!----------------------------------------------------------------------C 
! 
!     RDOUTP: ReaD OUTPut control input group. 
! 
!----------------------------------------------------------------------C  
      USE CONSTS
      USE DRLINK
      USE DRWFLD 
      USE DRSPFL
      USE FILINX
      USE FLDINX  
      USE LANDSF
      USE LOGCLS
      USE NUMBRS
      USE OUTPLT
      USE POINTS 
      USE REFERN
      USE SOLVAR 
      USE EP_MOD 
      use intersub 
      USE PVSPNL
!      USE GLOBAL
!----------------------------------------------------------------------C 
!     Implicit Double Precision. 
!----------------------------------------------------------------------C 
      IMPLICIT REAL*8 (A-H,O-Z) 
      IMPLICIT INTEGER*4 (I-N) 
!----------------------------------------------------------------------C 
!     Type Declarations. 
!----------------------------------------------------------------------C 
      CHARACTER*240 CHDUM 
      CHARACTER*80 UNITS,ADUM 
      CHARACTER*500 str,str1
      CHARACTER*80 solfmt
!----------------------------------------------------------------------C 
!     WRITE header to output file. 
!----------------------------------------------------------------------C 
      WRITE (IWR,'(/A)') ' Output Control' 
      WRITE (IWR,'(A )') ' --------------' 
!----------------------------------------------------------------------C 
!     Found 'Output Control Record' 
!----------------------------------------------------------------------C 
!----------------------------------------------------------------------!
!     Read Number of Output Solution Files.
!----------------------------------------------------------------------!
 write(IWR,'(A)') ' Solution Files Output:'
 read(IRD,'(A)') str 
 call lcase(str) 
 istart=1
 ires = RDI4(',',istart,str,NFOut)
 write(IWR,'(A,I0)') ' Number of output solution files:  ',NFOut
!----------------------------------------------------------------------!
!     Allocate arrays.
!----------------------------------------------------------------------!
 ALLOCATE(IOutType(NFOut),IOutMode(NFOut)); IOutType=0; IOutMode=0;
 ALLOCATE(IOutFile(NFOut)); IOutFile='';
 ALLOCATE(NSolFrames(NVARS,NFOut)); NSolFrames=0;
 ALLOCATE(TMIN_OUT(NFOut),TMAX_OUT(NFOut),DT_OUT(NFOut), TOuts(NFOut));
  TMIN_OUT=0.d0; TMAX_OUT=0.d0; DT_OUT=0.d0; TOuts=0.d0;
 ALLOCATE(fFirstSolWrite(NFOut), FOuts(NFOut)); fFirstSolWrite=.false.; FOuts=.false.;
 ALLOCATE(SolStartT(NVARS,NFOut),SolEndT(NVARS,NFOut)); SolStartT=0.d0; SolEndT=0.d0;
 ALLOCATE(IPLOTB(50,NFOut)); IPLOTB=.false.;
!----------------------------------------------------------------------!
!     Read parameters for each output solution file.
!----------------------------------------------------------------------!
 do k=1,NFOut
   read(IRD,'(A)') str
   istart=1
   !----------------------------------------------------------------------!
   !     Read type.
   !----------------------------------------------------------------------!
   ires = RDStr(',',istart,str,str1)
   if (str1(1:4)=='plot') then
     IOutType(k)=1
     solfmt='plot'
   elseif (str1(1:7)=='bouss2d') then
     IOutType(k)=2
     solfmt='Bouss2DBinary'
   elseif (str1(1:3)=='sms') then
     IOutType(k)=3
     solfmt='SMSGenericBinary'
   elseif (str1(1:5)=='immsp') then
     IOutType(k)=4
     solfmt='IMMSPBinary'
   else
     call Msg(0,IWR,'INPUT ERROR! - Unrecognized Output Solution Files format: '//trim(str1))
     STOP
   endif
   if (k>1) write(IWR,*)
   write(IWR,'(I4,2A)') k,':  File format - ',trim(solfmt)
   !----------------------------------------------------------------------!
   !     Read output frequency.
   !----------------------------------------------------------------------!
   ires = RDR8(',',istart,str,TMIN_OUT(k))
   ires = RDR8(',',istart,str,TMAX_OUT(k))
   if (ires==1) TMAX_OUT(k)=BIG !--Assume infinite if empty or blank string
   ires = RDR8(',',istart,str,DT_OUT(k))
   ires = RDStr(',',istart,str,Units)
   rU=1.d0; call ConvToSI(rU,Units)
   TMIN_OUT(k)=TMIN_OUT(k)*rU
   TMAX_OUT(k)=TMAX_OUT(k)*rU
   DT_OUT(k)=DT_OUT(k)*rU
   if (DT_OUT(k)<SMALL) DT_OUT(k)=BIG
   tt = TMIN_OUT(k) + dint((TMAX_OUT(k)-TMIN_OUT(k))/DT_OUT(k))*DT_OUT(k)
   if ((tt+DT_OUT(k))-TMAX_OUT(k)>1.d-10) TMAX_OUT(k)=tt !--Check for precise timing (otherwise can give minus one time step bcz of stupid FPU errors)
   fFirstSolWrite(k)=.true. !--For solution files output
   write(IWR,'(A,3(1PG11.5,A))') '       Start Time - ',TMIN_OUT(k)/rU,',  End Time - ',TMAX_OUT(k)/rU,&
         ',  Interval - ',DT_OUT(k)/rU,',  Units - '//trim(Units)
 enddo

!----------------------------------------------------------------------!
!     Read extraction points from xy file.
!----------------------------------------------------------------------!
 read(IRD,'(A)') str
 istart=1
 ires = RDStr(',',istart,str,EPFName,1)
 EPFName = adjustl(EPFName)
 write(IWR,'(/2A)') ' Name of xy file containing extraction points:  ',trim(EPFName)
!----------------------------------------------------------------------!
!     Read extraction points output frequency and units.
!----------------------------------------------------------------------!
 call ReadFrequency(str,istart,IFQEPType,IFQEP,DTOUT_EP,'  Extraction Points')
!----------------------------------------------------------------------!
!     Read extraction points from XY file.
!----------------------------------------------------------------------!
 if (EPFName/='null') then
   call CalcTextFileLength(EPFName,EPN)
   EPN=max(0,EPN-1)
   if (EPN>0) then
     ALLOCATE(EPXY(EPN,2)); EPXY=0.d0;
     ALLOCATE(EPIJ(EPN,2)); EPIJ=0;
     ALLOCATE(EPName(EPN)); EPName='';
     call ReadConvTable(EPFName,' ',0,EPN,1,2, EPXY(:,1),EPXY(:,2))
     !-----Try to read gauge names-----!
     OPEN(UNIT=1001, FILE=EPFName)
      read(1001,*)
      do k=1,EPN
        read(1001,'(A)') str
        read(str,*,err=10,end=10) rr,rr,EPName(k)
        10 CONTINUE
      enddo
     CLOSE(1001)
     !-----Find the cells-----!
     do kl=1,EPN   
       EPIJ(kl,1) = EPXY(kl,1)   
       EPIJ(kl,2) = EPXY(kl,2)   
!!       call FindNearestGridNode(IFLD,JFLD,XBD,YBD,EPXY(kl,1),EPXY(kl,2),EPIJ(kl,1),EPIJ(kl,2))
     enddo 
     deallocate(EPXY);
   endif
   write(IWR,'(A,I0)')'  Number of entries in extraction points file:  ', EPN
!   write(IWR,'(A,A5)')'  Units of extraction points file:  ', OutlocUn
 endif

!----------------------------------------------------------------------C 
!     Read field variable output times. 
!     Read a record to see how many output times to read. 
!----------------------------------------------------------------------C 
      READ (IRD,*) NPRTM 
      ALLOCATE(PRTM(NPRTM)); PRTM=0.d0;
      write(IWR,'(/A,I0)') ' Number of output plot files:  ',NPRTM
      DO 210 II = 1,NPRTM 
        READ (IRD,'(A)') CHDUM 
        CALL LCASE( CHDUM ) 
        ISTART = 1 
        CALL RDDPR(ISTART,ICOMMA,CHDUM,PRTM(II)) 
        CALL RDCHR(ISTART,ICOMMA,CHDUM,UNITS) 
        write(IWR,'(A,1PG11.5,A)') '   At ',PRTM(ii),' '//trim(Units)
        CALL RDUNIT(UNITS,PRTM(II)) 
  210 CONTINUE 
!----------------------------------------------------------------------C 
!     Read reference node indices and output frequency. 
!----------------------------------------------------------------------C 
      READ (IRD,'(A)') CHDUM 
      CALL LCASE( CHDUM ) 
      ISTART = 1 
      NMREF = 0 
      DO 220 II = 1,4 
        CALL RDINT(ISTART,ICOMMA,CHDUM,IREF(II)) 
        CALL RDINT(ISTART,ICOMMA,CHDUM,JREF(II)) 
        IF( IREF(II) .NE. 0 .AND. JREF(II) .NE. 0 ) THEN 
          NMREF = NMREF +1 
          NREF(NMREF) = ND(IREF(II),JREF(II)) 
        ELSE 
          IREF(II) = 0 
          JREF(II) = 0 
          NREF(II) = 0 
        ENDIF 
  220 CONTINUE 
! 
      CALL RDINT(ISTART,ICOMMA,CHDUM,IFQREF) 
      IF( IFQREF .LT. 1 ) IFQREF = IBIG 
      write(IWR,'(A,I0,A)') ' Reference Node Output Frequency:  Every ', IFQREF,' Time Step(s)'
!      ires = RDI4(',',istart,str,IFQMV) 
!      if( IFQMV < 1 ) IFQMV = IBIG 
!      write(IWR,'(A,I0,A)') ' Max-Velocity Output Frequency  :  Every ', IFQMV,' Time Step(s)'
!      ires = RDI4(',',istart,str,IFQRST) 
      ires = RDI4(',',istart,chdum,IFQRST) 
      if( IFQRST < 1 ) IFQRST = IBIG 
      write(IWR,'(A,I0,A)') ' RESTART File Writing Frequency :  Every ', IFQRST,' Time Step(s)'
!----------------------------------------------------------------------C 
!     Initialize IPLOT(),IRNV() to false 
!----------------------------------------------------------------------C 
      IPLOTB(:,:) = .false.
      DO 230 I = 1,50 
        IPLOTEP(I) = .false. 
        IPLOT(I)   = .FALSE. 
        IRNV(I)    = .FALSE. 
  230 CONTINUE 
!----------------------------------------------------------------------!
!     Read output options for solution files.
!----------------------------------------------------------------------!
 do k=1,NFOut
   write(IWR,'(/A,I0,A)') ' Solution File ',k,' Variables'
   read(IRD,*) NLIN 
   do ii = 1,NLIN 
     read(IRD,'(A)') str 
     call lcase(str) 
     istart=1
     ires = RDStr(',',istart,str,str1)
     call FindVarName(str1,iv)
     if (iv==0) CYCLE
     IPLOTB(iv,k)=.true.
     write(IWR,'(3X,A)') trim(PNLU(iv))
   enddo
 enddo

!----------------------------------------------------------------------!
!     Read output options for extraction points.
!----------------------------------------------------------------------!
 write(IWR,'(/A)') ' Extraction Points Variables' 
 read(IRD,*) NLIN 
 do ii = 1,NLIN 
   read(IRD,'(A)') str 
   call lcase(str) 
   istart=1
   ires = RDStr(',',istart,str,str1)
   call FindVarName(str1,iv)
   if (iv==0) CYCLE
   IPLOTEP(iv)=.true.
   write(IWR,'(3X,A)') PNLU(iv)
 enddo

!----------------------------------------------------------------------C 
!     Read output options for plot file. 
!----------------------------------------------------------------------C 
      WRITE (IWR,'(/A/)') ' Plot File Variables' 
      READ (IRD,*) NLIN 
      DO 270 II = 1,NLIN 
        READ (IRD,'(A)') CHDUM 
        CALL LCASE( CHDUM ) 
        ISTART = 1 
        CALL RDCHR(ISTART,ICOMMA,CHDUM,ADUM) 
        WRITE (IWR,'(3X,A)' ) ADUM 
        IF( ADUM(1:13) .EQ. 'cross-section'   ) IPLOT( 1) = .TRUE. 
        IF( ADUM(1:12) .EQ. 'manning coef'    ) IPLOT( 2) = .TRUE. 
        IF( ADUM(1:12) .EQ. 'thalweg elev'    ) IPLOT( 3) = .TRUE. 
        IF( ADUM(1:13) .EQ. 'water surface'   ) IPLOT( 4) = .TRUE. 
        IF( ADUM(1:11) .EQ. 'water depth'     ) IPLOT( 5) = .TRUE. 
        IF( ADUM(1:14) .EQ. 'flow discharge'  ) IPLOT( 6) = .TRUE. 
        IF( ADUM(1:15) .EQ. 'top-layer depth' ) IPLOT( 7) = .TRUE. 
        IF( ADUM(1:14) .EQ. 'sediment conce'  ) IPLOT( 8) = .TRUE. 
        IF( ADUM(1:14) .EQ. 'sediment equil'  ) IPLOT( 9) = .TRUE. 
        IF( ADUM(1:14) .EQ. 'sol species co'  ) IPLOT(10) = .TRUE. 
        IF( ADUM(1:15) .EQ. 'part species co' ) IPLOT(11) = .TRUE. 
        IF( ADUM(1:14) .EQ. 'bot species co'  ) IPLOT(12) = .TRUE. 
        IF( ADUM(1:15) .EQ. 'deposition rate' ) IPLOT(13) = .TRUE. 
        IF( ADUM(1:12) .EQ. 'erosion rate'    ) IPLOT(14) = .TRUE. 
!        IF( ADUM(1:14) .EQ. 'x-dir. dischar'  ) IPLOT(13) = .TRUE. 
!        IF( ADUM(1:14) .EQ. 'y-dir. dischar'  ) IPLOT(14) = .TRUE. 
!        IF( ADUM(1:15) .EQ. 'y-dir. velocity' ) IPLOT( 7) = .TRUE. 
  270 CONTINUE 
!----------------------------------------------------------------------C 
!     Read reference node options for screen output. 
!----------------------------------------------------------------------C 
      WRITE (IWR,'(/A/)') ' Reference Node File Variables' 
      READ (IRD,*) NLIN 
      DO 310 II = 1,NLIN 
        READ (IRD,'(A)') CHDUM 
        CALL LCASE( CHDUM ) 
        ISTART = 1 
        CALL RDCHR(ISTART,ICOMMA,CHDUM,ADUM) 
        WRITE (IWR,'(3X,A)') ADUM  
        IF( ADUM(1:13) .EQ. 'cross-section'   ) IRNV( 1) = .TRUE. 
        IF( ADUM(1:12) .EQ. 'manning coef'    ) IRNV( 2) = .TRUE. 
        IF( ADUM(1:12) .EQ. 'thalweg elev'    ) IRNV( 3) = .TRUE. 
        IF( ADUM(1:13) .EQ. 'water surface'   ) IRNV( 4) = .TRUE. 
        IF( ADUM(1:11) .EQ. 'water depth'     ) IRNV( 5) = .TRUE. 
        IF( ADUM(1:14) .EQ. 'flow discharge'  ) IRNV( 6) = .TRUE. 
        IF( ADUM(1:15) .EQ. 'top-layer depth' ) IRNV( 7) = .TRUE. 
        IF( ADUM(1:14) .EQ. 'sediment conce'  ) IRNV( 8) = .TRUE. 
        IF( ADUM(1:14) .EQ. 'sediment equil'  ) IRNV( 9) = .TRUE. 
        IF( ADUM(1:14) .EQ. 'sol species co'  ) IRNV(10) = .TRUE. 
        IF( ADUM(1:15) .EQ. 'part species co' ) IRNV(11) = .TRUE. 
        IF( ADUM(1:14) .EQ. 'bot species co'  ) IRNV(12) = .TRUE. 
        IF( ADUM(1:15) .EQ. 'deposition rate' ) IRNV(13) = .TRUE. 
        IF( ADUM(1:12) .EQ. 'erosion rate'    ) IRNV(14) = .TRUE. 
!        IF( ADUM(1:14) .EQ. 'x-dir. dischar'  ) IRNV(13) = .TRUE. 
!        IF( ADUM(1:14) .EQ. 'y-dir. dischar'  ) IRNV(14) = .TRUE. 
!        IF( ADUM(1:15) .EQ. 'y-dir. velocity' ) IRNV( 7) = .TRUE. 
  310 CONTINUE 
!----------------------------------------------------------------------C 
!     Write output control data to output file. 
!----------------------------------------------------------------------C 
      WRITE (IWR,'(/A)') ' Field Variable Output Frequency:'   
      
      DO 400 II = 1,NPRTM 
        WRITE (IWR,'(A,E9.2,A)') '  At ',PRTM(II),' sec.' 
  400 CONTINUE 
!----------------------------------------------------------------------C 
!     Write reference node control data to output file. 
!----------------------------------------------------------------------C 
      DO 500 II = 1,NMREF 
        WRITE (IWR,'(3(A,I6))') ' Reference Node No.',II,': L = ',      &
         IREF(II),', I = ',JREF(II) 
  500 CONTINUE 
      WRITE (IWR,'(A,I6,A)') ' Reference Node Output Frequency: Every ',& 
       IFQREF,' Time Step(s)' 
      WRITE (IWR,'(//)') 
!----------------------------------------------------------------------C 
!     End of RDOUTP group. 
!----------------------------------------------------------------------C 
! 
      RETURN 
    END 
 
!======================================================================C 

