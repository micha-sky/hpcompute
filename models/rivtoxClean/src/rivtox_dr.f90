!======================================================================C 
! 
	PROGRAM RIVTOX_DR
!
!----------------------------------------------------------------------C
!
!     RIVTOX-DR: .
!
!     Version 3.0, IPMMS, Kiev, Ukraine (01.2016) 
!----------------------------------------------------------------------C  
      USE DRSPFL
      USE mGLOBAL 
      USE WTSOL  
      USE RW_MOD  
      USE REFERN 
      USE EP_MOD
      !USE IFPORT, ONLY : CLOCKX
!----------------------------------------------------------------------C
!     Implicit Double Precision.
!----------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
!----------------------------------------------------------------------C 
!     Type Declarations. 
!----------------------------------------------------------------------C 
!----------------------------------------------------------------------C 
!     Initialize convergence variables. 
!----------------------------------------------------------------------C
      CNVMX = 0.D+0 
      NFCNV = 0
      NRFW  = 0
      NRFS  = 0  
!----------------------------------------------------------------------C 
!     Define I/O file units and file names. 
!----------------------------------------------------------------------C
      ISF = 10
      IRD = 11 
      IWR = 12 
      IPL = 13 
      IRS = 14 
      FileName(ISF) = 'surface'
      FileName(IPL) = 'plot' 
      FileName(IRD) = 'input' 
      FileName(IWR) = 'output' 
      FileName(IRS) = 'restart' 

      OPEN( UNIT=IRD, FILE=FileName(IRD) ) 
      OPEN( UNIT=IWR, FILE=FileName(IWR),STATUS='UNKNOWN',           & 
                      FORM='FORMATTED') 
      CLOSE(UNIT=IWR, STATUS='DELETE') 
      OPEN( UNIT=IWR, FILE=FileName(IWR),STATUS='NEW',FORM='FORMATTED') 
      OPEN( UNIT=IRS, FILE=FileName(IRS) )
!----------------------------------------------------------------------C 
      iFUWTSol=21 !--Water solution files (21-30)
      iFUEP=2001 !--Extraction points (2001-3000)
!----------------------------------------------------------------------!
!     INITIALIZE ReadWriteRoutines.
!----------------------------------------------------------------------!
      RW_IFUTemp=1001
      RW_IFUEcho=IWR
!----------------------------------------------------------------------C 
!     Primary initializations. 
!----------------------------------------------------------------------C
      CALL INITAL 
!----------------------------------------------------------------------C 
!     Read input file. 
!----------------------------------------------------------------------C 
      CALL INPUT
!-----------------------------------------------------------------------! 
!     Find nodes for extraction points and create output files.
!-----------------------------------------------------------------------! 
      call InitExtractionPoints
!----------------------------------------------------------------------C
!     Compute cross-sectional area, width, and wetting perimeter.
!----------------------------------------------------------------------C
      CALL CRSSEC
!----------------------------------------------------------------------C
!     Load old time step arrays. 
!----------------------------------------------------------------------C 
      CALL LQPREV 
!----------------------------------------------------------------------C 
!     Reference node(s) output. 
!----------------------------------------------------------------------C  
      NITER = 0
      CALL REFNOD( NFCNV,NITER ) 
!----------------------------------------------------------------------C 
!     If initial conditions only requested, skip computations and quit. 
!----------------------------------------------------------------------C 
      IF( INITL ) GO TO 900 
!----------------------------------------------------------------------C 
!     Start time step loop. 
!----------------------------------------------------------------------C 
  100 CONTINUE 
!----------------------------------------------------------------------C 
!     Calculate time till next solution file output.
!----------------------------------------------------------------------C 
 TIMPR = BIG
 DTPL = BIG
 TOuts=BIG
 FOuts=.false.
 do k=1,NFOut
   if (TIME<TMIN_OUT(k)) then  
     TOuts(k) = TMIN_OUT(k)
   elseif (TIME<TMAX_OUT(k)) then  
     dNOUT = dint((TIME -TMIN_OUT(k))/DT_OUT(k)) +1
     TOuts(k) = TMIN_OUT(k) + dNOUT*DT_OUT(k)
     if (TOuts(k)-TIME<1.d-10) TOuts(k) = TOuts(k) + DT_OUT(k)
   endif
   TIMPR = min(TIMPR,TOuts(k))
 enddo
 do k=1,NFOut
   if (abs(TOuts(k)-TIMPR)<1.d-10) FOuts(k)=.true.
 enddo
 DTPL = TIMPR -TIME
!----------------------------------------------------------------------C 
!     Compute time remaining until next print output. 
!----------------------------------------------------------------------C 
      TIMPR = TIMAX  !!BIG 
      DO 200 NPR = 1,NPRTM 
        IF( TIME .LT. PRTM(NPR) ) TIMPR = MIN( TIMPR, PRTM(NPR) )  
  200 CONTINUE 
!!      DTPR = TIMPR - TIME 
      DTPL1 = TIMPR - TIME

 !-----Calc time till next extraction point output-----
 DTEP = BIG
 if (IFQEPType==2) then
   TOut = (dint(TIME/DTOUT_EP)+1.0)*DTOUT_EP
   if (TOut-TIME<1.d-10) TOut = TOut + DTOUT_EP !--Check for precise timing (otherwise can give TOut=TIME bcz of stupid FPU errors)
   DTEP = TOut - TIME
 endif

 !    Calculate time till next file output.
      DTPR = min(DTPL,DTPL1,DTEP)
      TIMPR = TIME +DTPR  
!----------------------------------------------------------------------C 
!     Store old time step size. 
!----------------------------------------------------------------------C 
      IF( DT*DTFAC .GT. DTPR ) DTOLD = DT
!----------------------------------------------------------------------C 
!     Compute time remaining until end of simulation. 
!----------------------------------------------------------------------C 
      DTQ = TIMAX -TIME 
!----------------------------------------------------------------------C 
!     If simulation time is completed, exit. 
!----------------------------------------------------------------------C 
      IF( DTQ .LE. 1.D-12 ) GO TO 900 
!----------------------------------------------------------------------C 
!     Select time step size. 
!----------------------------------------------------------------------C 
      DT = MIN( DTQ , DT*DTFAC, DTPR, DTMAX ) 
!----------------------------------------------------------------------C 
!     Advance simulation time by DT. 
!----------------------------------------------------------------------C 
      TIME = TIME +DT 
!----------------------------------------------------------------------C 
!     Increment time step counter. 
!----------------------------------------------------------------------C 
      NSTEP = NSTEP +1 
!      IF( NSTEP .GT. MXSTEP ) GO TO 900 
!----------------------------------------------------------------------!
!     If using Water Solutions then load next frames for each.
!     Then Calc water solutions for the moment.
!----------------------------------------------------------------------!
 if (NWTSol>0) then
   call LoadWaterSolutionsFrames
   call CalcWaterSolutionsAtTime
   CALL CRSSEC 
!----------------------------------------------------------------------!
!     If 1st step then init Before Vars.
!----------------------------------------------------------------------!
   if ((NSTEP==1).and.(ISOLVE(1)+ISOLVE(2) > 0)) call LQPREV
 endif
!----------------------------------------------------------------------C 
!     Restart at present time. 
!----------------------------------------------------------------------C 
  310 CONTINUE
!----------------------------------------------------------------------C 
!     Load old time step arrays. 
!----------------------------------------------------------------------C 
      IF( ISOLVE(1)+ISOLVE(2) > 0 ) CALL LQPREV 
!----------------------------------------------------------------------C 
!     Start of Saint Venant equations. 
!----------------------------------------------------------------------C 
      NDTR  = 0 
      NITER = 0 
!----------------------------------------------------------------------C 
!     Compute evapotranspiration rate. 
!----------------------------------------------------------------------C 
      CALL EVAPTR 
!----------------------------------------------------------------------C 
!     Compute precipitation fallout characteristics. 
!----------------------------------------------------------------------C 
      CALL RAINFL( 1 )
!----------------------------------------------------------------------C 
!     Start of iterations for transient problems. 
!----------------------------------------------------------------------C 
      IF( ISOLVE(1) .EQ. 0 ) GOTO 320 
  300 CONTINUE 
      NITER = NITER +1 
!----------------------------------------------------------------------C
!     Compute cross-sectional area, width, and wetting perimeter.
!----------------------------------------------------------------------C
      CALL CRSSEC 
!----------------------------------------------------------------------C 
!     Calculate lateral inflow. 
!----------------------------------------------------------------------C 
      IF( LQSOUR ) CALL LATINF
!----------------------------------------------------------------------C 
!     Proceeded forward sweep. 
!----------------------------------------------------------------------C 
      CALL COEFF  
      CALL MATRIX   
!----------------------------------------------------------------------C 
!     Insert nodal sources and sinks. 
!----------------------------------------------------------------------C 
      IF( LQSRND ) CALL NDLQSR
!----------------------------------------------------------------------C 
!     Modify the Jacobian matrix for non-zero flux boundary conditions. 
!----------------------------------------------------------------------C
      CALL BCLQ
!----------------------------------------------------------------------C 
!     Matrix Solver for Implicit Scheme. 
!----------------------------------------------------------------------C
      CALL SOLVL       
!----------------------------------------------------------------------C 
!     Proceeded backward sweep. 
!----------------------------------------------------------------------C 
      CALL BACKSW   
!----------------------------------------------------------------------C
!     Convergence test.
!----------------------------------------------------------------------C
      CALL CONVIT( NITER,NDTR,IRCODE,IERR )  
!            
      IF ( IRCODE .GE. 1) THEN 
        GOTO (300,310,900), IRCODE 
      ENDIF 
!
      IF( IBREAK .EQ. 1 ) THEN 
        IBREAK = 0 
      ELSEIF( IBREAK .GT. 1 ) THEN 
        IF( NITER .EQ. 1 ) THEN   
          IBREAK = IBREAK -1        
          DT = MIN(DTOLD,DT/DTR) 
        ELSE 
          IBREAK = IBREAK -1        
        ENDIF   
      ENDIF  
!----------------------------------------------------------------------C 
!     Set convergence indices and determine which equation has the
!     worst convergence. 
!----------------------------------------------------------------------C 
      CVW  = 0.D+0
      NRFW = 1
!----------------------------------------------------------------------C 
!     Calculate overland water discharges. 
!----------------------------------------------------------------------C
!      CALL DISCHG
!----------------------------------------------------------------------C 
!     Calculate liquid velocities. 
!----------------------------------------------------------------------C
!        IF( ISOLVE(1) .EQ. 1 ) CALL RNFLOW( RWORK(NLS),RWORK(NLS+LSX),  
!     &    RWORK(NLS+LSX+LSY),RWORK(NLS+2*LSX+LSY) ) 
  320 CONTINUE  
 www = hr(1) +hro(1) +qh(1) +qho(1) +tsf(1) +tsfo(1)
!----------------------------------------------------------------------C 
!     Runoff sediment transport equation. 
!----------------------------------------------------------------------C 
      IF( ISOLVE(2) == 1 ) THEN
!----------------------------------------------------------------------C 
!     Restore rainfall intensity. 
!----------------------------------------------------------------------C
!        CALL RAINFL( 2 ) 
!----------------------------------------------------------------------C 
!     Compute the equilibrium sediment concentration. 
!----------------------------------------------------------------------C 
        CALL SEDEQR 
!----------------------------------------------------------------------C 
!     Sediment transport equation. 
!----------------------------------------------------------------------C
        CALL SPEC( SH,SHO ) 
!----------------------------------------------------------------------C 
!     Incorporate the boundary conditions. 
!----------------------------------------------------------------------C
        CALL BCSEDM  
!----------------------------------------------------------------------C 
!     Incorporate sediment-bottom exchange. 
!----------------------------------------------------------------------C
        CALL SRCSED
!----------------------------------------------------------------------C 
!     Calculate lateral inflow. 
!----------------------------------------------------------------------C 
        IF( SDSOUR ) CALL SPLATR (3)
!----------------------------------------------------------------------C 
!     Insert nodal sources and sinks. 
!----------------------------------------------------------------------C 
        IF( SDSRND ) CALL NDSPSR (3)
!------------------------------------------------------------------------------C 
!     Store in ALC only non-zeros. 
!------------------------------------------------------------------------------C   
        CALL NOZERO( nPoint,nALC,NEL )    
!------------------------------------------------------------------------------C 
!     Sparse Matrix Solver. 
!------------------------------------------------------------------------------C  
        NVAR = nPoint
        if( .not.associated(RWORK) ) then   
          allocate( RWORK(nALC+9*NVAR) )  
        endif    
        if( .not.associated(IWORK) ) then   
          allocate( IWORK(12+nALC+5*NVAR) )  
        endif    
        CALL DSLUCS( NVAR,BLC,SH,NEL,IA,JA,ALC,0,2,1.D-13,500,ITER,               &
             ERR,IERR,0,RWORK,nALC+9*NVAR,IWORK,12+nALC+5*NVAR ) 
!      SUBROUTINE DSLUCS(N, B, X, NELT, IA, JA, A, ISYM, ITOL, TOL,
!     $     ITMAX, ITER, ERR, IERR, IUNIT, RWORK, LENW, IWORK, LENIW )
!----------------------------------------------------------------------C 
!     Land Surface Erosion/Deposition. 
!----------------------------------------------------------------------C  
        CALL EXNER 
      ENDIF  
!----------------------------------------------------------------------C 
!     Start of runoff sediment-species equations.  
!----------------------------------------------------------------------C  
      IF( ISOLVE(3)+ISOLVE(4) .GT. 0 ) THEN 
!----------------------------------------------------------------------C 
!     Load old time step arrays. 
!----------------------------------------------------------------------C 
        CALL SPPREV 
!----------------------------------------------------------------------C 
!     Runoff species transport equation. 
!----------------------------------------------------------------------C 
        IF( ISOLVE(3) .EQ. 1 ) THEN 
!----------------------------------------------------------------------C 
!     Species solute transport equation. 
!----------------------------------------------------------------------C 
          CALL SPEC( CH,CHO ) 
!----------------------------------------------------------------------C 
!     Incorporate the boundary conditions. 
!----------------------------------------------------------------------C 
          CALL BCSP
!----------------------------------------------------------------------C 
!     Calculate lateral inflow. 
!----------------------------------------------------------------------C 
          IF( CLSOUR ) CALL SPLATR (1)
!----------------------------------------------------------------------C 
!     Insert nodal sources and sinks. 
!----------------------------------------------------------------------C 
          IF( CLSRND ) CALL NDSPSR (1)
!------------------------------------------------------------------------------C 
!     Store in ALC only non-zeros. 
!------------------------------------------------------------------------------C   
          CALL NOZERO( nPoint,nALC,NEL )    
!----------------------------------------------------------------------C 
!     Sparse Matrix Solver. 
!----------------------------------------------------------------------C 
!          CALL DLUBCG( NFLD-NXP,IERR,2 ) 
          NVAR = nPoint
          if( .not.associated(RWORK) ) then   
            allocate( RWORK(nALC+9*NVAR) )  
          endif    
          if( .not.associated(IWORK) ) then   
            allocate( IWORK(12+nALC+5*NVAR) )  
          endif    
          CALL DSLUCS( NVAR,BLC,CH,NEL,IA,JA,ALC,0,2,1.D-13,500,ITER,             &
             ERR,IERR,0,RWORK,nALC+9*NVAR,IWORK,12+nALC+5*NVAR ) 
        ENDIF 
!----------------------------------------------------------------------C 
!     Particulate transport equation. 
!----------------------------------------------------------------------C  
        IF( ISOLVE(4) .EQ. 1 ) THEN  
          CALL SPEC( CSH,CSHO ) 
!----------------------------------------------------------------------C 
!     Incorporate the boundary conditions. 
!----------------------------------------------------------------------C 
          CALL BCPT   
!----------------------------------------------------------------------C 
!     Calculate lateral inflow. 
!----------------------------------------------------------------------C 
          IF( CPSOUR ) CALL SPLATR (2)
!----------------------------------------------------------------------C 
!     Insert nodal sources and sinks. 
!----------------------------------------------------------------------C 
          IF( CPSRND ) CALL NDSPSR (2)
!------------------------------------------------------------------------------C 
!     Store in ALC only non-zeros. 
!------------------------------------------------------------------------------C   
          CALL NOZERO( nPoint,nALC,NEL )    
!----------------------------------------------------------------------C 
!     Sparse Matrix Solver. 
!----------------------------------------------------------------------C 
!          CALL DLUBCG( NFLD-NXP,IERR,2 ) 
          NVAR = nPoint
          if( .not.associated(RWORK) ) then   
            allocate( RWORK(nALC+9*NVAR) )  
          endif    
          if( .not.associated(IWORK) ) then   
            allocate( IWORK(12+nALC+5*NVAR) )  
          endif    
          CALL DSLUCS( NVAR,BLC,CSH,NEL,IA,JA,ALC,0,2,1.D-13,500,ITER,            &
             ERR,IERR,0,RWORK,nALC+9*NVAR,IWORK,12+nALC+5*NVAR ) 
        ENDIF  
!----------------------------------------------------------------------C 
!     Upper Bottom Deposit Layer. 
!----------------------------------------------------------------------C 
        CALL UPLAYR
      ENDIF 
  800 CONTINUE 
!----------------------------------------------------------------------C
!     Call surface flux integrator.
!----------------------------------------------------------------------C
!      CALL SFINT
!----------------------------------------------------------------------C 
!     Diagnostic output at reference node(s). 
!----------------------------------------------------------------------C
      NFCNV = NRFW
      IF( MOD(NSTEP,IFQREF) .EQ. 0 ) CALL REFNOD( NFCNV,NITER ) 
!----------------------------------------------------------------------!
!     Write RESTART.
!----------------------------------------------------------------------!
      IF( MOD(NSTEP,IFQRST) .EQ. 0 ) THEN
       OPEN( UNIT=IRS, FILE=FileName(IRS) )
       CALL WRRST 
       CLOSE( IRS )
      ENDIF
!----------------------------------------------------------------------!
!     Write extraction points info.
!----------------------------------------------------------------------!
      if (IFQEPType==1) then
        if (mod(NSTEP,IFQEP)==0) call WriteExtractionPoints
      elseif (IFQEPType==2) then
        DTEP = DTEP-DT
        if (DTEP<=1.d-10) call WriteExtractionPoints
      endif
!----------------------------------------------------------------------C 
!     Write to output and plot files. 
!----------------------------------------------------------------------C 
      DTPR = TIMPR-TIME 
 DTPL = DTPL-DT
 DTPL1 = DTPL1-DT
 DTPR = min(DTPL,DTPL1)
 if (DTPR <= 1.D-12) then

   if (DTPL<=1.D-12) then 
       
     do k=1,NFOut
     if (FOuts(k)) then
       select case (IOutType(k))
       case (1)
         call WRPLOT(IPLOTB)
       case (2,3,4) !time-series files
         call WriteSolutionFiles(k)
       end select
     endif
     enddo
   endif

   if (DTPL1<=1.D-12) then
     call WRPLOT(IPLOT)
   endif

   DT = DTOLD
 endif
!      IF( DTPR .LE. 1.D-12 ) THEN 
!        CALL WRPLOT 
!        DT = DTOLD 
!      ENDIF 
!----------------------------------------------------------------------C 
!     New time step. 
!----------------------------------------------------------------------C 
      GO TO 100 
!----------------------------------------------------------------------C 
  900 CONTINUE   
!      CALL WRPLOT   
      call WRPLOT(IPLOT)
      WRITE (IWR,9000) 
      OPEN( UNIT=IRS, FILE=FileName(IRS) ) 
      CALL WRRST 
!----------------------------------------------------------------------C 
!     Close all files... 
!----------------------------------------------------------------------C 
      DO 910 N = 10,14 
        CLOSE (N) 
  910 CONTINUE  
      do i=1,NWTSol
        CLOSE(IWTSolIFU(i))
      enddo
      do i=1,EPN
        CLOSE(iFUEP+i-1)
      enddo
!----------------------------------------------------------------------C 
!     Terminate run here... 
!----------------------------------------------------------------------C 
      STOP '*** END OF SIMULATION ***' 
!----------------------------------------------------------------------C 
!     Format Statements. 
!----------------------------------------------------------------------C 
 9000 FORMAT(/,' END OF SIMULATION',/,' -----------------') 
!----------------------------------------------------------------------C 
!     End of RIVTOX_DR program. 
!----------------------------------------------------------------------C 
! 
      END 
 
!======================================================================C 
! 
      SUBROUTINE RDCHR( ISTART,ICOMMA,CHDUM,ADUM ) 
! 
!----------------------------------------------------------------------C 
! 
!     RDCHR: ReaD CHaRacter data. 
! 
!     MSTS (Multiphase Subsurface Transport Simulator) code.  
!     M.D. White, and W.E. Nichols. Pacific Northwest National Laboratory,  
!     Richland, WA, 1993. 
!----------------------------------------------------------------------C 
      USE FILINX
!----------------------------------------------------------------------C 
!     Implicit Double Precision. 
!----------------------------------------------------------------------C 
      IMPLICIT REAL*8 (A-H,O-Z) 
      IMPLICIT INTEGER*4 (I-N) 
!----------------------------------------------------------------------C 
!     Type Declarations. 
!----------------------------------------------------------------------C 
      CHARACTER*(*) ADUM 
      CHARACTER*(*) CHDUM 
!----------------------------------------------------------------------C 
!     Read characters between commas. 
!     Return 'null' for adjacent commas and for missing data. 
!     Issue a warning for missing data to alert user to use of defaults. 
!----------------------------------------------------------------------C 
      ICOMMA = INDEX (CHDUM(ISTART:), ',') + ISTART - 1 
      IF( ICOMMA .LT. ISTART ) THEN 
        WRITE (*,'(A)') ' Input Warning! Missing Character String Data' 
        WRITE (IWR,'(A)')' Input Warning! Missing Character String Data' 
        ADUM = 'null' 
        RETURN 
      ELSEIF( ICOMMA .EQ. ISTART ) THEN 
        ADUM = 'null' 
      ELSE 
   10   IF( CHDUM(ISTART:ISTART) .EQ. ' ' ) THEN 
          ISTART = ISTART + 1 
          GOTO 10 
        ENDIF    
        ADUM = ' ' 
        READ (CHDUM(ISTART:ICOMMA-1),'(A)') ADUM(1:ICOMMA-ISTART) 
      ENDIF 
      ISTART = ICOMMA +1 
!----------------------------------------------------------------------C 
!     End of RDCHR group. 
!----------------------------------------------------------------------C 
! 
      RETURN 
      END 
 
!======================================================================C 
! 
      SUBROUTINE RDDPR( ISTART,ICOMMA,CHDUM,DUMY ) 
! 
!----------------------------------------------------------------------C 
! 
!     RDDPR: ReaD Double Precision Real data. 
! 
!     MSTS (Multiphase Subsurface Transport Simulator) code.  
!     M.D. White, and W.E. Nichols. Pacific Northwest National Laboratory,  
!     Richland, WA, 1993. 
!----------------------------------------------------------------------C 
      USE FILINX
!----------------------------------------------------------------------C 
!     Implicit Double Precision. 
!----------------------------------------------------------------------C 
      IMPLICIT REAL*8 (A-H,O-Z) 
      IMPLICIT INTEGER*4 (I-N) 
!----------------------------------------------------------------------C 
!     Type Declarations. 
!----------------------------------------------------------------------C 
      CHARACTER*240 CHDUM 
      CHARACTER*6 FORM1 
      CHARACTER*7 FORM2 
!----------------------------------------------------------------------C 
!     Data Statements. 
!----------------------------------------------------------------------C 
      !!! SAVE FORM1,FORM2 
      !DATA FORM1 /'(D .0)'/ 
      !DATA FORM2 /'(D  .0)'/ 
      
      FORM1 = '(D .0)'
      FORM2 = '(D  .0)'
      
!----------------------------------------------------------------------C 
!     Read numbers between commas. 
!     Return '0.D+0' for adjacent commas and missing values. 
!     Issue a warning for missing data to alert user to use of defaults. 
!----------------------------------------------------------------------C 
      ICOMMA = INDEX (CHDUM(ISTART:140), ',') +ISTART -1 
      IF( ICOMMA .LT. ISTART ) THEN 
        WRITE (*,'(2A)') ' Input Warning! Missing Data: ',            &
         'Double Precision Real - Substituted a value of zero.' 
        WRITE (IWR,'(2A)') ' Input Warning! Missing Data: ',          &
         'Double Precision Real - Substituted a value of zero.' 
        DUMY = 0.D+0 
        RETURN 
      ELSEIF( ICOMMA .EQ. ISTART ) THEN 
        DUMY = 0.D+0 
      ELSE 
        NCHR = ICOMMA-ISTART 
        IF( NCHR .LT. 10 ) THEN 
          WRITE (FORM1(3:3),'(I1)') NCHR 
          READ (CHDUM(ISTART:ICOMMA-1), FORM1 ) DUMY 
        ELSEIF( NCHR .LT. 100 ) THEN 
          WRITE (FORM2(3:4),'(I2)') NCHR 
          READ (CHDUM(ISTART:ICOMMA-1), FORM2 ) DUMY 
        ELSE 
          WRITE (*,'(2A)') ' Input Error: Excessive Length: ',        &
           'Double Precision Real' 
          WRITE (IWR,'(2A)') ' Input Error: Excessive Length: ',      &
           'Double Precision Real' 
        ENDIF 
      ENDIF 
      ISTART = ICOMMA +1 
!----------------------------------------------------------------------C 
!     End of RDDPR group. 
!----------------------------------------------------------------------C 
! 
      RETURN 
      END 
 
!======================================================================C 
! 
      SUBROUTINE INPUT 
! 
!----------------------------------------------------------------------C 
! 
!     INPUT: read INPUT file. 
! 
!----------------------------------------------------------------------C 
      USE FILINX
!----------------------------------------------------------------------C 
!     Implicit Double Precision. 
!----------------------------------------------------------------------C 
      IMPLICIT REAL*8 (A-H,O-Z) 
      IMPLICIT INTEGER*4 (I-N) 
!----------------------------------------------------------------------C 
!      Type Declarations. 
!----------------------------------------------------------------------C 
      CHARACTER*240 CHDUM 
!----------------------------------------------------------------------C 
!     Loop over all input cards; logic will return to this location to 
!     read another group identification record. 
!----------------------------------------------------------------------C 
  100 CONTINUE  
! 
      READ (IRD,'(A)',END=500) CHDUM
      CALL LCASE( CHDUM ) 
! 
      IF( CHDUM( 1: 1) .NE. '~'         ) GOTO 100 
! 
      IF( CHDUM( 2:27) .EQ. 'simulation title and notes' ) CALL RDSIMU 
      IF( CHDUM( 2:17) .EQ. 'solution schemes'           ) CALL RDSOLU 
      IF( CHDUM( 2:13) .EQ. 'default para'               ) CALL RDDFPR
      IF( CHDUM( 2:19) .EQ. 'link topology data'         ) CALL RDLITO 
      IF( CHDUM( 2:19) .EQ. 'node topology data'         ) CALL RDNDTO 
      IF( CHDUM( 2:12) .EQ. 'points data'                ) CALL RDPOIN 
      IF( CHDUM( 2:14) .EQ. 'cross-section'              ) CALL RDSECT 
!      IF( CHDUM( 2: 15) .EQ. 'grid from file'      ) CALL RDGridFromFile
      IF( CHDUM( 2: 10) .EQ. 'dam input'                 ) CALL RDDAMS
!      IF( CHDUM( 2:13) .EQ. 'land surface'               ) CALL RDLSF 
      IF( CHDUM( 2:19) .EQ. 'upper bottom-layer'          ) CALL RDUpperSLDepth 
      IF( CHDUM( 2:18) .EQ. 'bottom rock types'          ) CALL RDROCK 
      IF( CHDUM( 2:13) .EQ. 'bottom types'               ) CALL RDBTTP 
      IF( CHDUM( 2: 5) .EQ. 'mech'                       ) CALL RDMECH 
      IF( CHDUM( 2:11) .EQ. 'bottom pro'                 ) CALL RDBTPR 
      IF( CHDUM( 2:12) .EQ. 'species pro'                ) CALL RDSPPR 
      IF( CHDUM( 2:13) .EQ. 'sediment pro'               ) CALL RDSEPR 
      IF( CHDUM( 2:13) .EQ. 'sediment bou'               ) CALL RDSEBC 
      IF( CHDUM( 2:12) .EQ. 'species bou'                ) CALL RDSPBC 
      IF( CHDUM( 2:14) .EQ. 'particulate b'              ) CALL RDPTBC 
      IF( CHDUM( 2:11) .EQ. 'node bound'                 ) CALL RDBCND 
      IF( CHDUM( 2:11) .EQ. 'link sourc'                 ) CALL RDLKSR 
      IF( CHDUM( 2:11) .EQ. 'node sourc'                 ) CALL RDNDSR 
      IF( CHDUM( 2: 5) .EQ. 'init'                       ) CALL RDINIT 
      if (chdum( 2:10)  ==  'water sol'                  ) call RDWaterSolutions
      IF( CHDUM( 2: 9) .EQ. 'rainfall'                   ) CALL RDRAIN 
      IF( CHDUM( 2:19) .EQ. 'evapotranspiration'         ) CALL RDEVAP 
!      IF( CHDUM( 2: 5) .EQ. 'wind'                       ) CALL RDWIND 
      IF( CHDUM( 2: 5) .EQ. 'outp'                       ) CALL RDOUTP 
!      IF( CHDUM( 2: 5) .EQ. 'surf'                       ) CALL RDSF
!----------------------------------------------------------------------C 
      GOTO 100 
  500 CONTINUE 
!----------------------------------------------------------------------C 
!     End of INPUT group. 
!----------------------------------------------------------------------C 
! 
      RETURN 
      END 
 
!======================================================================C 
! 
      SUBROUTINE RDINT( ISTART,ICOMMA,CHDUM,IDUMY ) 
! 
!----------------------------------------------------------------------C 
! 
!     RDINT: ReaD INTeger data. 
! 
!     MSTS (Multiphase Subsurface Transport Simulator) code.  
!     M.D. White, and W.E. Nichols. Pacific Northwest National Laboratory,  
!     Richland, WA, 1993. 
!----------------------------------------------------------------------C 
      USE FILINX
!----------------------------------------------------------------------C 
!     Implicit Double Precision. 
!----------------------------------------------------------------------C 
      IMPLICIT REAL*8 (A-H,O-Z) 
      IMPLICIT INTEGER*4 (I-N) 
!----------------------------------------------------------------------C 
!     Type Declarations. 
!----------------------------------------------------------------------C 
      CHARACTER*240 CHDUM 
      CHARACTER*4 FORM1 
!----------------------------------------------------------------------C 
!     Data Statements. 
!----------------------------------------------------------------------C 
      SAVE FORM1 
      DATA FORM1 /'(I )'/ 
!----------------------------------------------------------------------C 
!     Read numbers between commas. 
!     Return '0' for adjacent commas and for missing data. 
!     Issue a warning for missing data to alert user to use of defaults. 
!----------------------------------------------------------------------C 
      ICOMMA = INDEX (CHDUM(ISTART:), ',') + ISTART -1 
      IF( ICOMMA .LT. ISTART ) THEN 
        WRITE (*,'(2A)') ' Input Warning! Missing Data: ',            &
          'Integer- Zero Assumed.' 
        WRITE (IWR,'(2A)') ' Input Warning! Missing Data: ',          &
          'Integer- Zero Assumed.' 
        IDUMY = 0 
        RETURN 
      ELSEIF( ICOMMA .EQ. ISTART ) THEN 
        IDUMY = 0 
      ELSE 
        NCHR = ICOMMA - ISTART 
        IF( NCHR .LT. 10 ) THEN 
          WRITE (FORM1(3:3),'(I1)') NCHR 
          READ (CHDUM(ISTART:ICOMMA-1), FORM1 ) IDUMY 
        ELSE 
          WRITE (IWR,'(A)') ' Input Error: Excessive Length: Integer' 
          WRITE (*,'(A)') ' Input Error: Excessive Length: Integer' 
        ENDIF 
      ENDIF 
      ISTART = ICOMMA +1 
!----------------------------------------------------------------------C 
!     End of RDINT group. 
!----------------------------------------------------------------------C 
! 
      RETURN 
      END 
 
!======================================================================C 
! 
      SUBROUTINE RDSIMU 
! 
!----------------------------------------------------------------------C 
! 
!     RDSIMU: ReaD SIMUlation title input group. 
! 
!----------------------------------------------------------------------C 
      USE FILINX  
      USE HEADNG
      USE intersub
!----------------------------------------------------------------------C 
!     Implicit Double Precision. 
!----------------------------------------------------------------------C 
      IMPLICIT REAL*8 (A-H,O-Z) 
      IMPLICIT INTEGER*4 (I-N) 
!----------------------------------------------------------------------C 
!     Type Declarations. 
!----------------------------------------------------------------------C 
      CHARACTER*70 H(5) 
      INTEGER IARR(3) 
!----------------------------------------------------------------------C 
!     Data Statements. 
!----------------------------------------------------------------------C 
      SAVE H 
      DATA H /                                                        &
      'RRRR   III  V   V  TTTTT   OOO   X   X     DDDD   RRRR ',      &
      'R   R   I   V   V    T    O   O   X X      D   D  R   R',      & 
      'RRRR    I   V   V    T    O   O    X   --  D   D  RRRR ',      & 
      'R  R    I    V V     T    O   O   X X      D   D  R  R ',      &
      'R   R  III    V      T     OOO   X   X     DDDD   R   R'/ 
!----------------------------------------------------------------------C 
      ISTART = 1 
      IARR(1) = 4 
      IARR(2) = 18 
      IARR(3) = 1997 
      CALL IDATE (IARR) 
      IMONTH = IARR(1) 
      IDAY   = IARR(2) 
      IYEAR  = IARR(3) 
!----------------------------------------------------------------------C 
      WRITE (NDATE(1:2),'(I2)') IDAY 
      WRITE (NDATE(4:5),'(I2)') IMONTH 
      WRITE (NDATE(7:10),'(I4)') IYEAR 
      NDATE(3:3) = '/' 
      NDATE(6:6) = '/'  
      WRITE(NTIME(1:2),'(I2)')22 
      WRITE(NTIME(4:5),'(I2)')18 
      WRITE(NTIME(7:8),'(I2)')33 
      !CALL TIME( NTIME ) 
      NTIME(3:3) = ':' 
      NTIME(6:6) = ':' 
      NTIME(9:10) = '  '  
!----------------------------------------------------------------------C 
!     Write credits and timing information to the output file. 
!----------------------------------------------------------------------C 
      WRITE (IWR,920) 
      WRITE (IWR,940) NDATE, NTIME !NTIME GOES TO IWR NO-ONE CARES
!----------------------------------------------------------------------C 
!     Write credits and timing information to the 'screen' file. 
!----------------------------------------------------------------------C 
      WRITE (*,900) 
      DO 100 K=1,5 
        WRITE (*,910) H(K) 
  100 CONTINUE 
      WRITE (*,900) 
      WRITE (*,950)  
      WRITE (*,940) NDATE, NTIME 
      WRITE (*,900) 
!----------------------------------------------------------------------C 
!     WRITE header to output file. 
!----------------------------------------------------------------------C 
      WRITE (IWR,960) 
      WRITE (IWR,'(/A)') ' Simulation Title and Notes' 
      WRITE (IWR,'(A/)') ' --------------------------' 
!----------------------------------------------------------------------C 
!     Found 'Simulation Title Record' 
!----------------------------------------------------------------------C 
      READ (IRD,'(A)') TITLE 
      READ (IRD,'(A)') USER 
      READ (IRD,'(A)') COMP 
      READ (IRD,'(A)') DATVAL 
      READ (IRD,'(A)') TIMVAL  
!  
      WRITE (IWR,'(2A)') '   Simulation Title    : ', TITLE   
      WRITE (IWR,'(2A)') '   User Name           : ', USER  
      WRITE (IWR,'(2A)') '   Company Name        : ', COMP   
      WRITE (IWR,'(2A)') '   Input Creation Date : ', DATVAL  
      WRITE (IWR,'(2A)') '   Input Creation Time : ', TIMVAL  
      WRITE (IWR,'(/)') 
!----------------------------------------------------------------------C 
!     Read simulation notes. 
!----------------------------------------------------------------------C 
      READ (IRD,*) NUMLN 
      WRITE (IWR,'(A/)') ' Simulation notes: ' 
      DO 200 N=1, NUMLN 
        READ(IRD,'(A)') NOTES 
        WRITE(IWR,'(A)') NOTES 
  200 CONTINUE        
!----------------------------------------------------------------------C 
!     Format Statements. 
!----------------------------------------------------------------------C 
  900 FORMAT(/) 
  910 FORMAT(9X,70A) 
  920 FORMAT(                                                          &
         10X,'                 RIVTOX-DR Output File                '  &  
        /10X,'------------------ R I V T O X - D R ------------------' & 
        /) 
  940 FORMAT(16X,'Execution Date: ',A10,                               &
            /16X,'Execution Time: ',A10,/) 
  950 FORMAT(5X,                                                        &
      '------------------------- R I V T O X - D R ---------------------& 
      ----'      /) 
  960 FORMAT(//' Record of Input Data',                                &
             /' --------------------') 
!----------------------------------------------------------------------C 
!     End of RDSIMU group. 
!----------------------------------------------------------------------C 
! 
      RETURN 
      END 
 
!======================================================================C 
! 
      SUBROUTINE RDSOLU 
! 
!----------------------------------------------------------------------C 
! 
!     RDSOLU: ReaD SOLUtions schemes and numerical control input groups. 
! 
!----------------------------------------------------------------------C 
      USE FILINX 
      USE LOGCLS
      USE NUMBRS
      USE RESIDL
      USE SOLVAR
!----------------------------------------------------------------------C 
!     Implicit Double Precision. 
!----------------------------------------------------------------------C 
      IMPLICIT REAL*8 (A-H,O-Z) 
      IMPLICIT INTEGER*4 (I-N) 
!----------------------------------------------------------------------C 
!     Type Declarations. 
!----------------------------------------------------------------------C 
      CHARACTER*20 CHMEAN(10) 
      CHARACTER*80 ADUM,UNITS 
      CHARACTER*240 CHDUM 
!----------------------------------------------------------------------C 
!     WRITE header to output file. 
!----------------------------------------------------------------------C 
      WRITE (IWR,'(/A)') ' Solution Schemes and Options' 
      WRITE (IWR,'(A/)') ' ----------------------------' 
!----------------------------------------------------------------------C 
!     Define the maximum number of time step reductions and the time 
!     step reduction factor. 
!----------------------------------------------------------------------C 
      MXDTR = 30 
      DTR = 0.2 
!----------------------------------------------------------------------C 
!     Found 'Solution Schemes Record' 
!----------------------------------------------------------------------C 
      READ (IRD,'(A)') CHDUM 
      CALL LCASE( CHDUM ) 
!----------------------------------------------------------------------C 
!     Initial Condition and Restart option. 
!----------------------------------------------------------------------C 
      INITL   = .FALSE. 
      RESTART = .FALSE. 
      ISTART = 1 
      CALL RDCHR(ISTART,ICOMMA,CHDUM,ADUM) 
      IF( ADUM(1:4) .EQ. 'true' ) THEN 
         INITL = .TRUE. 
      ENDIF 
      CALL RDCHR(ISTART,ICOMMA,CHDUM,ADUM) 
      IF( ADUM(1:4) .EQ. 'true' ) THEN 
         RESTART = .TRUE. 
      ENDIF 
      CALL RDCHR(ISTART,ICOMMA,CHDUM,ADUM) 
      IF( ADUM(1:4) .EQ. 'true' ) THEN 
         ZSTART = .TRUE. 
      ENDIF 
      CALL RDCHR(ISTART,ICOMMA,CHDUM,ADUM) 
      IF( ADUM(1:4) .EQ. 'true' ) THEN 
         STEADY = .TRUE. 
      ENDIF 
! 
      IF( INITL ) THEN 
        WRITE (IWR,'(2A)') ' Simulation Control:   ',                &
         'Initial Conditions Only' 
      ELSE 
        WRITE (IWR,'(A)') ' Simulation Control:   Full Execution' 
      ENDIF 
! 
      IF( RESTART ) THEN 
        WRITE (IWR,'(A)') ' Simulation Type   :   Restart Simulation' 
      ELSE 
        WRITE (IWR,'(A)') ' Simulation Type   :   New Simulation' 
      ENDIF 
! 
      IF( STEADY ) THEN 
        WRITE (IWR,'(A)') ' Simulation Mode   :   Steady Flow' 
      ELSE 
        WRITE (IWR,'(A)') ' Simulation Mode   :   Unsteady Flow' 
      ENDIF 
!----------------------------------------------------------------------C 
!     Equations-to-solve options (10 logical options). 
!----------------------------------------------------------------------C 
      READ (IRD,'(A)') CHDUM 
      CALL LCASE( CHDUM ) 
      ISTART = 1 
      DO 240 I = 1,4 
         ISOLVE(I) = 0 
         CALL RDCHR(ISTART,ICOMMA,CHDUM,ADUM) 
         IF( ADUM(1:4) .EQ. 'true' ) ISOLVE(I) = 1 
  240 CONTINUE  
!      ISOLVE(1) = ISOLVE(2) 
!      ISOLVE(2) = 0 
! 
      WRITE(IWR,'(/)')  
      WRITE(IWR,'(A)') 'Solved Equations:' 
      IF( ISOLVE(1) .GE. 1 ) THEN 
        WRITE (IWR,'(A)') ' River Water Equations  :   solution on' 
      ELSE 
        WRITE (IWR,'(A)') ' River Water Equations  :   solution off' 
      ENDIF 
! 
      IF( ISOLVE(2) .GE. 1 ) THEN 
        WRITE (IWR,'(A)') ' Sediment Transport       :   solution on' 
      ELSE 
        WRITE (IWR,'(A)') ' Sediment Transport       :   solution off' 
      ENDIF 
! 
      IF( ISOLVE(3) .GE. 1 ) THEN 
        WRITE (IWR,'(A)') ' Solute Transport Equation:   solution on' 
      ELSE 
        WRITE (IWR,'(A)') ' Solute Transport Equation:   solution off' 
      ENDIF 
! 
      IF( ISOLVE(4) .GE. 1 ) THEN 
        WRITE (IWR,'(A)') ' Particulate Transport    :   solution on' 
      ELSE 
        WRITE (IWR,'(A)') ' Particulate Transport    :   solution off' 
      ENDIF 
!----------------------------------------------------------------------C 
!     Maximum time to simulate. 
!----------------------------------------------------------------------C 
      READ (IRD,'(A)') CHDUM 
      CALL LCASE( CHDUM ) 
      ISTART = 1 
      CALL RDDPR(ISTART,ICOMMA,CHDUM,TIMAX) 
      CALL RDCHR(ISTART,ICOMMA,CHDUM,UNITS) 
      CALL RDUNIT(UNITS,TIMAX) 

!----------------------------------------------------------------------C 
!     Time scale for output. 
!----------------------------------------------------------------------C 
      TSCL = 1.D+0 
      CALL RDUNIT(UNITS,TSCL) 
!----------------------------------------------------------------------C 
!     Initial time step size. 
!----------------------------------------------------------------------C 
      READ (IRD,'(A)') CHDUM 
      CALL LCASE( CHDUM ) 
      ISTART = 1 
      CALL RDDPR(ISTART,ICOMMA,CHDUM,DT) 
      CALL RDCHR(ISTART,ICOMMA,CHDUM,UNITS) 
      CALL RDUNIT(UNITS,DT) 
!----------------------------------------------------------------------C 
!     Time step size acceleration factor. 
!----------------------------------------------------------------------C 
      READ (IRD,*) DTFAC 
!----------------------------------------------------------------------C 
!     Maximum time step size. 
!----------------------------------------------------------------------C 
      READ (IRD,'(A)') CHDUM 
      CALL LCASE( CHDUM ) 
      ISTART = 1 
      CALL RDDPR(ISTART,ICOMMA,CHDUM,DTMAX) 
      CALL RDCHR(ISTART,ICOMMA,CHDUM,UNITS) 
      CALL RDUNIT(UNITS,DTMAX)  
!----------------------------------------------------------------------C 
!     Echo time controls. 
!----------------------------------------------------------------------C 
      WRITE (IWR,'(/A)') ' Time Advancement and Limits' 
      WRITE (IWR,'(A,1PE11.4)')  '   Maximum Time:                  ', & 
        TIMAX 
      IF( .NOT. RESTART )                                              &
        WRITE (IWR,'(A,1PE11.4)')'   Initial Time Step:             ', &
           DT 
        WRITE (IWR,'(A,1PE11.4)')'   Time Step Acceleration Factor: ', &
           DTFAC  

      MNSTEP = 1 
      IF( TIMAX/DT .LE. 32000. ) MXSTEP = 5*(INT(TIMAX/DT) + 1) 
!----------------------------------------------------------------------C 
!     Adjust DT by DTFAC. 
!----------------------------------------------------------------------C 
      IF( DTFAC .NE. ZERO ) DT = DT / DTFAC 
!----------------------------------------------------------------------C 
!     Decipher other options. 
!----------------------------------------------------------------------C 
      READ (IRD,'(A)') CHDUM 
      CALL LCASE( CHDUM ) 
      SCREEN  = .FALSE. 
! 
      ISTART = 1 
      CALL RDCHR(ISTART,ICOMMA,CHDUM,ADUM) 
      IF( ADUM(1:4) .EQ. 'true' ) SCREEN = .TRUE. 
! 
      WRITE (IWR,'(/A)')  ' Solution Options' 
! 
      IF( SCREEN ) THEN 
        WRITE (IWR,'(A)') '   Screen Printing           : On' 
      ELSE 
        WRITE (IWR,'(A)') '   Screen Printing           : Off' 
      ENDIF  
!----------------------------------------------------------------------C 
!     Read Numerical Control card header. 
!----------------------------------------------------------------------C 
      READ (IRD,'(A)') CHDUM 
!----------------------------------------------------------------------C 
!     WRITE header to output file. 
!----------------------------------------------------------------------C 
      WRITE (IWR,'(/A)') ' Numerical Control' 
      WRITE (IWR,'(A )') ' -----------------' 
!----------------------------------------------------------------------C 
!     Read Newton-Raphson Iteration Limit. 
!----------------------------------------------------------------------C 
      READ (IRD,'(A)') CHDUM 
      CALL LCASE(CHDUM) 
      ISTART = 1 
      CALL RDINT(ISTART,ICOMMA,CHDUM,NRSD) 
      CALL RDCHR(ISTART,ICOMMA,CHDUM,ADUM) 
      IF( ADUM(1:5) .EQ. 'least' .OR. ADUM(1:1) .EQ. '1' ) THEN 
        IRSD = 1 
      ELSEIF( ADUM(1:7) .EQ. 'maximum' .OR. ADUM(1:1) .EQ. '2'  ) THEN 
        IRSD = 2 
      ENDIF 
!----------------------------------------------------------------------C 
!     Read convergence test mode and convergence limits. 
!----------------------------------------------------------------------C 
      READ (IRD,'(A)') CHDUM 
      ISTART = 1 
      DO I=1,4 
        RSD(I)   = 0.D+0 
        CALL RDDPR(ISTART,ICOMMA,CHDUM,RSDMX(I)) 
      ENDDO 
!----------------------------------------------------------------------C 
!     Write convergence information. 
!----------------------------------------------------------------------C 
      WRITE (IWR,'(/A)')   ' Convergence criteria' 
      WRITE (IWR,'(A,I4)') '  Maximum Number of Iterations per Step : ',& 
        NRSD 
      IF( IRSD .EQ. 1 ) THEN 
        WRITE (IWR,'(2A)') '  Convergence Mode Index: ',                &
          'Least Squares Summation' 
      ELSE 
        WRITE (IWR,'(2A)') '  Convergence Mode Index: ',                &
          'Maximum Residual' 
      ENDIF 
!      WRITE (IWR,'(A,1PE11.4)' )  
!     &           ' X-Dir. Water Velocity Residual  : ',RSDMX(1) 
!      WRITE (IWR,'(A,1PE11.4)' )  
!     &           ' Y-Dir. Water Velocity Residual  : ',RSDMX(2) 
      WRITE (IWR,'(A,1PE11.4)' )                                       &
                '  Water Depth Residual            : ',RSDMX(1) 
      WRITE (IWR,'(A,1PE11.4)' )                                       &
                '  Sediment Transport Residual     : ',RSDMX(2) 
      WRITE (IWR,'(A,1PE11.4)' )                                       &
                '  Species Transport Residual      : ',RSDMX(3) 
      WRITE (IWR,'(A,1PE11.4)' )                                       &
                '  Particulate Transport Residual  : ',RSDMX(4) 
!----------------------------------------------------------------------C 
!   Read 'mean' options 
!----------------------------------------------------------------------C 
!      READ (IRD,'(A)') CHDUM 
!      CALL LCASE( CHDUM ) 
!      ISTART = 1 
!      DO 380 I = 1, 4 
!        CALL RDCHR(ISTART,ICOMMA,CHDUM,ADUM) 
!        IF( ADUM(1:4) .EQ. 'harm' ) THEN 
!          IDMEAN(I) = 1 
!        ELSEIF( ADUM(1:4) .EQ. 'geom' ) THEN 
!          IDMEAN(I) = 2 
!        ELSEIF( ADUM(1:4) .EQ. 'arit' ) THEN 
!          IDMEAN(I) = 3 
!        ELSEIF( ADUM(1:4) .EQ. 'upwi' ) THEN 
!          IDMEAN(I) = 4 
!        ENDIF 
!        CHMEAN(I) = ADUM(1:20) 
!        ISTART = ICOMMA +1 
!  380 CONTINUE 
!      WRITE (IWR,'(/A)') ' Node Boundary Averaging Schemes' 
!      WRITE (IWR,'(2A)') ' Shallow Water Equations     :   ',CHMEAN(1) 
!      WRITE (IWR,'(2A)') ' Sediment Transport Equation :   ',CHMEAN(2) 
!      WRITE (IWR,'(2A)') ' Species Transport Equation  :   ',CHMEAN(3) 
!      WRITE (IWR,'(2A)') ' Particulate Species Equation:   ',CHMEAN(4) 
!----------------------------------------------------------------------C 
!     End of RDSOLU group. 
!----------------------------------------------------------------------C 
! 
      RETURN 
      END 
 
!======================================================================C 
! 
      SUBROUTINE RDUNIT( UNITS,VALUE ) 
! 
!----------------------------------------------------------------------C 
! 
!     RDUNIT: ReaD UNITs and convert value passed in VALUE from input 
!             units to SI base units. 
! 
!     Rules for expressing units: 
!       Units may only contain one divisor. 
!       Units may only contain one colon. 
!       Permeability units expressed as hydraulic conductivity (L/T) 
!         instead of intrinsic permeability (1/L^2) may be preceded 
!         with an "hc" indicator and a colon separator. 
!       Components within a units character string must be separated 
!         with a blank or a divisor (or a colon for hydraulic  
!         conductivity indications). 
!       No spaces between component name and divisor. 
!       Units raised to powers are indicated with a '^' symbol. 
! 
!       Examples: 
!         'Btu in/h ft^2 F'  converts to SI units of  'W/m k' 
!         'lbm/h ft'         converts to SI units of  'Pa s' 
!         'gm/l'             converts to SI units of  'kg/m^3' 
!         'hc:m/s'           converts to SI units of  'm^2' 
! 
!     MSTS (Multiphase Subsurface Transport Simulator) code.  
!     M.D. White, and W.E. Nichols. Pacific Northwest National Laboratory,  
!     Richland, WA, 1993. 
!----------------------------------------------------------------------C 
      USE FILINX
!----------------------------------------------------------------------C 
!     Implicit Double Precision. 
!----------------------------------------------------------------------C 
      IMPLICIT REAL*8 (A-H,O-Z) 
      IMPLICIT INTEGER*4 (I-N) 
!----------------------------------------------------------------------C 
!     Parameter Statements. 
!----------------------------------------------------------------------C 
      PARAMETER (LUNS=55) 
!----------------------------------------------------------------------C 
!     Type Declarations. 
!----------------------------------------------------------------------C 
      CHARACTER*80 UNITS 
      CHARACTER*10 CHS(LUNS),CHD 
      REAL*8 CF(LUNS)
      CHARACTER FORM2*6 
!----------------------------------------------------------------------C 
!     Data Statements. 
!----------------------------------------------------------------------C 
      SAVE CHS,CF,FORM2
!      SAVE FORM1
      DATA CHS /'m','kg','s','j','celsius','pa','w','kgmol','rad',    &
               'liquid','liq.','gas','solid','sol.','soil',           &
               'mm','cm','in','ft','yd','km',                         &
               'liter','l','gal',                                     &
               'gm','lbm','slug',                                     &
               'min','h','day','wk','fortnight','yr',                 &
               'btu','cal','hp',                                      &
               'kelvin','fahrenheit','rankine',                       & 
               'k', 'f', 'r',                                         &
               'psi','bar','atm',                                     &
               'degrees','deg',                                       &
               'cp','p','hc','1','mol','lbmol','debyes',              &
               '(absolute)'/ 
      DATA CF  /1.D+0,1.D+0,1.D+0,1.D+0,1.D+0,1.D+0,1.D+0,1.D+0,1.D+0, & 
               1.D+0,1.D+0,1.D+0,1.D+0,1.D+0,1.D+0,                    &
               1.D-3,1.D-2,2.54D-2,0.3048D+0,0.9144D+0,1.D+3,          &
               1.D-3,1.D-3,3.7854D-3,                                  &
               1.D-3,4.5359D-1,1.4594D+1,                              &
               6.D+1,3.6D+3,8.64D+4,6.048D+5,1.2096D+6,3.15576D+7,     &
               1.0544D+3,4.184D+0,7.457D+2,                            &
               1.D+0,5.555556D-1,5.555556D-1,                          &
               1.D+0,5.555556D-1,5.555556D-1,                          &
               6894.8D+0,1.D+5,1.01325D+5,                             &
               0.01745,0.0174533,                                      &
               1.D-3,1.D-1,1.02286D-7,1.D+0,1.D-3,4.5359D-1,1.D+0,     &
               1.D+0/ 
!      DATA FORM1 /'(I )'/ 
      DATA FORM2 /'(F .0)'/ 
      
!----------------------------------------------------------------------C 
      IF( TRIM(UNITS) .EQ. 'null' .OR. TRIM(UNITS) .EQ. 'none' ) RETURN 
!----------------------------------------------------------------------C 
!     Decompose the units into components and convert individual  
!     components. 
!----------------------------------------------------------------------C 
      IS = 1 
      IDV = INDEX( UNITS(1:),'/'  ) -1 
      IE  = INDEX( UNITS(1:),'  ' ) -1 
!----------------------------------------------------------------------C 
!     Units without a divisor. 
!----------------------------------------------------------------------C 
      IF( IDV .EQ. -1 ) THEN 
  100   CONTINUE 
        ISP = INDEX( UNITS(IS:),' ' ) +IS -2 
        IB = MIN( IE,ISP ) 
        CHD = UNITS(IS:IB) 
        IC = INDEX( CHD(1:),'^' ) 
        IF( IC .EQ. 0 ) THEN 
          IP = 1
          PW = 1.D+0  
        ELSE 
          I1 = IC +1 
          I2 = IB -IS +1 
          I3 = I2 -I1 +1 
          WRITE (FORM2(3:3),'(I1)') I3 
          READ (CHD(I1:I2),FORM2) PW 
!          WRITE (FORM1(3:3),'(I1)') I3 
!          READ (CHD(I1:I2),FORM1) IP 
          I2 = IC -1 
          CHD = CHD(1:I2) 
        ENDIF 
        DO 110 N = 1,LUNS 
          IF( CHS(N) .EQ. CHD ) THEN 
            VALUE = VALUE*(CF(N)**PW) 
!            VALUE = VALUE*(CF(N)**IP) 
            GOTO 120 
          ENDIF 
  110   CONTINUE 
        WRITE (IWR,'(/2A)') ' ERROR: Unrecognized Input Units: ',UNITS 
        WRITE (*,'(/2A)') ' ERROR: Unrecognized Input Units: ',UNITS 
        STOP 
  120   CONTINUE 
        IF( IB .LT. IE ) THEN 
          IS = IB +2 
          GOTO 100 
        ENDIF 
!----------------------------------------------------------------------C 
!     Units with a divisor. 
!----------------------------------------------------------------------C 
      ELSE 
!----------------------------------------------------------------------C 
!     Components before the divisor. 
!----------------------------------------------------------------------C 
  200   CONTINUE  
        ISP = INDEX( UNITS(IS:),' ' ) +IS -2 
        ICO = INDEX( UNITS(IS:),':' ) +IS -2 
        IF( (ICO .GT. 0) .AND. (ICO .GT. IS) ) THEN 
          IB = MIN( IDV,ISP,ICO ) 
        ELSE 
          IB = MIN( IDV,ISP ) 
        ENDIF 
        CHD = UNITS(IS:IB) 
        IC = INDEX( CHD(1:),'^' ) 
        IF( IC .EQ. 0 ) THEN 
          IP = 1
          PW = 1.D+0 
        ELSE 
          I1 = IC +1 
          I2 = IB -IS +1 
          I3 = I2 -I1 +1 
          WRITE (FORM2(3:3),'(I1)') I3 
          READ (CHD(I1:I2),FORM2) PW 
!          WRITE (FORM1(3:3),'(I1)') I3 
!          READ (CHD(I1:I2),FORM1) IP 
          I2 = IC -1 
          CHD = CHD(1:I2) 
        ENDIF 
        DO 210 N = 1,LUNS 
          IF( CHS(N) .EQ. CHD ) THEN 
            VALUE = VALUE*(CF(N)**PW) 
!            VALUE = VALUE*(CF(N)**IP) 
            GOTO 220 
          ENDIF 
  210   CONTINUE 
        WRITE (IWR,'(/2A)') ' ERROR: Unrecognized Input Units: ',UNITS 
        WRITE (*,'(/2A)') ' ERROR: Unrecognized Input Units: ',UNITS 
        STOP 
  220   CONTINUE 
        IF( IB .LT. IDV ) THEN 
          IS = IB +2 
          GOTO 200 
        ELSE 
          IS = IB +2 
          GOTO 300 
        ENDIF 
!----------------------------------------------------------------------C 
!     Components after the divisor. 
!----------------------------------------------------------------------C 
  300   CONTINUE 
        ISP = INDEX( UNITS(IS:),' ' ) +IS -2 
        IB = MIN( IE,ISP ) 
        CHD = UNITS(IS:IB) 
        IC = INDEX( CHD(1:),'^' ) 
        IF( IC .EQ. 0 ) THEN 
          IP = 1
          PW = 1.D+0 
        ELSE 
          I1 = IC+1 
          I2 = IB-IS+1 
          I3 = I2-I1+1 
          WRITE (FORM2(3:3),'(I1)') I3 
          READ (CHD(I1:I2),FORM2) PW 
!          WRITE (FORM1(3:3),'(I1)') I3 
!          READ (CHD(I1:I2),FORM1) IP 
          I2 = IC-1 
          CHD = CHD(1:I2) 
        ENDIF 
        DO 310 N = 1,LUNS 
          IF( CHS(N) .EQ. CHD ) THEN 
            VALUE = VALUE/(CF(N)**PW) 
!            VALUE = VALUE/(CF(N)**IP) 
            GOTO 320 
          ENDIF 
  310   CONTINUE 
        WRITE (IWR,'(/2A)') ' ERROR: Unrecognized Input Units: ',UNITS 
        WRITE (*,'(/2A)') ' ERROR: Unrecognized Input Units: ',UNITS 
        STOP 
  320   CONTINUE 
        IF( IB .LT. IE ) THEN 
          IS = IB +2 
          GOTO 300 
        ENDIF 
      ENDIF 
!----------------------------------------------------------------------C 
!     End of RDUNIT group. 
!----------------------------------------------------------------------C 
! 
      RETURN 
      END 
!======================================================================C 
! 
      SUBROUTINE LCASE( CHDUM ) 
! 
!----------------------------------------------------------------------C 
! 
!     LCASE: Lower CASE text string converter. 
! 
!     Generic subroutine to convert all upper-case characters in a 
!     variable-length string variable to lower case.  This subroutine 
!     does not disturb non-alphabetic characters; only captial letters 
!     (ASCII 65 through 90) are modified. 
! 
!     MSTS (Multiphase Subsurface Transport Simulator) code.  
!     M.D. White, and W.E. Nichols. Pacific Northwest National Laboratory,  
!     Richland, WA, 1993. 
!----------------------------------------------------------------------C 
!     Implicit Double Precision. 
!----------------------------------------------------------------------C 
      IMPLICIT REAL*8 (A-H,O-Z) 
      IMPLICIT INTEGER*4 (I-N) 
!----------------------------------------------------------------------C 
!     Type Declarations. 
!----------------------------------------------------------------------C 
      CHARACTER*(*) CHDUM  
!----------------------------------------------------------------------C 
!     Determine length of the character string. 
!----------------------------------------------------------------------C 
      LENGTH = LEN(CHDUM) 
!----------------------------------------------------------------------C 
!     Convert each upper-case alphabetic character to lower case. 
!----------------------------------------------------------------------C 
      DO 10 M = 1,LENGTH 
        IDUM = ICHAR(CHDUM(M:M)) 
        IF( IDUM .GE. 65 .AND. IDUM .LE. 90 ) THEN 
          IDUM = IDUM +32 
          CHDUM(M:M) = CHAR(IDUM) 
        ENDIF 
   10 CONTINUE 
!----------------------------------------------------------------------C 
!     End of LCASE group. 
!----------------------------------------------------------------------C 
! 
      RETURN 
      END 
 
!======================================================================C 
! 
      SUBROUTINE INITAL 
! 
!----------------------------------------------------------------------C 
! 
!     INITAL: variable INITiALization. 
! 
!----------------------------------------------------------------------C 
      USE LOGCLS
      USE mGLOBAL 
      USE DRCATC 
      USE NUMBRS 
      USE WORK
!----------------------------------------------------------------------C 
!     Implicit Double Precision. 
!----------------------------------------------------------------------C 
      IMPLICIT REAL*8 (A-H,O-Z) 
      IMPLICIT INTEGER*4 (I-N) 
!----------------------------------------------------------------------C 
!     Character and logical declarations. 
!----------------------------------------------------------------------C 
!----------------------------------------------------------------------C 
!     Small and big numbers. 
!----------------------------------------------------------------------C 
      SMALL  = 1.D-20 
      BIG    = 1.D+20 
      ISMALL = -32000 
      IBIG   = 32000 
      ZERO   = 0.D+0 
      IBREAK = 0   
      HSmall = 1.d-4
!----------------------------------------------------------------------C 
      EPSMIN = 1.D+0 
  100 EPSMIN = 0.5D+0*EPSMIN 
      IF( 1.D+0+EPSMIN .NE. 1.D+0 ) GOTO 100 
      EPSMIN = 2.D+0*EPSMIN 
!----------------------------------------------------------------------C 
!     Reference properties of water: units are [ kg, m, s, K, J ] 
!     Water values suitable for 20 deg C. 
!     kg/m**3, W/m K, J/kg K 
!----------------------------------------------------------------------C 
      RHOLQ  = 9.9832238583783692D+2 
      SPHF   = 4.182D+3 
      VISLQ  = 1.0017364807268300D-3 
      GRAV   = 9.81D+0 
      PATM   = 1.013255D+5 
      PCR    = 2.212D+7 
      PMN    = 6.1125D+2  
      PI     = 3.14159265358979D+0  
!----------------------------------------------------------------------C 
!     Time step factors etc. 
!----------------------------------------------------------------------C 
      NSTEP  = 0 
      TIME   = 0.D+0 
      DT     = 1.0D+0 
      DTOVER = 0.D+0 
      DTFAC  = 1. 
      DTMAX  = BIG 
      DTSEDMIN=0.0
      NHDS=0
!----------------------------------------------------------------------C 
!     Sources indices. 
!----------------------------------------------------------------------C  
      LQSOUR = .false.
      SDSOUR = .false.
      CLSOUR = .false.
      CPSOUR = .false.
      LQSRND = .false.
      SDSRND = .false.
      CLSRND = .false.
      CPSRND = .false.
!----------------------------------------------------------------------C 
!     Other indices and variables. 
!----------------------------------------------------------------------C 
      NFRSD = 0 
!----------------------------------------------------------------------C 
!     Solution defaults 
!----------------------------------------------------------------------C 
      DO 280 N = 1,10 
        ISOLVE(N) = 0 
  280 CONTINUE 
      DO 290 N = 1,10 
        IDMEAN(N) = 1 
  290 CONTINUE 
      NRSD = 1 
      IRSD = 2 
      DO 300 N=1,4 
        RSDMX(N)  = 0.0001 
        RSD(N)    = 0.001 
  300 CONTINUE 
!----------------------------------------------------------------------C 
!     Reference node for diagnostics 
!----------------------------------------------------------------------C 
      IHREF   = 25 
      IFQREF  = IBIG 
      NMREF   = 0 
!----------------------------------------------------------------------C 
!     Iteration variables. 
!----------------------------------------------------------------------C 
      NSTEP  = 0 
      MNSTEP = 10 
!----------------------------------------------------------------------C 
!     Default parameters. 
!----------------------------------------------------------------------C  
      ALPHA = 1.0D+0 
      BETA  = 1.D+0  
      PHI   = 0.5D+0  
      THETA = 0.55D+0  
      iIFRIC = 3  
      rFRIC  = 0.2D-1  
      SIGMA = 1.d0
!----------------------------------------------------------------------C 
!     MXSTEP:  maximum number of time steps or maximum of iterations 
!----------------------------------------------------------------------C 
      MXSTEP = 32000. *300
      SCREEN = .FALSE. 
!----------------------------------------------------------------------C 
!     Time-history plots and tables 
!----------------------------------------------------------------------C 
      DTP       = BIG 
      DTC       = BIG 
      DTT       = BIG 
      DTPG      = BIG 
      NAMEQN(1) = 'WATER CONTINUITY' 
      NAMEQN(2) = 'SEDIMENT TRANSP.' 
      NAMEQN(3) = 'SOLUTE SPECIES  ' 
      NAMEQN(4) = 'PARTIC. SPECIES ' 
!----------------------------------------------------------------------C 
!     Print times. 
!----------------------------------------------------------------------C 
      NPRTM = 0 
!      DO 390 N = 1,40 
!        PRTM(N) = 0.D+0 
!  390 CONTINUE 
!----------------------------------------------------------------------!
!     Variable names.
!----------------------------------------------------------------------!
      call InitVariables
!----------------------------------------------------------------------!
!     Temporary arrays Sizes.
!----------------------------------------------------------------------!
      NRWORK=0
      NIWORK=0
!----------------------------------------------------------------------C 
!     End of INITAL group. 
!----------------------------------------------------------------------C 
! 
      RETURN 
      END 
 
!======================================================================C 

!======================================================================C 
! 
      SUBROUTINE RDINIT 
! 
!----------------------------------------------------------------------C 
! 
!     RDINIT: Read INITial conditions input group. 
! 
!----------------------------------------------------------------------C  
      USE CONSTS
      USE DRLINK 
      USE DRSPFL
      USE DRWFLD
      USE FILINX
      USE FLDINX
      USE HEADNG
      USE LOGCLS 
      USE LANDSF
      USE NUMBRS
      USE POINTS
      USE SOLVAR
      USE intersub
!----------------------------------------------------------------------C 
!     Implicit Double Precision. 
!----------------------------------------------------------------------C 
      IMPLICIT REAL*8 (A-H,O-Z) 
      IMPLICIT INTEGER*4 (I-N) 
!----------------------------------------------------------------------C 
!     Type Declarations. 
!----------------------------------------------------------------------C 
      CHARACTER*240 CHDUM 
      CHARACTER*80 ADUM,UNITS(4) 
      REAL*8 VALUE(4) 
      INTEGER*4 IRANGE(4) 
!----------------------------------------------------------------------C 
!     Write header to output file. 
!----------------------------------------------------------------------C 
      WRITE (IWR,'(/A)') ' Initial Conditions' 
      WRITE (IWR,'(A/)') ' ------------------' 
!----------------------------------------------------------------------C 
!     Restart file will be read for initial conditions. 
!----------------------------------------------------------------------C 
      IF( RESTART ) THEN 
        CALL RDRST 
        WRITE (IWR,'(A)') ' Restart Conditions ' 
        WRITE (IWR,'(A,1PE11.4)') '  Restart Time,       s: ',TIME 
        WRITE (IWR,'(A,1PE11.4)') '  Restart Time Step,  s: ',DT 
        GOTO 900
      ENDIF 
!----------------------------------------------------------------------C 
!     Found 'Initial Conditions' data. 
!     Read logic for computing the initial liquid pressure. 
!----------------------------------------------------------------------C 
!----------------------------------------------------------------------C 
!     Read the number of lines of initial conditions. 
!----------------------------------------------------------------------C 
      READ (IRD,*) NLIN 
      DO 300 NL = 1,NLIN 
        READ (IRD,'(A)') CHDUM 
        CALL LCASE( CHDUM ) 
        ISTART = 1 
        CALL RDCHR(ISTART,ICOMMA,CHDUM,ADUM) 
!----------------------------------------------------------------------C 
!     Read 2 real numbers (and their units) and a range for the data. 
!----------------------------------------------------------------------C 
        DO 310 I = 1,2 
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VALUE(I)) 
          CALL RDCHR(ISTART,ICOMMA,CHDUM,UNITS(I)) 
          IF( I .EQ. 1 ) THEN 
            CONV = 1.D+0 
            CALL RDUNIT(UNITS(I),CONV) 
            VALUE(I) = VALUE(I)*CONV 
          ELSE 
            CALL RDUNIT(UNITS(I),VALUE(I)) 
            VALUE(I) = VALUE(I)*CONV 
          ENDIF 
  310   CONTINUE
!----------------------------------------------------------------------C 
!     Read 4 integer numbers. 
!----------------------------------------------------------------------C 
        DO 320 I = 1,4 
          CALL RDINT(ISTART,ICOMMA,CHDUM,IRANGE(I)) 
  320   CONTINUE 
!----------------------------------------------------------------------C 
!     Store Flow Discharge data. 
!     Write initial conditions to output file. 
!----------------------------------------------------------------------C 
        IF( ADUM(1:14) .EQ. 'flow discharge' ) THEN 
          IF( NL .NE. 1) WRITE (IWR,'(/)') 
          WRITE (IWR,'(A)') '   Flow Discharge, m^3/s' 
!
          WRITE (IWR,'(2(A,I6),A,1PE11.4)') '  QH(',IRANGE(1),',',     &
           IRANGE(3), '): ', VALUE(1) 
          WRITE (IWR,'(A,1PE11.4)') '   X-direction Gradient: ',       &
           VALUE(2) 
          WRITE (IWR,'(2(A,I4))') '   Link : L = ',IRANGE(1),' to ',   &
           IRANGE(2) 
          WRITE (IWR,'(2(A,I4))') '   Point: I = ',IRANGE(3),' to ',   &
           IRANGE(4) 
! 
          call AssignGradientVariable(QH,Value(1),Value(2),iRange)
!----------------------------------------------------------------------C 
!     Store Water Surface data. 
!     Write initial conditions to output file. 
!----------------------------------------------------------------------C 
        ELSEIF( ADUM(1:13) .EQ. 'thalweg eleva' ) THEN 
          IF( NL .NE. 1) WRITE (IWR,'(/)') 
          WRITE (IWR,'(A)') '   Thalweg Elevation, m' 
          WRITE (IWR,'(2(A,I6),A,1PE11.4)') '  TSF(',IRANGE(1),',',   &
            IRANGE(3), '): ', VALUE(1) 
          WRITE (IWR,'(A,1PE11.4)') '   X-direction Gradient: ',      &
            VALUE(2) 
          WRITE (IWR,'(2(A,I4))') '   Link : L = ',IRANGE(1),' to ',  &
            IRANGE(2) 
          WRITE (IWR,'(2(A,I4))') '   Point: I = ',IRANGE(3),' to ',  &
            IRANGE(4) 
! 
          call AssignGradientVariable(TSF,Value(1),Value(2),iRange)
!----------------------------------------------------------------------C 
!     Store Water Surface data. 
!     Write initial conditions to output file. 
!----------------------------------------------------------------------C 
        ELSEIF( ADUM(1:13) .EQ. 'water surface' ) THEN 
          IF( NL .NE. 1) WRITE (IWR,'(/)') 
          WRITE (IWR,'(A)') '   Water Surface Elevation, m' 
          WRITE (IWR,'(2(A,I6),A,1PE11.4)') '  HR(',IRANGE(1),',',     &
           IRANGE(3), '): ', VALUE(1) 
          WRITE (IWR,'(A,1PE11.4)') '   X-direction Gradient: ',       &
           VALUE(2) 
          WRITE (IWR,'(2(A,I4))') '   Link : L = ',IRANGE(1),' to ',   &
           IRANGE(2) 
          WRITE (IWR,'(2(A,I4))') '   Point: I = ',IRANGE(3),' to ',   &
           IRANGE(4) 
! 
          call AssignGradientVariable(HR,Value(1),Value(2),iRange) 
          call AddVariableRegion(HR,HR,-TSF,iRange)
          HR = MAX(0.D+0,HR)  
!----------------------------------------------------------------------C 
!     Store Water Depth data. 
!     Write initial conditions to output file. 
!----------------------------------------------------------------------C 
        ELSEIF( ADUM(1:14) .EQ. 'water depth' ) THEN 
          IF( NL .NE. 1) WRITE (IWR,'(/)') 
          WRITE (IWR,'(A)') '   Water Depth, m' 
          WRITE (IWR,'(2(A,I6),A,1PE11.4)') '  HR(',IRANGE(1),',',    &
            IRANGE(3), '): ', VALUE(1) 
          WRITE (IWR,'(A,1PE11.4)') '   X-direction Gradient: ',      &
            VALUE(2) 
          WRITE (IWR,'(2(A,I4))') '   Link : L = ',IRANGE(1),' to ',  &
            IRANGE(2) 
          WRITE (IWR,'(2(A,I4))') '   Point: I = ',IRANGE(3),' to ',  &
            IRANGE(4) 
! 
          call AssignGradientVariable(HR,Value(1),Value(2),iRange) 
!----------------------------------------------------------------------C 
!     Store Sediment Concentration data. 
!     Write initial conditions to output file. 
!----------------------------------------------------------------------C 
        ELSEIF( ADUM(1:13) .EQ. 'sediment conc' ) THEN   
          IF( NL .NE. 1 ) WRITE (IWR,'(/)') 
          WRITE (IWR,'(A)') '   Sediment Concentration, 1/m^3' 
          WRITE (IWR,'(2(A,I6),A,1PE11.4)') '   SH(',IRANGE(1),',',    &
            IRANGE(3), '): ', VALUE(1) 
          WRITE (IWR,'(A,1PE11.4)') '   X-direction Gradient: ',       &
            VALUE(2) 
          WRITE (IWR,'(2(A,I4))') '   Link : L = ',IRANGE(1),' to ',   & 
            IRANGE(2) 
          WRITE (IWR,'(2(A,I4))') '   Point: I = ',IRANGE(3),' to ',   &
            IRANGE(4) 
 
          call AssignGradientVariable(SH,Value(1),Value(2),iRange) 
!----------------------------------------------------------------------C 
!     Store Species Concentration data. 
!     Write initial conditions to output file. 
!----------------------------------------------------------------------C 
        ELSEIF( ADUM(1:12) .EQ. 'species conc' ) THEN   
          IF( NL .NE. 1 ) WRITE (IWR,'(/)') 
          WRITE (IWR,'(A)') '   Species Concentration, 1/m^3' 
          WRITE (IWR,'(2(A,I6),A,1PE11.4)') '   CH(',IRANGE(1),',',   &
            IRANGE(3), '): ', VALUE(1) 
          WRITE (IWR,'(A,1PE11.4)') '   X-direction Gradient: ',      &
            VALUE(2) 
          WRITE (IWR,'(2(A,I4))') '   Link : L = ',IRANGE(1),' to ',  &
            IRANGE(2) 
          WRITE (IWR,'(2(A,I4))') '   Point: I = ',IRANGE(3),' to ',  &
            IRANGE(4) 
! 
          call AssignGradientVariable(CH,Value(1),Value(2),iRange) 
!----------------------------------------------------------------------C 
!     Store Particulate Concentration data. 
!     Write initial conditions to output file. 
!----------------------------------------------------------------------C 
        ELSEIF( ADUM(1:16) .EQ. 'particulate conc' ) THEN   
          IF( NL .NE. 1 ) WRITE (IWR,'(/)') 
          WRITE (IWR,'(A)') '   Particulate Concentration, 1/m^3' 
          WRITE (IWR,'(2(A,I6),A,1PE11.4)') '   CSH(',IRANGE(1),',',  &
            IRANGE(3), '): ', VALUE(1) 
          WRITE (IWR,'(A,1PE11.4)') '   X-direction Gradient: ',      &
            VALUE(2) 
          WRITE (IWR,'(2(A,I4))') '   Link : L = ',IRANGE(1),' to ',  &
            IRANGE(2) 
          WRITE (IWR,'(2(A,I4))') '   Point: I = ',IRANGE(3),' to ',  &
            IRANGE(4) 
! 
          call AssignGradientVariable(CSH,Value(1),Value(2),iRange) 
!----------------------------------------------------------------------C 
!     Store Bottom Species Concentration data. 
!     Write initial conditions to output file. 
!----------------------------------------------------------------------C 
        ELSEIF( ADUM(1:14) .EQ. 'bottom species' ) THEN   
          IF( NL .NE. 1 ) WRITE (IWR,'(/)') 
          WRITE (IWR,'(A)') '   Bottom Species Concentration, 1/m^3' 
          WRITE (IWR,'(2(A,I6),A,1PE11.4)') '   CB(',IRANGE(1),',',    &
            IRANGE(3), '): ', VALUE(1) 
          WRITE (IWR,'(A,1PE11.4)') '   X-direction Gradient: ',       &
            VALUE(2) 
          WRITE (IWR,'(2(A,I4))') '   Link : L = ',IRANGE(1),' to ',   &
            IRANGE(2) 
          WRITE (IWR,'(2(A,I4))') '   Point: I = ',IRANGE(3),' to ',   &
            IRANGE(4) 
! 
          call AssignGradientVariable(CB,Value(1),Value(2),iRange) 
        ENDIF 
  300 CONTINUE 
  900 CONTINUE 
!----------------------------------------------------------------------C 
!     End of RDINIT group. 
!----------------------------------------------------------------------C 
! 
      RETURN 
      END 
 
!======================================================================C 
!======================================================================C 
!
   SUBROUTINE ReallocateR8 (N0,N1,V)
!    
!----------------------------------------------------------------------C 
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
!----------------------------------------------------------------------C 
      INTEGER*4 N0,N1
      REAL*8, POINTER, intent(inout) :: V(:) 
!----------------------------------------------------------------------C 
      REAL*8, POINTER :: Vt(:) !--New array
!----------------------------------------------------------------------C 
        ALLOCATE(Vt(N1)); Vt=0.d0;
        if (ASSOCIATED(V)) then  
          N2 = N0  
          if( N0 == 0 ) N2 = 1   
          N2 = min(N1,N2)
            IF (N2 > 0) Vt(1:N2) = V(1:N2)
          DEALLOCATE(V)
        endif
        V=>Vt
!----------------------------------------------------------------------C 
!     End of ReallocateR8 group. 
!----------------------------------------------------------------------C 
! 
      RETURN
    END SUBROUTINE
!======================================================================C 
!
    SUBROUTINE ReallocateI4(N0,N1,V)
!
!----------------------------------------------------------------------C 
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
!----------------------------------------------------------------------C 
      INTEGER*4 N0,N1
      INTEGER*4, POINTER :: V(:)
!----------------------------------------------------------------------C 
      INTEGER*4, POINTER :: Vt(:) !--New array
!----------------------------------------------------------------------C 
      ALLOCATE(Vt(N1)); Vt=0;
      if (ASSOCIATED(V)) then 
        N2 = min(N0,N1)  
        Vt(1:N2) = V(1:N2)
        DEALLOCATE(V)
      endif
      V=>Vt
!----------------------------------------------------------------------C 
!     End of ReallocateI4 group. 
!----------------------------------------------------------------------C 
! 
      RETURN
    END SUBROUTINE
!======================================================================C 
!
    SUBROUTINE Reallocate2DI4(NX0,NY0,NX1,NY1,V)
!
!----------------------------------------------------------------------C 
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
!----------------------------------------------------------------------C 
      INTEGER*4 NX0,NY0,NX1,NY1
      INTEGER*4, POINTER :: V(:,:)
!----------------------------------------------------------------------C 
      INTEGER*4, POINTER :: Vt(:,:) !--New array
!----------------------------------------------------------------------C 
      ALLOCATE(Vt(NX1,NY1)); Vt=0;
      if (ASSOCIATED(V)) then 
        NX2 = MIN(NX0,NX1); 
        NY2 = MIN(NY0,NY1);
        Vt(1:NX2,1:NY2) = V(1:NX2,1:NY2)
        DEALLOCATE(V)
      endif
      V=>Vt
!----------------------------------------------------------------------C 
!     End of Reallocate2DI4 group. 
!----------------------------------------------------------------------C 
! 
      RETURN
    END SUBROUTINE
!======================================================================C 
!
    SUBROUTINE Reallocate2DR8(NX0,NY0,NX1,NY1,V)
!
!----------------------------------------------------------------------C 
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
!----------------------------------------------------------------------C 
      INTEGER*4 NX0,NY0,NX1,NY1
      REAL*8, POINTER :: V(:,:)
!----------------------------------------------------------------------C 
      REAL*8, POINTER :: Vt(:,:) !--New array
!----------------------------------------------------------------------C 
      ALLOCATE(Vt(NX1,NY1)); Vt=0.d0;
      if (ASSOCIATED(V)) then 
        NX2 = MIN(NX0,NX1); 
        NY2 = MIN(NY0,NY1);
        Vt(1:NX2,1:NY2) = V(1:NX2,1:NY2)
        DEALLOCATE(V)
      endif
      V=>Vt
!----------------------------------------------------------------------C 
!     End of Reallocate2DR8 group. 
!----------------------------------------------------------------------C 
! 
      RETURN
    END SUBROUTINE
!======================================================================C 
!
    SUBROUTINE Reallocate3DI4(NX0,NY0,NZ0,NX1,NY1,NZ1,V)
!
!----------------------------------------------------------------------C 
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
!----------------------------------------------------------------------C 
      INTEGER*4 NX0,NY0,NZ0,NX1,NY1,NZ1
      INTEGER*4, POINTER :: V(:,:,:)
!----------------------------------------------------------------------C 
      INTEGER*4, POINTER :: Vt(:,:,:) !--New array
!----------------------------------------------------------------------C 
      ALLOCATE(Vt(NX1,NY1,NZ1)); Vt=0;
      if (ASSOCIATED(V)) then 
        NX2 = MIN(NX0,NX1); 
        NY2 = MIN(NY0,NY1);
        NZ2 = MIN(NZ0,NZ1);
        Vt(1:NX2,1:NY2,1:NZ2) = V(1:NX2,1:NY2,1:NZ2)
        DEALLOCATE(V)
      endif
      V=>Vt
!----------------------------------------------------------------------C 
!     End of Reallocate3DI4 group. 
!----------------------------------------------------------------------C 
! 
      RETURN
    END SUBROUTINE
!======================================================================C 
!======================================================================C 
!======================================================================C 

