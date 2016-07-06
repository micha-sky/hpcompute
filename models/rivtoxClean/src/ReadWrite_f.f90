!==============================================================================C 
! 
      SUBROUTINE RDSECT 
! 
!------------------------------------------------------------------------------C 
! 
!     RDSECT: ReaD and process cross-SECTional data input group. 
! 
!----------------------------------------------------------------------C  
      USE FILINX 
      USE FLDINX 
      USE DRLINK 
      USE DRNAME 
      USE CRSECT 
      USE NUMBRS
!------------------------------------------------------------------------------C 
!     Implicit Double Precision. 
!------------------------------------------------------------------------------C 
      IMPLICIT REAL*8 (A-H,O-Z) 
      IMPLICIT INTEGER*4 (I-N) 
!------------------------------------------------------------------------------C 
!     Type Declarations. 
!------------------------------------------------------------------------------C 
      CHARACTER*240 CHDUM 
      CHARACTER*120 NOTES 
      CHARACTER ATYPE*80,PEROPT*80,UNITS*80, fname*80 
      INTEGER*4 NPZN(4) 

      REAL*8,DIMENSION(:),POINTER:: xTab,yTab,hTab
      REAL*8, POINTER :: rTmp(:) 
      CHARACTER*10, POINTER :: tmpSEC(:)  
      INTEGER*4, POINTER :: tmpNSCT(:,:)
!------------------------------------------------------------------------------C 
!     WRITE header to output file. 
!------------------------------------------------------------------------------C 
      WRITE (IWR,'(/A)') ' Cross-Sectional Data' 
      WRITE (IWR,'(A)')  ' --------------------' 
!------------------------------------------------------------------------------C 
!     Loop over cross-sectional data records. 
!     Read record containing number of sets of cross-sectional data.
!     Set starting location for tables to zero.
!------------------------------------------------------------------------------C 
      NSEC = 0  
      NCOUNT = 0 
      READ (IRD,*) NSETT
!----------------------------------------------------------------------C  
      ALLOCATE(TSEC(NSETT),NSCT(2,NSETT))  
      TSEC = ' '; NSCT = 0; IRDf = IRD; nSize = NSETT;
!
!      DO 600 N = 1,NSETT   
      N = 1 
10    continue      
      if( N > NSETT ) goto 600
!------------------------------------------------------------------------------C 
!     Read record into character buffer and convert to lower case.
!------------------------------------------------------------------------------C 
        READ (IRDf,'(A)') CHDUM 
        CALL LCASE( CHDUM )
!------------------------------------------------------------------------------C 
!     Read table name of cross-section.
!------------------------------------------------------------------------------C 
        ISTART = 1
        CALL RDCHR(ISTART,ICOMMA,CHDUM,ATYPE)  
        
        if( ATYPE(1:8) == 'fromfile' ) then   
           Nold = N; NSETTold = NSETT;  
           CALL RDCHR(ISTART,ICOMMA,CHDUM,ATYPE)  
           fname = trim(ATYPE);
           write(IWR,*)  
           write(IWR,*) 'Read cross-sections from file ',trim(fname); 
           write(IWR,*) '------------------------------',('-',i=1,len_trim(fname))
           
           open(270,file=trim(ATYPE));  
           IRDf = 270; 
           READ (IRDf,*) NSETT  
           if( associated( Tsec ) ) then  
             nSize1 = nSize +NSETT -1;  
             ALLOCATE(tmpSEC(nSize1),tmpNSCT(2,nSize1));   
             tmpSEC = ' '; tmpNSCT = 0; 
             tmpSEC(1:nSize) = TSEC; tmpNSCT(1:2,1:nSize) = NSCT; 
             deallocate( TSEC,NSCT ); nSize = nSize1; 
             TSEC => tmpSEC; NSCT => tmpNSCT;
           endif  
           N = 1;
           READ (IRDf,'(A)') CHDUM 
           CALL LCASE( CHDUM )
           ISTART = 1
           CALL RDCHR(ISTART,ICOMMA,CHDUM,ATYPE)  
        endif
            
        DO 100 I = 1,NSEC
          IF( ATYPE .EQ. TSEC(I) ) THEN
            WRITE (*,'(2A)') ATYPE,                                   &
                   '- Cross-section with such name already exists!'
            WRITE (IWR,'(2A)') ATYPE,                                 &
                   '- Cross-section with such name already exists!'
            WRITE (*,'(1A)') ' Run aborting...'
            STOP
          ENDIF
  100   CONTINUE
        NSEC = NSEC +1   
        TSEC(NSEC) = ATYPE 
!------------------------------------------------------------------------------C 
!     Read which type of Level-Width or Distance-Level pairs will be input.
!------------------------------------------------------------------------------C 
        CALL RDCHR(ISTART,ICOMMA,CHDUM,ATYPE) 
        iSec = 1 
        IF( ATYPE(1:5) .EQ. 'level' ) iSec = 2 
!------------------------------------------------------------------------------C 
        CALL RDCHR(ISTART,ICOMMA,CHDUM,PEROPT)
!------------------------------------------------------------------------------C 
!     Tabular.
!------------------------------------------------------------------------------C 
        IF( PEROPT(1:4) .EQ. 'tabu' ) THEN
!------------------------------------------------------------------------------C 
!     Read tabular data for cross-section.
!------------------------------------------------------------------------------C 
          READ (IRDf,*) NLIN    
          ALLOCATE( xTab(1:NLIN),yTab(1:NLIN) ); 
          xTab = 0.d0; yTab = 0.d0; iRec = 0; 
!
          DO 300 I = 1,NLIN
            READ (IRDf,'(A)') CHDUM
            CALL LCASE( CHDUM )
            ISTART = 1
!
            CALL RDDPR(ISTART,ICOMMA,CHDUM,RLX)
            CALL RDCHR(ISTART,ICOMMA,CHDUM,UNITS)
            CALL RDUNIT(UNITS,RLX)
!
            CALL RDDPR(ISTART,ICOMMA,CHDUM,RLY)
            CALL RDCHR(ISTART,ICOMMA,CHDUM,UNITS)
            CALL RDUNIT(UNITS,RLY)
!------------------------------------------------------------------------------C 
!     Sort table in ascending water depth order.
!------------------------------------------------------------------------------C  
            iRec = iRec +1  
            DO 130 J = 1, iRec -1
              IF( RLX .LE.  xTab(J) ) THEN
                DO 120 K = iRec,J+1,-1
                  xTab(K) = xTab(K-1)
                  yTab(K) = yTab(K-1)
  120           CONTINUE
                xTab(J) = RLX
                yTab(J) = RLY
                GOTO 140
              ENDIF
  130       CONTINUE
            xTab(iRec) = RLX
            yTab(iRec) = RLY
  140       CONTINUE
  300     CONTINUE
!------------------------------------------------------------------------------C 
!     Write tabular data for cross-section.
!------------------------------------------------------------------------------C 
          WRITE (IWR,'(/2A)') ' Table name: ', TSEC(NSEC) 
          WRITE (IWR,'(2A,I6,A)') '   Cross-sectional Function: ',       &
             'Tabular with ',NLIN,' entries'
          IF( iSec == 1 ) THEN 
            WRITE (IWR,'(A)') '   Level-Width pairs: ' 
          ELSE  
            WRITE (IWR,'(A)') '   Distance-Level pairs: ' 
          ENDIF   
!
          ICNT   = 0
          DO 310 II = 1,NLIN
            ICNT = ICNT +1
            IF( iSec == 1 ) THEN 
              WRITE (IWR,'(A,I4,2(A,1PE11.4))') '   Entry No.',ICNT,     &
                 '   Depth: ', xTab(II),                                 &
                 '   Width: ', yTab(II)
            ELSE  
              WRITE (IWR,'(A,I4,2(A,1PE11.4))') '   Entry No.',ICNT,     &
                 '   X-coord.: ', xTab(II),                              &
                 '   Level   : ', yTab(II)
            ENDIF   
  310     CONTINUE
!------------------------------------------------------------------------------C 
!     Compute cross-section area, width and wetted perimeter.
!------------------------------------------------------------------------------C 
          if( nCount == 0 ) then  
            ALLOCATE(SWD(NLIN),SHG(NLIN),SAR(NLIN),SPer(NLIN))  
            SWD = 0.d0; SHG = 0.d0; SAR = 0.d0; SPer = 0.d0;
          else  
            ALLOCATE(rTmp(nCount)); 
            
            rTmp(1:nCount) = SWD(1:nCount); 
            DEALLOCATE(SWD); ALLOCATE(SWD(nCount+NLIN)); 
            SWD(1:nCount) = rTmp; 
            
            rTmp(1:nCount) = SAR(1:nCount);
            DEALLOCATE(SAR); ALLOCATE(SAR(nCount+NLIN)); 
            SAR(1:nCount) = rTmp; 
            
            rTmp(1:nCount) = SPer(1:nCount);
            DEALLOCATE(SPer); ALLOCATE(SPer(nCount+NLIN)); 
            SPer(1:nCount) = rTmp; 
            
            rTmp(1:nCount) = SHG(1:nCount);
            DEALLOCATE(SHG); ALLOCATE(SHG(nCount+NLIN)); 
            SHG(1:nCount) = rTmp; 
            
            DEALLOCATE(rTmp);
          endif    

          IF( iSec == 1 ) THEN 
!     Level-Width pairs. 
            NCOUNT = NCOUNT +1
            NSCT(1,NSEC) = NCOUNT
            DO i=1,NLIN 
              SHG(NCOUNT) = xTab(i) 
              SWD(NCOUNT) = yTab(i) 
              If( i == 1) Then 
                Sar(NCOUNT)  = 0.d0; 
                SPer(NCOUNT) = yTab(i);   
              Else  
                Sar(NCOUNT) = Sar(NCOUNT-1)                          &
                   +0.5d0*(yTab(i) +yTab(i-1))*(xTab(i) -xTab(i-1))   
                SPer(NCOUNT) = SPer(NCOUNT-1)                        &
                   +2.d0*SQRT((xTab(i) -xTab(i-1))**2                &
                   +0.25d0*(yTab(i) -yTab(i-1))**2)    
              EndIf 
              NCOUNT = NCOUNT +1       
            ENDDO  
            nCount = nCount -1
          ELSE    ! iSec == 2
!     Distance-Level pairs.            
            ALLOCATE( hTab(1:NLIN) ); hTab = yTab; 
            Do i=1,NLIN-1 
              Do j=i+1,NLIN    
                If( hTab(j) < hTab(i) ) Then 
                  wr = hTab(i); hTab(i) = hTab(j); hTab(j) = wr; 
                EndIf  
              EndDo      
            EndDo  
                   
            NCOUNT = NCOUNT +1
            NSCT(1,NSEC) = NCOUNT 
            Do i=1,NLIN-1  
              SHG(NCOUNT) = hTab(i) -hTab(1)
              SWD(NCOUNT) = 0.d0; 
              If( i == 1 ) Then 
                Sar(NCOUNT) = 0.d0; SPer(NCOUNT) = 0.d0;  
              Else 
                Sar(NCOUNT)  = Sar(NCOUNT-1)   
                SPer(NCOUNT) = SPer(NCOUNT-1)  
              EndIf    

              Ist = 0; xb = xTab(NLIN); xe = xTab(1); 
              Do j=2,NLIN 
                If( hTab(i) <= yTab(j-1) .and. hTab(i) >= yTab(j) ) Then  ! begin of subcross-section  
                  if( hTab(i) == yTab(j) .and. i /= 1 ) cycle
                  tg_b = (yTab(j) -yTab(j-1))/(xTab(j) -xTab(j-1)) 
                  If( abs(tg_b) <= small ) Then  
                    xb = xTab(j-1)   
                  Else  
                    xb = xTab(j-1) +(hTab(i) -yTab(j-1))/tg_b  
                  EndIf  
                  Ist = j;
!                  If( tg_a < 0.d0 .and. Ist == 0 ) then 
!                    xb = xs; Ist = j; tg_s = tg_a;  
!                  ElseIf( tg_a >= 0.d0 .and. Ist /= 0 ) Then 
!                    xe = xs; 
!                    If( tg_a > Small ) Ist = 0;   
!                  Endif               
                ElseIf( hTab(i) >= yTab(j-1) .and. hTab(i) <= yTab(j) ) Then  ! End of subcross-section
                  if( hTab(i) == yTab(j-1) .and. i /= 1 ) cycle
                  
                  If( abs( xTab(j)-xTab(j-1) ) <= small ) Then
                     xe = xTab(j-1)
                  else
                    tg_a = (yTab(j) -yTab(j-1))/(xTab(j) -xTab(j-1)) 
                    If( abs(tg_a) <= small ) Then  
                      xe = xTab(j-1)   
                    Else  
                      xe = xTab(j-1) +(hTab(i) -yTab(j-1))/tg_a  
                    EndIf 
                  Endif
!     Compute cross-section area, width and wetted perimeter. 
                  SWD(NCOUNT) = SWD(NCOUNT) +xe -xb  
                  If( tg_a <= small ) Then  
                    SPer(NCOUNT) = SPer(NCOUNT) +xe -xb 
                  ElseIf( i > 1 ) Then   
!                    xb1 = xTab(Ist-1) +(hTab(i-1) -yTab(Ist-1))/tg_b   
!                    xe1 = xTab(j-1) +(hTab(i-1) -yTab(j-1))/tg_a  

                    dHr = hTab(i) -hTab(i-1);
                    Sar(NCOUNT)  = Sar(NCOUNT)                          &
                         +dHr*(xe -xb -0.5d0*dHr*(-1.d0/tg_b +1.d0/tg_a))
                    SPer(NCOUNT) = SPer(NCOUNT)                         &
                         +dHr*(sqrt(1.d0 +1.d0/tg_b**2) +sqrt(1.d0 +1.d0/tg_a**2))   
                  EndIf                        
                EndIf   
              EndDo  
              NCOUNT = NCOUNT +1       
            EndDo  
            nCount = nCount -1
            DeAllocate( hTab );
          ENDIF
!          
          NSCT(2,NSEC) = NCOUNT  
          DeAllocate( xTab,yTab ); 
!------------------------------------------------------------------------------C 
        ENDIF
!------------------------------------------------------------------------------C    
        if( IRDf == 270 .and. N == NSETT ) then    
!       Reading from file is complited  
          write(IWR,*)  
          write(IWR,*) 'Reading cross-sections from file ',trim(fname),' is finished'  
          write(IWR,*) '---------------------------------',('-',i=1,len_trim(fname)), &
                       '------------'
          IRDf = IRD; N = Nold; NSETT = NSETTold;  
        endif    
        N = N +1  
        goto 10
  600 CONTINUE
!------------------------------------------------------------------------------C 
!     Format Statements. 
!------------------------------------------------------------------------------C 
 9007 FORMAT(/' Input or Compilation Error: The number of table ',     &
              'entries specified [',I6,'] for input of the ',          &
              'cross-sectional name ', A,                              &
              ' exceeds the maximum set in parameter LSCR [',I6,'].')
 9008 FORMAT(/' Input or Compilation Error: The number of table',      &
              ' entries for cross-sectional data exceeds',/,           &
              ' the maximum set in',                                   &
              ' parameter LSCT [',I6,'].') 
!------------------------------------------------------------------------------C 
!     End of RDSECT group. 
!------------------------------------------------------------------------------C 
! 
      RETURN 
      END 
 
!==============================================================================C 
