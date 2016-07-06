C======================================================================C
C
      SUBROUTINE RDBCND
C
C----------------------------------------------------------------------C
C
C     RDBCND: ReaD Boundary Conditions for NoDe input group.  
C
C----------------------------------------------------------------------C  
      USE FLDINX
      USE DRLINK
      USE DRNAME  
      USE DRNODE
      USE DRBCLQ 
      USE FILINX
C----------------------------------------------------------------------C
C     Implicit Double Precision.
C----------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C----------------------------------------------------------------------C
      INCLUDE 'utils.inc' ! subroutine interfaces
C----------------------------------------------------------------------C
C     Type Declarations.
C----------------------------------------------------------------------C
      CHARACTER*240 CHDUM
      CHARACTER*80 ADUM,UNITS,TYPE
      CHARACTER*20 TABLE  
      CHARACTER*20, POINTER :: TBNAME(:)
      REAL*8 DUMY(6),VALUE
C----------------------------------------------------------------------C 
C     For i-th boundary node (1 <= i <= NBCL):
C        IBCND(i)  - number of i-th node;
C        IBCTY(i)  - type of boundary condition;
C        MBCND(i)  - table number of boundary values; 
C     For j-th table of boundary values:
C        NBCT(1,j) - first line in the table;
C        NBCT(2,j) - last line in the table.
C     For k-th line (1 <= k <= NBTU) in table of boundary values:
C        BCV(k)  -  boundary value at
C        BCT(k)  -  time moment.
C----------------------------------------------------------------------C
C     Write header to output file.
C----------------------------------------------------------------------C
      WRITE (IWR,'(/A)') ' Node Water Boundary Conditions'
      WRITE (IWR,'(A )') ' ------------------------------'
C----------------------------------------------------------------------C
C     Found 'Water Boundaries' data.
C     Read record containing the number of sets of boundary
C     condition data to read.
C----------------------------------------------------------------------C
      READ (IRD,*) NSETT  
C----------------------------------------------------------------------C 
      nk = NSETT +1 
      ALLOCATE(IBCND(nBound),IBCTY(nBound),MBCND(nBound),NBCT(2,nk))  
      ALLOCATE(TBNAME(nk),BCV(nk),BCT(nk)); TBNAME = ' ';
C----------------------------------------------------------------------C
C     Initialize boundary condition variables.
C----------------------------------------------------------------------C
      NBCL = 0 
      do n=1, nBound  
        IBCND(N)  = 1
        IBCTY(N)  = 2
        MBCND(N)  = 0
      enddo    
      DO 100 N = 1,nk
        NBCT(1,N) = 0
        NBCT(2,N) = 0  
        BCT(N) = 0.D+0
        BCV(N) = 0.D+0
  100 CONTINUE   
      NBTU = 0
      ITABLE = 0
      ILINE = 0
C----------------------------------------------------------------------C
      DO 600 NS = 1,NSETT
C----------------------------------------------------------------------C
C     Read record into character buffer and convert to lower case.
C----------------------------------------------------------------------C
        READ (IRD,'(A)') CHDUM
        CALL LCASE( CHDUM )
C----------------------------------------------------------------------C
C     Read 'boundary type'
C----------------------------------------------------------------------C 
        ISTART = 1
        CALL RDCHR(ISTART,ICOMMA,CHDUM,ADUM)
        IF( ADUM(1:20) .EQ. 'discharge @ boundary' ) THEN
          IBTYPE = 1
        ELSEIF( ADUM(1:18) .EQ. 'surface @ boundary' ) THEN
          IBTYPE = 2
        ENDIF
C----------------------------------------------------------------------C
C     Read boundary temporal variation.
C----------------------------------------------------------------------C
        CALL RDCHR(ISTART,ICOMMA,CHDUM,TYPE)
        IF( TYPE(1:8) .EQ. 'constant' ) THEN
C----------------------------------------------------------------------C
C     Read boundary values and units.
C----------------------------------------------------------------------C 
          CALL RDDPR(ISTART,ICOMMA,CHDUM,VALUE)
          CALL RDCHR(ISTART,ICOMMA,CHDUM,UNITS)
          CALL RDUNIT(UNITS,VALUE)

          ITEMP = 1 
          NBTU = NBTU +1
          ITAB = NBTU
          TBNAME(NBTU) = 'constant'
          ILINE = ILINE +1 
          NBCT(1,NBTU) = ILINE
          NBCT(2,NBTU) = ILINE
          BCV(ILINE) = VALUE 
        ELSE
C----------------------------------------------------------------------C
C     Read table name of boundary values.
C----------------------------------------------------------------------C
          CALL RDCHR(ISTART,ICOMMA,CHDUM,ADUM)
          TABLE = ADUM
          CALL RDCHR(ISTART,ICOMMA,CHDUM,ADUM)
          ITEMP = 2
          DO 200 N = 1,NBTU
            IF( TBNAME(N) .EQ. TABLE ) THEN
              ITAB = N
              GOTO 250
            ENDIF  
  200     CONTINUE
          ITABLE = ITABLE +1
          NBTU = NBTU +1
          TBNAME(NBTU) = TABLE
          ITAB = NBTU
  250     CONTINUE  
        ENDIF
C----------------------------------------------------------------------C
C     Read nodal range of boundary (2 integers).
C----------------------------------------------------------------------C
        CALL RDINT(ISTART,ICOMMA,CHDUM,IBS)
        CALL RDINT(ISTART,ICOMMA,CHDUM,IBE)
C----------------------------------------------------------------------C
C     Write input boundary condition data to output file.
C----------------------------------------------------------------------C
        WRITE (IWR,'(2(A,I4))') 
     &             ' Boundary Domain (Nodes): I = ',IBS,' to ',IBE

        IF( IBTYPE .EQ. 1 ) THEN
          WRITE (IWR,'(A)') ' Boundary Condition: Discharge @ Boundary'
          IF( ITEMP .EQ. 1 ) THEN
            WRITE (IWR,'(A)')     ' Temporal Variation: Constant'
            WRITE (IWR,'(A,1PE11.4)') 
     &             ' Water Discharge, m^3/s: ',VALUE
          ELSE
            WRITE (IWR,'(A)')     ' Temporal Variation: Tabular'
            WRITE (IWR,'(2A)')    ' Boundary Table Name: ',TBNAME(NS)
          ENDIF 
        ELSEIF( IBTYPE .EQ. 2 ) THEN
          WRITE (IWR,'(A)') ' Boundary Condition: Surface @ Boundary'
          IF( ITEMP .EQ. 1 ) THEN
            WRITE (IWR,'(A)')     ' Temporal Variation: Constant'
            WRITE (IWR,'(A,1PE11.4)') 
     &             ' Water Surface Elevation,    m: ',VALUE
          ELSE
            WRITE (IWR,'(A)')     ' Temporal Variation: Tabular'
            WRITE (IWR,'(2A)')    ' Boundary Table Name: ',TBNAME(NS)
          ENDIF 
        ENDIF
C----------------------------------------------------------------------C
C     Assign values to boundary variables.
C----------------------------------------------------------------------C
        ISTOP = 0
        DO 300 I = IBS,IBE
C----------------------------------------------------------------------C
C     Check for boundary values applied to interior surfaces.
C----------------------------------------------------------------------C
          IERROR = 0
          IF( (IBTYPE .EQ. 2) .AND. (KK(I) > 1) ) THEN
             IERROR = 1
          ENDIF
          IF( IERROR .EQ. 1 ) THEN
            WRITE (IWR,9001) I,NodeName(I)
            ISTOP = 1
          ENDIF
          NBCL = NBCL +1
          IBCND(NBCL) = I 
          IBCTY(NBCL) = IBTYPE
          MBCND(NBCL) = ITAB
  300   CONTINUE
  600 CONTINUE
      IF( ISTOP .EQ. 1 ) THEN
        WRITE (*,'(/2A)') ' INPUT OR COMPILATION ERROR: ',
     &    'Water Boundary Conditions Input (check output file).' 
        WRITE (IWR,'(/2A)') ' INPUT OR COMPILATION ERROR: ',
     &    'Water Boundary Conditions Input.'
        STOP
      ENDIF  
C----------------------------------------------------------------------C
C     Read boundary tables.
C----------------------------------------------------------------------C
      IF( ITABLE .EQ. 0 ) GOTO 900 
C----------------------------------------------------------------------C
C     Length of table.
C----------------------------------------------------------------------C 
      READ (IRD,*) NLIN 
!----------------------------------------------------------------------!
!       REALLOCATE  BCVS,BCTS arrays.
!----------------------------------------------------------------------!
      nk = iLine +NLIN +1;  
      nkk = iLine
      call ReallocateR8 (nkk,nk,BCV)
      call ReallocateR8 (nkk,nk,BCT)
C----------------------------------------------------------------------C
C     Read a table ...
C     Read the abcissa value, the ordinate value and their units.
C----------------------------------------------------------------------C
        DO 700 I = 1,NLIN
          READ (IRD,'(A)') CHDUM
          CALL LCASE( CHDUM )
          ISTART = 1
          CALL RDCHR(ISTART,ICOMMA,CHDUM,ADUM)
          TABLE = ADUM
C----------------------------------------------------------------------C 
C     Read a table ...
C     Read the abcissa value and units.
C----------------------------------------------------------------------C
          CALL RDDPR(ISTART,ICOMMA,CHDUM,DUMY(1))
          CALL RDCHR(ISTART,ICOMMA,CHDUM,UNITS)
          CALL RDUNIT(UNITS,DUMY(1))
C----------------------------------------------------------------------C
C     Read a table ...
C     Read the ordinate values and units.
C----------------------------------------------------------------------C 
          DO J=2,2
            CALL RDDPR(ISTART,ICOMMA,CHDUM,DUMY(J))
            CALL RDCHR(ISTART,ICOMMA,CHDUM,UNITS)
            CALL RDUNIT(UNITS,DUMY(J))
          ENDDO
C----------------------------------------------------------------------C
C     Check for matching table name.
C----------------------------------------------------------------------C
          DO 710 N = 1, NBTU
            IF( TABLE .EQ. TBNAME(N)) THEN
              NTAB = N
              GOTO 720
            ENDIF
  710     CONTINUE
          WRITE (*,'(1A)') ' No matching table name!'
          WRITE (*,'(2A)') ' Type we''re looking for is: ',TABLE
          WRITE (*,'(1A)') ' Run aborting...'
          STOP
  720     CONTINUE
C----------------------------------------------------------------------C
C     Assign boundary variables.
C----------------------------------------------------------------------C
          ILINE = ILINE +1
          IF( NBCT(1,NTAB) .EQ. 0 ) NBCT(1,NTAB) = ILINE
          NBCT(2,NTAB) = ILINE
          BCT(ILINE) = DUMY(1)
          BCV(ILINE) = DUMY(2)
  700   CONTINUE
C----------------------------------------------------------------------C
C     Write boundary tables.
C----------------------------------------------------------------------C
        DO 800 N = 1,NBTU
          IF( NBCT(1,N) .EQ. NBCT(2,N) ) GOTO 800
          DO 820 I = 1,NBCL
            IF( MBCND(I) .EQ. N ) THEN
              IBTYPE = IBCTY(I)
              GOTO 830
            ENDIF  
  820     CONTINUE        
  830     CONTINUE        
          WRITE (IWR,'(/2A)') ' Boundary Table: ',TBNAME(N)
          NLIN = NBCT(2,N) -NBCT(1,N) +1
          DO 810 I = 1,NLIN
            J = NBCT(1,N) + I - 1
            IF( IBTYPE .EQ. 1 ) THEN
              WRITE (IWR,'(A,I4,3(A,1PE11.4))') ' Entry No.',I,
     &          '  Time, s: ',BCT(J),
     &          ' Water Discharge, m^3/s: ',BCV(J)
            ELSEIF( IBTYPE .EQ. 2 ) THEN
              WRITE (IWR,'(A,I4,3(A,1PE11.4))') ' Entry No.',I,
     &          '  Time, s: ',BCT(J),
     &          ' Water Surface Elevation,m: ',BCV(J)
            ENDIF
  810     CONTINUE 
  800   CONTINUE
C----------------------------------------------------------------------C
  900 CONTINUE
C----------------------------------------------------------------------C
 9001 FORMAT(/' Input Error: A water boundary condition has ',
     &       'been specified for an interior node.',/,' Node: ',I6,
     &       '  Node Name: ',A)
C----------------------------------------------------------------------C
C     End of RDBCND group.
C----------------------------------------------------------------------C
C
      RETURN
      END

C======================================================================C 
C 
      SUBROUTINE RDRAIN1 
C 
C----------------------------------------------------------------------C 
C 
C     RDRAIN: ReaD RAINfall input group. 
C 
C----------------------------------------------------------------------C 
      USE FILINX 
      USE DRRAIN 
C----------------------------------------------------------------------C 
C     Implicit Double Precision. 
C----------------------------------------------------------------------C 
      IMPLICIT REAL*8 (A-H,O-Z) 
      IMPLICIT INTEGER*4 (I-N) 
C----------------------------------------------------------------------C 
C     Type Declarations. 
C----------------------------------------------------------------------C 
      CHARACTER*240 CHDUM 
      CHARACTER*80 UNITS,TYPE 
      REAL*8 DUMY(5) 
C----------------------------------------------------------------------C 
C     Set table counters to zero. 
C----------------------------------------------------------------------C 
      NSPREC = 0 
      NRAIN  = 0 
C----------------------------------------------------------------------C 
C     Write header to output file. 
C----------------------------------------------------------------------C 
      WRITE (IWR,'(/A)') ' Rainfall Rate' 
      WRITE (IWR,'(A )') ' -------------' 
C----------------------------------------------------------------------C 
C     Found 'Rainfall rate' data 
C----------------------------------------------------------------------C 
      READ (IRD,*) NSZN 
!----------------------------------------------------------------------!
!     ALLOCATING ARRAYS.
!----------------------------------------------------------------------!
      nk = NSZN
      ALLOCATE(IPRCV(nk),KSPREC(2,nk),ISPREC(4,nk)); IPRCV=0; KSPREC=0;
C----------------------------------------------------------------------C 
C     Loop over the number of rainfall rate records. 
C----------------------------------------------------------------------C 
      DO 900 N = 1,NSZN 
        READ (IRD,'(A)') CHDUM 
        CALL LCASE( CHDUM ) 
C----------------------------------------------------------------------C 
C     Read source type. 
C----------------------------------------------------------------------C 
        ITAB = 0 
        ISTART = 1 
        NSPREC = NSPREC + 1 
C----------------------------------------------------------------------C 
C     Read rainfall temporal variation. 
C----------------------------------------------------------------------C  
        CALL RDCHR(ISTART,ICOMMA,CHDUM,TYPE) 
        IF( TYPE(1:8) .EQ. 'constant') THEN  
          IPRCV(NSPREC) = 1               
C----------------------------------------------------------------------C  
C     Read 3 parameter and unit combinations. 
C----------------------------------------------------------------------C  
          DO 100 I = 1, 3 
            CALL RDDPR(ISTART,ICOMMA,CHDUM,DUMY(I))  
            CALL RDCHR(ISTART,ICOMMA,CHDUM,UNITS) 
            CALL RDUNIT(UNITS,DUMY(I)) 
  100     CONTINUE 
C----------------------------------------------------------------------C 
C     Assign rainfall rate variables. 
C----------------------------------------------------------------------C 
          NRAIN = NRAIN + 1  
          IF( NRAIN+1 .GT. LTRAIN ) THEN 
            WRITE (IWR,9020) LTRAIN 
            WRITE (*,'(/2A)') ' INPUT OR COMPILATION ERROR: ',  
     &          'Rainfall Rate Input (check output file).'  
            STOP 
          ENDIF                    
          KSPREC(1,NSPREC) = NRAIN  
          KSPREC(2,NSPREC) = NRAIN+1 
          SPREC(NRAIN)   = DUMY(1) 
          SPREC_TM(NRAIN)  = DUMY(2) 
          NRAIN = NRAIN + 1 
          SPREC_TM(NRAIN)  = DUMY(3) 
        ELSEIF( TYPE(1:7) .EQ. 'tabular' ) THEN 
          ITAB = 1 
C----------------------------------------------------------------------C 
C     Read 5 blank data entries. 
C----------------------------------------------------------------------C 
          DO 200 I = 1,3 
            CALL RDDPR(ISTART,ICOMMA,CHDUM,DUMY(I)) 
            CALL RDCHR(ISTART,ICOMMA,CHDUM,UNITS) 
            CALL RDUNIT(UNITS,DUMY(I)) 
  200     CONTINUE 
C   
          IPRCV(NSPREC) = 2 
C          NRAIN = NRAIN +1 
C          SPREC(NRAIN,1) = DUMY(2) 
        ENDIF 
C----------------------------------------------------------------------C 
C     Read 4 integers (range of rainfall rate) 
C     The rainfall domain KSPREC refers to the field domain directly. 
C----------------------------------------------------------------------C 
        DO 300 I = 1,4 
          CALL RDINT(ISTART,ICOMMA,CHDUM,ISPREC(I,NSPREC)) 
  300   CONTINUE 
C----------------------------------------------------------------------C 
C     Write rainfall rate data to output file. 
C----------------------------------------------------------------------C 
        IF( N .NE. 1 ) WRITE (IWR,'(/)') 
C----------------------------------------------------------------------C 
C     Write rainfall rate records. 
C----------------------------------------------------------------------C 
        WRITE (IWR,'(A)') ' Precipitation Fallout' 
        IF( IPRCV(NSPREC) .EQ. 1 ) THEN 
          WRITE (IWR,'(A)') ' Temporal Variation: Constant' 
          WRITE (IWR,'(A,1PE11.4)') ' Rainfall Rate,      m/s: ', 
     &        SPREC(KSPREC(1,NSPREC)) 
          WRITE (IWR,'(A,1PE11.4)') ' Rainfall Start Time,  s: ', 
     &        SPREC_TM(KSPREC(1,NSPREC)) 
          WRITE (IWR,'(A,1PE11.4)') ' Rainfall Stop Time,   s: ', 
     &        SPREC_TM(KSPREC(2,NSPREC)) 
        ELSEIF( IPRCV(NSPREC) .EQ. 2 ) THEN 
          WRITE (IWR,'(A)') ' Temporal Variation: Tabular' 
        ENDIF   
C----------------------------------------------------------------------C 
C     Write rainfall domain information. 
C----------------------------------------------------------------------C 
        WRITE (IWR,'(A,I4,A,I4)') '   Links: L = ',ISPREC(1,NSPREC), 
     &      ' to ',ISPREC(2,NSPREC) 
        WRITE (IWR,'(A,I4,A,I4)') '  Points: I = ',ISPREC(3,NSPREC), 
     &      ' to ',ISPREC(4,NSPREC) 
C----------------------------------------------------------------------C 
C     Read source tables. 
C----------------------------------------------------------------------C 
        IF( ITAB .EQ. 0 ) GOTO 900  
C----------------------------------------------------------------------C 
C     Length of table. 
C----------------------------------------------------------------------C 
        READ (IRD,*) NLIN 
C----------------------------------------------------------------------C 
C     Table start index. 
C----------------------------------------------------------------------C 
        KSPREC(1,NSPREC) = NRAIN +1 
        KSPREC(2,NSPREC) = KSPREC(1,NSPREC) + NLIN -1 
C----------------------------------------------------------------------C 
C     Read a table ... 
C     Read the abcissa value, the ordinate value and their units. 
C----------------------------------------------------------------------C 
        DO 400 I = 1,NLIN 
          READ (IRD,'(A)') CHDUM 
          CALL LCASE( CHDUM ) 
          ISTART = 1 
          CALL RDDPR(ISTART,ICOMMA,CHDUM,DUMY(1)) 
          CALL RDCHR(ISTART,ICOMMA,CHDUM,UNITS) 
          CALL RDUNIT(UNITS,DUMY(1)) 
C----------------------------------------------------------------------C 
C     Read a table ... 
C     Read the ordinate value and units. 
C----------------------------------------------------------------------C
          DO 450 J=2,2 
            CALL RDDPR(ISTART,ICOMMA,CHDUM,DUMY(J)) 
            CALL RDCHR(ISTART,ICOMMA,CHDUM,UNITS) 
            CALL RDUNIT(UNITS,DUMY(J))
  450     CONTINUE
C----------------------------------------------------------------------C 
C     Assign source variables. 
C----------------------------------------------------------------------C 
          NRAIN = NRAIN +1 
          IF( NRAIN .GT. LTRAIN ) THEN 
            WRITE (IWR,9020) LTRAIN 
            WRITE (*,'(/2A)') ' INPUT OR COMPILATION ERROR: ', 
     &          'Rainfall Rate Input (check output file).' 
            STOP 
          ENDIF 
          SPREC(NRAIN)    = DUMY(2) 
          SPREC_TM(NRAIN) = DUMY(1) 
  400   CONTINUE 
C----------------------------------------------------------------------C 
C     Write source tables. 
C----------------------------------------------------------------------C 
        WRITE (IWR,'(A)') ' Precipitation Fallout Table ' 
        DO 500 I = 1,NLIN 
          J = KSPREC(1,NSPREC) + I - 1 
          WRITE (IWR,'(A,I4,4(A,1PE11.4))') ' Entry No.',I, 
     &          '  Time, s: ',SPREC_TM(J), 
     &          '  Rainfall Rate, m/s: ',SPREC(J) 
  500   CONTINUE 
C----------------------------------------------------------------------C 
  900 CONTINUE 
C----------------------------------------------------------------------C 
C     Format Statements. 
C----------------------------------------------------------------------C 
 9020 FORMAT(/' Input or Compilation Error: The number of rainfall ', 
     &        'rate table entries plus ',/,' doubled constant ', 
     &        'rainfall specified exceeds the maximum set in parameter ' 
     &        'LTRAIN [',I6,'].') 
C----------------------------------------------------------------------C 
C     End of RDRAIN group. 
C----------------------------------------------------------------------C 
C 
      RETURN 
      END 
 
!= =====================================================================C
C
      SUBROUTINE RAINFL( IR )
C
C----------------------------------------------------------------------C
C
C     RAINFL: RAINFaLl rate for shallow water equations.
C
C----------------------------------------------------------------------C 
      USE DRRAIN
      USE DRWFLD
      USE FLDINX 
      USE SOLVAR
C----------------------------------------------------------------------C
C     Implicit Double Precision.
C----------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C----------------------------------------------------------------------C
C     Loop over sources.
C----------------------------------------------------------------------C
      DO 100 NS = 1,NSPREC
        NSB = KSPREC(1,NS)
        NSE = KSPREC(2,NS)
C----------------------------------------------------------------------C
C     Constant source.
C----------------------------------------------------------------------C
        IF( IPRCV(NS) .EQ. 1 ) THEN
          IF( TIME .GT. SPREC_TM(NSB) .AND. 
     &        TIME-DT .LT. SPREC_TM(NSE) ) THEN
            TMIN = MAX( TIME-DT, SPREC_TM(NSB) )
            TMAX = MIN( TIME, SPREC_TM(NSE) )
            DTSRC = TMAX -TMIN 
            SRCX = SPREC(NSB)*DTSRC/DT
          ELSE
            SRCX = 0.D+0
          ENDIF
C----------------------------------------------------------------------C
C     Tabular source.
C----------------------------------------------------------------------C
        ELSEIF( IPRCV(NS) .EQ. 2 ) THEN
          IF( TIME .GT. SPREC_TM(NSB) .AND.
     &        TIME-DT .LT. SPREC_TM(NSE)) THEN
            SRCX = 0.D+0
            DO 140 I = NSB+1,NSE
              IF( TIME .GT. SPREC_TM(I-1) .AND.
     &            TIME-DT .LT. SPREC_TM(I)) THEN
                TMIN = MAX( TIME-DT, SPREC_TM(I-1) )
                TMAX = MIN( TIME, SPREC_TM(I) )
                DTSRC = TMAX -TMIN
                TMID = 0.5*(TMIN +TMAX)
                NLEN = 2
                SRV = SPREC(I)
C                CALL LININT(TMID,SRV,ZZ,SPREC_TM(I-1),SPREC(I-1),NLEN)
                SRCX = SRCX + SRV*DTSRC/DT 
              ENDIF
  140       CONTINUE   
          ELSE
            SRCX = 0.D+0
          ENDIF
        ENDIF
C----------------------------------------------------------------------C 
C     Loop over source domain.
C----------------------------------------------------------------------C
        DO 110 I = ISPREC(1,NS),ISPREC(2,NS)
          DO 120 J = ISPREC(3,NS),ISPREC(4,NS) 
            N = ND(I,J)
            RAIN(N) = SRCX
  120     CONTINUE
  110   CONTINUE
  100 CONTINUE
C----------------------------------------------------------------------C 
C     End of RAINFL group.
C----------------------------------------------------------------------C
C
      RETURN
      END

!= =====================================================================C
C 
      SUBROUTINE RDEVAP1 
C 
C----------------------------------------------------------------------C 
C 
C     RDEVAP: ReaD EVAPotranspiration input group. 
C 
C----------------------------------------------------------------------C 
      USE DREVAP
      USE FILINX 
C----------------------------------------------------------------------C 
C     Implicit Double Precision. 
C----------------------------------------------------------------------C 
      IMPLICIT REAL*8 (A-H,O-Z) 
      IMPLICIT INTEGER*4 (I-N) 
C----------------------------------------------------------------------C 
C     Type Declarations. 
C----------------------------------------------------------------------C 
      CHARACTER*240 CHDUM 
      CHARACTER*80 UNITS,TYPE 
      REAL*8 DUMY(5) 
C----------------------------------------------------------------------C 
C     Set table counters to zero. 
C----------------------------------------------------------------------C 
      NSEVAP = 0 
      NEVAP  = 0 
C----------------------------------------------------------------------C 
C     Write header to output file. 
C----------------------------------------------------------------------C 
      WRITE (IWR,'(/A)') ' Evapotranspiration Rate' 
      WRITE (IWR,'(A )') ' -----------------------' 
C----------------------------------------------------------------------C 
C     Found 'Evapotranspiration rate' data 
C----------------------------------------------------------------------C 
      READ (IRD,*) NSZN 
C----------------------------------------------------------------------C 
C     Loop over the number of evapotranspiration rate records. 
C----------------------------------------------------------------------C 
      DO 900 N = 1,NSZN 
        READ (IRD,'(A)') CHDUM 
        CALL LCASE( CHDUM ) 
C----------------------------------------------------------------------C 
C     Read source type. 
C----------------------------------------------------------------------C 
        ITAB = 0 
        ISTART = 1 
        NSEVAP = NSEVAP + 1 
C----------------------------------------------------------------------C 
C     Check number of evapotranspiration zones specified against 
C     parameter limit. 
C----------------------------------------------------------------------C 
!        IF( NSEVAP .GT. LDEVAP ) THEN 
!          WRITE (IWR,9000) NSEVAP,LDEVAP 
!          WRITE (*,'(/2A)') ' INPUT OR COMPILATION ERROR: ', 
!     &      'Evapotranspiration Rate Input (check output file).' 
!          STOP 
!        ENDIF  
C----------------------------------------------------------------------C 
C     Read evapotranspiration temporal variation. 
C----------------------------------------------------------------------C 
        CALL RDCHR(ISTART,ICOMMA,CHDUM,TYPE) 
        IF( TYPE(1:8) .EQ. 'constant') THEN 
          IEVPV(NSEVAP) = 1              
C----------------------------------------------------------------------C 
C     Read 3 parameter and unit combinations. 
C----------------------------------------------------------------------C  
          DO 100 I = 1, 3  
            CALL RDDPR(ISTART,ICOMMA,CHDUM,DUMY(I)) 
            CALL RDCHR(ISTART,ICOMMA,CHDUM,UNITS) 
            CALL RDUNIT(UNITS,DUMY(I)) 
  100     CONTINUE 
C----------------------------------------------------------------------C 
C     Assign evapotranspiration rate variables. 
C----------------------------------------------------------------------C 
          NEVAP = NEVAP + 1  
!          IF( NEVAP+1 .GT. LTEVAP ) THEN 
!            WRITE (IWR,9020) LTEVAP 
!            WRITE (*,'(/2A)') ' INPUT OR COMPILATION ERROR: ', 
!     &          'Evapotranspiration Rate Input (check output file).' 
!            STOP 
!          ENDIF                    
          KSEVAP(1,NSEVAP) = NEVAP 
          KSEVAP(2,NSEVAP) = NEVAP+1 
          SEVAP(NEVAP)     = DUMY(1) 
          SEVAP_TM(NEVAP)  = DUMY(2) 
          NEVAP = NEVAP + 1 
          SEVAP_TM(NEVAP)  = DUMY(3) 
        ELSEIF( TYPE(1:7) .EQ. 'tabular' ) THEN 
          ITAB = 1 
C----------------------------------------------------------------------C 
C     Read 3 blank data entries. 
C----------------------------------------------------------------------C 
          DO 200 I = 1,3 
            CALL RDDPR(ISTART,ICOMMA,CHDUM,DUMY(I)) 
            CALL RDCHR(ISTART,ICOMMA,CHDUM,UNITS) 
            CALL RDUNIT(UNITS,DUMY(I)) 
  200     CONTINUE 
C   
          IEVPV(NSEVAP) = 2 
C          NEVAP = NEVAP +1 
C          SEVAP(NEVAP,1) = DUMY(2) 
        ENDIF 
C----------------------------------------------------------------------C 
C     Read 4 integers (range of evapotranspiration rate) 
C     The evapotranspiration domain KSEVAP refers to the field 
C     domain directly. 
C----------------------------------------------------------------------C 
        DO 300 I = 1,4 
          CALL RDINT(ISTART,ICOMMA,CHDUM,ISEVAP(I,NSEVAP)) 
  300   CONTINUE 
C----------------------------------------------------------------------C 
C     Write evapotranspiration rate data to output file. 
C----------------------------------------------------------------------C 
        IF( N .NE. 1 ) WRITE (IWR,'(/)') 
C----------------------------------------------------------------------C 
C     Write evapotranspiration rate records. 
C----------------------------------------------------------------------C 
        WRITE (IWR,'(A)') ' Evapotranspiration' 
        IF( IEVPV(NSEVAP) .EQ. 1 ) THEN 
          WRITE (IWR,'(A)') ' Temporal Variation: Constant' 
          WRITE (IWR,'(A,1PE11.4)')  
     &      ' Evapotranspiration Rate,     m/s: ',
     &        SEVAP(KSEVAP(1,NSEVAP)) 
          WRITE (IWR,'(A,1PE11.4)') 
     &      ' Evapotranspiration Start Time, s: ', 
     &        SEVAP_TM(KSEVAP(1,NSEVAP)) 
          WRITE (IWR,'(A,1PE11.4)') 
     &      ' Evapotranspiration Stop Time,  s: ', 
     &        SEVAP_TM(KSEVAP(2,NSEVAP)) 
        ELSEIF( IEVPV(NSEVAP) .EQ. 2 ) THEN 
          WRITE (IWR,'(A)') ' Temporal Variation: Tabular' 
        ENDIF   
C----------------------------------------------------------------------C 
C     Write evapotranspiration domain information. 
C----------------------------------------------------------------------C 
        WRITE (IWR,'(A,I4,A,I4)') '   Links : L = ',ISEVAP(1,NSEVAP), 
     &      ' to ',ISEVAP(2,NSEVAP) 
        WRITE (IWR,'(A,I4,A,I4)') '   Points: I = ',ISEVAP(3,NSEVAP), 
     &      ' to ',ISEVAP(4,NSEVAP) 
C----------------------------------------------------------------------C 
C     Read source tables. 
C----------------------------------------------------------------------C 
        IF( ITAB .EQ. 0 ) GOTO 900  
C----------------------------------------------------------------------C 
C     Length of table. 
C----------------------------------------------------------------------C 
        READ (IRD,*) NLIN 
C----------------------------------------------------------------------C 
C     Table start index. 
C----------------------------------------------------------------------C 
        KSEVAP(1,NSEVAP) = NEVAP +1 
        KSEVAP(2,NSEVAP) = KSEVAP(1,NSEVAP) + NLIN -1 
C----------------------------------------------------------------------C 
C     Read a table ... 
C     Read the abcissa value, the ordinate value and their units. 
C----------------------------------------------------------------------C 
        DO 400 I = 1,NLIN 
          READ (IRD,'(A)') CHDUM 
          CALL LCASE( CHDUM ) 
          ISTART = 1 
          CALL RDDPR(ISTART,ICOMMA,CHDUM,DUMY(1)) 
          CALL RDCHR(ISTART,ICOMMA,CHDUM,UNITS) 
          CALL RDUNIT(UNITS,DUMY(1)) 
C----------------------------------------------------------------------C 
C     Read a table ... 
C     Read the ordinate value and units. 
C----------------------------------------------------------------------C
          DO 450 J=2,4 
            CALL RDDPR(ISTART,ICOMMA,CHDUM,DUMY(J)) 
            CALL RDCHR(ISTART,ICOMMA,CHDUM,UNITS) 
            CALL RDUNIT(UNITS,DUMY(J))
  450     CONTINUE
C----------------------------------------------------------------------C 
C     Assign source variables. 
C----------------------------------------------------------------------C 
          NEVAP = NEVAP +1 
!          IF( NEVAP .GT. LTEVAP ) THEN 
!            WRITE (IWR,9020) LTEVAP 
!            WRITE (*,'(/2A)') ' INPUT OR COMPILATION ERROR: ', 
!     &          'Evapotranspiration Rate Input (check output file).' 
!            STOP 
!          ENDIF 
          SEVAP(NEVAP)     = DUMY(2) 
          SEVAP_TM(NEVAP)  = DUMY(1) 
  400   CONTINUE 
C----------------------------------------------------------------------C 
C     Write source tables. 
C----------------------------------------------------------------------C 
        WRITE (IWR,'(A)') ' Evapotranspiration Table ' 
        DO 500 I = 1,NLIN 
          J = KSEVAP(1,NSEVAP) + I - 1 
          WRITE (IWR,'(A,I4,2(A,1PE11.4))') ' Entry No.',I, 
     &          '  Time, s: ',SEVAP_TM(J), 
     &          '  Evapotranspiration Rate,      m/s: ',SEVAP(J) 
  500   CONTINUE 
C----------------------------------------------------------------------C 
  900 CONTINUE 
C----------------------------------------------------------------------C 
C     End of RDEVAP group. 
C----------------------------------------------------------------------C 
C 
      RETURN 
      END 
 
C======================================================================C
C
      SUBROUTINE EVAPTR
C
C----------------------------------------------------------------------C
C
C     EVAPTR: EVAPoTRanspiration rate for shallow water equations.
C
C----------------------------------------------------------------------C 
      USE DREVAP
      USE DRWFLD
      USE FLDINX 
      USE SOLVAR
C----------------------------------------------------------------------C
C     Implicit Double Precision.
C----------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C----------------------------------------------------------------------C
C     Loop over sources.
C----------------------------------------------------------------------C
      SCL = 0.D+0
      DO 100 NS = 1,NSEVAP
        NSB = KSEVAP(1,NS)
        NSE = KSEVAP(2,NS)
C----------------------------------------------------------------------C
C     Constant source.
C----------------------------------------------------------------------C
        IF( IEVPV(NS) .EQ. 1 ) THEN
          IF( TIME .GT. SEVAP_TM(NSB) .AND. 
     &        TIME-DT .LT. SEVAP_TM(NSE) ) THEN
            TMIN = MAX( TIME-DT, SEVAP_TM(NSB) )
            TMAX = MIN( TIME, SEVAP_TM(NSE) )
            DTSRC = TMAX -TMIN 
            SRCX = SEVAP(NSB)*DTSRC/DT
          ELSE
            SRCX = 0.D+0
            SCL  = 0.D+0
          ENDIF
C----------------------------------------------------------------------C
C     Tabular source.
C----------------------------------------------------------------------C
        ELSEIF( IEVPV(NS) .EQ. 2 ) THEN
          IF( TIME .GT. SEVAP_TM(NSB) .AND.
     &        TIME-DT .LT. SEVAP_TM(NSE)) THEN
            SRCX = 0.D+0
            DO 140 I = NSB+1,NSE
              IF( TIME .GT. SEVAP_TM(I-1) .AND.
     &            TIME-DT .LT. SEVAP_TM(I)) THEN
                TMIN = MAX( TIME-DT, SEVAP_TM(I-1) )
                TMAX = MIN( TIME, SEVAP_TM(I) )
                DTSRC = TMAX -TMIN
                TMID = 0.5*(TMIN +TMAX)
                NLEN = 2
                SRV = SEVAP(I)
C                CALL LININT(TMID,SRV,ZZ,SEVAP_TM(I-1),SEVAP(I-1),NLEN)
                SRCX = SRCX + SRV*DTSRC/DT 
              ENDIF
  140       CONTINUE   
          ELSE
            SRCX = 0.D+0
            SCL  = 0.D+0
          ENDIF
        ENDIF
C----------------------------------------------------------------------C
C     Loop over source domain.
C----------------------------------------------------------------------C
        DO 110 I = ISEVAP(1,NS),ISEVAP(2,NS)
          DO 120 J = ISEVAP(3,NS),ISEVAP(4,NS)
            N = ND(I,J)
            EVAP(N) = 0.D+0
            IF( IXP(N) .LE. 0 ) GOTO 120
            EVAP(N) = SRCX
  120     CONTINUE
  110   CONTINUE
  100 CONTINUE
C----------------------------------------------------------------------C
C     End of EVAPTR group.
C----------------------------------------------------------------------C
C
      RETURN
      END

C======================================================================C
C
      SUBROUTINE RDDFPR
C
C----------------------------------------------------------------------C
C
C     RDDFPR: ReaD DeFault PaRameters input group.
C
C----------------------------------------------------------------------C  
      USE DRSOLV 
      USE FILINX
C----------------------------------------------------------------------C
C     Implicit Double Precision.
C----------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C----------------------------------------------------------------------C
C     Type Declarations.
C----------------------------------------------------------------------C
      CHARACTER*240 CHDUM
      CHARACTER*80 UNITS,ATYPE,PEROPT
C----------------------------------------------------------------------C
C     Write header to output file.
C----------------------------------------------------------------------C
      WRITE (IWR,'(/A)') ' Default Parameters'
      WRITE (IWR,'(A )') ' ------------------'
C----------------------------------------------------------------------C 
C     Read Weighting Factors for St.Venant Equations. 
C----------------------------------------------------------------------C 
      READ (IRD,'(A)') CHDUM 
      CALL LCASE(CHDUM) 
      ISTART = 1 
      CALL RDDPR(ISTART,ICOMMA,CHDUM,ALPHA) 
      CALL RDDPR(ISTART,ICOMMA,CHDUM,BETA) 
      CALL RDDPR(ISTART,ICOMMA,CHDUM,PHI) 
      CALL RDDPR(ISTART,ICOMMA,CHDUM,THETA)
      IF( ALPHA .LT. 0.D+0 .OR. ALPHA .GT. 1.D+0 ) ALPHA = 1.D+0  
      IF( BETA .LT. 0.D+0 .OR. BETA .GT. 1.D+0 )    BETA = 1.D+0  
      IF( PHI .LT. 0.D+0 .OR. PHI .GT. 1.D+0 )       PHI = 0.5D+0  
      IF( THETA .LT. 0.D+0 .OR. THETA .GT. 1.D+0 ) THETA = 0.55D+0  
C----------------------------------------------------------------------C 
C     Write Weighting Factors for St.Venant Equations. 
C----------------------------------------------------------------------C 
      WRITE (IWR,'(/A)')   ' Weighting Factors for St.Venant Equations' 
      WRITE (IWR,'(A,1PE11.4)')
     &             '   Momentum Correction Factor   : ',ALPHA 
      WRITE (IWR,'(A,1PE11.4)')'   Energy Slope Weighting Factor: ',BETA 
      WRITE (IWR,'(A,1PE11.4)')'   Space Weighting Factor       : ',PHI 
      WRITE (IWR,'(A,1PE11.4)')
     &      '   Time-Weighting Parameter     : ',THETA 
C----------------------------------------------------------------------C
C     Read Hydraulic Resistance Properties.
C----------------------------------------------------------------------C
      CALL RDCHR(ISTART,ICOMMA,CHDUM,PEROPT)
C----------------------------------------------------------------------C
C     Darcy-Weisbach and Manning.
C----------------------------------------------------------------------C
      IF( PEROPT(1:14) .EQ. 'darcy-weisbach') THEN
        iIFRIC = 1
        CALL RDDPR(ISTART,ICOMMA,CHDUM,rFRIC)
        CALL RDCHR(ISTART,ICOMMA,CHDUM,UNITS) 
        CALL RDUNIT(UNITS,rFRIC) 
        WRITE (IWR,'(/2A)') ' Hydraulic Resistance ',
     &     'Function: Darcy-Weisbach Equation'
        WRITE (IWR,'(A,1PE11.4)') 
     &    '   Surface Friction Coefficient             : ',rFRIC
      ELSEIF( PEROPT(1:5) .EQ. 'chezy') THEN
        iIFRIC = 2
        CALL RDDPR(ISTART,ICOMMA,CHDUM,rFRIC)
        CALL RDCHR(ISTART,ICOMMA,CHDUM,UNITS) 
        CALL RDUNIT(UNITS,rFRIC) 
        WRITE (IWR,'(/2A)') ' Hydraulic Resistance ',
     &     'Function: Chezy Law'
        WRITE (IWR,'(A,1PE11.4)') 
     &    '   Chezy Friction Coefficient,        m^0.5/s: ',rFRIC
      ELSEIF( PEROPT(1:7) .EQ. 'manning') THEN
        iIFRIC = 3
        CALL RDDPR(ISTART,ICOMMA,CHDUM,rFRIC)
        CALL RDCHR(ISTART,ICOMMA,CHDUM,UNITS) 
        CALL RDUNIT(UNITS,rFRIC) 
        WRITE (IWR,'(/2A)') ' Hydraulic Resistance ',
     &     'Function: Manning Equation'
        WRITE (IWR,'(A,1PE11.4)') 
     &    '   Manning''s Roughness Coefficient, s/m^0.33: ',rFRIC 
      ELSE 
        iIFRIC = 3  
        rFRIC = 0.2D-1 
        WRITE (IWR,'(/2A)') ' Hydraulic Resistance ',
     &     'Function: Manning Equation'
        WRITE (IWR,'(A,1PE11.4)') 
     &    '   Manning''s Roughness Coefficient, s/m^0.33: ',rFRIC 
      ENDIF
C----------------------------------------------------------------------C 
C     Read depth-averaged longitudinal and transverse  
C     turbulent diffusivities.  
C----------------------------------------------------------------------C 
C        CALL RDDPR(ISTART,ICOMMA,CHDUM,DFSL(IPLANT)) 
C        CALL RDDPR(ISTART,ICOMMA,CHDUM,DFST(IPLANT)) 
C        CALL RDCHR(ISTART,ICOMMA,CHDUM,UNITS) 
C        CALL RDUNIT(UNITS,DFSL(IPLANT)) 
C        CALL RDUNIT(UNITS,DFST(IPLANT)) 
C        WRITE (IWR,'(2A,1PE11.4)') '   Longitudinal Turbulent ', 
C     &    'Diffusion Coefficient, m^2/s: ',DFSL(IPLANT) 
C        WRITE (IWR,'(2A,1PE11.4)') '   Transverse Turbulent ', 
C     &    'Diffusion Coefficient,   m^2/s: ',DFST(IPLANT) 
C----------------------------------------------------------------------C
  900 CONTINUE
C----------------------------------------------------------------------C
C     End of RDDFPR group.
C----------------------------------------------------------------------C
C
      RETURN
      END

C======================================================================C
C
      SUBROUTINE RDLKSR
C
C----------------------------------------------------------------------C
C
C     RDLKSR: ReaD LinK SouRces and sinks input group.
C
C----------------------------------------------------------------------C 
      USE FILINX 
      USE LOGCLS
      USE LQSORC 
      USE SPSOUR 
      use intersub
C----------------------------------------------------------------------C
C     Implicit Double Precision.
C----------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

      INCLUDE 'utils.inc' ! subroutine interfaces
C----------------------------------------------------------------------C
C     Type Declarations.
C----------------------------------------------------------------------C
      CHARACTER*240 CHDUM
      CHARACTER*80 ADUM,UNITS,TYPE  
      CHARACTER STR*35
      REAL*8 DUMY(3) 
      integer*4 NSTSP(3)  
C----------------------------------------------------------------------C
C     Set table counters to zero.
C----------------------------------------------------------------------C
      NSRCLQ = 0
      NSTLQ  = 0
      NSRCSP = 0
      NSTSP  = 0
C----------------------------------------------------------------------C
C     Write header to output file.
C----------------------------------------------------------------------C
      WRITE (IWR,'(/A)') ' Link Sources and Sinks'
      WRITE (IWR,'(A )') ' ----------------------'
C----------------------------------------------------------------------C
C     Found 'Sources and Sinks' data
C----------------------------------------------------------------------C
      READ (IRD,*) NSZN
C----------------------------------------------------------------------C
C     Allocate sources data
C----------------------------------------------------------------------C 
      Allocate(ISLQV(NSZN),ISSPV(NSZN,3),ISPTP(NSZN,3)) 
      ISLQV = 0; ISSPV = 0; ISPTP = 0; 
      Allocate(ISRCLQ(4,NSZN),KSRCLQ(2,NSZN))  
      ISRCLQ = 0; KSRCLQ = 0; 
      Allocate(ISRCSP(4,NSZN,3),KSRCSP(2,NSZN,3))  
      ISRCSP = 0; KSRCSP = 0; 
      Allocate(SRCLQ(NSZN),SRCLQ_TM(2*NSZN)); 
      SRCLQ = 0.d0; SRCLQ_TM = 0.d0;   
      Allocate(SRCSP1(NSZN),SRCSP_TM1(2*NSZN))  
      SRCSP1 = 0.d0; SRCSP_TM1 = 0.d0;
      Allocate(SRCSP2(NSZN),SRCSP_TM2(2*NSZN))  
      SRCSP2 = 0.d0; SRCSP_TM2 = 0.d0;
      Allocate(SRCSP3(NSZN),SRCSP_TM3(2*NSZN))  
      SRCSP3 = 0.d0; SRCSP_TM3 = 0.d0;
C----------------------------------------------------------------------C
C     Loop over the number of sources and/or sinks.
C----------------------------------------------------------------------C
      DO 900 N = 1,NSZN
        READ (IRD,'(A)') CHDUM
        CALL LCASE( CHDUM )
C----------------------------------------------------------------------C
C     Read source type.
C----------------------------------------------------------------------C
        ITAB = 0
        ISTYPE = 2
        ISTART = 1
        CALL RDCHR(ISTART,ICOMMA,CHDUM,ADUM)
        IF( ADUM(1:11) .EQ. 'liquid disc') THEN  
          LQSOUR = .true.  
          ISTYPE = 1
          NSRCLQ = NSRCLQ + 1
        ELSEIF( ADUM(1:16) .EQ. 'species liq conc' ) THEN
          CLSOUR = .true.  
          NSRCSP(1) = NSRCSP(1) + 1
          ISPTP(NSRCSP(1),1) = 1 
          is = 1
        ELSEIF( ADUM(1:16) .EQ. 'species liq disc' ) THEN
          CLSOUR = .true.  
          NSRCSP(1) = NSRCSP(1) + 1
          ISPTP(NSRCSP(1),1) = 2 
          is = 1
        ELSEIF( ADUM(1:16) .EQ. 'species pat conc' ) THEN
          CPSOUR = .true.  
          NSRCSP(2) = NSRCSP(2) + 1
          ISPTP(NSRCSP(2),2) = 1  
          is = 2
        ELSEIF( ADUM(1:16) .EQ. 'species pat disc' ) THEN
          CPSOUR = .true.  
          NSRCSP(2) = NSRCSP(2) + 1
          ISPTP(NSRCSP(2),2) = 2  
          is = 2
        ELSEIF( ADUM(1:13) .EQ. 'sediment conc' ) THEN
          SDSOUR = .true.  
          NSRCSP(3) = NSRCSP(3) + 1
          ISPTP(NSRCSP(3),3) = 1   
          is = 3
        ELSEIF( ADUM(1:13) .EQ. 'sediment disc' ) THEN
          SDSOUR = .true.  
          NSRCSP(3) = NSRCSP(3) + 1
          ISPTP(NSRCSP(3),3) = 2   
          is = 3
        ENDIF
C----------------------------------------------------------------------C
C     Read source temporal variation.
C----------------------------------------------------------------------C
        CALL RDCHR(ISTART,ICOMMA,CHDUM,TYPE)
        IF( TYPE(1:8) .EQ. 'constant') THEN
          IF( ISTYPE .EQ. 1 ) THEN
            ISLQV(NSRCLQ) = 1             
          ELSE
            ISSPV(NSRCSP(is),is) = 1
          ENDIF
C----------------------------------------------------------------------C
C     Read 3 parameter and unit combinations.
C----------------------------------------------------------------------C
          DO 100 I = 1, 3
            CALL RDDPR(ISTART,ICOMMA,CHDUM,DUMY(I))
            CALL RDCHR(ISTART,ICOMMA,CHDUM,UNITS)
            CALL RDUNIT(UNITS,DUMY(I))
  100     CONTINUE
C----------------------------------------------------------------------C
C     Assign source variables.
C----------------------------------------------------------------------C
          IF( ISTYPE .EQ. 1 ) THEN 
            NSTLQ = NSTLQ + 1 
            KSRCLQ(1,NSRCLQ) = NSTLQ
            KSRCLQ(2,NSRCLQ) = NSTLQ+1
            SRCLQ(NSTLQ)    = DUMY(1)
            SRCLQ_TM(NSTLQ) = DUMY(2)
            NSTLQ = NSTLQ + 1
            SRCLQ_TM(NSTLQ) = DUMY(3)
          ELSE
            NSTSP(is) = NSTSP(is) + 1 
            KSRCSP(1,NSRCSP(is),is) = NSTSP(is)
            KSRCSP(2,NSRCSP(is),is) = NSTSP(is) +1 
            if( is == 1 ) then   
              SRCSP1(NSTSP(is))    = DUMY(1)
              SRCSP_TM1(NSTSP(is)) = DUMY(2)
              NSTSP(is) = NSTSP(is) + 1
              SRCSP_TM1(NSTSP(is)) = DUMY(3)
            elseif( is == 2 ) then   
              SRCSP2(NSTSP(is))    = DUMY(1)
              SRCSP_TM2(NSTSP(is)) = DUMY(2)
              NSTSP(is) = NSTSP(is) + 1
              SRCSP_TM2(NSTSP(is)) = DUMY(3)
            elseif( is == 3 ) then
              SRCSP3(NSTSP(is))    = DUMY(1)
              SRCSP_TM3(NSTSP(is)) = DUMY(2)
              NSTSP(is) = NSTSP(is) + 1
              SRCSP_TM3(NSTSP(is)) = DUMY(3)
            endif
            
          ENDIF
       ELSEIF( TYPE(1:7) .EQ. 'tabular' ) THEN
          ITAB = 1
          IF( ISTYPE .EQ. 1 ) THEN
            ISLQV(NSRCLQ) = 2             
          ELSE
            ISSPV(NSRCSP(is),is) = 2
          ENDIF
C----------------------------------------------------------------------C
C     Read 6 blank data entries.
C----------------------------------------------------------------------C
          DO 200 I = 1,6
            CALL RDCHR(ISTART,ICOMMA,CHDUM,ADUM)
  200     CONTINUE
        ENDIF
C----------------------------------------------------------------------C
C     Read 4 integers (range of source/sink)
C     The source domain KSRCLQ and KSRCSP refer to the field domain
C     directly.
C----------------------------------------------------------------------C
        DO 300 I = 1,4
          IF( ISTYPE .EQ. 1 ) THEN
            CALL RDINT(ISTART,ICOMMA,CHDUM,ISRCLQ(I,NSRCLQ))
          ELSE
            CALL RDINT(ISTART,ICOMMA,CHDUM,ISRCSP(I,NSRCSP(is),is))
          ENDIF  
  300   CONTINUE
C----------------------------------------------------------------------C
C     Write source and sink data to output file.
C----------------------------------------------------------------------C
        IF( N .NE. 1 ) WRITE (IWR,'(/)')
C----------------------------------------------------------------------C
C     Write liquid sources.
C----------------------------------------------------------------------C
        IF( ISTYPE .EQ. 1 ) THEN
          WRITE (IWR,'(A)') ' Liquid Discharge Source'
          IF( ISLQV(NSRCLQ) .EQ. 1 ) THEN
            WRITE (IWR,'(A)') ' Temporal Variation: Constant'
            WRITE (IWR,'(A,1PE11.4)') ' Source Rate,          m^3/s: ',
     &        SRCLQ(KSRCLQ(1,NSRCLQ))
            WRITE (IWR,'(A,1PE11.4)') ' Source Start Time,        s: ',
     &        SRCLQ_TM(KSRCLQ(1,NSRCLQ))
            WRITE (IWR,'(A,1PE11.4)') ' Source Stop Time,         s: ',
     &        SRCLQ_TM(KSRCLQ(2,NSRCLQ))
          ELSEIF( ISLQV(NSRCLQ) .EQ. 2 ) THEN
            WRITE (IWR,'(A)') ' Temporal Variation: Tabular'
          ENDIF
C----------------------------------------------------------------------C
C     Write species sources.
C----------------------------------------------------------------------C
        ELSE  
          if( is == 1 ) then  
            if( isptp(nsrcsp(is),is) == 1 ) then  
              WRITE (IWR,'(A)') ' Species Liquid Concentration Source'  `
              STR = ' Species Concentration,   1/m^3: '  
            else   
              WRITE (IWR,'(A)') ' Species Liquid Discharge Source' 
              STR = ' Species Discharge,         1/s: '  
            endif     
          elseif( is == 2 ) then  
            if( isptp(nsrcsp(is),is) == 1 ) then  
              WRITE (IWR,'(A)')   
     &          ' Species Particulate Concentration Source' 
              STR = ' Species Concentration,   1/m^3: '  
            else   
              WRITE (IWR,'(A)') ' Species Particulate Discharge Source' 
              STR = ' Species Discharge,         1/s: '  
            endif     
          elseif( is == 3 ) then  
            if( isptp(nsrcsp(is),is) == 1 ) then  
              WRITE (IWR,'(A)') ' Sediment Concentration Source' 
              STR = ' Sediment Concentration, kg/m^3: '  
            else   
              WRITE (IWR,'(A)') ' Sediment Discharge Source' 
              STR = ' Sediment Discharge,       kg/s: '  
            endif     
          endif      
            
          IF( ISSPV(NSRCSP(is),is) == 1 ) THEN
            WRITE (IWR,'(A)') ' Temporal Variation: Constant'
            WRITE (IWR,'(A,1PE11.4)') STR,DUMY(1)
!     &        SRCSP(KSRCSP(1,NSRCSP(is),is),is)
            WRITE (IWR,'(A,1PE11.4)') ' Source Start Time,        s: ',
     &        DUMY(2)
            WRITE (IWR,'(A,1PE11.4)') ' Source Stop Time,         s: ',
     &        DUMY(3)
          ELSEIF( ISSPV(NSRCSP(is),is) .EQ. 2 ) THEN
            WRITE (IWR,'(A)') ' Temporal Variation: Tabular'
          ENDIF
        ENDIF  
C----------------------------------------------------------------------C
C     Write source domain information.
C----------------------------------------------------------------------C
        IF( ISTYPE .EQ. 1 ) THEN
          WRITE (IWR,'(A,I4,A,I4)') '   Link : L = ',ISRCLQ(1,NSRCLQ),
     &      ' to ',ISRCLQ(2,NSRCLQ)
          WRITE (IWR,'(A,I4,A,I4)') '   Point: I = ',ISRCLQ(3,NSRCLQ),
     &      ' to ',ISRCLQ(4,NSRCLQ)
        ELSE
          WRITE (IWR,'(A,I4,A,I4)') '   Link : L = ',
     &      ISRCSP(1,NSRCSP(is),is),' to ',ISRCSP(2,NSRCSP(is),is)
          WRITE (IWR,'(A,I4,A,I4)') '   Point: I = ',
     &      ISRCSP(3,NSRCSP(is),is),' to ',ISRCSP(4,NSRCSP(is),is)
        ENDIF
C----------------------------------------------------------------------C
C     Read source tables.
C----------------------------------------------------------------------C
        IF( ITAB .EQ. 0 ) GOTO 900 
C----------------------------------------------------------------------C
C     Length of table.
C----------------------------------------------------------------------C
        READ (IRD,*) NLIN
!------------------------------------------------------------------------------C
!     REALLOCATE  BCVS,BCTS arrays.
!------------------------------------------------------------------------------C
        IF( ISTYPE .EQ. 1 ) THEN  
          nk = NSTLQ +NLIN;
          call ReallocateR8(NSTLQ,nk,SRCLQ)
          call ReallocateR8(NSTLQ,nk,SRCLQ_TM)
        ELSE   
          nk = NSTSP(is) +NLIN;  
          if( is == 1 ) then 
            call ReallocateR8(NSTSP(is),nk,SRCSP1)
            call ReallocateR8(NSTSP(is),nk,SRCSP_TM1)  
          elseif( is == 2 ) then  
            call ReallocateR8(NSTSP(is),nk,SRCSP2)
            call ReallocateR8(NSTSP(is),nk,SRCSP_TM2)  
          elseif( is == 3 ) then   
            call ReallocateR8(NSTSP(is),nk,SRCSP3)
            call ReallocateR8(NSTSP(is),nk,SRCSP_TM3)  
          endif    
        ENDIF    
C----------------------------------------------------------------------C
C     Table start index.
C----------------------------------------------------------------------C 
        IF( ISTYPE .EQ. 1 ) THEN
          KSRCLQ(1,NSRCLQ) = NSTLQ +1
          KSRCLQ(2,NSRCLQ) = KSRCLQ(1,NSRCLQ) + NLIN -1 
        ELSE
          KSRCSP(1,NSRCSP(is),is) = NSTSP(is) +1
          KSRCSP(2,NSRCSP(is),is) = KSRCSP(1,NSRCSP(is),is) + NLIN -1 
        ENDIF 
C----------------------------------------------------------------------C
C     Read a table ... 
C     Read the abcissa value, the ordinate value and their units.
C----------------------------------------------------------------------C
        DO 400 I = 1,NLIN
          READ (IRD,'(A)') CHDUM
          CALL LCASE( CHDUM )
          ISTART = 1
          CALL RDDPR(ISTART,ICOMMA,CHDUM,DUMY(1))
          CALL RDCHR(ISTART,ICOMMA,CHDUM,UNITS)
          CALL RDUNIT(UNITS,DUMY(1))
C----------------------------------------------------------------------C
C     Read a table ...
C     Read the ordinate value and units.
C----------------------------------------------------------------------C
          CALL RDDPR(ISTART,ICOMMA,CHDUM,DUMY(2))
          CALL RDCHR(ISTART,ICOMMA,CHDUM,UNITS)
          CALL RDUNIT(UNITS,DUMY(2))
C----------------------------------------------------------------------C
C     Assign source variables.
C----------------------------------------------------------------------C
          IF( ISTYPE .EQ. 1 ) THEN
            NSTLQ = NSTLQ +1
            SRCLQ(NSTLQ)     = DUMY(2)
            SRCLQ_TM(NSTLQ)  = DUMY(1)
          ELSE
            NSTSP(is) = NSTSP(is) +1
            if( is == 1 ) then  
              SRCSP1(NSTSP(is))     = DUMY(2)
              SRCSP_TM1(NSTSP(is))  = DUMY(1)
            elseif( is == 2 ) then   
              SRCSP2(NSTSP(is))     = DUMY(2)
              SRCSP_TM2(NSTSP(is))  = DUMY(1)
            elseif( is == 3 ) then  
              SRCSP3(NSTSP(is))     = DUMY(2)
              SRCSP_TM3(NSTSP(is))  = DUMY(1)
            endif    
          ENDIF
  400   CONTINUE
C----------------------------------------------------------------------C
C     Write source tables.
C----------------------------------------------------------------------C
        WRITE (IWR,'(A)') ' Source Table '
        IF( ISTYPE .EQ. 1 ) THEN
            DO 500 I = 1,NLIN
              J = KSRCLQ(1,NSRCLQ) + I - 1
              WRITE (IWR,'(A,I4,4(A,1PE11.4))') ' Entry No.',I,
     &          '  Time, s: ',SRCLQ_TM(J),
     &          '  Source Rate, m^3/s: ',SRCLQ(J)
  500       CONTINUE
        ELSE  
          if( is == 1 ) then   
            if( isptp(nsrcsp(is),is) == 1 ) then   
              STR = ' Liq Species Concentration, 1/m^3: '   
            else   
              STR = ' Liq Species Discharge,     m^3/s: '   
            endif    
          elseif( is == 2 ) then   
            if( isptp(nsrcsp(is),is) == 1 ) then  
              STR = ' Pat Species Concentration, 1/m^3: '   
            else   
              STR = ' Pat Species Discharge,     m^3/s: '   
            endif    
          elseif( is == 3 ) then   
            if( isptp(nsrcsp(is),is) == 1 ) then  
              STR = ' Sediment Concentration, kg/m^3: '   
            else   
              STR = ' Sediment Discharge,       kg/s: '   
            endif    
          endif    
          DO 600 I = 1,NLIN
            J = KSRCSP(1,NSRCSP(is),is) +I -1  
            if( is == 1 ) then   
              WRITE (IWR,'(A,I4,2(A,1PE11.4))') ' Entry No.',I,
     &          '  Time, s: ',SRCSP_TM1(J),STR,SRCSP1(J)
            elseif( is == 2 ) then   
              WRITE (IWR,'(A,I4,2(A,1PE11.4))') ' Entry No.',I,
     &          '  Time, s: ',SRCSP_TM2(J),STR,SRCSP2(J)
            elseif( is == 3 ) then  
              WRITE (IWR,'(A,I4,2(A,1PE11.4))') ' Entry No.',I,
     &          '  Time, s: ',SRCSP_TM3(J),STR,SRCSP3(J)
            endif    
  600     CONTINUE
        ENDIF 
C----------------------------------------------------------------------C
  900 CONTINUE
C----------------------------------------------------------------------C
C     End of RDLKSR group.
C----------------------------------------------------------------------C
C
      RETURN
      END

C======================================================================C
C
      SUBROUTINE BCLQ
C
C----------------------------------------------------------------------C
C
C     BCLQ: Boundary Condition for LiQuid transport equation.
C
C----------------------------------------------------------------------C 
      USE DRBCLQ 
      USE DRLINK  
      USE DRNODE  
      USE DRWFLD
      USE FILINX  
      USE FLDINX  
      USE LANDSF
      USE LINEQN  
      USE LNDAMS
      USE SOLVAR
C----------------------------------------------------------------------C
C     Implicit Double Precision.
C----------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C----------------------------------------------------------------------C
C     Type Declarations.
C----------------------------------------------------------------------C
      LOGICAL SCREEN,INITL,RESTART,ZSTART,IPLOT,IRNV
      REAL*8 RU(4),RV(4),RW(4)
      REAL*8, POINTER :: VBC(:) => null()
C----------------------------------------------------------------------C 
      if( .not.associated(VBC) ) then  
        allocate(VBC(nBound)) 
        VBC = 0.d0; 
      endif    
C----------------------------------------------------------------------C
C     Loop over dams. (Temporal implementation.)
C----------------------------------------------------------------------C
      IF( NDAM .GE. 1 ) THEN
        DO 100 i=1,NDAM-1
	    jup  = IUPDAM(i)  ! number of upper dam node.
	    jlow = IDNDAM(i)  ! number of lower dam node. 
          DO K=1,KK(i) 
            IF( LK(jup,K) .LT. 0 ) ivt = LK(jup,K) 
          ENDDO    

C     (predpolagaetsa kol-vo vtek. bran. = 1 ) 
 
          N = ND(ivt,IL(ivt))
          IF( QH(N) .GT. 0 ) THEN 
            BLU(JLOW) = BLU(JLOW) -QH(N)  
CC            QGRNO(jlow)   = Q(2,ivt,IL(ivt))
C            C_GRAN(jlow)  = C(2,ivt,IL(ivt))
C            CS_GRAN(jlow) = CS(2,ivt,IL(ivt))
C            S_GRAN(jlow)  = S(2,ivt,IL(ivt))
          ENDIF
  100   ENDDO   ! (po i, i=1,NDAM)	
	ENDIF
C----------------------------------------------------------------------C
C     Loop over boundary value tables.
C----------------------------------------------------------------------C
      DO 200 NS = 1,NBTU
        NSB = NBCT(1,NS)
        NSE = NBCT(2,NS)
C----------------------------------------------------------------------C
C     Constant source.
C----------------------------------------------------------------------C
        IF( NSB .EQ. NSE ) THEN
          VBC(NS) = BCV(NSB)
C----------------------------------------------------------------------C
C     Tabular source.
C----------------------------------------------------------------------C
        ELSE
          VBC(NS) = 0.D+0
          IF( TIME .GT. BCT(NSB) .AND.
     &        TIME-DT .LT. BCT(NSE)) THEN
            DO 210 I = NSB+1,NSE
              IF( TIME .GT. BCT(I-1) .AND.
     &            TIME-DT .LT. BCT(I)) THEN 
                TMIN = MAX( TIME-DT, BCT(I-1) )
                TMAX = MIN( TIME, BCT(I) )
                DTBND = TMAX -TMIN
                TMID = 0.5*(TMIN +TMAX)
                NLEN = 2
c                CALL LININT( TIME,VALUE,ZZ,BCT(I-1),BCV(I-1),NLEN)
                CALL LININT( TMID,VALUE,ZZ,BCT(I-1:),BCV(I-1:),NLEN) 
C                VALUE = BCV(I-1) 
C                VBC(NS) = VBC(NS) + VALUE
                VBC(NS) = VBC(NS) + VALUE*DTBND/DT
              ENDIF  
  210       CONTINUE    
          ENDIF
        ENDIF
  200 CONTINUE         
C----------------------------------------------------------------------C
C     Shallow Liquid Boundary Conditions.
C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C     Jacobian Matrix Modifier.
C----------------------------------------------------------------------C
      DO 300 NBC = 1,NBCL
        N = IBCND(NBC) 
        IBC = IBCTY(NBC) 
        id = 1 
        if( LK(n,1) > 0 ) id = -1   ! inlet boundary
C----------------------------------------------------------------------C
C     Discharge @ boundary.
C----------------------------------------------------------------------C
        IF( IBC .EQ. 1 ) THEN 
          BLU(N) = BLU(N) +id*VBC(MBCND(NBC)) 
C----------------------------------------------------------------------C
C     Surface @ boundary.
C----------------------------------------------------------------------C
        ELSEIF( IBC .EQ. 2 ) THEN 
          L = LK(n,1)   
          if( L < 0 ) then  
            L = -L  
            np = ND(L,IL(L))  
          else  
            np = ND(L,1)  
          endif    
!          BLU(N) = BLU(N) +ALU(n,n)*(VBC(MBCND(NBC)) -HR(np) -TSF(np))  
          do k=1,Nnode  
            ALU(n,k) = 0.d0  
          enddo   
          ALU(n,n) = 1.d0; BLU(n) = VBC(MBCND(NBC)) -HR(np) -TSF(np)
        ENDIF
  300 CONTINUE
C----------------------------------------------------------------------C
C     End of BCRNLQ group.
C----------------------------------------------------------------------C
C
      RETURN
      END

C======================================================================C
C
      SUBROUTINE RDDAMS
C
C----------------------------------------------------------------------C
C
C     RDDAMS: ReaD DAMS input group.
C
C----------------------------------------------------------------------C 
      USE DRNODE  
      USE DRSOLV 
      USE FILINX  
      USE LNDAMS
C----------------------------------------------------------------------C
C     Implicit Double Precision.
C----------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C----------------------------------------------------------------------C
C     Type Declarations.
C----------------------------------------------------------------------C
      CHARACTER*240 CHDUM
      CHARACTER*80 UNITS,ATYPE,PEROPT
C----------------------------------------------------------------------C
C     Write header to output file.
C----------------------------------------------------------------------C
      WRITE (IWR,'(/A)') ' Dam Parameters'
      WRITE (IWR,'(A )') ' --------------'
C----------------------------------------------------------------------C 
C     Read number of dams. 
C----------------------------------------------------------------------C 
      READ (IRD,*) NDAM
C----------------------------------------------------------------------C 
C     Read Dam Input group. 
C----------------------------------------------------------------------C 
      DO 900 N=1,NDAM  
        READ (IRD,'(A)') CHDUM 
        CALL LCASE(CHDUM) 
        ISTART = 1 
        CALL RDINT(ISTART,ICOMMA,CHDUM,IUPDAM(N))
        CALL RDINT(ISTART,ICOMMA,CHDUM,IDNDAM(N))
        CALL RDINT(ISTART,ICOMMA,CHDUM,NGSIZE(N+1))
C----------------------------------------------------------------------C 
C     Write Dam parameters. 
C----------------------------------------------------------------------C 
        WRITE (IWR,'(/A,I5)') '   Number of Dams                   : ',
     &       NDAM 
        WRITE (IWR,'(/A,I5)') '   Number of Upper Dam Node         : ',
     &       IUPDAM(N)
        WRITE (IWR,'(A,I5)')  '   Number of Lower Dam Node         : ',
     &       IDNDAM(N) 
        WRITE (IWR,'(A,I5)')  '   Number of Points in the Dam Group: ',
     &       NGSIZE(N+1)    
C----------------------------------------------------------------------C
  900 CONTINUE
C----------------------------------------------------------------------C
C     End of RDDAMS group.
C----------------------------------------------------------------------C
C
      RETURN
      END

C======================================================================C
C======================================================================C
C
      SUBROUTINE LATINF
C
C----------------------------------------------------------------------C
C
C     LATINF: LATeral INFlow to water flow equation.
C
C----------------------------------------------------------------------C   
      USE DRLINK
      USE FLDINX
      USE LQSORC  
      USE NUMBRS
      USE POINTS 
      USE SOLVAR
C----------------------------------------------------------------------C
C     Implicit Double Precision.
C----------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
!----------------------------------------------------------------------C 
!     Type Declarations. 
!----------------------------------------------------------------------C  
      INTEGER*4 iRange(4)
C----------------------------------------------------------------------C
C     Loop over sources.
C----------------------------------------------------------------------C
      DO 100 NS = 1,NSRCLQ
        NSB = KSRCLQ(1,NS)
        NSE = KSRCLQ(2,NS)
C----------------------------------------------------------------------C
C     Constant source.
C----------------------------------------------------------------------C
        IF( ISLQV(NS) .EQ. 1 ) THEN
          IF( TIME .GT. SRCLQ_TM(NSB) .AND. 
     &        TIME-DT .LT. SRCLQ_TM(NSE) ) THEN
            TMIN = MAX( TIME-DT, SRCLQ_TM(NSB) )
            TMAX = MIN( TIME, SRCLQ_TM(NSE) )
            DTSRC = TMAX -TMIN 
            SRCX = SRCLQ(NSB)*DTSRC/DT
          ELSE
            SRCX = 0.D+0
          ENDIF
C----------------------------------------------------------------------C
C     Tabular source.
C----------------------------------------------------------------------C
        ELSEIF( ISLQV(NS) .EQ. 2 ) THEN
          IF( TIME .GT. SRCLQ_TM(NSB) .AND.
     &        TIME-DT .LT. SRCLQ_TM(NSE)) THEN
            SRCX = 0.D+0  
!            HR = 0.D+0 
!            PART = 0.D+0 
            DO 140 I = NSB+1,NSE
              IF( TIME .GT. SRCLQ_TM(I-1) .AND.
     &            TIME-DT .LT. SRCLQ_TM(I)) THEN
                TMIN = MAX( TIME-DT, SRCLQ_TM(I-1) )
                TMAX = MIN( TIME, SRCLQ_TM(I) )
                DTSRC = TMAX -TMIN
                TMID = 0.5*(TMIN +TMAX)
                NLEN = 2  
!                CALL LININT( TMID,SRCY,ZZ,SRCLQ_TM(I-1),SRCLQ(I-1),NLEN) 
                SRCY = SRCLQ(I-1)
                SRCX = SRCX + SRCY*DTSRC/DT 
              ENDIF
  140       CONTINUE   
          ELSE
            SRCX = 0.D+0
          ENDIF
        ENDIF
C----------------------------------------------------------------------C
C     Loop over source domain.
C----------------------------------------------------------------------C  
        iRange(1) = ISRCLQ(1,NS); iRange(2) = ISRCLQ(2,NS);
        if( iRange(1) == 0 ) iRange(1) = 1; 
        if( iRange(2) == 0 ) iRange(2) = nLink; 
        Link1 = max(1,iRange(1)); Link2 = min(nLink,iRange(2));  
        iRange(3) = ISRCLQ(3,NS); iRange(4) = ISRCLQ(4,NS);
        
        DO 110 L = Link1,Link2
          if( iRange(3) == 0 ) then 
            iP1 = 1;  
          else  
            iP1 = max(1,iRange(3));   
          endif    
          if( iRange(4) == 0 ) then 
            iP2 = IL(L); 
          else  
            iP2 = min(IL(L),iRange(4)); 
          endif 
          DO 120 I = iP1, iP2
            N = ND(L,I)  
            FLQ(N) = SRCX
  120     CONTINUE
  110   CONTINUE
  100 CONTINUE
C----------------------------------------------------------------------C
C     End of LATINF group.
C----------------------------------------------------------------------C
C
      RETURN
      END

C======================================================================C
C
      SUBROUTINE SPLATR (is)
C
C----------------------------------------------------------------------C
C
C     SPLATR: SPecies sources in LATeRal inflow.
C
C     is = 1 - liquid species source; 
C     is = 2 - particulate species source;
C     is = 3 - sediment source.
C----------------------------------------------------------------------C
      USE DRLINK
      USE FLDINX
      USE LQSORC  
      USE LINEQN
      USE NUMBRS
      USE POINTS 
      USE SOLVAR 
      USE SPSOUR
C----------------------------------------------------------------------C
C     Implicit Double Precision.
C----------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
!----------------------------------------------------------------------C 
!     Type Declarations. 
!----------------------------------------------------------------------C  
      INTEGER*4 iRange(4) 
      Real*8 FSP(nPoint)
C----------------------------------------------------------------------C
C     Loop over species sources related to lateral inflow. 
C----------------------------------------------------------------------C
      DO 100 NS = 1,NSRCLQ
        NSB = KSRCLQ(1,NS)
        NSE = KSRCLQ(2,NS)
C----------------------------------------------------------------------C
C     Constant source.
C----------------------------------------------------------------------C
        IF( ISLQV(NS) .EQ. 1 ) THEN
          IF( TIME .GT. SRCLQ_TM(NSB) .AND. 
     &        TIME-DT .LT. SRCLQ_TM(NSE) ) THEN
            TMIN = MAX( TIME-DT, SRCLQ_TM(NSB) )
            TMAX = MIN( TIME, SRCLQ_TM(NSE) )
            DTSRC = TMAX -TMIN 
            SRCW = SRCLQ(NSB)*DTSRC/DT
          ELSE
            SRCW = 0.D+0
          ENDIF
C----------------------------------------------------------------------C
C     Tabular source.
C----------------------------------------------------------------------C
        ELSEIF( ISLQV(NS) .EQ. 2 ) THEN
          IF( TIME .GT. SRCLQ_TM(NSB) .AND.
     &        TIME-DT .LT. SRCLQ_TM(NSE)) THEN
            SRCW = 0.D+0
            DO 140 I = NSB+1,NSE
              IF( TIME .GT. SRCLQ_TM(I-1) .AND.
     &            TIME-DT .LT. SRCLQ_TM(I)) THEN
                TMIN = MAX( TIME-DT, SRCLQ_TM(I-1) )
                TMAX = MIN( TIME, SRCLQ_TM(I) )
                DTSRC = TMAX -TMIN
                TMID = 0.5*(TMIN +TMAX)
                NLEN = 2
!                CALL LININT( TMID,SRCY,ZZ,SRCLQ_TM(I-1),SRCLQ(I-1),NLEN)
!                CALL LININT( TMID,SRCP,ZZ,SRCLQ_TM(I-1),SPCLQ(I-1),NLEN)
                CALL LININT( TMID,SRCY,ZZ,SRCLQ_TM,SRCLQ,NLEN)
!                CALL LININT( TMID,SRCP,ZZ,SRCLQ_TM,SPCLQ,NLEN)
                SRCW = SRCW + SRCY*DTSRC/DT 
              ENDIF
  140       CONTINUE   
          ELSE
            SRCW = 0.D+0
          ENDIF
        ENDIF
C----------------------------------------------------------------------C
C     Loop over infiltration/leakage source domain.
C----------------------------------------------------------------------C
        iRange(1) = ISRCLQ(1,NS); iRange(2) = ISRCLQ(2,NS);
        if( iRange(1) == 0 ) iRange(1) = 1; 
        if( iRange(2) == 0 ) iRange(2) = nLink; 
        Link1 = max(1,iRange(1)); Link2 = min(nLink,iRange(2));  
        iRange(3) = ISRCLQ(3,NS); iRange(4) = ISRCLQ(4,NS);
        
        DO 110 L = Link1,Link2
          if( iRange(3) == 0 ) then 
            iP1 = 1;  
          else  
            iP1 = max(1,iRange(3));   
          endif    
          if( iRange(4) == 0 ) then 
            iP2 = IL(L); 
          else  
            iP2 = min(IL(L),iRange(4)); 
          endif 
          DO 120 I = iP1, iP2
            N = ND(L,I)  
            FLQ(N) = SRCW
  120     CONTINUE
  110   CONTINUE
  100 CONTINUE
C----------------------------------------------------------------------C
C     Loop over species sources.
C----------------------------------------------------------------------C
      DO 500 NS = 1,NSRCSP(is)
        NSB = KSRCSP(1,NS,is)
        NSE = KSRCSP(2,NS,is)
        SRCX = 0.D+0  
C----------------------------------------------------------------------C
        if( is == 1 ) then  
          tBeg = SRCSP_TM1(NSB);  
          tEnd = SRCSP_TM1(NSE);   
          tSr = SRCSP1(NSB)
        elseif( is == 2 ) then  
          tBeg = SRCSP_TM2(NSB);  
          tEnd = SRCSP_TM2(NSE);  
          tSr = SRCSP2(NSB)
        elseif( is == 3 ) then  
          tBeg = SRCSP_TM3(NSB);  
          tEnd = SRCSP_TM3(NSE);  
          tSr = SRCSP3(NSB)
        endif  
C----------------------------------------------------------------------C
C     Constant source.
C----------------------------------------------------------------------C
        IF( ISSPV(NS,is) .EQ. 1 ) THEN 
          IF( TIME .GT. tBeg .AND. 
     &        TIME-DT .LT. tEnd ) THEN
            TMIN = MAX( TIME-DT, tBeg )
            TMAX = MIN( TIME, tEnd )
            DTSRC = TMAX -TMIN
            SRCX = tSr*DTSRC/DT
          ENDIF  
C----------------------------------------------------------------------C
C     Tabular source.
C----------------------------------------------------------------------C
        ELSEIF( ISSPV(NS,is) .EQ. 2 ) THEN
          IF( TIME .GT. tBeg .AND. 
     &        TIME-DT .LT. tEnd ) THEN
            DO 400 I = NSB+1,NSE 
              if( is == 1 ) then  
                tBeg = SRCSP_TM1(i-1);  
                tEnd = SRCSP_TM1(i);   
                tSr = SRCSP1(i-1)
              elseif( is == 2 ) then  
                tBeg = SRCSP_TM2(i-1);  
                tEnd = SRCSP_TM2(i);  
                tSr = SRCSP2(i-1)
              elseif( is == 3 ) then  
                tBeg = SRCSP_TM3(i-1);  
                tEnd = SRCSP_TM3(i);  
                tSr = SRCSP3(i-1)
              endif  
                
              IF( TIME .GT. tBeg .AND. 
     &            TIME-DT .LT. tEnd ) THEN
                TMIN = MAX( TIME-DT, tBeg )
                TMAX = MIN( TIME, tEnd )
                DTSRC = TMAX -TMIN
                TMID = 0.5*(TMIN +TMAX)
                NLEN = 2
!                CALL LININT( TMID,SRCY,ZZ,SRCSP_TM(I-1),SRCSP(I-1),NLEN)
                SRCX = SRCX + tSr*DTSRC/DT
              ENDIF  
  400       CONTINUE   
          ENDIF
        ENDIF  
C----------------------------------------------------------------------C
C     Loop over source domain.
C----------------------------------------------------------------------C
        iRange(1) = ISRCSP(1,NS,is); iRange(2) = ISRCSP(2,NS,is);
        if( iRange(1) == 0 ) iRange(1) = 1; 
        if( iRange(2) == 0 ) iRange(2) = nLink; 
        Link1 = max(1,iRange(1)); Link2 = min(nLink,iRange(2));  
        iRange(3) = ISRCSP(3,NS,is); iRange(4) = ISRCSP(4,NS,is); 
!        
        FSP = 0.d0;
        DO 330 L = Link1,Link2
          if( iRange(3) == 0 ) then 
            iP1 = 1;  
          else  
            iP1 = max(1,iRange(3));   
          endif    
          if( iRange(4) == 0 ) then 
            iP2 = IL(L); 
          else  
            iP2 = min(IL(L),iRange(4)); 
          endif 
          DO 320 I = iP1, iP2
            N = ND(L,I)  
            FSP(N) = SRCX
  320     CONTINUE
  330   CONTINUE
C----------------------------------------------------------------------C
C     Incorporate Species sources into transport equation.
C----------------------------------------------------------------------C  
      do L=1, NLINK  
        do i=1,IL(L) 
          N = ND(L,I)  
C----------------------------------------------------------------------C
C     Species concentration source.
C----------------------------------------------------------------------C 
          NEL = iEq(N)
          IF( ISPTP(NS,is) == 1 ) THEN  
            ALC(NEL) = ALC(NEL) +FLQ(N)  
            BLC(N) = BLC(N) + FSP(N)*FLQ(N)
C----------------------------------------------------------------------C
C     Species source.
C----------------------------------------------------------------------C
          ELSEIF( ISPTP(NS,is) == 2 ) THEN
            ALC(NEL) = ALC(NEL) +FLQ(N)  
            BLC(N) = BLC(N) + FSP(N)
          ENDIF 
        enddo  
      enddo  
  500 CONTINUE
C----------------------------------------------------------------------C
C     End of SPLATR group.
C----------------------------------------------------------------------C
C
      RETURN
      END

C======================================================================C
C======================================================================C
C
      SUBROUTINE RDNDSR
C
C----------------------------------------------------------------------C
C
C     RDNDSR: ReaD NoDe SouRces and sinks input group.  
C
C----------------------------------------------------------------------C  
      USE FILINX   
      use FLDINX
      USE LOGCLS
      USE LQNODE 
      USE SPNODE 
      use intersub
C----------------------------------------------------------------------C
C     Implicit Double Precision.
C----------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C----------------------------------------------------------------------C
      INCLUDE 'utils.inc' ! subroutine interfaces
C----------------------------------------------------------------------C
C     Type Declarations.
C----------------------------------------------------------------------C
      CHARACTER*240 CHDUM
      CHARACTER*80 ADUM,UNITS,TYPE
      CHARACTER*20 TABLE  
      CHARACTER*20, POINTER :: TBNAME(:)
      REAL*8 DUMY(6),VALUE
      CHARACTER STR*35
      integer*4 NSTSP(3)  
      
      REAL*8, POINTER, SAVE :: tempArray(:)
      
C----------------------------------------------------------------------C 
C     For i-th boundary node (1 <= i <= NBCL):
C        IBCND(i)  - number of i-th node;
C        IBCTY(i)  - type of boundary condition;
C        MBCND(i)  - table number of boundary values; 
C     For j-th table of boundary values:
C        NBCT(1,j) - first line in the table;
C        NBCT(2,j) - last line in the table.
C     For k-th line (1 <= k <= NBTU) in table of boundary values:
C        BCV(k)  -  boundary value at
C        BCT(k)  -  time moment.
C----------------------------------------------------------------------C
C     Write header to output file.
C----------------------------------------------------------------------C
      WRITE (IWR,'(/A)') ' Node Sources and Sinks'
      WRITE (IWR,'(A )') ' ----------------------'
C----------------------------------------------------------------------C
C     Set table counters to zero.
C----------------------------------------------------------------------C
      NDSRLQ = 0
      NSTLQ  = 0
      NDSRSP = 0
      NSTSP  = 0
C----------------------------------------------------------------------C
C     Found 'Sources and Sinks' data
C----------------------------------------------------------------------C
      READ (IRD,*) NSZN !NSZN - number of sources and sinks data
C----------------------------------------------------------------------C
C     Allocate sources data
C----------------------------------------------------------------------C 
      Allocate(INDLQV(NSZN),INDSPV(NSZN,3),INDTP(NSZN,3)) 
      INDLQV = 0; INDSPV = 0; INDTP = 0; 
      Allocate(INDSLQ(2,NSZN),KNDSLQ(2,NSZN))  
      INDSLQ = 0; KNDSLQ = 0; 
      Allocate(SNDLQ(2*NSZN),SNDLQ_TM(2*NSZN)); 
      SNDLQ = 0.d0; SNDLQ_TM = 0.d0;   
      
      Allocate(INDSSP(2,NSZN,3),KNDSSP(2,NSZN,3))  
      INDSSP = 0; KNDSSP = 0; 
      Allocate(SNDSP1(2*NSZN),SNDSP_TM1(2*NSZN))  
      SNDSP1 = 0.d0; SNDSP_TM1 = 0.d0;
      Allocate(SNDSP2(2*NSZN),SNDSP_TM2(2*NSZN))  
      SNDSP2 = 0.d0; SNDSP_TM2 = 0.d0;
      Allocate(SNDSP3(2*NSZN),SNDSP_TM3(2*NSZN))  
      SNDSP3 = 0.d0; SNDSP_TM3 = 0.d0;
C----------------------------------------------------------------------C
C     Loop over the number of sources and/or sinks.
C----------------------------------------------------------------------C
      DO 900 N = 1,NSZN
        READ (IRD,'(A)') CHDUM
        CALL LCASE( CHDUM )
C----------------------------------------------------------------------C
C     Read source type.
C----------------------------------------------------------------------C
        ITAB = 0
        ISTYPE = 2
        ISTART = 1
        CALL RDCHR(ISTART,ICOMMA,CHDUM,ADUM)
        IF( ADUM(1:11) .EQ. 'liquid disc') THEN  
          LQSRND = .true.  
          ISTYPE = 1
          NDSRLQ = NDSRLQ + 1
        ELSEIF( ADUM(1:16) .EQ. 'species liq conc' ) THEN
          CLSRND = .true.  
          NDSRSP(1) = NDSRSP(1) + 1
          INDTP(NDSRSP(1),1) = 1 
          is = 1
        ELSEIF( ADUM(1:16) .EQ. 'species liq disc' ) THEN
          CLSRND = .true.  
          NDSRSP(1) = NDSRSP(1) + 1
          INDTP(NDSRSP(1),1) = 2 
          is = 1
        ELSEIF( ADUM(1:16) .EQ. 'species pat conc' ) THEN
          CPSRND = .true.  
          NDSRSP(2) = NDSRSP(2) + 1
          INDTP(NDSRSP(2),2) = 1  
          is = 2
        ELSEIF( ADUM(1:16) .EQ. 'species pat disc' ) THEN
          CPSRND = .true.  
          NDSRSP(2) = NDSRSP(2) + 1
          INDTP(NDSRSP(2),2) = 2  
          is = 2
        ELSEIF( ADUM(1:13) .EQ. 'sediment conc' ) THEN
          SDSRND = .true.  
          NDSRSP(3) = NDSRSP(3) + 1
          INDTP(NDSRSP(3),3) = 1   
          is = 3
        ELSEIF( ADUM(1:13) .EQ. 'sediment disc' ) THEN
          SDSRND = .true.  
          NDSRSP(3) = NDSRSP(3) + 1
          INDTP(NDSRSP(3),3) = 2   
          is = 3
        ENDIF
C----------------------------------------------------------------------C
C     Read source temporal variation.
C----------------------------------------------------------------------C
        CALL RDCHR(ISTART,ICOMMA,CHDUM,TYPE)
        IF( TYPE(1:8) .EQ. 'constant') THEN
          IF( ISTYPE .EQ. 1 ) THEN
            INDLQV(NDSRLQ) = 1             
          ELSE
            INDSPV(NDSRSP(is),is) = 1
          ENDIF
C----------------------------------------------------------------------C
C     Read 3 parameter and unit combinations.
C----------------------------------------------------------------------C
          DO 100 I = 1, 3
            CALL RDDPR(ISTART,ICOMMA,CHDUM,DUMY(I))
            CALL RDCHR(ISTART,ICOMMA,CHDUM,UNITS)
            CALL RDUNIT(UNITS,DUMY(I))
  100     CONTINUE
C----------------------------------------------------------------------C
C     Assign source variables.
C----------------------------------------------------------------------C
          IF( ISTYPE .EQ. 1 ) THEN 
            NSTLQ = NSTLQ + 1 
            KNDSLQ(1,NDSRLQ) = NSTLQ
            KNDSLQ(2,NDSRLQ) = NSTLQ+1
            SNDLQ(NSTLQ)    = DUMY(1)
            SNDLQ_TM(NSTLQ) = DUMY(2)
            NSTLQ = NSTLQ + 1
            SNDLQ_TM(NSTLQ) = DUMY(3)
          ELSE
            NSTSP(is) = NSTSP(is) + 1 
            KNDSSP(1,NDSRSP(is),is) = NSTSP(is)
            KNDSSP(2,NDSRSP(is),is) = NSTSP(is) +1 
            if( is == 1 ) then   
              SNDSP1(NSTSP(is))    = DUMY(1)
              SNDSP_TM1(NSTSP(is)) = DUMY(2)
              NSTSP(is) = NSTSP(is) + 1
              SNDSP_TM1(NSTSP(is)) = DUMY(3)
            elseif( is == 2 ) then   
              SNDSP2(NSTSP(is))    = DUMY(1)
              SNDSP_TM2(NSTSP(is)) = DUMY(2)
              NSTSP(is) = NSTSP(is) + 1
              SNDSP_TM2(NSTSP(is)) = DUMY(3)
            elseif( is == 3 ) then
              SNDSP3(NSTSP(is))    = DUMY(1)
              SNDSP_TM3(NSTSP(is)) = DUMY(2)
              NSTSP(is) = NSTSP(is) + 1
              SNDSP_TM3(NSTSP(is)) = DUMY(3)
            endif
            
          ENDIF
       ELSEIF( TYPE(1:7) .EQ. 'tabular' ) THEN
          ITAB = 1
          IF( ISTYPE .EQ. 1 ) THEN
            INDLQV(NDSRLQ) = 2             
          ELSE
            INDSPV(NDSRSP(is),is) = 2
          ENDIF
C----------------------------------------------------------------------C
C     Read table name.
C----------------------------------------------------------------------C
          CALL RDCHR(ISTART,ICOMMA,CHDUM,TABLE)
C----------------------------------------------------------------------C
C     Read 5 blank data entries.
C----------------------------------------------------------------------C
          DO 200 I = 1,5
            CALL RDCHR(ISTART,ICOMMA,CHDUM,ADUM)
  200     CONTINUE
        ENDIF
C----------------------------------------------------------------------C
C     Read 2 NODAL integers (range of source/sink)
C     The source domain KDSRLQ and KDSRSP refer to the NODE
C     directly.
C----------------------------------------------------------------------C
        DO 300 I = 1,2
          IF( ISTYPE .EQ. 1 ) THEN
            CALL RDINT(ISTART,ICOMMA,CHDUM,INDSLQ(I,NDSRLQ))
          ELSE
            CALL RDINT(ISTART,ICOMMA,CHDUM,INDSSP(I,NDSRSP(is),is))
          ENDIF  
  300   CONTINUE
C----------------------------------------------------------------------C
C     Write source and sink data to output file.
C----------------------------------------------------------------------C
        IF( N .NE. 1 ) WRITE (IWR,'(/)')
C----------------------------------------------------------------------C
C     Write liquid sources.
C----------------------------------------------------------------------C
        IF( ISTYPE .EQ. 1 ) THEN
          WRITE (IWR,'(A)') ' Liquid Discharge Source'
          IF( INDLQV(NDSRLQ) .EQ. 1 ) THEN
            WRITE (IWR,'(A)') ' Temporal Variation: Constant'
            WRITE (IWR,'(A,1PE11.4)') ' Source Rate,          m^3/s: ',
     &        SNDLQ(KNDSLQ(1,NDSRLQ))
            WRITE (IWR,'(A,1PE11.4)') ' Source Start Time,        s: ',
     &        SNDLQ_TM(KNDSLQ(1,NDSRLQ))
            WRITE (IWR,'(A,1PE11.4)') ' Source Stop Time,         s: ',
     &        SNDLQ_TM(KNDSLQ(2,NDSRLQ))
          ELSEIF( INDLQV(NDSRLQ) .EQ. 2 ) THEN
            WRITE (IWR,'(A)') ' Temporal Variation: Tabular'
            WRITE (IWR,'(2A)') ' Table Name        : ',TABLE
          ENDIF
C----------------------------------------------------------------------C
C     Write species sources.
C----------------------------------------------------------------------C
        ELSE  
          if( is == 1 ) then  
            if( iNDtp(nDsrsp(is),is) == 1 ) then  
              WRITE (IWR,'(A)') ' Species Liquid Concentration Source'  `
              STR = ' Species Concentration,   1/m^3: '  
            else   
              WRITE (IWR,'(A)') ' Species Liquid Discharge Source' 
              STR = ' Species Discharge,         1/s: '  
            endif     
          elseif( is == 2 ) then  
            if( iNDtp(nDsrsp(is),is) == 1 ) then  
              WRITE (IWR,'(A)')   
     &          ' Species Particulate Concentration Source' 
              STR = ' Species Concentration,   1/m^3: '  
            else   
              WRITE (IWR,'(A)') ' Species Particulate Discharge Source' 
              STR = ' Species Discharge,         1/s: '  
            endif     
          elseif( is == 3 ) then  
            if( iNDtp(nDsrsp(is),is) == 1 ) then  
              WRITE (IWR,'(A)') ' Sediment Concentration Source' 
              STR = ' Sediment Concentration, kg/m^3: '  
            else   
              WRITE (IWR,'(A)') ' Sediment Discharge Source' 
              STR = ' Sediment Discharge,       kg/s: '  
            endif     
          endif      
            
          IF( INDSPV(NDSRSP(is),is) == 1 ) THEN
            WRITE (IWR,'(A)') ' Temporal Variation: Constant'
            WRITE (IWR,'(A,1PE11.4)') STR,DUMY(1)
!     &        SRCSP(KSRCSP(1,NSRCSP(is),is),is)
            WRITE (IWR,'(A,1PE11.4)') ' Source Start Time,        s: ',
     &        DUMY(2)
            WRITE (IWR,'(A,1PE11.4)') ' Source Stop Time,         s: ',
     &        DUMY(3)
          ELSEIF( INDSPV(NDSRSP(is),is) .EQ. 2 ) THEN
            WRITE (IWR,'(A)') ' Temporal Variation: Tabular'
          ENDIF
        ENDIF  
C----------------------------------------------------------------------C
C     Write source domain information.
C----------------------------------------------------------------------C
        IF( ISTYPE .EQ. 1 ) THEN
          WRITE (IWR,'(A,I4,A,I4)') '   Node : N = ',INDSLQ(1,NDSRLQ),
     &      ' to ',INDSLQ(2,NDSRLQ)
        ELSE
          WRITE (IWR,'(A,I4,A,I4)') '   Node : N = ',
     &      INDSSP(1,NDSRSP(is),is),' to ',INDSSP(2,NDSRSP(is),is)
        ENDIF
C----------------------------------------------------------------------C
C     Read source tables.
C----------------------------------------------------------------------C
        IF( ITAB .EQ. 0 ) GOTO 900 
C----------------------------------------------------------------------C
C     Length of table.
C----------------------------------------------------------------------C
        READ (IRD,*) NLIN
!------------------------------------------------------------------------------C
!     REALLOCATE  BCVS,BCTS arrays.
!------------------------------------------------------------------------------C
        IF( ISTYPE .EQ. 1 ) THEN  
          nk = NSTLQ +NLIN;
          CALL ReallocateR8(NSTLQ,nk,SNDLQ)
          call ReallocateR8(NSTLQ,nk,SNDLQ_TM)
        ELSE   
          if( is == 1 ) then 
            call ReallocateR8(NSTSP(is),nk,SNDSP1)
            call ReallocateR8(NSTSP(is),nk,SNDSP_TM1)  
          elseif( is == 2 ) then  
            call ReallocateR8(NSTSP(is),nk,SNDSP2)
            call ReallocateR8(NSTSP(is),nk,SNDSP_TM2)  
          elseif( is == 3 ) then   
            call ReallocateR8(NSTSP(is),nk,SNDSP3)
            call ReallocateR8(NSTSP(is),nk,SNDSP_TM3)  
          endif    
        ENDIF    
C----------------------------------------------------------------------C
C     Table start index.
C----------------------------------------------------------------------C 
        IF( ISTYPE .EQ. 1 ) THEN
          KNDSLQ(1,NDSRLQ) = NSTLQ +1
          KNDSLQ(2,NDSRLQ) = KNDSLQ(1,NDSRLQ) + NLIN -1 
        ELSE
          KNDSSP(1,NDSRSP(is),is) = NSTSP(is) +1
          KNDSSP(2,NDSRSP(is),is) = KNDSSP(1,NDSRSP(is),is) + NLIN -1 
        ENDIF 
C----------------------------------------------------------------------C
C     Read a table ... 
C     Read the abcissa value, the ordinate(next block) value and their units.
C----------------------------------------------------------------------C
        DO 400 I = 1,NLIN
          READ (IRD,'(A)') CHDUM
          CALL LCASE( CHDUM )
          ISTART = 1
          CALL RDDPR(ISTART,ICOMMA,CHDUM,DUMY(1))
          CALL RDCHR(ISTART,ICOMMA,CHDUM,UNITS)
          CALL RDUNIT(UNITS,DUMY(1))
C----------------------------------------------------------------------C
C     Read a table ...
C     Read the ordinate value and units.
C----------------------------------------------------------------------C
          CALL RDDPR(ISTART,ICOMMA,CHDUM,DUMY(2))
          CALL RDCHR(ISTART,ICOMMA,CHDUM,UNITS)
          CALL RDUNIT(UNITS,DUMY(2))
C----------------------------------------------------------------------C
C     Assign source variables.
C----------------------------------------------------------------------C
          IF( ISTYPE .EQ. 1 ) THEN
            NSTLQ = NSTLQ +1
            SNDLQ(NSTLQ)     = DUMY(2)
            SNDLQ_TM(NSTLQ)  = DUMY(1)
          ELSE
            NSTSP(is) = NSTSP(is) +1
            if( is == 1 ) then  
              SNDSP1(NSTSP(is))     = DUMY(2)
              SNDSP_TM1(NSTSP(is))  = DUMY(1)
            elseif( is == 2 ) then   
              SNDSP2(NSTSP(is))     = DUMY(2)
              SNDSP_TM2(NSTSP(is))  = DUMY(1)
            elseif( is == 3 ) then  
              SNDSP3(NSTSP(is))     = DUMY(2)
              SNDSP_TM3(NSTSP(is))  = DUMY(1)
            endif    
          ENDIF
  400   CONTINUE
C----------------------------------------------------------------------C
C     Write source tables.
C----------------------------------------------------------------------C
        WRITE (IWR,'(A)') ' Nodal Source Table '
        IF( ISTYPE .EQ. 1 ) THEN
            DO 500 I = 1,NLIN
              J = KNDSLQ(1,NDSRLQ) + I - 1
              WRITE (IWR,'(A,I4,4(A,1PE11.4))') ' Entry No.',I,
     &          '  Time, s: ',SNDLQ_TM(J),
     &          '  Source Rate, m^3/s: ',SNDLQ(J)
  500       CONTINUE
        ELSE  
          if( is == 1 ) then   
            if( indtp(ndsrsp(is),is) == 1 ) then   
              STR = ' Liq Species Concentration, 1/m^3: '   
            else   
              STR = ' Liq Species Discharge,     m^3/s: '   
            endif    
          elseif( is == 2 ) then   
            if( indtp(ndsrsp(is),is) == 1 ) then  
              STR = ' Pat Species Concentration, 1/m^3: '   
            else   
              STR = ' Pat Species Discharge,     m^3/s: '   
            endif    
          elseif( is == 3 ) then   
            if( indtp(ndsrsp(is),is) == 1 ) then  
              STR = ' Sediment Concentration, kg/m^3: '   
            else   
              STR = ' Sediment Discharge,       kg/s: '   
            endif    
          endif    
          DO 600 I = 1,NLIN
            J = KNDSSP(1,NDSRSP(is),is) +I -1  
            if( is == 1 ) then   
              WRITE (IWR,'(A,I4,2(A,1PE11.4))') ' Entry No.',I,
     &          '  Time, s: ',SNDSP_TM1(J),STR,SNDSP1(J)
            elseif( is == 2 ) then   
              WRITE (IWR,'(A,I4,2(A,1PE11.4))') ' Entry No.',I,
     &          '  Time, s: ',SNDSP_TM2(J),STR,SNDSP2(J)
            elseif( is == 3 ) then  
              WRITE (IWR,'(A,I4,2(A,1PE11.4))') ' Entry No.',I,
     &          '  Time, s: ',SNDSP_TM3(J),STR,SNDSP3(J)
            endif    
  600     CONTINUE
        ENDIF 
C----------------------------------------------------------------------C
  900 CONTINUE  
C----------------------------------------------------------------------C
C     Reallocate memory.
C----------------------------------------------------------------------C
      if( LQSRND ) then  
        call ReallocateI4 ( NSZN,NDSRLQ,INDLQV ); 
        call Reallocate2DI4 ( 2,NSZN,2,NDSRLQ,INDSLQ );
        call Reallocate2DI4 ( 2,NSZN,2,NDSRLQ,KNDSLQ );
        if( NSTLQ < 2*NSZN ) then   
           call ReallocateR8( 2*NSZN,NSTLQ,SNDLQ ); 
           call ReallocateR8( 2*NSZN,NSTLQ,SNDLQ_TM ); 
        endif   
        Allocate (SRCndLQ(NNode)); SRCndLQ = 0.0d0;
      else   
        deallocate( INDLQV,KNDSLQ,INDSLQ,SNDLQ,SNDLQ_TM );    
      endif   
      
      NNsp = max( NDSRSP(1),NDSRSP(2),NDSRSP(3) );
      if( SDSRND .or. CLSRND .or. CPSRND ) then  
        call Reallocate2DI4 ( NSZN,3,NNsp,3,INDSPV );   
        call Reallocate2DI4 ( NSZN,3,NNsp,3,INDTP );    
        
        call Reallocate3DI4 ( 2,NSZN,3,2,NNsp,3,INDSSP );
        call Reallocate3DI4 ( 2,NSZN,3,2,NNsp,3,KNDSSP );
        Allocate (SRCndSP(NNode),SPndConc(NNode)); 
        SRCndSP = 0.0d0; SPndConc = 0.0d0;  
        if( .not. LQSRND ) then  
          Allocate (SRCndLQ(NNode)); 
          SRCndLQ = 0.0d0;   
        endif  
        
!        if( .not. LQSRND ) then   
!          do i=1,NNsp  
!             do j=1,3    
!                if( INDTP(i,j) == 1 ) then    
!                   write(*,*) 
!     &  'Warning: Species node source type is concentration type but',
!     &  'source water discharge at the node equals to zero.',
!     &  'Node number is ',INDSSP(1,NDSRSP(j),j)              
!                   write(IWR,*) 
!     &  'Warning: Species node source type is concentration type but',
!     &  'source water discharge at the node equals to zero.',
!     &  'Node number is ',INDSSP(1,NDSRSP(j),j)              
!                endif    
!             enddo    
!          enddo    
!        endif    
      else   
        deallocate( INDSPV,INDTP,INDSSP,KNDSSP );    
      endif  
      
      if( SDSRND ) then  
        if( NSTSP(3) < 2*NSZN ) then   
           call ReallocateR8( 2*NSZN,NSTSP(3),SNDSP3 ); 
           call ReallocateR8( 2*NSZN,NSTSP(3),SNDSP_TM3 ); 
        endif  
      else  
        deallocate( SNDSP3,SNDSP_TM3 );  
      endif  
      
      if( CLSRND ) then  
        if( NSTSP(1) < 2*NSZN ) then   
           call ReallocateR8( 2*NSZN,NSTSP(1),SNDSP1 ); 
           call ReallocateR8( 2*NSZN,NSTSP(1),SNDSP_TM1 ); 
        endif  
      else  
        deallocate( SNDSP1,SNDSP_TM1 );  
      endif    
      
      if( CPSRND ) then  
        if( NSTSP(2) < 2*NSZN ) then   
           call ReallocateR8( 2*NSZN,NSTSP(2),SNDSP2 ); 
           call ReallocateR8( 2*NSZN,NSTSP(2),SNDSP_TM2 ); 
        endif  
      else  
        deallocate( SNDSP2,SNDSP_TM2 );  
      endif    
C----------------------------------------------------------------------C
C     End of RDNDSR group.
C----------------------------------------------------------------------C
C
      RETURN
      END

C======================================================================C 
C
      SUBROUTINE NDLQSR
C
C----------------------------------------------------------------------C
C
C     NDLQSR: inserting NoDal LiQuid sources or sinks to water flow 
!             equation.
C
C----------------------------------------------------------------------C   
      USE DRLINK
      USE FLDINX
      USE LQNODE 
      USE NUMBRS
      USE POINTS 
      USE SOLVAR
      USE LINEQN 
C----------------------------------------------------------------------C
C     Implicit Double Precision.
C----------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
!----------------------------------------------------------------------C 
!     Type Declarations. 
!----------------------------------------------------------------------C  
      INTEGER*4 iRange(2)
C----------------------------------------------------------------------C
C     Loop over sources.
C----------------------------------------------------------------------C
      DO 100 NS = 1,NDSRLQ
        NSB = KNDSLQ(1,NS)
        NSE = KNDSLQ(2,NS)
C----------------------------------------------------------------------C
C     Constant source.
C----------------------------------------------------------------------C
        IF( INDLQV(NS) .EQ. 1 ) THEN
          IF( TIME .GT. SNDLQ_TM(NSB) .AND. 
     &        TIME-DT .LT. SNDLQ_TM(NSE) ) THEN
            TMIN = MAX( TIME-DT, SNDLQ_TM(NSB) )
            TMAX = MIN( TIME, SNDLQ_TM(NSE) )
            DTSRC = TMAX -TMIN 
            SRCX = SNDLQ(NSB)*DTSRC/DT
          ELSE
            SRCX = 0.D+0
          ENDIF
C----------------------------------------------------------------------C
C     Tabular source.
C----------------------------------------------------------------------C
        ELSEIF( INDLQV(NS) .EQ. 2 ) THEN
          IF( TIME .GT. SNDLQ_TM(NSB) .AND.
     &        TIME-DT .LT. SNDLQ_TM(NSE)) THEN
            SRCX = 0.D+0  
!            HR = 0.D+0 
!            PART = 0.D+0 
            DO 140 I = NSB+1,NSE
              IF( TIME .GT. SNDLQ_TM(I-1) .AND.
     &            TIME-DT .LT. SNDLQ_TM(I)) THEN
                TMIN = MAX( TIME-DT, SNDLQ_TM(I-1) )
                TMAX = MIN( TIME, SNDLQ_TM(I) )
                DTSRC = TMAX -TMIN
                TMID = 0.5*(TMIN +TMAX)
                NLEN = 2  
!                CALL LININT( TMID,SRCY,ZZ,SRCLQ_TM(I-1),SRCLQ(I-1),NLEN) 
                SRCY = SNDLQ(I-1)
                SRCX = SRCX + SRCY*DTSRC/DT 
              ENDIF
  140       CONTINUE   
          ELSE
            SRCX = 0.D+0
          ENDIF
        ENDIF
C----------------------------------------------------------------------C
C     Loop over source domain.
C----------------------------------------------------------------------C  
        iRange(1) = INDSLQ(1,NS); iRange(2) = INDSLQ(2,NS);
        if( iRange(1) == 0 ) iRange(1) = 1; 
        if( iRange(2) == 0 ) iRange(2) = nNode; 
        Node1 = max(1,iRange(1)); Node2 = min(nNode,iRange(2));  
        
        DO 110 in = Node1,Node2
          BLU(in) = BLU(in) -SRCX  
          SRCndLQ(in) = SRCX
  110   CONTINUE
  100 CONTINUE
C----------------------------------------------------------------------C
C     End of NDLQSR group.
C----------------------------------------------------------------------C
C
      RETURN
      END

C======================================================================C
C
      SUBROUTINE NDSPSR (is)
C
C----------------------------------------------------------------------C
C
C     NDSPSR: inserting NoDal SPecies sources or sinks to water flow 
!             equation.
C
C----------------------------------------------------------------------C   
      USE DRLINK  
      USE DRNODE  
      USE DRWFLD
      USE FLDINX
      USE LQNODE 
      USE NUMBRS
      USE POINTS 
      USE SOLVAR
      USE LINEQN  
      USE SPNODE
C----------------------------------------------------------------------C
C     Implicit Double Precision.
C----------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
!----------------------------------------------------------------------C 
!     Type Declarations. 
!----------------------------------------------------------------------C  
      INTEGER*4 iRange(2)
C----------------------------------------------------------------------C
C     Loop over sources.
C----------------------------------------------------------------------C 
      SRCndSP = 0.0d0;
      DO 100 NS = 1,NDSRSP(is)
        NSB = KNDSSP(1,NS,is)
        NSE = KNDSSP(2,NS,is)
C----------------------------------------------------------------------C
C     Constant source.
C----------------------------------------------------------------------C
        IF( INDSPV(NS,is) .EQ. 1 ) THEN  
          select case (is)  
            case (1)  
              IF( TIME .GT. SNDSP_TM1(NSB) .AND. 
     &            TIME-DT .LT. SNDSP_TM1(NSE) ) THEN
                TMIN = MAX( TIME-DT, SNDSP_TM1(NSB) )
                TMAX = MIN( TIME, SNDSP_TM1(NSE) )
                DTSRC = TMAX -TMIN 
                SRCX = SNDSP1(NSB)*DTSRC/DT
              ELSE
                SRCX = 0.D+0
              ENDIF 
            case (2)  
              IF( TIME .GT. SNDSP_TM2(NSB) .AND. 
     &            TIME-DT .LT. SNDSP_TM2(NSE) ) THEN
                TMIN = MAX( TIME-DT, SNDSP_TM2(NSB) )
                TMAX = MIN( TIME, SNDSP_TM2(NSE) )
                DTSRC = TMAX -TMIN 
                SRCX = SNDSP2(NSB)*DTSRC/DT
              ELSE
                SRCX = 0.D+0
              ENDIF 
            case (3)  
              IF( TIME .GT. SNDSP_TM3(NSB) .AND. 
     &            TIME-DT .LT. SNDSP_TM3(NSE) ) THEN
                TMIN = MAX( TIME-DT, SNDSP_TM3(NSB) )
                TMAX = MIN( TIME, SNDSP_TM3(NSE) )
                DTSRC = TMAX -TMIN 
                SRCX = SNDSP3(NSB)*DTSRC/DT
              ELSE
                SRCX = 0.D+0
              ENDIF 
          end select
C----------------------------------------------------------------------C
C     Tabular source.
C----------------------------------------------------------------------C
        ELSEIF( INDSPV(NS,is) .EQ. 2 ) THEN  
          select case (is)    
            case (1)  
              IF( TIME .GT. SNDSP_TM1(NSB) .AND.
     &            TIME-DT .LT. SNDSP_TM1(NSE)) THEN
                SRCX = 0.D+0  
!                HR = 0.D+0 
!                PART = 0.D+0 
                DO 140 I = NSB+1,NSE
                  IF( TIME .GT. SNDSP_TM1(I-1) .AND.
     &                TIME-DT .LT. SNDSP_TM1(I)) THEN
                    TMIN = MAX( TIME-DT, SNDSP_TM1(I-1) )
                    TMAX = MIN( TIME, SNDSP_TM1(I) )
                    DTSRC = TMAX -TMIN
                    TMID = 0.5*(TMIN +TMAX)
                    NLEN = 2  
!                    CALL LININT( TMID,SRCY,ZZ,SRCSP_TM(I-1),SRCSP(I-1),NLEN) 
                    SRCY = SNDSP1(I-1)
                    SRCX = SRCX + SRCY*DTSRC/DT 
                  ENDIF  
  140           CONTINUE   
              ELSE
                SRCX = 0.D+0
              ENDIF
            case (2)  
              IF( TIME .GT. SNDSP_TM2(NSB) .AND.
     &            TIME-DT .LT. SNDSP_TM2(NSE)) THEN
                SRCX = 0.D+0  
!                HR = 0.D+0 
!                PART = 0.D+0 
                DO 150 I = NSB+1,NSE
                  IF( TIME .GT. SNDSP_TM2(I-1) .AND.
     &                TIME-DT .LT. SNDSP_TM2(I)) THEN
                    TMIN = MAX( TIME-DT, SNDSP_TM2(I-1) )
                    TMAX = MIN( TIME, SNDSP_TM2(I) )
                    DTSRC = TMAX -TMIN
                    TMID = 0.5*(TMIN +TMAX)
                    NLEN = 2  
!                    CALL LININT( TMID,SRCY,ZZ,SRCSP_TM2(I-1),SRCSP2(I-1),NLEN) 
                    SRCY = SNDSP2(I-1)
                    SRCX = SRCX + SRCY*DTSRC/DT 
                  ENDIF  
  150           CONTINUE   
              ELSE
                SRCX = 0.D+0
              ENDIF
            case (3)  
              IF( TIME .GT. SNDSP_TM3(NSB) .AND.
     &            TIME-DT .LT. SNDSP_TM3(NSE)) THEN
                SRCX = 0.D+0  
 !               HR = 0.D+0 
 !               PART = 0.D+0 
                DO 160 I = NSB+1,NSE
                  IF( TIME .GT. SNDSP_TM3(I-1) .AND.
     &                TIME-DT .LT. SNDSP_TM3(I)) THEN
                    TMIN = MAX( TIME-DT, SNDSP_TM3(I-1) )
                    TMAX = MIN( TIME, SNDSP_TM3(I) )
                    DTSRC = TMAX -TMIN
                    TMID = 0.5*(TMIN +TMAX)
                    NLEN = 2  
!                    CALL LININT( TMID,SRCY,ZZ,SRCSP_TM3(I-1),SRCSP3(I-1),NLEN) 
                    SRCY = SNDSP3(I-1)
                    SRCX = SRCX + SRCY*DTSRC/DT 
                  ENDIF  
  160           CONTINUE   
              ELSE
                SRCX = 0.D+0
              ENDIF
          end select
        ENDIF
C----------------------------------------------------------------------C
C     Loop over source domain.
C----------------------------------------------------------------------C  
        iRange(1) = INDSSP(1,NS,is); iRange(2) = INDSSP(2,NS,is);
        if( iRange(1) == 0 ) iRange(1) = 1; 
        if( iRange(2) == 0 ) iRange(2) = nNode; 
        Node1 = max(1,iRange(1)); Node2 = min(nNode,iRange(2));  
        
        DO 110 in = Node1,Node2  
          if( INDTP(NS,is) == 1 ) then     
             SRCndSP(in) = SRCX *SRCndLQ(in)  
          else  
             SRCndSP(in) = SRCX   
          endif    
  110   CONTINUE
  100 CONTINUE
C----------------------------------------------------------------------C
C     Modify right-hand part of the linear equation system.
C----------------------------------------------------------------------C  
      do inode=1, Nnode
         KKK = KK(inode)  
         if( KKK == 1 ) cycle;
         Qout = 0.d0; 
         do j=1,KKK  
            Lin = -LK(inode,j) 
            if( Lin > 0 ) then
               nIn = nd(Lin,IL(Lin))  
            else  
               nIn = nd(-Lin,1)  
            endif   
!                
            if( Lin*QH(nIn) < 0 ) then  
               Qout = Qout +abs(QH(nIn))  
            endif    
         enddo  
!              
         do j=1,KKK 
           Lin = -LK(inode,j) 
           if( Lin > 0 ) then
             nIn = nd(Lin,IL(Lin))  
           else  
             nIn = nd(-Lin,1)  
           endif  
!                
           if( Lin*QH(nIn) < 0 ) then  
             BLC(nIn) = BLC(nIn) +SRCndSP(inode)/Qout *abs(QH(nIn))    
           endif   
         enddo  
      enddo                
C----------------------------------------------------------------------C
C     End of NDSPSR group.
C----------------------------------------------------------------------C
C
      RETURN
      END

C======================================================================C
