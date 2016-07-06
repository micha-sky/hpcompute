!==============================================================================C 
  SUBROUTINE SPPREV
!------------------------------------------------------------------------------C
!
!     SPPREV: save SPecies transport variables at PREVious
!             time step.
!
!------------------------------------------------------------------------------C  
  USE FILINX 
  USE FLDINX
  USE DRSPFL 
  USE DRLINK 
  USE SOLVAR
!------------------------------------------------------------------------------C
  IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
!------------------------------------------------------------------------------C
!------------------------------------------------------------------------------C
!     Assign the primary variables at previous time step.
!------------------------------------------------------------------------------C 
    do 100 L = 1,NLINK  
      do I =1, IL(L)  
        N = ND(L,I)  
        CHO(N)  = CH(N)
        CSHO(N) = CSH(N)
        CBO(N)  = CB(N)  
      enddo
100 CONTINUE
!------------------------------------------------------------------------------C
!     End of SPPREV group. 
!----------------------------------------------------------------------C 
! 
  RETURN
  END
!==============================================================================C 

!==============================================================================C 
  SUBROUTINE RDSEBC
!------------------------------------------------------------------------------C
!
!     RDSEBC: ReaD SEdiment Boundary Conditions input group.
!
!------------------------------------------------------------------------------C
  USE DRBCSE 
  USE DRLINK 
  USE DRNAME
  USE DRNODE
  USE FILINX
  USE FLDINX 
  use intersub 
  USE SOLVAR
!------------------------------------------------------------------------------C
  IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
      INCLUDE 'utils.inc' ! subroutine interfaces
!------------------------------------------------------------------------------C
!     For i-th boundary node (1 <= i <= NBCS):
!        IBCNS(i)  - number of i-th node;
!        IBCDS(i)  - direction for i-th node;
!        IBCTS(i)  - Type of boundary condition;
!        MBCS(i)   - table number of boundary values;
!     For j-th table of boundary values:
!        NBRNS(1,j) - first line in the table;
!        NBRNS(2,j) - last line in the table.
!     For k-th line (1 <= k <= NBTLS) in table of boundary values:
!        BCVS(k)  -  boundary value for X-direction water-flow velocity
!        BCTS(k)  -  time moment.
!------------------------------------------------------------------------------C 
!------------------------------------------------------------------------------C
  CHARACTER*500 str,str1
  CHARACTER*80 Units,Type
  CHARACTER*150 TABLE,TABLEFN
  REAL*8 Value

  CHARACTER*150 vStr
  INTEGER*4 lvS
!---------------------TEMPORARY**ARRAYS-----------------------------------
  CHARACTER*150, POINTER ::  TBNAME(:)=>NULL() !--NBT
!------------------------------------------------------------------------------C
  if( iSolve(2) /= 1 ) return
!------------------------------------------------------------------------------C
!     Write header to output file.
!------------------------------------------------------------------------------C
  write(IWR,'(/A)') ' Sediment Boundary Conditions'
  write(IWR,'(A)')  ' ----------------------------'
!------------------------------------------------------------------------------C
!     Initialize boundary condition variables.
!------------------------------------------------------------------------------C
  NBCS = 0
  NBTS = 0
  NBTLS = 0 !--Current index in all table lines
  NVTABLE = 0 !--Number of non constant tables

!------------------------------------------------------------------------------C
!     Read the number of BCS input lines.
!------------------------------------------------------------------------------C
  read(IRD,*) nLNS
!------------------------------------------------------------------------------C
!     ALLOCATING ARRAYS.
!------------------------------------------------------------------------------C
   nk = nBound;
   ALLOCATE(IBCNS(nk), IBCTS(nk), MBCS(nk));  
   IBCNS=0; IBCTS=0; MBCS=0;
!------------------------------------------------------------------------------C
!     ALLOCATING TEMPORARY ARRAYS.
!------------------------------------------------------------------------------C
   nk=nLNS+1; !--(min size)
   ALLOCATE(NBRNS(2,nk)); NBRNS=0;
   ALLOCATE(TBNAME(nk)); TBNAME=' ';
   ALLOCATE(BCVS(nk), BCTS(nk));  BCVS=0.d0; BCTS=0.d0;
!------------------------------------------------------------------------------C
!     Read BCS input lines.
!------------------------------------------------------------------------------C
  do iLNS = 1,nLNS
    read(IRD,'(A)') str
    istart=1
    if (iLNS/=1) write(IWR,*)
!------------------------------------------------------------------------------C
!     Read 'boundary type'.
!------------------------------------------------------------------------------C
    ires = RDStr(',',istart,str,str1)
    IBType=0
    if (str1(1:9) == 'dirichlet') IBType = 1
    if (str1(1:7) == 'neumann')   IBType = 2
    if (str1(1:7) == 'outflow')   IBType = 3
    if (IBType==0) then
      call Msg(0,IWR,'INPUT ERROR! - Unrecognized Sediment Boundary Condition Type: '//trim(str1))
      STOP
    endif
!------------------------------------------------------------------------------C
!     Read boundary temporal variation.
!------------------------------------------------------------------------------C
    ITEMP = 0 !--Current temporal variation
    ITAB = 0 !--Current table number
    ires = RDStr(',',istart,str,Type)
    if (Type(1:8) == 'constant') then
      ITEMP = 1
!------------------------------------------------------------------------------C
!     Read boundary values and units.
!------------------------------------------------------------------------------C
      ires = RDR8(',',istart,str,Value)
      ires = RDStr(',',istart,str,Units)
      call ConvToSI(Value,Units)
      NBTS = NBTS + 1
      ITAB = NBTS
      TBNAME(NBTS) = 'constant' !--Current table name
      NBTLS = NBTLS + 1
      NBRNS(1,NBTS) = NBTLS
      NBRNS(2,NBTS) = NBTLS
      BCVS(NBTLS) = Value

    elseif (Type(1:7) == 'tabular') then
      ITEMP = 2
!------------------------------------------------------------------------------C
!     Read table name of boundary values.
!------------------------------------------------------------------------------C
      ires = RDStr(',',istart,str,TABLE)
!      ires = RDStr(',',istart,str,str1) !--Unused value
      do N = 1,NBTS
        if (TBNAME(N)==TABLE) then !--Already mentioned table name
          ITAB = N
          GOTO 250
        endif
      enddo
      NVTABLE = NVTABLE +1
      NBTS = NBTS +1
      TBNAME(NBTS) = TABLE
      ITAB = NBTS
 250  CONTINUE
    endif
!------------------------------------------------------------------------------C
!     Read nodal range of boundary (2 integers).
!------------------------------------------------------------------------------C
    ires = RDI4(',',istart,str,IBS)
    ires = RDI4(',',istart,str,IBE)
!------------------------------------------------------------------------------C
!     Write input boundary condition data to output file.
!------------------------------------------------------------------------------C
    WRITE (IWR,'(2(A,I4))') ' Boundary Domain (Nodes): I = ',IBS,' to ',IBE
!------------------------------------------------------------------------------C
!     Write input boundary condition data to output file.
!------------------------------------------------------------------------------C
    vStr=''
    if (IBType==1) then
      write(IWR,'(A)') ' Boundary Condition: Dirichlet @ Boundary'
      vStr='   Sediment Concentration,    kg/m^3: '
    elseif (IBType == 2 ) then
      write(IWR,'(A)') ' Boundary Condition: Neumann @ Boundary'
      vStr='   Sediment Flux,               kg/s: '
    elseif (IBType == 3 ) then
      write(IWR,'(A)') ' Boundary Condition: Outflow @ Boundary'
    endif
    lvS=len_trim(vStr)

    if (ITEMP==1) then
      write(IWR,'(A)')     ' Temporal Variation: Constant'
      if (lvS/=0) write(IWR,'(A,1PG12.5)') vStr(1:lvS), Value
    elseif (ITEMP == 2) then
      write(IWR,'(2A)')     ' Temporal Variation: Tabular - ', trim(TBNAME(iLNS))
    endif
!------------------------------------------------------------------------------C
!     Assign values to boundary variables.
!------------------------------------------------------------------------------C
    ISTOP = 0
    DO 300 I = IBS,IBE
!------------------------------------------------------------------------------C
!     Check for boundary values applied to interior surfaces.
!------------------------------------------------------------------------------C
      IERROR = 0
      IF( (IBTYPE .EQ. 2) .AND. (KK(I) > 1) ) THEN
         IERROR = 1
      ENDIF
      IF( IERROR .EQ. 1 ) THEN
        WRITE (IWR,9001) I,NodeName(I)
        ISTOP = 1
      ENDIF
      NBCS = NBCS +1
      IBCNS(NBCS) = I 
      IBCTS(NBCS) = IBTYPE
      MBCS(NBCS)  = ITAB
300   CONTINUE
!------------------------------------------------------------------------------C
  enddo
  IF( ISTOP .EQ. 1 ) THEN
    WRITE (*,'(/2A)') ' INPUT OR COMPILATION ERROR: ',           &
       'Sediment Boundary Conditions Input (check output file).' 
    WRITE (IWR,'(/2A)') ' INPUT OR COMPILATION ERROR: ',         &
       'Sediment Boundary Conditions Input.'
      STOP
  ENDIF  
!------------------------------------------------------------------------------C
!     If no temporal boundary tables were mentioned above then finish.
!------------------------------------------------------------------------------C
  if (NVTABLE==0) GOTO 900

!------------------------------------------------------------------------------C
!     Read boundary tables.
!     Read the number of BC table input lines.
!------------------------------------------------------------------------------C
  read(IRD,*) NLIN
!------------------------------------------------------------------------------C
!     REALLOCATE  BCVS,BCTS arrays.
!------------------------------------------------------------------------------C
  nk=NBTLS+NLIN+1;
  call ReallocateR8(NBTLS,nk,BCVS)
  call ReallocateR8(NBTLS,nk,BCTS)
!------------------------------------------------------------------------------C
!     Reading table.
!------------------------------------------------------------------------------C
  call ReadTableCOASTOX(2,IRD,NLIN,NBTS,NBTLS,NBRNS,TBNAME,BCTS,1,BCVS,BCVS,BCVS,.true.)


 900 CONTINUE
!------------------------------------------------------------------------------C
!     Write boundary tables.
!------------------------------------------------------------------------------C
  do N = 1,NBTS
    if (TBNAME(N)(1:8) == 'constant') CYCLE
    do i = 1,NBCS
      if (MBCS(i)==N) then
        IBType = IBCTS(i)
        EXIT
      endif
    enddo

    write(IWR,'(/2A)') ' Boundary Table: ', TBNAME(N)
    NLIN = NBRNS(2,N) - NBRNS(1,N) + 1
    do i = 1,NLIN
      j = NBRNS(1,N) + i - 1
      if (IBType==1) then
        write(IWR,'(A,I4,4(A,1PG12.5))') ' Entry No.',i,'  Time, s: ',BCTS(j),' Sediment Concentration, kg/m^3 liquid: ',BCVS(j)
      elseif (IBType==2) then
        write(IWR,'(A,I4,3(A,1PG12.5))') ' Entry No.',i,'  Time, s: ',BCTS(j),' Sediment Flux,                   kg/s: ',BCVS(j)
      elseif (IBType == 3 ) then
        write(IWR,'(A,I4,3(A,1PG12.5))') ' Entry No.',i,'  Time, s: ',BCTS(j),&
              ' Sediment Concentration Derivation at Outflow Boundary: ',BCVS(j)
      endif
    enddo
  enddo
!------------------------------------------------------------------------------C
!     Deallocating Temporary arrays.
!------------------------------------------------------------------------------C
  if (ASSOCIATED(TBNAME)) DEALLOCATE(TBNAME)

!------------------------------------------------------------------------------C
!     Format Statements.
!------------------------------------------------------------------------------C
9001 FORMAT(/' INPUT ERROR! - A Sediment boundary condition has been specified for an interior node.',/, &
  ' Node: ',I0, ',  Boundary Direction Index: ',I0)
!------------------------------------------------------------------------------C
!     End of RDSEBC group. 
!----------------------------------------------------------------------C 
! 
  RETURN
  END
!==============================================================================C 

!==============================================================================C 
  SUBROUTINE RDSPBC
!------------------------------------------------------------------------------C
!
!     RDSPBC: ReaD SPecies Boundary Conditions input group.
!
!------------------------------------------------------------------------------C
  USE DRBCCL 
  USE DRNAME
  USE DRNODE
  USE DRLINK
  USE FILINX
  USE FLDINX
  use intersub
  USE NUMBRS
  USE SOLVAR
!------------------------------------------------------------------------------C
  IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
      INCLUDE 'utils.inc' ! subroutine interfaces
!------------------------------------------------------------------------------C
!     For i-th boundary node (1 <= i <= NBCS):
!        IBCNS(i)  - number of i-th node;
!        IBCDS(i)  - direction for i-th node;
!        IBCTS(i)  - Type of boundary condition;
!        MBCS(i)   - table number of boundary values;
!     For j-th table of boundary values:
!        NBRNS(1,j) - first line in the table;
!        NBRNS(2,j) - last line in the table.
!     For k-th line (1 <= k <= NBTLS) in table of boundary values:
!        BCVS(k)  -  boundary value for X-direction water-flow velocity
!        BCTS(k)  -  time moment.
!------------------------------------------------------------------------------C 
!------------------------------------------------------------------------------C
  CHARACTER*500 str,str1
  CHARACTER*80 Units,Type
  CHARACTER*150 TABLE,TABLEFN
  REAL*8 Value

  CHARACTER*150 vStr
  INTEGER*4 lvS
!---------------------TEMPORARY**ARRAYS-----------------------------------
  CHARACTER*150, POINTER ::  TBNAME(:)=>NULL() !--NBT
!------------------------------------------------------------------------------C
  IF( ISOLVE(3) /= 1 ) RETURN
!------------------------------------------------------------------------------C
!     Write header to output file.
!------------------------------------------------------------------------------C
  write(IWR,'(/A)') ' Species Boundary Conditions'
  write(IWR,'(A)')  ' ---------------------------'
!------------------------------------------------------------------------------C
!     Initialize boundary condition variables.
!------------------------------------------------------------------------------C
  NBCT = 0
  NBTT = 0
  NBTLT = 0 !--Current index in all table lines
  NVTABLE = 0 !--Number of non constant tables

!------------------------------------------------------------------------------C
!     Read the number of BCS input lines.
!------------------------------------------------------------------------------C
  read(IRD,*) nLNS
!------------------------------------------------------------------------------C
!     ALLOCATING ARRAYS.
!------------------------------------------------------------------------------C
   nk = nBound;
   ALLOCATE(IBCNT(nk), IBCTT(nk), MBCT(nk));  
   IBCNT=0; IBCTT=0; MBCT=0;
!------------------------------------------------------------------------------C
!     ALLOCATING TEMPORARY ARRAYS.
!------------------------------------------------------------------------------C
   nk=nLNS+1; !--(min size)
   ALLOCATE(NBTST(2,nk)); NBTST=0;
   ALLOCATE(TBNAME(nk)); TBNAME=' ';
   ALLOCATE(BCVT(nk), BCTT(nk));  BCVT=0.d0; BCTT=0.d0;
!------------------------------------------------------------------------------C
!     Read BCS input lines.
!------------------------------------------------------------------------------C
  do iLNS = 1,nLNS
    read(IRD,'(A)') str
    istart=1
    if (iLNS/=1) write(IWR,*)
!------------------------------------------------------------------------------C
!     Read 'boundary type'.
!------------------------------------------------------------------------------C
    ires = RDStr(',',istart,str,str1)
    IBType=0
    if (str1(1:9) == 'dirichlet') IBType = 1
    if (str1(1:7) == 'neumann')   IBType = 2
    if (str1(1:7) == 'outflow')   IBType = 3
    if (IBType==0) then
      call Msg(0,IWR,'INPUT ERROR! - Unrecognized Species Boundary Condition Type: '//trim(str1))
      STOP
    endif
!------------------------------------------------------------------------------C
!     Read boundary temporal variation.
!------------------------------------------------------------------------------C
    ITEMP = 0 !--Current temporal variation
    ITAB = 0 !--Current table number
    ires = RDStr(',',istart,str,Type)
    if (Type(1:8) == 'constant') then
      ITEMP = 1
!------------------------------------------------------------------------------C
!     Read boundary values and units.
!------------------------------------------------------------------------------C
      ires = RDR8(',',istart,str,Value)
      ires = RDStr(',',istart,str,Units)
      call ConvToSI(Value,Units)
      NBTT = NBTT + 1
      ITAB = NBTT
      TBNAME(NBTT) = 'constant' !--Current table name
      NBTLT = NBTLT + 1
      NBTST(1,NBTT) = NBTLT
      NBTST(2,NBTT) = NBTLT
      BCVT(NBTLT) = Value

    elseif (Type(1:7) == 'tabular') then
      ITEMP = 2
!------------------------------------------------------------------------------C
!     Read table name of boundary values.
!------------------------------------------------------------------------------C
      ires = RDStr(',',istart,str,TABLE)
!      ires = RDStr(',',istart,str,str1) !--Unused value
      do N = 1,NBTT
        if (TBNAME(N)==TABLE) then !--Already mentioned table name
          ITAB = N
          GOTO 250
        endif
      enddo
      NVTABLE = NVTABLE +1
      NBTT = NBTT +1
      TBNAME(NBTT) = TABLE
      ITAB = NBTT
 250  CONTINUE
    endif
!------------------------------------------------------------------------------C
!     Read nodal range of boundary (2 integers).
!------------------------------------------------------------------------------C
    ires = RDI4(',',istart,str,IBS)
    ires = RDI4(',',istart,str,IBE)
!------------------------------------------------------------------------------C
!     Write input boundary condition data to output file.
!------------------------------------------------------------------------------C
    WRITE (IWR,'(2(A,I4))') ' Boundary Domain (Nodes): I = ',IBS,' to ',IBE
!------------------------------------------------------------------------------C
!     Write input boundary condition data to output file.
!------------------------------------------------------------------------------C
    vStr=''
    if (IBType==1) then
      write(IWR,'(A)') ' Boundary Condition: Dirichlet @ Boundary'
      vStr='   Species Concentration,    1/m^3: '
    elseif (IBType == 2 ) then
      write(IWR,'(A)') ' Boundary Condition: Neumann @ Boundary'
      vStr='   Species Flux,               1/s: '
    elseif (IBType == 3 ) then
      write(IWR,'(A)') ' Boundary Condition: Outflow @ Boundary'
    endif
    lvS=len_trim(vStr)

    if (ITEMP==1) then
      write(IWR,'(A)')     ' Temporal Variation: Constant'
      if (lvS/=0) write(IWR,'(A,1PG12.5)') vStr(1:lvS), Value
    elseif (ITEMP == 2) then
      write(IWR,'(2A)')     ' Temporal Variation: Tabular - ', trim(TBNAME(iLNS))
    endif
!------------------------------------------------------------------------------C
!     Assign values to boundary variables.
!------------------------------------------------------------------------------C
    ISTOP = 0
    DO 300 I = IBS,IBE
!------------------------------------------------------------------------------C
!     Check for boundary values applied to interior surfaces.
!------------------------------------------------------------------------------C
      IERROR = 0
      IF( (IBTYPE .EQ. 2) .AND. (KK(I) > 1) ) THEN
         IERROR = 1
      ENDIF
      IF( IERROR .EQ. 1 ) THEN
        WRITE (IWR,9001) I,NodeName(I)
        ISTOP = 1
      ENDIF
      NBCT = NBCT +1
      IBCNT(NBCT) = I 
      IBCTT(NBCT) = IBTYPE
      MBCT(NBCT)  = ITAB
300   CONTINUE
!------------------------------------------------------------------------------C
  enddo
  IF( ISTOP .EQ. 1 ) THEN
    WRITE (*,'(/2A)') ' INPUT OR COMPILATION ERROR: ',           &
       'Species Boundary Conditions Input (check output file).' 
    WRITE (IWR,'(/2A)') ' INPUT OR COMPILATION ERROR: ',         &
       'Species Boundary Conditions Input.'
      STOP
  ENDIF  
!------------------------------------------------------------------------------C
!     If no temporal boundary tables were mentioned above then finish.
!------------------------------------------------------------------------------C
  if (NVTABLE==0) GOTO 900
!------------------------------------------------------------------------------C
!     Read boundary tables.
!     Read the number of BC table input lines.
!------------------------------------------------------------------------------C
  read(IRD,*) NLIN
!------------------------------------------------------------------------------C
!     REALLOCATE  BCVS,BCTS arrays.
!------------------------------------------------------------------------------C
  nk=NBTLT+NLIN+1;
  call ReallocateR8(NBTLT,nk,BCVT)
  call ReallocateR8(NBTLT,nk,BCTT)
!------------------------------------------------------------------------------C
!     Reading table.
!------------------------------------------------------------------------------C
  call ReadTableCOASTOX(2,IRD,NLIN,NBTT,NBTLT,NBTST,TBNAME,BCTT,1,BCVT,BCVT,BCVT,.true.)
 900 CONTINUE
!------------------------------------------------------------------------------C
!     Write boundary tables.
!------------------------------------------------------------------------------C
  do N = 1,NBTT
    if (TBNAME(N)(1:8) == 'constant') CYCLE
    do i = 1,NBCT
      if (MBCT(i)==N) then
        IBType = IBCTT(i)
        EXIT
      endif
    enddo

    write(IWR,'(/2A)') ' Boundary Table: ', TBNAME(N)
    NLIN = NBTST(2,N) - NBTST(1,N) + 1
    do i = 1,NLIN
      j = NBTST(1,N) + i - 1
      if (IBType==1) then
        write(IWR,'(A,I4,4(A,1PG12.5))') ' Entry No.',i,'  Time, s: ',BCTT(j),' Species Concentration, kg/m^3 liquid: ',BCVT(j)
      elseif (IBType==2) then
        write(IWR,'(A,I4,3(A,1PG12.5))') ' Entry No.',i,'  Time, s: ',BCTT(j),' Species Flux,                   kg/s: ',BCVT(j)
      elseif (IBType == 3 ) then
        write(IWR,'(A,I4,3(A,1PG12.5))') ' Entry No.',i,'  Time, s: ',BCTT(j),&
            ' Species Concentration Derivation at Outflow Boundary: ',BCVT(j)
      endif
    enddo
  enddo
!------------------------------------------------------------------------------C
!     Deallocating Temporary arrays.
!------------------------------------------------------------------------------C
  if (ASSOCIATED(TBNAME)) DEALLOCATE(TBNAME)
!------------------------------------------------------------------------------C
!     Format Statements.
!------------------------------------------------------------------------------C
9001 FORMAT(/' INPUT ERROR! - A Species boundary condition has been specified for an interior node.',/, &
  ' Node: ',I0, ',  Boundary Direction Index: ',I0)
!------------------------------------------------------------------------------C
!     End of RDSPBC group. 
!----------------------------------------------------------------------C 
! 
  RETURN
  END
!==============================================================================C 

!==============================================================================C 
  SUBROUTINE RDPTBC
!------------------------------------------------------------------------------C
!
!     RDPTBC: ReaD ParTiculate Boundary Conditions input group.
!
!------------------------------------------------------------------------------C
  USE DRBCCP 
  USE DRLINK
  USE DRNAME
  USE DRNODE
  USE FILINX
  USE FLDINX
  use intersub 
  USE NUMBRS
  USE SOLVAR
!------------------------------------------------------------------------------C
  IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
      INCLUDE 'utils.inc' ! subroutine interfaces
!------------------------------------------------------------------------------C
!     For i-th boundary node (1 <= i <= NBCS):
!        IBCNS(i)  - number of i-th node;
!        IBCDS(i)  - direction for i-th node;
!        IBCTS(i)  - Type of boundary condition;
!        MBCS(i)   - table number of boundary values;
!     For j-th table of boundary values:
!        NBRNS(1,j) - first line in the table;
!        NBRNS(2,j) - last line in the table.
!     For k-th line (1 <= k <= NBTLS) in table of boundary values:
!        BCVS(k)  -  boundary value for X-direction water-flow velocity
!        BCTS(k)  -  time moment.
!------------------------------------------------------------------------------C
!------------------------------------------------------------------------------C 
!----------------------------------------------------------------------C 
!     Type Declarations. 
!----------------------------------------------------------------------C 
  CHARACTER*500 str,str1
  CHARACTER*80 Units,Type
  CHARACTER*150 TABLE,TABLEFN
  REAL*8 Value

  CHARACTER*150 vStr
  INTEGER*4 lvS
!---------------------TEMPORARY**ARRAYS-----------------------------------
  CHARACTER*150, POINTER ::  TBNAME(:)=>NULL() !--NBT
!------------------------------------------------------------------------------C
  if( iSolve(4) /= 1 ) return
!------------------------------------------------------------------------------C
!     Write header to output file.
!------------------------------------------------------------------------------C
  write(IWR,'(/A)') ' Particulate Boundary Conditions'
  write(IWR,'(A)')  ' -------------------------------'
!------------------------------------------------------------------------------C
!     Initialize boundary condition variables.
!------------------------------------------------------------------------------C
  NBCR = 0
  NBTR = 0
  NBTLR = 0 !--Current index in all table lines
  NVTABLE = 0 !--Number of non constant tables
!------------------------------------------------------------------------------C
!     Read the number of BCS input lines.
!------------------------------------------------------------------------------C
  read(IRD,*) nLNS
!------------------------------------------------------------------------------C
!     ALLOCATING ARRAYS.
!------------------------------------------------------------------------------C
   nk = nBound;
   ALLOCATE(IBCNR(nk), IBCTR(nk), MBCR(nk));  
   IBCNR=0; IBCTR=0; MBCR=0;
!------------------------------------------------------------------------------C
!     ALLOCATING TEMPORARY ARRAYS.
!------------------------------------------------------------------------------C
   nk=nLNS+1; !--(min size)
   ALLOCATE(NBTSR(2,nk)); NBTSR=0;
   ALLOCATE(TBNAME(nk)); TBNAME=' ';
   ALLOCATE(BCVR(nk), BCTR(nk));  BCVR=0.d0; BCTR=0.d0;
!------------------------------------------------------------------------------C
!     Read BCS input lines.
!------------------------------------------------------------------------------C
  do iLNS = 1,nLNS
    read(IRD,'(A)') str
    istart=1
    if (iLNS/=1) write(IWR,*)
!------------------------------------------------------------------------------C
!     Read 'boundary type'.
!------------------------------------------------------------------------------C
    ires = RDStr(',',istart,str,str1)
    IBType=0
    if (str1(1:9) == 'dirichlet') IBType = 1
    if (str1(1:7) == 'neumann')   IBType = 2
    if (str1(1:7) == 'outflow')   IBType = 3
    if (IBType==0) then
      call Msg(0,IWR,'INPUT ERROR! - Unrecognized Sediment Boundary Condition Type: '//trim(str1))
      STOP
    endif
!------------------------------------------------------------------------------C
!     Read boundary temporal variation.
!------------------------------------------------------------------------------C
    ITEMP = 0 !--Current temporal variation
    ITAB = 0 !--Current table number
    ires = RDStr(',',istart,str,Type)
    if (Type(1:8) == 'constant') then
      ITEMP = 1
!------------------------------------------------------------------------------C
!     Read boundary values and units.
!------------------------------------------------------------------------------C
      ires = RDR8(',',istart,str,Value)
      ires = RDStr(',',istart,str,Units)
      call ConvToSI(Value,Units)
      NBTR = NBTR + 1
      ITAB = NBTR
      TBNAME(NBTR) = 'constant' !--Current table name
      NBTLR = NBTLR + 1
      NBTSR(1,NBTR) = NBTLR
      NBTSR(2,NBTR) = NBTLR
      BCVR(NBTLR) = Value

    elseif (Type(1:7) == 'tabular') then
      ITEMP = 2
!------------------------------------------------------------------------------C
!     Read table name of boundary values.
!------------------------------------------------------------------------------C
      ires = RDStr(',',istart,str,TABLE)
!      ires = RDStr(',',istart,str,str1) !--Unused value
      do N = 1,NBTR
        if (TBNAME(N)==TABLE) then !--Already mentioned table name
          ITAB = N
          GOTO 250
        endif
      enddo
      NVTABLE = NVTABLE +1
      NBTR = NBTR +1
      TBNAME(NBTR) = TABLE
      ITAB = NBTR
 250  CONTINUE
    endif
!------------------------------------------------------------------------------C
!     Read nodal range of boundary (2 integers).
!------------------------------------------------------------------------------C
    ires = RDI4(',',istart,str,IBS)
    ires = RDI4(',',istart,str,IBE)
!------------------------------------------------------------------------------C
!     Write input boundary condition data to output file.
!------------------------------------------------------------------------------C
    WRITE (IWR,'(2(A,I4))') ' Boundary Domain (Nodes): I = ',IBS,' to ',IBE
!------------------------------------------------------------------------------C
!     Write input boundary condition data to output file.
!------------------------------------------------------------------------------C
    vStr=''
    if (IBType==1) then
      write(IWR,'(A)') ' Boundary Condition: Dirichlet @ Boundary'
      vStr='   Particulate Concentration,    1/m^3: '
    elseif (IBType == 2 ) then
      write(IWR,'(A)') ' Boundary Condition: Neumann @ Boundary'
      vStr='   Particulate Flux,               1/s: '
    elseif (IBType == 3 ) then
      write(IWR,'(A)') ' Boundary Condition: Outflow @ Boundary'
    endif
    lvS=len_trim(vStr)

    if (ITEMP==1) then
      write(IWR,'(A)')     ' Temporal Variation: Constant'
      if (lvS/=0) write(IWR,'(A,1PG12.5)') vStr(1:lvS), Value
    elseif (ITEMP == 2) then
      write(IWR,'(2A)')     ' Temporal Variation: Tabular - ', trim(TBNAME(iLNS))
    endif
!------------------------------------------------------------------------------C
!     Assign values to boundary variables.
!------------------------------------------------------------------------------C
    ISTOP = 0
    DO 300 I = IBS,IBE
!------------------------------------------------------------------------------C
!     Check for boundary values applied to interior surfaces.
!------------------------------------------------------------------------------C
      IERROR = 0
      IF( (IBTYPE .EQ. 2) .AND. (KK(I) > 1) ) THEN
         IERROR = 1
      ENDIF
      IF( IERROR .EQ. 1 ) THEN
        WRITE (IWR,9001) I,NodeName(I)
        ISTOP = 1
      ENDIF
      NBCR = NBCR +1
      IBCNR(NBCR) = I 
      IBCTR(NBCR) = IBTYPE
      MBCR(NBCR)  = ITAB
300   CONTINUE
!------------------------------------------------------------------------------C
  enddo
  IF( ISTOP .EQ. 1 ) THEN
    WRITE (*,'(/2A)') ' INPUT OR COMPILATION ERROR: ',           &
       'Particulate Boundary Conditions Input (check output file).' 
    WRITE (IWR,'(/2A)') ' INPUT OR COMPILATION ERROR: ',         &
       'Particulate Boundary Conditions Input.'
      STOP
  ENDIF  
!------------------------------------------------------------------------------C
!     If no temporal boundary tables were mentioned above then finish.
!------------------------------------------------------------------------------C
  if (NVTABLE==0) GOTO 900
!------------------------------------------------------------------------------C
!     Read boundary tables.
!     Read the number of BC table input lines.
!------------------------------------------------------------------------------C
  read(IRD,*) NLIN
!------------------------------------------------------------------------------C
!     REALLOCATE  BCVS,BCTS arrays.
!------------------------------------------------------------------------------C
  nk=NBTLR+NLIN+1;
  call ReallocateR8(NBTLR,nk,BCVR)
  call ReallocateR8(NBTLR,nk,BCTR)
!------------------------------------------------------------------------------C
!     Reading table.
!------------------------------------------------------------------------------C
  call ReadTableCOASTOX(2,IRD,NLIN,NBTR,NBTLR,NBTSR,TBNAME,BCTR,1,BCVR,BCVR,BCVR,.true.)

 900 CONTINUE
!------------------------------------------------------------------------------C
!     Write boundary tables.
!------------------------------------------------------------------------------C
  do N = 1,NBTR
    if (TBNAME(N)(1:8) == 'constant') CYCLE
    do i = 1,NBCR
      if (MBCR(i)==N) then
        IBType = IBCTR(i)
        EXIT
      endif
    enddo

    write(IWR,'(/2A)') ' Boundary Table: ', TBNAME(N)
    NLIN = NBTSR(2,N) - NBTSR(1,N) + 1
    do i = 1,NLIN
      j = NBTSR(1,N) + i - 1
      if (IBType==1) then
        write(IWR,'(A,I4,4(A,1PG12.5))') ' Entry No.',i,'  Time, s: ',BCTR(j),' Particulate Concentration, kg/m^3 liquid: ',BCVR(j)
      elseif (IBType==2) then
        write(IWR,'(A,I4,3(A,1PG12.5))') ' Entry No.',i,'  Time, s: ',BCTR(j),' Particulate Flux,                   kg/s: ',BCVR(j)
      elseif (IBType == 3 ) then
        write(IWR,'(A,I4,3(A,1PG12.5))') ' Entry No.',i,'  Time, s: ',BCTR(j),&
             ' Particulate Concentration Derivation at Outflow Boundary: ',BCVR(j)
      endif
    enddo
  enddo
!------------------------------------------------------------------------------C
!     Deallocating Temporary arrays.
!------------------------------------------------------------------------------C
  if (ASSOCIATED(TBNAME)) DEALLOCATE(TBNAME)
!------------------------------------------------------------------------------C
!     Format Statements.
!------------------------------------------------------------------------------C
9001 FORMAT(/' INPUT ERROR! - A Particulate boundary condition has been specified for an interior node.',/, &
  ' Node: ',I0, ',  Boundary Direction Index: ',I0)
!------------------------------------------------------------------------------C
!     End of RDPTBC group. 
!------------------------------------------------------------------------------C
! 
  RETURN
  END
!==============================================================================C 

!==============================================================================C 
SUBROUTINE BCSP
!------------------------------------------------------------------------------C
!
!     BCSP: Boundary Condition for the SPecies
!             transport equation.
!
!------------------------------------------------------------------------------C
USE DRBCCL 
USE DRLINK 
USE DRSPFL
USE DRWFLD  
USE DRNODE
USE FLDINX 
USE LINEQN
USE NUMBRS
USE SOLVAR  
!------------------------------------------------------------------------------C
IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)  
!------------------------------------------------------------------------------C
!     Type Declarations. 
!------------------------------------------------------------------------------C
 REAL*8, POINTER :: wBCV(:) => null()
!------------------------------------------------------------------------------C
!      IUPW = 4 
 if( .not.associated(wBCV) ) allocate(wBCV(NBTT));
!------------------------------------------------------------------------------C
!     Loop over boundary value tables.
!------------------------------------------------------------------------------C
 do NS = 1,NBTT
   NSB = NBTST(1,NS)
   NSE = NBTST(2,NS)
   !------------------------------------------------------------------------------C
   !     Constant source.
   !------------------------------------------------------------------------------C
   if (NSB == NSE) then
     wBCV(NS) = BCVT(NSB)
   !------------------------------------------------------------------------------C
   !     Tabular source.
   !------------------------------------------------------------------------------C
   else
     wBCV(NS) = 0.d0
     if (TIME > BCTT(NSB) .and.         TIME-DT < BCTT(NSE)) then
       do I = NSB+1,NSE
         if (TIME > BCTT(I-1) .and.             TIME-DT < BCTT(I)) then
           TMIN = max(TIME-DT, BCTT(I-1))
           TMAX = min(TIME, BCTT(I))
           DTBND = TMAX -TMIN
           TMID = 0.5*(TMIN +TMAX)
           NLEN = 2
           call LININT(TMID,Value,ZZ,BCTT(I-1:),BCVT(I-1:),NLEN)
           wBCV(NS) = wBCV(NS) + Value*DTBND/DT
         endif
       enddo
     endif
   endif
 enddo
!------------------------------------------------------------------------------C
!     Species Boundary Conditions.
!------------------------------------------------------------------------------C
 do 200 iBC = 1,NBCT
   N = IBCNT(iBC) 
!   IZN = IZ(N)
!   IZPLN = IZPL(N)
   IBCT = IBCTT(iBC)
   BCS = wBCV(MBCT(iBC))  
   
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
   
   CP = SIGMA*CH(i) + (1.d0 -SIGMA)*CHO(i)
   HP = SIGMA*HR(i) + (1.d0 -SIGMA)*HRO(i)
   QP = SIGMA*QH(i) + (1.d0 -SIGMA)*QHO(i)  
!------------------------------------------------------------------------------C
!     Boundary direction controller.
!------------------------------------------------------------------------------C 
   AB = 0.d0; AP = 0.d0;
   IF( L > 0 ) GOTO 210
   IF( L < 0 ) GOTO 220  
   GOTO 230
!------------------------------------------------------------------------------C
!    Left Link Boundary.
!------------------------------------------------------------------------------C
210   CONTINUE
!     Diffusion terms
!   DLP   = tDFSL*DXIW*HP 
   DLP = 0.d0
   if (IBCT == 3) then
     DLP = 0.d0
   endif 

   AB = DLP + max(QP,ZERO)
   AP = DLP - min(QP,ZERO)
   GOTO 230
!------------------------------------------------------------------------------C
!     Right Link Boundary.
!------------------------------------------------------------------------------C
220   CONTINUE
!     Diffusion terms
!   DLP   = tDFSL*DXIE*HP 
   DLP = 0.d0
   if (IBCT == 3) then
     DLP = 0.d0
   endif

   AB = DLP - min(QP,ZERO)
   AP = DLP + max(QP,ZERO)
230 CONTINUE
!------------------------------------------------------------------------------C
!     Modify the Equation Matrix and solution vector according to the
!     boundary conditions.
!------------------------------------------------------------------------------C 
   NEL = iEq(i) 
   if ((IBCT == 1) .or. (IBCT == 3)) then
     ALC(NEL) = ALC(NEL) + AB
     BLC(i)   = BLC(i) + AB*BCS
   elseif (IBCT == 2) then
     BLC(i) = BLC(i) - BCS
   endif
200 CONTINUE
!------------------------------------------------------------------------------C
!     End of BCSP group. 
!------------------------------------------------------------------------------C
! 
RETURN
END
!==============================================================================C 

!==============================================================================C 
SUBROUTINE BCPT
!------------------------------------------------------------------------------C
!
!     BCPT: Boundary Condition for the ParTiculate
!             transport equation.
!
!------------------------------------------------------------------------------C   
USE DRBCCP 
USE DRLINK 
USE DRNODE
USE DRSPFL
USE DRWFLD 
USE FLDINX   
USE LINEQN
USE NUMBRS
USE SOLVAR  
!------------------------------------------------------------------------------C
IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
!------------------------------------------------------------------------------C
!     Type Declarations. 
!------------------------------------------------------------------------------C
 REAL*8, POINTER :: wBCV(:) => null()
!------------------------------------------------------------------------------C
!      IUPW = 4
 if( .not.associated(wBCV) ) allocate(wBCV(NBTR));
!------------------------------------------------------------------------------C
!     Loop over boundary value tables.
!------------------------------------------------------------------------------C
 do 100 NS = 1,NBTR
   NSB = NBTSR(1,NS)
   NSE = NBTSR(2,NS)
!------------------------------------------------------------------------------C
!     Constant source.
!------------------------------------------------------------------------------C
   if (NSB == NSE) then
     wBCV(NS) = BCVR(NSB)
!------------------------------------------------------------------------------C
!     Tabular source.
!------------------------------------------------------------------------------C
   else
     wBCV(NS) = 0.d0
     if (TIME > BCTR(NSB) .and.         TIME-DT < BCTR(NSE)) then
       do 110 I = NSB+1,NSE
         if (TIME > BCTR(I-1) .and.             TIME-DT < BCTR(I)) then
           TMIN = max(TIME-DT, BCTR(I-1))
           TMAX = min(TIME, BCTR(I))
           DTBND = TMAX -TMIN
           TMID = 0.5*(TMIN +TMAX)
           NLEN = 2
           call LININT(TMID,Value,ZZ,BCTR(I-1:),BCVR(I-1:),NLEN)
           wBCV(NS) = wBCV(NS) + Value*DTBND/DT
         endif
110      CONTINUE
     endif
   endif
100 CONTINUE
!------------------------------------------------------------------------------C
!     Particulate Boundary Conditions.
!------------------------------------------------------------------------------C
 do 200 iBC = 1,NBCR
   N = IBCNR(iBC)
!   IZN = IZ(N)
!   IZPLN = IZPL(N)
   IBCT = IBCTR(iBC)
   BCS = wBCV(MBCR(iBC))
   
   KKK = KK(N) 
   if( KKK > 1 ) then   
     write(*,'(a,i5,a)') ' Node ',N,' is not a boundary node'  
     write(IWR,'(a,i5,a)') ' Node ',N,' is not a boundary node'  
     stop
   endif    
!------------------------------------------------------------------------------C
!     Local grid node indices.
!------------------------------------------------------------------------------C
!------------------------------------------------------------------------------C
!     Local grid node indices.
!------------------------------------------------------------------------------C
   L = LK(N,1) 
   if( L > 0 ) then 
     i = ND(L,1); ! Number of computational node (first in link). 
   else  
     i = ND(-L,IL(-L)); ! Number of computational node (last in link).
   endif    
   
   CP = SIGMA*CSH(i) + (1.d0 -SIGMA)*CSHO(i)
   HP = SIGMA*HR(i) + (1.d0 -SIGMA)*HRO(i)
   QP = SIGMA*QH(i) + (1.d0 -SIGMA)*QHO(i)
!------------------------------------------------------------------------------C
!     Boundary direction controller.
!------------------------------------------------------------------------------C
   AB = 0.d0; AP = 0.d0;
   IF( L > 0 ) GOTO 210
   IF( L < 0 ) GOTO 220  
   GOTO 230
!------------------------------------------------------------------------------C
!    Left Link Boundary.
!------------------------------------------------------------------------------C
210   CONTINUE
!     Diffusion terms
!   DLP   = tDFSL*DYIS*HP 
   DLP = 0.d0;
   if (IBCT == 3) then
     DLP = 0.d0
   endif

   AB = DLP + max(QP,ZERO)
   AP = DLP - min(QP,ZERO)
   GOTO 230
!------------------------------------------------------------------------------C
!     Right Link Boundary.
!------------------------------------------------------------------------------C
220   CONTINUE
!     Diffusion terms
!   DLP   = tDFSL*DXIE*HP 
   DLP = 0.0d0
   if (IBCT == 3) then
     DLP = 0.d0
   endif

   AB = DLP - min(QP,ZERO)
   AP = DLP + max(QP,ZERO)
230   CONTINUE
!------------------------------------------------------------------------------C
!     Modify the Equation Matrix and solution vector according to the
!     boundary conditions.
!------------------------------------------------------------------------------C
   NEL = iEq(i) 
   if ((IBCT == 1) .or. (IBCT == 3)) then
     ALC(NEL) = ALC(NEL) + AB
     BLC(i)   = BLC(i) + AB*BCS
   elseif (IBCT == 2) then
     BLC(i) = BLC(i) - BCS
   endif
200 CONTINUE
!------------------------------------------------------------------------------C
!     End of BCPT group. 
!----------------------------------------------------------------------C 
! 
RETURN
END
!==============================================================================C 

!==============================================================================C 
  SUBROUTINE UPLAYR1(IRET)
!------------------------------------------------------------------------------C
!
!     UPLAYR: transport equations for UPper LAYeR of bottom deposit.
!
!------------------------------------------------------------------------------C   
  USE CONSTS 
  USE FLDINX
  USE DRCATC     
  USE DRLINK  
  USE DRSEDT
  USE DRSPFL 
  USE DRWFLD
  USE NUMBRS  
  USE PLZONE 
  USE PROPER
  USE SOLVAR 
  USE SPECIE 
  USE LANDSF
!------------------------------------------------------------------------------C
  IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
!------------------------------------------------------------------------------C
!------------------------------------------------------------------------------C
!     Compute thickness of upper bottom deposit layer and radionuclide
!     concentration in it.
!------------------------------------------------------------------------------C
 RSB = 0
 do 200 L = 1,NLINK
   do 100 I = 1,IL(L)
     N = ND(L,I)
     IZN = IZ(N)
!          SP  = SIGMA*SH(N)  + (1.d0 -SIGMA)*SH0(N)
     HP  = SIGMA*HR(N)  + (1.d0 -SIGMA)*HRO(N)
     CP  = SIGMA*CH(N)  + (1.d0 -SIGMA)*CHO(N)
     CSP = SIGMA*CSH(N)/max(SH(N),SMALL) + (1.d0 -SIGMA)*CSHO(N)/max(SHO(N),SMALL)
!------------------------------------------------------------------------------C
!     Thickness of upper bottom deposit.
!------------------------------------------------------------------------------C
     QS = SDOWN(N)/(1.d0 -POR(IZN))
     QB = RSUP(N)/(1.d0 -POR(IZN))/RHOS(IZN)
!------------------------------------------------------------------------------C
!     Exchangeable phase radionuclide activity of upper soil layer.
!------------------------------------------------------------------------------C
     ZP  = SIGMA*ZSD(N) + (1.d0 -SIGMA)*ZSDO(N)
     ZP  = max(0.d0,ZP)
     ERH = HP/max(HP,SMALL)
     if (ISOLVE(3)+ISOLVE(4) > 0) then
       PCB = EXSB(IZN)*PCSB(IZN)*RHOS(IZN)/RHOLQ
       ACB = ZP*(1.d0/DT/SIGMA + ERH*EXSB(IZN) + HFLF(IZN)) + QB
       FCB = max(0.d0,ZSDO(N))*CBO(N)/DT/SIGMA + ERH*PCB*ZP*CP + QS*CSP
       if (ACB <= SMALL) then
         CBS = 0.d0
       else
         CBS = FCB/ACB
       endif

       CNEW = (CBS - (1.d0 -SIGMA)*CBO(N))/SIGMA
       RES = ABS(CBS -CB(N))/max(ABS(CBS),SMALL)
       CB(N) = (CBS - (1.d0 -SIGMA)*CBO(N))/SIGMA
       RSB = max(RSB,RES)
     endif
100  CONTINUE
200 CONTINUE
  IRET = 0
!  if (RSB > RSDMX(4)) then
    IRET = 1
!  endif
!------------------------------------------------------------------------------C
!     End of UPLAYR group. 
!----------------------------------------------------------------------C 
! 
  RETURN
    END
!==============================================================================C 

!==============================================================================C 
SUBROUTINE UPLAYR
!------------------------------------------------------------------------------C
!
!     UPLAYR: transport equations for UPper LAYeR of bottom deposit.
!
!------------------------------------------------------------------------------C
  USE CONSTS 
  USE FLDINX
  USE DRCATC     
  USE DRLINK  
  USE DRSEDT
  USE DRSPFL 
  USE DRWFLD
  USE NUMBRS  
  USE PLZONE 
  USE PROPER
  USE SOLVAR 
  USE SPECIE 
  USE LANDSF 
  USE POINTS
!  USE SOLVERSED
!------------------------------------------------------------------------------C
IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
!------------------------------------------------------------------------------C
REAL*8 aCf(3,4), Kds,Kdb,Lambda,As,Ab
!------------------------------------------------------------------------------C 
 if( iSolve(2) +iSolve(3) +iSolve(4) == 0 ) return
!------------------------------------------------------------------------------C
!     At this moment HR/=HR0, SH/=SH0, CH/=CH0, CSH/=CSH0,
!     only CB = CB0.
!------------------------------------------------------------------------------C
 do L=1,NLINK  
   do i=1,IL(L)   
     ic = nd(L,i)  
     IPLANT = IZPL(ic)
     IZN = IZ(ic)

     aCf = 0.d0
!------------------------------------------------------------------------------C
!     Calculate Qs (deposition) and Qb (erosion).
!     Use already updated SH for simplicity.
!------------------------------------------------------------------------------C
     Qs = 0.d0;  Qb = 0.d0; 
     if( ISOLVE(2) == 0 ) goto 100
     if (ILOAD/=9) then !--Non-cohesive
       tErod = EROD(IPLANT); 
       tt = SEQ(ic) -SH(ic)
       if (tt>=0.d0) then
         Qb = VSD*tErod*tt !--Erosion
       else
         Qs = -VSD*tt !--Deposition
       endif
     else !--Cohesive
       UP = sqrt(2.d0)*HR(ic)*QH(ic)/sqrt(HR(ic)**4 +max(HR(ic)**4,small))
       TAY = UP*UP*RHOLQ
       if (TAY < TAY_D) Qs = VSD*(1.d0 -TAY/TAY_D) !--Deposition
       if (TAY > TAY_E) Qb = WCM*(TAY/TAY_E -1.d0) !--Erosion
     endif
100  continue
!------------------------------------------------------------------------------C
!     Calculate params.
!------------------------------------------------------------------------------C
     zp = (1.d0 -POR(IZN))*ZSD(ic)
     Kds  = PCSD(IZN)/RHOLQ
     Kdb  = PCSB(IZN)*RHOS(IZN)/RHOLQ
     As = EXSD(IZN)
     Ab = EXSB(IZN)
     Lambda = HFLF(IZN)

!------------------------------------------------------------------------------C
!     Fill equation coefs.
!------------------------------------------------------------------------------C 
     if( ISOLVE(2) == 0 ) then  
       SHic = 0.d0;  
     else 
       SHic = SH(ic)  
     endif    
     if (HR(ic) >= HSMALL .or. AW(ic) >= HSMALL ) then
!------------------------------------------------------------------------------C
!     Fill sediments in aqueous phase coefs.
!------------------------------------------------------------------------------C
!       aCf(1,1) = HR(ic)*(1.d0/DT + Lambda + As*Kds*SH(ic)/RHOSE) + zp*Ab*Kdb
       aCf(1,1) = AW(ic)*(1.d0/DT + Lambda + As*Kds*SHic) + zp*Ab*Kdb*BW(ic)
       aCf(1,2) = -As*AW(ic)
       aCf(1,3) = -zp*Ab*BW(ic)
       aCf(1,4)  = AW(ic)*CH(ic)/DT
!------------------------------------------------------------------------------C
!     Fill suspended sediments coefs.
!------------------------------------------------------------------------------C
!       aCf(2,1) = -HR(ic)*As*Kds*SH(ic)/RHOSE
!       aCf(2,2) = HR(ic)*(1.d0/DT +Lambda +As) +Qs/max(small,SH(ic))
       aCf(2,1) = -AW(ic)*As*Kds*SHic
       aCf(2,2) = AW(ic)*(1.d0/DT +Lambda +As) +BW(ic)*Qs/max(small,SHic)
       aCf(2,3) = -BW(ic)*Qb/RHOS(IZN)
       aCf(2,4) = AW(ic)*CSH(ic)/DT
     else
       aCf(1,1) = 1.d0
       aCf(2,2) = 1.d0
     endif
!------------------------------------------------------------------------------C
!     Fill bottom sediments coefs.
!------------------------------------------------------------------------------C
     if( zp >= SMALL ) then
       aCf(3,1) = -zp*Ab*Kdb
!       aCf(3,2) = -Qs/max(small,SH(ic))
       aCf(3,2) = -Qs/max(small,SHic)
       aCf(3,3) = zp*(1.0/DT +Lambda +Ab) +Qb/RHOS(IZN)
       aCf(3,4) = (1.d0-POR(IZN))*ZSD(ic)*CB(ic)/DT
     else
       aCf(3,3) = 1.d0
     endif

!------------------------------------------------------------------------------C
!     Solve 3 species equations.
!------------------------------------------------------------------------------C
     call Solve3x3(aCf, CH(ic),CSH(ic),CB(ic))
   enddo
 enddo
!------------------------------------------------------------------------------C
!     End of UPLAYR group. 
!----------------------------------------------------------------------C 
! 
RETURN
END
!==============================================================================C 

!==============================================================================C 
SUBROUTINE Solve3x3 (aCf, x1,x2,x3)
!------------------------------------------------------------------------------C
!
!     Solve3x3: Solution of 3x3 linear equations system.
!
!------------------------------------------------------------------------------C
!USE GLOBAL
!------------------------------------------------------------------------------C
IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
!------------------------------------------------------------------------------C
REAL*8 aCf(3,4), x1,x2,x3
!------------------------------------------------------------------------------C
!------------------------------------------------------------------------------C
!     Calculate determinants.
!------------------------------------------------------------------------------C
  d  = aCf(1,1) *(aCf(2,2)*aCf(3,3) -aCf(2,3)*aCf(3,2))     &
      -aCf(2,1) *(aCf(1,2)*aCf(3,3) -aCf(1,3)*aCf(3,2))     &
      +aCf(3,1) *(aCf(1,2)*aCf(2,3) -aCf(1,3)*aCf(2,2))

  d1 = aCf(1,4) *(aCf(2,2)*aCf(3,3) -aCf(2,3)*aCf(3,2))     &
      -aCf(2,4) *(aCf(1,2)*aCf(3,3) -aCf(1,3)*aCf(3,2))     &
      +aCf(3,4) *(aCf(1,2)*aCf(2,3) -aCf(1,3)*aCf(2,2))

  d2 = aCf(1,1) *(aCf(2,4)*aCf(3,3) -aCf(2,3)*aCf(3,4))     &
      -aCf(2,1) *(aCf(1,4)*aCf(3,3) -aCf(1,3)*aCf(3,4))     &
      +aCf(3,1) *(aCf(1,4)*aCf(2,3) -aCf(1,3)*aCf(2,4))

  d3 = aCf(1,1) *(aCf(2,2)*aCf(3,4) -aCf(2,4)*aCf(3,2))     &
      -aCf(2,1) *(aCf(1,2)*aCf(3,4) -aCf(1,4)*aCf(3,2))     &
      +aCf(3,1) *(aCf(1,2)*aCf(2,4) -aCf(1,4)*aCf(2,2))
!------------------------------------------------------------------------------C
!     Calculate solution.
!------------------------------------------------------------------------------C
  x1 = d1/d; x2 = d2/d; x3 = d3/d;
!------------------------------------------------------------------------------C
!     End of Solve3x3 group. 
!----------------------------------------------------------------------C 
! 
RETURN
END
!==============================================================================C 

    
