!= =====================================================================C 
C 
      SUBROUTINE RDLITO 
C 
C----------------------------------------------------------------------C 
C 
C     RDLITO: ReaD LInk TOpology data input group. 
C 
C----------------------------------------------------------------------C  
      USE FILINX 
      USE FLDINX 
      USE DRLINK 
      USE DRNAME 
      USE DRWFLD  
      USE DRSPFL
      USE LANDSF 
      USE LINEQN 
      USE PLZONE
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
      CHARACTER*240 CHDUM 
      CHARACTER*120 NOTES 
      CHARACTER*80  ATYPE  
      INTEGER*4 NPZN(4) 
      
      CHARACTER*2 is1, is2, is3, is4, is5, is6
      
C----------------------------------------------------------------------C 
C     WRITE header to output file. 
C----------------------------------------------------------------------C 
      WRITE (IWR,'(/A)') ' Link Topological Data' 
      WRITE (IWR,'(A/)') ' ---------------------' 
C----------------------------------------------------------------------C 
C     Loop over link topology data records. 
C----------------------------------------------------------------------C 
      READ (IRD,*) NLIN    
C----------------------------------------------------------------------C 
C     Allocate memory for link variables. 
C----------------------------------------------------------------------C 
      ALLOCATE(IL(NLIN),MU(NLIN),MD(NLIN)); IL = 0; MU = 0; MD = 0; 
      ALLOCATE(LinkName(NLIN),NMU(NLIN),NMD(NLIN));   
      LinkName = ' '; NMU = ' '; NMD = ' '; 
      MaxLPoints = 0; mPoints = 0;
C----------------------------------------------------------------------C 
      IS = LOG10(REAL(NLIN)) 
      IF( IS .LE. 2 ) IS = 2  
      IS = 9 -IS   
      
      is1 = "  "; is2 = "  "; is3 = "  "; is4 = "  "; is5 = "  "
      is6 = "  "
      
      write(is1, '(I0)') 9-IS
      write(is2, '(I0)') 17-IS
      write(is3, '(I0)') 29-IS 
      write(is4, '(I0)') 39-IS 
      write(is5, '(I0)') 49-IS 
      write(is6, '(I0)') 8-IS
     
      WRITE(IWR, '(T' // TRIM(is1) // ',A,T' // TRIM(is2) // 
     &    ',A,T' // TRIM(is3) // ',A,T' // TRIM(is4) // 
     &    ',A,T' // TRIM(is5) // ',A)')     
     &            'No','LINK','IL','MU','MD'  !!,'LTYPE','MSEQ'
     
      WRITE(IWR, '(T' // TRIM(is6) // ',45(1H-))')
      WRITE(IWR, '(T' // TRIM(is6) //',64(1H-))')
         
!      WRITE (IWR,'(T<9-IS>,A,T<17-IS>,A,T<29-IS>,A,T<39-IS>,A,T<49-IS>, 
!     &       A,T<56-IS>,A,T<67-IS>,A)') 
!     &       'No','LINK','IL','MU','MD'  !!,'LTYPE','MSEQ'
!      WRITE (IWR,'(T<8-IS>,45(1H-))') 
      
!      WRITE (IWR,'(T<8-IS>,64(1H-))') 
C
      NLINK = 0 
      DO 600 NL = 1,NLIN 
C----------------------------------------------------------------------C 
C     Read a set of link topology data. 
C----------------------------------------------------------------------C 
        READ (IRD,'(A)') CHDUM 
        CALL LCASE( CHDUM ) 
C----------------------------------------------------------------------C 
C     Read sequential link number. 
C----------------------------------------------------------------------C 
        ISTART = 1 
        CALL RDINT(ISTART,ICOMMA,CHDUM,Lnum) 
C----------------------------------------------------------------------C 
C     Read link name. 
C----------------------------------------------------------------------C 
        ATYPE(1:) = ' ' 
        CALL RDCHR(ISTART,ICOMMA,CHDUM,ATYPE) 
        DO 300 M = 1,NLINK 
          IF( LinkName(M) .EQ. ATYPE ) THEN 
            ILINK = M 
            GOTO 301 
          ENDIF 
  300   CONTINUE 
C----------------------------------------------------------------------C 
C     If a previously undefined link name, increment the number 
C     of link name (NLINK) and check to ensure parameter limit 
C     is not exceeded. 
C----------------------------------------------------------------------C 
        NLINK = NLINK +1  
        if( nLink /= Lnum ) then  
           write(*,*) ' Warning: Incorrect sequential number of link ',
     &     Lnum,' under reading link topology'   
           write(IWR,*)' Warning: Incorrect sequential number of link ',
     &     Lnum,' under reading link topology'   
        endif    
        LinkName(NLINK) = ATYPE 
        ILINK = NLINK 
  301   CONTINUE 
C----------------------------------------------------------------------C 
C     Read link topology data. 
C----------------------------------------------------------------------C 
        CALL RDINT(ISTART,ICOMMA,CHDUM,IL(ILINK)) 
        CALL RDCHR(ISTART,ICOMMA,CHDUM,NMU(ILINK)) 
        CALL RDCHR(ISTART,ICOMMA,CHDUM,NMD(ILINK))  
        MaxLPoints = max(MaxLPoints,IL(iLink)) 
        mPoints = mPoints +IL(iLink)
C----------------------------------------------------------------------C 
C     Write link topology information to output file. 
C----------------------------------------------------------------------C 
        LN  = LEN_TRIM(LinkName(ILINK))    
        LN1 = LEN_TRIM(NMU(ILINK))    
        LN2 = LEN_TRIM(NMD(ILINK)) 
      
      is1 = '  '; is2 = '  '; is3 = '  '; is4 = '  '; is5 = '  '
      is6 = '  '  
        
      write(is1, '(I0)') 10-IS 
      write(is2, '(I0)') 21-LN-IS 
      write(is3, '(I0)') 21-IS 
      write(is4, '(I0)') 41-LN1-IS 
      write(is5, '(I0)') 51-LN2-IS 
      write(is6, '(I0)') 51-IS
     
      WRITE(IWR, '(I' // TRIM(is1) // ',T' // TRIM(is2) // 
     &    ',A,T' // TRIM(is3) // ',I10,T' // TRIM(is4) // 
     &    ',A,T' // TRIM(is5) // ',A,T' // TRIM(is6) // ',2I10)')
     &       ILINK,LinkName(ILINK),IL(ILINK),NMU(ILINK),
     &       NMD(ILINK)     
     
        
   !     WRITE (IWR,'(I<10-IS>,T<21-LN-IS>,A,T<21-IS>,I10,T<41-LN1-IS>,
   !  &         A,T<51-LN2-IS>,A,T<51-IS>,2I10)') 
   !  &         ILINK,LinkName(ILINK),IL(ILINK),NMU(ILINK),
   !  &         NMD(ILINK)
!     &         ,LTYPE(ILINK),MSEQ(ILINK)   

  550   CONTINUE 
  600 CONTINUE  
C----------------------------------------------------------------------C 
      ALLOCATE(ILK(mPoints),JPT(mPoints),IXP(mPoints)); 
      ILK = 0; JPT = 0; IXP = 0; 
      ALLOCATE(NSX(mPoints),NSY(mPoints),ND(NLIN,MaxLPoints));
      NSX = 0; NSY = 0; ND = 0;
C----------------------------------------------------------------------C 
C     Define the domain translator arrays. 
C----------------------------------------------------------------------C 
      NPOINT = 0 
	DO 710 I = 1,NLINK 
	  DO 700 J = 1,IL(I) 
	    NPOINT = NPOINT +1 
		ND(I,J) = NPOINT
		ILK(NPOINT) = I 
		JPT(NPOINT) = J  
  700   CONTINUE 
  710 CONTINUE   		      
C----------------------------------------------------------------------C 
C     Allocate memory for computational point variables. 
C----------------------------------------------------------------------C  
!----------------------------------------------------------------------C 
!     Primary variables. 
!----------------------------------------------------------------------C   
      ALLOCATE(X(nPoint),TSF(nPoint),TSFO(nPoint),ROU(nPoint))
      ALLOCATE(MSEC(nPoint),AW(nPoint),AWO(nPoint),FLQ(nPoint))
      ALLOCATE(PER(nPoint),PERO(nPoint),DAW(nPoint),DPER(nPoint)); 
      ALLOCATE(HR(nPoint),HRO(nPoint),QH(nPoint),QHO(nPoint))  
      allocate(BW(nPoint)); BW = 0.d0;
      if( iSolve(1) > 0 ) then
        ALLOCATE(YH(nPoint)); YH = 0.d0;  
        allocate(CRA(nPoint,10),CRH(nPoint,6)); 
        CRA = 0.d0; CRH = 0.d0;
      endif
        
      do i=1, nPoint
        X(i) = 0.0d0; TSF(i) = 0.d0; ROU(i) = 0.d0; MSEC(i) = 0.d0; 
        AW(i) = 0.d0; AWO(i) = 0.d0; DAW(i) = 0.d0; FLQ(i) = 0.0d0;
        PER(i) = 0.d0; PERO(i) = 0.d0; dPER(i) = 0.d0; 
        HR(i) = 0.d0; HRO(i) = 0.d0; QH(i) = 0.d0; QHO(i) = 0.d0;
      enddo 
      ALLOCATE(DH(NPoint),DQ(NPoint));  
      allocate(RH(NPoint),RQ(NPoint));  
!      
      if( iSolve(2) == 1 ) then  
        allocate(SH(NPoint),SHO(NPoint)) 
        allocate(SEQ(NPoint))  
        SH = 0.d0; SHO = 0.d0; SEQ = 0.d0; 
      endif  
!   
      if( iSolve(3) +iSolve(4) > 0 ) then  
        allocate(CH(NPoint),CHO(NPoint)) 
        CH = 0.d0; CHO = 0.d0;
        allocate(CSH(NPoint),CSHO(NPoint)) 
        CSH = 0.d0; CSHO = 0.d0;
        allocate(CB(NPoint),CBO(NPoint))  
        CB = 0.d0; CBO = 0.d0;
      endif  
!
      if( iSolve(2) +iSolve(3) +iSolve(4) > 0 ) then  
        allocate(ZSD(NPoint),ZSDO(NPoint))  
        ZSD = 0.d0; ZSDO = 0.d0;  
        allocate(IZ(NPoint),IZPL(NPoint)) 
        IZ = 0; IZPL = 0;
      endif  
!      
C----------------------------------------------------------------------C 
C     End of RDLITO group. 
C----------------------------------------------------------------------C 
C 
      RETURN 
      END 
 
!= =====================================================================C 
C 
      SUBROUTINE RDNDTO 
C 
C----------------------------------------------------------------------C 
C 
C     RDNDTO: ReaD NoDe TOpology data input group. 
C 
C----------------------------------------------------------------------C  
      USE FILINX 
      USE FLDINX 
      USE DRLINK 
      USE DRNAME 
      USE DRNODE 
      USE LINEQN 
      USE SOLVAR
C----------------------------------------------------------------------C 
C     Implicit Double Precision. 
C----------------------------------------------------------------------C 
      IMPLICIT REAL*8 (A-H,O-Z) 
      IMPLICIT INTEGER*4 (I-N) 
C----------------------------------------------------------------------C 
C     Type Declarations. 
C----------------------------------------------------------------------C 
      CHARACTER*240 CHDUM 
      CHARACTER*120 NOTES,FORM,FORM1 
      CHARACTER*80  ATYPE,STR  
      INTEGER*4 NPZN(4) 
      CHARACTER NodeLink*10     
      CHARACTER*6, POINTER :: tAr(:,:) 
      INTEGER*4, POINTER :: iTemp(:)
C----------------------------------------------------------------------C 
C     WRITE header to output file. 
C----------------------------------------------------------------------C 
      WRITE (IWR,'(/A)') ' Node Topological Data' 
      WRITE (IWR,'(A/)') ' ---------------------' 
C----------------------------------------------------------------------C 
C     Loop over node topology data records. 
C----------------------------------------------------------------------C 
      READ (IRD,*) NLIN 
C----------------------------------------------------------------------C 
C     Allocate memory for node variables. 
C----------------------------------------------------------------------C 
      ALLOCATE(NodeName(NLIN),KK(NLIN)); 
      NodeName = ' '; KK = 0; 
      ALLOCATE(NGRP(NLIN)); 
      NGRP = 1; ! NPOS = 0; 
      MaxLinks = 0;  
      lMaxL = 5; lGroup = 5; 
      ALLOCATE(NLK(NLIN,lMaxL)); NLK = ' '; 
      ALLOCATE(nGSize(lGroup)); nGSize = 0;
C----------------------------------------------------------------------C 
!      WRITE (IWR,'(T5,A,T10,A,T16,A,T24,A,T33,A,T38,A,T46,A,T54,A)') 
!    &       'No','NODE','NGROUP','NPOS','KK','NL(1)','NL(2)','.....'
!      WRITE (IWR,'(T4,65(1H-))') 
!------------------------------------------------------------------------------C 
!     Output head of the node topology table. 
!------------------------------------------------------------------------------C 
      WRITE (IWR,'(T5,A,T10,A,T19,A,T25,A,T33,A,T41,A)')          
     &        'No','NODE','KK','NL(1)','NL(2)','.....'
      WRITE (IWR,'(T4,65(1H-))') 
C
      NNODE = 0; NGroup = 0; nBound = 0;
      DO 400 NL = 1,NLIN 
C----------------------------------------------------------------------C 
C     Read a set of node topology data. 
C----------------------------------------------------------------------C 
        READ (IRD,'(A)') CHDUM 
        CALL LCASE( CHDUM ) 
C----------------------------------------------------------------------C 
C     Read node sequential number. 
C----------------------------------------------------------------------C 
        ISTART = 1 
        CALL RDINT(ISTART,ICOMMA,CHDUM,NDnum) 
C----------------------------------------------------------------------C 
C     Read node name. 
C----------------------------------------------------------------------C 
        ATYPE(1:) = ' ' 
        CALL RDCHR(ISTART,ICOMMA,CHDUM,ATYPE) 
        DO 300 M = 1,NNODE 
          IF( NodeName(M) .EQ. ATYPE ) THEN 
            INODE = M 
            GOTO 301 
          ENDIF 
  300   CONTINUE 
C----------------------------------------------------------------------C 
C     If a previously undefined node name, increment the number 
C     of node name (NNODE) and check to ensure parameter limit 
C     is not exceeded. 
C----------------------------------------------------------------------C 
        NNODE = NNODE +1 
        if( nNode /= NDnum ) then  
           write(*,*) ' Warning: Incorrect sequential number of node ',
     &     NDnum,' under reading node topology'   
           write(IWR,*)' Warning: Incorrect sequential number of node ',
     &     NDnum,' under reading node topology'   
        endif    
        NodeName(NNODE) = ATYPE 
        INODE = NNODE 
  301   CONTINUE 
C----------------------------------------------------------------------C 
C     Read node topology data. 
C----------------------------------------------------------------------C 
!        CALL RDINT(ISTART,ICOMMA,CHDUM,NGRP(INODE)) 
!        CALL RDINT(ISTART,ICOMMA,CHDUM,NPOS(INODE)) 
        CALL RDINT(ISTART,ICOMMA,CHDUM,KKK) 
        MaxLinks = max(MaxLinks,KKK) 
!        NGroup = max(NGroup,NGRP(INODE))
c        NGROUP = MAX(NGROUP,NGRP(INODE)) 
        if( KKK == 1 ) nBound = nBound +1
C
        IF( KKK .GT. lMaxL ) THEN 
          ALLOCATE( tAr(NLIN,lMaxL) ); tAr = NLK; 
          DEALLOCATE(NLK);  
          ALLOCATE(NLK(NLIN,lMaxL+5)); NLK(1:NLIN,1:lMaxL) = tAr;
          lMaxL = lMaxL +5; DEALLOCATE(tAr);
        ENDIF
C
        DO K=1,KKK  
          CALL RDCHR(ISTART,ICOMMA,CHDUM,NLK(INODE,K)) 
        ENDDO   
        KK(INODE) = KKK 
C
        IF( NGRP(INODE) .GT. lGroup ) THEN 
          ALLOCATE(iTemp(lGroup)); iTemp = nGSize; 
          DEALLOCATE(nGSize); ALLOCATE(nGSize(lGroup+5)); 
          nGSize(1:lGroup) = iTemp; DEALLOCATE(iTemp); 
          lGroup = lGroup+5;
        ENDIF
        nGSize(NGRP(INODE)) = nGSize(NGRP(INODE)) +1  
C----------------------------------------------------------------------C 
C     Write link topology information to output file. 
C----------------------------------------------------------------------C 
        LN = LEN_TRIM(NodeName(INODE))  
        FORM = '(T2,I5,T  ,A,T16,I5'  
        WRITE(FORM(9:10),'(I2)') 14-LN  
        DO K=1,KKK 
          STR = ',T'  
          LP = 22+8*K -LEN_TRIM(NLK(INODE,K))  
          IF( LP .LT. 100 ) THEN 
            WRITE(STR(3:4),'(I2)') LP
          ELSE    
            WRITE(STR(3:5),'(I3)') LP
          ENDIF 
          STR = TRIM(STR) // ',A'   
          FORM = TRIM(FORM) // TRIM(STR)    
        ENDDO  
        FORM = TRIM(FORM) // ')'  
        WRITE (IWR,TRIM(FORM)) 
     &       INODE,NodeName(INODE),KK(INODE),
     &       (TRIM(NLK(INODE,K)),K=1,KKK)    
  400 CONTINUE  
C----------------------------------------------------------------------C 
      ALLOCATE(LK(NNODE,MaxLinks)); LK = 0;
      ALLOCATE( tAr(NNODE,MaxLinks) ); tAr = NLK(1:NNODE,1:MaxLinks); 
      DEALLOCATE(NLK); ALLOCATE(NLK(NNODE,MaxLinks)); NLK = tAr;
      DEALLOCATE(tAr);
C----------------------------------------------------------------------C
C     Verify that all MU,MD correspond to defined nodes and 
C     represent legal node group combinations. 
C----------------------------------------------------------------------C 
      DO 560 L=1,NLINK
        MUNODE = 0
        MDNODE = 0
        DO 510 I=1,NNODE
          IF( NMU(L) .EQ. NodeName(I) ) MUNODE = I
          IF( NMD(L) .EQ. NodeName(I) ) MDNODE = I
  510   CONTINUE
C
        IF( MUNODE*MDNODE .EQ. 0 ) THEN  
          WRITE (IWR,9012) LinkName(L) 
          WRITE (IWR,'(/2A)') ' INPUT OR COMPILATION ERROR: ', 
     &      'Node Topology Data Input.' 
          WRITE (*,'(/2A)') ' INPUT OR COMPILATION ERROR: ', 
     &      'Node Topology Data Input (check output file).' 
          STOP 
        ENDIF 
           
        IF( MUNODE*MDNODE .NE. 0 .AND. 
     &      IABS(NGRP(MUNODE)-NGRP(MDNODE)) .GT. 1 ) THEN
          WRITE (IWR,9013) LinkName(L) 
          WRITE (IWR,'(/2A)') ' INPUT OR COMPILATION ERROR: ', 
     &      'Node Topology Data Input.' 
          WRITE (*,'(/2A)') ' INPUT OR COMPILATION ERROR: ', 
     &      'Node Topology Data Input (check output file).' 
          STOP 
        ENDIF                         

        IF( MUNODE .LE. 0 ) GOTO 530
        KKK = KK(MUNODE)
        DO 520 K=1,KKK 
          IF( NLK(MUNODE,K)(1:1) .EQ. '-' ) THEN 
            NodeLink = NLK(MUNODE,K)(2:)  
          ELSE  
            NodeLink = NLK(MUNODE,K)  
          ENDIF     
C          IF( LinkName(L) .EQ. NLK(MUNODE,K) ) GOTO 530
          IF( LinkName(L) .EQ. NodeLink ) GOTO 530
  520   CONTINUE
C     Error message.
        WRITE (IWR,9014) LinkName(L),MUNODE 
        WRITE (IWR,'(/2A)') ' INPUT OR COMPILATION ERROR: ', 
     &    'Node Topology Data Input.' 
        WRITE (*,'(/2A)') ' INPUT OR COMPILATION ERROR: ', 
     &    'Node Topology Data Input (check output file).' 
        STOP 
  530   IF( MDNODE .LE. 0 ) GOTO 560
C
        KKK = KK(MDNODE)
        DO 540 K=1,KKK
          IF( NLK(MDNODE,K)(1:1) .EQ. '-' ) THEN 
            NodeLink = NLK(MDNODE,K)(2:)  
          ELSE  
            NodeLink = NLK(MDNODE,K)  
          ENDIF     
          IF( LinkName(L) .EQ. NodeLink ) GOTO 560
  540   CONTINUE
C     Error message.
        WRITE (IWR,9014) LinkName(L),MDNODE 
        WRITE (IWR,'(/2A)') ' INPUT OR COMPILATION ERROR: ', 
     &    'Node Topology Data Input.' 
        WRITE (*,'(/2A)') ' INPUT OR COMPILATION ERROR: ', 
     &    'Node Topology Data Input (check output file).' 
        STOP 
  560 CONTINUE
C----------------------------------------------------------------------C
C     Verify that all node definitions refer to valid associated links. 
C----------------------------------------------------------------------C 
      DO 640 I=1,NNODE
        KKK = KK(I)
        DO 630 K=1,KKK
          IF( NLK(I,K)(1:1) .EQ. '-' ) THEN 
            NodeLink = NLK(I,K)(2:)  
          ELSE  
            NodeLink = NLK(I,K)  
          ENDIF     
          DO 610 L=1,NLINK
            IF( NodeLink .EQ. LinkName(L) ) GOTO 620
  610     CONTINUE
C     Error message.
          WRITE (IWR,9015) NodeName(I),NodeLink 
          WRITE (IWR,'(/2A)') ' INPUT OR COMPILATION ERROR: ', 
     &      'Node Topology Data Input.' 
          WRITE (*,'(/2A)') ' INPUT OR COMPILATION ERROR: ', 
     &      'Node Topology Data Input (check output file).' 
          STOP 
          GOTO 630
  620     IF( NodeName(I) .NE. NMU(L).AND.NodeName(I) .NE. NMD(L) ) THEN 
            WRITE (IWR,9016) NodeName(I),NLK(I,K) 
            WRITE (IWR,'(/2A)') ' INPUT OR COMPILATION ERROR: ', 
     &        'Node Topology Data Input.' 
            WRITE (*,'(/2A)') ' INPUT OR COMPILATION ERROR: ', 
     &        'Node Topology Data Input (check output file).' 
            STOP 
          ENDIF   
          IF( (NLK(I,K)(1:1) .EQ. '-'.AND. NodeName(I) .NE. NMD(L)) .OR.
     &      (NLK(I,K)(1:1) .NE. '-'.AND. NodeName(I) .NE. NMU(L)) ) THEN 
            WRITE (IWR,9016) NodeName(I),NLK(I,K) 
            WRITE (IWR,'(/2A)') ' INPUT OR COMPILATION ERROR: ', 
     &        'Node Topology Data Input.' 
            WRITE (*,'(/2A)') ' INPUT OR COMPILATION ERROR: ', 
     &        'Node Topology Data Input (check output file).' 
            STOP 
          ENDIF  
  630   CONTINUE
  640 CONTINUE
C----------------------------------------------------------------------C
C     Replace references to node names by references to node sequence. 
C----------------------------------------------------------------------C   
      nALC = 0
      DO 660 L=1,NLINK
        DO 650 I=1,NNODE
          IF( NMU(L) .EQ. NodeName(I) ) MUNODE = I
          IF( NMD(L) .EQ. NodeName(I) ) MDNODE = I
  650   CONTINUE
        MU(L)=MUNODE
        MD(L)=MDNODE 
        nALC = nALC +3*IL(L) 
  660 CONTINUE
C----------------------------------------------------------------------C
C     Replace references to link names by references to link sequence. 
C----------------------------------------------------------------------C 
      DO 680 I=1,NNODE
        KKK = KK(I) 
        nALC = nALC +KKK*(KKK-1) 
        DO 680 K=1,KKK
          DO 670 L=1,NLINK
            IF( NLK(I,K)(1:1) .EQ. '-' ) THEN 
              NodeLink = NLK(I,K)(2:) 
              LSG = -1   
            ELSE  
              NodeLink = NLK(I,K) 
              LSG = 1   
            ENDIF     
            IF( NodeLink .EQ. LinkName(L) ) LK(I,K) = L*LSG
  670     CONTINUE
  680 CONTINUE
!----------------------------------------------------------------------C 
!     Lower-Upper decomposition initializations. 
!----------------------------------------------------------------------C  
      if( iSolve(1) == 1 ) then
        allocate(ALU(NNODE,NNODE),BLU(NNODE))
        DO 750 N = 1,NNODE 
          DO 740 M = 1,NNODE 
            ALU(M,N) = 0.D+0 
  740     CONTINUE 
          BLU(N) = 0.D+0 
  750   CONTINUE   
      endif  
      if( iSolve(2) +iSolve(3) +iSolve(4) > 0 ) then  
        allocate(ALC(nALC),BLC(nPoint),iEq(nPoint)) 
        allocate(IA(nALC),JA(nALC))
        nEq = 1  
        DO L=1,NLINK  
          DO I=1,IL(L)    
            N = ND(L,I)    
!            nEq = nEq +1
            iEq(N) = nEq   
            IF( I == 1 ) THEN  
              inode = mu(L)  
              nEq = nEq +KK(inode) +2
            ELSEIF( I == IL(L) ) THEN  
              inode = md(L)  
              nEq = nEq +KK(inode) +2
            ELSE   
!              iEq(N) = nEq 
              nEq = nEq +3
            ENDIF    
          ENDDO  
        ENDDO  
      endif    
C----------------------------------------------------------------------C 
C     Format Statements. 
C----------------------------------------------------------------------C 
 9012 FORMAT(/' Input or Compilation Error: Incorrect definition of MU', 
     &        ',MD for link name [',A,'].') 
 9013 FORMAT(/' Input or Compilation Error: Illegal connection between', 
     &        ' node groups as defined',/,' by MU,MD for link name [',
     &        A,'].') 
 9014 FORMAT(/' Input or Compilation Error: Link [',A,'] has an MU or', 
     &        ' MD [',I6,'] which is',/,' inconsistent with the ',
     &        'corresponding node definitions.') 
 9015 FORMAT(/' Input or Compilation Error: Definition of node [',A, 
     &        '] refers to a link name [',A,'] which is ',/,
     &        'not defined.') 
 9016 FORMAT(/' Input or Compilation Error: Definition of node [',A, 
     &        '] refers to a link name [',A,'] neither of whose ',/,
     &        'MU or MD refer to this node.') 
C----------------------------------------------------------------------C 
C     End of RDNDTO group. 
C----------------------------------------------------------------------C 
C 
      RETURN 
      END 
 
!= =====================================================================C 
C 
      SUBROUTINE RDPOIN 
C 
C----------------------------------------------------------------------C 
C 
C     RDPOIN: ReaD POINts data input group. 
C 
C----------------------------------------------------------------------C  
      USE FILINX 
      USE FLDINX 
      USE CRSECT 
      USE DRLINK 
      USE DRNAME 
      USE DRNODE 
      USE DRSOLV 
      USE DRWFLD  
      USE LANDSF
      USE POINTS
C----------------------------------------------------------------------C 
C     Implicit Double Precision. 
C----------------------------------------------------------------------C 
      IMPLICIT REAL*8 (A-H,O-Z) 
      IMPLICIT INTEGER*4 (I-N) 
C----------------------------------------------------------------------C 
C     Type Declarations. 
C----------------------------------------------------------------------C 
      CHARACTER*240 CHDUM 
      CHARACTER*120 NOTES 
      CHARACTER*10  ATYPE,LINKNM 
      CHARACTER*80  UNITS 
      CHARACTER*2 is1,is2
      INTEGER*4 NPZN(4) 
C----------------------------------------------------------------------C 
C     WRITE header to output file. 
C----------------------------------------------------------------------C 
      WRITE (IWR,'(/A)') ' Computational Points Data' 
      WRITE (IWR,'(A/)') ' -------------------------' 
C----------------------------------------------------------------------C 
C     Loop over link topology data records. 
C----------------------------------------------------------------------C 
      WRITE (IWR,'(T02,A,T9,A,T17,A,2X,A,3X,A,6X,A,6X,A,
     &       3X,A,2X,A,2X,A,4X,A,4X,A)') 
     &  'LINK','Pt.','X(km)','TALW(m)','RC','SEC','Y(m)','Q(m^3/s)',     
     &  ' APHA',' BETA','FX','SED'   
      WRITE (IWR,'(T2,88(1H-))')  
C
      READ (IRD,*) NLIN  
      NPOIN = 0; eps_mach = sqrt(1.d-17); 
      DO 100 NR = 1,NLIN 
C----------------------------------------------------------------------C 
C     Read a set of link topology data. 
C----------------------------------------------------------------------C 
        READ (IRD,'(A)') CHDUM 
        CALL LCASE( CHDUM ) 
C----------------------------------------------------------------------C 
C     Read link name. 
C----------------------------------------------------------------------C 
        ISTART = 1 
        LINKNM(1:) = ' ' 
        CALL RDCHR(ISTART,ICOMMA,CHDUM,LINKNM)  
C----------------------------------------------------------------------C 
C     Define number of link.
C----------------------------------------------------------------------C 
        DO 50 L=1,NLINK
          IF( LINKNM .EQ. LinkName(L) ) GOTO 60
   50   CONTINUE
        WRITE (IWR,9008) LINKNM,IPT 
        WRITE (IWR,'(/2A)') ' INPUT OR COMPILATION ERROR: ', 
     &    'Computational Points Data Input.' 
        WRITE (*,'(/2A)') ' INPUT OR COMPILATION ERROR: ', 
     &    'Computational Points Data Input (check output file).'  
        STOP  
C----------------------------------------------------------------------C 
C     Read computational points data. 
C----------------------------------------------------------------------C 
   60   CALL RDINT(ISTART,ICOMMA,CHDUM,IPT) 
        NP = ND(L,IPT)  
C     Read point distance.
        CALL RDDPR(ISTART,ICOMMA,CHDUM,X(NP)) 
        CALL RDCHR(ISTART,ICOMMA,CHDUM,UNITS) 
        CALL RDUNIT(UNITS,X(NP)) 
C     Read elevation of lowest point in cross-section.
        CALL RDDPR(ISTART,ICOMMA,CHDUM,TSF(NP)) 
        CALL RDCHR(ISTART,ICOMMA,CHDUM,UNITS) 
        CALL RDUNIT(UNITS,TSF(NP))  
C     Read roughness coefficient.          
        CALL RDDPR(ISTART,ICOMMA,CHDUM,ROU(NP)) 
C        CALL RDCHR(ISTART,ICOMMA,CHDUM,UNITS) 
C        CALL RDUNIT(UNITS,ROU(NP))  
        IF( ROU(NP) .LE. 0.D+0 ) ROU(NP) = rFRIC   
!      ROU(NP) =1.d-5;  
C     Read number of cross-section.           
C        CALL RDINT(ISTART,ICOMMA,CHDUM,ISEC) 
        CALL RDCHR(ISTART,ICOMMA,CHDUM,ATYPE)  
        DO 70 I = 1,NSEC
          IF( ATYPE .EQ. TSEC(I) ) THEN
            ISEC = I
            GOTO 80
          ENDIF
   70   CONTINUE
        WRITE (*,'(1A)') ' No cross-sectional table!'
        WRITE (*,'(2A)') ' Table we''re looking for is: ',ATYPE
        WRITE (*,'(1A)') ' Run aborting...'
        STOP
   80   CONTINUE 
        MSEC(NP) = ISEC 
C     Read water surface elevation.           
        CALL RDDPR(ISTART,ICOMMA,CHDUM,HR(NP)) 
        CALL RDCHR(ISTART,ICOMMA,CHDUM,UNITS) 
        CALL RDUNIT(UNITS,HR(NP)) 
C     Read water discharge.           
        CALL RDDPR(ISTART,ICOMMA,CHDUM,QH(NP)) 
        CALL RDCHR(ISTART,ICOMMA,CHDUM,UNITS) 
        CALL RDUNIT(UNITS,QH(NP))  
C     Read Alpha and Beta.          
        CALL RDDPR(ISTART,ICOMMA,CHDUM,ALPHAP) 
        CALL RDDPR(ISTART,ICOMMA,CHDUM,BETAP)  
!        IF( ALPHAP(NP) .LE. 0.D+0 ) ALPHAP(NP) = ALPHA 
!        IF( BETAP(NP) .LE. 0.D+0 ) BETAP(NP) = BETA 
C----------------------------------------------------------------------C 
C     Write points data.        
C----------------------------------------------------------------------C  
        L1 = LEN_TRIM(LINKNM)  
        L2 = LEN_TRIM(ATYPE)  
       
      is1 = '  '; is2 = '  '  
        
      write(is1, '(I0)') 6-L1 
      write(is2, '(I0)') 48-L2 
     
      WRITE(IWR, '(T' // TRIM(is1) // ',A,I5,F11.3,F8.2,F6.2,T' 
     &    // TRIM(is2) //',A,F7.2,F10.3,F8.2,F7.2)')  
     &        TRIM(LINKNM),IPT,X(NP)/1.D+3,TSF(NP),ROU(NP),
     &        TRIM(ATYPE),HR(NP),QH(NP),ALPHAP,BETAP
C     &         ,NSED(I,L)
        HR(NP) = MAX(0.D+0,HR(NP)-TSF(NP))  
  100 CONTINUE         
C----------------------------------------------------------------------C 
C     Check water surface level at node.        
C----------------------------------------------------------------------C 
      DO 350 I=1,NNODE
        KKK = KK(I)
        YSUM = 0.0
        DO 310 K=1,KKK
          L = IABS(LK(I,K))
C          J = IL(L)
C          IF( LK(I,K) .LT. 0) J = 1  
          J = 1
          IF( LK(I,K) .LT. 0) J = IL(L)  
C          
          NP = ND(L,J)          
          IF( K .EQ. 1) YYY = HR(NP) +TSF(NP)
          IF( abs(YYY-HR(NP)-TSF(NP)) > eps_mach*abs(YYY) ) THEN   
            WRITE (IWR,9009) NodeName(I),LinkName(L) 
            WRITE (IWR,'(/2A)') ' INPUT OR COMPILATION ERROR: ', 
     &        'Computational Points Data Input.' 
            WRITE (*,'(/2A)') ' INPUT OR COMPILATION ERROR: ', 
     &        'Computational Points Data Input (check output file).'  
          ENDIF  
          YSUM = YSUM +TSF(NP) +HR(NP)
  310   CONTINUE
        YAVG = YSUM/KKK
        DO 320 K=1,KKK
          L = IABS(LK(I,K))
C          J = IL(L)
C          IF( LK(I,K) .LT. 0 ) J=1    
          J = 1
          IF( LK(I,K) .LT. 0) J = IL(L)  
          NP = ND(L,J)          
          HR(NP) = YAVG -TSF(NP)
  320   CONTINUE
C----------------------------------------------------------------------C 
C     Detect the nodes with two links having the different 
C     cross-sections.        
C----------------------------------------------------------------------C 
        IF( KKK .NE. 2 ) GOTO 350
        LK1 = IABS(LK(I,1))
        LK2 = IABS(LK(I,2))
C        NSEC1 = NSEC(1,LK1)
C        IF( LK(I,1) .GT. 0 ) NSEC1 = NSEC(II(LK(I,1)),LK(I,1))
C        NSEC2 = NSEC(1,LK2)
C        IF( LK(I,2) .GT. 0 ) NSEC2 = NSEC(II(LK(I,2)),LK(I,2))
C        IF( LTYPE(LK1) .EQ. 0 .AND. LTYPE(LK2) .EQ. 0 .AND.
C     &      NSEC1 .NE. NSEC2 ) CALL ERRWAR(1,2,0,NODNAM(M),0,0.0)
  350 CONTINUE    

15    CONTINUE
C----------------------------------------------------------------------C 
C     Format Statements. 
C----------------------------------------------------------------------C 
 9008 FORMAT(/' Input or Compilation Error: A non-existent link [',A, 
     &        '] for ',/,'the point No [',I6,'].') 
 9009 FORMAT(/' Input or Compilation Error: Water surface elevation ',  
     &        'for link name [',A,'] differs ',/,'from water surface ',
     &        ' elevation at the node [',A,'].') 
C----------------------------------------------------------------------C 
C     End of RDPOIN group. 
C----------------------------------------------------------------------C 
C 
      RETURN 
      END 
 
!= =====================================================================C 
C 
      SUBROUTINE RDSECT1 
C 
C----------------------------------------------------------------------C 
C 
C     RDSECT: ReaD cross-SECTional data input group. 
C 
C----------------------------------------------------------------------C  
      USE FILINX 
      USE FLDINX 
      USE DRLINK 
      USE DRNAME 
      USE CRSECT
C----------------------------------------------------------------------C 
C     Implicit Double Precision. 
C----------------------------------------------------------------------C 
      IMPLICIT REAL*8 (A-H,O-Z) 
      IMPLICIT INTEGER*4 (I-N) 
C----------------------------------------------------------------------C 
C     Type Declarations. 
C----------------------------------------------------------------------C 
      CHARACTER*240 CHDUM 
      CHARACTER*120 NOTES 
      CHARACTER  ATYPE*10,PEROPT*80,UNITS*80 
      INTEGER*4 NPZN(4) 
      REAL*8, POINTER :: rTmp(:)
C----------------------------------------------------------------------C 
C     WRITE header to output file. 
C----------------------------------------------------------------------C 
      WRITE (IWR,'(/A)') ' Cross-Sectional Data' 
      WRITE (IWR,'(A)') ' --------------------' 
C----------------------------------------------------------------------C 
C     Loop over cross-sectional data records. 
C     Read record containing number of sets of cross-sectional data.
C     Set starting location for tables to zero.
C----------------------------------------------------------------------C 
      NSEC = 0  
	NCOUNT = 0 
      READ (IRD,*) NSETT
C
C----------------------------------------------------------------------C  
      ALLOCATE(TSEC(NSETT),NSCT(2,NSETT))  
      TSEC = ' '; NSCT = 0;
C
      DO 600 N = 1,NSETT 
C----------------------------------------------------------------------C 
C     Read record into character buffer and convert to lower case.
C----------------------------------------------------------------------C 
        READ (IRD,'(A)') CHDUM 
        CALL LCASE( CHDUM )
        ISTART = 1
        CALL RDCHR(ISTART,ICOMMA,CHDUM,ATYPE)
        DO 100 I = 1,NSEC
          IF( ATYPE .EQ. TSEC(I) ) THEN
            WRITE (*,'(2A)') ATYPE,
     &             '- Cross-section with such name already exists!'
            WRITE (IWR,'(2A)') ATYPE,
     &             '- Cross-section with such name already exists!'
            WRITE (*,'(1A)') ' Run aborting...'
            STOP
          ENDIF
  100   CONTINUE
C----------------------------------------------------------------------C 
        NSEC = NSEC +1   
        TSEC(NSEC) = ATYPE 
        CALL RDCHR(ISTART,ICOMMA,CHDUM,PEROPT)
C----------------------------------------------------------------------C
C     Tabular.
C----------------------------------------------------------------------C
        IF( PEROPT(1:4) .EQ. 'tabu' ) THEN
C----------------------------------------------------------------------C
C     Read tabular data for cross-section.
C----------------------------------------------------------------------C
          READ (IRD,*) NLIN 
          if( nCount == 0 ) then  
            ALLOCATE(SWD(NLIN),SHG(NLIN))  
            SWD = 0.d0; SHG = 0.d0;
          else  
            ALLOCATE(rTmp(nCount)); rTmp = SWD; 
            DEALLOCATE(SWD); ALLOCATE(SWD(nCount+NLIN)); 
            SWD(1:nCount) = rTmp; rTmp = SHG;
            DEALLOCATE(SHG); ALLOCATE(SHG(nCount+NLIN)); 
            SHG(1:nCount) = rTmp; nCount = nCount +NLIN; 
            DEALLOCATE(rTmp);
          endif    
          NCOUNT = NCOUNT +1
          NSCT(1,NSEC) = NCOUNT
          DO 300 I = 1,NLIN
            READ (IRD,'(A)') CHDUM
            CALL LCASE( CHDUM )
            ISTART = 1
C
            CALL RDDPR(ISTART,ICOMMA,CHDUM,RLX)
            CALL RDCHR(ISTART,ICOMMA,CHDUM,UNITS)
            CALL RDUNIT(UNITS,RLX)
C
            CALL RDDPR(ISTART,ICOMMA,CHDUM,RLY)
            CALL RDCHR(ISTART,ICOMMA,CHDUM,UNITS)
            CALL RDUNIT(UNITS,RLY)
C----------------------------------------------------------------------C
C     Sort table in ascending water depth order.
C----------------------------------------------------------------------C
            IF( I .GT. 1 ) THEN
              DO 130 J = NSCT(1,NSEC), NCOUNT -1
                IF( RLX .LE. SHG(J) ) THEN
                  DO 120 K = NCOUNT,J+1,-1
                    SHG(K) = SHG(K-1)
                    SWD(K) = SWD(K-1)
  120             CONTINUE
                  SHG(J) = RLX
                  SWD(J) = RLY
                  GOTO 140
                ENDIF
  130         CONTINUE
              SHG(NCOUNT) = RLX
              SWD(NCOUNT) = RLY
  140         CONTINUE
            ELSE
              SHG(NCOUNT) = RLX
              SWD(NCOUNT) = RLY
            ENDIF
            IF( I .NE. NLIN ) NCOUNT = NCOUNT +1
  300     CONTINUE
C
          WRITE (IWR,'(/2A)') ' Table name: ', ATYPE
          WRITE (IWR,'(2A,I6,A)') '   Cross-sectional Function: ',
     &        'Tabular with ',NLIN,' entries'
C
          NSCT(2,NSEC) = NCOUNT
C
          ISTART = NSCT(1,NSEC)
          IFIN   = NSCT(2,NSEC)
          ICNT   = 0
C
          DO 310 II = ISTART,IFIN
            ICNT = ICNT +1
            WRITE (IWR,'(A,I4,2(A,1PE11.4))') '   Entry No.',ICNT,
     &          '   Depth: ', SHG(II),
     &          '   Width: ', SWD(II)
  310     CONTINUE
C----------------------------------------------------------------------C
        ENDIF
C----------------------------------------------------------------------C
  600 CONTINUE
C----------------------------------------------------------------------C 
C     Format Statements. 
C----------------------------------------------------------------------C 
 9007 FORMAT(/' Input or Compilation Error: The number of table ',
     &        'entries specified [',I6,'] for input of the ',
     &        'cross-sectional name ', A,
     &        ' exceeds the maximum set in parameter LSCR [',I6,'].')
 9008 FORMAT(/' Input or Compilation Error: The number of table', 
     &        ' entries for cross-sectional data exceeds',/,
     &        ' the maximum set in', 
     &        ' parameter LSCT [',I6,'].') 
C----------------------------------------------------------------------C 
C     End of RDSECT group. 
C----------------------------------------------------------------------C 
C 
      RETURN 
      END 
 
!= =====================================================================C 