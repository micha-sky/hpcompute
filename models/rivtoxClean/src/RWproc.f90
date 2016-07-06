!==============================================================================!
SUBROUTINE ReadRange(delim,istart,str, iRange)
!------------------------------------------------------------------------------!
!     ReadRange: Read link name and points range values.
!------------------------------------------------------------------------------!
  USE FILINX 
  USE DRNAME  
  USE DRLINK
!------------------------------------------------------------------------------!
  IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
!------------------------------------------------------------------------------! 
!  interface   
!    integer*4 function RDStr(delim,istart,str, substr,ikeep)   
!      integer*4 istart 
!      character*(*) str
!      character*(*) substr 
!      character delim 
!      integer*4, optional :: ikeep
!    end function
!  end interface
  CHARACTER delim
  CHARACTER*(*) str
  INTEGER*4 istart
  INTEGER*4 iRange(4)  
  CHARACTER*80 Lname
!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
!     Read Link Name.
!------------------------------------------------------------------------------! 
  Lname = ' ';
!  call RDCHR(istart,icomma,str,Lname) 
  ires = RDStr(delim,istart,str, Lname)
  if( Lname == 'null' .or. Lname == 'all' ) then   
    iRange(1) = 0; iRange(2) = 0;  
    goto 301;
  endif    
  do m = 1,NLink 
    IF( LinkName(M) .EQ. Lname ) THEN 
      ILink = m  
      iRange(1) = m; iRange(2) = m;
      goto 301 
    ENDIF 
  enddo 
  WRITE (IWR,'(2a)')                                                &
      ' INPUT OR COMPILATION ERROR: There is not link with name ', Lname 
  WRITE (IWR,'(/2A)') ' INPUT OR COMPILATION ERROR: ',              &
       'Upper Bottom-Layer Data Input.' 
  WRITE (*,'(/2A)') ' INPUT OR COMPILATION ERROR: ',                & 
       'Upper Bottom-Layer Data Input (check output file).'  
  STOP  
 301 continue 
  !----------------------------------------------------------------------C
  !     Read 2 integers (range of upper layer bottom depth)
  !----------------------------------------------------------------------C
!  do i = 3,4
!    call RDINT(istart,icomma,str,iRange(i))
!  enddo
  
  do i = 3,4
    ires = RDI4(delim,istart,str, iRange(i))
  enddo
  !----------------------------------------------------------------------C
  !     Assign upper layer bottom depths.
  !----------------------------------------------------------------------C
 
!  call WholeRANGE(iRange)  !-----If all zeroes then whole zone-----

!------------------------------------------------------------------------------!
!     End of ReadRange group. 
!------------------------------------------------------------------------------!
  RETURN
END
!==============================================================================!

!==============================================================================!
SUBROUTINE WriteRange(iRange)
!------------------------------------------------------------------------------!
!     WriteRange: Write range to output.
!------------------------------------------------------------------------------!
USE FILINX   
USE FLDINX
USE DRLINK
IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
!------------------------------------------------------------------------------!
INTEGER*4 iRange(4)
!------------------------------------------------------------------------------!
  Link1 = max(1,iRange(1)); Link2 = min(nLink,iRange(2)); 
  if( Link2 == 0 ) Link2 = nLink; 
  iP1 = max(1,iRange(3)); iP2 = min(MaxLPoints,iRange(4));  
  if( iRange(4) == 0 ) iP2 = MaxLPoints; 
  if( Link1 == Link2 ) ip2 = min(IL(Link1),iRange(4))
  write(IWR,'(2(A,I4))') '   Domain: Links  = ',Link1,' to ',Link2
  write(IWR,'(2(A,I4))') '           Points = ',ip1,' to ',ip2
!------------------------------------------------------------------------------!
!     End of WriteRange group. 
!------------------------------------------------------------------------------!
RETURN
END
!==============================================================================!

!==============================================================================!
SUBROUTINE WholeRANGE(iRange)
!------------------------------------------------------------------------------!
!     WholeRANGE: Changes the Range to the whole domain if it was zero.
!------------------------------------------------------------------------------!
USE mGLOBAL
IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
!------------------------------------------------------------------------------!
INTEGER*4 iRange(4)
!------------------------------------------------------------------------------!
 !-----If all zeroes then whole zone-----
 if ((iRange(1)==0).and.(iRange(2)==0)) then
   iRange(1)=1
   iRange(2)=IFLD
 endif
 if ((iRange(3)==0).and.(iRange(4)==0)) then
   iRange(3)=1
   iRange(4)=JFLD
 endif
!------------------------------------------------------------------------------!
!     End of WholeRANGE group. 
!------------------------------------------------------------------------------!
RETURN
END
!==============================================================================!

!==============================================================================!
SUBROUTINE IsWholeRANGE(ARange,Flag)
!------------------------------------------------------------------------------!
!     IsWholeRANGE: Check if Range is the whole domain.
!------------------------------------------------------------------------------!
USE mGLOBAL
IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
!------------------------------------------------------------------------------!
INTEGER*4 ARange(4)
LOGICAL*4 Flag
!------------------------------------------------------------------------------!
 Flag=.false.
 if ((ARange(1)==1).and.(ARange(2)==IFLD).and.(ARange(3)==1).and.(ARange(4)==JFLD)) Flag=.true.
!------------------------------------------------------------------------------!
!     End of IsWholeRANGE group. 
!------------------------------------------------------------------------------!
  RETURN
END
!==============================================================================!

!==============================================================================!
SUBROUTINE AssignVariableInRegionI4(Value,V)
!------------------------------------------------------------------------------!
!     AssignVariableInRegionI4: Assign value to array v in region using
!       GRange global variable.
!       I4 variant.
!------------------------------------------------------------------------------!
  USE DRLINK
  USE FLDINX
!------------------------------------------------------------------------------!
  IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
!------------------------------------------------------------------------------!
  INTEGER*4 Value
  INTEGER*4, POINTER :: V(:)
!------------------------------------------------------------------------------!
  Link1 = max(1,GRange(1)); Link2 = min(nLink,GRange(2)); 
  if( Link2 == 0 ) Link2 = nLink; 
  do i = Link1,Link2 
    iP1 = max(1,GRange(3)); iP2 = min(IL(i),GRange(4));
    if( GRange(4) == 0 ) iP2 = IL(i); 
    do j = iP1,iP2
       N = ND(i,j)
       V(N) = Value
    enddo
  enddo
!------------------------------------------------------------------------------!
!     End of AssignVariableInRegionI4 group. 
!------------------------------------------------------------------------------!
  RETURN
END
!==============================================================================!

!==============================================================================!
SUBROUTINE AssignVariableInRegionR8(Value,V)
!------------------------------------------------------------------------------!
!     AssignVariableInRegionR8: Assign value to field array v in region using
!       GRange global variable.
!       R8 version.
!------------------------------------------------------------------------------! 
  USE FLDINX
  USE DRLINK
!------------------------------------------------------------------------------!
  IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
!------------------------------------------------------------------------------!
  REAL*8 Value,V(*)
!------------------------------------------------------------------------------! 
  Link1 = max(1,GRange(1)); Link2 = min(nLink,GRange(2)); 
  if( Link2 == 0 ) Link2 = nLink; 
  do i = Link1,Link2 
    iP1 = max(1,GRange(3)); iP2 = min(IL(i),GRange(4));
    if( GRange(4) == 0 ) iP2 = IL(i); 
    do j = iP1,iP2
      N = ND(i,j)
      V(N) = Value
    enddo
  enddo
!------------------------------------------------------------------------------!
!     End of AssignVariableInRegionR8 group. 
!------------------------------------------------------------------------------!
  RETURN
END
!==============================================================================!

!==============================================================================!
INTEGER*4 FUNCTION RDI4(delim,istart,str, Value) RESULT(iresult)
!------------------------------------------------------------------------------!
! RDI4: ReaD Integer*4.
! From string STR starting from istart position to the next DELIM.
! After EXIT istart points to the next position after the DELIM => usable for sequential read.
!------------------------------------------------------------------------------! 
  USE FILINX  
!------------------------------------------------------------------------------!
  IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
!------------------------------------------------------------------------------!
  CHARACTER delim
  INTEGER*4 istart
  CHARACTER*(*) Str
  INTEGER*4 Value
!------------------------------------------------------------------------------!
  INTEGER*4 icomma,ln
!------------------------------------------------------------------------------!
  iresult=0 !--OK
  Value=0
  ln=len_trim(Str)
  icomma = index(Str(istart:), delim) + istart -1
  if (icomma < istart) then
    if (istart<=ln) then !--Read evrth after the delim (non-blank)
      read(Str(istart:ln),*) Value
      istart=ln+1
    else !--Missing value
      call Msg(0,IWR,'INPUT WARNING! - Missing Integer Value - Zero assumed.')
      iresult=2 !--Missing value
    endif
  elseif (len_trim(Str(istart:icomma-1))==0) then !--If empty or blank string
    istart=icomma+1
    iresult=1 !--Empty or blank string
  else
!   read(Str(istart:icomma-1), '(I100)') Value !--For reading empty strings
    read(Str(istart:icomma-1), *) Value !--Can read " NN ..." strings
    istart=icomma+1
  endif
!------------------------------------------------------------------------------!
!     End of RDI4 group. 
!------------------------------------------------------------------------------!
  RETURN
END FUNCTION
!==============================================================================!

!==============================================================================!
INTEGER*4 FUNCTION RDR8(delim,istart,str, Value) RESULT(iresult)
!------------------------------------------------------------------------------!
! RDR8: ReaD Real*8.
! From string STR starting from istart position to the next DELIM.
! After EXIT istart points to the next position after the DELIM => usable for sequential read.
!------------------------------------------------------------------------------!
  USE FILINX
!------------------------------------------------------------------------------!
  IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
!------------------------------------------------------------------------------!
  CHARACTER delim
  INTEGER*4 istart
  CHARACTER*(*) Str
  REAL*8 Value
!------------------------------------------------------------------------------!
  INTEGER*4 icomma,ln
!------------------------------------------------------------------------------!
  iresult=0 !--OK
  Value=0.d0
  ln=len_trim(Str)
  icomma = index(Str(istart:), delim) + istart -1
  if (icomma < istart) then
    if (istart<=ln) then !--Read evrth after the delim (non-blank)
      read(Str(istart:ln),*) Value
      istart=ln+1
    else !--Missing value
      call Msg(0,IWR,'INPUT WARNING! - Missing Real Double Precision Value - Zero assumed.')
      iresult=2 !--Missing value
    endif
  elseif (len_trim(Str(istart:icomma-1))==0) then !--If empty or blank string
    istart=icomma+1
    iresult=1 !--Empty or blank string
  else
    read(Str(istart:icomma-1),*) Value
    istart=icomma+1
  endif
!------------------------------------------------------------------------------!
!     End of RDR8 group. 
!------------------------------------------------------------------------------!
  RETURN
END FUNCTION
!==============================================================================!

!==============================================================================!
INTEGER*4 FUNCTION RDStr(delim,istart,str, StrValue,iKeepCase) RESULT(iresult)
!------------------------------------------------------------------------------!
! RDStr: ReaD a String.
! From string STR starting from istart position to the next DELIM.
! After EXIT istart points to the next position after the DELIM => usable for sequential read.
!------------------------------------------------------------------------------!
  USE FILINX
!------------------------------------------------------------------------------!
  IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
!------------------------------------------------------------------------------!
  CHARACTER delim
  INTEGER*4 istart
  CHARACTER*(*) Str, StrValue
  INTEGER*4, OPTIONAL :: iKeepCase !--optional argument: 0 or absent - convert to lower case, otherwise - keep case
!------------------------------------------------------------------------------!
  CHARACTER*255 str1
  INTEGER*4 icomm,ln
!------------------------------------------------------------------------------!
  iresult=0 !--OK
  ln=len_trim(Str)
  icomma = index(Str(istart:), delim) + istart -1
  if (icomma < istart) then
    if (istart<=ln) then !--Read evrth after the delim (non-blank)
      read(Str(istart:ln), '(A)') StrValue
      StrValue=adjustl(StrValue)
      istart=ln+1
    else !--Missing value
      call Msg(0,IWR,'INPUT WARNING! - Missing Character String Value - "null" assumed.')
      StrValue = 'null'
      iresult=2 !--Missing value
    endif
  elseif (len_trim(Str(istart:icomma-1))==0) then !--If empty or blank string
    StrValue = 'null'
    istart = icomma +1
    iresult = 1 !--Empty or blank string
  else
    read(Str(istart:icomma-1), '(A)') StrValue
    StrValue=adjustl(StrValue)
    istart=icomma+1
  endif
 !-----Remove leading and trailing quotes-----!
  ln=len_trim(StrValue)
  if (StrValue(ln:ln)=='"') StrValue(ln:ln)=' '
  if (StrValue(1:1)=='"') StrValue(1:1)=' '
  StrValue=adjustl(StrValue)
 !-----Optionally convert to lower case-----!
  if (PRESENT(iKeepCase)) then
    if (iKeepCase==0) call lcase(StrValue)
  else
    call lcase(StrValue)
  endif
!------------------------------------------------------------------------------!
!     End of RDStr group. 
!------------------------------------------------------------------------------!
  RETURN
END FUNCTION
!==============================================================================!

!==============================================================================!
SUBROUTINE ConvToSI(Value, Units)
!------------------------------------------------------------------------------!
!     ConvToSI: Convert Value from Units To SI base units.
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
!------------------------------------------------------------------------------!
  USE FILINX
!------------------------------------------------------------------------------!
  IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
!------------------------------------------------------------------------------!
  REAL*8 Value
  CHARACTER*80 Units
!------------------------------------------------------------------------------!
!     Parameter Statements.
!------------------------------------------------------------------------------!
  PARAMETER (LUNS=56)
!------------------------------------------------------------------------------!
!  CHARACTER*4 FORM1
  CHARACTER*10 CHS(LUNS),CHD
  REAL*8 CF(LUNS)
  CHARACTER FORM2*6
!------------------------------------------------------------------------------!
!     Data Statements.
!------------------------------------------------------------------------------!
  SAVE CHS,CF,FORM2
!  SAVE FORM1
  DATA CHS /'m','kg','s','j','celsius','pa','w','kgmol','rad','n', &
 'liquid','liq.','gas','solid','sol.','soil', &
 'mm','cm','in','ft','yd','km', &
 'liter','l','gal', &
 'gm','lbm','slug', &
 'min','h','day','wk','fortnight','yr', &
 'btu','cal','hp', &
 'kelvin','fahrenheit','rankine', &
 'k', 'f', 'r', &
 'psi','bar','atm', &
 'degrees','deg', &
 'cp','p','hc','1','mol','lbmol','debyes', &
 '(absolute)'/
 DATA CF  /1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0, &
  1.d0,1.d0,1.d0,1.d0,1.d0,1.d0, &
  1.D-3,1.D-2,2.54D-2,0.3048D+0,0.9144D+0,1.d3, &
  1.D-3,1.D-3,3.7854D-3, &
  1.D-3,4.5359D-1,1.4594D+1, &
  6.d1,3.6D+3,8.64D+4,6.048D+5,1.2096D+6,3.15576D+7, &
  1.0544D+3,4.184D+0,7.457D+2, &
  1.d0,5.555556D-1,5.555556D-1, &
  1.d0,5.555556D-1,5.555556D-1, &
  6894.8D+0,1.d5,1.01325D+5, &
  1.745329251994328D-002,1.745329251994328D-002, &
  1.D-3,1.D-1,1.02286D-7,1.d0,1.D-3,4.5359D-1,1.d0, &
  1.d0/
!  DATA FORM1 /'(I)'/
  DATA FORM2 /'(F .0)'/
!------------------------------------------------------------------------------!
  call lcase(Units)
  if (Units=='null' .or. Units=='none') RETURN
!----------------------------------------------------------------------!
!     Decompose the units into components and convert individual
!     components.
!----------------------------------------------------------------------!
  IS = 1
  IDV = index(Units(1:),'/') -1
  IE  = index(Units(1:),'  ') -1
!----------------------------------------------------------------------!
!     Units without a divisor.
!----------------------------------------------------------------------!
  if (IDV==-1) then
   100 CONTINUE

    ISP = index(Units(IS:),' ') +IS -2
    IB = min(IE,ISP)
    CHD = Units(IS:IB)
    IC = index(CHD(1:),'^')

    if (IC == 0) then
      IP = 1
      PW = 1.d0
    else
      I1 = IC +1
      I2 = IB -IS +1
      I3 = I2 -I1 +1
      write(FORM2(3:3),'(I1)') I3
      read(CHD(I1:I2),FORM2) PW
!      write(FORM1(3:3),'(I1)') I3
!      read(CHD(I1:I2),FORM1) IP
      I2 = IC -1
      CHD = CHD(1:I2)
    endif

    do N = 1,LUNS
      if (CHS(N) == CHD) then
        Value = Value*(CF(N)**PW)
!        Value = Value*(CF(N)**IP)
        GOTO 120
      endif
    enddo

    call Msg(0,IWR,'ERROR! - Unrecognized Input Units: '//trim(Units))
    STOP

   120 CONTINUE
    if (IB<IE) then
      IS = IB +2
      GOTO 100
    endif
!----------------------------------------------------------------------!
!     Units with a divisor.
!----------------------------------------------------------------------!
  else
   !----------------------------------------------------------------------!
   !     Components before the divisor.
   !----------------------------------------------------------------------!
   200 CONTINUE
    ISP = index(Units(IS:),' ') +IS -2
    ICO = index(Units(IS:),':') +IS -2

    if ((ICO>0).and.(ICO>IS)) then
      IB = min(IDV,ISP,ICO)
    else
      IB = min(IDV,ISP)
    endif

    CHD = Units(IS:IB)
    IC = index(CHD(1:),'^')

    if (IC == 0) then
      IP = 1
      PW = 1.d0
    else
      I1 = IC +1
      I2 = IB -IS +1
      I3 = I2 -I1 +1
      write(FORM2(3:3),'(I1)') I3
      read(CHD(I1:I2),FORM2) PW
!      write(FORM1(3:3),'(I1)') I3
!      read(CHD(I1:I2),FORM1) IP
      I2 = IC -1
      CHD = CHD(1:I2)
    endif

    do N = 1,LUNS
      if (CHS(N) == CHD) then
        Value = Value*(CF(N)**PW)
!        Value = Value*(CF(N)**IP)
        GOTO 220
      endif
    enddo

    call Msg(0,IWR,'ERROR! - Unrecognized Input Units: '//trim(Units))
    STOP

   220 CONTINUE
    if (IB<IDV) then
      IS = IB +2
      GOTO 200
    else
      IS = IB +2
      GOTO 300
    endif
   !----------------------------------------------------------------------!
   !     Components after the divisor.
   !----------------------------------------------------------------------!
   300 CONTINUE

    ISP = index(Units(IS:),' ') +IS -2
    IB = min(IE,ISP)
    CHD = Units(IS:IB)
    IC = index(CHD(1:),'^')

    if (IC == 0) then
      IP = 1
      PW = 1.d0
    else
      I1 = IC+1
      I2 = IB-IS+1
      I3 = I2-I1+1
      write(FORM2(3:3),'(I1)') I3
      read(CHD(I1:I2),FORM2) PW
!      write(FORM1(3:3),'(I1)') I3
!      read(CHD(I1:I2),FORM1) IP
      I2 = IC-1
      CHD = CHD(1:I2)
    endif

    do N = 1,LUNS
      if (CHS(N) == CHD) then
        Value = Value/(CF(N)**PW)
!        Value = Value/(CF(N)**IP)
        GOTO 320
      endif
    enddo

    call Msg(0,IWR,'ERROR! - Unrecognized Input Units: '//trim(Units))
    STOP

   320 CONTINUE
    if (IB<IE) then
      IS = IB +2
      GOTO 300
    endif
  endif
!------------------------------------------------------------------------------!
!     End of ConvToSI group. 
!------------------------------------------------------------------------------!
  RETURN
END SUBROUTINE
!==============================================================================!
    
!==============================================================================!
SUBROUTINE Msg(IFU1,IFU2,str,iCon)
!------------------------------------------------------------------------------!
  IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
!------------------------------------------------------------------------------!
  INTEGER*4 IFU1,IFU2
  CHARACTER*(*) str
  INTEGER*4, OPTIONAL :: iCon !--optional argument: 0 or absent - break line, otherwise - continue line
!------------------------------------------------------------------------------!
  CHARACTER fmt*4
!------------------------------------------------------------------------------!
  ln=len_trim(str)
  if (PRESENT(iCon)) then
    if (iCon==0) then
      fmt='(A)'
    else
      fmt='(A$)'
    endif
  else
    fmt='(A)'
  endif
  if (IFU1>0) then
    write(IFU1,fmt) str(1:ln)
  elseif (IFU1==0) then
    write(*,fmt) str(1:ln)
  endif
  if (IFU2/=IFU1) then !--No double echo
    if (IFU2>0)then
      write(IFU2,fmt) str(1:ln)
    elseif (IFU2==0) then
      write(*,fmt) str(1:ln)
    endif
  endif
!------------------------------------------------------------------------------!
!     End of Msg group. 
!------------------------------------------------------------------------------!
  RETURN
    END SUBROUTINE
!==============================================================================!
  
!==============================================================================!
SUBROUTINE AssignGradientVariable(V,Value,GValue,iRange)
!------------------------------------------------------------------------------!
!     AssignGradientVariable: Assign value according to gradient to field 
!                             array V in region using iRange.
!------------------------------------------------------------------------------! 
  USE FLDINX
  USE DRLINK 
  USE POINTS
!------------------------------------------------------------------------------!
  IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
!------------------------------------------------------------------------------!
  REAL*8 Value,GValue,V(*)   
  INTEGER*4 iRange(4)
!------------------------------------------------------------------------------! 
  if( iRange(1) == 0 ) iRange(1) = 1; 
  if( iRange(2) == 0 ) iRange(2) = nLink; 
  Link1 = max(1,iRange(1)); Link2 = min(nLink,iRange(2)); 
  do i = Link1,Link2 
    if( iRange(3) == 0 ) then 
      iP1 = 1;  
    else  
      iP1 = max(1,iRange(3));   
    endif    
    if( iRange(4) == 0 ) then 
      iP2 = IL(i); 
    else  
      iP2 = min(IL(i),iRange(4)); 
    endif 
    x1 = X(nd(i,iP1))
    do j = iP1,iP2
      N = ND(i,j)
      V(N) = Value +GValue*(X(N) -x1)
    enddo
  enddo
!------------------------------------------------------------------------------!
!     End of AssignGradientVariable group. 
!------------------------------------------------------------------------------!
  RETURN
END
!==============================================================================!

  
!==============================================================================!
SUBROUTINE AddVariableRegion(R,V,W,iRange)
!------------------------------------------------------------------------------!
!     PlusVariableRegion: Result (R) is addition of two 
!                         arrays V and W in region using iRange.
!------------------------------------------------------------------------------! 
  USE FLDINX
  USE DRLINK 
  USE POINTS
!------------------------------------------------------------------------------!
  IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
!------------------------------------------------------------------------------!
  REAL*8 W(*),V(*),R(*)   
  INTEGER*4 iRange(4)
!------------------------------------------------------------------------------! 
  if( iRange(1) == 0 ) iRange(1) = 1; 
  if( iRange(2) == 0 ) iRange(2) = nLink; 
  Link1 = max(1,iRange(1)); Link2 = min(nLink,iRange(2)); 
  do i = Link1,Link2 
    if( GRange(3) == 0 ) then 
      iP1 = 1;  
    else  
      iP1 = max(1,GRange(3));   
    endif    
    if( GRange(4) == 0 ) then 
      iP2 = IL(i); 
    else  
      iP2 = min(IL(i),GRange(4)); 
    endif 
    x1 = X(nd(i,iP1))
    do j = iP1,iP2
      N = ND(i,j)
      R(N) = V(N) +W(N)
    enddo
  enddo
!------------------------------------------------------------------------------!
!     End of AddVariableRegion group. 
!------------------------------------------------------------------------------!
  RETURN
END
!==============================================================================!

    

    
