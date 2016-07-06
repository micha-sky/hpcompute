!======================================================================C 
! 
INTEGER*4 FUNCTION ReadInteger(IFU,iSize,x) RESULT(result)
!------------------------------------------------------------------------------C
IMPLICIT NONE
!------------------------------------------------------------------------------C
INTEGER*4 IFU,iSize,x
!------------------------------------------------------------------------------C
INTEGER*1 i1wr
INTEGER*2 i2wr
INTEGER*4 i4wr
INTEGER*8 i8wr
!------------------------------------------------------------------------------C
 select case (iSize)
 case(1)
   read(IFU,err=1,end=2) i1wr
   x=i1wr
 case(2)
   read(IFU,err=1,end=2) i2wr
   x=i2wr
 case(4)
   read(IFU,err=1,end=2) i4wr
   x=i4wr
 case(8)
   read(IFU,err=1,end=2) i8wr
   x=i8wr
 end select
 result=0 !--OK
 RETURN
 1 CONTINUE !--Error reading
 result=1
 RETURN
 2 CONTINUE !--End-of-file
 result=2
 RETURN
!------------------------------------------------------------------------------C
END
!======================================================================C 

!======================================================================C 
! 
INTEGER*4 FUNCTION ReadIntegerPos(IFU,iSize,IPos,x) RESULT(result)
!------------------------------------------------------------------------------C
IMPLICIT NONE
!------------------------------------------------------------------------------C
INTEGER*4 IFU,iSize,x
INTEGER*8 IPos
!------------------------------------------------------------------------------C
INTEGER*1 i1wr
INTEGER*2 i2wr
INTEGER*4 i4wr
INTEGER*8 i8wr
!------------------------------------------------------------------------------C
 select case (iSize)
 case(1)
   read(IFU,POS=IPos,err=1,end=2) i1wr
   x=i1wr
 case(2)
   read(IFU,POS=IPos,err=1,end=2) i2wr
   x=i2wr
 case(4)
   read(IFU,POS=IPos,err=1,end=2) i4wr
   x=i4wr
 case(8)
   read(IFU,POS=IPos,err=1,end=2) i8wr
   x=i8wr
 end select
 result=0 !--OK
 RETURN
 1 CONTINUE !--Error reading
 result=1
 RETURN
 2 CONTINUE !--End-of-file
 result=2
 RETURN
!------------------------------------------------------------------------------C
END
!======================================================================C 

!======================================================================C 
! 
INTEGER*4 FUNCTION ReadReal(IFU,iSize,x) RESULT(result)
!------------------------------------------------------------------------------C
IMPLICIT NONE
!------------------------------------------------------------------------------C
INTEGER*4 IFU,iSize
REAL*8 x
!------------------------------------------------------------------------------C
REAL*4 r4wr
REAL*8 r8wr
!------------------------------------------------------------------------------C
 select case (iSize)
 case(4)
   read(IFU,err=1,end=2) r4wr
   x=r4wr
 case(8)
   read(IFU,err=1,end=2) r8wr
   x=r8wr
 end select
 result=0 !--OK
 RETURN
 1 CONTINUE !--Error reading
 result=1
 RETURN
 2 CONTINUE !--End-of-file
 result=2
 RETURN
!------------------------------------------------------------------------------C
END
!======================================================================C 

!======================================================================C 
! 
INTEGER*4 FUNCTION ReadRealPos(IFU,iSize,IPos,x) RESULT(result)
!------------------------------------------------------------------------------C
IMPLICIT NONE
!------------------------------------------------------------------------------C
INTEGER*4 IFU,iSize
INTEGER*8 IPos
REAL*8 x
!------------------------------------------------------------------------------C
REAL*4 r4wr
REAL*8 r8wr
!------------------------------------------------------------------------------C
 select case (iSize)
 case(4)
   read(IFU,POS=IPos,err=1,end=2) r4wr
   x=r4wr
 case(8)
   read(IFU,POS=IPos,err=1,end=2) r8wr
   x=r8wr
 end select
 result=0 !--OK
 RETURN
 1 CONTINUE !--Error reading
 result=1
 RETURN
 2 CONTINUE !--End-of-file
 result=2
 RETURN
!------------------------------------------------------------------------------C
END
!======================================================================C 

!======================================================================C 
! 
SUBROUTINE ReadIntegerA(IFU,iSize,n,x)
!------------------------------------------------------------------------------C
IMPLICIT NONE
!------------------------------------------------------------------------------C
INTEGER*4 IFU,iSize,n
INTEGER*4 x(n)
!------------------------------------------------------------------------------C
INTEGER*4 i
INTEGER*1 iwa1(n)
INTEGER*2 iwa2(n)
INTEGER*8 iwa8(n)
!------------------------------------------------------------------------------C
 select case (iSize)
 case(1)
   read(IFU) iwa1(1:n)
   x(1:n)=iwa1(1:n)
 case(2)
   read(IFU) iwa2(1:n)
   x(1:n)=iwa2(1:n)
 case(4)
   read(IFU) x(1:n)
 case(8)
   read(IFU) iwa8(1:n)
   x(1:n)=iwa8(1:n)
 end select
!------------------------------------------------------------------------------C
END
!======================================================================C 

!======================================================================C 
! 
SUBROUTINE ReadRealA(IFU,iSize,n,x)
!------------------------------------------------------------------------------C
IMPLICIT NONE
!------------------------------------------------------------------------------C
INTEGER*4 IFU,iSize,n
REAL*8 x(n)
!------------------------------------------------------------------------------C
INTEGER*4 i
REAL*4 rwa4(n)
!------------------------------------------------------------------------------C
 select case (iSize)
 case(4)
   read(IFU) rwa4(1:n)
   x(1:n)=rwa4(1:n)
 case(8)
   read(IFU) x(1:n)
 end select
!------------------------------------------------------------------------------C
END
!======================================================================C 

!======================================================================C 
! 
SUBROUTINE ReadSolutionHeader_SMSGenericBinary(IFU, ID,IObT,NC,NFr,NDim,NN,ITU,StartT,EndT,iPar)
!------------------------------------------------------------------------------C
!
!     PURPOSE.
!     Read solution file header.
!     Calculate the actual number of data sets in the solution file,
!     read start and end times,
!     and position the file at the end of file header.
!
!------------------------------------------------------------------------------C
!USE IFPORT
IMPLICIT NONE
!------------------------------------------------------------------------------C
INTEGER*4 IFU
CHARACTER*(*) ID
INTEGER*4 IObT,NC,NFr,NDim,NN,ITU,iPar(10)
REAL*8 StartT,EndT

INTENT(IN) IFU
INTENT(INOUT) ITU
INTENT(OUT) ID,IObT,NC,NFr,NDim,NN,StartT,EndT,iPar
!------------------------------------------------------------------------------C
INTEGER*4 IDCard, IDDIM

INTEGER*4 iwVersion,iwObjectType,iwSFLT,iwSFLG,iwVECTYPE
INTEGER*4 iwOBJID,iwNUMDATA,iwNUMCELLS,iwISTAT,iwTimeUnits
REAL*8 rwDateTimeJulian
CHARACTER*40 iwNAME

REAL*8 rwTIME,rw

INTEGER*8 IPos
INTEGER*4 iw,k,i,iresult
!------------------------------------------------------------------------------C
INTEGER*4 ReadInteger,ReadIntegerPos,ReadReal,ReadRealPos
!------------------------------------------------------------------------------C
 NFr=0
!----------------------------------------------------------------------!
!     Loop through cards.
!----------------------------------------------------------------------!
 k=0  !--header size in 4 bytes
 do while (.true.)

   read(IFU,err=1,end=1) IDCard
   k=k+1

   select case (IDCard)
   case(3000)
   case(100)
     read(IFU) iwObjectType
     k=k+1
   case(110)
     read(IFU) iwSFLT
     k=k+1
   case(120)
     read(IFU) iwSFLG
     k=k+1
   case(130)
     NDim=1
   case(140)
     NDim=2
   case(150)
     read(IFU) iwVECTYPE
     k=k+1
   case(160)
     read(IFU) iwOBJID
     k=k+1
   case(170)
     read(IFU) iwNUMDATA
     k=k+1
   case(180)
     read(IFU) iwNUMCELLS
     k=k+1
   case(190)
     read(IFU) iwNAME
     k=k+10
   case(240)
     read(IFU) rwDateTimeJulian
     k=k+2
   case(250)
     read(IFU) iwTimeUnits
     k=k+1
     ITU=iwTimeUnits
   case(200)
     k=k-1
     iresult = ReadInteger(IFU,iwSFLG,iwISTAT)
     if (iresult/=0) GOTO 1
     iresult = ReadReal(IFU,iwSFLT,rwTIME)
     if (iresult/=0) GOTO 1
     !-----Skip frame-----!
     !IPos=FTELLI8(IFU)+1 !--Add 1 because FTELL=IPos-1
     IPos=FTELL(IFU) + 1
     !replace by ftell
     if (iwISTAT>0) IPos=IPos+iwNUMCELLS*iwSFLG !--Skip status flags !! if (iwISTAT>0) ???
     IPos=IPos+(NDim*iwNUMDATA-1)*iwSFLT !--Skip data values except the last
     iresult = ReadRealPos(IFU,iwSFLT,IPos,rw) !--Try to read the last data value of the frame, and position the file index at the beginning of the next frame
     if (iresult/=0) GOTO 1
     !----------------------!
     NFr=NFr+1
     if (NFr==1) StartT=rwTIME
     EndT=rwTIME
   case(210)
     k=k-1
     GOTO 2
   case default
     k=k-1
     write(*,'(A,I0,A)') 'WARNING! - Unidentified card: ',IDCard,'.'
     GOTO 2
   end select

 enddo

 1 CONTINUE
 write(*,'(A)') 'WARNING! - File is corrupted.'
 2 CONTINUE
 ID= iwNAME
 IObT=iwObjectType
 NC=iwNUMCELLS
 NN=iwNUMDATA
 iPar(1)=iwSFLT
 iPar(2)=iwSFLG
 iPar(3)=iwNUMCELLS
 !----------

 read(IFU,POS=4*k+1) !--Position the file after the end of file header (FTELL=4*k)
!------------------------------------------------------------------------------C
END
!======================================================================C 

!======================================================================C 
! 
SUBROUTINE ReadSolutionHeaderTimes_SMSGenericBinary(IFU, ID,IObT,NC,NFr,NDim,NN,ITU,Times,iPar)
!------------------------------------------------------------------------------C
!
!     PURPOSE.
!     Read solution file header.
!     Calculate the actual number of data sets in the solution file,
!     read all frame times,
!     and position the file at the end of file header.
!
!------------------------------------------------------------------------------C
!USE IFPORT
IMPLICIT NONE
!------------------------------------------------------------------------------C
INTEGER*4 IFU
CHARACTER*(*) ID
INTEGER*4 IObT,NC,NFr,NDim,NN,ITU,iPar(10)
REAL*8 Times(*)

INTENT(IN) IFU
INTENT(INOUT) ITU
INTENT(OUT) ID,IObT,NC,NFr,NDim,NN,Times,iPar
!------------------------------------------------------------------------------C
INTEGER*4 IDCard, IDDIM

INTEGER*4 iwVersion,iwObjectType,iwSFLT,iwSFLG,iwVECTYPE
INTEGER*4 iwOBJID,iwNUMDATA,iwNUMCELLS,iwISTAT,iwTimeUnits
REAL*8 rwDateTimeJulian
CHARACTER*40 iwNAME

REAL*8 rwTIME,rw

INTEGER*8 IPos
INTEGER*4 iw,k,i,iresult
!------------------------------------------------------------------------------C
INTEGER*4 ReadInteger,ReadIntegerPos,ReadReal,ReadRealPos
!------------------------------------------------------------------------------C
 NFr=0
!----------------------------------------------------------------------!
!     Loop through cards.
!----------------------------------------------------------------------!
 k=0  !--header size in 4 bytes
 do while (.true.)

   read(IFU,err=1,end=1) IDCard
   k=k+1

   select case (IDCard)
   case(3000)
   case(100)
     read(IFU) iwObjectType
     k=k+1
   case(110)
     read(IFU) iwSFLT
     k=k+1
   case(120)
     read(IFU) iwSFLG
     k=k+1
   case(130)
     NDim=1
   case(140)
     NDim=2
   case(150)
     read(IFU) iwVECTYPE
     k=k+1
   case(160)
     read(IFU) iwOBJID
     k=k+1
   case(170)
     read(IFU) iwNUMDATA
     k=k+1
   case(180)
     read(IFU) iwNUMCELLS
     k=k+1
   case(190)
     read(IFU) iwNAME
     k=k+10
   case(240)
     read(IFU) rwDateTimeJulian
     k=k+2
   case(250)
     read(IFU) iwTimeUnits
     k=k+1
     ITU=iwTimeUnits
   case(200)
     k=k-1
     iresult = ReadInteger(IFU,iwSFLG,iwISTAT)
     if (iresult/=0) GOTO 1
     iresult = ReadReal(IFU,iwSFLT,rwTIME)
     if (iresult/=0) GOTO 1
     !-----Skip frame-----!
     IPos=FTELL(IFU)+1 !--Add 1 because FTELL=IPos-1
     if (iwISTAT>0) IPos=IPos+iwNUMCELLS*iwSFLG !--Skip status flags
     IPos=IPos+(NDim*iwNUMDATA-1)*iwSFLT !--Skip data values except the last
     iresult = ReadRealPos(IFU,iwSFLT,IPos,rw) !--Try to read the last data value of the frame, and position the file index at the beginning of the next frame
     if (iresult/=0) GOTO 1
     !----------------------!
     NFr=NFr+1
     Times(NFr)=rwTIME
   case(210)
     k=k-1
     GOTO 2
   case default
     k=k-1
     write(*,'(A,I0,A)') 'WARNING! - Unidentified card: ',IDCard,'.'
     GOTO 2
   end select

 enddo

 1 CONTINUE
 write(*,'(A)') 'WARNING! - File is corrupted.'
 2 CONTINUE
 ID= iwNAME
 IObT=iwObjectType
 NC=iwNUMCELLS
 NN=iwNUMDATA
 iPar(1)=iwSFLT
 iPar(2)=iwSFLG
 iPar(3)=iwNUMCELLS
 !----------

 read(IFU,POS=4*k+1) !--Position the file after the end of file header (FTELL=4*k)
!------------------------------------------------------------------------------C
END
!======================================================================C 

!======================================================================C 
! 
SUBROUTINE ReadSolutionFrame_SMSGenericBinary(IFU,NDim,NN,iPar, frTime,v1,v2)
!------------------------------------------------------------------------------C
!USE IFPORT
IMPLICIT NONE
!------------------------------------------------------------------------------C
INTEGER*4 IFU
INTEGER*4 NDim,NN,iPar(10)
REAL*8 frTime
REAL*8 v1(NN),v2(NN)

INTENT(IN) IFU,NDim,NN,iPar
INTENT(OUT) frTime,v1,v2
!------------------------------------------------------------------------------C
INTEGER*4 IDCard, IDDIM
INTEGER*4 iwVersion,iwObjectType,iwSFLT,iwSFLG,iwVECTYPE
INTEGER*4 iwOBJID,iwNUMDATA,iwNUMCELLS,iwISTAT

INTEGER*8 IPos
INTEGER*4 iw,i,iresult
REAL*8 rwa(2*NN)
!------------------------------------------------------------------------------C
INTEGER*4 ReadInteger,ReadIntegerPos,ReadReal,ReadRealPos
!------------------------------------------------------------------------------C
 iwSFLT=iPar(1)
 iwSFLG=iPar(2)
 iwNUMCELLS=iPar(3)

 read(IFU) IDCard
 iresult = ReadInteger(IFU,iwSFLG,iwISTAT)
 iresult = ReadReal(IFU,iwSFLT,frTime)

 if (iwISTAT>0) then !--Skip status flags
   IPos=FTELL(IFU)+1
   IPos=IPos+iwNUMCELLS*iwSFLG
   read(IFU,POS=IPos) !--Position the file after the end of this section (FTELL=IPos-1)
 endif

 if (NDim==1) then
   call ReadRealA(IFU,iwSFLT,NN,v1)
   do i=1,NN
     if (v1(i)==-99999) v1(i)=0.d0 !--Set water surface at SMS dry node to zero
   enddo
 else
   call ReadRealA(IFU,iwSFLT,2*NN,rwa)
   do i=1,NN
     v1(i)=rwa(2*i-1)
     v2(i)=rwa(2*i)
   enddo
 endif
!------------------------------------------------------------------------------C
END
!======================================================================C 

!======================================================================C 
! 
SUBROUTINE ReadSolutionFrameR4_SMSGenericBinary(IFU,NDim,NN,iPar, frTime,v1,v2)
!------------------------------------------------------------------------------C
!USE IFPORT
IMPLICIT NONE
!------------------------------------------------------------------------------C
INTEGER*4 IFU
INTEGER*4 NDim,NN,iPar(10)
REAL*4 frTime
REAL*4 v1(NN),v2(NN)

INTENT(IN) IFU,NDim,NN,iPar
INTENT(OUT) frTime,v1,v2
!------------------------------------------------------------------------------C
INTEGER*4 IDCard, IDDIM
INTEGER*4 iwVersion,iwObjectType,iwSFLT,iwSFLG,iwVECTYPE
INTEGER*4 iwOBJID,iwNUMDATA,iwNUMCELLS,iwISTAT

INTEGER*8 IPos
INTEGER*4 i,iw,iresult
REAL*4 rwa(2*NN)
!------------------------------------------------------------------------------C
INTEGER*4 ReadInteger,ReadIntegerPos,ReadReal,ReadRealPos
!------------------------------------------------------------------------------C
 iwSFLT=iPar(1)
 iwSFLG=iPar(2)
 iwNUMCELLS=iPar(3)

 read(IFU) IDCard
 iresult = ReadInteger(IFU,iwSFLG,iwISTAT)
 read(IFU) frTime

 if (iwISTAT>0) then !--Skip status flags
   IPos=FTELL(IFU)+1
   IPos=IPos+iwNUMCELLS*iwSFLG
   read(IFU,POS=IPos) !--Position the file after the end of this section (FTELL=IPos-1)
 endif

 if (NDim==1) then
   read(IFU) v1(1:NN)
   do i=1,NN
     if (v1(i)==-99999) v1(i)=0.0 !--Set water surface at SMS dry node to zero
   enddo
 else
   read(IFU) rwa(1:2*NN)
   do i=1,NN
     v1(i)=rwa(2*i-1)
     v2(i)=rwa(2*i)
   enddo
 endif
!------------------------------------------------------------------------------C
END
!======================================================================C 


!======================================================================C 
! 
SUBROUTINE SkipSolutionFrame_SMSGenericBinary(IFU,NDim,NN,iPar)
!------------------------------------------------------------------------------C
!USE IFPORT
IMPLICIT NONE
!------------------------------------------------------------------------------C
INTEGER*4 IFU
INTEGER*4 NDim,NN,iPar(10)

INTENT(IN) IFU,NDim,NN,iPar
!------------------------------------------------------------------------------C
INTEGER*4 IDCard, IDDIM
INTEGER*4 iwVersion,iwObjectType,iwSFLT,iwSFLG,iwVECTYPE
INTEGER*4 iwOBJID,iwNUMDATA,iwNUMCELLS,iwISTAT

INTEGER*8 IPos
INTEGER*4 iresult
!------------------------------------------------------------------------------C
INTEGER*4 ReadInteger,ReadIntegerPos,ReadReal,ReadRealPos
!------------------------------------------------------------------------------C
 iwSFLT=iPar(1)
 iwSFLG=iPar(2)
 iwNUMCELLS=iPar(3)

 IPos=FTELL(IFU)+1
 iresult = ReadIntegerPos(IFU,iwSFLG,IPos+4,iwISTAT)
 IPos=IPos+4+iwSFLG+iwSFLT !--Skip frame header
 if (iwISTAT>0) IPos=IPos+iwNUMCELLS*iwSFLG !--Skip status flags
 IPos=IPos+NDim*NN*iwSFLT !--Skip data values
 read(IFU,POS=IPos) !--Position the file index at the beginning of the next frame
!------------------------------------------------------------------------------C
END
!======================================================================C 

!======================================================================C 
! 
SUBROUTINE WriteSolutionHeader_SMSGenericBinary(IFU,ID,IObT,NC,NFr,NDim,NN,ITU,Times)
!------------------------------------------------------------------------------C
IMPLICIT NONE
!------------------------------------------------------------------------------C
INTEGER*4 IFU
CHARACTER*(*) ID
INTEGER*4 IObT,NC,NFr,NDim,NN,ITU
REAL*8 Times(*)

INTENT(IN) IFU,ID,IObT,NC,NFr,NDim,NN,ITU,Times
!------------------------------------------------------------------------------C
INTEGER*4 ID100,ID110,ID120,ID130,ID140,ID150,ID160,ID170,ID180
INTEGER*4 ID190,ID200,ID210,IDDim,ID240,ID250

INTEGER*4 iwVersion,iwObjectType,iwSFLT,iwSFLG,iwVECTYPE
INTEGER*4 iwOBJID,iwNUMDATA,iwNUMCELLS,iwTimeUnits
REAL*8 rwDateTimeJulian
CHARACTER*40 iwNAME
!------------------------------------------------------------------------------C
 iwVersion=3000
 iwObjectType=IObT !3-mesh, 5-scatter points
 iwSFLT=4
 iwSFLG=4
 iwOBJID=0
 iwNUMDATA=NN
 if (IObT==3) then !--mesh
   iwNUMCELLS=NC
 elseif (IObT==5) then !--scatter set
   iwNUMCELLS=iwNUMDATA
 endif
 iwNAME=ID
 iwTimeUnits=ITU !--0-hours, 1-minutes, 2-seconds, 4-days

 ID100=100
 ID110=110
 ID120=120
 ID130=130
 ID140=140
 ID150=150
 ID160=160
 ID170=170
 ID180=180
 ID190=190
 ID200=200
 ID210=210
 ID240=240
 ID250=250

 write(IFU) iwVersion, ID100,iwObjectType, ID110,iwSFLT, ID120,iwSFLG

 IDDim=120+NDim*10 !-- 130-scalar, 140-vector --
 write(IFU) IDDim

 write(IFU) ID160,iwOBJID !--necessary for scatter points

 write(IFU) ID170,iwNUMDATA, ID180,iwNUMCELLS, ID190,iwNAME, ID250,iwTimeUnits
!------------------------------------------------------------------------------C
END
!======================================================================C 

!======================================================================C 
! 
SUBROUTINE WriteSolutionFrame_SMSGenericBinary(IFU,NDim,NN,frTime,v1,v2)
!------------------------------------------------------------------------------C
IMPLICIT NONE
!------------------------------------------------------------------------------C
INTEGER*4 IFU
INTEGER*4 NDim,NN
REAL*8 frTime
REAL*8 v1(NN),v2(NN)

INTENT(IN) IFU,frTime,NDim,NN,v1,v2
!------------------------------------------------------------------------------C
INTEGER*4 ID200
INTEGER*4 iwISTAT
REAL*4 rwTime
REAL*4 rwa(2*NN)
INTEGER*4 i
!------------------------------------------------------------------------------C
 iwISTAT=0
 rwTIME=frTime
 ID200=200

 write(IFU) ID200,iwISTAT,rwTIME

 if (NDim==1) then
   rwa(1:NN)=v1(1:NN)
   write(IFU) rwa(1:NN)
 else
   do i=1,NN
     rwa(2*i-1)=v1(i)
     rwa(2*i)=v2(i)
   enddo
   write(IFU) rwa(1:2*NN)
 endif
!------------------------------------------------------------------------------C
END
!======================================================================C 

!======================================================================C 
! 
SUBROUTINE wSMSGenericBinary(FName, fRewrite,fFirstWrite, ParNameL, frTime, NDim, NN, V1,V2, IObT,NC, ITU, StartT,EndT,NFrames)
!------------------------------------------------------------------------------C
!
!     wSMSGenericBinary:  All inclusive writing routine.
!     Can create or append.
!     If appending then gets the number of frames and start,end time,
!     and checks against model parameters (NDim,NN,etc.).
!     And never rewrites already written frames.
!     This routine needs to be called once for every frame,
!     with callee updating StartT,EndT,NFrames and fFirstWrite.
!     Then it always updates file header with correct info.
!
!------------------------------------------------------------------------------C
IMPLICIT NONE
!------------------------------------------------------------------------------C
CHARACTER*(*) FName
LOGICAL*4 fRewrite,fFirstWrite
CHARACTER*(*) ParNameL
REAL*8 frTime
INTEGER*4 NDim, NN
REAL*8 V1(NN),V2(NN)
INTEGER*4 IObT,NC, ITU
REAL*8 StartT,EndT
INTEGER*4 NFrames

INTENT(IN) FName,fRewrite,ParNameL,frTime,NDim, NN, V1,V2, IObT,NC, ITU
INTENT(INOUT) fFirstWrite, StartT,EndT,NFrames
!------------------------------------------------------------------------------C
INTEGER*4 IFU
CHARACTER*20 FStatus
INTEGER*8 FileHeaderSize,FrameHeaderSize
INTEGER*8 IPos,FrameSize

INTEGER*4 IDCard,IDDim,iFDim
INTEGER*4 ID100,ID110,ID120,ID130,ID140,ID150,ID160,ID170,ID180
INTEGER*4 ID190,ID200,ID210,ID240,ID250

INTEGER*4 iwCardID
INTEGER*4 iwVersion,iwObjectType,iwSFLT,iwSFLG,iwVECTYPE
INTEGER*4 iwOBJID,iwNUMDATA,iwNUMCELLS,iwISTAT,iwTimeUnits
REAL*8 rwDateTimeJulian
REAL*4 rwTIME
CHARACTER*40 iwNAME

REAL*4 rwa(2*NN)

LOGICAL*4 fFExist,fFlag

INTEGER*4, PARAMETER :: nbyte=4

INTEGER*4 i
!------------------------------------------------------------------------------C
!----------------------------------------------------------------------!
!   Create or open.
!----------------------------------------------------------------------!
 fFExist=.false.
 INQUIRE(FILE=FName, EXIST = fFExist)

 fFlag=(fRewrite .and. fFirstWrite) .or. (.not. fFExist)

 if (fFlag) then
   FStatus='REPLACE'
   StartT=frTime
   EndT=-1.d30
   NFrames=0
 else
   FStatus='UNKNOWN'
 endif

 IFU=1001

 OPEN(UNIT=IFU, FILE=FName, ACCESS='STREAM', FORM='UNFORMATTED', STATUS=FStatus) !rewrite to gfortran standart - remove buffered?

!----------------------------------------------------------------------!
!   Define header sizes.
!----------------------------------------------------------------------!
 FileHeaderSize=(11+18)*nbyte
 FrameHeaderSize=3*nbyte
 FrameSize=NN*NDim*nbyte

!----------------------------------------------------------------------!
!     If appending, get number of frames and start,end time.
!     Check if frame is already written.
!----------------------------------------------------------------------!
 if (.not. fFlag) then
   !----------------------------------------------------------------------!
   !     For the 1st write get number of frames and start,end time.
   !     Check against model dimensions.
   !----------------------------------------------------------------------!
   if (fFirstWrite) then
     read(IFU) iwVersion, ID100,iwObjectType, &
         ID110,iwSFLT, ID120,iwSFLG, IDDim, ID150,iwVECTYPE, &
         ID160,iwOBJID, ID170,iwNUMDATA, ID180,iwNUMCELLS, ID190,iwNAME, ID250,iwTimeUnits
     !-----Check-----!
     iFDim=(IDDim-120)/10  !--find file dimension--130-scalar,140-vector
     if (iFDim/=NDim) then
       write(*,'(A)') 'ERROR!'
       write(*,'(A,I0)') 'Dimension of file '//trim(FName)//' = ', iFDim
       write(*,'(A,I0)') 'Requested dimension = ', NDim
       GOTO 99
     endif
     if (iwNUMDATA/=NN) then
       write(*,'(A)') 'ERROR!'
       write(*,'(A,I0)') 'Number of nodes of file '//trim(FName)//' = ', iwNUMDATA
       write(*,'(A,I0)') 'Requested number of nodes = ', NN
       GOTO 99
     endif
     !-----Calculate number of frames and start,end times-----!
     StartT=1.d+20
     NFrames=0
     IPos=1 !--Starting position {1-Intel, 0-Portland}
     IPos=IPos+FileHeaderSize
     do
       read(IFU,POS=IPos,err=10,end=10) IDCard, iwISTAT, rwTIME
       if (IDCard/=200) GOTO 10 !--Finish if not new frame
       IPos=IPos+(FrameHeaderSize+FrameSize)
       NFrames=NFrames+1
       StartT=min(StartT,rwTIME)
       EndT=rwTIME
     enddo
     10 CONTINUE
   endif
   !----------------------------------------------------------------------!
   !     Always skip frames already written.
   !----------------------------------------------------------------------!
   if (frTime<=EndT) GOTO 99
 endif

!----------------------------------------------------------------------!
!     Update variables. 
!     (22.03.2012) - Now it's done outside.
!----------------------------------------------------------------------!
! EndT=frTime
! NFrames=NFrames+1

!----------------------------------------------------------------------!
!     Prepare for writing headers.
!----------------------------------------------------------------------!
 iwVersion=3000
 iwObjectType=IObT !3-mesh, 5-scatter points
 iwSFLT=nbyte
 iwSFLG=nbyte
 iwVECTYPE=0
 iwOBJID=0
 iwNUMDATA=NN
 if (IObT==3) then !--mesh
   iwNUMCELLS=NC
 elseif (IObT==5) then !--scatter set
   iwNUMCELLS=iwNUMDATA
 endif
 iwNAME=trim(ParNameL)
 iwISTAT=0
! rwDateTimeJulian=(2000+4713)*365.25 !--01.01.2000 12:00:00
 iwTimeUnits=ITU !--0-hours, 1-minutes, 2-seconds, 4-days

 rwTIME=frTime

 ID100=100
 ID110=110
 ID120=120
 IDDim=120+NDim*10 !-- 130-scalar, 140-vector --
 ID150=150
 ID160=160
 ID170=170
 ID180=180
 ID190=190
 ID200=200
 ID210=210
 ID240=240
 ID250=250

!----------------------------------------------------------------------!
!     Write file header.
!----------------------------------------------------------------------!
 IPos=1 !--Starting position {1-Intel, 0-Portland}
 write(IFU,POS=IPos) iwVersion, ID100,iwObjectType, &
     ID110,iwSFLT, ID120,iwSFLG, IDDim, ID150,iwVECTYPE, &
     ID160,iwOBJID, ID170,iwNUMDATA, ID180,iwNUMCELLS, ID190,iwNAME, ID250,iwTimeUnits

!----------------------------------------------------------------------!
!     Write last frame header.
!----------------------------------------------------------------------!
 IPos=IPos+FileHeaderSize
! IPos=IPos+(NFrames-1)*(FrameHeaderSize+FrameSize) !--Skip all previous frames
 IPos=IPos+NFrames*(FrameHeaderSize+FrameSize) !--Skip all previous frames (NFrames is updated outside now)
 write(IFU,POS=IPos) ID200, iwISTAT, rwTIME

!----------------------------------------------------------------------!
!     Write frame.
!----------------------------------------------------------------------!
 IPos=IPos+FrameHeaderSize
 if (NDim==1) then
   do i=1,NN
     rwa(i)=V1(i)
   enddo
   write(IFU,POS=IPos) rwa(1:NN)
 else
   do i=1,NN
     rwa(2*i-1)=V1(i)
     rwa(2*i)=V2(i)
   enddo
   write(IFU,POS=IPos) rwa(1:2*NN)
 endif

!----------------------------------------------------------------------!
!     Write end of file.
!----------------------------------------------------------------------!
 IPos=IPos+FrameSize
 write(IFU,POS=IPos) ID210


 99 CONTINUE
!----------------------------------------------------------------------!
!     Close file.
!----------------------------------------------------------------------!
 CLOSE(IFU)
!------------------------------------------------------------------------------C
 RETURN
 END
!======================================================================C 
