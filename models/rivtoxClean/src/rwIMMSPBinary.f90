!======================================================================C 
! 
SUBROUTINE ReadSolutionHeader_IMMSPBinary(IFU, ID,NFr,NDim,NN,StartT,EndT)
!------------------------------------------------------------------------------C
!
!     PURPOSE.
!     Read solution file header.
!     Read start and end times,
!     and position the file at the end of file header.
!
!------------------------------------------------------------------------------C
IMPLICIT NONE
!------------------------------------------------------------------------------C
INTEGER*4 IFU
CHARACTER*(*) ID
INTEGER*4 NFr,NDim,NN
REAL*8 StartT,EndT

INTENT(IN) IFU
INTENT(OUT) ID,NFr,NDim,NN,StartT,EndT
!------------------------------------------------------------------------------C
CHARACTER iwFTID*12, iwNAME*40
INTEGER*4 iwNFrames,iwNDim,iwNN
REAL*4 rwStartTime,rwEndTime
!------------------------------------------------------------------------------C
 read(IFU) iwFTID,iwNAME,iwNFrames,iwNDim,iwNN,rwStartTime,rwEndTime
 ID=iwNAME
 NFr=iwNFrames
 NDim=iwNDim
 NN=iwNN
 StartT=rwStartTime
 EndT=rwEndTime
!------------------------------------------------------------------------------C
END
!======================================================================C 

!======================================================================C 
! 
SUBROUTINE ReadSolutionHeaderTimes_IMMSPBinary(IFU, ID,NFr,NDim,NN,Times)
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
INTEGER*4 NFr,NDim,NN
REAL*8 Times(*)

INTENT(IN) IFU
INTENT(OUT) ID,NFr,NDim,NN,Times
!------------------------------------------------------------------------------C
CHARACTER iwFTID*12, iwNAME*40
INTEGER*4 iwNFrames,iwNDim,iwNN
REAL*4 rwStartTime,rwEndTime

REAL*4 rwTime

INTEGER*8 IPos
INTEGER*4 iw
REAL*4 rw
!------------------------------------------------------------------------------C
!----------------------------------------------------------------------!
!     Read solution file header.
!----------------------------------------------------------------------!
 read(IFU) iwFTID,iwNAME,iwNFrames,iwNDim,iwNN,rwStartTime,rwEndTime
 ID=iwNAME
 NFr=iwNFrames
 NDim=iwNDim
 NN=iwNN
!----------------------------------------------------------------------!
!     Loop through frames.
!----------------------------------------------------------------------!
 NFr=0
 do while (.true.)
   read(IFU,end=1) rwTime
   !-----Skip frame-----!
   !IPos=FTELLI8(IFU)+1
   IPos=FTELL(IFU)+1
   IPos=IPos+NDim*NN*4
   read(IFU,POS=IPos) !--Position the file after the end of this frame (FTELL=IPos-1)
   !----------------------!
   NFr=NFr+1
   Times(NFr)=rwTime
 enddo
 1 CONTINUE
!----------------------------------------------------------------------!
!     Position the file after the end of file header.
!----------------------------------------------------------------------!
 read(IFU,POS=4*18+1) !--Position the file after the end of file header (FTELL=4*18)
!------------------------------------------------------------------------------C
END
!======================================================================C 

!======================================================================C 
! 
SUBROUTINE ReadSolutionFrame_IMMSPBinary(IFU,NDim,NN, frTime,v1,v2)
!------------------------------------------------------------------------------C
IMPLICIT NONE
!------------------------------------------------------------------------------C
INTEGER*4 IFU
INTEGER*4 NDim,NN
REAL*8 frTime
REAL*8 v1(NN),v2(NN)

INTENT(IN) IFU,NDim,NN
INTENT(OUT) frTime,v1,v2
!------------------------------------------------------------------------------C
REAL*4 rwTime
REAL*4 rwa(2*NN)
INTEGER*4 i
!------------------------------------------------------------------------------C
 read(IFU) rwTime
 frTime=rwTime

 if (NDim==1) then
   read(IFU) rwa(1:NN)
   do i=1,NN
     v1(i)=rwa(i)
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
SUBROUTINE ReadSolutionFrameR4_IMMSPBinary(IFU,NDim,NN, frTime,v1,v2)
!------------------------------------------------------------------------------C
IMPLICIT NONE
!------------------------------------------------------------------------------C
INTEGER*4 IFU
INTEGER*4 NDim,NN
REAL*4 frTime
REAL*4 v1(NN),v2(NN)

INTENT(IN) IFU,NDim,NN
INTENT(OUT) frTime,v1,v2
!------------------------------------------------------------------------------C
REAL*4 rwa(2*NN)
INTEGER*4 i
!------------------------------------------------------------------------------C
 read(IFU) frTime

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
SUBROUTINE SkipSolutionFrame_IMMSPBinary(IFU,NDim,NN)
!------------------------------------------------------------------------------C
!USE IFPORT
IMPLICIT NONE
!------------------------------------------------------------------------------C
INTEGER*4 IFU
INTEGER*4 NDim,NN

INTENT(IN) IFU,NDim,NN
!------------------------------------------------------------------------------C
INTEGER*8 IPos
!------------------------------------------------------------------------------C
 !IPos=FTELLI8(IFU)+1
 IPos=FTELL(IFU)+1
 IPos=IPos+4 !--Skip frame header
 IPos=IPos+NDim*NN*4 !--Skip data values
 read(IFU,POS=IPos) !--Position the file index at the beginning of the next frame
!------------------------------------------------------------------------------C
END
!======================================================================C 

!======================================================================C 
! 
SUBROUTINE WriteSolutionHeader_IMMSPBinary(IFU,ID,NFr,NDim,NN,Times)
!------------------------------------------------------------------------------C
IMPLICIT NONE
!------------------------------------------------------------------------------C
INTEGER*4 IFU
CHARACTER*(*) ID
INTEGER*4 NFr,NDim,NN
REAL*8 Times(*)

INTENT(IN) IFU,ID,NFr,NDim,NN,Times
!------------------------------------------------------------------------------C
CHARACTER iwFTID*12, iwNAME*40
INTEGER*4 iwNFrames,iwNDim,iwNN
REAL*4 rwStartTime,rwEndTime
!------------------------------------------------------------------------------C
 iwFTID='IMMSPBinary '
 iwNAME=ID
 iwNFrames=NFr
 iwNDim=NDim
 iwNN=NN
 rwStartTime=Times(1)
 rwEndTime=Times(NFr)

 write(IFU) iwFTID,iwNAME,iwNFrames,iwNDim,iwNN,rwStartTime,rwEndTime
!------------------------------------------------------------------------------C
END
!======================================================================C 

!======================================================================C 
! 
SUBROUTINE WriteSolutionFrame_IMMSPBinary(IFU,NDim,NN,frTime,v1,v2)
!------------------------------------------------------------------------------C
IMPLICIT NONE
!------------------------------------------------------------------------------C
INTEGER*4 IFU
INTEGER*4 NDim,NN
REAL*8 frTime
REAL*8 v1(NN),v2(NN)

INTENT(IN) IFU,frTime,NDim,NN,v1,v2
!------------------------------------------------------------------------------C
REAL*4 rwTime
REAL*4 rwa(2*NN)
INTEGER*4 i
!------------------------------------------------------------------------------C
 rwTime=frTime
 write(IFU) rwTime

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




!//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
SUBROUTINE wIMMSPBinary(FName, fRewrite,fFirstWrite, ParNameL, frTime, NDim, NN, V1,V2, StartT,EndT,NFrames)
!------------------------------------------------------------------------------C
!
!     wIMMSPBinary:  All inclusive writing routine.
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
REAL*8 StartT,EndT
INTEGER*4 NFrames

INTENT(IN) FName,fRewrite,ParNameL,frTime,NDim, NN, V1,V2
INTENT(INOUT) fFirstWrite, StartT,EndT,NFrames
!------------------------------------------------------------------------------C
INTEGER*4 IFU
CHARACTER*20 FStatus
INTEGER*8 FileHeaderSize,FrameHeaderSize
INTEGER*8 IPos,FrameSize

CHARACTER iwFTID*12,iwNAME*40
INTEGER*4 iwNFrames,iwNDim,iwNN
REAL*4 rwStartTime,rwEndTime,rwTime

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

 !OPEN(UNIT=IFU, FILE=FName, ACCESS='STREAM', FORM='BINARY', STATUS=FStatus, BUFFERED='YES') ! Buffered, Form = 'UNFORMATTED'
 OPEN(UNIT=IFU, FILE=FName, ACCESS='STREAM', FORM='UNFORMATTED', STATUS=FStatus) 
!----------------------------------------------------------------------!
!   Define header sizes.
!----------------------------------------------------------------------!
 FileHeaderSize=(3+10+5)*nbyte
 FrameHeaderSize=1*nbyte
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
     read(IFU) iwFTID,iwNAME,iwNFrames,iwNDim,iwNN,rwStartTime,rwEndTime
     !-----Save-----!
     StartT=rwStartTime
     EndT=rwEndTime
     NFrames=iwNFrames
     !-----Check-----!
     if (iwNDim/=NDim) then
       write(*,'(A)') 'ERROR!'
       write(*,'(A,I0)') 'Dimension of file '//trim(FName)//' = ', iwNDim
       write(*,'(A,I0)') 'Requested dimension = ', NDim
       GOTO 99
     endif
     if (iwNN/=NN) then
       write(*,'(A)') 'ERROR!'
       write(*,'(A,I0)') 'Number of nodes of file '//trim(FName)//' = ', iwNN
       write(*,'(A,I0)') 'Requested number of nodes = ', NN
       GOTO 99
     endif
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
 iwFTID='IMMSPBinary '
 iwNAME=trim(ParNameL)
! iwNFrames=NFrames
 iwNFrames=NFrames+1 !--NFrames is updated outside now
 iwNDim=NDim
 iwNN=NN
 rwStartTime=StartT
 rwEndTime=frTime
 rwTime=frTime

!----------------------------------------------------------------------!
!     Write file header.
!----------------------------------------------------------------------!
 IPos=1 !--Starting position {1-Intel, 0-Portland}
 write(IFU,POS=IPos) iwFTID,iwNAME,iwNFrames,iwNDim,iwNN,rwStartTime,rwEndTime

!----------------------------------------------------------------------!
!     Write last frame header.
!----------------------------------------------------------------------!
 IPos=IPos+FileHeaderSize
! IPos=IPos+(NFrames-1)*(FrameHeaderSize+FrameSize) !--Skip all previous frames
 IPos=IPos+NFrames*(FrameHeaderSize+FrameSize) !--Skip all previous frames (NFrames is updated outside now)
 write(IFU,POS=IPos) rwTime

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


 99 CONTINUE
!----------------------------------------------------------------------!
!     Close file.
!----------------------------------------------------------------------!
 CLOSE(IFU)
!------------------------------------------------------------------------------C
RETURN
END
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

