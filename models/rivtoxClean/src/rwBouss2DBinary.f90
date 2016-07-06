!======================================================================C 
! 
SUBROUTINE wBouss2DBinary(FName, fRewrite,fFirstWrite, ParNameL, frTime, NDim, NX,NY, DelX, DelY, Grid_Orientation, Xmin, Ymin, & 
    DelT, V1,V2, StartT,EndT,NFrames)
!------------------------------------------------------------------------------C
IMPLICIT NONE
!------------------------------------------------------------------------------C
CHARACTER*(*) FName
LOGICAL*4 fRewrite,fFirstWrite
CHARACTER*(*) ParNameL
REAL*8 frTime
INTEGER*4 NDim, NX,NY
REAL*8 DelX,DelY, Grid_Orientation, Xmin,Ymin, DelT
REAL*8 V1(NX*NY),V2(NX*NY)
REAL*8 StartT,EndT
INTEGER*4 NFrames

INTENT(IN) FName,fRewrite,ParNameL,frTime,NDim, NX,NY, DelX, DelY, Grid_Orientation, Xmin, Ymin, DelT, V1,V2
INTENT(INOUT) fFirstWrite, StartT,EndT,NFrames
!------------------------------------------------------------------------------C
INTEGER*4 IFU
CHARACTER*20 FStatus
INTEGER*8 FileHeaderSize,FrameHeaderSize
INTEGER*8 IPos,FrameSize

INTEGER*4 iwNFrames,iwNX,iwNY,iwIType
REAL*4 rwDelX,rwDelY, rwGrid_Orientation, rwXmin,rwYmin, rwStartTime,rwEndTime, rwDelT, rwTime

REAL*4 rwa(NX*NY)

LOGICAL*4 fFExist,fFlag

INTEGER*4, PARAMETER :: nbyte=4
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
   NFrames=0
 else
   FStatus='UNKNOWN'
 endif

 IFU=1001

 OPEN(UNIT=IFU, FILE=FName, ACCESS='STREAM', FORM='UNFORMATTED', STATUS=FStatus) !OPENMARK

!----------------------------------------------------------------------!
!   Define header sizes.
!----------------------------------------------------------------------!
 FileHeaderSize=11*nbyte
 FrameHeaderSize=3*nbyte
 FrameSize=NX*NY*nbyte

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
     read(IFU) iwNX,iwNY, rwDelX, rwDelY, rwGrid_Orientation, iwIType, rwXmin, rwYmin, rwStartTime, rwEndTime, rwDelT
     !-----Save-----!
     StartT=rwStartTime
     EndT=rwEndTime
     !-----Check-----!
     if (iwIType/=NDim) then
       write(*,'(A)') 'ERROR!'
       write(*,'(A,I0)') 'Dimension of file '//trim(FName)//' = ', iwIType
       write(*,'(A,I0)') 'Requested dimension = ', NDim
       GOTO 99
     endif
     !-----Calculate number of frames-----!
     NFrames=0
     IPos=1 !--Starting position {1-Intel, 0-Portland}
     IPos=IPos+FileHeaderSize
     do
       read(IFU,POS=IPos,err=10,end=10)
       IPos=IPos+NDim*(FrameHeaderSize+FrameSize)
       NFrames=NFrames+1
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
 iwNX=NX
 iwNY=NY
 rwDelX=DelX
 rwDelY=DelY
 rwGrid_Orientation=Grid_Orientation
 iwIType=NDim
 rwXmin=Xmin
 rwYmin=Ymin
 rwStartTime=StartT
 rwEndTime=frTime
 rwDelT=DelT

! iwNFrames=NFrames
 iwNFrames=NFrames+1 !--NFrames is updated outside now
 rwTime=frTime

!----------------------------------------------------------------------!
!     Write file header.
!----------------------------------------------------------------------!
 IPos=1 !--Starting position {1-Intel, 0-Portland}
 write(IFU,POS=IPos) iwNX,iwNY, rwDelX, rwDelY, rwGrid_Orientation, iwIType, rwXmin, rwYmin, rwStartTime, rwEndTime, rwDelT

!----------------------------------------------------------------------!
!     Write last frame header.
!----------------------------------------------------------------------!
 IPos=IPos+FileHeaderSize
! IPos=IPos+NDim*(NFrames-1)*(FrameHeaderSize+FrameSize) !--Skip all previous frames
 IPos=IPos+NDim*NFrames*(FrameHeaderSize+FrameSize) !--Skip all previous frames (NFrames is updated outside now)
 write(IFU,POS=IPos) iwNFrames, iwNFrames, rwTime

!----------------------------------------------------------------------!
!     Write frame.
!----------------------------------------------------------------------!
 IPos=IPos+FrameHeaderSize
 rwa(1:NX*NY)=V1(1:NX*NY)
 write(IFU,POS=IPos) rwa(1:NX*NY)

!----------------------------------------------------------------------!
!     For 2D write last header and frame of 2nd variable.
!----------------------------------------------------------------------!
 if (NDim==2) then
   !-----Write frame header-----!
   IPos=IPos+FrameSize
   write(IFU,POS=IPos) iwNFrames, iwNFrames, rwTime

   !-----Write frame-----!
   IPos=IPos+FrameHeaderSize
   rwa(1:NX*NY)=V2(1:NX*NY)
   write(IFU,POS=IPos) rwa(1:NX*NY)
 endif


 99 CONTINUE
!----------------------------------------------------------------------!
!     Close file.
!----------------------------------------------------------------------!
 CLOSE(IFU)
!------------------------------------------------------------------------------C
RETURN
END
!======================================================================C 
