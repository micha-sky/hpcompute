
!==============================================================================!
SUBROUTINE RDUpperSLDepth
!------------------------------------------------------------------------------!
!
!     RDUpperSLDepth: ReaD Upper Surface Layer Depth.
!
!------------------------------------------------------------------------------!
  USE FILINX
  USE DRLINK 
  USE LANDSF  
  use intersub 
  use solvar
!------------------------------------------------------------------------------!
  IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
!------------------------------------------------------------------------------!
  CHARACTER str*500, Units*80, Lname*80  
  integer ivl(4)
!------------------------------------------------------------------------------!
   if( iSolve(2) +iSolve(3) +iSolve(4) == 0 ) return  
!------------------------------------------------------------------------------!
!     Write header to output file.
!------------------------------------------------------------------------------!
   write(IWR,'(/A)') ' Upper Sediment Layer Depths'
   write(IWR,'(A)')  ' ---------------------------'
!------------------------------------------------------------------------------!
!     Found 'Bottom Surface of Upper Layer' record.
!------------------------------------------------------------------------------!
   read(IRD,*) nLNS
!------------------------------------------------------------------------------!
!     Loop over the number of bottom layer records sets.
!------------------------------------------------------------------------------!
   do iLNS = 1,nLNS
     read(IRD,'(A)') str
     call LCASE( str )
     istart=1
!------------------------------------------------------------------------------!
!     Read the upper layer bottom depth.
!------------------------------------------------------------------------------!
     call RDDPR(istart,icomma,str,Value)
     call RDCHR(istart,icomma,str,Units)
     call RDUNIT(Units,Value)
!------------------------------------------------------------------------------!
!     Read range of upper layer bottom depth.
!------------------------------------------------------------------------------!  
     call ReadRange(',',istart,str, GRange)
!---------------------------------------------------------------------------C
!     Assign upper layer bottom depths.
!---------------------------------------------------------------------------C 
     call AssignVariableInRegionR8(Value,ZSD);
!-----------------------------------------------------------------------------C
!     Write upper layer bottom depths to output file.
!-----------------------------------------------------------------------------C
     WRITE (IWR,'(A,1PE11.4)')                                             & 
            ' Depth of Upper Sediment Layer,             m: ',VALUE  
     call WriteRange(GRange)
!------------------------------------------------------------------------------!
   enddo
!------------------------------------------------------------------------------!
!     End of RDUpperSLDepth group. 
!------------------------------------------------------------------------------!
  RETURN
END
!==============================================================================!
    
!==============================================================================!
SUBROUTINE RDROCK
!------------------------------------------------------------------------------!
!
!     RDROCK: ReaD bottom ROCK types input group.
!
!------------------------------------------------------------------------------!
  USE DRCATC  
  USE DRLINK
  USE DRSEDT
  USE FILINX 
  USE INPVAR 
  USE PLZONE
  USE PROPER 
  USE SOLVAR
  USE SPECIE
  USE intersub
   ! subroutine AssignVariableInRegionI4(Value,V)  
   !   integer*4 Value 
   !   integer*4, Pointer :: V(:)
   ! end subroutine
  !end interface
!------------------------------------------------------------------------------!
  CHARACTER*500 str
  CHARACTER*120 NOTES
  CHARACTER*80 AType
!------------------------------------------------------------------------------!
   if( iSolve(2) +iSolve(3) +iSolve(4) == 0 ) return  
!----------------------------------------------------------------------!
!     Write header to output file.
!----------------------------------------------------------------------!
   write(IWR,'(/A)') ' Bottom Rock Types'
   write(IWR,'(A)')  ' -----------------'
!----------------------------------------------------------------------!
!     Read Number of Note Lines.
!----------------------------------------------------------------------!
   read(IRD,*) NUMLN
!----------------------------------------------------------------------!
!     Read Notes.
!----------------------------------------------------------------------!
   write(IWR,'(A)') ' Rock Notes: '
   do N=1, NUMLN
     read(IRD,'(A)') NOTES
     write(IWR,'(A)') NOTES
   enddo

!----------------------------------------------------------------------!
!     Read Number of Lines.
!----------------------------------------------------------------------!
   read(IRD,*) nLNS
!----------------------------------------------------------------------!
!     ALLOCATING ARRAYS.
!----------------------------------------------------------------------!
   ALLOCATE(ROCTYP(nLNS)); ROCTYP=' ';
   ALLOCATE(RHOS(nLNS),POR(nLNS)); RHOS = 2650.d0; POR = 0.3d0;

   if( ISOLVE(3)+ISOLVE(4) > 0 ) then
     ALLOCATE(PCSD(nLNS),PCSB(nLNS),EXSD(nLNS),EXSB(nLNS));
     PCSD=0.d0; PCSB=0.d0; EXSD=0.d0; EXSB=0.d0;

     ALLOCATE(SMDB(nLNS),SMDL(nLNS),TPKD(nLNS),PCKD(nLNS),HFLF(nLNS));
     SMDB=0.d0; SMDL=0.d0; TPKD=0.d0; PCKD=0.d0; HFLF=1.d20;
   endif
!----------------------------------------------------------------------!

   NROCK = 0
   write(IWR,'(/A)') ' Bottom Rock Types and Domains'
!----------------------------------------------------------------------!
!     Loop over Bottom Rock Type records.
!----------------------------------------------------------------------!
   do iLNS = 1,nLNS
     read(IRD,'(A)') str
     call LCASE( str ) 
     istart=1
     if (iLNS/=1) write(IWR,*) 
   !----------------------------------------------------------------------!
   !     Read Bottom Rock Type.
   !----------------------------------------------------------------------!
     AType(1:) = ' ' 
     call RDCHR(istart,icomma,str,AType) 
   !----------------------------------------------------------------------!
   !     Look if previously defined.
   !----------------------------------------------------------------------!
     do k=1,NROCK
       if (AType==ROCTYP(k)) then
         IROCK = k
         GOTO 301
       endif
     enddo
   !----------------------------------------------------------------------!
   !     If a previously undefined bottom rock type, increment the number
   !     of bottom rock types (NROCK) and check to ensure parameter
   !     limit is not exceeded.
   !----------------------------------------------------------------------!
     NROCK = NROCK +1
     ROCTYP(NROCK) = AType
     IROCK = NROCK
 301 CONTINUE

   !----------------------------------------------------------------------!
   !     Read Range.
   !----------------------------------------------------------------------!
     call ReadRange(',',istart,str,GRange)
   !----------------------------------------------------------------------!
   !     Write Bottom Rock Type information to output file.
   !----------------------------------------------------------------------!
     write(IWR,'(A,I0)') '   Rock Number:  ',IROCK
     write(IWR,'(2A)')   '   Name: ', trim(ROCTYP(IROCK))
   !----------------------------------------------------------------------!
   !     Write range.
   !----------------------------------------------------------------------!
     call WriteRange(GRange)
   !----------------------------------------------------------------------!
   !     Assign Bottom Rock Types to cells.
   !----------------------------------------------------------------------!
     call AssignVariableInRegionI4(IROCK,IZ)

!----------------------------------------------------------------------!
   enddo
!----------------------------------------------------------------------C 
!     End of RDROCK group. 
!----------------------------------------------------------------------C 
RETURN
END
!==============================================================================!

!==============================================================================!
SUBROUTINE RDBTTP
!----------------------------------------------------------------------!
!
!     RDBTTP: ReaD BoTtom Types input group.
!
!----------------------------------------------------------------------! 
USE DRLINK
USE FILINX 
USE INPVAR  
USE intersub 
USE PLZONE  
USE SOLVAR
!----------------------------------------------------------------------C 
IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
!----------------------------------------------------------------------!
!  interface  
!    subroutine AssignVariableInRegionI4(Value,V)  
!      integer*4 Value 
!      integer*4, Pointer :: V(:)
!    end subroutine
!  end interface
!------------------------------------------------------------------------------!
CHARACTER*500 str
CHARACTER*120 NOTES
CHARACTER*80 AType
!----------------------------------------------------------------------!
!------------------------------------------------------------------------------!
   if( iSolve(2) +iSolve(3) +iSolve(4) == 0 ) return  
!----------------------------------------------------------------------!
!     Write header to output file.
!----------------------------------------------------------------------!
 write(IWR,'(/A)') ' Bottom Types'
 write(IWR,'(A)')  ' ------------'
!----------------------------------------------------------------------!
!     Read Number of Note Lines.
!----------------------------------------------------------------------!
 read(IRD,*) NUMLN
!----------------------------------------------------------------------!
!     Read Notes.
!----------------------------------------------------------------------!
 write(IWR,'(A)') ' Bottom Type Notes: '
 do N=1, NUMLN
   read(IRD,'(A)') NOTES
   write(IWR,'(A)') NOTES
 enddo

!----------------------------------------------------------------------!
!     Read Number of Lines.
!----------------------------------------------------------------------!
 read(IRD,*) nLNS
!----------------------------------------------------------------------!
!     ALLOCATING ARRAYS.
!----------------------------------------------------------------------!
! call AllocateBottomTypesArrays(nLNS)
 ALLOCATE(PLNTYP(nLNS)); PLNTYP=' ';

! nk = nLNS;
! ALLOCATE(ROCTYP(nk)); ROCTYP=' ';
! ALLOCATE(RHOS(nk),POR(nk)); RHOS=2650.d0; POR=0.d0;

! if (ISOLVE(4)+ISOLVE(5)>0) then
!   ALLOCATE(PCSD(nk),PCSB(nk),EXSD(nk),EXSB(nk));
!   PCSD=0;PCSB=0.d0;EXSD=0.d0;EXSB=0.d0;

!   ALLOCATE(SMDB(nk),SMDL(nk),TPKD(nk),PCKD(nk),HFLF(nk));
!   SMDB=0.d0;SMDL=0.d0;TPKD=0.d0;PCKD=0.d0;HFLF=1.d20;
! endif
!----------------------------------------------------------------------C 

 NPLANT = 0
 write(IWR,'(/A)') ' Bottom Types and Domains'
!----------------------------------------------------------------------!
!     Loop over Bottom Type records.
!----------------------------------------------------------------------!
 do iLNS = 1,nLNS
   read(IRD,'(A)') str
   istart=1
   if (iLNS/=1) write(IWR,*)
   !----------------------------------------------------------------------!
   !     Read Bottom Type.
   !----------------------------------------------------------------------!
   ires = RDStr(',',istart,str,AType)
   !----------------------------------------------------------------------!
   !     Look if previously defined.
   !----------------------------------------------------------------------!
   do k=1,NPLANT
     if (AType==PLNTYP(k)) then
       IPLANT = k
       GOTO 301
     endif
   enddo
   !----------------------------------------------------------------------!
   !     If a previously undefined bottom type, increment the number
   !     of bottom types (NPLANT) and check to ensure parameter
   !     limit is not exceeded.
   !----------------------------------------------------------------------!
   NPLANT = NPLANT +1
   PLNTYP(NPLANT) = AType
   IPLANT = NPLANT
   301 CONTINUE
   !----------------------------------------------------------------------!
   !     Read Range.
   !----------------------------------------------------------------------!
   call ReadRange(',',istart,str,GRange)
   !----------------------------------------------------------------------!
   !     Write Bottom Type information to output file.
   !----------------------------------------------------------------------!
   write(IWR,'(A,I0)') '   Bottom Type Number:  ',IPLANT
   write(IWR,'(2A)')   '   Name: ', trim(PLNTYP(IPLANT))
   !----------------------------------------------------------------------!
   !     Write range.
   !----------------------------------------------------------------------!
   call WriteRange(GRange)
   !----------------------------------------------------------------------!
   !     Assign Bottom Types to cells.
   !----------------------------------------------------------------------!
   call AssignVariableInRegionI4(IPLANT,IZPL)

!----------------------------------------------------------------------!
 enddo
!----------------------------------------------------------------------!
!----------------------------------------------------------------------C 
!     End of RDBTTP group. 
!----------------------------------------------------------------------C 
RETURN
END
!==============================================================================!

!==============================================================================!
SUBROUTINE RDMECH
!----------------------------------------------------------------------!
!
!     RDMECH: ReaD MECHanical properties input group.
!
!----------------------------------------------------------------------!
USE DRLINK
USE FILINX 
USE INPVAR
USE intersub  
USE PROPER
USE SOLVAR
!----------------------------------------------------------------------!
IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
!----------------------------------------------------------------------!
CHARACTER*500 str,str1
CHARACTER*80 Units,AType
!------------------------------------------------------------------------------!
   if( iSolve(2) +iSolve(3) +iSolve(4) == 0 ) return  
!----------------------------------------------------------------------!
!     Write header to output file.
!----------------------------------------------------------------------!
 write(IWR,'(/A)') ' Mechanical Properties'
 write(IWR,'(A)')  ' ---------------------'
!----------------------------------------------------------------------!
!     Read Number of Lines.
!----------------------------------------------------------------------!
 read(IRD,*) nLNS

!----------------------------------------------------------------------!
!     Loop over Mechanical Properties records.
!----------------------------------------------------------------------!
 do iLNS = 1,nLNS
   read(IRD,'(A)') str
   istart=1
   if (iLNS/=1) write(IWR,*)
   !----------------------------------------------------------------------!
   !     Read rock name.
   !----------------------------------------------------------------------!
   ires = RDStr(',',istart,str,AType)
   do k = 1, NROCK !--Looking for a corresponding Bottom Rock Type name
     if (AType==ROCTYP(k)) then
       IROCK = k
       GOTO 860
     endif
   enddo
   call Msg(0,IWR,' INPUT ERROR! - No matching Rock Type!')
   call Msg(0,IWR,' Type we''re looking for is: '//trim(AType))
   call Msg(0,IWR,' Run aborting...')
   STOP
   860 CONTINUE
   write(IWR,'(2A)') ' Rock type: ', trim(AType)
   !----------------------------------------------------------------------!
   !     Read the porosity and storativity.
   !----------------------------------------------------------------------!
   ires = RDR8(',',istart,str,POR(IROCK))
   !----------------------------------------------------------------------!
   !     Read 'rock density'.
   !----------------------------------------------------------------------!
   ires = RDR8(',',istart,str,RHOS(IROCK))
   ires = RDStr(',',istart,str,Units)
   call ConvToSI(RHOS(IROCK),Units)
   !----------------------------------------------------------------------!
   !     Output information.
   !----------------------------------------------------------------------!
   write(IWR,'(A,1PG12.5)') '   Porosity            : ', POR(IROCK)
   write(IWR,'(A,1PG12.5)') '   Rock Density, kg/m^3: ', RHOS(IROCK)

 !----------------------------------------------------------------------!
 enddo
!----------------------------------------------------------------------!
!----------------------------------------------------------------------C 
!     End of RDMECH group. 
!----------------------------------------------------------------------C 
RETURN
END
!==============================================================================!

!==============================================================================!
SUBROUTINE RDBTPR
!----------------------------------------------------------------------!
!
!     RDBTPR: ReaD BoTtom PRoperties input group.
!
!----------------------------------------------------------------------!
USE DRCATC  
USE FILINX 
USE INPVAR 
USE intersub  

USE DIFF
USE SOLVAR
!----------------------------------------------------------------------!
IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
!----------------------------------------------------------------------!
LOGICAL fDIFF1,fDIFF2,fDIFF3
CHARACTER*500 str,str1
CHARACTER*80 Units,AType,PEROPT
CHARACTER*80 cfStr !--string defining which variables are read from files and which are constant
CHARACTER*150 zFile,zFile1,zFile2
!----------------------------------------------------------------------!
!------------------------------------------------------------------------------!
 if( iSolve(2) +iSolve(3) +iSolve(4) == 0 ) return  
!----------------------------------------------------------------------!
!     Write header to output file.
!----------------------------------------------------------------------!
 write(IWR,'(/A)') ' Bottom Properties'
 write(IWR,'(A)')  ' -----------------'

!----------------------------------------------------------------------!
!     ALLOCATING ARRAYS.
!----------------------------------------------------------------------!
 ALLOCATE(Erod(NPLANT),DFSL(NPLANT),Fric(NPLANT),Rough(NPLANT));
 Erod=0d0; DFSL=0d0; Fric=0d0; Rough=0d0;
 allocate(IDiff(NPLANT),IFRIC(NPLANT)); IDiff = 0; IFRIC = 0;
!----------------------------------------------------------------------!
!     Read Number of Lines.
!----------------------------------------------------------------------!
 read(IRD,*) nLNS

 fDIFF3=.false.
!----------------------------------------------------------------------!
!     Loop over Bottom Type records.
!----------------------------------------------------------------------!
 do iLNS = 1,nLNS
   read(IRD,'(A)') str
   istart=1
   if (iLNS/=1) write(IWR,*)
   !----------------------------------------------------------------------!
   !     Read bottom type name.
   !----------------------------------------------------------------------!
   ires = RDStr(',',istart,str,AType)
   do k = 1, NPLANT !--Looking for a corresponding Bottom Type name
     if (AType==PLNTYP(k)) then
       IPLANT = k
       GOTO 860
     endif
   enddo
   call Msg(0,IWR,' INPUT ERROR! - No matching Bottom Type!')
   call Msg(0,IWR,' Type we''re looking for is: '//trim(AType))
   call Msg(0,IWR,' Run aborting...')
   STOP
   860 CONTINUE
   write(IWR,'(2A)') ' Bottom Type: ', trim(AType)
   !----------------------------------------------------------------------!
   !     Read erodibility coefficient.
   !----------------------------------------------------------------------!
     ires = RDR8(',',istart,str,EROD(IPLANT))
     write(IWR,'(A,1PG12.5)') '   Erodibility Coefficient: ', EROD(IPLANT)

   !----------------------------------------------------------------------!
   !     Read Mixing Parameters in the Water Column (Eddy Viscosity
   !     Coefficient),
   !     or
   !     Read depth-averaged turbulent diffusivity.
   !----------------------------------------------------------------------!
   fDIFF1=.false.
   fDIFF2=.false.
     IDIFF(IPLANT)=2
       ires = RDR8(',',istart,str,DFSL(IPLANT))
       ires = RDStr(',',istart,str,Units)
       call ConvToSI(DFSL(IPLANT),Units)
       fDIFF2 = (DFSL(IPLANT) > SMALL)
   !----------------------------------------------------------------------!
   !     Output Information.
   !----------------------------------------------------------------------!
     write(IWR,'(A)')'   Diffusion: Fixed '
     write(IWR,'(A,1PG12.5)') '   Turbulent Diffusion Coefficient, m^2/s: ', DFSL(IPLANT)
   !----------------------------------------------------------------------!
   !     Read Hydraulic Resistance Function.
   !----------------------------------------------------------------------!
   ires = RDStr(',',istart,str,PEROPT)
   !----------------------------------------------------------------------!
   !     Read Friction Coefficient.
   !----------------------------------------------------------------------!
     ires = RDR8(',',istart,str,FRIC(IPLANT))
     ires = RDStr(',',istart,str,Units)
     call ConvToSI(FRIC(IPLANT),Units)

   if (PEROPT(1:5)=='darcy') then !--Cb=rFr
     IFRIC(IPLANT) = 1
     write(IWR,'(A)') '   Hydraulic Resistance Function: Darcy-Weisbach Equation'
     write(IWR,'(A,1PG12.5)') '   Darcy-Weisbach Friction Coefficient: ', FRIC(IPLANT)
   elseif (PEROPT(1:5)=='chezy') then
     IFRIC(IPLANT) = 2
     write(IWR,'(A)') '   Hydraulic Resistance Function: Chezy Law'
     write(IWR,'(A,1PG12.5)') '   Chezy Friction Coefficient, m^0.5/s: ', FRIC(IPLANT)
     FRIC(IPLANT)=GRAV/max(FRIC(IPLANT)*FRIC(IPLANT),SMALL)     !--Prepare for Cb calculation, now Cb=rFr
   elseif (PEROPT(1:7)=='manning') then
     IFRIC(IPLANT) = 3
     write(IWR,'(A)') '   Hydraulic Resistance Function: Manning Equation'
     write(IWR,'(A,1PG12.5)') '   Manning''s Roughness Coefficient, s/m^0.33: ', FRIC(IPLANT)
     FRIC(IPLANT)=GRAV*FRIC(IPLANT)*FRIC(IPLANT) !--Prepare for Cb calculation, now Cb=rFr/h**(1/3)
   endif
   !----------------------------------------------------------------------!
   !     Read roughness parameter.
   !----------------------------------------------------------------------!
     ires = RDR8(',',istart,str,ROUGH(IPLANT))
     ires = RDStr(',',istart,str,Units)
     call ConvToSI(ROUGH(IPLANT),Units)
     write(IWR,'(1A,1PG12.5)')'   Roughness parameter, m: ',ROUGH(IPLANT)

 !----------------------------------------------------------------------!
 enddo


 fDIFF=fDIFF .and. fDIFF3
!----------------------------------------------------------------------!
!----------------------------------------------------------------------C 
!     End of RDBTPR group. 
!----------------------------------------------------------------------C 
RETURN
END
!==============================================================================!

!==============================================================================!
SUBROUTINE RDSEPR
!----------------------------------------------------------------------C 
!
!     RDSEPR: ReaD SEdiment PRoperties input group.
!
!----------------------------------------------------------------------C  
USE DRSEDT
USE FILINX
USE intersub  
USE SOLVAR
!----------------------------------------------------------------------C 
IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
!----------------------------------------------------------------------C 
CHARACTER*500 str,str1
CHARACTER*80 Units,Type
CHARACTER*150 TABLEFN
INTEGER*4 KSMAF(2,1)
!----------------------------------------------------------------------C 
!------------------------------------------------------------------------------!
   if( iSolve(2) +iSolve(3) +iSolve(4) == 0 ) return  
!----------------------------------------------------------------------!
!     Write header to output file.
!----------------------------------------------------------------------!
 write(IWR,'(/A)') ' Sediment Properties'
 write(IWR,'(A)')  ' -------------------'
!----------------------------------------------------------------------!
!     Read Number of Lines.
!----------------------------------------------------------------------!
 read(IRD,*) nLNS

!----------------------------------------------------------------------!
!     Loop over Sediment Properties records.
!----------------------------------------------------------------------!
 do iLNS = 1,nLNS
   read(IRD,'(A)') str
   istart=1
   !----------------------------------------------------------------------!
   !     Read sediment diameter.
   !----------------------------------------------------------------------!
   ires = RDR8(',',istart,str,DSE)
   ires = RDR8(',',istart,str,D90)
   ires = RDStr(',',istart,str,Units)
   call ConvToSI(DSE,Units)
   call ConvToSI(D90,Units)
   !----------------------------------------------------------------------!
   !     Read density of sediment particle.
   !----------------------------------------------------------------------!
   ires = RDR8(',',istart,str,RHOSE)
   ires = RDStr(',',istart,str,Units)
   call ConvToSI(RHOSE,Units)

   !----------------------------------------------------------------------!
   !     Read 'sediment transport capacity option'.
   !----------------------------------------------------------------------!
   ires = RDStr(',',istart,str,str1)

   !----------------------------------------------------------------------!
   !     No equation, no erosion\deposition (dye for example).
   !----------------------------------------------------------------------!
   if (str1(1:4)=='none') then
     ILOAD = 0
     write(IWR,'(A)') ' Transport Capacity Equation: None'
   !----------------------------------------------------------------------!
   !     Yalin bedload equation.
   !----------------------------------------------------------------------!
   elseif (str1(1:5)=='yalin') then
     ILOAD = 1
     write(IWR,'(A)') ' Transport Capacity Equation: Yalin'
   !----------------------------------------------------------------------!
   !     Modified Yalin bedload equation.
   !----------------------------------------------------------------------!
   elseif (str1(1:14)=='modified yalin') then
     ILOAD = 2
     write(IWR,'(A)') ' Transport Capacity Equation: modified Yalin'
   !----------------------------------------------------------------------!
   !     Engelund-Hansen total load equation.
   !----------------------------------------------------------------------!
   elseif (str1(1:15)=='engelund-hansen') then
     ILOAD = 3
     write(IWR,'(A)') ' Transport Capacity Equation: Engelund-Hansen'
   !----------------------------------------------------------------------!
   !     Einstein-Brown bedload equation.
   !----------------------------------------------------------------------!
   elseif (str1(1:14)=='einstein-brown') then
     ILOAD = 4
     write(IWR,'(A)') ' Transport Capacity Equation: Einstein-Brown'
   !----------------------------------------------------------------------!
   !     Bagnold total load equation.
   !----------------------------------------------------------------------!
   elseif (str1(1:7)=='bagnold') then
     ILOAD = 5
     write(IWR,'(A)') ' Transport Capacity Equation: Bagnold'
   !----------------------------------------------------------------------!
   !     Ackers-White total load equation.
   !----------------------------------------------------------------------!
   elseif (str1(1:12)=='ackers-white') then
     ILOAD = 6
     write(IWR,'(A)') ' Transport Capacity Equation: Ackers-White'
   !----------------------------------------------------------------------!
   !     Rijn total load equation.
   !----------------------------------------------------------------------!
   elseif (str1(1:9)=='rijn-1993') then
     ILOAD = 7
     write(IWR,'(A)') ' Transport Capacity Equation: van Rijn (1993)'
   !----------------------------------------------------------------------!
   !     Rijn total load equation.
   !----------------------------------------------------------------------!
   elseif (str1(1:15)=='simplified rijn') then
     ILOAD = 8
     write(IWR,'(A)') ' Transport Capacity Equation: simplified van Rijn'
   !----------------------------------------------------------------------!
   !     Ariathurai-Krone load equation (cohesive sediment).
   !----------------------------------------------------------------------!
   elseif (str1(1:15)=='cohesive') then
     ILOAD = 9
     write(IWR,'(A)') ' Transport Capacity Equation: Ariathurai-Krone (cohesive sediment)'
     !----------------------------------------------------------------------!
     !     Read cohesive params.
     !----------------------------------------------------------------------!
     ires = RDR8(',',istart,str,TAY_D); if (ires>0) TAY_D=0.1; !--Read or set default value if missing
     ires = RDR8(',',istart,str,TAY_E); if (ires>0) TAY_E=2.5; !--Read or set default value if missing
     ires = RDStr(',',istart,str,Units)
     call ConvToSI(TAY_D,Units)
     call ConvToSI(TAY_E,Units)
     TAY_D = max(TAY_D,0.d0) !--TAY_D must be non-negative
     TAY_E = max(TAY_E,1.d-10) !--TAY_E must be positive
     ires = RDR8(',',istart,str,WCM); if (ires>0) WCM=3.d-5; !--Read or set default value if missing
     write(IWR,'(A,1PG12.5)') '  Critical shear stress for deposition, N/m^2: ',TAY_D
     write(IWR,'(A,1PG12.5)') '  Critical shear stress for erosion,    N/m^2: ',TAY_E
     write(IWR,'(A,1PG12.5)') '  Erosion constant (M): ',WCM
   !----------------------------------------------------------------------!
   !     Bijker total load equation (W&C).
   !----------------------------------------------------------------------!
   elseif (str1(1:6)=='bijker') then
     ILOAD = 10
     write(IWR,'(A)') ' Transport Capacity Equation: Bijker (W&C)'
   !----------------------------------------------------------------------!
   !     Van Rijn total load equation (W&C).
   !----------------------------------------------------------------------!
   elseif (str1(1:11)=='wave & rijn') then
     ILOAD = 11
     write(IWR,'(A)') ' Transport Capacity Equation: Van Rijn (W&C)'
   !----------------------------------------------------------------------!
   !     Dibajnia-Watanabe total load equation (W&C).
   !----------------------------------------------------------------------!
   elseif (str1(1:17)=='dibajnia-watanabe') then
     ILOAD = 12
     write(IWR,'(A)') ' Transport Capacity Equation: Dibajnia-Watanabe (W&C)'
   !----------------------------------------------------------------------!
   !     Meuer-Peter & Muler bed load equation (W&C).
   !----------------------------------------------------------------------!
   elseif (str1(1:11)=='meyer-peter') then
     ILOAD = 13
     write(IWR,'(A)') ' Transport Capacity Equation: Meyer-Peter & Muler (W&C)'
   !----------------------------------------------------------------------!
   !     Nielson bed load equation (W&C).
   !----------------------------------------------------------------------!
   elseif (str1(1:7)=='nielson') then
     ILOAD = 14
     write(IWR,'(A)') ' Transport Capacity Equation: Nielson (W&C)'
   !----------------------------------------------------------------------!
   !     Soulsby bed load equation (W&C).
   !----------------------------------------------------------------------!
   elseif (str1(1:7)=='soulsby') then
     ILOAD = 15
     write(IWR,'(A)') ' Transport Capacity Equation: Soulsby (W&C)'
   !----------------------------------------------------------------------!
   !     Camenen & Larson bed load equation (W&C).
   !----------------------------------------------------------------------!
   elseif (str1(1:16)=='camenen & larson') then
     ILOAD = 16
     write(IWR,'(A)') ' Transport Capacity Equation: Camenen & Larson (W&C)'
   !----------------------------------------------------------------------!
   !     Yang's bed load equation (W&C).
   !----------------------------------------------------------------------!
   elseif (str1(1:4)=='yang') then
     ILOAD = 17
     write(IWR,'(A)') ' Transport Capacity Equation: Yang (W&C)'
   !----------------------------------------------------------------------!
   !     Yang's bed load equation (W&C).
   !----------------------------------------------------------------------!
   elseif (str1(1:10)=='river yang') then
     ILOAD = 18
     write(IWR,'(A)') ' Transport Capacity Equation: Yang for river'
   !----------------------------------------------------------------------!
   !     Van Rijn total load equation similar mike21.
   !----------------------------------------------------------------------!
   elseif (str1(1:13)=='mike21 & rijn') then
     ILOAD = 19
     write(IWR,'(A)') ' Transport Capacity Equation: Van Rijn similar mike21'
   !----------------------------------------------------------------------!
   !     Van Rijn total load equation for silt&sand.
   !----------------------------------------------------------------------!
   elseif (str1(1:9)=='rijn-2004') then
     ILOAD = 20
     write(IWR,'(A)') ' Transport Capacity Equation: van Rijn (2004)'
   else
     ILOAD = 1
     write(IWR,'(A)') ' Transport Capacity Equation: Yalin'
   endif

   ires = RDStr(',',istart,str,str1)
   !----------------------------------------------------------------------!
   !     Read 'Shear stress mode' if using 'wave & rijn' formula.
   !----------------------------------------------------------------------!
   if (ILOAD==11) then
     !----------------------------------------------------------------------!
     !     Log-Velocity Profile shear stress mode.
     !----------------------------------------------------------------------!
     if (str1(1:12)=='log-velocity') then
       ISHMODE = 1
       write(IWR,'(A)') '  Shear stress Mode: Log-Velocity Profile'
     !----------------------------------------------------------------------!
     !     Manning shear stress mode.
     !----------------------------------------------------------------------!
     elseif (str1(1:7)=='manning') then
       ISHMODE = 2
       write(IWR,'(A)') '  Shear stress Mode: Manning'
     !----------------------------------------------------------------------!
     !     Bijker shear stress mode.
     !----------------------------------------------------------------------!
     elseif (str1(1:11)=='bijker-type') then
       ISHMODE = 3
       write(IWR,'(A)') '  Shear stress Mode: Bijker-type'
     else
       ISHMODE = 3
       write(IWR,'(A)') '  Shear stress Mode: Bijker-type'
     endif
   endif

   !----------------------------------------------------------------------!
   !     Output Information.
   !----------------------------------------------------------------------!
   if (ILOAD>0) then
     write(IWR,'(A,1PG12.5)') ' Diameter of Sediment Particle (D50), m: ',DSE
     write(IWR,'(A,1PG12.5)') ' Diameter of Sediment Particle (D90), m: ',D90
     write(IWR,'(A,1PG12.5)') ' Density of Sediment Particle,   kg/m^3: ',RHOSE
   endif

   !----------------------------------------------------------------------!
   !     Read Morphological Acceleration Factor if needed.
   !----------------------------------------------------------------------!
!   NMAF=0
!   MAF=1.d0
!   InvMAF=1.d0
!   if (ILOAD==0) CYCLE
   !----------------------------------------------------------------------!
   !     Read temporal variation.
   !----------------------------------------------------------------------!
!   ITEMP = 0 !--Current temporal variation
!   ires = RDStr(',',istart,str,Type)
!   if (Type(1:8)=='constant') then
!     ITEMP = 1
!     ires = RDR8(',',istart,str,MAF)
!     InvMAF=0.d0; if (MAF>SMALL) InvMAF=1.d0/MAF
!   elseif (Type(1:11) == 'tabularfile' ) then
!     ITEMP = 21
     !----------------------------------------------------------------------!
     !     Read external file containing table of values.
     !----------------------------------------------------------------------!
!     ires = RDStr(',',istart,str,TABLEFN,1)

!     call CalcTextFileLength(TABLEFN,NLIN)
     !----------------------------------------------------------------------!
     !     ALLOCATING ARRAYS.
     !----------------------------------------------------------------------!
!     nk=NLIN
!     ALLOCATE(SMAF(nk), SInvMAF(nk), SMAF_TM(nk)); SMAF=1.d0; SInvMAF=1.d0;  SMAF_TM=0.d0;
     !----------------------------------------------------------------------!
     !     Reading table.
     !----------------------------------------------------------------------!
!     IFU=1001
!     OPEN(UNIT=IFU, FILE=TABLEFN)
!     call ReadTableCOASTOX(1,IFU,NLIN,1,NMAF,KSMAF,'',SMAF_TM,1,SMAF,SMAF,SMAF,.false.)
!     CLOSE(IFU)
!   elseif (Type(1:7)=='tabular') then
!     ITEMP = 2
!   endif

   !----------------------------------------------------------------------!
   !     Output Information.
   !----------------------------------------------------------------------!
!   if (ILOAD>0) then
!     write(IWR,'(A)') ' Morphology Acceleration Factor'
!     if (ITEMP==1) then
!       write(IWR,'(A)') '   Temporal Variation: Constant'
!       write(IWR,'(A,1PG12.5)') '   MAF value: ',MAF
!     elseif (ITEMP==21) then
!       write(IWR,'(2A)') '   Temporal Variation: Tabular from File - ', trim(TABLEFN)
!     elseif (ITEMP==2) then
!       write(IWR,'(A)') '   Temporal Variation: Tabular'
!     endif
!   endif
   !----------------------------------------------------------------------!
   !     Read source tables.
   !----------------------------------------------------------------------!
!   if (ITEMP/=2) GOTO 900

   !----------------------------------------------------------------------!
   !     Length of table.
   !----------------------------------------------------------------------!
!   read(IRD,*) NLIN
   !----------------------------------------------------------------------!
   !     ALLOCATING ARRAYS.
   !----------------------------------------------------------------------!
!   nk=NLIN
!   ALLOCATE(SMAF(nk), SInvMAF(nk), SMAF_TM(nk)); SMAF=1.d0; SInvMAF=1.d0;  SMAF_TM=0.d0;
   !----------------------------------------------------------------------!
   !     Reading table.
   !----------------------------------------------------------------------!
!   call ReadTableCOASTOX(1,IRD,NLIN,1,NMAF,KSMAF,'',SMAF_TM,1,SMAF,SMAF,SMAF,.false.)


!   900 CONTINUE

!   if (ITEMP<=1) CYCLE
   !----------------------------------------------------------------------!
   !     Write source table header to output file.
   !----------------------------------------------------------------------!
!   write(IWR,'(A)') '   MAF Table'
   !----------------------------------------------------------------------!
   !     Write source table to output file.
   !----------------------------------------------------------------------!
!   do i = 1,NMAF
!     write(IWR,'(A,I4,2(A,1PG12.5))') '   Entry No.',i, ';  Time, s: ',SMAF_TM(i), ';  MAF value: ',SMAF(i)
!   enddo

 !----------------------------------------------------------------------!
 enddo
!----------------------------------------------------------------------C 
!     End of RDSEPR group. 
!----------------------------------------------------------------------C 
RETURN
END
!==============================================================================!

!==============================================================================!
SUBROUTINE RDSPPR
!----------------------------------------------------------------------C 
!
!     RDSPPR: ReaD SPecies PRoperties input group.
!
!----------------------------------------------------------------------C 
USE DRCATC 
USE FILINX
USE intersub  
USE INPVAR
USE SOLVAR 
USE SPECIE
!----------------------------------------------------------------------C 
IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
!----------------------------------------------------------------------C 
CHARACTER*500 str
CHARACTER*80 AType,SPCNAM,Units
!----------------------------------------------------------------------C 
!------------------------------------------------------------------------------!
   if( iSolve(3) +iSolve(4) == 0 ) return  
!----------------------------------------------------------------------!
!     Write header to output file.
!----------------------------------------------------------------------!
 write(IWR,'(/A)') ' Species Properties'
 write(IWR,'(A)')  ' ------------------'
!----------------------------------------------------------------------!
!     Check for matching rock type.
!----------------------------------------------------------------------!
 read(IRD,*) nLNS
 do iLNS = 1,nLNS
   read(IRD,'(A)') str
   istart=1
   if (iLNS/=1) write(IWR,*)
   !----------------------------------------------------------------------!
   !     Read rock type name.
   !----------------------------------------------------------------------!
   ires = RDStr(',',istart,str,AType)
   do k = 1, NROCK !--Looking for a corresponding Rock Type name
     if (AType == ROCTYP(k)) then
       IROCK = k
       GOTO 860
     endif
   enddo
   call Msg(0,IWR,' INPUT ERROR! - No matching Rock Type!')
   call Msg(0,IWR,' Type we''re looking for is: '//trim(AType))
   call Msg(0,IWR,' Run aborting...')
   STOP
   860 CONTINUE
   write(IWR,'(2A)') ' Rock Type: ', trim(AType)
   !----------------------------------------------------------------------!
   !     Read species name.
   !----------------------------------------------------------------------!
   ires = RDStr(',',istart,str,SPCNAM)
   write(IWR,'(2A)') '   Species: ', trim(SPCNAM)
   !----------------------------------------------------------------------!
   !     Read sediment-water and sediment-bottom partition coefficient.
   !----------------------------------------------------------------------!
   ires = RDR8(',',istart,str,PCSD(IROCK))
   ires = RDR8(',',istart,str,PCSB(IROCK))
   ires = RDStr(',',istart,str,Units)
   call ConvToSI(PCSD(IROCK),Units)
   call ConvToSI(PCSB(IROCK),Units)
   write(IWR,'(A,1PG12.5)') '   Sediment-Water Partition Coefficient, kg liq./kg sol.: ', PCSD(IROCK)
   write(IWR,'(A,1PG12.5)') '   Water-Bottom Partition Coefficient,   kg liq./kg sol.: ', PCSB(IROCK)
   !----------------------------------------------------------------------!
   !     Read rates of sediment-water and sediment-bottom exchange.
   !----------------------------------------------------------------------!
   ires = RDR8(',',istart,str,EXSD(IROCK))
   ires = RDR8(',',istart,str,EXSB(IROCK))
   ires = RDStr(',',istart,str,Units)
   call ConvToSI(EXSD(IROCK),Units)
   call ConvToSI(EXSB(IROCK),Units)
   write(IWR,'(A,1PG12.5)') '   Sediment-Water Exchange Rate, 1/s: ', EXSD(IROCK)
   write(IWR,'(A,1PG12.5)') '   Water-Bottom Exchange Rate,   1/s: ', EXSB(IROCK)
   !----------------------------------------------------------------------!
   !     Read species half life.
   !----------------------------------------------------------------------!
   ires = RDR8(',',istart,str,HFLF(IROCK))
   ires = RDStr(',',istart,str,Units)
   call ConvToSI(HFLF(IROCK),Units)
   if (HFLF(IROCK) <= SMALL) then
     HFLF(IROCK) = ZERO
   else
     HFLF(IROCK) = 0.6931D+0/HFLF(IROCK)
   endif
   write(IWR,'(A,1PG12.5)') '   Species Decay Constant, 1/s: ', HFLF(IROCK)
!----------------------------------------------------------------------!
 enddo
!
!----------------------------------------------------------------------C 
!     End of RDSPPR group. 
!----------------------------------------------------------------------C 
RETURN
END
!==============================================================================!

!==============================================================================!
SUBROUTINE RDRAIN
!----------------------------------------------------------------------C 
!
!     RDRAIN: ReaD RAIN input group.
!
!----------------------------------------------------------------------C 
USE FILINX 
USE DRRAIN 
USE intersub  
!----------------------------------------------------------------------C 
IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
      INCLUDE 'utils.inc' ! subroutine interfaces
!----------------------------------------------------------------------C 
CHARACTER*500 str,str1
CHARACTER*80 Units,Type
CHARACTER*150 TABLEFN
REAL*8 Value(5)
!----------------------------------------------------------------------C 
!----------------------------------------------------------------------!
!     Write header to output file.
!----------------------------------------------------------------------!
 write(IWR,'(/A)') ' Rainfall Rate'
 write(IWR,'(A)')  ' -------------'
!----------------------------------------------------------------------!
!     Set table counters to zero.
!----------------------------------------------------------------------!
 NSPREC = 0
 NRAIN  = 0

!----------------------------------------------------------------------!
!     Read the number of RAIN input lines.
!----------------------------------------------------------------------!
 read(IRD,*) nLNS

!----------------------------------------------------------------------!
!     ALLOCATING ARRAYS.
!----------------------------------------------------------------------!
 nk=nLNS
 ALLOCATE(IPRCV(nk), KSPREC(2,nk), ISPREC(4,nk)); IPRCV=0; KSPREC=0;

!----------------------------------------------------------------------!
!     Loop over the number of lines.
!----------------------------------------------------------------------!
 do iLNS = 1,nLNS
   read(IRD,'(A)') str
   istart=1
   NSPREC = NSPREC + 1
   !----------------------------------------------------------------------!
   !     Read temporal variation.
   !----------------------------------------------------------------------!
   ITEMP = 0 !--Current temporal variation
   ires = RDStr(',',istart,str,Type)
   if (Type(1:8)=='constant') then
     ITEMP = 1
     !----------------------------------------------------------------------!
     !     REALLOCATING ARRAYS.
     !----------------------------------------------------------------------!
     nk=NRAIN + 2
!     call Reallocate2DR8(NRAIN,3,nk,3,SPREC)  !!#
     call ReallocateR8(NRAIN,nk,SPREC)  !!#
     call ReallocateR8(NRAIN,nk,SPREC_TM)
     !----------------------------------------------------------------------!
     !     Read Rainfall rate, starting and ending time.
     !----------------------------------------------------------------------!
     do i=1,3
       ires = RDR8(',',istart,str,Value(i))
       ires = RDStr(',',istart,str,Units)
       call ConvToSI(Value(i),Units)
     enddo
     !----------------------------------------------------------------------!
     !     Assign Rainfall variables.
     !----------------------------------------------------------------------!
     NRAIN = NRAIN + 1
     KSPREC(1,NSPREC) = NRAIN
     KSPREC(2,NSPREC) = NRAIN+1
     SPREC(NRAIN)   = Value(1)
!     SPREC(NRAIN,1)   = Value(1)
!     SRAIN(NRAIN,2)   = Value(2)
!     SRAIN(NRAIN,3)   = Value(3)
!     SRAIN_TM(NRAIN)  = Value(4)
     NRAIN = NRAIN + 1
     SPREC_TM(NRAIN)  = Value(2)

   elseif (Type(1:11) == 'tabularfile') then
     ITEMP = 21
     !----------------------------------------------------------------------!
     !     Read external file containing table of values.
     !----------------------------------------------------------------------!
     ires = RDStr(',',istart,str,TABLEFN,1)

     call CalcTextFileLength(TABLEFN,NLIN)
     !----------------------------------------------------------------------!
     !     REALLOCATING ARRAYS.
     !----------------------------------------------------------------------!
     nk=NRAIN + NLIN
!     call Reallocate2DR8(NRAIN,3,nk,3,SRAIN)
     call ReallocateR8(NRAIN,nk,SPREC)
     call ReallocateR8(NRAIN,nk,SPREC_TM)
     !----------------------------------------------------------------------!
     !     Reading table.
     !----------------------------------------------------------------------!
     IFU=1001
     OPEN(UNIT=IFU, FILE=TABLEFN)
     call ReadTableCOASTOX(1,IFU,NLIN,NSPREC,NRAIN,KSPREC,'',SPREC_TM,1,SPREC,SPREC,SPREC,.true.)
     CLOSE(IFU)

   elseif (Type(1:7)=='tabular') then
     ITEMP = 2
   endif
   IPRCV(NSPREC) = ITEMP

   !----------------------------------------------------------------------!
   !     Read Range.
   !----------------------------------------------------------------------!
   call ReadRange(',',istart,str,ISPREC(:,NSPREC))

   !----------------------------------------------------------------------!
   !     Write Rain data to output file.
   !----------------------------------------------------------------------!
   if (iLNS/=1) write(IWR,*)
   write(IWR,'(A)') ' Rainfall Fallout'
   if (ITEMP==1) then
     write(IWR,'(A)') ' Temporal Variation: Constant'
     write(IWR,'(A,1PG12.5)') ' Rainfall Rate,       m/s: ', Value(1)
     write(IWR,'(A,1PG12.5)') ' Rainfall Start Time,   s: ', Value(2)
     write(IWR,'(A,1PG12.5)') ' Rainfall Stop Time,    s: ', Value(3)
   elseif (ITEMP==21) then
     write(IWR,'(2A)') ' Temporal Variation: Tabular from File - ', trim(TABLEFN)
   elseif (ITEMP==2) then
     write(IWR,'(A)') ' Temporal Variation: Tabular'
   endif
   !----------------------------------------------------------------------!
   !     Write range.
   !----------------------------------------------------------------------!
   call WriteRange(ISPREC(:,NSPREC))

   !----------------------------------------------------------------------!
   !     Read source tables.
   !----------------------------------------------------------------------!
   if (ITEMP/=2) GOTO 900

   !----------------------------------------------------------------------!
   !     Length of table.
   !----------------------------------------------------------------------!
   read(IRD,*) NLIN
   !----------------------------------------------------------------------!
   !     REALLOCATING ARRAYS.
   !----------------------------------------------------------------------!
   nk=NRAIN + NLIN
!   call Reallocate2DR8(NRAIN,3,nk,3,SRAIN)
   call ReallocateR8(NRAIN,nk,SPREC)
   call ReallocateR8(NRAIN,nk,SPREC_TM)
   !----------------------------------------------------------------------!
   !     Reading table.
   !----------------------------------------------------------------------!
   call ReadTableCOASTOX(1,IRD,NLIN,NSPREC,NRAIN,KSPREC,'',SPREC_TM,1,SPREC,SPREC,SPREC,.true.)


   900 CONTINUE

   if (ITEMP<=1) CYCLE
   !----------------------------------------------------------------------!
   !     Write source table header to output file.
   !----------------------------------------------------------------------!
   write(IWR,'(A)') ' Rainfall Table'
   !----------------------------------------------------------------------!
   !     Write source table to output file.
   !----------------------------------------------------------------------!
   do i = KSPREC(1,NSPREC),KSPREC(2,NSPREC)
     write(IWR,'(A,I4,4(A,1PG12.5))') ' Entry No.',i-KSPREC(1,NSPREC)+1,   &
               ';  Time, s: ',SPREC_TM(i), ';  Rainfall Rate, m/s: ',SPREC(i) 
   enddo
   !----------------------------------------------------------------------!
   !     Change tabularfile to tabular type for farther.
   !----------------------------------------------------------------------!
   if( IPRCV(NSPREC)==21 ) IPRCV(NSPREC)=2
!----------------------------------------------------------------------!
 enddo
!----------------------------------------------------------------------!
!----------------------------------------------------------------------C 
!     End of RDRAIN group. 
!----------------------------------------------------------------------C 
RETURN
END
!==============================================================================!

!==============================================================================!
SUBROUTINE RDEVAP
!----------------------------------------------------------------------C 
!
!     RDEVAP: ReaD EVAP input group.
!
!----------------------------------------------------------------------C 
USE DREVAP
USE FILINX 
USE intersub  
!----------------------------------------------------------------------C 
IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
      INCLUDE 'utils.inc' ! subroutine interfaces
!----------------------------------------------------------------------C 
CHARACTER*500 str,str1
CHARACTER*80 Units,Type
CHARACTER*150 TABLEFN
REAL*8 Value(3)
!----------------------------------------------------------------------C 
!----------------------------------------------------------------------!
!     Write header to output file.
!----------------------------------------------------------------------!
 write(IWR,'(/A)') ' Evapotranspiration Rate'
 write(IWR,'(A)')  ' -----------------------'
!----------------------------------------------------------------------!
!     Set table counters to zero.
!----------------------------------------------------------------------!
 NSEVAP = 0
 NEVAP  = 0

!----------------------------------------------------------------------!
!     Read the number of EVAP input lines.
!----------------------------------------------------------------------!
 read(IRD,*) nLNS

!----------------------------------------------------------------------!
!     ALLOCATING ARRAYS.
!----------------------------------------------------------------------!
 nk=nLNS
 ALLOCATE(IEVPV(nk), KSEVAP(2,nk), ISEVAP(4,nk)); IEVPV=0; KSEVAP=0; ISEVAP=0.d0;

!----------------------------------------------------------------------!
!     Loop over the number of lines.
!----------------------------------------------------------------------!
 do iLNS = 1,nLNS
   read(IRD,'(A)') str
   istart=1
   NSEVAP = NSEVAP + 1
   !----------------------------------------------------------------------!
   !     Read temporal variation.
   !----------------------------------------------------------------------!
   ITEMP = 0 !--Current temporal variation
   ires = RDStr(',',istart,str,Type)
   if (Type(1:8)=='constant') then
     ITEMP = 1
     !----------------------------------------------------------------------!
     !     REALLOCATING ARRAYS.
     !----------------------------------------------------------------------!
     nk=NEVAP + 2
     call ReallocateR8(NEVAP,nk,SEVAP)
     call ReallocateR8(NEVAP,nk,SEVAP_TM)
     !----------------------------------------------------------------------!
     !     Read Evap rate, starting and ending time.
     !----------------------------------------------------------------------!
     do i=1,3
       ires = RDR8(',',istart,str,Value(i))
       ires = RDStr(',',istart,str,Units)
       call ConvToSI(Value(i),Units)
     enddo
     !----------------------------------------------------------------------!
     !     Assign Evap variables.
     !----------------------------------------------------------------------!
     NEVAP = NEVAP + 1
     KSEVAP(1,NSEVAP) = NEVAP
     KSEVAP(2,NSEVAP) = NEVAP+1
     SEVAP(NEVAP)   = Value(1)
     SEVAP_TM(NEVAP)  = Value(2)
     NEVAP = NEVAP + 1
     SEVAP_TM(NEVAP)  = Value(3)

   elseif (Type(1:11) == 'tabularfile' ) then
     ITEMP = 21
     !----------------------------------------------------------------------!
     !     Read external file containing table of values.
     !----------------------------------------------------------------------!
     ires = RDStr(',',istart,str,TABLEFN,1)

     call CalcTextFileLength(TABLEFN,NLIN)
     !----------------------------------------------------------------------!
     !     REALLOCATING ARRAYS.
     !----------------------------------------------------------------------!
     nk=NEVAP + NLIN
     call ReallocateR8(NEVAP,nk,SEVAP)
     call ReallocateR8(NEVAP,nk,SEVAP_TM)
     !----------------------------------------------------------------------!
     !     Reading table.
     !----------------------------------------------------------------------!
     IFU=1001
     OPEN(UNIT=IFU, FILE=TABLEFN)
     call ReadTableCOASTOX(1,IFU,NLIN,NSEVAP,NEVAP,KSEVAP,'',SEVAP_TM,1,SEVAP,SEVAP,SEVAP,.true.)
     CLOSE(IFU)

   elseif (Type(1:7)=='tabular' ) then
     ITEMP = 2
   endif
   IEVPV(NSEVAP) = ITEMP

   !----------------------------------------------------------------------!
   !     Read Range.
   !----------------------------------------------------------------------!
   call ReadRange(',',istart,str,ISEVAP(:,NSEVAP))

   !----------------------------------------------------------------------!
   !     Write Evap data to output file.
   !----------------------------------------------------------------------!
   if (iLNS/=1) write(IWR,*)
   write(IWR,'(A)') ' Evapotranspiration'
   if (ITEMP==1) then
     write(IWR,'(A)') ' Temporal Variation: Constant'
     write(IWR,'(A,1PG12.5)') ' Evapotranspiration Rate,       m/s: ', Value(1)
     write(IWR,'(A,1PG12.5)') ' Evapotranspiration Start Time,   s: ', Value(2)
     write(IWR,'(A,1PG12.5)') ' Evapotranspiration Stop Time,    s: ', Value(3)
   elseif (ITEMP==21) then
     write(IWR,'(2A)') ' Temporal Variation: Tabular from File - ', trim(TABLEFN)
   elseif (ITEMP==2) then
     write(IWR,'(A)') ' Temporal Variation: Tabular'
   endif
   !----------------------------------------------------------------------!
   !     Write range.
   !----------------------------------------------------------------------!
   call WriteRange(ISEVAP(:,NSEVAP))

   !----------------------------------------------------------------------!
   !     Read source tables.
   !----------------------------------------------------------------------!
   if (ITEMP/=2) GOTO 900

   !----------------------------------------------------------------------!
   !     Length of table.
   !----------------------------------------------------------------------!
   read(IRD,*) NLIN
   !----------------------------------------------------------------------!
   !     REALLOCATING ARRAYS.
   !----------------------------------------------------------------------!
   nk=NEVAP + NLIN
   call ReallocateR8(NEVAP,nk,SEVAP)
   call ReallocateR8(NEVAP,nk,SEVAP_TM)
   !----------------------------------------------------------------------!
   !     Reading table.
   !----------------------------------------------------------------------!
   call ReadTableCOASTOX(1,IRD,NLIN,NSEVAP,NEVAP,KSEVAP,'',SEVAP_TM,1,SEVAP,SEVAP,SEVAP,.true.)


   900 CONTINUE

   if (ITEMP<=1) CYCLE
   !----------------------------------------------------------------------!
   !     Write source table to output file.
   !----------------------------------------------------------------------!
   write(IWR,'(A)') ' Evapotranspiration Table'
   do i = KSEVAP(1,NSEVAP),KSEVAP(2,NSEVAP)
     write(IWR,'(A,I4,2(A,1PG12.5))') ' Entry No.',i-KSEVAP(1,NSEVAP)+1,  &
       ';  Time, s: ',SEVAP_TM(i), ';  Evapotranspiration Rate, m/s: ',SEVAP(i)
   enddo
   !----------------------------------------------------------------------!
   !     Change tabularfile to tabular type for farther.
   !----------------------------------------------------------------------!
   if (IEVPV(NSEVAP)==21) IEVPV(NSEVAP)=2
!----------------------------------------------------------------------!
 enddo
!----------------------------------------------------------------------!
!----------------------------------------------------------------------C 
!     End of RDEVAP group. 
!----------------------------------------------------------------------C 
RETURN
END
!==============================================================================!


!==============================================================================!
SUBROUTINE RDWIND
!----------------------------------------------------------------------C 
!
!     RDWIND: ReaD WIND input group.
!
!----------------------------------------------------------------------C 
USE FILINX
USE WIND
USE intersub  
!----------------------------------------------------------------------C 
IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
      INCLUDE 'utils.inc' ! subroutine interfaces
!----------------------------------------------------------------------C 
CHARACTER*500 str,str1
CHARACTER*80 Units,Type
CHARACTER*150 TABLEFN
REAL*8 Value(4)
!----------------------------------------------------------------------C 
!----------------------------------------------------------------------!
!     Write header to output file.
!----------------------------------------------------------------------!
 write(IWR,'(/A)') ' Wind Speed'
 write(IWR,'(A)')  ' ----------'
!----------------------------------------------------------------------!
!     Set table counters to zero.
!----------------------------------------------------------------------!
 NSWIND = 0
 NWIND  = 0

 read(IRD,'(A)') str
 istart=1
!----------------------------------------------------------------------!
!     Read type of formula for wind drag coefficient.
!----------------------------------------------------------------------!
 ires = RDStr(',',istart,str,str1)
 if (str1(1:3)=='hsu') then !--Hsu
   iWDCType=1
   write(IWR,'(A)') ' Wind drag coefficient formula: Hsu, 1988'
 elseif (str1(1:4)=='garr') then !--Garratt
   iWDCType=2
   write(IWR,'(A)') ' Wind drag coefficient formula: Garrat, 1977'
 else !--Hsu
   iWDCType=1
   write(IWR,'(A)') ' Wind drag coefficient formula: Hsu, 1988'
 endif
!----------------------------------------------------------------------!
!     Read air density.
!----------------------------------------------------------------------!
 ires = RDR8(',',istart,str,RHAIR)
 ires = RDStr(',',istart,str,Units)
 call ConvToSI(RHAIR,Units)
 write(IWR,'(A,1PG12.5)') ' Air Density, kg/m^3: ', RHAIR
!----------------------------------------------------------------------!
!     Read heigth of anemometer.
!----------------------------------------------------------------------!
 ires = RDR8(',',istart,str,ANEMHEIGHT)
 ires = RDStr(',',istart,str,Units)
 call ConvToSI(ANEMHEIGHT,Units)
 if (ANEMHEIGHT<1.d0) ANEMHEIGHT=10.d0 !--If height is too small set it to default
 PreW10 = (10.d0/ANEMHEIGHT)**(1.d0/7.d0) !--Calc PreW10 value
 WindPreCoef = RHAIR/RHOLQ*PreW10
 write(IWR,'(A,1PG12.5)') ' Anemometer height, m: ', ANEMHEIGHT

!----------------------------------------------------------------------!
!     Read the number of WIND input lines.
!----------------------------------------------------------------------!
 read(IRD,*) nLNS
!----------------------------------------------------------------------!
!     ALLOCATING ARRAYS.
!----------------------------------------------------------------------!
 nk=nLNS
 ALLOCATE(IWINDV(nk), KSWIND(2,nk), WINDR(4,nk)); IWINDV=0; KSWIND=0; WINDR=0.d0;

!----------------------------------------------------------------------!
!     Loop over the number of lines.
!----------------------------------------------------------------------!
 do iLNS = 1,nLNS
   read(IRD,'(A)') str
   istart=1
   NSWIND = NSWIND + 1
   !----------------------------------------------------------------------!
   !     Read wind velocity temporal variation.
   !----------------------------------------------------------------------!
   ITEMP = 0 !--Current temporal variation
   ires = RDStr(',',istart,str,Type)
   if (Type(1:8)=='constant') then
     ITEMP = 1
     !----------------------------------------------------------------------!
     !     REALLOCATING ARRAYS.
     !----------------------------------------------------------------------!
     nk=NWIND + 2
!     call Reallocate2DR8(NWIND,2,nk,2,SWIND)
     call ReallocateR8(NWIND,nk,SWIND_TM)
     !----------------------------------------------------------------------!
     !     Read Wind Velocity,angle, starting and ending time.
     !----------------------------------------------------------------------!
     do i = 1, 4
       ires = RDR8(',',istart,str,Value(i))
       ires = RDStr(',',istart,str,Units)
       call ConvToSI(Value(i),Units)
     enddo
     !----------------------------------------------------------------------!
     !     Assign wind speed variables.
     !----------------------------------------------------------------------!
     NWIND = NWIND + 1
     KSWIND(1,NSWIND) = NWIND
     KSWIND(2,NSWIND) = NWIND+1

     SWIND(NWIND,1) = Value(1)*cos(Value(2))
     SWIND(NWIND,2) = Value(1)*sin(Value(2))
     SWIND_TM(NWIND) = Value(3)
     NWIND = NWIND + 1
     SWIND_TM(NWIND) = Value(4)

   elseif (Type(1:11) == 'tabularfile' ) then
     ITEMP = 21
     !----------------------------------------------------------------------!
     !     Read external file containing table of values.
     !----------------------------------------------------------------------!
     ires = RDStr(',',istart,str,TABLEFN,1)

     call CalcTextFileLength(TABLEFN,NLIN)
     !----------------------------------------------------------------------!
     !     REALLOCATING ARRAYS.
     !----------------------------------------------------------------------!
     nk=NWIND + NLIN
!     call Reallocate2DR8(NWIND,2,nk,2,SWIND)
     call ReallocateR8(NWIND,nk,SWIND_TM)
     !----------------------------------------------------------------------!
     !     Reading table.
     !----------------------------------------------------------------------!
     IFU=1001
     OPEN(UNIT=IFU, FILE=TABLEFN)
     call ReadTableCOASTOX(1,IFU,NLIN,NSWIND,NWIND,KSWIND,'',SWIND_TM,2,SWIND(:,1),SWIND(:,2),SWIND,.true.)
     CLOSE(IFU)

   elseif (Type(1:7)=='tabular' ) then
     ITEMP = 2
   endif
   IWINDV(NSWIND) = ITEMP

   !----------------------------------------------------------------------!
   !     Read Range.
   !----------------------------------------------------------------------!
   call ReadRange(',',istart,str,WINDR(:,NSWIND))

   !----------------------------------------------------------------------!
   !     Write wind data to output file.
   !----------------------------------------------------------------------!
   write(IWR,*)
   if (ITEMP==1) then
     write(IWR,'(A)') ' Temporal Variation: Constant'
     write(IWR,'(A,1PG12.5)') ' Wind Speed at anemometer height, m/s: ', Value(1)
     write(IWR,'(A,1PG12.5)') ' Wind Speed at 10m height, W10,   m/s: ', Value(1)*PreW10
     write(IWR,'(A,1PG12.5)') ' Wind Direction,  rad: ', Value(2)
     write(IWR,'(A,1PG12.5)') ' Wind Start Time,   s: ', Value(3)
     write(IWR,'(A,1PG12.5)') ' Wind Stop Time,    s: ', Value(4)
   elseif (ITEMP==21) then
     write(IWR,'(2A)') ' Temporal Variation: Tabular from File - ', trim(TABLEFN)
   elseif (ITEMP==2) then
     write(IWR,'(A)') ' Temporal Variation: Tabular'
   endif
   !----------------------------------------------------------------------!
   !     Write range.
   !----------------------------------------------------------------------!
   call WriteRange(WINDR(:,NSWIND))

   !----------------------------------------------------------------------!
   !     Read source tables.
   !----------------------------------------------------------------------!
   if (ITEMP/=2) GOTO 900

   !----------------------------------------------------------------------!
   !     Length of table.
   !----------------------------------------------------------------------!
   read(IRD,*) NLIN
   !----------------------------------------------------------------------!
   !     REALLOCATING ARRAYS.
   !----------------------------------------------------------------------!
   nk=NWIND+NLIN
!   call Reallocate2DR8(NWIND,2,nk,2,SWIND)
   call ReallocateR8(NWIND,nk,SWIND_TM)
   !----------------------------------------------------------------------!
   !     Reading table.
   !----------------------------------------------------------------------!
   call ReadTableCOASTOX(1,IRD,NLIN,NSWIND,NWIND,KSWIND,'',SWIND_TM,2,SWIND(:,1),SWIND(:,2),SWIND,.true.)


   900 CONTINUE

   if (ITEMP<=1) CYCLE
   !----------------------------------------------------------------------!
   !     Write source table header to output file.
   !----------------------------------------------------------------------!
   write(IWR,'(A)') ' Wind Speed Table '
   !----------------------------------------------------------------------!
   !     Write source table to output file.
   !----------------------------------------------------------------------!
   do i = KSWIND(1,NSWIND),KSWIND(2,NSWIND)
     wMagn=SWIND(i,1)
     wDir=SWIND(i,2)
     write(IWR,'(A,I4,5(A,1PG12.5))') ' Entry No.',i-KSWIND(1,NSWIND)+1, ';  Time, s: ',SWIND_TM(i), &
                 ';  Wind Speed, m/s: ',wMagn, ';  Wind Direction, rad: ',wDir, ';  W10 speed, m/s: ',wMagn*PreW10
    SWIND(i,1)   = wMagn*cos(wDir)
    SWIND(i,2)   = wMagn*sin(wDir)
   enddo
   !----------------------------------------------------------------------!
   !     Change tabularfile to tabular type for farther (WINDCALC, etc.).
   !----------------------------------------------------------------------!
   if (IWINDV(NSWIND)==21) IWINDV(NSWIND)=2

!----------------------------------------------------------------------!
 enddo
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
!     ALLOCATING ARRAYS.
!----------------------------------------------------------------------!
 if (NSWIND>0) fWind=.true. 
! 
!----------------------------------------------------------------------C 
!     End of RDWIND group. 
!----------------------------------------------------------------------C 
RETURN
END
!==============================================================================!

!==============================================================================!
SUBROUTINE ReadFrequency(str,istart,IFQType,IFQ,DTOUT,PrefixStr)
!----------------------------------------------------------------------C 
!
!     ReadFrequency: READ coastox output frequency and units.
!
!----------------------------------------------------------------------C 
USE FILINX
USE NUMBRS
!USE RW_MOD
use intersub
!----------------------------------------------------------------------C 
IMPLICIT NONE
!----------------------------------------------------------------------C 
CHARACTER*(*) str,PrefixStr
INTEGER*4 istart,IFQType,IFQ
REAL*8 DTOUT
!----------------------------------------------------------------------C 
CHARACTER*80 Units
INTEGER*4 ires
REAL*8 tt,rU
!----------------------------------------------------------------------C 
 ires = RDR8(',',istart,str,tt)
 ires = RDStr(',',istart,str,Units)
 if ((Units(1:8)=='timestep').or.(Units(1:4)=='null')) then !--Timestep based
   IFQType = 1
   IFQ = nint(tt)
   if (IFQ<1) then
     IFQ = IBIG
     write(IWR,'(A)') trim(PrefixStr)//' Output Frequency:  Never'
   else
     write(IWR,'(A,I0,A)') trim(PrefixStr)//' Output Frequency:  Every ',IFQ,' Time Step(s)'
   endif
 else !--Time based
   IFQType = 2
   rU = 1.d0; call ConvToSI(rU,Units)
   DTOUT = tt*rU
   if (DTOUT<SMALL) then
     DTOUT = BIG
     write(IWR,'(A)') trim(PrefixStr)//' Output Frequency:  Never'
   else
     write(IWR,'(A,1PG12.5,A)') trim(PrefixStr)//' Output Frequency:  Every ',DTOUT/rU,' '//trim(Units)
   endif
 endif
!
!----------------------------------------------------------------------C 
!     End of ReadFrequency group. 
!----------------------------------------------------------------------C 
RETURN
END
!==============================================================================!

!==============================================================================!
SUBROUTINE ReadTableCOASTOX(iMode,IFU,NTLines,NTAB,iLine,TI,TBNAME,TIME,NV,V1,V2,V3,fUnits)
!----------------------------------------------------------------------C 
!
!     ReadTableCOASTOX: READ universal TABLE from file.
!
!----------------------------------------------------------------------C 
!USE RW_MOD
USE intersub  
!----------------------------------------------------------------------C 
IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
!----------------------------------------------------------------------C 
INTEGER*4 iMode,IFU,NTLines,NTAB,iLine,NV
REAL*8 TIME(*),V1(*),V2(*),V3(*)
INTEGER*4 TI(2,*)
CHARACTER*(*) TBNAME(*)
LOGICAL*4 fUnits !--Read or not V1,V2,V3 units
!----------------------------------------------------------------------C 
CHARACTER*500 str,TABLE
CHARACTER*80 str1,Units
REAL*8 Value(6)
!----------------------------------------------------------------------C 
ITAB=NTAB
!----------------------------------------------------------------------!
!     Read a table ...
!----------------------------------------------------------------------!
 do i = 1,NTLines
   read(IFU,'(A)') str
   istart=1
   !----------------------------------------------------------------------!
   !     Read the table name.
   !----------------------------------------------------------------------!
   if (iMode==2) then
     ires = RDStr(',',istart,str,TABLE)
     !----------------------------------------------------------------------!
     !     Check for matching table name.
     !----------------------------------------------------------------------!
     do ITAB = 1, NTAB
       if (TABLE == TBNAME(ITAB)) GOTO 10
     enddo
     write(*,'(1A)') ' INPUT ERROR! - No matching table name!'
     write(*,'(2A)') ' INPUT ERROR! - Name we''re looking for is: ',TABLE
     STOP
     10 CONTINUE
   endif
   !----------------------------------------------------------------------!
   !     Read time with units.
   !----------------------------------------------------------------------!
   ires = RDR8(',',istart,str,Value(1))
   ires = RDStr(',',istart,str,Units)
   call ConvToSI(Value(1),Units)
   !----------------------------------------------------------------------!
   !     Read NV variables with or without units.
   !----------------------------------------------------------------------!
   do j = 2, NV+1
     ires = RDR8(',',istart,str,Value(j))
     if (fUnits) then
       ires = RDStr(',',istart,str,Units)
       call ConvToSI(Value(j),Units)
     endif
   enddo
   !----------------------------------------------------------------------!
   !     Assign source variables.
   !----------------------------------------------------------------------!
   iLine = iLine +1
   TIME(iLine) = Value(1)
   V1(iLine)   = Value(2)
   if (NV>=2) V2(iLine) = Value(3)
   if (NV>=3) V3(iLine) = Value(4)

   !----------------------------------------------------------------------!
   !     Table start\end indices.
   !----------------------------------------------------------------------!
   if (TI(1,ITAB) == 0) TI(1,ITAB) = iLine !--If first entry of this table then set the starting position
   TI(2,ITAB) = iLine
 enddo
!
!----------------------------------------------------------------------C 
!     End of ReadTableCOASTOX group. 
!----------------------------------------------------------------------C 
RETURN
END
!==============================================================================! 

!==============================================================================! 
SUBROUTINE CalcTextFileLength(FName,NN)
!----------------------------------------------------------------------C 
USE RW_MOD
!----------------------------------------------------------------------C 
IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
!----------------------------------------------------------------------C 
CHARACTER*(*) FName
INTEGER*4 NN
REAL*8 x,y
CHARACTER*1000 str
!----------------------------------------------------------------------C 
 write(*,'(3A$)') 'Calculating number of rows in text file:  ', trim(FName),' .....'

 OPEN(UNIT=RW_IFUTemp, FILE=FName, STATUS='OLD')
  NN=0
  ne=0
  do while(.true.)
    read(RW_IFUTemp,'(A)',end=10,err=10) str
    ln=len_trim(str)
    if (ln==0) then !--Don't count last blank lines only (inner are counted)
      ne=ne+1
    else
      ne=0
    endif
    NN=NN+1
  enddo
  10 CONTINUE
  NN=NN-ne
 CLOSE(RW_IFUTemp)

 write(*,'(A,I0,A)') 'Done (',NN,' rows)'
!
!----------------------------------------------------------------------C 
!     End of CalcTextFileLength group. 
!----------------------------------------------------------------------C 
RETURN
END SUBROUTINE
!==============================================================================! 
    
