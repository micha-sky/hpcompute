      Interface   
        SUBROUTINE ReallocateR8 (N0,N1,V)
          INTEGER*4 N0,N1
          REAL*8, POINTER, intent(inout) :: V(:) 
        END SUBROUTINE

        SUBROUTINE ReallocateI4(N0,N1,V)
          INTEGER*4 N0,N1
          INTEGER*4, POINTER :: V(:)
        END SUBROUTINE

        SUBROUTINE Reallocate2DI4(NX0,NY0,NX1,NY1,V)
          INTEGER*4 NX0,NY0,NX1,NY1
          INTEGER*4, POINTER :: V(:,:)
        END SUBROUTINE

        SUBROUTINE Reallocate2DR8(NX0,NY0,NX1,NY1,V)
          INTEGER*4 NX0,NY0,NX1,NY1
          REAL*8, POINTER :: V(:,:)
        END SUBROUTINE

        SUBROUTINE Reallocate3DI4(NX0,NY0,NZ0,NX1,NY1,NZ1,V)
          INTEGER*4 NX0,NY0,NZ0,NX1,NY1,NZ1
          INTEGER*4, POINTER :: V(:,:,:)
        END SUBROUTINE
      end Interface
