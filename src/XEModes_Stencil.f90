!
! Adapted from https://github.com/lanl/feots (POP_Stencil_Class)
!
! ////////////////////////////////////////////////////////////////////////////////////////////////

 
MODULE XEModes_Stencil


IMPLICIT NONE

   TYPE, ABSTRACT :: Stencil
      INTEGER              :: nPoints ! Number of points in the stencil
      INTEGER, ALLOCATABLE :: neighbors(:,:) ! indices of the neighbors relative to the center

      CONTAINS 
  
      PROCEDURE(Build_Stencil), DEFERRED :: Build
      PROCEDURE :: Trash => Trash_Stencil

   END TYPE Stencil

   TYPE, EXTENDS(Stencil) :: Laplacian5OStencil
     !! Class for a 5 point laplacian overlap stencil in 2-D
     !! This stencil is equivalent to the stencil obtained by
     !! applying the 5-point laplacian operator twice on
     !! an impulse field

     CONTAINS

     PROCEDURE :: Build => Build_Laplacian5OStencil

   END TYPE Laplacian5OStencil


   INTERFACE 
    SUBROUTINE Build_Stencil( this )
      IMPORT Stencil
      IMPLICIT NONE
      CLASS(Stencil),INTENT(out) :: this
    END SUBROUTINE Build_Stencil
   END INTERFACE


CONTAINS

 SUBROUTINE Trash_Stencil( this )

   IMPLICIT NONE
   CLASS( Stencil ),INTENT(inout) :: this

   DEALLOCATE( this % neighbors )
   this % nPoints = -1

 END SUBROUTINE Trash_Stencil


 SUBROUTINE Build_Laplacian5OStencil( this ) 
 !! Constructor for the Laplacian5OStencil
   IMPLICIT NONE
   CLASS( Laplacian5OStencil ), INTENT(out) :: this
   ! Local
   INTEGER :: i, j, nid

     this % nPoints = 13 
     ALLOCATE( this % neighbors(1:3,13) )

     nid = 1
     this % neighbors(1:3,nid) = (/0,-2,0/)

     j = -1
     DO i = -1, 1
       nid = nid + 1
       this % neighbors(1:3,nid) = (/ i,j,0 /)
     ENDDO

     j = 0
     DO i = -2, 2
       nid = nid + 1
       this % neighbors(1:3,nid) = (/ i,j,0 /)
     ENDDO

     j = 1
     DO i = -1, 1
       nid = nid + 1
       this % neighbors(1:3,nid) = (/ i,j,0 /)
     ENDDO

     nid = nid + 1
     this % neighbors(1:3,nid) = (/0,2,0/)


 END SUBROUTINE Build_Laplacian5OStencil

END MODULE XEModes_Stencil
