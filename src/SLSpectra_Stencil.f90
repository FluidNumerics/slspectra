!
! Adapted from https://github.com/lanl/feots (POP_Stencil_Class)
!
! ////////////////////////////////////////////////////////////////////////////////////////////////

 
MODULE SLSpectra_Stencil


IMPLICIT NONE


   TYPE, ABSTRACT :: Stencil
     !! The Stencil class is an abstract class that lays out the structure for defining a 
     !! discrete computational stencil. The Build method for this abstract class is deferred
     !! intentionally. The Free method is provided so that you don't have to add your own
     !! Free method for a type extended class.
     !!
     !! Adding a specific stencil requires that you declare a type extension of
     !! the Stencil class and then concretize the Build type bound procedure.
     !! When using a specific stencil, you simply have to declare your object as one of the
     !! concretized types.
     !!
     INTEGER              :: nPoints ! Number of points in the stencil
     INTEGER, ALLOCATABLE :: neighbors(:,:) ! indices of the neighbors relative to the center

     CONTAINS 
  
     PROCEDURE(Build_Stencil), DEFERRED :: Build
     PROCEDURE :: Free => Free_Stencil

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

 SUBROUTINE Free_Stencil( this )

   IMPLICIT NONE
   CLASS( Stencil ),INTENT(inout) :: this

   DEALLOCATE( this % neighbors )
   this % nPoints = -1

 END SUBROUTINE Free_Stencil


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

END MODULE SLSpectra_Stencil
