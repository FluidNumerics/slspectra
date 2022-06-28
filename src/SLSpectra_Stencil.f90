!
! Adapted from https://github.com/lanl/feots (POP_Stencil_Class)
!
! ////////////////////////////////////////////////////////////////////////////////////////////////

 
MODULE SLSpectra_Stencil


IMPLICIT NONE


   TYPE:: Stencil
     !! The Stencil class is a class that lays out the structure for defining a 
     !! discrete computational stencil. The Build method for this class is provided for
     !! memory allocation only. The Free method is provided so that you don't have to add your own
     !! Free method for a type extended class.
     !!
     INTEGER              :: nPoints ! Number of points in the stencil
     INTEGER, ALLOCATABLE :: neighbors(:,:) ! indices of the neighbors relative to the center

     CONTAINS 
  
     PROCEDURE :: Build => Build_Stencil
     PROCEDURE :: Free => Free_Stencil

   END TYPE Stencil

   TYPE :: Laplacian5Stencil
     !! Class for a 5 point laplacian overlap stencil in 2-D
     !! This stencil is equivalent to the stencil obtained by
     !! applying the 5-point laplacian operator twice on
     !! an impulse field
    
       TYPE(Stencil) :: firstOrder
       TYPE(Stencil) :: secondOrder
       
     CONTAINS

     PROCEDURE :: Build => Build_Laplacian5Stencil
     PROCEDURE :: Free => Free_Laplacian5Stencil

   END TYPE Laplacian5Stencil



CONTAINS

 SUBROUTINE Build_Stencil( this, nPoints ) 
 !! Constructor 
   IMPLICIT NONE
   CLASS( Stencil ), INTENT(out) :: this
   INTEGER,  INTENT(in) :: nPoints


     this % nPoints = nPoints
     
     ALLOCATE( this % neighbors(1:3, nPoints) )
     
 END SUBROUTINE Build_Stencil

 SUBROUTINE Free_Stencil( this )

   IMPLICIT NONE
   CLASS( Stencil ),INTENT(inout) :: this

   DEALLOCATE( this % neighbors )
   this % nPoints = -1

 END SUBROUTINE Free_Stencil


 SUBROUTINE Build_Laplacian5Stencil( this ) 
 !! Constructor for the Laplacian5Stencil
   IMPLICIT NONE
   CLASS( Laplacian5Stencil ), INTENT(out) :: this
   ! Local
   INTEGER :: i, j, nid

     ! Allocate space for the 5-point Laplacian stencil
     CALL this % firstOrder % Build( 5 )
     
     ! Allocate space for the overlap stencil associated with
     ! the 5 point Laplacian stencil
     CALL this % secondOrder % Build( 13 )
     
     ! Create the 5-point Laplacian stencil
     nid = 1
     this % firstOrder % neighbors(1:3,nid) = (/0,-1,0/)
     
     j = 0
     DO i = -1, 1
       nid = nid + 1
       this % firstOrder % neighbors(1:3,nid) = (/ i,j,0 /)
     ENDDO
     
     nid = nid + 1
     this % firstOrder % neighbors(1:3,nid) = (/0,1,0/)

     ! Create the overlap stencil (equivalent to applying the Laplacian twice) 
     nid = 1
     this % secondOrder % neighbors(1:3,nid) = (/0,-2,0/)

     j = -1
     DO i = -1, 1
       nid = nid + 1
       this % secondOrder % neighbors(1:3,nid) = (/ i,j,0 /)
     ENDDO

     j = 0
     DO i = -2, 2
       nid = nid + 1
       this % secondOrder % neighbors(1:3,nid) = (/ i,j,0 /)
     ENDDO

     j = 1
     DO i = -1, 1
       nid = nid + 1
       this % secondOrder % neighbors(1:3,nid) = (/ i,j,0 /)
     ENDDO

     nid = nid + 1
     this % secondOrder % neighbors(1:3,nid) = (/0,2,0/)


 END SUBROUTINE Build_Laplacian5Stencil
 
 SUBROUTINE Free_Laplacian5Stencil( this ) 
 !! Constructor for the Laplacian5Stencil
   IMPLICIT NONE
   CLASS( Laplacian5Stencil ), INTENT(inout) :: this
   
     CALL this % firstOrder % Free()
     CALL this % secondOrder % Free()
   
 END SUBROUTINE Free_Laplacian5Stencil

END MODULE SLSpectra_Stencil
