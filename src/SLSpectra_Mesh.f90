MODULE SLSpectra_Mesh

USE SLSpectra_Precision
 IMPLICIT NONE

   TYPE Mesh
      INTEGER                 :: nX, nY, nDOF

      ! Derived quantities
      REAL(prec), ALLOCATABLE :: tracermask(:,:)
      REAL(prec), ALLOCATABLE :: x(:,:)
      REAL(prec), ALLOCATABLE :: y(:,:)
      REAL(prec) :: dx, dy
      INTEGER, ALLOCATABLE    :: DOFtoIJ(:,:), IJtoDOF(:,:)

       
       CONTAINS

          PROCEDURE :: Build => Build_Mesh
          PROCEDURE :: Free => Free_Mesh
    
          PROCEDURE :: ConstructDirichletCube
          PROCEDURE :: ConstructWetPointMap
          PROCEDURE :: FlatMap
          PROCEDURE :: GridMap

    END TYPE Mesh

    !TYPE, EXTENDS(Mesh) :: MITgcmMesh
    !  INTEGER :: filePrecision
    !  REAL(prec), ALLOCATABLE :: 
    !END TYPE MITgcmMesh
    
    REAL(prec), PARAMETER :: SLSpectra_wetValue = 1.0_prec
    REAL(prec), PARAMETER :: SLSpectra_dryValue = 0.0_prec
    REAL(prec), PARAMETER :: SLSpectra_fillValue = -99999.0_prec


 CONTAINS
!
!==================================================================================================!
!----------------------------- Manual Constructor/Destructor --------------------------------------!
!==================================================================================================!
!
 SUBROUTINE Build_Mesh( this, nX, nY )
 ! S/R Build
 !  
 !    
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Mesh), INTENT(out) :: this
   INTEGER, INTENT(in)      :: nX, nY

      this % nX = nX
      this % nY = nY
      ! Tracer mesh
      ALLOCATE( this % tracermask(1:nX,1:nY), &
                this % IJtoDOF(1:nX,1:nY), &
                this % x(1:nX,1:nY), &
                this % y(1:nX,1:nY))

      
 END SUBROUTINE Build_Mesh
!
 SUBROUTINE Free_Mesh( this )
 ! S/R Free
 !  
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Mesh), INTENT(inout) :: this

      DEALLOCATE( this % tracermask )
      DEALLOCATE( this % x )
      DEALLOCATE( this % y )

      IF( ALLOCATED( this % DOFtoIJ ) )THEN
         DEALLOCATE( this % DOFtoIJ )
      ENDIF
      IF( ALLOCATED( this % IJtoDOF ) )THEN
         DEALLOCATE( this % IJtoDOF )
      ENDIF

 END SUBROUTINE Free_Mesh
!
 SUBROUTINE ConstructDirichletCube( this )
 !! Constructs a cube domain on [0,1]^2 with the tracermask
 !! set to 0 on the south, east, north, and west boundaries.
 !!
 !! This mesh construction is meant to be used for testing purposes
 !!
   IMPLICIT NONE
   CLASS(Mesh), INTENT(inout) :: this
   ! Local
   INTEGER    :: i,j
   REAL(prec) :: dx, dy

     dx = 1.0_prec/REAL(this % nX,prec)
     dy = 1.0_prec/REAL(this % nY,prec)
     this % dx = dx
     this % dy = dy

     DO j = 1, this % nY
       DO i = 1, this % nX
         this % x(i,j) = (REAL(i,prec) + 0.5_prec)*dx
         this % y(i,j) = (REAL(j,prec) + 0.5_prec)*dy
       ENDDO
     ENDDO

     this % tracerMask = SLSpectra_wetValue
     ! South and North Boundaries
     DO i = 1, this % nX
       this % tracermask(i,1) = SLSpectra_dryValue
       this % tracermask(i,this % nY) = SLSpectra_dryValue
     ENDDO

     ! East and West Boundaries
     DO j = 1, this % nY
       this % tracermask(1,j) = SLSpectra_dryValue
       this % tracermask(this % nX,j) = SLSpectra_dryValue
     ENDDO

     CALL this % ConstructWetPointMap()

 END SUBROUTINE ConstructDirichletCube

 SUBROUTINE ConstructWetPointMap( this )
   IMPLICIT NONE
   CLASS( Mesh ), INTENT(inout) :: this
   ! LOCAL
   INTEGER :: i, j, nDOF


     ! Count the number of wet-points and set the IJK to DOF mapping
     nDOF = 0
     DO j = 1, this % nY
       DO i = 1, this % nX
         IF( this % tracerMask(i,j) == SLSpectra_wetValue )THEN
           nDOF = nDOF + 1
           this % IJtoDOF(i,j) = nDOF
         ENDIF
       ENDDO
     ENDDO

     PRINT*, "Number of DOF : ", nDOF

     ALLOCATE( this % DOFtoIJ(1:2,1:nDOF) ) 
     this % nDOF = nDOF

     ! Now we can set the DOF to IJ mapping
     nDOF = 0
     DO j = 1, this % nY
       DO i = 1, this % nX
         IF( this % tracerMask(i,j) == SLSpectra_wetValue )THEN
           nDOF = nDOF + 1
           this % DOFtoIJ(1,nDOF) = i
           this % DOFtoIJ(2,nDOF) = j
         ENDIF
       ENDDO
     ENDDO

 END SUBROUTINE ConstructWetPointMap

 FUNCTION GridMap( this, dofArray ) RESULT( ijArray )
   IMPLICIT NONE
   CLASS( Mesh ) :: this
   REAL(prec)    :: dofArray(1:this % nDOF)
   REAL(prec)    :: ijArray(1:this % nX, 1:this % nY)
   ! Local
   INTEGER :: i, j, l


      ijArray = 0.0_prec
      DO l = 1, this % nDof

         i = this % DOFtoIJ(1,l)
         j = this % DOFtoIJ(2,l)

         ijArray(i,j) = dofArray(l)

      ENDDO

 END FUNCTION GridMap

 FUNCTION FlatMap( this, ijArray ) RESULT( dofArray )
   IMPLICIT NONE
   CLASS( Mesh ) :: this
   REAL(prec)    :: ijArray(1:this % nX, 1:this % nY)
   REAL(prec)    :: dofArray(1:this % nDOF)
   ! Local
   INTEGER :: i, j, l

      DO l = 1, this % nDof

         i = this % DOFtoIJ(1,l)
         j = this % DOFtoIJ(2,l)

         dofArray(l) = ijArray(i,j)

      ENDDO

 END FUNCTION FlatMap

END MODULE SLSpectra_Mesh
