MODULE SLSpectra_Mesh

USE SLSpectra_Precision
 IMPLICIT NONE

   TYPE :: Mesh
      INTEGER                 :: nX, nY, nZ, nDOF

      ! Derived quantities
      REAL(prec), ALLOCATABLE :: tracermask(:,:,:)
      REAL(prec), ALLOCATABLE :: x(:,:)
      REAL(prec), ALLOCATABLE :: y(:,:)
      REAL(prec), ALLOCATABLE :: z(:)
      INTEGER, ALLOCATABLE    :: DOFtoIJK(:,:), IJKtoDOF(:,:,:)

       
       CONTAINS

          PROCEDURE :: Build => Build_Mesh
          PROCEDURE :: Free => Free_Mesh
    
          PROCEDURE :: ConstructDirichletCube
          PROCEDURE :: ConstructWetPointMap
          PROCEDURE :: FlatMap
          PROCEDURE :: GridMap

    END TYPE Mesh

    REAL(prec), PARAMETER, PRIVATE :: wetValue = 1.0_prec
    REAL(prec), PARAMETER, PRIVATE :: dryValue = 0.0_prec
    REAL(prec), PARAMETER, PRIVATE :: fillValue = -99999.0_prec


 CONTAINS
!
!==================================================================================================!
!----------------------------- Manual Constructor/Destructor --------------------------------------!
!==================================================================================================!
!
 SUBROUTINE Build_Mesh( this, nX, nY, nZ )
 ! S/R Build
 !  
 !    
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(Mesh), INTENT(out) :: this
   INTEGER, INTENT(in)      :: nX, nY, nZ

      this % nX = nX
      this % nY = nY
      this % nZ = nZ 
      ! Tracer mesh
      ALLOCATE( this % tracermask(1:nX,1:nY,1:nZ), &
                this % IJKtoDOF(1:nX,1:nY,1:nZ), &
                this % x(1:nX,1:nY), &
                this % y(1:nX,1:nY), &
                this % z(1:nZ))

      
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
      DEALLOCATE( this % z )

      IF( ALLOCATED( this % DOFtoIJK ) )THEN
         DEALLOCATE( this % DOFtoIJK )
      ENDIF
      IF( ALLOCATED( this % IJKtoDOF ) )THEN
         DEALLOCATE( this % IJKtoDOF )
      ENDIF

 END SUBROUTINE Free_Mesh
!
 SUBROUTINE ConstructDirichletCube( this )
 !! Constructs a cube domain on [0,1]^3 with the tracermask
 !! set to 0 on the south, east, north, and west boundaries.
 !!
 !! This mesh construction is meant to be used for testing purposes
 !!
   IMPLICIT NONE
   CLASS(Mesh), INTENT(inout) :: this
   ! Local
   INTEGER    :: i,j,k
   REAL(prec) :: dx, dy, dz

     dx = 1.0_prec/REAL(this % nX,prec)
     dy = 1.0_prec/REAL(this % nY,prec)
     dz = 1.0_prec/REAL(this % nZ,prec)

     DO j = 1, this % nY
       DO i = 1, this % nX
         this % x(i,j) = (REAL(i,prec) + 0.5_prec)*dx
         this % y(i,j) = (REAL(j,prec) + 0.5_prec)*dy
       ENDDO
     ENDDO

     DO k = 1, this % nZ
       this % z(k) = (REAL(k,prec) + 0.5_prec)*dz
     ENDDO

     this % tracerMask = wetValue
     ! South and North Boundaries
     DO k = 1, this % nZ
       DO i = 1, this % nX
         this % tracermask(i,1,k) = dryValue
         this % tracermask(i,this % nY,k) = dryValue
       ENDDO
     ENDDO

     ! East and West Boundaries
     DO k = 1, this % nZ
       DO j = 1, this % nY
         this % tracermask(1,j,k) = dryValue
         this % tracermask(this % nX,j,k) = dryValue
       ENDDO
     ENDDO

     CALL this % ConstructWetPointMap()

 END SUBROUTINE ConstructDirichletCube

 SUBROUTINE ConstructWetPointMap( this )
   IMPLICIT NONE
   CLASS( Mesh ), INTENT(inout) :: this
   ! LOCAL
   INTEGER :: i, j, k, nDOF


     ! Count the number of wet-points and set the IJK to DOF mapping
     nDOF = 0
     DO k = 1, this % nZ
       DO j = 1, this % nY
         DO i = 1, this % nX
           IF( this % tracerMask(i,j,k) == wetValue )THEN
              nDOF = nDOF + 1
              this % IJKtoDOF(i,j,k) = nDOF
           ENDIF
         ENDDO
       ENDDO
     ENDDO

     PRINT*, "Number of DOF : ", nDOF

     ALLOCATE( this % DOFtoIJK(1:3,1:nDOF) ) 
     this % nDOF = nDOF

     ! Now we can set the DOF to IJK mapping
     nDOF = 0
     DO j = 1, this % nY
        DO i = 1, this % nX
           DO k = 1, this % nZ
              IF( this % tracerMask(i,j,k) == wetValue )THEN
                 nDOF = nDOF + 1
                 this % DOFtoIJK(1,nDOF) = i
                 this % DOFtoIJK(2,nDOF) = j
                 this % DOFtoIJK(3,nDOF) = k
              ENDIF
           ENDDO
        ENDDO
     ENDDO

 END SUBROUTINE ConstructWetPointMap

 FUNCTION GridMap( this, dofArray ) RESULT( ijkArray )
   IMPLICIT NONE
   CLASS( Mesh ) :: this
   REAL(prec)    :: dofArray(1:this % nDOF)
   REAL(prec)    :: ijkArray(1:this % nX, 1:this % nY, 1:this % nZ)
   ! Local
   INTEGER :: i, j, k, l


      ijkArray = fillValue
      DO l = 1, this % nDof

         i = this % DOFtoIJK(1,l)
         j = this % DOFtoIJK(2,l)
         k = this % DOFtoIJK(3,l)

         ijkArray(i,j,k) = dofArray(l)

      ENDDO

 END FUNCTION GridMap

 FUNCTION FlatMap( this, ijkArray ) RESULT( dofArray )
   IMPLICIT NONE
   CLASS( Mesh ) :: this
   REAL(prec)    :: ijkArray(1:this % nX, 1:this % nY, 1:this % nZ)
   REAL(prec)    :: dofArray(1:this % nDOF)
   ! Local
   INTEGER :: i, j, k, l

      DO l = 1, this % nDof

         i = this % DOFtoIJK(1,l)
         j = this % DOFtoIJK(2,l)
         k = this % DOFtoIJK(3,l)

         dofArray(l) = ijkArray(i,j,k)

      ENDDO

 END FUNCTION FlatMap

END MODULE SLSpectra_Mesh
