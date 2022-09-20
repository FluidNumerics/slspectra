MODULE SLSpectra_MITgcmMesh

USE SLSpectra_Precision
USE SLSpectra_SupportRoutines
USE SLSpectra_Mesh


  TYPE, EXTENDS(Mesh) :: MITgcmMesh
    INTEGER :: nZ
    INTEGER :: kLev
    INTEGER :: filePrecision
    REAL(prec), ALLOCATABLE :: dxg(:,:)
    REAL(prec), ALLOCATABLE :: dyg(:,:)
    REAL(prec), ALLOCATABLE :: rac(:,:)
    REAL(prec), ALLOCATABLE :: hFacC(:,:,:)
    REAL(prec), ALLOCATABLE :: drf(:)
 
    CONTAINS

      GENERIC :: Build => Build_MITgcmMesh
      PROCEDURE, PRIVATE :: Build_MITgcmMesh
      PROCEDURE :: Free => Free_MITgcmMesh

      PROCEDURE :: Load => Load_MITgcmMesh
      PROCEDURE :: GenerateMaskFromhFac

  END TYPE MITgcmMesh

  CONTAINS
  
  SUBROUTINE Build_MITgcmMesh( this, nX, nY, nZ )

   IMPLICIT NONE
   CLASS(MITgcmMesh), INTENT(out) :: this
   INTEGER, INTENT(in)      :: nX, nY, nZ

      this % nX = nX
      this % nY = nY
      this % nZ = nZ
      
      ! Tracer mesh
      ALLOCATE( this % tracermask(1:nX,1:nY), &
                this % IJtoDOF(1:nX,1:nY), &
                this % x(1:nX,1:nY), &
                this % y(1:nX,1:nY), &
                this % dxg(1:nX,1:nY), &
                this % dyg(1:nX,1:nY), &
                this % rac(1:nX,1:nY), &
                this % hFacC(1:nX,1:nY,1:nZ), &
                this % drf(1:nZ))
      
 END SUBROUTINE Build_MITgcmMesh

 SUBROUTINE Free_MITgcmMesh( this )

   IMPLICIT NONE
   CLASS(MITgcmMesh), INTENT(inout) :: this

      DEALLOCATE( this % tracermask )
      DEALLOCATE( this % x )
      DEALLOCATE( this % y )
      DEALLOCATE( this % dxg )
      DEALLOCATE( this % dyg )
      DEALLOCATE( this % rac )
      DEALLOCATE( this % hFacC )
      DEALLOCATE( this % drf )

      IF( ALLOCATED( this % DOFtoIJ ) )THEN
         DEALLOCATE( this % DOFtoIJ )
      ENDIF
      IF( ALLOCATED( this % IJtoDOF ) )THEN
         DEALLOCATE( this % IJtoDOF )
      ENDIF

 END SUBROUTINE Free_MITgcmMesh
 
 SUBROUTINE Load_MITgcmMesh( this, mdpath )
 !!
 !! Loads the MITgcm tracer mesh. 
 !!
 !! Input:
 !! 
 !!   mdpath - the full path to a directory containing MITgcm mesh files
 !!
 !! The following files must be found in the mdpath
 !!
 !!   XC.data 
 !!   XC.meta
 !!   YC.data
 !!   YC.meta
 !!   DXG.data
 !!   DXG.meta
 !!   DYG.data
 !!   DYG.meta
 !!   RAC.data
 !!   RAC.meta
 !!   hFacC.data
 !!   hFacC.meta
 !!
   IMPLICIT NONE
   CLASS(MITgcmMesh), INTENT(out) :: this
   CHARACTER(*), INTENT(in) :: mdpath
   ! Local
   INTEGER :: fUnit
   INTEGER :: nDims
   INTEGER :: io
!   INTEGER :: ioprec
   INTEGER :: i, m, n
   INTEGER :: dims(1:3)
   CHARACTER(50) :: line
  
   ! Load meta file for nx, ny
    OPEN( UNIT=NEWUNIT(fUnit), &
      FILE= TRIM(mdpath)//'/hFacC.meta', &
      ACTION= 'READ', &
      FORM='formatted', &
      STATUS='old')
    
    READ(fUnit,'(A)') line
    IF( line(1:6) == ' nDims' )THEN
      READ(line(12:15),'(I4)') nDims
    ENDIF
    
    DO
      READ(fUnit,'(A)',iostat=io) line
      
      IF( io == 0 )THEN
      
        IF(  line(1:8) == ' dimList' )THEN
          DO i = 1, nDims
            READ(fUnit, *, iostat=io) dims(i), m, n
          ENDDO
          EXIT
        !ELSEIF(  line(1:9) == ' dataprec' )THEN

        !  IF( line(14:24) == ' float32 ' )THEN
        !    ioprec = 4
        !  ELSE
        !    ioprec = 8
        !  ENDIF
          
        ENDIF
        
      ELSE
      
        EXIT
        
      ENDIF
      
    ENDDO
    CLOSE(fUnit)
    
    CALL this % Build( dims(1), dims(2), dims(3) )
    
    ! Read longitude
    OPEN( UNIT=NEWUNIT(fUnit), &
      FILE= TRIM(mdpath)//'/XC.data', &
      ACTION= 'READ', &
      FORM='UNFORMATTED', &
      STATUS='OLD', &
      ACCESS='DIRECT', &
      CONVERT='BIG_ENDIAN', &
      RECL=prec*this % nX*this % nY)
    READ(fUnit,REC=1) this % x
    CLOSE(fUnit)
    
    ! Read latitude
    OPEN( UNIT=NEWUNIT(fUnit), &
      FILE= TRIM(mdpath)//'/YC.data', &
      ACTION= 'READ', &
      FORM='UNFORMATTED', &
      STATUS='OLD', &
      ACCESS='DIRECT', &
      CONVERT='BIG_ENDIAN', &
      RECL=prec*this % nX*this % nY)
    READ(fUnit,REC=1) this % y
    CLOSE(fUnit)
    
    ! Read x grid spacing
    OPEN( UNIT=NEWUNIT(fUnit), &
      FILE= TRIM(mdpath)//'/DXG.data', &
      ACTION= 'READ', &
      FORM='UNFORMATTED', &
      STATUS='OLD', &
      ACCESS='DIRECT', &
      CONVERT='BIG_ENDIAN', &
      RECL=prec*this % nX*this % nY)
    READ(fUnit,REC=1) this % dxg
    CLOSE(fUnit)
    
    ! Read y grid spacing
    OPEN( UNIT=NEWUNIT(fUnit), &
      FILE= TRIM(mdpath)//'/DYG.data', &
      ACTION= 'READ', &
      FORM='UNFORMATTED', &
      STATUS='OLD', &
      ACCESS='DIRECT', &
      CONVERT='BIG_ENDIAN', &
      RECL=prec*this % nX*this % nY)
    READ(fUnit,REC=1) this % dyg
    CLOSE(fUnit)
    
    ! Read hFacC
    OPEN( UNIT=NEWUNIT(fUnit), &
      FILE= TRIM(mdpath)//'/hFacC.data', &
      ACTION= 'READ', &
      FORM='UNFORMATTED', &
      STATUS='OLD', &
      ACCESS='DIRECT', &
      CONVERT='BIG_ENDIAN', &
      RECL=prec*this % nX*this % nY*this % nZ)
    READ(fUnit,REC=1) this % hFacC
    CLOSE(fUnit)
    
    ! Read cell area
    OPEN( UNIT=NEWUNIT(fUnit), &
      FILE= TRIM(mdpath)//'/RAC.data', &
      ACTION= 'READ', &
      FORM='UNFORMATTED', &
      STATUS='OLD', &
      ACCESS='DIRECT', &
      CONVERT='BIG_ENDIAN', &
      RECL=prec*this % nX*this % nY)
    READ(fUnit,REC=1) this % rac
    CLOSE(fUnit)
    
    ! Set the scalar dx and dy (for testing purposes only)
    this % dx = this % dxg(1,1)
    this % dy = this % dyg(1,1)
    
   
  END SUBROUTINE Load_MITgcmMesh
  
  SUBROUTINE GenerateMaskFromhFac( this, klev )
    IMPLICIT NONE
    CLASS(MITgcmMesh),  INTENT(inout) :: this
    INTEGER, INTENT(in) :: klev
    ! Local
    INTEGER :: i, j
    
      this % tracerMask = SLSpectra_wetValue
      DO j = 1, this % nY
        DO i = 1, this % nX
        
          IF( this % hFacC(i,j,klev) == 0.0_prec ) THEN
            this % tracerMask(i,j) = SLSpectra_dryValue
          ENDIF
          
        ENDDO
      ENDDO
      
      ! Mask boundary cells for homogeneous Dirichlet BC
      this % tracerMask(:,1) = SLSpectra_dryValue
      this % tracerMask(:,this % nY) = SLSpectra_dryValue
      this % tracerMask(1,:) = SLSpectra_dryValue
      this % tracerMask(this % nX,:) = SLSpectra_dryValue
  
  END SUBROUTINE GenerateMaskFromhFac
  
  
   
END MODULE SLSpectra_MITgcmMesh
