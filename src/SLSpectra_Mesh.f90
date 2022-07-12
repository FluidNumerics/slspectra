MODULE SLSpectra_Mesh

USE SLSpectra_Precision
 IMPLICIT NONE

   TYPE Mesh
      INTEGER                 :: nX, nY, nZ, nDOF
      INTEGER :: kLev
      ! Derived quantities
      REAL(prec), ALLOCATABLE :: tracermask(:,:)
      REAL(prec), ALLOCATABLE :: x(:,:)
      REAL(prec), ALLOCATABLE :: y(:,:)
      REAL(prec) :: dx, dy
      INTEGER, ALLOCATABLE    :: DOFtoIJ(:,:), IJtoDOF(:,:)

      ! MITgcm Data
      REAL(prec), ALLOCATABLE :: dxg(:,:)
      REAL(prec), ALLOCATABLE :: dyg(:,:)
      REAL(prec), ALLOCATABLE :: rac(:,:)
      REAL(prec), ALLOCATABLE :: hFacC(:,:,:)
      REAL(prec), ALLOCATABLE :: drf(:)
 
       
       CONTAINS

          GENERIC :: Build => Build_Mesh
          PROCEDURE, PRIVATE :: Build_Mesh
          PROCEDURE :: Free => Free_Mesh
    
          PROCEDURE :: ConstructDirichletCube
          PROCEDURE :: ConstructIrregularGeometryDemo
          PROCEDURE :: ConstructCircularGeometryDemo
          PROCEDURE :: ConstructWetPointMap
          PROCEDURE :: FlatMap
          PROCEDURE :: GridMap
          
          PROCEDURE :: Load => Load_Mesh
         PROCEDURE :: GenerateMaskFromhFac

          PROCEDURE :: WriteTecplot
          
    END TYPE Mesh


    
    REAL(prec), PARAMETER :: SLSpectra_wetValue = 1.0_prec
    REAL(prec), PARAMETER :: SLSpectra_dryValue = 0.0_prec
    REAL(prec), PARAMETER :: SLSpectra_fillValue = -99999.0_prec


 CONTAINS
 
  SUBROUTINE Build_Mesh( this, nX, nY, nZ )

   IMPLICIT NONE
   CLASS(Mesh), INTENT(out) :: this
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
      
 END SUBROUTINE Build_Mesh

 SUBROUTINE Free_Mesh( this )

   IMPLICIT NONE
   CLASS(Mesh), INTENT(inout) :: this

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

 END SUBROUTINE Free_Mesh
 
 SUBROUTINE Load_Mesh( this, mdpath )
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
   CLASS(Mesh), INTENT(out) :: this
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
    
   
  END SUBROUTINE Load_Mesh
  
  SUBROUTINE GenerateMaskFromhFac( this, x1, x2, y1, y2, klev )
    IMPLICIT NONE
    CLASS(Mesh),  INTENT(inout) :: this
    REAL(prec), INTENT(in) :: x1, x2, y1, y2
    INTEGER, INTENT(in) :: klev
    ! Local
    INTEGER :: i, j
    
      this % tracerMask = SLSpectra_dryValue
      DO j = 1, this % nY
        DO i = 1, this % nX
        
          IF( this % x(i,j) > x1 .AND. this % x(i,j) < x2 .AND. &
              this % y(i,j) > y1 .AND. this % y(i,j) < y2 .AND. &
              this % hFacC(i,j,klev) /= 0.0_prec )THEN
            this % tracerMask(i,j) = SLSpectra_wetValue
          ENDIF
          
        ENDDO
      ENDDO
  
  
      CALL this % ConstructWetPointMap()
      
  END SUBROUTINE GenerateMaskFromhFac
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
 
 SUBROUTINE ConstructIrregularGeometryDemo( this )
 !! Constructs a square domain on [0,1]^2 with the tracermask
 !! set to 0 on the south, east, north, and west boundaries.
 !! On the west boundary we create an irregular coastline with the
 !! tracermask.
 !!
 !! This mesh construction is meant to be used for testing purposes
 !!
   IMPLICIT NONE
   CLASS(Mesh), INTENT(inout) :: this
   ! Local
   INTEGER    :: i,j
   REAL(prec) :: dx, dy, coastX

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
       DO i = 1, this % nX
       
         ! West boundary
         coastX = 0.25_prec*exp( -(this % y(i,j) - 0.5_prec)**2/(2.0_prec*0.25_prec*0.25_prec) )
         IF( this % x(i,j) <= coastX )THEN
           this % tracermask(i,j) = SLSpectra_dryValue
         ELSE
           EXIT
         ENDIF
       ENDDO
       
       ! East Boundary
       this % tracermask(this % nX,j) = SLSpectra_dryValue
       
     ENDDO

     CALL this % ConstructWetPointMap()

 END SUBROUTINE ConstructIrregularGeometryDemo
 
  SUBROUTINE ConstructCircularGeometryDemo( this )
 !! Constructs a square domain on [0,1]^2 with the tracermask
 !! set to 0 on the south, east, north, and west boundaries.
 
 !!
 !! This mesh construction is meant to be used for testing purposes
 !!
   IMPLICIT NONE
   CLASS(Mesh), INTENT(inout) :: this
   ! Local
   INTEGER    :: i,j
   REAL(prec) :: dx, dy, r

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

     DO j = 1, this % nY
       DO i = 1, this % nX
       
         r = sqrt( (this % x(i,j) - 0.5_prec)**2 + (this % y(i,j) - 0.5_prec)**2 )
        
         IF( r > 0.4_prec )THEN
           this % tracermask(i,j) = SLSpectra_dryValue
         ENDIF
         
       ENDDO
     ENDDO

     CALL this % ConstructWetPointMap()

 END SUBROUTINE ConstructCircularGeometryDemo

 SUBROUTINE ConstructWetPointMap( this )
   IMPLICIT NONE
   CLASS( Mesh ), INTENT(inout) :: this
   ! LOCAL
   INTEGER :: i, j, nDOF


     PRINT*, "nX, nY (nCells): ", this % nX, this % nY, this % nX*this % nY

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
 
 SUBROUTINE WriteTecplot( this, ijArray, arrayName, filename )
   IMPLICIT NONE
   CLASS( Mesh ), INTENT(in) :: this
   REAL(prec), INTENT(in) :: ijArray(1:this % nX, 1:this % nY)
   CHARACTER(*), INTENT(in) :: arrayName
   CHARACTER(*), INTENT(in) :: filename
   ! Local
   INTEGER :: fUnit, i, j 
  
   
   OPEN( UNIT=NEWUNIT(fUnit), &
      FILE= TRIM(filename), &
      FORM='formatted', &
      STATUS='replace')

    WRITE(fUnit,*) 'VARIABLES = "X", "Y", "'//TRIM(arrayName)//'"'

      
    WRITE(fUnit,*) 'ZONE T="el00", I=',this % nX,', J=',this % nY

    DO j = 1, this % nY
      DO i = 1, this % nX

        WRITE(fUnit,'(3(ES16.7E3,1x))') this % x(i,j), this % y(i,j), ijArray(i,j)

      ENDDO
    ENDDO

    CLOSE(UNIT=fUnit)
   
   
 END SUBROUTINE WriteTecplot
 
 INTEGER FUNCTION NewUnit(thisunit)

    IMPLICIT NONE
    INTEGER,INTENT(out),OPTIONAL :: thisunit
    ! Local
    INTEGER,PARAMETER :: unitMin = 100,unitMax = 1000
    LOGICAL :: isopened
    INTEGER :: iUnit

    newunit = -1

    DO iUnit = unitMin,unitMax
      ! Check to see IF this UNIT is opened
      INQUIRE (UNIT=iUnit,opened=isopened)
      IF (.not. isopened) THEN
        newunit = iUnit
        EXIT
      END IF
    END DO

    IF (PRESENT(thisunit)) thisunit = newunit

  END FUNCTION NewUnit

END MODULE SLSpectra_Mesh
