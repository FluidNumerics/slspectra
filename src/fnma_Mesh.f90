MODULE fnma_Mesh

USE fnma_Precision

 IMPLICIT NONE

   TYPE fnmaMesh
      INTEGER :: nX, nY, nZ, nDOF
      INTEGER :: kLev
      REAL(prec), ALLOCATABLE :: mask(:,:)

      REAL(prec), ALLOCATABLE :: xc(:,:)
      REAL(prec), ALLOCATABLE :: yc(:,:)
      REAL(prec), ALLOCATABLE :: xg(:,:)
      REAL(prec), ALLOCATABLE :: yg(:,:)

      REAL(prec), ALLOCATABLE :: dxc(:,:)
      REAL(prec), ALLOCATABLE :: dyc(:,:)
      REAL(prec), ALLOCATABLE :: dxg(:,:)
      REAL(prec), ALLOCATABLE :: dyg(:,:)
      REAL(prec), ALLOCATABLE :: rac(:,:)
      REAL(prec), ALLOCATABLE :: raw(:,:)
      REAL(prec), ALLOCATABLE :: ras(:,:)
      REAL(prec), ALLOCATABLE :: hFacC(:,:,:)
      REAL(prec), ALLOCATABLE :: hFacS(:,:,:)
      REAL(prec), ALLOCATABLE :: hFacW(:,:,:)
      REAL(prec), ALLOCATABLE :: drf(:)
       
       CONTAINS

          GENERIC :: Build => Build_fnmaMesh
          PROCEDURE, PRIVATE :: Build_fnmaMesh
          PROCEDURE :: Free => Free_fnmaMesh
    
          PROCEDURE :: Demo_DirichletCube
!          PROCEDURE :: Demo_IrregularGeometry
!          PROCEDURE :: Demo_CircularGeometry
          
          PROCEDURE :: Load => Load_fnmaMesh
          PROCEDURE :: GenerateMaskFromhFac

!          PROCEDURE :: WriteTecplot
          
    END TYPE fnmaMesh


    
    REAL(prec), PARAMETER :: fnma_wetValue = 1.0_prec
    REAL(prec), PARAMETER :: fnma_dryValue = 0.0_prec
    REAL(prec), PARAMETER :: fnma_fillValue = -99999.0_prec


 CONTAINS
 
  SUBROUTINE Build_fnmaMesh( this, nX, nY, nZ )

   IMPLICIT NONE
   CLASS(fnmaMesh), INTENT(out) :: this
   INTEGER, INTENT(in)      :: nX, nY, nZ

      this % nX = nX
      this % nY = nY
      this % nZ = nZ
      
      ! Tracer mesh
      ALLOCATE( this % mask(1:nX,1:nY), &
                this % xc(1:nX,1:nY), &
                this % yc(1:nX,1:nY), &
                this % xg(1:nX,1:nY), &
                this % yg(1:nX,1:nY), &
                this % dxg(1:nX,1:nY), &
                this % dyg(1:nX,1:nY), &
                this % dxc(1:nX,1:nY), &
                this % dyc(1:nX,1:nY), &
                this % rac(1:nX,1:nY), &
                this % raw(1:nX,1:nY), &
                this % ras(1:nX,1:nY), &
                this % hFacC(1:nX,1:nY,1:nZ), &
                this % hFacW(1:nX,1:nY,1:nZ), &
                this % hFacS(1:nX,1:nY,1:nZ), &
                this % drf(1:nZ))
      
 END SUBROUTINE Build_fnmaMesh

 SUBROUTINE Free_fnmaMesh( this )

   IMPLICIT NONE
   CLASS(fnmaMesh), INTENT(inout) :: this

      DEALLOCATE( this % mask )
      DEALLOCATE( this % xc )
      DEALLOCATE( this % yc )
      DEALLOCATE( this % xg )
      DEALLOCATE( this % yg )
      DEALLOCATE( this % dxg )
      DEALLOCATE( this % dyg )
      DEALLOCATE( this % dxc )
      DEALLOCATE( this % dyc )
      DEALLOCATE( this % rac )
      DEALLOCATE( this % ras )
      DEALLOCATE( this % raw )
      DEALLOCATE( this % hFacC )
      DEALLOCATE( this % hFacS )
      DEALLOCATE( this % hFacW )
      DEALLOCATE( this % drf )

 END SUBROUTINE Free_fnmaMesh
 
 SUBROUTINE Load_fnmaMesh( this, mdpath )
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
 !!   XG.data 
 !!   XG.meta
 !!   YG.data
 !!   YG.meta
 !!   DXG.data
 !!   DXG.meta
 !!   DYG.data
 !!   DYG.meta
 !!   DRF.data
 !!   DRF.meta
 !!   RAC.data
 !!   RAC.meta
 !!   RAW.data
 !!   RAW.meta
 !!   RAS.data
 !!   RAS.meta
 !!   hFacC.data
 !!   hFacC.meta
 !!   hFacW.data
 !!   hFacW.meta
 !!   hFacS.data
 !!   hFacS.meta
 !!
   IMPLICIT NONE
   CLASS(fnmaMesh), INTENT(out) :: this
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
    READ(fUnit,REC=1) this % xc
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
    READ(fUnit,REC=1) this % yc
    CLOSE(fUnit)

    ! Read longitude
    OPEN( UNIT=NEWUNIT(fUnit), &
      FILE= TRIM(mdpath)//'/XG.data', &
      ACTION= 'READ', &
      FORM='UNFORMATTED', &
      STATUS='OLD', &
      ACCESS='DIRECT', &
      CONVERT='BIG_ENDIAN', &
      RECL=prec*this % nX*this % nY)
    READ(fUnit,REC=1) this % xg
    CLOSE(fUnit)
    
    ! Read latitude
    OPEN( UNIT=NEWUNIT(fUnit), &
      FILE= TRIM(mdpath)//'/YG.data', &
      ACTION= 'READ', &
      FORM='UNFORMATTED', &
      STATUS='OLD', &
      ACCESS='DIRECT', &
      CONVERT='BIG_ENDIAN', &
      RECL=prec*this % nX*this % nY)
    READ(fUnit,REC=1) this % yg
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
    
    ! Read x grid spacing
    OPEN( UNIT=NEWUNIT(fUnit), &
      FILE= TRIM(mdpath)//'/DXC.data', &
      ACTION= 'READ', &
      FORM='UNFORMATTED', &
      STATUS='OLD', &
      ACCESS='DIRECT', &
      CONVERT='BIG_ENDIAN', &
      RECL=prec*this % nX*this % nY)
    READ(fUnit,REC=1) this % dxc
    CLOSE(fUnit)
    PRINT*, "DXC MAX : ", MAXVAL(this % dxc)
    PRINT*, "DXC MIN : ", MINVAL(this % dxc)
    
    ! Read y grid spacing
    OPEN( UNIT=NEWUNIT(fUnit), &
      FILE= TRIM(mdpath)//'/DYC.data', &
      ACTION= 'READ', &
      FORM='UNFORMATTED', &
      STATUS='OLD', &
      ACCESS='DIRECT', &
      CONVERT='BIG_ENDIAN', &
      RECL=prec*this % nX*this % nY)
    READ(fUnit,REC=1) this % dyc
    CLOSE(fUnit)
    PRINT*, "DYC MAX : ", MAXVAL(this % dyc)
    PRINT*, "DYC MIN : ", MINVAL(this % dyc)
    
    ! Read z grid spacing
    OPEN( UNIT=NEWUNIT(fUnit), &
      FILE= TRIM(mdpath)//'/DRF.data', &
      ACTION= 'READ', &
      FORM='UNFORMATTED', &
      STATUS='OLD', &
      ACCESS='DIRECT', &
      CONVERT='BIG_ENDIAN', &
      RECL=prec*this % nZ)
    READ(fUnit,REC=1) this % drf
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

    ! Read hFacS
    OPEN( UNIT=NEWUNIT(fUnit), &
      FILE= TRIM(mdpath)//'/hFacS.data', &
      ACTION= 'READ', &
      FORM='UNFORMATTED', &
      STATUS='OLD', &
      ACCESS='DIRECT', &
      CONVERT='BIG_ENDIAN', &
      RECL=prec*this % nX*this % nY*this % nZ)
    READ(fUnit,REC=1) this % hFacS
    CLOSE(fUnit)

    ! Read hFacW
    OPEN( UNIT=NEWUNIT(fUnit), &
      FILE= TRIM(mdpath)//'/hFacW.data', &
      ACTION= 'READ', &
      FORM='UNFORMATTED', &
      STATUS='OLD', &
      ACCESS='DIRECT', &
      CONVERT='BIG_ENDIAN', &
      RECL=prec*this % nX*this % nY*this % nZ)
    READ(fUnit,REC=1) this % hFacW
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
    PRINT*, "RAC MAX : ", MAXVAL(this % rac)
    PRINT*, "RAC MIN : ", MINVAL(this % rac)
    
    ! Read cell area
    OPEN( UNIT=NEWUNIT(fUnit), &
      FILE= TRIM(mdpath)//'/RAW.data', &
      ACTION= 'READ', &
      FORM='UNFORMATTED', &
      STATUS='OLD', &
      ACCESS='DIRECT', &
      CONVERT='BIG_ENDIAN', &
      RECL=prec*this % nX*this % nY)
    READ(fUnit,REC=1) this % raw
    CLOSE(fUnit)
    PRINT*, "RAW MAX : ", MAXVAL(this % raw)
    PRINT*, "RAW MIN : ", MINVAL(this % raw)
    
    ! Read cell area
    OPEN( UNIT=NEWUNIT(fUnit), &
      FILE= TRIM(mdpath)//'/RAS.data', &
      ACTION= 'READ', &
      FORM='UNFORMATTED', &
      STATUS='OLD', &
      ACCESS='DIRECT', &
      CONVERT='BIG_ENDIAN', &
      RECL=prec*this % nX*this % nY)
    READ(fUnit,REC=1) this % ras
    CLOSE(fUnit)
    PRINT*, "RAS MAX : ", MAXVAL(this % ras)
    PRINT*, "RAS MIN : ", MINVAL(this % ras)
    
  END SUBROUTINE Load_fnmaMesh
  
  SUBROUTINE GenerateMaskFromhFac( this, x1, x2, y1, y2, klev )
    IMPLICIT NONE
    CLASS(fnmaMesh),  INTENT(inout) :: this
    REAL(prec), INTENT(in) :: x1, x2, y1, y2
    INTEGER, INTENT(in) :: klev
    ! Local
    INTEGER :: i, j
    
      this % mask = fnma_dryValue
      DO j = 1, this % nY
        DO i = 1, this % nX
        
          IF( this % xc(i,j) > x1 .AND. this % xc(i,j) < x2 .AND. &
              this % yc(i,j) > y1 .AND. this % yc(i,j) < y2 .AND. &
              this % hFacC(i,j,klev) /= 0.0_prec )THEN
            this % mask(i,j) = fnma_wetValue
          ENDIF
          
        ENDDO
      ENDDO
      
  END SUBROUTINE GenerateMaskFromhFac
!
 SUBROUTINE Demo_DirichletCube( this )
 !! Constructs a cube domain on [0,1]^2 with the mask
 !! set to 0 on the south, east, north, and west boundaries.
 !!
 !! This mesh construction is meant to be used for testing purposes
 !!
   IMPLICIT NONE
   CLASS(fnmaMesh), INTENT(inout) :: this
   ! Local
   INTEGER    :: i,j,k
   REAL(prec) :: dx, dy

     dx = 1.0_prec/REAL(this % nX,prec)
     dy = 1.0_prec/REAL(this % nY,prec)

     DO j = 1, this % nY
       DO i = 1, this % nX
         this % xg(i,j) = (REAL(i-1,prec))*dx
         this % yg(i,j) = (REAL(j-1,prec))*dy
         this % xc(i,j) = this % xg(i,j) + 0.5_prec*dx
         this % yc(i,j) = this % yg(i,j) + 0.5_prec*dy
       ENDDO
     ENDDO

     this % mask = fnma_wetValue
     this % hFacC = 1.0_prec
     this % hFacS = 1.0_prec
     this % hFacW = 1.0_prec

     ! South and North Boundaries
     DO k = 1, this % nZ
       DO i = 1, this % nX
         this % hFacC(i,1,k) = 0.0_prec
         this % hFacC(i,this % nY,k) = 0.0_prec
         this % hFacW(i,1,k) = 0.0_prec
         this % hFacW(i,this % nY,k) = 0.0_prec
         this % hFacS(i,1,k) = 0.0_prec
         this % hFacS(i,2,k) = 0.0_prec
         this % hFacS(i,this % nY,k) = 0.0_prec
       ENDDO
     ENDDO
     
     DO i = 1, this % nX
       this % mask(i,1) = fnma_dryValue
       this % mask(i,this % nY) = fnma_dryValue
     ENDDO

     ! East and West Boundaries
     DO k = 1, this % nZ
       DO j = 1, this % nY
         this % hFacC(1,j,k) = 0.0_prec
         this % hFacC(this % nX,j,k) = 0.0_prec
         this % hFacS(1,j,k) = 0.0_prec
         this % hFacS(this % nX,j,k) = 0.0_prec

         this % hFacW(1,j,k) = 0.0_prec
         this % hFacW(2,j,k) = 0.0_prec
         this % hFacW(this % nX,j,k) = 0.0_prec
       ENDDO
     ENDDO

     DO j = 1, this % nY
       this % mask(1,j) = fnma_dryValue
       this % mask(this % nX,j) = fnma_dryValue
     ENDDO

 END SUBROUTINE Demo_DirichletCube
 
 !SUBROUTINE Demo_IrregularGeometry( this )
 !!! Constructs a square domain on [0,1]^2 with the mask
 !!! set to 0 on the south, east, north, and west boundaries.
 !!! On the west boundary we create an irregular coastline with the
 !!! mask.
 !!!
 !!! This mesh construction is meant to be used for testing purposes
 !!!
 !  IMPLICIT NONE
 !  CLASS(fnmaMesh), INTENT(inout) :: this
 !  ! Local
 !  INTEGER    :: i,j
 !  REAL(prec) :: dx, dy, coastX

 !    dx = 1.0_prec/REAL(this % nX,prec)
 !    dy = 1.0_prec/REAL(this % nY,prec)
 !    this % dx = dx
 !    this % dy = dy

 !    DO j = 1, this % nY
 !      DO i = 1, this % nX
 !        this % x(i,j) = (REAL(i,prec) + 0.5_prec)*dx
 !        this % y(i,j) = (REAL(j,prec) + 0.5_prec)*dy
 !      ENDDO
 !    ENDDO

 !    this % tracerMask = fnma_wetValue
 !    ! South and North Boundaries
 !    DO i = 1, this % nX
 !      this % mask(i,1) = fnma_dryValue
 !      this % mask(i,this % nY) = fnma_dryValue
 !    ENDDO

 !    ! East and West Boundaries
 !    DO j = 1, this % nY
 !      DO i = 1, this % nX
 !      
 !        ! West boundary
 !        coastX = 0.25_prec*exp( -(this % y(i,j) - 0.5_prec)**2/(2.0_prec*0.25_prec*0.25_prec) )
 !        IF( this % x(i,j) <= coastX )THEN
 !          this % mask(i,j) = fnma_dryValue
 !        ELSE
 !          EXIT
 !        ENDIF
 !      ENDDO
 !      
 !      ! East Boundary
 !      this % mask(this % nX,j) = fnma_dryValue
 !      
 !    ENDDO

 !END SUBROUTINE Demo_IrregularGeometry
 !
 ! SUBROUTINE Demo_CircularGeometry( this )
 !!! Constructs a square domain on [0,1]^2 with the mask
 !!! set to 0 on the south, east, north, and west boundaries.
 !
 !!!
 !!! This mesh construction is meant to be used for testing purposes
 !!!
 !  IMPLICIT NONE
 !  CLASS(fnmaMesh), INTENT(inout) :: this
 !  ! Local
 !  INTEGER    :: i,j
 !  REAL(prec) :: dx, dy, r

 !    dx = 1.0_prec/REAL(this % nX,prec)
 !    dy = 1.0_prec/REAL(this % nY,prec)
 !    this % dx = dx
 !    this % dy = dy

 !    DO j = 1, this % nY
 !      DO i = 1, this % nX
 !        this % x(i,j) = (REAL(i,prec) + 0.5_prec)*dx
 !        this % y(i,j) = (REAL(j,prec) + 0.5_prec)*dy
 !      ENDDO
 !    ENDDO

 !    this % tracerMask = fnma_wetValue

 !    DO j = 1, this % nY
 !      DO i = 1, this % nX
 !      
 !        r = sqrt( (this % x(i,j) - 0.5_prec)**2 + (this % y(i,j) - 0.5_prec)**2 )
 !       
 !        IF( r > 0.4_prec )THEN
 !          this % mask(i,j) = fnma_dryValue
 !        ENDIF
 !        
 !      ENDDO
 !    ENDDO

 !END SUBROUTINE Demo_CircularGeometry

 !SUBROUTINE WriteTecplot( this, ijArray, arrayName, filename )
 !  IMPLICIT NONE
 !  CLASS( fnmaMesh ), INTENT(in) :: this
 !  REAL(prec), INTENT(in) :: ijArray(1:this % nX, 1:this % nY)
 !  CHARACTER(*), INTENT(in) :: arrayName
 !  CHARACTER(*), INTENT(in) :: filename
 !  ! Local
 !  INTEGER :: fUnit, i, j 
 ! 
 !  
 !  OPEN( UNIT=NEWUNIT(fUnit), &
 !     FILE= TRIM(filename), &
 !     FORM='formatted', &
 !     STATUS='replace')

 !   WRITE(fUnit,*) 'VARIABLES = "X", "Y", "'//TRIM(arrayName)//'"'

 !     
 !   WRITE(fUnit,*) 'ZONE T="el00", I=',this % nX,', J=',this % nY

 !   DO j = 1, this % nY
 !     DO i = 1, this % nX

 !       WRITE(fUnit,'(3(ES16.7E3,1x))') this % x(i,j), this % y(i,j), ijArray(i,j)

 !     ENDDO
 !   ENDDO

 !   CLOSE(UNIT=fUnit)
 !  
 !  
 !END SUBROUTINE WriteTecplot
 
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

END MODULE fnma_Mesh
