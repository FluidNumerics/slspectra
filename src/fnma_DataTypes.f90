MODULE fnma_DataTypes

! The fnma_dataObj class is adapted from https://github.com/FluidNumerics/SELF/blob/main/src/SELF_Data.f90

!
!
! Indexing layout relative to tracer cell (t(i,j))
! 
!  > Vorticity point (i,j) is at the south-west corner
!  > U-velocity point (i,j) is centered on the west edge
!  > V-velocity point (i,j) is centered on the south edge
!  
! z(i,j+1) -- v(i,j+1) -- z(i+1,j+1) 
!    |                         |
!    |                         |
!  u(i,j)     t(i,j)       u(i+1,j) 
!    |                         |
!    |                         |
!  z(i,j) --- v(i,j) --- z(i+1,j) 
        
USE fnma_Precision
USE fnma_Mesh
USE fnma_Metadata

  TYPE,PUBLIC :: fnmaDataObj
    INTEGER :: nVar, nX, nY
    TYPE(fnmaMesh), POINTER :: mesh
    TYPE(fnmaMetadata), ALLOCATABLE :: meta(:)
    REAL(prec), ALLOCATABLE :: var(:,:,:)
    
    CONTAINS
    PROCEDURE,PUBLIC :: Init => Init_DataObj
    PROCEDURE,PUBLIC :: Free => Free_DataObj

    ! Procedures for setting metadata for
    PROCEDURE,PUBLIC :: SetName => SetName_DataObj
    PROCEDURE,PUBLIC :: SetDescription => SetDescription_DataObj
    PROCEDURE,PUBLIC :: SetUnits => SetUnits_DataObj


  END TYPE fnmaDataObj

  TYPE,EXTENDS(fnmaDataObj) :: fnmaArCT
    !! *Ar*akawa *C* *T*racer point (fnmaArCT)

    CONTAINS
    PROCEDURE, PUBLIC :: diffX => diffX_fnmaArCT
    PROCEDURE, PUBLIC :: diffY => diffY_fnmaArCT

    ! Write to tecplot
    PROCEDURE,PUBLIC :: WriteTecplot => WriteTecplot_fnmaArCT

  END TYPE fnmaArCT

  TYPE,EXTENDS(fnmaDataObj) :: fnmaArCU
    !! *Ar*akawa *C* *U*-velocity point (fnmaArCU)

    CONTAINS
    PROCEDURE, PUBLIC :: diffX => diffX_fnmaArCU
    PROCEDURE, PUBLIC :: diffY => diffY_fnmaArCU

    ! Write to tecplot
    PROCEDURE,PUBLIC :: WriteTecplot => WriteTecplot_fnmaArCU

  END TYPE fnmaArCU

  TYPE,EXTENDS(fnmaDataObj) :: fnmaArCV
     !! *Ar*akawa *C* *V*-velocity point (fnmaArCV)

     CONTAINS
     PROCEDURE, PUBLIC :: diffX => diffX_fnmaArCV
     PROCEDURE, PUBLIC :: diffY => diffY_fnmaArCV

    ! Write to tecplot
    PROCEDURE,PUBLIC :: WriteTecplot => WriteTecplot_fnmaArCV

  END TYPE fnmaArCV

  TYPE,EXTENDS(fnmaDataObj) :: fnmaArCZ
     !! *Ar*akawa *C* *Z*eta (voriticity) point (fnmaArCZ)

     CONTAINS
     PROCEDURE, PUBLIC :: diffX => diffX_fnmaArCZ
     PROCEDURE, PUBLIC :: diffY => diffY_fnmaArCZ

     ! Write to tecplot
     PROCEDURE,PUBLIC :: WriteTecplot => WriteTecplot_fnmaArCZ

  END TYPE fnmaArCZ

CONTAINS

! -- DataObj -- !

  SUBROUTINE Init_DataObj(fnmaStorage,mesh,nX,nY,nVar)
    IMPLICIT NONE
    CLASS(fnmaDataObj),INTENT(out)  :: fnmaStorage
    TYPE(fnmaMesh),INTENT(in),TARGET :: mesh
    INTEGER,INTENT(in) :: nX
    INTEGER,INTENT(in) :: nY
    INTEGER,INTENT(in) :: nVar

    fnmaStorage % nX = nX
    fnmaStorage % nY = nY
    fnmaStorage % nVar = nVar
    fnmaStorage % mesh => mesh
    ALLOCATE( fnmaStorage % meta(1:nVar) )
    ALLOCATE( fnmaStorage % var(0:nX+1,0:nY+1,1:nVar) )

  END SUBROUTINE Init_DataObj

  SUBROUTINE Free_DataObj(fnmaStorage)
    IMPLICIT NONE
    CLASS(fnmaDataObj),INTENT(inout) :: fnmaStorage

    fnmaStorage % mesh => NULL()
    DEALLOCATE( fnmaStorage % meta )

  END SUBROUTINE Free_DataObj

  SUBROUTINE SetName_DataObj(fnmaStorage,ivar,name)
    !! Set the name of the `ivar-th` variable
    IMPLICIT NONE
    CLASS(fnmaDataObj),INTENT(inout) :: fnmaStorage
    INTEGER,INTENT(in) :: ivar
    CHARACTER(*),INTENT(in) :: name

    CALL fnmaStorage % meta(ivar) % SetName(name) 

  END SUBROUTINE SetName_DataObj

  SUBROUTINE SetDescription_DataObj(fnmaStorage,ivar,description)
    !! Set the description of the `ivar-th` variable
    IMPLICIT NONE
    CLASS(fnmaDataObj),INTENT(inout) :: fnmaStorage
    INTEGER,INTENT(in) :: ivar
    CHARACTER(*),INTENT(in) :: description

    CALL fnmaStorage % meta(ivar) % SetDescription(description) 

  END SUBROUTINE SetDescription_DataObj

  SUBROUTINE SetUnits_DataObj(fnmaStorage,ivar,units)
    !! Set the units of the `ivar-th` variable
    IMPLICIT NONE
    CLASS(fnmaDataObj),INTENT(inout) :: fnmaStorage
    INTEGER,INTENT(in) :: ivar
    CHARACTER(*),INTENT(in) :: units

    CALL fnmaStorage % meta(ivar) % SetUnits(units) 

  END SUBROUTINE SetUnits_DataObj

  !! Derivative Operators !!

  ! > Tracer Points
  SUBROUTINE DiffX_fnmaArCT( f, dfdx )
    !! Calculates the derivative of "f", stored
    !! at Arakawa C Tracer points, and returns
    !! "dfdx", stored at Arakawa C U-velocity points.
    IMPLICIT NONE
    CLASS(fnmaArCT), INTENT(in)   :: f
    TYPE(fnmaArCU), INTENT(inout) :: dfdx
    ! Local
    INTEGER :: i, j, iVar

    DO iVar = 1, f % nVar
      DO j = 1, f % nY
        DO i = 1, f % nX
          dfdx % var(i,j,iVar) = ( f % var(i,j,iVar) - &
                  f % var(i-1,j,iVar) ) / f % mesh % dxc(i,j)
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE DiffX_fnmaArCT

  SUBROUTINE DiffY_fnmaArCT( f, dfdy )
    !! Calculates the derivative of "f", stored
    !! at Arakawa C Tracer points, and returns
    !! "dfdy", stored at Arakawa C V-velocity points.
    IMPLICIT NONE
    CLASS(fnmaArCT), INTENT(in)   :: f
    TYPE(fnmaArCV), INTENT(inout) :: dfdy
    ! Local
    INTEGER :: i, j, iVar

    DO iVar = 1, f % nVar
      DO j = 1, f % nY
        DO i = 1, f % nX
          dfdy % var(i,j,iVar) = ( f % var(i,j,iVar) - &
                  f % var(i,j-1,iVar) ) / f % mesh % dyc(i,j)
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE DiffY_fnmaArCT

  ! > U-velocity Points
  SUBROUTINE DiffX_fnmaArCU( f, dfdx )
    !! Calculates the derivative of "f", stored
    !! at Arakawa U-velocity points, and returns
    !! "dfdx", stored at Arakawa C tracer points.
    IMPLICIT NONE
    CLASS(fnmaArCU), INTENT(in)   :: f
    TYPE(fnmaArCT), INTENT(inout) :: dfdx
    ! Local
    INTEGER :: i, j, iVar

    DO iVar = 1, f % nVar
      DO j = 1, f % nY
        DO i = 1, f % nX
          dfdx % var(i,j,iVar) = ( f % var(i+1,j,iVar) - &
                  f % var(i,j,iVar) ) / f % mesh % dxg(i,j)
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE DiffX_fnmaArCU

  SUBROUTINE DiffY_fnmaArCU( f, dfdy )
    !! Calculates the derivative of "f", stored
    !! at Arakawa C U-velocity points, and returns
    !! "dfdy", stored at Arakawa C Z (vorticity) points.
    IMPLICIT NONE
    CLASS(fnmaArCU), INTENT(in)   :: f
    TYPE(fnmaArCZ), INTENT(inout) :: dfdy
    ! Local
    INTEGER :: i, j, iVar

    DO iVar = 1, f % nVar
      DO j = 1, f % nY
        DO i = 1, f % nX
          dfdy % var(i,j,iVar) = ( f % var(i,j,iVar) - &
                  f % var(i,j-1,iVar) ) / f % mesh % dyc(i,j)
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE DiffY_fnmaArCU

  ! > V-velocity Points
  SUBROUTINE DiffX_fnmaArCV( f, dfdx )
    !! Calculates the derivative of "f", stored
    !! at Arakawa V-velocity points, and returns
    !! "dfdx", stored at Arakawa C Z (vorticity) points.
    IMPLICIT NONE
    CLASS(fnmaArCV), INTENT(in)   :: f
    TYPE(fnmaArCZ), INTENT(inout) :: dfdx
    ! Local
    INTEGER :: i, j, iVar

    DO iVar = 1, f % nVar
      DO j = 1, f % nY
        DO i = 1, f % nX
          dfdx % var(i,j,iVar) = ( f % var(i,j,iVar) - &
                  f % var(i-1,j,iVar) ) / f % mesh % dxc(i,j)
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE DiffX_fnmaArCV

  SUBROUTINE DiffY_fnmaArCV( f, dfdy )
    !! Calculates the derivative of "f", stored
    !! at Arakawa C V-velocity points, and returns
    !! "dfdy", stored at Arakawa C Tracer points.
    IMPLICIT NONE
    CLASS(fnmaArCV), INTENT(in)   :: f
    TYPE(fnmaArCT), INTENT(inout) :: dfdy
    ! Local
    INTEGER :: i, j, iVar

    DO iVar = 1, f % nVar
      DO j = 1, f % nY
        DO i = 1, f % nX
          dfdy % var(i,j,iVar) = ( f % var(i,j+1,iVar) - &
                  f % var(i,j,iVar) ) / f % mesh % dyg(i,j)
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE DiffY_fnmaArCV

  ! > Z vorticity Points
  SUBROUTINE DiffX_fnmaArCZ( f, dfdx )
    !! Calculates the derivative of "f", stored
    !! at Arakawa C Z (vorticity) points, and returns
    !! "dfdx", stored at Arakawa C V-velocity points.
    IMPLICIT NONE
    CLASS(fnmaArCZ), INTENT(in)   :: f
    TYPE(fnmaArCV), INTENT(inout) :: dfdx
    ! Local
    INTEGER :: i, j, iVar

    DO iVar = 1, f % nVar
      DO j = 1, f % nY
        DO i = 1, f % nX
          dfdx % var(i,j,iVar) = ( f % var(i+1,j,iVar) - &
                  f % var(i,j,iVar) ) / f % mesh % dxg(i,j)
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE DiffX_fnmaArCZ

  SUBROUTINE DiffY_fnmaArCZ( f, dfdy )
    !! Calculates the derivative of "f", stored
    !! at Arakawa C Z (vorticity) points, and returns
    !! "dfdy", stored at Arakawa C U-velocity points.
    IMPLICIT NONE
    CLASS(fnmaArCZ), INTENT(in)   :: f
    TYPE(fnmaArCU), INTENT(inout) :: dfdy
    ! Local
    INTEGER :: i, j, iVar

    DO iVar = 1, f % nVar
      DO j = 1, f % nY
        DO i = 1, f % nX
          dfdy % var(i,j,iVar) = ( f % var(i,j+1,iVar) - &
                  f % var(i,j,iVar) ) / f % mesh % dyg(i,j)
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE DiffY_fnmaArCZ

  ! Tecplot IO

  SUBROUTINE WriteTecplot_fnmaArCT(this, filename)
    IMPLICIT NONE
    CLASS(fnmaArCT), INTENT(inout) :: this
    CHARACTER(*), INTENT(in) :: filename
    ! Local
    INTEGER :: fUnit
    INTEGER :: i, j, iVar
    CHARACTER(LEN=256) :: tecHeader
    CHARACTER(LEN=256) :: fmat

    OPEN( UNIT=NEWUNIT(fUnit), &
     FILE= TRIM(filename), &
     FORM='formatted', &
     STATUS='replace')

    tecHeader = 'VARIABLES = "X", "Y"'
    DO iVar = 1, this % nVar
      tecHeader = TRIM(tecHeader)//', "'//TRIM(this % meta(iVar) % name)//'"'
    ENDDO

    WRITE(fUnit,*) TRIM(tecHeader)

    ! Create format statement
    WRITE(fmat,*) this % nvar+2
    fmat = '('//TRIM(fmat)//'(ES16.7E3,1x))'

    WRITE(fUnit,*) 'ZONE T="el0", I=',this % mesh % nx,&
                               ', J=',this % mesh % ny
    DO j = 1, this % mesh % ny
      DO i = 1, this % mesh % nx

          WRITE(fUnit,fmat) this % mesh % xc(i,j),&
                            this % mesh % yc(i,j),&
                            this % var(i,j,1:this % nvar) 

      ENDDO
    ENDDO

    CLOSE(UNIT=fUnit)

  END SUBROUTINE WriteTecplot_fnmaArCT

  SUBROUTINE WriteTecplot_fnmaArCU(this, filename)
    IMPLICIT NONE
    CLASS(fnmaArCU), INTENT(inout) :: this
    CHARACTER(*), INTENT(in) :: filename
    ! Local
    INTEGER :: fUnit
    INTEGER :: i, j, iVar
    CHARACTER(LEN=256) :: tecHeader
    CHARACTER(LEN=256) :: fmat

    OPEN( UNIT=NEWUNIT(fUnit), &
     FILE= TRIM(filename), &
     FORM='formatted', &
     STATUS='replace')

    tecHeader = 'VARIABLES = "X", "Y"'
    DO iVar = 1, this % nVar
      tecHeader = TRIM(tecHeader)//', "'//TRIM(this % meta(iVar) % name)//'"'
    ENDDO

    WRITE(fUnit,*) TRIM(tecHeader)

    ! Create format statement
    WRITE(fmat,*) this % nvar+2
    fmat = '('//TRIM(fmat)//'(ES16.7E3,1x))'

    WRITE(fUnit,*) 'ZONE T="el0", I=',this % mesh % nx,&
                               ', J=',this % mesh % ny
    DO j = 1, this % mesh % ny
      DO i = 1, this % mesh % nx

          WRITE(fUnit,fmat) this % mesh % xg(i,j),&
                            this % mesh % yc(i,j),&
                            this % var(i,j,1:this % nvar) 

      ENDDO
    ENDDO

    CLOSE(UNIT=fUnit)

  END SUBROUTINE WriteTecplot_fnmaArCU

  SUBROUTINE WriteTecplot_fnmaArCV(this, filename)
    IMPLICIT NONE
    CLASS(fnmaArCV), INTENT(inout) :: this
    CHARACTER(*), INTENT(in) :: filename
    ! Local
    INTEGER :: fUnit
    INTEGER :: i, j, iVar
    CHARACTER(LEN=256) :: tecHeader
    CHARACTER(LEN=256) :: fmat

    OPEN( UNIT=NEWUNIT(fUnit), &
     FILE= TRIM(filename), &
     FORM='formatted', &
     STATUS='replace')

    tecHeader = 'VARIABLES = "X", "Y"'
    DO iVar = 1, this % nVar
      tecHeader = TRIM(tecHeader)//', "'//TRIM(this % meta(iVar) % name)//'"'
    ENDDO

    WRITE(fUnit,*) TRIM(tecHeader)

    ! Create format statement
    WRITE(fmat,*) this % nvar+2
    fmat = '('//TRIM(fmat)//'(ES16.7E3,1x))'

    WRITE(fUnit,*) 'ZONE T="el0", I=',this % mesh % nx,&
                               ', J=',this % mesh % ny
    DO j = 1, this % mesh % ny
      DO i = 1, this % mesh % nx

          WRITE(fUnit,fmat) this % mesh % xc(i,j),&
                            this % mesh % yg(i,j),&
                            this % var(i,j,1:this % nvar) 

      ENDDO
    ENDDO

    CLOSE(UNIT=fUnit)

  END SUBROUTINE WriteTecplot_fnmaArCV

  SUBROUTINE WriteTecplot_fnmaArCZ(this, filename)
    IMPLICIT NONE
    CLASS(fnmaArCZ), INTENT(inout) :: this
    CHARACTER(*), INTENT(in) :: filename
    ! Local
    INTEGER :: fUnit
    INTEGER :: i, j, iVar
    CHARACTER(LEN=256) :: tecHeader
    CHARACTER(LEN=256) :: fmat

    OPEN( UNIT=NEWUNIT(fUnit), &
     FILE= TRIM(filename), &
     FORM='formatted', &
     STATUS='replace')

    tecHeader = 'VARIABLES = "X", "Y"'
    DO iVar = 1, this % nVar
      tecHeader = TRIM(tecHeader)//', "'//TRIM(this % meta(iVar) % name)//'"'
    ENDDO

    WRITE(fUnit,*) TRIM(tecHeader)

    ! Create format statement
    WRITE(fmat,*) this % nvar+2
    fmat = '('//TRIM(fmat)//'(ES16.7E3,1x))'

    WRITE(fUnit,*) 'ZONE T="el0", I=',this % mesh % nx,&
                               ', J=',this % mesh % ny
    DO j = 1, this % mesh % ny
      DO i = 1, this % mesh % nx

          WRITE(fUnit,fmat) this % mesh % xg(i,j),&
                            this % mesh % yg(i,j),&
                            this % var(i,j,1:this % nvar) 

      ENDDO
    ENDDO

    CLOSE(UNIT=fUnit)

  END SUBROUTINE WriteTecplot_fnmaArCZ


END MODULE fnma_DataTypes
