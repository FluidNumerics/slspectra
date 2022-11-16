PROGRAM MITgcmGeometry

USE ISO_FORTRAN_ENV
USE SLSpectra_Precision
USE SLSpectra_SupportRoutines
USE SLSpectra_Mesh
USE SLSpectra_Stencil
USE SLSpectra_AdjacencyGraph
USE SLSpectra_Generators
USE SLSpectra_Model

IMPLICIT NONE

  
  INTEGER, PARAMETER :: klev =  25 ! Vertical level of interest in the MITgcm Mesh
  INTEGER, PARAMETER :: maxValence = 5
  INTEGER, PARAMETER :: maxEvalWrite = 100
  REAL(prec), PARAMETER :: x1 = 5.5_prec
  REAL(prec), PARAMETER :: x2 = 7.0_prec
  REAL(prec), PARAMETER :: y1 = 34.5_prec
  REAL(prec), PARAMETER :: y2 = 36.0_prec
  REAL(prec), PARAMETER :: pi = 4.0_prec*atan(1.0_prec)
  TYPE(Laplacian5Stencil) :: modelStencil
  TYPE(Mesh), TARGET :: modelMesh
  TYPE(Generator) :: laplacian
  TYPE(AdjacencyGraph), TARGET :: diagnosisGraph
  TYPE(AdjacencyGraph), TARGET :: overlapGraph
  TYPE(Model) :: slModel
  CHARACTER(LEN=255) :: MDPath
  REAL(prec), ALLOCATABLE :: v(:,:)
  REAL(prec), ALLOCATABLE :: projection(:,:)
  REAL(prec), ALLOCATABLE :: residual(:,:)
  REAL(prec), ALLOCATABLE :: modalCoeffs(:)
  INTEGER :: fUnit, i, j 
  
     
  ! Create a stencil
  CALL modelStencil % Build()
  
  CALL get_environment_variable("SLS_MDPATH", MDPath)

  ! Create a mesh
  CALL modelMesh % Load(MDPath)
  
  CALL modelMesh % GenerateMaskFromhFac( x1, x2, y1, y2, klev )
  
  CALL modelMesh % WriteTecplot( modelMesh % tracerMask, 'mask', 'mask.tec' )

  ! Associate the mesh with the Generator
  CALL laplacian % AssociateMesh( modelMesh )

  ! Point operator to the 5-point laplacian with neumann boundary condition
  laplacian % SLOperator => SLOperator_MITgcmNeumann
  ! Point the inner product to the MITgcm inner product
  laplacian % SLInnerProduct => SLInnerProduct_MITgcm
  
  ! Create an adjacency graph using the mesh and the stencil
  CALL diagnosisGraph % Build( modelMesh, modelStencil % firstOrder )
  CALL overlapGraph % Build( modelMesh, modelStencil % secondOrder )
  
  ! Create the model
  CALL slModel % Init( diagnosisGraph, overlapGraph, laplacian )

  CALL slModel % CalculateEigenModes( )

!  CALL slModel % VerifyMutualOrthogonality( )

  ! Write the Eigenvectors to file (first 100 modes)
!  CALL slModel % WriteTecplot( 100 )
  
  ALLOCATE( v(1:modelMesh % nX, modelMesh % nY),&
            projection(1:modelMesh % nX, modelMesh % nY),&
            residual(1:modelMesh % nX, modelMesh % nY),&
            modalCoeffs(1:modelMesh % nDOF) )

  ! Project example data onto the eigenvectors
  ! > Example - create a sinusoid with known wavelength
  !v = slModel % slGenerator % modelMesh % GridMap( slModel % evectors(:,modelMesh % nDOF) )
  DO j = 1, modelMesh % nY
    DO i = 1, modelMesh % nX
      v(i,j) = sin( 2.0*pi*modelMesh % x(i,j)/10.0_prec )*sin( 2.0*pi*10.0*modelMesh % y(i,j)/10.0_prec )
    ENDDO
  ENDDO

  ! Get the modal projections
  CALL slModel % ModalProjection( v, modalCoeffs, projection, residual )

  ! Write the data, projection, and residual to tecplot
  CALL slModel % slGenerator % modelMesh % WriteTecplot(v, 'data', 'sinusoid.tec' )
  CALL slModel % slGenerator % modelMesh % WriteTecplot(projection, 'data', 'sinusoid-projection.tec' )
  CALL slModel % slGenerator % modelMesh % WriteTecplot(residual, 'data', 'sinusoid-residual.tec' )


  ! Write the  modal coefficients
  OPEN( UNIT=NEWUNIT(fUnit), &
    FILE= 'sin-10-modal.curve', &
    FORM='formatted', &
    STATUS='replace')
  
  DO i = 1, modelMesh % nDOF
    WRITE(fUnit,*) slModel % eValues(i), modalCoeffs(i)
  ENDDO
  CLOSE(UNIT=fUnit)
    

  DEALLOCATE(v, projection, residual, modalCoeffs )

  CALL modelMesh % Free()
  CALL modelStencil % Free()
  CALL diagnosisGraph % Free()
  CALL overlapGraph % Free()
  CALL slModel % Free()
  
END PROGRAM MITgcmGeometry
