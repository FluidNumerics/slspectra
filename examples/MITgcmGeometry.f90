PROGRAM IrregularGeometry

USE ISO_FORTRAN_ENV
USE SLSpectra_Precision
USE SLSpectra_SupportRoutines
USE SLSpectra_Mesh
USE SLSpectra_Stencil
USE SLSpectra_AdjacencyGraph
USE SLSpectra_Generators

IMPLICIT NONE

  
  INTEGER, PARAMETER :: klev =  25 ! Vertical level of interest in the MITgcm Mesh
  INTEGER, PARAMETER :: maxValence = 5
  INTEGER, PARAMETER :: maxEvalWrite = 100
  REAL(prec), PARAMETER :: x1 = 5.5_prec
  REAL(prec), PARAMETER :: x2 = 8.0_prec
  REAL(prec), PARAMETER :: y1 = 34.5_prec
  REAL(prec), PARAMETER :: y2 = 36.0_prec
  REAL(prec), PARAMETER :: pi = 4.0_prec*atan(1.0_prec)
  TYPE(Laplacian5Stencil) :: modelStencil
  TYPE(Mesh), TARGET :: modelMesh
  TYPE(Generator) :: laplacian
  TYPE(AdjacencyGraph) :: diagnosisGraph
  TYPE(AdjacencyGraph) :: overlapGraph
  REAL(prec), ALLOCATABLE :: impulseDof(:,:), irfDof(:,:)
  REAL(prec), ALLOCATABLE :: irf(:,:)
  REAL(prec), ALLOCATABLE :: Lmatrix(:,:)
  REAL(prec), ALLOCATABLE :: v(:,:), vDof(:)
  REAL(prec), ALLOCATABLE :: Av(:,:), AvDof(:), AvCheck(:,:)
  REAL(prec), ALLOCATABLE :: uniqueEval(:)
  INTEGER, ALLOCATABLE :: uniqueIndices(:)
  REAL(prec) :: absmaxDiff
  INTEGER :: row, i, j, nUnique
  INTEGER :: N, lwork, liwork, info
  REAL(prec), ALLOCATABLE :: w(:) ! Eigenvalues
  REAL(prec), ALLOCATABLE :: work(:)
  INTEGER, ALLOCATABLE :: iwork(:)
  REAL(prec), ALLOCATABLE :: modalCoeffs(:)
  CHARACTER(8) :: evec
  REAL(prec) :: vDotL, lNorm
  INTEGER :: fUnit
  CHARACTER(LEN=255) :: MDPath
  
     
  ! Create a stencil
  CALL modelStencil % Build()
  
  CALL get_environment_variable("SLS_MDPATH", MDPath)

  ! Create a mesh
  CALL modelMesh % Load(MDPath)
  
  CALL modelMesh % GenerateMaskFromhFac( x1, x2, y1, y2, klev )
  
  CALL modelMesh % WriteTecplot( modelMesh % tracerMask, 'mask', 'mask.tec' )

 ! Associate the mesh with the Generator
  CALL laplacian % AssociateMesh( modelMesh )
  
  ! Create an adjacency graph using the mesh and the stencil
  CALL diagnosisGraph % Build( modelMesh, modelStencil % firstOrder )
  CALL overlapGraph % Build( modelMesh, modelStencil % secondOrder )
  
  ! Color the graph
  CALL overlapGraph % GreedyColoring()
  
  ! Allocate space for the impulse and impulse response functions
  ALLOCATE(impulseDof(1:modelMesh % nDOF, overlapGraph % nColors), &
           irfDof(1:modelMesh % nDOF, overlapGraph % nColors), &
           irf(1:modelMesh % nX, 1:modelMesh % nY), &
           LMatrix(1:modelMesh % nDOF,1:modelMesh % nDOF), &
           v(1:modelMesh % nX,1:modelMesh % nY), &
           vDof(1:modelMesh % nDOF), &
           Av(1:modelMesh % nX,1:modelMesh % nY), &
           AvDof(1:modelMesh % nDOF), &
           modalCoeffs(1:modelMesh % nDOF), &
           uniqueEval(1:modelMesh % nDOF), &
           uniqueIndices(1:modelMesh % nDOF))
  
 ! ! Initialize all values to 0
  impulseDof = 0.0_prec
  irfDof = 0.0_prec
  irf = 0.0_prec
  
  ! Create the impulse fields from the colored graph
  impulseDof = overlapGraph % ImpulseFields()
  
  DO i = 1, overlapGraph % nColors
    ! Convert the impulse fields to ij format from dof format
    v = modelMesh % gridMap(impulseDof(:,i))
    ! Calculate the impulse response
    irf = laplacian % SLOperator( v )  
    ! Convert the impulse response field to dof format from ij format
    irfDof(:,i) = modelMesh % FlatMap( irf )
  ENDDO
  
  ! Create a matrix from the graph and the irf
  LMatrix = overlapGraph % DenseMatrix(diagnosisGraph, irfDof)

 ! Get Eigenvalues and Eigenvectors
  N = modelMesh % nDOF
  lwork = 1+6*N + 2*N*N
  liwork = 3+5*N
  ALLOCATE( w(1:N), work(1:lwork), iwork(1:liwork) )
  PRINT*, 'Finding eigenvalues'
  
  IF( prec == real32 )THEN
    CALL SSYEVD( 'V', 'U', N, &
                Lmatrix, N, &
                w, work, lwork, &
                iwork, liwork, info )
  ELSE
    CALL DSYEVD( 'V', 'U', N, &
                Lmatrix, N, &
                w, work, lwork, &
                iwork, liwork, info )
  ENDIF
  
  IF( info == 0 )THEN
    PRINT*,'Eigenvalues + Eigenvectors found!'
    
    DO i = 1, maxEvalWrite
      WRITE(evec,'(I8.8)') i
      j = modelMesh % nDOF + 1 - i
      v = modelMesh % gridMap(Lmatrix(:,j))
      CALL modelMesh % WriteTecplot(v, 'evec', 'evec.'//evec//'.tec' )
    ENDDO
    
    ! Write the eigenvalues
    OPEN( UNIT=NEWUNIT(fUnit), &
      FILE= 'evalues.curve', &
      FORM='formatted', &
      STATUS='replace')
    DO i = 1, modelMesh % nDOF
      j = modelMesh % nDOF + 1 - i
      WRITE(fUnit,*) i, ABS(w(j))
    ENDDO
    CLOSE(UNIT=fUnit)
    
    ! Find degenerate eigenvalues
    N = modelMesh % nDOF
    j = 1
    uniqueEval = 0.0_prec
    uniqueEval(1) = w(1)
    uniqueIndices(1) = 1
    DO i = 2, N
      IF( ANY( uniqueEval == w(i) ) ) THEN
         uniqueIndices(i) = j
        CYCLE
      ENDIF
      j = j + 1
      uniqueIndices(i) = j
      uniqueEval(j) = w(i)
    ENDDO
    nUnique = j
  
    PRINT*, 'Number of unique eigenvalues : ', j
    
    ! Project example data onto the eigenvectors
    
 !   ! > Example - create a sinusoid with known wavelength
 !   DO j = 1, modelMesh % nY
 !     DO i = 1, modelMesh % nY
 !       v(i,j) = sin( 2.0*pi*10.0*modelMesh % x(i,j) )*sin( 2.0*pi*10.0*modelMesh % y(i,j) )
 !     ENDDO
 !   ENDDO
 !   
 !   vDof = modelMesh % FlatMap( v )
 !   
 !   modalCoeffs = 0.0_prec
 !   DO j = 1, N
 !     i = uniqueIndices(j)
 !     vDotL = DOT_PRODUCT(vDOF, Lmatrix(:,j))
 !     Lnorm = SQRT(DOT_PRODUCT(Lmatrix(:,j),Lmatrix(:,j)))
 !     modalCoeffs(i) = modalCoeffs(i) + (vDotL/Lnorm)**2
 !   ENDDO
 !       
 !   ! Write the  modal coefficients
 !   OPEN( UNIT=NEWUNIT(fUnit), &
 !     FILE= 'sin-10-modal.curve', &
 !     FORM='formatted', &
 !     STATUS='replace')
 !   
 !   DO i = 1, nUnique
 !     WRITE(fUnit,*) uniqueEval(i), modalCoeffs(i)
 !   ENDDO
 !   CLOSE(UNIT=fUnit)
 !   
 !   ! > Example - create a gaussian with known halfwidth
 !   DO j = 1, modelMesh % nY
 !     DO i = 1, modelMesh % nY
 !       v(i,j) = exp( -( (modelMesh % x(i,j)-0.5_prec)**2 + (modelMesh % y(i,j)-0.5_prec)**2 )/(2.0_prec*0.01_prec) )
 !     ENDDO
 !   ENDDO
 !   
 !   vDof = modelMesh % FlatMap( v )
 !   
 !   modalCoeffs = 0.0_prec
 !   DO j = 1, N
 !     i = uniqueIndices(j)
 !     vDotL = DOT_PRODUCT(vDOF, Lmatrix(:,j))
 !     Lnorm = SQRT(DOT_PRODUCT(Lmatrix(:,j),Lmatrix(:,j)))
 !     modalCoeffs(i) = modalCoeffs(i) + (vDotL/Lnorm)**2
 !   ENDDO
 !       
 !   ! Write the  modal coefficients
 !   OPEN( UNIT=NEWUNIT(fUnit), &
 !     FILE= 'guassian-modal.curve', &
 !     FORM='formatted', &
 !     STATUS='replace')
 !   
 !   DO i = 1, nUnique
 !     WRITE(fUnit,*) uniqueEval(i), modalCoeffs(i)
 !   ENDDO
 !   CLOSE(UNIT=fUnit)
 !     
 !   
  ELSEIF( info < 0 )THEN
    PRINT*, 'Illegal value in argument : ',ABS(info)
  ELSE
    PRINT*, 'Algorithm failed!'
  ENDIF
 ! 
 ! ! Clean up memory
 ! DEALLOCATE(impulseDof, &
 !            irfDof, &
 !            impulse, &
 !            irf, &
 !            LMatrix, &
 !            v, Av, &
 !            vDof, AvDof, &
 !            AvCheck, &
 !            modalCoeffs, &
 !            uniqueEval, &
 !            uniqueIndices)
 ! 
 ! DEALLOCATE( w, work, iwork )

 ! CALL laplacian % DisassociateMesh( )
  CALL modelMesh % Free()
  CALL modelStencil % Free()
 ! CALL diagnosisGraph % Free()
 ! CALL overlapGraph % Free()
  

END PROGRAM IrregularGeometry
