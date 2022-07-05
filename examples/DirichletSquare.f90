PROGRAM DirichletSquare

USE ISO_FORTRAN_ENV
USE SLSpectra_Precision
USE SLSpectra_SupportRoutines
USE SLSpectra_Mesh
USE SLSpectra_Stencil
USE SLSpectra_AdjacencyGraph
USE SLSpectra_Generators

IMPLICIT NONE

  INTEGER, PARAMETER :: nX = 100
  INTEGER, PARAMETER :: nY = 100
  REAL(prec), PARAMETER :: pi = 4.0_prec*atan(1.0_prec)
  TYPE(Laplacian5Stencil) :: modelStencil
  TYPE(Mesh), TARGET :: modelMesh
  TYPE(Generator) :: laplacian
  TYPE(AdjacencyGraph) :: diagnosisGraph
  TYPE(AdjacencyGraph) :: overlapGraph
  REAL(prec), ALLOCATABLE :: impulseDof(:,:), irfDof(:,:)
  REAL(prec), ALLOCATABLE :: impulse(:,:), irf(:,:)
  REAL(prec), ALLOCATABLE :: colors(:,:)
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
  
     
  ! Create a stencil
  CALL modelStencil % Build()
  
  ! Create a mesh
  CALL modelMesh % Build(nX,nY)
  CALL modelMesh % ConstructDirichletCube()
  
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
           impulse(1:nX, 1:nY), &
           irf(1:nX, 1:nY), &
           colors(1:nX, 1:nY), &
           LMatrix(1:modelMesh % nDOF,1:modelMesh % nDOF), &
           v(1:nX,1:nY), &
           vDof(1:modelMesh % nDOF), &
           Av(1:nX,1:nY), &
           AvDof(1:modelMesh % nDOF), &
           AvCheck(1:nX,1:nY),&
           modalCoeffs(1:modelMesh % nDOF), &
           uniqueEval(1:modelMesh % nDOF), &
           uniqueIndices(1:modelMesh % nDOF))
  
  ! Initialize all values to 0
  impulseDof = 0.0_prec
  irfDof = 0.0_prec
  impulse = 0.0_prec
  irf = 0.0_prec
  
  ! Create the impulse fields from the colored graph
  impulseDof = overlapGraph % ImpulseFields()
  
  colors = 0.0_prec
  DO i = 1, overlapGraph % nColors
    ! Convert the impulse fields to ij format from dof format
    impulse = modelMesh % gridMap(impulseDof(:,i))
    colors = colors+impulse*REAL(i,prec)
    
    WRITE(evec,'(I8.8)') i
    CALL modelMesh % WriteTecplot( impulse, 'impulse', 'impulsefield.'//evec//'.tec' )
    
    ! Calculate the impulse response
    irf = laplacian % SLOperator( impulse )
   
    ! Convert the impulse response field to dof format from ij format
    irfDof(:,i) = modelMesh % FlatMap( irf )
    
  ENDDO
  CALL modelMesh % WriteTecplot( colors, 'colors', 'colors.tec' )
  
  ! Create a matrix from the graph and the irf
  LMatrix = overlapGraph % DenseMatrix(diagnosisGraph, irfDof)

  ! Create a step function where v = 1 in the upper right corner of the domain
  PRINT*, 'Verifying generated matrix'
  v = 0.0_prec
  DO j = 1, nY
    DO i = 1, nX
      IF( modelMesh % x(i,j) > 0.5_prec .AND. &
          modelMesh % y(i,j) > 0.5_prec .AND. &
          modelMesh % tracerMask(i,j) == SLSpectra_wetValue )THEN
        v(i,j) = 1.0_prec
      ENDIF
    ENDDO
  ENDDO
  
  ! Create an equivalent array for "v" in DOF format
  vDof = modelMesh % flatMap( v )
  
  ! Calculate matrix action using SLOperator
  Av = laplacian % SLOperator( v )
  
  ! Calculate matrix action using diagnosed matrix
  DO row = 1, modelMesh % nDOF
    AvDof(row) = 0.0_prec
    DO i = 1, modelMesh % nDOF
        AvDof(row) = AvDof(row) + LMatrix(row,i)*vDof(i)
    ENDDO
  ENDDO
      
  ! Compare results and verify consistency
  AvCheck = modelMesh % gridMap(AvDof)
  absmaxDiff = 0.0_prec
  DO j = 1, nY
    DO i = 1, nX
      IF( modelMesh % tracerMask(i,j) == SLSpectra_wetValue )THEN
        absmaxDiff = MAX(absmaxDiff, ABS(AvCheck(i,j) - Av(i,j)))
      ENDIF
    ENDDO
  ENDDO
  
  PRINT*, 'Absolute max difference : ', absmaxDiff
  
  ! To do : Get Eigenvalues and Eigenvectors
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
    
    ! > Example - create a sinusoid with known wavelength
    DO j = 1, modelMesh % nY
      DO i = 1, modelMesh % nY
        v(i,j) = sin( 2.0*pi*10.0*modelMesh % x(i,j) )*sin( 2.0*pi*10.0*modelMesh % y(i,j) )
      ENDDO
    ENDDO
    
    vDof = modelMesh % FlatMap( v )
    
    modalCoeffs = 0.0_prec
    DO j = 1, N
      i = uniqueIndices(j)
      vDotL = DOT_PRODUCT(vDOF, Lmatrix(:,j))
      Lnorm = SQRT(DOT_PRODUCT(Lmatrix(:,j),Lmatrix(:,j)))
      modalCoeffs(i) = modalCoeffs(i) + (vDotL/Lnorm)**2
    ENDDO
        
    ! Write the  modal coefficients
    OPEN( UNIT=NEWUNIT(fUnit), &
      FILE= 'sin-10-modal.curve', &
      FORM='formatted', &
      STATUS='replace')
    
    DO i = 1, nUnique
      WRITE(fUnit,*) uniqueEval(i), modalCoeffs(i)
    ENDDO
    CLOSE(UNIT=fUnit)
    
    ! > Example - create a gaussian with known halfwidth
    DO j = 1, modelMesh % nY
      DO i = 1, modelMesh % nY
        v(i,j) = exp( -( (modelMesh % x(i,j)-0.5_prec)**2 + (modelMesh % y(i,j)-0.5_prec)**2 )/(2.0_prec*100.0_prec) )
      ENDDO
    ENDDO
    
    vDof = modelMesh % FlatMap( v )
    
    modalCoeffs = 0.0_prec
    DO j = 1, N
      i = uniqueIndices(j)
      vDotL = DOT_PRODUCT(vDOF, Lmatrix(:,j))
      Lnorm = SQRT(DOT_PRODUCT(Lmatrix(:,j),Lmatrix(:,j)))
      modalCoeffs(i) = modalCoeffs(i) + (vDotL/Lnorm)**2
    ENDDO
        
    ! Write the  modal coefficients
    OPEN( UNIT=NEWUNIT(fUnit), &
      FILE= 'guassian-modal.curve', &
      FORM='formatted', &
      STATUS='replace')
    
    DO i = 1, nUnique
      WRITE(fUnit,*) uniqueEval(i), modalCoeffs(i)
    ENDDO
    CLOSE(UNIT=fUnit)
      
    
  ELSEIF( info < 0 )THEN
    PRINT*, 'Illegal value in argument : ',ABS(info)
  ELSE
    PRINT*, 'Algorithm failed!'
  ENDIF
  
  ! Clean up memory
  DEALLOCATE(impulseDof, &
             irfDof, &
             impulse, &
             irf, &
             LMatrix, &
             v, Av, &
             vDof, AvDof, &
             AvCheck, &
             modalCoeffs, &
             uniqueEval, &
             uniqueIndices)
  
  DEALLOCATE( w, work, iwork )

  CALL laplacian % DisassociateMesh( )
  CALL modelMesh % Free()
  CALL modelStencil % Free()
  CALL diagnosisGraph % Free()
  CALL overlapGraph % Free()
  

END PROGRAM DirichletSquare
