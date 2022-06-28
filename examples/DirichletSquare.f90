PROGRAM DirichletSquare

USE SLSpectra_Precision
USE SLSpectra_Mesh
USE SLSpectra_Stencil
USE SLSpectra_AdjacencyGraph
USE SLSpectra_Generators

IMPLICIT NONE

  INTEGER, PARAMETER :: nX = 5
  INTEGER, PARAMETER :: nY = 5
  
  TYPE(Laplacian5Stencil) :: modelStencil
  TYPE(Mesh), TARGET :: modelMesh
  TYPE(Generator) :: laplacian
  TYPE(AdjacencyGraph) :: diagnosisGraph
  TYPE(AdjacencyGraph) :: overlapGraph
  REAL(prec), ALLOCATABLE :: impulseDof(:,:), irfDof(:,:)
  REAL(prec), ALLOCATABLE :: impulse(:,:), irf(:,:)
  REAL(prec), ALLOCATABLE :: Lmatrix(:,:)
  REAL(prec), ALLOCATABLE :: v(:,:), vDof(:)
  REAL(prec), ALLOCATABLE :: Av(:,:), AvDof(:), AvCheck(:,:)
  REAL(prec) :: absmaxDiff
  INTEGER :: row, i, j 
  
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
           LMatrix(1:modelMesh % nDOF,1:modelMesh % nDOF), &
           v(1:nX,1:nY), &
           vDof(1:modelMesh % nDOF), &
           Av(1:nX,1:nY), &
           AvDof(1:modelMesh % nDOF), &
           AvCheck(1:nX,1:nY))
  
  ! Initialize all values to 0
  impulseDof = 0.0_prec
  irfDof = 0.0_prec
  impulse = 0.0_prec
  irf = 0.0_prec
  
  ! Create the impulse fields from the colored graph
  impulseDof = overlapGraph % ImpulseFields()
  
  DO i = 1, overlapGraph % nColors
    ! Convert the impulse fields to ij format from dof format
    impulse = modelMesh % gridMap(impulseDof(:,i))
    
    ! Calculate the impulse response
    irf = laplacian % SLOperator( impulse )
   
    ! Convert the impulse response field to dof format from ij format
    irfDof(:,i) = modelMesh % FlatMap( irf )
    
  ENDDO
  
  ! Create a matrix from the graph and the irf
  LMatrix = overlapGraph % DenseMatrix(diagnosisGraph, irfDof)

  ! Create a step function where v = 1 in the upper right corner of the domain
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
  
  ! To do : Calculate matrix action using SLOperator
  Av = laplacian % SLOperator( v )
  
  ! To do : Calculate matrix action using diagnosed matrix
  DO row = 1, modelMesh % nDOF
    AvDof(row) = 0.0_prec
    DO i = 1, modelMesh % nDOF
        AvDof(row) = AvDof(row) + LMatrix(row,i)*vDof(i)
    ENDDO
  ENDDO
      
  ! To do : Compare results and verify consistency
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
  
  ! Clean up memory
  DEALLOCATE(impulseDof, &
             irfDof, &
             impulse, &
             irf, &
             LMatrix, &
             v, Av, &
             vDof, AvDof, &
             AvCheck)
  CALL laplacian % DisassociateMesh( )
  CALL modelMesh % Free()
  CALL modelStencil % Free()
  CALL diagnosisGraph % Free()
  CALL overlapGraph % Free()
  

END PROGRAM DirichletSquare
