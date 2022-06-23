PROGRAM DirichletSquare

USE SLSpectra_Precision
USE SLSpectra_Mesh
USE SLSpectra_Stencil
USE SLSpectra_AdjacencyGraph
USE SLSpectra_Generators

IMPLICIT NONE

  INTEGER, PARAMETER :: nX = 5
  INTEGER, PARAMETER :: nY = 5
  
  TYPE(Laplacian5OStencil) :: modelStencil
  TYPE(Mesh), TARGET :: modelMesh
  TYPE(Generator) :: laplacian
  TYPE(AdjacencyGraph) :: graph
  REAL(prec), ALLOCATABLE :: impulseDof(:,:), irfDof(:,:)
  REAL(prec), ALLOCATABLE :: impulse(:,:), irf(:,:)
  REAL(prec), ALLOCATABLE :: Lmatrix(:,:)
  INTEGER :: i, j 
  
  ! Create a stencil
  CALL modelStencil % Build()
  
  ! Create a mesh
  CALL modelMesh % Build(nX,nY)
  CALL modelMesh % ConstructDirichletCube()
  
 ! Associate the mesh with the Generator
  CALL laplacian % AssociateMesh( modelMesh )
  
  ! Create an adjacency graph using the mesh and the stencil
  CALL graph % Build( modelMesh, modelStencil )
  
  ! Color the graph
  CALL graph % GreedyColoring()
  
  ! Allocate space for the impulse and impulse response functions
  ALLOCATE(impulseDof(1:modelMesh % nDOF, graph % nColors), &
           irfDof(1:modelMesh % nDOF, graph % nColors), &
           impulse(1:nX, 1:nY), &
           irf(1:nX, 1:nY), &
           LMatrix(1:modelMesh % nDOF,1:modelMesh % nDOF))
  
  ! Initialize all values to 0
  impulseDof = 0.0_prec
  irfDof = 0.0_prec
  impulse = 0.0_prec
  irf = 0.0_prec
  
  ! Create the impulse fields from the colored graph
  impulseDof = graph % ImpulseFields()
  
  DO i = 1, graph % nColors
    ! Convert the impulse fields to ij format from dof format
    impulse = modelMesh % gridMap(impulseDof(:,i))
    
    ! Calculate the impulse response
    irf = laplacian % SLOperator( impulse )
   
    ! Convert the impulse response field to dof format from ij format
    irfDof(:,i) = modelMesh % FlatMap( irf )
    
  ENDDO
  

  ! Create a matrix from the graph and the irf
  LMatrix = graph % DenseMatrix(irfDof)
  DO j = 1, modelMesh % nDOF
    DO i = 1, modelMesh % nDOF
      PRINT*, LMatrix(i,j) - LMatrix(j,i)
    ENDDO
  ENDDO
  ! To do : Calculate matrix action using SLOperator
  ! To do : Calculate matrix action using diagnosed matrix
  ! To do : Compare results and verify consistency
  
  ! Clean up memory
  DEALLOCATE(impulseDof, &
             irfDof, &
             impulse, &
             irf, &
             LMatrix)
  CALL laplacian % DisassociateMesh( )
  CALL modelMesh % Free()
  CALL modelStencil % Free()
  CALL graph % Free()
  

END PROGRAM DirichletSquare
