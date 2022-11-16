MODULE SLSpectra_Model

USE ISO_FORTRAN_ENV
USE SLSpectra_Precision
USE SLSpectra_SupportRoutines
USE SLSpectra_Mesh
USE SLSpectra_Stencil
USE SLSpectra_AdjacencyGraph
USE SLSpectra_Generators

IMPLICIT NONE

  TYPE Model
    INTEGER :: nDOF
    INTEGER :: nx, ny ! Convenience - copied from slGenerator % modelMesh attributes
    TYPE(AdjacencyGraph), POINTER :: diagGraph
    TYPE(AdjacencyGraph), POINTER :: overlapGraph
    TYPE(Generator), POINTER :: slGenerator !** Mesh inherited from generator
    REAL(prec), ALLOCATABLE :: eValues(:)
    REAL(prec), ALLOCATABLE :: eVectors(:,:)
    INTEGER, ALLOCATABLE :: uniqueEvalIndices(:)

    CONTAINS

    PROCEDURE :: Init => Init_Model
    PROCEDURE :: Free => Free_Model

    PROCEDURE :: SetDiagGraph
    PROCEDURE :: SetOverlapGraph
    PROCEDURE :: SetGenerator

    PROCEDURE :: CalculateEigenModes

    PROCEDURE, PRIVATE :: Lapack_SYEVDWrapper
    PROCEDURE, PRIVATE :: SLOpMatrix

    PROCEDURE :: ModalProjection

    PROCEDURE :: WriteTecplot => WriteTecplot_Model

    PROCEDURE :: VerifyMutualOrthogonality

  END TYPE Model

  CONTAINS

  SUBROUTINE Init_Model( this, dGraph, oGraph, slGenerator )
    IMPLICIT NONE
    CLASS(Model), INTENT(out) :: this
    TYPE(AdjacencyGraph), INTENT(in), TARGET :: dGraph
    TYPE(AdjacencyGraph), INTENT(in), TARGET :: oGraph
    TYPE(Generator), INTENT(in), TARGET :: slGenerator
    INTEGER :: nDOF

      CALL this % SetDiagGraph( dGraph )
      CALL this % SetOverlapGraph( oGraph )
      CALL this % SetGenerator( slGenerator )

      nDOF = this % slGenerator % modelMesh % nDOF
      ALLOCATE( this % eValues(1:nDOF), &
                this % eVectors(1:nDOF,1:nDOF), &
                this % uniqueEvalIndices(1:nDOF) )
      this % nDOF = nDOF
      this % nx = this % slGenerator % modelMesh % nX
      this % ny = this % slGenerator % modelMesh % nY

      this % eValues = 0.0_prec
      this % eVectors = 0.0_prec

  END SUBROUTINE Init_Model

  SUBROUTINE Free_Model( this )
    IMPLICIT NONE
    CLASS(Model), INTENT(inout) :: this

      this % diagGraph => NULL()
      this % overlapGraph => NULL()
      this % slGenerator => NULL()

      DEALLOCATE( this % eValues, &
                  this % eVectors, &
                  this % uniqueEvalIndices )

  END SUBROUTINE Free_Model

  SUBROUTINE SetDiagGraph( this, graph )
    IMPLICIT NONE
    CLASS(Model), INTENT(inout) :: this
    TYPE(AdjacencyGraph), INTENT(in), TARGET :: graph
    
      this % diagGraph => graph
      
  END SUBROUTINE SetDiagGraph

  SUBROUTINE SetOverlapGraph( this, graph )
    IMPLICIT NONE
    CLASS(Model), INTENT(inout) :: this
    TYPE(AdjacencyGraph), INTENT(in), TARGET :: graph
    
      this % overlapGraph => graph
      
  END SUBROUTINE SetOverlapGraph

  SUBROUTINE SetGenerator( this, slGenerator )
    IMPLICIT NONE
    CLASS(Model), INTENT(inout) :: this
    TYPE(Generator), INTENT(in), TARGET :: slGenerator
    
      this % slGenerator => slGenerator
      
  END SUBROUTINE SetGenerator

  SUBROUTINE CalculateEigenModes( this )
    IMPLICIT NONE
    CLASS(Model), INTENT(inout) :: this
    ! Local 
    REAL(prec), ALLOCATABLE :: Lmatrix(:,:)
    INTEGER :: info

    ! Color the graph
    CALL this % overlapGraph % GreedyColoring()

    ! Allocate space for the impulse and impulse response functions
    ALLOCATE(LMatrix(1:this % nDOF,1:this % nDOF))
    
    ! Generate a Dense Matrix consistent with the Sturm-Liouville Operator
    LMatrix = this % SLOpMatrix( )

    ! Call *SYEVD from Lapack/OpenBLAS and return eigenvalues/eigenvectors
    ! with unique eigenvalues detected (degenerate modes are identified)
    CALL this % Lapack_SYEVDWrapper( LMatrix, info ) 

    DEALLOCATE( Lmatrix )

  END SUBROUTINE CalculateEigenModes

  SUBROUTINE Lapack_SYEVDWrapper( this, L, info) 
    IMPLICIT NONE
    CLASS(Model), INTENT(inout) :: this
    REAL(prec), INTENT(inout) :: L(1:this % nDOF,1:this % nDOF)
    INTEGER, INTENT(out) :: info
    ! Local
    INTEGER :: i, j, nUnique
    INTEGER :: N, lwork, liwork
    REAL(prec), ALLOCATABLE :: uniqueEval(:)
    INTEGER, ALLOCATABLE :: uniqueIndices(:)
    REAL(prec), ALLOCATABLE :: w(:) ! Eigenvalues
    REAL(prec), ALLOCATABLE :: work(:)
    INTEGER, ALLOCATABLE :: iwork(:)
    REAL(prec) :: evec(1:this % nx, 1:this % ny)
    REAL(preC) :: LDotL

    ALLOCATE(&
      uniqueEval(1:this % nDOF), &
      uniqueIndices(1:this % nDOF))
    
    N = this % slGenerator % modelMesh % nDOF
    lwork = 1+6*N + 2*N*N
    liwork = 3+5*N

    ALLOCATE( w(1:N), work(1:lwork), iwork(1:liwork) )

    PRINT*, 'Finding eigenvalues'
    
    IF( prec == real32 )THEN
      CALL SSYEVD( 'V', 'U', N, &
                  L, N, &
                  w, work, lwork, &
                  iwork, liwork, info )
    ELSE
      CALL DSYEVD( 'V', 'U', N, &
                  L, N, &
                  w, work, lwork, &
                  iwork, liwork, info )
    ENDIF
    
    IF( info == 0 )THEN
      PRINT*,'Eigenvalues + Eigenvectors found!'
      
      ! Find degenerate eigenvalues
      N = this % slGenerator % modelMesh % nDOF
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
    
      ! Save the eigenvalues
      this % uniqueEvalIndices = uniqueIndices
      this % eValues = w
      
      PRINT*, "Normalizing eigenvectors"
      ! Normalize the eigenvectors wrt inner product
      DO i = 1, N
        evec = this % slGenerator % modelMesh % GridMap(L(:,i))
        LDotL = this % slGenerator % SLInnerProduct( evec, evec )
        this % eVectors(:,i) = L(:,i)/SQRT(LDotL)
      ENDDO


    ELSEIF( info < 0 )THEN
      PRINT*, 'Illegal value in argument : ',ABS(info)
    ELSE
      PRINT*, 'Algorithm failed!'
    ENDIF


    DEALLOCATE(&
      uniqueEval, &
      uniqueIndices, &
      w, work, iwork)
    
  END SUBROUTINE Lapack_SYEVDWrapper

  FUNCTION SLOpMatrix( this ) RESULT( L )
    IMPLICIT NONE
    CLASS(Model) :: this
    REAL(prec) :: L(1:this % nDOF,1:this % nDOF)
    ! Local
    REAL(prec), ALLOCATABLE :: impulseDof(:,:), irfDof(:,:)
    REAL(prec), ALLOCATABLE :: irf(:,:), v(:,:)
    INTEGER :: i

    ! Allocate space for the impulse and impulse response functions
    ALLOCATE(&
      impulseDof(1:this % nDOF, this % overlapGraph % nColors), &
      irfDof(1:this % nDOF, this % overlapGraph % nColors), &
      irf(1:this % slGenerator % modelMesh % nX, 1:this % slGenerator % modelMesh % nY), &
      v(1:this % slGenerator % modelMesh % nX,1:this % slGenerator % modelMesh % nY) )

 !   ! Initialize all values to 0
    impulseDof = 0.0_prec
    irfDof = 0.0_prec
    irf = 0.0_prec
    
    ! Create the impulse fields from the colored graph
    impulseDof = this % overlapGraph % ImpulseFields()
    
    DO i = 1, this % overlapGraph % nColors
      ! Convert the impulse fields to ij format from dof format
      v = this % slGenerator % modelMesh % gridMap(impulseDof(:,i))
      ! Calculate the impulse response
      irf = this % slGenerator % SLOperator( v )  
      ! Convert the impulse response field to dof format from ij format
      irfDof(:,i) = this % slGenerator % modelMesh % FlatMap( irf )
    ENDDO
    
    ! Create a matrix from the graph and the irf
    L = this % overlapGraph % DenseMatrix(this % diagGraph, irfDof)

    DEALLOCATE( impulseDof, irfDof, irf, v )

  END FUNCTION SLOpMatrix

  SUBROUTINE ModalProjection( this, vData, modalCoeffs, projection, residual )
  !!  Calculates the modal projection of vData onto the eigenvectors
  !!  Input 
  !!    this - Model
  !!    vData - 2-D array of data; size is consistent with model's mesh
  !!
  !!  Output
  !!    vHat - modal coefficients
  !!    residual - 2-D array showing the residual from the projection
    IMPLICIT NONE
    CLASS(Model), INTENT(in) :: this
    REAL(prec), INTENT(in) :: vData(1:this % nx, 1:this % ny)
    REAL(prec), INTENT(out) :: modalCoeffs(1:this % ndof)
    REAL(prec), INTENT(out) :: projection(1:this % nx, 1:this % ny)
    REAL(prec), INTENT(out) :: residual(1:this % nx, 1:this % ny)
    ! Local
    INTEGER :: i, j, k
    REAL(prec) :: evec(1:this % nx, 1:this % ny)
    REAL(prec) :: vhat, enorm
    
    PRINT*, "Calculating projection"
    modalCoeffs = 0.0_prec
    projection = 0.0_prec
    DO i = this % nDOF, 1, -1
      
      evec = this % slGenerator % modelMesh % GridMap(this % eVectors(:,i))

      vhat = this % slGenerator % SLInnerProduct( vData, evec )
      !enorm = this % slGenerator % SLInnerProduct( evec, evec )

      ! Save the coefficients
      modalCoeffs(i) = vhat

      ! Calculate the projection
      DO k = 1, this % ny
        DO j = 1, this % nx
          projection(j,k) = projection(j,k) + vhat*evec(j,k)*&
                this % slGenerator % modelMesh % tracerMask(j,k)
        ENDDO
      ENDDO

    ENDDO

    ! Calculate the residual
    DO k = 1, this % ny
      DO j = 1, this % nx
        residual(j,k) = (vData(j,k) - projection(j,k))*&
                this % slGenerator % modelMesh % tracerMask(j,k)
      ENDDO
    ENDDO

    PRINT*, "Absolute Max Residual", MAXVAL(ABS(residual))

  END SUBROUTINE ModalProjection

  SUBROUTINE VerifyMutualOrthogonality( this )
    IMPLICIT NONE
    CLASS(Model), INTENT(in) :: this
    ! Local
    INTEGER :: i, j, k
    REAL(prec) :: ei(1:this % nx, 1:this % ny)
    REAL(prec) :: ej(1:this % nx, 1:this % ny)
    REAL(prec) :: ref, eiDotej
    
    PRINT*, "Verifying mutual orthogonality"

    DO i = 1, this % nDOF
      ei = this % slGenerator % modelMesh % GridMap(this % eVectors(:,i))
      DO j = 1, i
        ej = this % slGenerator % modelMesh % GridMap(this % eVectors(:,j))
        eiDotej = this % slGenerator % SLInnerProduct( ei, ej )
        IF( i == j )THEN
          ref = 1.0_prec
        ELSE
          ref = 0.0_prec
        ENDIF

        IF( .NOT. AlmostEqual(eiDotej,ref) )THEN 
          PRINT*, "Possible Orthogonality Violation at (i,j) : ", i, j, eiDotej 
        ENDIF

      ENDDO
    ENDDO
      
  END SUBROUTINE VerifyMutualOrthogonality

  SUBROUTINE WriteTecplot_Model( this, nModes ) 
  !! Writes the eigenvectors to tecplot
    IMPLICIT NONE
    CLASS(Model), INTENT(in) :: this 
    INTEGER, INTENT(in) :: nModes
    ! Local
    REAL(prec):: v(1:this % nx,1:this % ny)
    INTEGER :: i, j
    CHARACTER(8) :: evec

    DO i = 1, nModes
      WRITE(evec,'(I8.8)') i
      j = this % slGenerator % modelMesh % nDOF + 1 - i
      v = this % slGenerator % modelMesh % gridMap(this % eVectors(:,j))
      CALL this % slGenerator % modelMesh % WriteTecplot(v, 'evec', 'evec.'//evec//'.tec' )
    ENDDO

  END SUBROUTINE WriteTecplot_Model

END MODULE SLSpectra_Model
