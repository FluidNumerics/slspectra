MODULE SLSpectra_Generators
!! This module is used to specify routines that implement the action of Sturm-Liouville operators
!! The name "Generator" refers to the fact that the Sturm-Liouville boundary value problems
!! provide a means to generate a basis

USE SLSpectra_Precision
USE SLSpectra_Mesh

IMPLICIT NONE


  TYPE Generator
    TYPE(Mesh), POINTER :: modelMesh => NULL()
    PROCEDURE(SLSOp), POINTER :: SLOperator => SLOperator_Default
    CONTAINS
    
      PROCEDURE :: AssociateMesh
      
      PROCEDURE :: DisassociateMesh
      PROCEDURE :: SLOperator_Default
      PROCEDURE :: SLOperator_Neumann
      
      PROCEDURE :: ApplyNeumannCondition
      
  END TYPE Generator
  
  INTERFACE
    FUNCTION SLSOp( this, v ) RESULT(Av)
      USE SLSpectra_Precision, ONLY : prec
      IMPORT Generator
      IMPLICIT NONE
      CLASS(Generator) :: this
      REAL(prec) :: v(1:this % modelMesh % nX, 1:this % modelMesh % nY)
      REAL(prec) :: Av(1:this % modelMesh % nX, 1:this % modelMesh % nY)
    END FUNCTION SLSOp
  END INTERFACE 
CONTAINS

  SUBROUTINE AssociateMesh( this, modelMesh )
    IMPLICIT NONE
    CLASS(Generator), INTENT(inout) :: this
    TYPE(Mesh), INTENT(in), TARGET :: modelMesh
    
      this % modelMesh => modelMesh
      
  END SUBROUTINE AssociateMesh
  
  SUBROUTINE DisassociateMesh( this )
    IMPLICIT NONE
    CLASS(Generator), INTENT(inout) :: this
    
      this % modelMesh => NULL()
      
  END SUBROUTINE DisassociateMesh
  
  SUBROUTINE ApplyNeumannCondition( this, v )
  !! 
    IMPLICIT NONE
    CLASS(Generator), INTENT(in) :: this
    REAL(prec), INTENT(inout) :: v(1:this % modelMesh % nX, 1:this % modelMesh % nY)
    ! Local
    REAL(prec) :: wetV
    INTEGER :: i, j
    
    DO j = 2, this % modelMesh % nY-1
        DO i = 2, this % modelMesh % nX-1
          
          IF ( this % modelMesh % tracerMask(i,j) == SLSpectra_wetValue ) THEN
          
            wetV = v(i,j)
                        
            IF( this % modelMesh % tracerMask(i-1,j) == SLSpectra_dryValue )THEN
              v(i-1,j) = wetV
            ENDIF
            
            IF( this % modelMesh % tracerMask(i+1,j) == SLSpectra_dryValue )THEN
              v(i+1,j) = wetV
            ENDIF
            
            IF( this % modelMesh % tracerMask(i,j-1) == SLSpectra_dryValue )THEN
              v(i,j-1) = wetV
            ENDIF
            
            IF( this % modelMesh % tracerMask(i,j+1) == SLSpectra_dryValue )THEN
              v(i,j+1) = wetV
            ENDIF
              
            
          ENDIF
          
        ENDDO
      ENDDO
    
  END SUBROUTINE ApplyNeumannCondition

  FUNCTION SLOperator_Default( this, v ) RESULT(Av)
  !! The default SLOperator takes in a 2-D array "v" (in ij format) and returns
  !! the action of the default Sturm-Liouville operator in "Av" (in ij format)
  !! The default Sturm-Liouville operator is a 2nd order finite difference 
  !! approximation for a laplacian on a uniform grid
    IMPLICIT NONE
    CLASS(Generator) :: this
    REAL(prec) :: v(1:this % modelMesh % nX, 1:this % modelMesh % nY)
    REAL(prec) :: Av(1:this % modelMesh % nX, 1:this % modelMesh % nY)
    ! Local
    INTEGER :: i, j
    REAL(prec) :: dv2dx2, dv2dy2
    
      
      Av = 0.0_prec
      
      DO j = 1, this % modelMesh % nY
        DO i = 1, this % modelMesh % nX
           
          IF ( this % modelMesh % tracerMask(i,j) == SLSpectra_wetValue ) THEN
          
            dv2dx2 = ( v(i+1,j) - 2.0_prec*v(i,j) + v(i-1,j) )/(this % modelMesh % dx**2)
            dv2dy2 = ( v(i,j+1) - 2.0_prec*v(i,j) + v(i,j-1) )/(this % modelMesh % dy**2)
            Av(i,j) = dv2dx2 + dv2dy2 
            
          ENDIF
          
        ENDDO
      ENDDO
      
  END FUNCTION SLOperator_Default
  
  FUNCTION SLOperator_Neumann( this, v ) RESULT(Av)
  !! The default SLOperator takes in a 2-D array "v" (in ij format) and returns
  !! the action of the default Sturm-Liouville operator in "Av" (in ij format)
  !! The default Sturm-Liouville operator is a 2nd order finite difference 
  !! approximation for a laplacian on a uniform grid
    IMPLICIT NONE
    CLASS(Generator) :: this
    REAL(prec) :: v(1:this % modelMesh % nX, 1:this % modelMesh % nY)
    REAL(prec) :: Av(1:this % modelMesh % nX, 1:this % modelMesh % nY)
    ! Local
    INTEGER :: i, j
    REAL(prec) :: dv2dx2, dv2dy2
    REAL(prec) :: fWest, fEast, fSouth, fNorth
    
      
      Av = 0.0_prec
      
      DO j = 1, this % modelMesh % nY
        DO i = 1, this % modelMesh % nX
           
          IF ( this % modelMesh % tracerMask(i,j) == SLSpectra_wetValue ) THEN
          
            ! West
            IF( this % modelMesh % tracerMask(i-1,j) == SLSpectra_dryValue ) THEN
              fWest = 0.0_prec
            ELSE
              fWest = (v(i,j) - v(i-1,j))/this % modelMesh % dx
            ENDIF
            
            ! East
            IF( this % modelMesh % tracerMask(i+1,j) == SLSpectra_dryValue ) THEN
              fEast = 0.0_prec
            ELSE
              fEast = (v(i+1,j) - v(i,j))/this % modelMesh % dx
            ENDIF
            
            ! South
            IF( this % modelMesh % tracerMask(i,j-1) == SLSpectra_dryValue ) THEN
              fSouth = 0.0_prec
            ELSE
              fSouth = (v(i,j) - v(i,j-1))/this % modelMesh % dy
            ENDIF
            
            ! North
            IF( this % modelMesh % tracerMask(i,j+1) == SLSpectra_dryValue ) THEN
              fNorth = 0.0_prec
            ELSE
              fNorth = (v(i,j+1) - v(i,j))/this % modelMesh % dy
            ENDIF
            
            dv2dx2 = ( fEast - fWest )/(this % modelMesh % dx)
            dv2dy2 = ( fNorth - fSouth )/(this % modelMesh % dy)
            Av(i,j) = dv2dx2 + dv2dy2 
            
          ENDIF
          
        ENDDO
      ENDDO
      
  END FUNCTION SLOperator_Neumann
  
  !SUBROUTINE Lanczos( this, nEval, v )
  !! This method applies the Lancsoz algorithm to obtain the 
  
END MODULE SLSpectra_Generators
