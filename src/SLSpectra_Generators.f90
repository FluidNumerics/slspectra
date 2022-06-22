MODULE SLSpectra_Generators
!! This module is used to specify routines that implement the action of Sturm-Liouville operators
!! The name "Generator" refers to the fact that the Sturm-Liouville boundary value problems
!! provide a means to generate a basis

USE SLSpectra_Precision
USE SLSpectra_Mesh

IMPLICIT NONE

  TYPE Generator
    TYPE(Mesh), POINTER :: modelMesh => NULL()
    
    CONTAINS
    
      PROCEDURE :: AssociateMesh
      PROCEDURE :: SLOperator => SLOperator_Default
      
  END TYPE Generator
  
CONTAINS

  SUBROUTINE AssociateMesh( this, modelMesh )
    IMPLICIT NONE
    CLASS(Generator), INTENT(inout) :: this
    TYPE(Mesh), INTENT(in), TARGET :: modelMesh
    
      this % modelMesh => modelMesh
      
  END SUBROUTINE AssociateMesh
  
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
  
END MODULE SLSpectra_Generators
