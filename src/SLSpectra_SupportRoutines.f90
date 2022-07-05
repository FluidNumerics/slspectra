MODULE SLSpectra_SupportRoutines

USE ISO_FORTRAN_ENV

IMPLICIT NONE

  INTERFACE AlmostEqual
    MODULE PROCEDURE AlmostEqual_r64, &
                     AlmostEqual_r32
  END INTERFACE AlmostEqual

CONTAINS

  FUNCTION AlmostEqual_r64(a,b) RESULT(AisB)
  
    IMPLICIT NONE
    REAL(real64) :: a,b
    LOGICAL :: AisB
  
    IF (a == 0.0_real64 .OR. b == 0.0_real64) THEN
      IF (ABS(a - b) <= EPSILON(1.0_real64)) THEN
        AisB = .TRUE.
      ELSE
        AisB = .FALSE.
      END IF
    ELSE
      IF ((abs(a - b) <= EPSILON(1.0_real64)*abs(a)) .OR. (abs(a - b) <= EPSILON(1.0_real64)*abs(b))) THEN
        AisB = .TRUE.
      ELSE
        AisB = .FALSE.
      END IF
    END IF
  
  END FUNCTION AlmostEqual_r64
  
  FUNCTION AlmostEqual_r32(a,b) RESULT(AisB)
  
    IMPLICIT NONE
    REAL(real32) :: a,b
    LOGICAL :: AisB
  
    IF (a == 0.0_real32 .OR. b == 0.0_real32) THEN
      IF (ABS(a - b) <= EPSILON(1.0_real32)) THEN
        AisB = .TRUE.
      ELSE
        AisB = .FALSE.
      END IF
    ELSE
      IF ((abs(a - b) <= EPSILON(1.0_real32)*abs(a)) .OR. (abs(a - b) <= EPSILON(1.0_real32)*abs(b))) THEN
        AisB = .TRUE.
      ELSE
        AisB = .FALSE.
      END IF
    END IF
  
  END FUNCTION AlmostEqual_r32

END MODULE SLSpectra_SupportRoutines
