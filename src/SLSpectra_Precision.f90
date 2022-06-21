MODULE SLSpectra_Precision

USE ISO_FORTRAN_ENV

IMPLICIT NONE

#ifdef DOUBLE_PRECISION
  INTEGER,PARAMETER :: prec = real64
#else
  INTEGER,PARAMETER :: prec = real32
#endif

END MODULE SLSpectra_Precision
