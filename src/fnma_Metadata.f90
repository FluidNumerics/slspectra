! fnma_Metadata.f90
! Adapted from SELF_Metadata.F90
!
! Copyright 2020-2023 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

MODULE fnma_Metadata

  INTEGER,PARAMETER,PUBLIC :: fnma_MTD_NameLength = 250
  INTEGER,PARAMETER,PUBLIC :: fnma_MTD_DescriptionLength = 1000
  INTEGER,PARAMETER,PUBLIC :: fnma_MTD_UnitsLength = 20

  ! A class for storing metadata information, intended for file IO
  TYPE fnmaMetadata
    CHARACTER(fnma_MTD_NameLength) :: name
    CHARACTER(fnma_MTD_DescriptionLength) :: description
    CHARACTER(fnma_MTD_UnitsLength) :: units

  CONTAINS

    PROCEDURE,PUBLIC :: SetName => SetName_fnmaMetadata
    PROCEDURE,PUBLIC :: SetDescription => SetDescription_fnmaMetadata
    PROCEDURE,PUBLIC :: SetUnits => SetUnits_fnmaMetadata

  END TYPE fnmaMetadata

CONTAINS

  SUBROUTINE SetName_fnmaMetadata(mtd,name)
    IMPLICIT NONE
    CLASS(fnmaMetadata),INTENT(inout) :: mtd
    CHARACTER(*),INTENT(in) :: name

    mtd % name = name

  END SUBROUTINE SetName_fnmaMetadata

  SUBROUTINE SetDescription_fnmaMetadata(mtd,description)
    IMPLICIT NONE
    CLASS(fnmaMetadata),INTENT(inout) :: mtd
    CHARACTER(*),INTENT(in) :: description

    mtd % description = description

  END SUBROUTINE SetDescription_fnmaMetadata

  SUBROUTINE SetUnits_fnmaMetadata(mtd,units)
    IMPLICIT NONE
    CLASS(fnmaMetadata),INTENT(inout) :: mtd
    CHARACTER(*),INTENT(in) :: units

    mtd % units = units

  END SUBROUTINE SetUnits_fnmaMetadata

END MODULE fnma_Metadata
