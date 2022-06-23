!
! Adapted from https://github.com/lanl/feots (POP_AdjacencyGraph_Class)
!
MODULE SLSpectra_AdjacencyGraph

USE SLSpectra_Stencil
USE SLSpectra_Mesh

IMPLICIT NONE

   TYPE AdjacencyGraph
      LOGICAL              :: colored
      INTEGER              :: nDOF
      INTEGER              :: maxValence
      INTEGER              :: nColors
      INTEGER, ALLOCATABLE :: valence(:)
      INTEGER, ALLOCATABLE :: color(:)
      INTEGER, ALLOCATABLE :: neighbors(:,:)

      CONTAINS

      GENERIC,PUBLIC :: Build => Build_AdjacencyGraph, &
                                 BuildFromMeshAndL5OStencil_AdjacencyGraph 
                                 
      PROCEDURE, PRIVATE :: Build_AdjacencyGraph
      PROCEDURE, PRIVATE :: BuildFromMeshAndL5OStencil_AdjacencyGraph
      
      PROCEDURE :: Free => Free_AdjacencyGraph      

      PROCEDURE :: GreedyColoring => GreedyColoring_AdjacencyGraph
      
      PROCEDURE :: ImpulseFields
      PROCEDURE :: DenseMatrix
    
!      PROCEDURE :: Write_HDF5
!      PROCEDURE :: Read_HDF5

   END TYPE AdjacencyGraph

CONTAINS

 SUBROUTINE Build_AdjacencyGraph( this, nDOF, maxValence )
 ! This subroutine allocates space for the adjacency graph and initializes
 ! all attributes to 0.
 ! INPUT: 
 !     nDOF - is the number of Degrees Of Freedom; ie, the number of wet
 !     gridpoints in the POP mesh.
 !     
 !     maxValence - the maximum number of neighbors that a grid point will
 !     have. This number depends on the stencil that is being used to construct
 !     the adjacency graph.
 ! ============================================================================ !
   IMPLICIT NONE
   CLASS( AdjacencyGraph ), INTENT(inout) :: this
   INTEGER, INTENT(in)                    :: nDOF
   INTEGER, INTENT(in)                    :: maxValence
   
      this % nDOF       = nDOF
      this % maxValence = maxValence
      this % nColors    = 0
      this % colored = .FALSE.

      ALLOCATE( this % valence(1:nDOF), & 
                this % color(1:nDOF), &
                this % neighbors(1:maxValence,1:nDOF) )

      this % valence    = 0
      this % color      = 0
      this % neighbors = 0

 END SUBROUTINE Build_AdjacencyGraph

 SUBROUTINE BuildFromMeshAndL5OStencil_AdjacencyGraph( this, modelMesh, relStencil )
 ! Given a Mesh and a 5-point Laplacian overlap stencil, this routine constructs an adjacency
 ! graph using only the wet points in the mesh.
    IMPLICIT NONE
    CLASS( AdjacencyGraph ), INTENT(inout) :: this
    TYPE( Mesh ), INTENT(in)  :: modelMesh
    TYPE( Laplacian5OStencil ), INTENT(in) :: relStencil
    ! Local
    INTEGER :: i, j,  e1, e2, m
    INTEGER :: ni, nj
    REAL(prec) :: t1, t2

      PRINT*, ' S/R BuildFromMeshAndStencil : Start! '

      CALL CPU_TIME( t1 )

      
      this % nDOF       = modelMesh % nDOF
      this % maxValence = relStencil % nPoints
      this % nColors    = 0
      this % colored = .FALSE.

      ALLOCATE( this % valence(1:modelMesh % nDOF), & 
                this % color(1:modelMesh % nDOF), &
                this % neighbors(1:this % maxValence,1:modelMesh % nDOF) )

      this % valence = 0
      this % color = 0
      this % neighbors = 0
      
       DO e1 = 1, modelMesh % nDOF

          i = modelMesh % DOFtoIJ(1,e1)
          j = modelMesh % DOFtoIJ(2,e1)

          DO m = 1, relStencil % nPoints ! Loop over the stencil points
   
             ! Find the i,j indices for this point in the stencil around e1.
             ni = i + relStencil % neighbors(1,m)
             nj = j + relStencil % neighbors(2,m)
 
             IF( ni > 0 .AND. ni <= modelMesh % nX .AND. &
                 nj > 0 .AND. nj <= modelMesh % nY ) THEN

                ! Obtain the "degree of freedom" index for this neighbor 
                e2 = modelMesh % ijToDOF( ni, nj )

                IF( e2 /= e1 .AND. e2 /= 0 )THEN
                   this % valence(e1) = this % valence(e1) + 1
                   this % neighbors(this % valence(e1), e1) = e2
                ENDIF

             ENDIF

          ENDDO

      ENDDO

      CALL CPU_TIME( t2 )
      PRINT *,' S/R BuildFromMeshAndStencil : Completed in ', t2-t1, ' seconds.'

 END SUBROUTINE BuildFromMeshAndL5OStencil_AdjacencyGraph

 SUBROUTINE Free_AdjacencyGraph( this )
 ! This subroutine frees memory held by the allocatable attributes of the
 ! AdjacencyGraph data structure.
 ! ======================================================================= !
   IMPLICIT NONE
   CLASS( AdjacencyGraph ), INTENT(inout) :: this

      DEALLOCATE( this % valence, &
                  this % color, &
                  this % neighbors )

 END SUBROUTINE Free_AdjacencyGraph
!
!
 SUBROUTINE GreedyColoring_AdjacencyGraph( this )
   ! Executes the Greedy coloring algorithm to fill in the "color" attribute
   ! given knowledge of the adjacency graph specified in the "neighbors"
   ! attribute.
   ! ====================================================================== !
   IMPLICIT NONE
   CLASS( AdjacencyGraph ), INTENT(inout) :: this
   ! Local
   INTEGER :: e1, e2
   INTEGER :: i, j
   INTEGER :: c2(1:this % maxValence)
   LOGICAL :: neednewcolor, colorfound
   REAL :: t1, t2

      PRINT*, ' S/R GreedyColoring : Start! '
     
      CALL CPU_TIME( t1 )

      this % ncolors = 1
      this % color(1) = 1

      DO e1 = 2, this % nDOF

         ! First build a list of all the colors of this element's neighbors
         c2 = 0
         DO i = 1, this % valence(e1) 
            e2    = this % neighbors(i,e1)
            c2(i) = this % color(e2)
         ENDDO

         neednewcolor = .TRUE.

         DO i = 1, this % nColors  ! only search through the colors already in use in the graph

           colorfound = .FALSE.
            DO j = 1, this % valence(e1) ! Search through the neighbor-color list           
               IF( c2(j) == i )THEN         ! and determine if this color is already in use by any neighbors
                  colorfound = .TRUE.
               ENDIF
            ENDDO
            ! If we do not find color "i" then we can assign element "e1" the
            ! color "i"
            IF( .NOT. colorfound )THEN
               this % color(e1) = i
               neednewcolor = .FALSE.
               EXIT
            ENDIF

         ENDDO

         ! In the event all of the possible colors are found in e1's neighbors,
         ! e1 must be colored with a new color      
         IF( neednewcolor )THEN
            this % nColors   = this % nColors + 1
            this % color(e1) = this % nColors
         ENDIF

      ENDDO
      
      CALL CPU_TIME( t2 )
      this % colored = .TRUE.
      
      PRINT *,' S/R GreedyColoring : Coloring completed with ', this % nColors, ' colors.'
      PRINT *,' S/R GreedyColoring : Coloring completed in ', t2-t1, ' seconds.'
      
 END SUBROUTINE GreedyColoring_AdjacencyGraph

 FUNCTION ImpulseFields( this ) RESULT( dofArray )
 !! Uses the colored graph to create arrays of impules field in DOF coordinates
 !!  
   IMPLICIT NONE
   CLASS( AdjacencyGraph ) :: this
   REAL(prec) :: dofArray(1:this % nDOF, 1:this % nColors)
   ! Local
   INTEGER :: i, j
   
     dofArray = 0.0_prec
     
     DO i = 1, this % nDOF
       j = this % color(i)
       dofArray(i,j) = 1.0_prec
     ENDDO
 
 END FUNCTION ImpulseFields
 
 FUNCTION DenseMatrix( this, irfDof ) RESULT( A )
 !! Uses the adjacency graph and the impulse response function to create a 
 !! dense matrix
   IMPLICIT NONE
   CLASS( AdjacencyGraph ) :: this
   REAL(prec) :: irfDof(1:this % nDOF, 1:this % nColors)
   REAL(prec) :: A(1:this % nDOF, 1:this % nDOF)
   ! Local
   INTEGER :: row, col, color, nid
   
     A = 0.0_prec
     DO col = 1, this % nDOF
       color = this % color(col)
       A(col,col) = irfDof(col,color)
       DO nid = 1, this % valence(col) 
         row = this % neighbors(nid,col)
         A(row,col) = irfDof(row,color)
       ENDDO
     ENDDO
     
 END FUNCTION DenseMatrix
       

! SUBROUTINE Read_HDF5( this, filename )
!   IMPLICIT NONE
!   CLASS( AdjacencyGraph ), INTENT(out) :: this
!   CHARACTER(*), INTENT(in)                 :: filename
!   ! Local
!   INTEGER(HID_T) :: file_id
!   INTEGER(HSIZE_T) :: v2dim(1:2), vdim(1:1)
!   INTEGER(HID_T)   :: dataset_id, filespace
!   INTEGER          :: error
!
!     CALL h5open_f(error)
!     CALL h5fopen_f(TRIM(filename), H5F_ACC_RDWR_F, file_id, error)
!
!     CALL Get_HDF5_Obj_Dimensions( file_id,'/graph/cell_valence', 1, vdim )
!     CALL Get_HDF5_Obj_Dimensions( file_id,'/graph/cell_neighbors', 2, v2dim )
!
!     CALL this % Build( nDOF = INT(vdim(1),4), &
!                           maxValence = INT(v2dim(1),4) )
!
!     ! Get the cell valence 
!     CALL h5dopen_f(file_id, '/graph/cell_valence', dataset_id, error)
!     CALL h5dget_space_f( dataset_id, filespace, error )
!     CALL h5dread_f( dataset_id, H5T_STD_I32LE, &
!                     this % valence, &
!                     vdim, error)
!     CALL h5dclose_f(dataset_id, error)
!     CALL h5sclose_f(filespace, error)
!
!     ! Get the cell color 
!     CALL h5dopen_f(file_id, '/graph/cell_color', dataset_id, error)
!     CALL h5dget_space_f( dataset_id, filespace, error )
!     CALL h5dread_f( dataset_id, H5T_STD_I32LE, &
!                     this % color, &
!                     vdim, error)
!     CALL h5dclose_f(dataset_id, error)
!     CALL h5sclose_f(filespace, error)
!
!     ! Get the cell neighbors 
!     CALL h5dopen_f(file_id, '/graph/cell_neighbors', dataset_id, error)
!     CALL h5dget_space_f( dataset_id, filespace, error )
!     CALL h5dread_f( dataset_id, H5T_STD_I32LE, &
!                     this % neighbors, &
!                     vdim, error)
!     CALL h5dclose_f(dataset_id, error)
!     CALL h5sclose_f(filespace, error)
!
!     CALL h5fclose_f(file_id, error)
!     CALL h5close_f(error)
!
!     this % nColors = MAXVAL(this % color)
!
! END SUBROUTINE Read_HDF5
! 
! SUBROUTINE Write_HDF5( this, filename )
!   IMPLICIT NONE
!   CLASS( AdjacencyGraph ), INTENT(in) :: this
!   CHARACTER(*), INTENT(in)                :: filename
!   ! Local
!   INTEGER(HID_T) :: file_id
!   INTEGER(HSIZE_T) :: vdim(1:1), v2dim(1:2)
!   INTEGER(HID_T)   :: group_id
!   INTEGER          :: error
!
!       CALL h5open_f(error)
!
!       CALL h5fcreate_f(TRIM(filename), H5F_ACC_TRUNC_F, file_id, error)
!
!       CALL h5gcreate_f( file_id, "/graph", group_id, error )
!       CALL h5gclose_f( group_id, error )
!
!       vdim(1:1) = (/this % nDOF/)
!       v2dim(1:2) = (/this % maxValence, this % nDOF/)
!
!       ! Write the cell valence
!       CALL Add_IntObj_to_HDF5( rank=1,&
!                                dimensions=vdim(1:1),&
!                                variable_name='/graph/cell_valence',&
!                                variable=this % valence,&
!                                file_id=file_id )
!       ! Write the cell color
!       CALL Add_IntObj_to_HDF5( rank=1,&
!                                dimensions=vdim(1:1),&
!                                variable_name='/graph/cell_color',&
!                                variable=this % color,&
!                                file_id=file_id )
!       ! Write the cell neighbors
!       CALL Add_IntObj_to_HDF5( rank=2,&
!                                dimensions=v2dim(1:2),&
!                                variable_name='/graph/cell_neighbors',&
!                                variable=this % neighbors,&
!                                file_id=file_id )
!
!       CALL h5fclose_f( file_id, error )
!       CALL h5close_f( error )
!
! END SUBROUTINE Write_HDF5

END MODULE SLSpectra_AdjacencyGraph
