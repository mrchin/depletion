module sparse_header

  use constants, only: NULL_COLUMN, ZERO

  implicit none

!===============================================================================
! SPARSECSR stores indices and values of a sparse matrix stored in
! compressed sparse row (CSR) format. This format has the non-zeros of the
! matrix stored in contiguous memory locations. The 'columns' array stores the
! column indices of the elements of the non-zero array and the 'row_ptr' array
! stores stores the indices in the columns and values arrays corresponding to
! beginning of each row. The total number of elements required is:
! 4 + 2*n_maximum + (m + 1).
!
! Types are extended for Real*8 and Complex*8 in the values array
!
! SPARSECSRREAL    defines  real(8)   values(:) array
! SPARSECSRCOMPLEX defines complex(8) values(:) array
!
! This pattern extends to init and insert subroutines
!===============================================================================

  type :: SparseCsr
    ! Information about size of matrix
    integer :: m         ! number of rows
    integer :: n         ! number of columns
    integer :: n_nonzero ! number of non-zero elements
    integer :: n_maximum ! actual maximum number preallocated 2x n_nonzero

    ! Storage of locations and values
    integer,    allocatable :: row_ptr(:) ! indices within values array
    integer,    allocatable :: columns(:) ! column indices with non-zeros
  end type

  type, extends(SparseCsr) :: SparseCsrComplex
    complex(8), allocatable :: values(:)  ! non-zero complex values in matrix
  contains
    procedure :: init => sparse_csr_init_complex
    procedure :: insert => sparse_csr_insert_complex
  end type

  type, extends(SparseCsr) :: SparseCsrReal
    real(8), allocatable :: values(:)     ! non-zero real values in matrix
  contains
    procedure :: init => sparse_csr_init_real
    procedure :: insert => sparse_csr_insert_real
    procedure :: fill_in => compute_fill_in_real
    procedure :: col_below_diag => sparse_column_below_diagonal_real
    procedure :: row_rgtof_diag => sparse_row_rgtof_diagonal_real
  end type

contains

!===============================================================================
! SPARSE_CSR_INIT allocates space for a compressed space row based on the size
! of the matrix and the expected number of non-zeros.
!===============================================================================

  subroutine sparse_csr_init(A, m, n, n_nonzero)

    class(SparseCsr) ::        A
    integer, intent(in)     :: m
    integer, intent(in)     :: n
    integer, intent(in)     :: n_nonzero

    ! Set integer constants
    A % m = m
    A % n = n
    A % n_nonzero = n_nonzero
    A % n_maximum = 2 * n_nonzero

    ! If any components are already allocated, remove that allocation
    if (allocated(A % row_ptr)) deallocate(A % row_ptr)
    if (allocated(A % columns)) deallocate(A % columns)

    ! Allocate components
    allocate(A % row_ptr(m + 1))
    allocate(A % columns(A % n_maximum))

    ! Initialize component arrays
    A % row_ptr = 0
    A % columns = NULL_COLUMN

  end subroutine sparse_csr_init

!===============================================================================
! SPARSE_CSR_INIT_COMPLEX allocates space for the values array using complex  
! values, and extends the allocation performed by SPARSE_CSR_INIT
!===============================================================================

  subroutine sparse_csr_init_complex(A, m, n, n_nonzero)

    class(SparseCsrComplex) :: A
    integer, intent(in)     :: m
    integer, intent(in)     :: n
    integer, intent(in)     :: n_nonzero

    call sparse_csr_init(A, m, n, n_nonzero)

    ! Allocate type specific components
    if (allocated(A % values)) deallocate(A % values)
    allocate(A % values(A % n_maximum))

    ! Type specific to complex
    A % values = ZERO

  end subroutine sparse_csr_init_complex

!===============================================================================
! SPARSE_CSR_INIT_REAL allocates space for the values array using real
! values, and extends the allocation performed by SPARSE_CSR_INIT
!===============================================================================

  subroutine sparse_csr_init_real(A, m, n, n_nonzero)

    class(SparseCsrReal) :: A
    integer, intent(in)     :: m
    integer, intent(in)     :: n
    integer, intent(in)     :: n_nonzero

    call sparse_csr_init(A, m, n, n_nonzero)

    ! Allocate type specific components
    if (allocated(A % values)) deallocate(A % values)
    allocate(A % values(A % n_maximum))

    ! Type specific to real(8)
    A % values = 0.  

  end subroutine sparse_csr_init_real   

!===============================================================================
! SPARSE_ROW_RGTOF_DIAGONAL_REAL
!===============================================================================

  subroutine sparse_row_rgtof_diagonal_real(A, k, values, nzcols, nnz)
    class(SparseCsrReal) :: A
    integer, intent(in)  :: k ! row/column diagonal in square matrix
    real(8), intent(out) :: values( A % m )
    integer, intent(out) :: nzcols( A % m ), nnz
    integer i, first_col_idx, last_col_idx, row_idx, nnz_count

    values = 0.D0
    nzcols = 0
    nnz_count = 0

    first_col_idx = A % row_ptr(k)
    last_col_idx = A % row_ptr(k + 1) - 1

    do i = first_col_idx, last_col_idx
      if (A % columns(i) > k) then
        nnz_count = nnz_count + 1
        values(nnz_count) = A % values(i)
        nzcols(nnz_count) = A % columns(i)
      end if
    end do

    nnz = nnz_count

  end subroutine sparse_row_rgtof_diagonal_real

!===============================================================================
! SPARSE_COLUMN_BELOW_DIAGONAL_REAL
!===============================================================================

  subroutine sparse_column_below_diagonal_real(A, k, values, nzrows, nnz)
    class(SparseCsrReal) :: A
    integer, intent(in)  :: k ! row/column diagonal in square matrix
    real(8), intent(out) :: values( A % m )
    integer, intent(out) :: nzrows( A % m ), nnz
    integer i, first_col_idx, last_col_idx, row_idx, nnz_count

    values = 0.D0
    nzrows = 0
    nnz_count = 0

    do row_idx = k+1, A % m

      first_col_idx = A % row_ptr(row_idx)
      last_col_idx = A % row_ptr(row_idx + 1) - 1

      do i = first_col_idx, last_col_idx
        if (k == A % columns(i)) then
          nnz_count = nnz_count + 1
          values(nnz_count) = A % values(i)
          nzrows(nnz_count) = row_idx
        end if
      end do

    end do

    nnz = nnz_count

  end subroutine sparse_column_below_diagonal_real

!===============================================================================
! SPARSE_MAX_BELOW_DIAGONAL
!===============================================================================

  subroutine sparse_max_below_diagonal_real(A, k, mmaxval, mmaxloc)
    class(SparseCsrReal) :: A
    integer, intent(in)  :: k ! column
    real(8), intent(out) :: mmaxval
    integer, intent(out) :: mmaxloc

    integer i, first_col_idx, last_col_idx, row_idx

    mmaxval = 0.D0
    mmaxloc = 0

    do row_idx = 1, A % m

      ! Evaluate the columns that are zero in the row of interest
      first_col_idx = A % row_ptr(row_idx)
      last_col_idx = A % row_ptr(row_idx + 1) - 1

      do i = first_col_idx, last_col_idx
      
        if (k == A % columns(i)) then
          if (abs(A % values(i)) > mmaxval) then
            mmaxval = abs(A % values(i))
            mmaxloc = row_idx
          end if
        end if

      end do

    end do

  end subroutine
!===============================================================================
! SPARSE_CSR_INSERT inserts a new value to csr sparse matrix. This is used when
! computing fill-in elements. A new value requires inserting to columns(:),
! values(:), and adjusting row_ptr(:), and increasing n_nonzero by one.
! The inputs are traditional row,col,val (regular matrix triplet)
!===============================================================================

  subroutine sparse_csr_insert(A, row, col, current_col_idx)

    class(SparseCsr) :: A     ! sparse matrix
    integer, intent(in)     :: row   ! row to insert using reg. matrix notation
    integer, intent(in)     :: col   ! col to insert using reg. matrix notation
    integer, intent(out)    :: current_col_idx
    integer i, first_col_idx, last_col_idx
    logical at_zero_position

    ! Evaluate the columns that are zero in the row of interest
    first_col_idx = A % row_ptr(row)
    last_col_idx = A % row_ptr(row + 1) - 1

    at_zero_position = .false.

    do i = first_col_idx, last_col_idx
      current_col_idx = i

      if (col == A % columns(current_col_idx)) then
        print *, 'invalid insert at non-zero location'
        stop

      ! If the col pos is <  'least' column
      else if (current_col_idx == first_col_idx .and. &
               col < A % columns(current_col_idx)) then

        at_zero_position = .true.; exit

      ! If the col pos is somewhere in between
      else if (current_col_idx > 1 .and. &
               col < A % columns(current_col_idx) .and. &
               col > A % columns(current_col_idx-1)) then

        at_zero_position = .true.; exit

      ! If the col pos is >  'most' column
      else if (current_col_idx == last_col_idx .and. &
               col > A % columns(current_col_idx)) then
        current_col_idx = current_col_idx + 1

        at_zero_position = .true.; exit

      end if
     
      if (at_zero_position) exit

    end do

    if ( at_zero_position ) then
      ! Increase nonzero size
      A % n_nonzero = A % n_nonzero + 1 

      ! Shift all row pointers that are +1 from the row insertion
      shift_row_ptr : do i = row + 1, A % m + 1
        A % row_ptr(i) = A % row_ptr(i) + 1    
      end do shift_row_ptr

      ! Shift all columns from the end to the current col idx pos + 1
      shift_col_pos: do i = A % n_nonzero, current_col_idx + 1, -1
        A % columns(i) = A % columns(i-1)
      end do shift_col_pos
      
      ! The insert goes here
      A % columns(current_col_idx) = col           

    else
      print *, 'invalid insert at non-zero location'
      stop
    end if

  end subroutine sparse_csr_insert

!===============================================================================
! SPARSE_CSR_INSERT_COMPLEX extends on SPARSE_CSR_INSERT
! The inputs are traditional row,col,val (regular matrix triplet)
!===============================================================================

  subroutine sparse_csr_insert_complex(A, row, col, val)

    class(SparseCsrComplex) :: A     ! sparse matrix
    integer, intent(in)     :: row   ! row to insert using reg. matrix notation
    integer, intent(in)     :: col   ! col to insert using reg. matrix notation
    complex, intent(in)     :: val   ! value of item         
    integer current_col_idx, i

    call sparse_csr_insert(A, row, col, current_col_idx)

    shift_val_pos: do i = A % n_nonzero, current_col_idx+1, -1
      A % values(i) = A % values(i-1)
    end do shift_val_pos
    A % values(current_col_idx) = val           

  end subroutine sparse_csr_insert_complex

!===============================================================================
! SPARSE_CSR_INSERT_REAL extends on SPARSE_CSR_INSERT
! The inputs are traditional row,col,val (regular matrix triplet)
!===============================================================================

  subroutine sparse_csr_insert_real(A, row, col, val)

    class(SparseCsrReal)    :: A     ! sparse matrix
    integer, intent(in)     :: row   ! row to insert using reg. matrix notation
    integer, intent(in)     :: col   ! col to insert using reg. matrix notation
    real(8), intent(in)     :: val   ! value of item         
    integer current_col_idx, i

    call sparse_csr_insert(A, row, col, current_col_idx)

    shift_val_pos: do i = A % n_nonzero, current_col_idx+1, -1
      A % values(i) = A % values(i-1)
    end do shift_val_pos
    A % values(current_col_idx) = val           

  end subroutine sparse_csr_insert_real   

!===============================================================================
! COMPUTE_FILL_IN analyzes a csr sparse matrix (real)
! Needs VALUE_IN_ARRAY function
!===============================================================================

  function value_in_array(val, array, array_size) result(bool)
    integer, intent(in) :: val, array_size
    integer, intent(in) :: array(array_size)
    logical :: bool
    integer :: i
    
    bool = .FALSE.
    do i = 1, array_size
      if (val == array(i)) bool = .TRUE.
    end do

  end function

  subroutine compute_fill_in_real(A)

    class(SparseCsrReal) :: A ! sparse matrix 
    integer :: i,j,k, threshold, first_col_idx, last_col_idx
    integer :: set_size, set(A % n) ! columns in A
    integer :: front_size, front(A % n) ! columns in A
    integer :: fill_size, fill(A % n)
    logical :: already_added, already_in_set

    ! Pre-initialization
    set = 0

    ! v -------> w is the directed graph
    ! Wish to find some intermediate u, 
    ! Such that  v -----> u -----> w
    ! And that the vertex value of v < both u and w

    ! Initialization let set = A(v), threshold = 1 and F(v) = null
    ! p. 182 Tarjan & Rose

    ! Define 'set' as the vertices that are in current row

    ! For threshold := 1 until n
    do threshold = 1, A % n  

      ! Clear set and fill
      set_size = 0; fill_size = 0
      set = 0     ; fill = 0

      ! Record all vertices w reachable from v

      ! Evaluate the columns that are zero in the row of interest
      first_col_idx = A % row_ptr(threshold)
      last_col_idx = A % row_ptr(threshold + 1) - 1

      ! print *, first_col_idx, last_col_idx
      do i = first_col_idx, last_col_idx
        set_size = set_size + 1
        set(set_size) = A % columns(i)
      end do

      ! j is the working threshold value, up to threshold
      ! i is the working element in set or fill arrays
      ! k is the evaluation of the front array

      ! The work loops are nearly identical, but the 2nd work loop moves through
      ! the fill array determined in the 1st work loop

      do j = 1, threshold 

        ! work through the set array
        WORK_LOOP_1: do i = 1, set_size
          if (j /= threshold) then
            if (j == set(i)) then

              ! Evaluate the columns that are zero in the row of interest
              first_col_idx = A % row_ptr(j)
              last_col_idx = A % row_ptr(j + 1) - 1

              front_size = 0
              front = 0
              do k = first_col_idx, last_col_idx
                front_size = front_size + 1
                front(front_size) = A % columns(k)
                if (front(front_size) > j .and. &
                    threshold /= front(front_size)) then

                  ! Add front value if not already added
                  already_added = &
                    value_in_array(front(front_size),fill,fill_size)

                  ! Ensure front value is not in set
                  already_in_set = &
                    value_in_array(front(front_size),set,set_size)

                  if (.not. already_added .and. .not. already_in_set) then
                    fill_size = fill_size + 1
                    fill(fill_size) = front(front_size) 
                  end if
                end if
              end do
            end if
          end if
        end do WORK_LOOP_1

        ! work through the fill array
        WORK_LOOP_2: do i = 1, fill_size
          if (j /= threshold) then
            if (j == fill(i)) then

              ! Evaluate the columns that are zero in the row of interest
              first_col_idx = A % row_ptr(j)
              last_col_idx = A % row_ptr(j + 1) - 1

              front_size = 0
              front = 0
              do k = first_col_idx, last_col_idx
                front_size = front_size + 1
                front(front_size) = A % columns(k)

                if (front(front_size) > j .and. &
                    threshold /= front(front_size)) then

                  ! Add front value if not already added
                  already_added = &
                    value_in_array(front(front_size),fill,fill_size)

                  ! Ensure front value is not in set
                  already_in_set = &
                    value_in_array(front(front_size),set,set_size)

                  if (.not. already_added .and. .not. already_in_set) then
                    fill_size = fill_size + 1
                    fill(fill_size) = front(front_size) 
                  end if

                end if
              end do

            end if
          end if
        end do WORK_LOOP_2

      end do

      if (fill_size > 0) print *, 'fill array:', fill, ' row:', threshold

    end do

  end subroutine compute_fill_in_real
 
end module sparse_header
