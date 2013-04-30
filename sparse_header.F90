module sparse_header

  use constants, only: NULL_COLUMN, ZERO

  implicit none

!===============================================================================
! SPARSECSR stores indices and values of a sparse matrix stored in
! compressed sparse row (CSR) format. This format has the non-zeros of the
! matrix stored in contiguous memory locations. The 'columns' array stores the
! column indices of the elements of the non-zero array and the 'row_ptr' array
! stores stores the indices in the columns and values arrays corresponding to
! beginning of each row. The total number of elements required is 2*n_nonzero +
! n + 1.
!
! Types are extended for Real*8 and Complex*8 in the values array
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

  subroutine sparse_csr_init_real(A, m, n, n_nonzero)

    class(SparseCsrReal) :: A
    integer, intent(in)     :: m
    integer, intent(in)     :: n
    integer, intent(in)     :: n_nonzero

    call sparse_csr_init(A, m, n, n_nonzero)

    ! Allocate type specific components
    if (allocated(A % values)) deallocate(A % values)
    allocate(A % values(A % n_maximum))

    ! type specific to real(8)
    A % values = 0.  

  end subroutine sparse_csr_init_real   

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
    !complex, intent(in)     :: val   ! value of item         
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

      ! if the col pos is <  'least' column
      else if (current_col_idx == first_col_idx .and. &
               col < A % columns(current_col_idx)) then

        at_zero_position = .true.; exit

      ! if the col pos is somewhere in between
      else if (current_col_idx > 1 .and. &
               col < A % columns(current_col_idx) .and. &
               col > A % columns(current_col_idx-1)) then

        at_zero_position = .true.; exit

      ! if the col pos is >  'most' column
      else if (current_col_idx == last_col_idx .and. &
               col > A % columns(current_col_idx)) then
        current_col_idx = current_col_idx + 1

        at_zero_position = .true.; exit

      end if
     
      if (at_zero_position) exit

    end do

    if ( at_zero_position ) then
      ! increase nonzero size
      A % n_nonzero = A % n_nonzero + 1 

      ! shift all row pointers that are +1 from the row insertion
      shift_row_ptr : do i = row + 1, A % m + 1
        A % row_ptr(i) = A % row_ptr(i) + 1    
      end do shift_row_ptr

      ! shift all columns from the end to the current col idx pos + 1
      shift_col_pos: do i = A % n_nonzero, current_col_idx + 1, -1
        A % columns(i) = A % columns(i-1)
      end do shift_col_pos
      
      ! the insert goes here
      A % columns(current_col_idx) = col           

    else
      print *, 'invalid insert at non-zero location'
      stop
    end if

  end subroutine sparse_csr_insert

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


end module sparse_header
