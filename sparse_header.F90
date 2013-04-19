module sparse_header

  use constants, only: NULL_COLUMN, ZERO

  implicit none

!===============================================================================
! SRARSECSRCOMPLEX stores indices and values of a sparse matrix stored in
! compressed sparse row (CSR) format. This format has the non-zeros of the
! matrix stored in contiguous memory locations. The 'columns' array stores the
! column indices of the elements of the non-zero array and the 'row_ptr' array
! stores stores the indices in the columns and values arrays corresponding to
! beginning of each row. The total number of elements required is 2*n_nonzero +
! n + 1.
!===============================================================================

  type SparseCsrComplex
    ! Information about size of matrix
    integer :: m         ! number of rows
    integer :: n         ! number of columns
    integer :: n_nonzero ! number of non-zero elements
    integer :: n_maximum ! actual maximum number preallocated 2x n_nonzero

    ! Storage of locations and values
    integer,    allocatable :: row_ptr(:) ! indices within values array
    integer,    allocatable :: columns(:) ! column indices with non-zeros
    complex(8), allocatable :: values(:)  ! non-zero values in matrix
  contains
    procedure :: init => sparse_csr_init_complex
    procedure :: expand => sparse_csr_expand_complex
    procedure :: print_values => sparse_csr_print_values  
    procedure :: print_dense => sparse_csr_print_dense   
    procedure :: print_graph => sparse_csr_print_graph   
  end type SparseCsrComplex

contains

!===============================================================================
! SPARSE_CSR_INIT allocates space for a compressed space row based on the size
! of the matrix and the expected number of non-zeros.
!===============================================================================

  subroutine sparse_csr_init_complex(A, m, n, n_nonzero)

    class(SparseCsrComplex) :: A
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
    if (allocated(A % values)) deallocate(A % values)

    ! Allocate components
    allocate(A % row_ptr(m + 1))
    allocate(A % columns(A % n_maximum))
    allocate(A % values(A % n_maximum))

    ! Initialize component arrays
    A % row_ptr = 0
    A % columns = NULL_COLUMN
    A % values = ZERO

  end subroutine sparse_csr_init_complex

!===============================================================================
! SPARSE_CSR_EXPAND adds extra space to a row in a sparse matrix. This is used
! in sparse factorization when fill-in elements are added.
!===============================================================================

  subroutine sparse_csr_expand_complex(A, i, space)

    class(SparseCsrComplex) :: A     ! sparse matrix
    integer, intent(in)     :: i     ! row to expand
    integer, intent(in)     :: space ! number of items to add

    integer :: n      ! original size of A%columns
    integer :: i_val  ! index in columns/values
    integer :: i_val2 ! another index in columns/values
    integer,    allocatable :: columns(:) ! expanded copy of columns
    complex(8), allocatable :: values(:)  ! expanded copy of values

    ! Get original size of columns
    n = size(A % columns)

    ! Allocate temporary arrays
    allocate(columns(n + space))
    allocate(values(n + space))

    ! Get pointers to start of rows i and i+1
    i_val = A % row_ptr(i)
    i_val2 = A % row_ptr(i+1)

    ! Copy rows 1, i-1
    columns(1 : i_val-1) = A % columns(1 : i_val-1)
    values(1 : i_val-1)  = A % values(1 : i_val-1)

    ! Copy row i
    columns(i_val : i_val2-1) = A % columns(i_val : i_val2-1)
    values(i_val : i_val2-1)  = A % values(i_val : i_val2-1)

    ! Initialize extra space in row i
    columns(i_val2 : i_val2+space-1) = -1
    values(i_val2 : i_val2+space-1)   = ZERO

    ! Copy rows i+1, n
    columns(i_val2+space:) = A % columns(i_val2:)
    values(i_val2+space:)  = A % values(i_val2:)

    ! Move allocation back
    call move_alloc(FROM=columns, TO=A%columns)
    call move_alloc(FROM=values, TO=A%values)

    ! Adjust row pointers
    A % row_ptr(i+1:) = A % row_ptr(i+1:) + space

    ! Increase number of non-zeros
    A % n_nonzero = A % n_nonzero + space

  end subroutine sparse_csr_expand_complex

!===============================================================================
! SPARSE_CSR_PRINT_VALUES prints the column_index and value(complex) in pairs 
! and also the related row pointers (not actual pointers)
! This is a routine for debugging.
!===============================================================================

  subroutine sparse_csr_print_values(A)
    class(SparseCsrComplex) :: A
    integer :: i

    do i = 1, A % n_nonzero
      write(*,*) A % columns(i), A % values(i)
    end do

    write(*,*)

    do i = 1, A % m + 1    
      write(*,*) A % row_ptr(i)
    end do

  end subroutine

!===============================================================================
! SPARSE_CSR_PRINT_DENSEF prints the dense matrix form of the sparse matrix   
! This is a routine for debugging.
!===============================================================================

  subroutine sparse_csr_print_dense(A)
    class(SparseCsrComplex) :: A
    integer :: i
    integer :: col, row
    logical :: col_nonzero_found

    ! for all columns in dense matrix
    
    row = 1

    SEARCH_ROW: do while (row <= A % n)

    col_nonzero_found = .FALSE.

        SEARCH_COL: do col = 1, A % n

        ! perform search
        do i = A % row_ptr(row) , A % row_ptr(row+1) - 1
            if (col == A % columns (i)) then
              print *, row, col, A % values (i), ' i =', i
              col_nonzero_found = .TRUE.
              exit 
            end if
        end do

        if (.not. col_nonzero_found) print *, row, col, 0
        col_nonzero_found = .FALSE.

        end do SEARCH_COL

    row = row + 1
    print *

    end do SEARCH_ROW

  end subroutine

!===============================================================================
! SPARSE_CSR_PRINT_GRAPH prints the directed graph for all elements in matrix 
! This is a routine for debugging.
!===============================================================================

  subroutine sparse_csr_print_graph(A)
    class(SparseCsrComplex) :: A
    integer :: i
    integer :: col, row
    logical :: col_nonzero_found

    ! for all columns in dense matrix
    
    row = 1

    SEARCH_ROW: do while (row <= A % n)

    col_nonzero_found = .FALSE.

        SEARCH_COL: do col = 1, A % n

        ! perform search
        do i = A % row_ptr(row) , A % row_ptr(row+1) - 1
            if (col == A % columns (i)) then
              print *, row, col , A % values (i), ' i =', i
              col_nonzero_found = .TRUE.
              exit 
            end if
        end do

        if (.not. col_nonzero_found) print *, row, col, 0
        col_nonzero_found = .FALSE.

        end do SEARCH_COL

    row = row + 1
    print *

    end do SEARCH_ROW

  end subroutine

end module sparse_header
