
module sparse_header_debug

  use sparse_header

  implicit none

  type, extends(SparseCsrComplex) :: SparseCsrComplexDbg
    contains
    procedure :: print_values => sparse_csr_print_complex_values  
  end type SparseCsrComplexDbg

  type, extends(SparseCsrReal   ) :: SparseCsrRealDbg
    contains
    procedure :: print_values => sparse_csr_print_real_values  
  end type SparseCsrRealDbg

  contains
    
!===============================================================================
! SPARSE_CSR_PRINT_VALUES prints the column_index and value(complex) in pairs 
! and also the related row pointers (not actual pointers)
! This is a routine for debugging.
!===============================================================================

  subroutine sparse_csr_print_complex_values(A)
    class(SparseCsrComplexDbg) :: A
    integer :: i

    !print *, 'made it here'
    do i = 1, A % n_nonzero
      write(*,*) A % columns(i), A % values(i)
    end do

    write(*,*)

    do i = 1, A % m + 1    
      write(*,*) A % row_ptr(i)
    end do

  end subroutine

  subroutine sparse_csr_print_real_values(A)
    class(SparseCsrRealDbg) :: A
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
! SPARSE_CSR_PRINT_DENSE prints the dense matrix form of the sparse matrix   
! This is a routine for debugging.
!===============================================================================

  subroutine sparse_csr_print_complex_dense(A)
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

  subroutine sparse_csr_print_real_dense(A)
    class(SparseCsrReal) :: A
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

end module
