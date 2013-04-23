
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

end module
