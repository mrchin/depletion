
program sparse
use sparse_header
use depletion_header
use depletion
type(SparseCsrComplex) :: A
type(SparseCsrComplex) :: fill  ! fill-in matrix

print *, 'sparse matrix tester'

  call A % init ( 7, 7, 19 ) 

  A % row_ptr(1) = 1
  A % row_ptr(2) = 4
  A % row_ptr(3) = 7
  A % row_ptr(4) = 9
  A % row_ptr(5) = 12
  A % row_ptr(6) = 14
  A % row_ptr(7) = 16
  A % row_ptr(8) = 20

  A % columns(1) = 1
  A % columns(2) = 2
  A % columns(3) = 4
  A % columns(4) = 2
  A % columns(5) = 5
  A % columns(6) = 7
  A % columns(7) = 3
  A % columns(8) = 6
  A % columns(9)  = 1
  A % columns(10) = 3
  A % columns(11) = 4
  A % columns(12) = 5
  A % columns(13) = 6
  A % columns(14) = 3
  A % columns(15) = 6
  A % columns(16) = 3
  A % columns(17) = 4
  A % columns(18) = 5
  A % columns(19) = 7

  A % values(1 ) = CMPLX( 1.00000 , 1.00000)
  A % values(2 ) = CMPLX( 1.00000 , 1.00000)
  A % values(3 ) = CMPLX( 1.00000 , 1.00000)
  A % values(4 ) = CMPLX( 1.00000 , 1.00000)
  A % values(5 ) = CMPLX( 1.00000 , 1.00000)
  A % values(6 ) = CMPLX( 1.00000 , 1.00000)
  A % values(7 ) = CMPLX( 1.00000 , 1.00000)
  A % values(8 ) = CMPLX( 1.00000 , 1.00000)
  A % values(9 ) = CMPLX( 1.00000 , 1.00000)
  A % values(10) = CMPLX( 1.00000 , 1.00000)
  A % values(11) = CMPLX( 1.00000 , 1.00000)
  A % values(12) = CMPLX( 1.00000 , 1.00000)
  A % values(13) = CMPLX( 1.00000 , 1.00000)
  A % values(14) = CMPLX( 1.00000 , 1.00000)
  A % values(15) = CMPLX( 1.00000 , 1.00000)
  A % values(16) = CMPLX( 1.00000 , 1.00000)
  A % values(17) = CMPLX( 1.00000 , 1.00000)
  A % values(18) = CMPLX( 1.00000 , 1.00000)
  A % values(19) = CMPLX( 1.00000 , 1.00000)

  call A % print_graph  
  stop
  call A % print_dense  
  call A % print_values 

  ! Perform symbolic factorization on the burnup matrix for LU decomposition
  print *, 'symbolic_factorization'
  call symbolic_factorization(A, fill)
  
end program
