
program sparse
!use constants, only : ZERO
use sparse_header_debug
use depletion_header
use depletion
!type(SparseCsrComplexDbg) :: ! A
type(SparseCsrRealDbg) :: B, C, A, R
real(8), allocatable :: rvalues(:), cvalues(:)
integer, allocatable :: nzrows(:), nzcols(:)
integer :: nnzr, nnzc
!type(SparseCsrComplex) :: fill  ! fill-in matrix

print *, 'sparse matrix tester'

!===============================================================================
! R Matrix Test
!===============================================================================

  call R % init (  6, 6,  15 )

  R % row_ptr(1) = 1;
  R % columns(1) = 1; R %  values(1) = 5;

  R % row_ptr(2) = 2;
  R % columns(2) = 2; R %  values(2) = 2;
  R % columns(3) = 3; R %  values(3) = -1;
  R % columns(4) = 4; R %  values(4) = 2;

  R % row_ptr(3) = 5;
  R % columns(5) = 3; R %  values(5) = 3;

  R % row_ptr(4) = 6;
  R % columns(6) = 1; R %  values(6) = -2;
  R % columns(7) = 4; R %  values(7) = 1;
  R % columns(8) = 5; R %  values(8) = 1;

  R % row_ptr(5) = 9;
  R % columns(9)  = 1; R %  values(9)  = -1;
  R % columns(10) = 4; R %  values(10) = -1;
  R % columns(11) = 5; R %  values(11) = 2;
  R % columns(12) = 6; R %  values(12) = -3;

  R % row_ptr(6) = 13;
  R % columns(13) = 1; R %  values(13) = -1;
  R % columns(14) = 2; R %  values(14) = -1;
  R % columns(15) = 6; R %  values(15) = 6;

  R % row_ptr(7) = 16;

  call R % fill_in
  allocate(rvalues(R % n), nzrows(R % n))
  allocate(cvalues(R % n), nzcols(R % n))
  call R % col_below_diag(1, rvalues, nzrows, nnzr) 
  call R % row_rgtof_diag(1, cvalues, nzcols, nnzc) 
  print*, '--- 1 ---'          
  print*, (rvalues(i),i=1,nnzr)
  print*, (nzrows(i),i=1,nnzr)
  print*, (cvalues(i),i=1,nnzc)
  print*, (nzcols(i),i=1,nnzc)
  call R % col_below_diag(2, rvalues, nzrows, nnzr) 
  call R % row_rgtof_diag(2, cvalues, nzcols, nnzc) 
  print*, '--- 2 ---'          
  print*, (rvalues(i),i=1,nnzr)
  print*, (nzrows(i),i=1,nnzr)
  print*, (cvalues(i),i=1,nnzc)
  print*, (nzcols(i),i=1,nnzc)
  call R % col_below_diag(3, rvalues, nzrows, nnzr) 
  call R % row_rgtof_diag(3, cvalues, nzcols, nnzc) 
  print*, '--- 3 ---'          
  print*, (rvalues(i),i=1,nnzr)
  print*, (nzrows(i),i=1,nnzr)
  print*, (cvalues(i),i=1,nnzc)
  print*, (nzcols(i),i=1,nnzc)
  call R % col_below_diag(4, rvalues, nzrows, nnzr) 
  call R % row_rgtof_diag(4, cvalues, nzcols, nnzc) 
  print*, '--- 4 ---'          
  print*, (rvalues(i),i=1,nnzr)
  print*, (nzrows(i),i=1,nnzr)
  print*, (cvalues(i),i=1,nnzc)
  print*, (nzcols(i),i=1,nnzc)
  call R % col_below_diag(5, rvalues, nzrows, nnzr) 
  call R % row_rgtof_diag(5, cvalues, nzcols, nnzc) 
  print*, '--- 5 ---'          
  print*, (rvalues(i),i=1,nnzr)
  print*, (nzrows(i),i=1,nnzr)
  print*, (cvalues(i),i=1,nnzc)
  print*, (nzcols(i),i=1,nnzc)
  call R % col_below_diag(6, rvalues, nzrows, nnzr) 
  call R % row_rgtof_diag(6, cvalues, nzcols, nnzc) 
  print*, '--- 6 ---'          
  print*, (rvalues(i),i=1,nnzr)
  print*, (nzrows(i),i=1,nnzr)
  print*, (cvalues(i),i=1,nnzc)
  print*, (nzcols(i),i=1,nnzc)
  stop

  call C % init ( 10, 10, 36 )

  C % row_ptr(1) = 1;
  C % columns(1) = 1; C % columns(2) = 3; C % columns(3) = 8

  C % row_ptr(2) = 4;
  C % columns(4) = 2; C % columns(5) = 5; C % columns(6) = 10

  C % row_ptr(3) = 7;
  C % columns(7) = 1; C % columns(8) = 3; C % columns(9) = 7

  C % row_ptr(4) = 10;
  C % columns(10) = 4; C % columns(11) = 9; C % columns(12) = 10

  C % row_ptr(5) = 13;
  C % columns(13) = 2; C % columns(14) = 5; 
  C % columns(15) = 9; C % columns(16) = 10;

  C % row_ptr(6) = 17;
  C % columns(17) = 6; C % columns(18) = 7; C % columns(19) = 8

  C % row_ptr(7) = 20;
  C % columns(20) = 3; C % columns(21) = 6; C % columns(22) = 7
  
  C % row_ptr(8) = 23;
  C % columns(23) = 1; C % columns(24) = 6; C % columns(25) = 8
  C % columns(26) = 9; C % columns(27) = 10; 

  C % row_ptr(9) = 28;
  C % columns(28) = 4; C % columns(29) = 5; 
  C % columns(30) = 8; C % columns(31) = 9;

  C % row_ptr(10) = 32;
  C % columns(32) = 2; C % columns(33) = 4; C % columns(34) = 5
  C % columns(35) = 8; C % columns(36) = 10; 

  C % row_ptr(11) = 37;

  C % values(1:36) = 1.00000

  !call C % print_values
  !call C % fill_in
  !stop

  call A % init ( 7, 7, 19 ) 

  A % row_ptr(1) = 1; A % row_ptr(2) = 4
  A % row_ptr(3) = 7; A % row_ptr(4) = 9
  A % row_ptr(5) = 12; A % row_ptr(6) = 14
  A % row_ptr(7) = 16; A % row_ptr(8) = 20

  A % columns(1) = 1; A % columns(2) = 2; A % columns(3) = 4
  A % columns(4) = 2; A % columns(5) = 5; A % columns(6) = 7
  A % columns(7) = 3; A % columns(8) = 6; A % columns(9)  = 1
  A % columns(10) = 3; A % columns(11) = 4; A % columns(12) = 5
  A % columns(13) = 6; A % columns(14) = 3; A % columns(15) = 6
  A % columns(16) = 3; A % columns(17) = 4; A % columns(18) = 5
  A % columns(19) = 7

  !A % values(1:19) = CMPLX( 1.00000 , 1.00000)
  A % values(1:19) = 1.00000 
  call A % print_values
  call A % fill_in
  stop
  
  !print *, 'printing values of case 1 matrix'
  !call A % print_values
  !print *
  !print *, 'printing values of case 1 matrix augment 1'
  !call A % insert( 1, 3, CMPLX( 100.00000 , 0.00000) )
  !call A % print_values
  !print *
  !print *, 'printing values of case 1 matrix augment 2a'
  !call A % insert( 1, 5, CMPLX( 200.00000 , 0.00000) )
  !call A % print_values
  !print *
  !print *, 'printing values of case 1 matrix augment 2b'
  !call A % insert( 1, 6, CMPLX( 210.00000 , 0.00000) )
  !call A % print_values
  !print *
  !print *, 'printing values of case 1 matrix augment 2c'
  !call A % insert( 1, 7, CMPLX( 220.00000 , 0.00000) )
  !call A % print_values
  !print *
  !print *, 'printing values of case 1 matrix augment 3a'
  !call A % insert( 2, 1, CMPLX( 300.00000 , 0.00000) )
  !call A % print_values
  !print *
  !print *, 'printing values of case 1 matrix augment 3b'
  !call A % insert( 2, 3, CMPLX( 310.00000 , 0.00000) )
  !call A % print_values
  !print *
  !print *, 'printing values of case 1 matrix augment 3c'
  !call A % insert( 2, 6, CMPLX( 320.00000 , 0.00000) )
  !call A % print_values
  !print *, 'printing values of case 1 matrix augment 3c'
  !call A % insert( 2, 6, CMPLX( 320.00000 , 0.00000) )
  !call A % print_values

  !stop

  call B % init ( 6, 6, 19 ) 

  B % row_ptr(1) = 1; B % row_ptr(2) = 3; B % row_ptr(3) = 6
  B % row_ptr(4) = 9; B % row_ptr(5) = 13; B % row_ptr(6) = 17
  B % row_ptr(7) = 20
  
  B % columns(1) = 1; B % columns(2) = 5; B % columns(3) = 1
  B % columns(4) = 2; B % columns(5) = 6; B % columns(6) = 2
  B % columns(7) = 3; B % columns(8) = 4; B % columns(9)  = 1
  B % columns(10) = 3; B % columns(11) = 4; B % columns(12) = 5
  B % columns(13) = 2; B % columns(14) = 4; B % columns(15) = 5 
  B % columns(16) = 6; B % columns(17) = 2; B % columns(18) = 5
  B % columns(19) = 6

  B % values(1 ) = 10; B % values(2 ) = -2; B % values(3 ) = 3
  B % values(4 ) = 9; B % values(5 ) = 3; B % values(6 ) = 7
  B % values(7 ) = 8; B % values(8 ) = 7; B % values(9 ) = 3
  B % values(10) = 8; B % values(11) = 7; B % values(12) = 5
  B % values(13) = 8; B % values(14) = 9; B % values(15) = 9
  B % values(16) = 13; B % values(17) = 4; B % values(18) = 2
  B % values(19) = -1


  call B % print_values
  print *
  call B % fill_in
  stop
  print *, 'printing values of case 2 matrix'
  print *
  print *, 'printing values of case 1 matrix augment 1'
  call B % insert( 1, 2, 100.00000_8 )
  call B % print_values
  print *
  print *, 'printing values of case 1 matrix augment 2'
  call B % insert( 1, 6, 200.00000_8 )
  call B % print_values
  print *
  print *, 'printing values of case 1 matrix augment 3'
  call B % insert( 3, 1, 300.00000_8 )
  call B % print_values
  print *
  print *, 'printing values of case 1 matrix augment 4'
  call B % insert( 3, 2, 300.00000_8 )
  call B % print_values

!  print *, 'printing dense values of case 2 matrix'
!  call B % print_dense 
  stop
!  call A % print_dense  
!  call A % print_values

!  call B % print_graph
!  call B % print_dense
!  call B % print_values 

  ! Perform symbolic factorization on the burnup matrix for LU decomposition
  !print *, 'symbolic_factorization'
  !call symbolic_factorization(A, fill)
  
end program
