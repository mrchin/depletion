module depletion

  use constants
  use depletion_header
  !use eigenvalue,       only: run_eigenvalue
  !use global
  !use output,           only: header
  use sparse_header,    only: SparseCsrComplex

contains

!===============================================================================
! RUN_DEPLETION
!===============================================================================

!  subroutine run_depletion()
!
!    integer :: i
!
!    !disable if (master) call header("DEPLETION SIMULATION", level=1)
!
!    do i = 1, n_depletion_steps
!      ! call run_eigenvalue()
!    end do
!
!  end subroutine run_depletion

!===============================================================================
! SOLVE_CRAM solves the matrix exponential using the Chebyshev rational
! approximation method. This method is described in M. Pusa and J. Leppanen,
! "Computing the Matrix Exponential in Burnup Calculations," Nucl. Sci. Eng.,
! 164, p. 140--150 (2010) as well as M. Pusa, "Rational Approximations to the
! Matrix Exponential in Burnup Calculations," Nucl. Sci. Eng., 169, p. 155--167
! (2011).
!===============================================================================

  !subroutine solve_cram(A, x0, x, t)
  subroutine solve_cram(A, x0, x)

    type(SparseCsrComplex), intent(in) :: A ! burnup matrix
    complex(8), intent(in)  :: x0(:)        ! starting concentration vector
    complex(8), intent(out) :: x(:)         ! end-of-step concentration vector
    !real(8),    intent(in)  :: t            ! time interval

    integer :: i     ! row index
    integer :: j     ! column index
    integer :: k     ! loop index for CRAM order
    integer :: i_val ! index in columns/values
    integer :: n     ! size of solution vector
    complex(8), allocatable :: b(:) ! right-hand side
    complex(8), allocatable :: y(:) ! solution to linear system
    type(SparseCsrComplex) :: fill  ! fill-in matrix
    type(SparseCsrComplex) :: lu    ! matrix used during LU decomposition

    ! Allocate arrays
    n = size(x0)
    allocate(b(n))
    allocate(y(n))

    ! Perform symbolic factorization on the burnup matrix for LU decomposition
    call symbolic_factorization(A, fill)

    ! Multiply elements in LU by time interval
    !fill % values = cmplx(real(fill % values) * t, aimag(fill % values))

    ! Compute first term of solution
    x = CRAM_ALPHA0 * x0

    do k = 1, CRAM_ORDER/2
      ! Copy fill matrix
      lu = fill

      ! Calculate right hand side
      b = CRAM_ALPHA(k) * x0

      ! Subtract theta on diagonal terms
      ROWS: do i = 1, n
        COLUMNS: do i_val = lu % row_ptr(i), lu % row_ptr(i+1) - 1
          ! Get index of column
          j = lu % columns(i_val)

          ! Subtract theta from diagonal term
          if (i == j) lu % values(i_val) = lu % values(i_val) - CRAM_THETA(k)
        end do COLUMNS
      end do ROWS

      ! Perform gaussian elimination
      call numerical_elimination(lu, b, y)

      ! Add to solution
      x = x + TWO * y
    end do

    ! Free space for temporary arrays
    deallocate(b)
    deallocate(y)

  end subroutine solve_cram

!===============================================================================
! SYMBOLIC_FACTORIZATION computes the non-zero structure of the fill-in matrix
! resulting from Gaussian elimination on a sparse matrix 'A'. The algorithm used
! here is the FILL2 algorithm from D. J. Rose and R. E. Tarjan, "Algorithmic
! aspects of vertex elimination on directed graphs," SIAM J. Appl. Math, 40,
! pp. 176--197 (1978).
!===============================================================================

  subroutine symbolic_factorization(matrix, fill)

    type(SparseCsrComplex), intent(in)    :: matrix
    type(SparseCsrComplex), intent(inout) :: fill

    integer :: i         ! row index
    integer :: i_val     ! index in columns/values
    integer :: i_val2    ! another index in columns/values
    integer :: j         ! column index
    integer :: k         ! another column index
    integer :: m         ! number of rows   
    integer :: n         ! number of columns
    integer :: n_nonzero ! number of non-zeros
    integer :: extra     ! number of extra elements on each row
    integer :: i_Omega   ! length of Omega
    integer :: tmpcol    ! saved column used during sorting
    complex(8) :: tmpval ! saved value used during sorting
    integer, allocatable :: Omega(:) ! list of rows to check
    integer, allocatable :: a(:)     ! temporary fill list (one row)
    integer, allocatable :: fill_array(:)     ! temporary fill list (one row)

    real(8), parameter :: EXTRA_SPACE_FRAC = 0.1

    ! ==========================================================================
    ! Initialize fill-in matrix

    n = matrix % n ! Copy size of matrix
    m = matrix % m ! Copy size of matrix

    ! Allocate Omega and fill_row
    allocate(a(n))
    allocate(Omega(n))
    allocate(fill_array(n))
    print *, 'allocating a and Omega with size:', n

    ! Allocate F with 10% more non-zero values
    extra = int(EXTRA_SPACE_FRAC*n)
    n_nonzero = matrix % n_nonzero + n*extra
    call fill % init(n, n, n_nonzero)
    print *, 'n_nonzero estimated to be:', n_nonzero, extra, EXTRA_SPACE_FRAC

    ! Copy existing values from A into F
    do i = 1, m
      ! Get original number of non-zeros in column
      n_nonzero = matrix % row_ptr(i+1) - matrix % row_ptr(i) 

      ! Set column pointers for F matrix (with extra space)
      fill % row_ptr(i) = matrix % row_ptr(i) + (i-1)*extra

      ! Copy non-zero rows and values from A into F
      i_val  = fill % row_ptr(i)
      i_val2 = matrix % row_ptr(i)
      fill % columns(i_val  : i_val+n_nonzero-1) = &
           matrix % columns(i_val2 : i_val2+n_nonzero-1)
      fill % values(i_val : i_val+n_nonzero-1) = &
           matrix % values(i_val2 : i_val2+n_nonzero-1)

      ! Initialize extra values that were allocated -- the extra non-zero rows
      ! are set to -1 to indicate that they haven't been used yet
      !fill % columns(i_val+n_nonzero : i_val+n_nonzero+extra-1) = -1
      !fill % values(i_val+n_nonzero : i_val+n_nonzero+extra-1) = ZERO
    end do

    ! Set final column ptr
    fill % row_ptr(m+1) = fill % n_nonzero + 1

    !call fill % print_values

    ! ==========================================================================
    ! SYMBOLIC DECOMPOSITION

    fill_array = 0
    ! Indicate that the Omega set is empty
    i_Omega = 0

    ! Set the vector 'a' to zero
    a = 0

    ROWS: do i = 2, n

      fill_array(i) = i

      ! ========================================================================
      ! Find non-zeros in row 

      COLUMNS_IN_ROW_I: do i_val = matrix % row_ptr(i), matrix % row_ptr(i+1) - 1
        ! Get index of column
        j = matrix % columns(i_val)

        ! Save positions of non-zero index
        a(j) = i

        i_Omega = i_Omega + 1
        Omega(i_Omega) = i

      end do COLUMNS_IN_ROW_I
    print *, 'i=', i, ' a=', a    
    print *, 'Omega=', (Omega(i_val),i_val=1,i_Omega)
    print *, 'hereka'

      ! ========================================================================
      ! Add fill-in where necessary

      ! Loop while the list of neighbor rows is empty
      NEIGHBOR_LIST: do while (i_Omega > 0)
        ! Get last item on Omega list
        j = Omega(i_Omega)
        i_Omega = i_Omega - 1

        ! Now loop over the columns in row j
        print *, 'herlawjerlas'
        !print *, fill % columns(fill % row_ptr(j)) , fill % columns(fill % row_ptr(j+1) -1)
        COLUMNS_IN_ROW_J: do i_val = fill % row_ptr(j), fill % row_ptr(j+1) - 1
        !COLUMNS_IN_ROW_J: do k = fill % columns(fill % row_ptr(j)) , fill % columns(fill % row_ptr(j+1) -1)
          ! Get index of column
          k = fill % columns(i_val)
        ! Now loop over the SPANNING columns in row j, k_min to k_max
        print *, 'graph j->k',j,k,', and a(k)', a(k)

          ! Check for end of non-zero columns in this row
          if (k == NULL_COLUMN) exit

          cycle

          ! If column k is in lower triangular part and the k-th column in row i
          ! is a zero entry, then (i,k) should be added to the fill-in
          if (k < j .and. a(k) == 0) then
            print *, 'j < k and a(k) = 0 (i,k)=', i, k
            ! Index in columns for new entry in row j
            i_val2 = fill % row_ptr(i) + n_nonzero

            ! Check if enough extra space exists -- if not, add extra space
            if (i_val2 >= fill % row_ptr(i+1)) then
                print *, 'adding extra space'
            !    call fill % expand(i, extra)
            end if

            ! Add (k,j) to fill-in matrix
            fill % columns(i_val2) = k
            n_nonzero = n_nonzero + 1

            ! Save k to temporary fill list
            a(k) = 1

            ! If j < k < i, there is a prospect for a longer path to node i
            ! through the nodes j and k, and therefore the node k should be
            ! added to Omega
            if (k < i) then
              print *, 'a node k should be added to Omega'
              i_Omega = i_Omega + 1
              Omega(i_Omega) = k
            end if
          end if
        end do COLUMNS_IN_ROW_J
      end do NEIGHBOR_LIST
    end do ROWS

    ! Free space from the arrays a and Omega
    deallocate(a)
    deallocate(Omega)

    ! ==========================================================================
    ! SORT LIST OF NON-ZERO COLUMNS FOR EACH ROW

    do i = 1, n
      do i_val = fill % row_ptr(i) + 1, fill % row_ptr(i+1) - 1
        ! Get index of column
        j = fill % columns(i_val)

        ! Check for end of non-zero columns in this row
        if (j == NULL_COLUMN) exit

        ! Save value to move
        k = i_val
        tmpcol = fill % columns(i_val)
        tmpval = fill % values(i_val)

        do
          ! Check if insertion value is greater than (k-1)th value
          if (tmpcol >= fill % columns(k-1)) exit

          ! Move values over until hitting one that's not larger
          fill % columns(k) = fill % columns(k-1)
          fill % values(k)  = fill % values(k-1)
          k = k - 1

          ! Exit if we've reached the beginning of this row
          if (k == fill % row_ptr(i)) exit
        end do

        ! Put the original value into its new position
        fill % columns(k) = tmpcol
        fill % values(k)  = tmpval
      end do
    end do

    !call fill % print_values

  end subroutine symbolic_factorization

!===============================================================================
! NUMERICAL_ELIMINATION zeros out entries in the matrix below the diagonal using
! Gaussian elimination. A recent paper -- M. Pusa and J. Leppanen, "An efficient
! implementation of the Chebyshev Rational Approximation Method for solving the
! burnup equations," Proc. PHYSOR 2012, Knoxville, TN, Apr. 15-20, 2012 --
! argues that partial pivoting is not required (based on properties of the
! burnup matrix). The actual algorithm for elimination on the sparse matrix is
! the CELIMINATE algorithm from R. E. Tarjan, "Graph Theory and Gaussian
! Elimination," Stanford STAN-CS-75-526, November 1975.
!===============================================================================

  subroutine numerical_elimination(fill, b, x)

    type(SparseCsrComplex), intent(inout) :: fill ! fill matrix
    complex(8),         intent(inout) :: b(:) ! right-hand side
    complex(8),         intent(out)   :: x(:) ! solution

    integer :: i          ! row index
    integer :: i_val      ! index in columns/values
    integer :: i_val2     ! another index in columns/values
    integer :: j          ! column index
    integer :: k          ! another column index
    integer :: n          ! number of columns
    complex(8) :: fill_ij ! (i,j) element in fill matrix
    complex(8) :: fill_jj ! (j,j) element in fill matrix
    complex(8) :: fill_jk ! (j,k) element in fill matrix
    complex(8) :: frac    ! F_ij / F_jj
    complex(8) :: sum     ! temporary sum for back substitution
    complex(8), allocatable :: v(:)
    complex(8), allocatable :: diag(:)

    ! Size of matrix
    n = fill % n

    ! Allocate arrays for storing a single row and the diagonal
    allocate(v(n))
    allocate(diag(n))

    ! Initialize v and diag
    v = ZERO
    diag = ZERO

    ROWS: do i = 2, n

      ! Copy row i to vector v and save diagonal
      do i_val = fill % row_ptr(i), fill % row_ptr(i+1) - 1
        j = fill % columns(i_val)

        ! Check for end of non-zero columns in this row
        if (j == NULL_COLUMN) exit

        ! Copy value into v vector
        v(j) = fill % values(i_val)

        ! Save diagonal
        if (i == j) diag(j) = v(j)
      end do

      COL_IN_ROW_I: do i_val = fill % row_ptr(i), fill % row_ptr(i+1) - 1
        j = fill % columns(i_val)

        ! Check for end of non-zero columns in this row
        if (j == NULL_COLUMN) exit

        if (j < i) then
          fill_ij = v(j)
          fill_jj = diag(j)

          ! TODO: Check norm of diagonal term

          ! Update right-hand side term in row i
          frac = fill_ij/fill_jj
          b(i) = b(i) - frac*b(j)

          ! Update matrix elements
          COL_IN_ROW_J: do i_val2 = fill % row_ptr(j), fill % row_ptr(j+1) - 1
            k = fill % columns(i_val2)

            ! Check for end of non-zero columns in this row
            if (k == NULL_COLUMN) exit

            ! We only need to update the terms
            if (j < k) then
              fill_jk = fill % columns(i_val2)
              v(k) = v(k) - frac*fill_jk
            end if             
          end do COL_IN_ROW_J
        end if

      end do COL_IN_ROW_I

      ! ========================================================================
      ! UPDATE ROW I

      do i_val = fill % row_ptr(i), fill % row_ptr(i+1) - 1
        j = fill % columns(i_val)

        ! Update diagonal
        if (i == j) diag(j) = v(j)

        ! Update element (i,j)
        fill % values(i_val) = v(j)
      end do

    end do ROWS

    !===========================================================================
    ! BACK_SUBSTITUTION

    do i = n, 1, -1
      ! Initialize sum
      sum = ZERO

      do i_val = fill % row_ptr(i), fill % row_ptr(i+1) - 1
        ! Get index of column
        j = fill % columns(i_val)

        ! Check for end of non-zero columns in this row
        if (j == NULL_COLUMN) exit

        ! Add A_ij * x_j to the sum -- only use upper triangular portion
        if (j > i) sum = sum + fill % values(i_val) * x(j)
      end do

      ! TODO: check norm of diagonal

      ! Calculate i-th element of solution
      x(i) = (b(i) - sum)/diag(i)
    end do

    ! Free up space from arrays
    deallocate(v)
    deallocate(diag)

  end subroutine numerical_elimination

end module depletion
