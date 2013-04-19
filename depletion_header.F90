module depletion_header

  use constants

  implicit none

!===============================================================================
! DEPLETIONSTEP
!===============================================================================

  type :: DepletionStep
    integer :: units        ! MWd/kgU or days
    real(8) :: value        ! user-provided value
    real(8) :: time  = NONE ! length of depletion step in seconds
    real(8) :: power = ONE  ! power fraction for this step
  end type DepletionStep

end module depletion_header
