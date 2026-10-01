!> @file run_controller.f90
!> @brief controls final run
!> @author ap
module run_controller_m
contains
!> @brief subroutine for final simulation
  subroutine simulate()
    use run_driver_m
    implicit none
    !
    call time_evol_transverse_ising_1d()
  end subroutine simulate
end module run_controller_m
