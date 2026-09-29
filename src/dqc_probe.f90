!> @file dqc_probe.f90
!> @brief module for providing quantities for probing dynamical quantum chaos
!> @author ap
module dqc_probe_m
contains
!> @brief function for calculating loschmidt echo
!> @param[in]     psi1      one state vector
!> @param[in]     psi2      another state vector
!> @return        val       loschmidt echo value
  function le(psi1, psi2)result(val)
    use qd_helper_m
    implicit none
    ! io variables
    complex(8), intent(in), dimension(:)    :: psi1
    complex(8), intent(in), dimension(:)    :: psi2
    double precision                        :: val
    ! internal variables
    complex(8)                                :: inner_pdt
    !
    inner_pdt = inner_product_dis(psi1, psi2)
    val = abs(conjg(inner_pdt)*inner_pdt)
  end function le
end module dqc_probe_m
