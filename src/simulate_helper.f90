!> @file simulate_helper.f90
!> @provides additional routines for simulation
!> @author ap
module simulate_helper_m
contains
!> @brief function for calculating optimized truncation point index for interaction potential
!> @param[in]       n           number of sites in spin chain
!> @param[in]       j           the spin-spin coupling constant
!> @param[in]       h           the spin-field coupling constant
!> @param[in]       s1          spin matrix along spin chain
!> @param[in]       s2          spin matrix parpendiculat to spin chain
!> @return          val         truncation index
  function vi_truncation_index(n, j, h, s1, s2)result(val)
    use global_m
    use qd_helper_m
    use ising_model_m
    use matrix_generator_m
    use linear_algebra_helper_m
    implicit none
    ! io variables
    integer, intent(in)                       :: n
    double precision, intent(in)              :: j
    complex(8), intent(in)                    :: h
    complex(8), dimension(:,:), intent(in)    :: s1
    complex(8), dimension(:,:), intent(in)    :: s2
    integer                                   :: val
    ! internal variables
    integer                                   :: dim
    integer                                   :: dimm
    complex(8), allocatable, dimension(:,:)   :: ham0
    complex(8), allocatable, dimension(:,:)   :: v_s
    complex(8), allocatable, dimension(:,:)   :: v_i1
    complex(8), allocatable, dimension(:,:)   :: v_i2
    complex(8), allocatable, dimension(:,:)   :: diff
    double precision, parameter               :: t=1.d0
    logical                                   :: stat
    integer                                   :: idx
    !
    dim = size(s1,1)
    dimm = dim**n
    !
    ham0 =  spin_spin_interaction(s1, n, "o", j)
    v_s = zeeman_term(s2, n, h)
    !
    stat = .false.
    allocate(diff(dimm,dimm))
    do idx = 1,20
      v_i1 = cmplx(0.d0,0.d0)
      v_i2 = cmplx(0.d0,0.d0)
      diff = cmplx(0.d0,0.d0)
      v_i1 = operator_interaction(v_s, ham0, t, idx)
      v_i2 = operator_interaction(v_s, ham0, t, idx+1)
      diff = v_i1 - v_i2
      call if_null_c(diff, stat)
      if (stat) then
        val = idx
        return
      end if
    end do
    deallocate(diff)
  end function vi_truncation_index

end module simulate_helper_m
