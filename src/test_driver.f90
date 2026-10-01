!> @file test_driver.f90
!> @brief provides unit test routines
!> @author ap
module test_driver_m
contains
!> @brief test for  if_null_c routine in linear_algebra_helper.f90
  subroutine test_if_null_c(test_stat)
    use global_m
    use linear_algebra_helper_m
    use matrix_generator_m
    implicit none
    ! io variables
    logical, intent(out)            :: test_stat
    ! internal variables
    integer, parameter              :: dim=5
    complex(8), dimension(dim,dim)  :: mat1
    complex(8), dimension(dim,dim)  :: mat2
    complex(8), dimension(dim,dim)  :: mat3
    logical                         :: stat1
    logical                         :: stat2
    logical                         :: stat3
    !
    mat1 = null_matrix_complex(dim)
    call if_null_c(mat1, stat1)
    !
    mat2 = identity_matrix_complex(dim)
    call if_null_c(mat2, stat2)
    !
    mat3 = random_complex_matrix(dim,dim)
    call if_null_c(mat3, stat3)
    !
    test_stat = .false.
    if (stat1 .and. .not.stat2 .and. .not.stat3) then
      test_stat = .true.
    end if
  end subroutine test_if_null_c
!> @brief subroutine for testing calculate_trace_c routine in linear_algebra_helper.f90
  subroutine test_calculate_trace_c(test_stat)
    use global_m
    use linear_algebra_helper_m
    use matrix_generator_m
    implicit none
    ! io variables
    logical, intent(out)             :: test_stat
    ! internal variables
    integer, parameter               :: dim=5
    complex(8), dimension(2,2)       :: matp1
    complex(8), dimension(2,2)       :: matp2
    complex(8), dimension(2,2)       :: matp3
    complex(8), dimension(dim,dim)   :: mat1
    complex(8), dimension(dim,dim)   :: mat2
    logical                          :: stat1
    logical                          :: stat2
    logical                          :: stat3
    complex(8)                       :: tracep1
    complex(8)                       :: tracep2
    complex(8)                       :: tracep3
    complex(8)                       :: trace1
    complex(8)                       :: trace2
    !
    matp1 = pauli_matrices(1)
    matp2 = pauli_matrices(2)
    matp3 = pauli_matrices(3)
    call calculate_trace_c(matp1, tracep1)
    call calculate_trace_c(matp2, tracep2)
    call calculate_trace_c(matp3, tracep3)
    stat1 = .false.
    if(tracep1.eq.cmplx(0.0d0,0.d0) .and. tracep2.eq.cmplx(0.d0,0.d0) .and. tracep3.eq.cmplx(0.d0,0.d0)) then
      stat1 = .true.
    end if
    !
    mat1 = null_matrix_complex(dim)
    call calculate_trace_c(mat1, trace1)
    stat2 = .false.
    if (trace1 .eq. cmplx(0.d0,0.d0)) then
      stat2 = .true.
    end if
    !
    mat2 = identity_matrix_complex(dim)
    call calculate_trace_c(mat2, trace2)
    stat3 = .false.
    if (trace2 .eq. cmplx(5.0d0,0.d0)) then
      stat3 = .true.
    end if
    !
    test_stat = .false.
    if (stat1 .and. stat2 .and. stat3) then
      test_stat = .true.
    end if
  end subroutine test_calculate_trace_c
!> @brief subroutine for testing kron_product routine in liner_algebra_helper.f90
!> @todo can add more unit tests
  subroutine test_kron_product(test_stat)
    use matrix_generator_m
    use linear_algebra_helper_m
    implicit none
    ! io variables
    logical, intent(out)                    :: test_stat
    ! internal variables
    integer, parameter                      :: dim = 10
    complex(8), dimension(dim,dim)          :: mat1
    complex(8), allocatable, dimension(:,:) :: mat1kron
    complex(8), dimension(dim**2,dim**2)    :: mat1ref
    complex(8), dimension(dim**2,dim**2)    :: mat1dif
    !
    mat1 = identity_matrix_complex(dim)
    call kron_product(mat1, mat1, mat1kron)
    mat1ref = identity_matrix_complex(dim**2)
    mat1dif = mat1kron-mat1ref
    test_stat = .false.
    call if_null_c(mat1dif,test_stat)
  end subroutine test_kron_product
!> @brief subroutine for testing diagonalize_matrix in linear_algebra_helper_m
  subroutine test_diagonalize_matrix(test_stat)
    use global_m
    use matrix_generator_m
    use linear_algebra_helper_m
    implicit none
    ! io variables
    logical, intent(out)                        :: test_stat
    ! internal variables
    complex(8), dimension(2,2)                  :: matp1
    complex(8), dimension(2,2)                  :: matp2
    complex(8), dimension(2,2)                  :: matp3
    double precision, dimension(2)              :: diagp1
    double precision, dimension(2)              :: diagp2
    double precision, dimension(2)              :: diagp3
    complex(8), dimension(2,2)                  :: eigp1
    complex(8), dimension(2,2)                  :: eigp2
    complex(8), dimension(2,2)                  :: eigp3
    double precision, dimension(2), parameter   :: val = (/-1.d0, 1.d0/)
    logical                                     :: stat1
    logical                                     :: stat2
    logical                                     :: stat3
    !
    matp1 = pauli_matrices(1)
    matp2 = pauli_matrices(2)
    matp3 = pauli_matrices(3)
    call diagonalize_matrix(2, matp1, eigp1, diagp1)
    call diagonalize_matrix(2, matp2, eigp2, diagp2)
    call diagonalize_matrix(2, matp3, eigp3, diagp3)
    !
    if (abs(diagp1(1)-val(1)).le.tol .and. abs(diagp1(2)-val(2)).le.tol) then
      stat1 = .true.
    else
      stat1 = .false.
    end if
    !
    if (abs(diagp2(1)-val(1)).le.tol .and. abs(diagp2(2)-val(2)).le.tol) then
      stat2 = .true.
    else
      stat2 = .false.
    end if
    !
    if (abs(diagp3(1)-val(1)).le.tol .and. abs(diagp3(2)-val(2)).le.tol) then
      stat3 = .true.
    else
      stat3 = .false.
    end if
    !
    test_stat = .false.
    if (stat1 .and. stat2 .and. stat3) then
      test_stat = .true.
    end if
  end subroutine test_diagonalize_matrix
!> @brief subroutine for testing calculate_commutator_c routine in linear_algebra_helper.f90
  subroutine test_calculate_commutator_c(test_stat)
    use global_m
    use linear_algebra_helper_m
    use matrix_generator_m
    use math_helper_m
    implicit none
    ! io variables
    logical, intent(out)          :: test_stat
    ! internal variable
    complex(8), dimension(2,2)    :: matp1
    complex(8), dimension(2,2)    :: matp2
    complex(8), dimension(2,2)    :: matp3
    complex(8), allocatable       :: comp12(:,:)
    complex(8), allocatable       :: comp23(:,:)
    complex(8), allocatable       :: comp31(:,:)
    complex(8), dimension(2,2)    :: difp12
    complex(8), dimension(2,2)    :: difp23
    complex(8), dimension(2,2)    :: difp31
    logical                       :: stat1
    logical                       :: stat2
    logical                       :: stat3
    !
    matp1 =  pauli_matrices(1)
    matp2 =  pauli_matrices(2)
    matp3 =  pauli_matrices(3)
    !
    call calculate_commutator_c(matp1, matp2, comp12)
    call calculate_commutator_c(matp2, matp3, comp23)
    call calculate_commutator_c(matp3, matp1, comp31)
    !
    difp12 = comp12 - 2.d0*iota*levi_civita(1,2,3)*matp3
    difp23 = comp23 - 2.d0*iota*levi_civita(2,3,1)*matp1
    difp31 = comp31 - 2.d0*iota*levi_civita(3,1,2)*matp2
    !
    call if_null_c(difp12, stat1)
    call if_null_c(difp23, stat2)
    call if_null_c(difp31, stat3)
    !
    test_stat = .false.
    if (stat1 .and. stat2 .and. stat3) then
      test_stat = .true.
    end if
  end subroutine test_calculate_commutator_c
!> @brief subroutine for testing  check_hermiticity routine in linear_algebra_helper.f90
  subroutine test_check_hermiticity(test_stat)
    use global_m
    use linear_algebra_helper_m
    use matrix_generator_m
    implicit none
    ! io variables
    logical, intent(out)            :: test_stat
    ! internal variables
    integer, parameter              :: dim=13
    complex(8), dimension(2,2)      :: mat1
    complex(8), dimension(2,2)      :: mat2
    complex(8), dimension(2,2)      :: mat3
    complex(8), dimension(dim,dim)  :: mat4
    complex(8), dimension(dim,dim)  :: mat5
    logical                         :: stat1
    logical                         :: stat2
    logical                         :: stat3
    logical                         :: stat4
    logical                         :: stat5
    !
    mat1 = pauli_matrices(1)
    mat2 = pauli_matrices(2)
    mat3 = pauli_matrices(3)
    mat4 = identity_matrix_complex(dim)
    mat5 = random_complex_matrix(dim,dim)
    !
    stat1 = .false.
    stat2 = .false.
    stat3 = .false.
    stat4 = .false.
    stat5 = .false.
    !
    call check_hermiticity(mat1, stat1)
    call check_hermiticity(mat2, stat2)
    call check_hermiticity(mat3, stat3)
    call check_hermiticity(mat4, stat4)
    call check_hermiticity(mat5, stat5)
    !
    test_stat = .false.
    if (stat1 .and. stat2 .and. stat3 .and. stat4 .and. .not.stat5) then
      test_stat = .true.
    end if
  end subroutine test_check_hermiticity
!> @brief subroutine for testing unitarity_check in linear_algebra_helper.f90
  subroutine test_unitarity_check(test_stat)
    use global_m
    use linear_algebra_helper_m
    use matrix_generator_m
    implicit none
    ! io variables
    logical, intent(out)            :: test_stat
    ! internal variables
    integer, parameter              :: dim=13
    complex(8), dimension(2,2)      :: mat1
    complex(8), dimension(2,2)      :: mat2
    complex(8), dimension(2,2)      :: mat3
    complex(8), dimension(dim,dim)  :: mat4
    complex(8), dimension(dim,dim)  :: mat5
    logical                         :: stat1
    logical                         :: stat2
    logical                         :: stat3
    logical                         :: stat4
    logical                         :: stat5
    !
    mat1 = pauli_matrices(1)
    mat2 = pauli_matrices(2)
    mat3 = pauli_matrices(3)
    mat4 = identity_matrix_complex(dim)
    mat5 = random_complex_matrix(dim,dim)
    !
    stat1 = .false.
    stat2 = .false.
    stat3 = .false.
    stat4 = .false.
    stat5 = .false.
    !
    call unitarity_check(mat1, stat1)
    call unitarity_check(mat2, stat2)
    call unitarity_check(mat3, stat3)
    call unitarity_check(mat4, stat4)
    call unitarity_check(mat5, stat5)
    !
    test_stat = .false.
    if (stat1 .and. stat2 .and. stat3 .and. stat4 .and. .not.stat5) then
      test_stat = .true.
    end if
  end subroutine test_unitarity_check
!> @brief subroutine for testing inner_product_dis routine in qd_helper.f90
  subroutine test_inner_product_dis(test_stat)
    use global_m
    use qd_helper_m
    use linear_algebra_helper_m
    use matrix_generator_m
    implicit none
    ! io variables
    logical, intent(out)          :: test_stat
    ! internal variables
    complex(8), dimension(2,2)        :: mat1
    complex(8), dimension(2,2)        :: mat2
    complex(8), dimension(2,2)        :: eig1
    complex(8), dimension(2,2)        :: eig2
    double precision, dimension(2)    :: eigval1
    double precision, dimension(2)    :: eigval2
    complex(8)                        :: val11
    complex(8)                        :: val12
    complex(8)                        :: val13
    complex(8)                        :: val21
    complex(8)                        :: val22
    complex(8)                        :: val23
    logical                           :: stat1
    logical                           :: stat2
    !
    mat1 = pauli_matrices(1)
    mat2 = pauli_matrices(2)
    !
    call diagonalize_matrix(2, mat1, eig1, eigval1)
    call diagonalize_matrix(2, mat2, eig2, eigval2)
    !
    val11 = inner_product_dis(eig1(:,1),eig1(:,1))
    val12 = inner_product_dis(eig1(:,2),eig1(:,2))
    val13 = inner_product_dis(eig1(:,2),eig1(:,1))
    val21 = inner_product_dis(eig2(:,1),eig2(:,1))
    val22 = inner_product_dis(eig2(:,2),eig2(:,2))
    val23 = inner_product_dis(eig2(:,2),eig2(:,1))
    !
    stat1 = .false.
    if (abs(abs(val11)-1.d0).le.tol .and. abs(abs(val12)-1.d0).le.tol .and. abs(val13).le.tol) then
      stat1 = .true.
    end if
    stat2 = .false.
    if (abs(abs(val21)-1.d0).le.tol .and. abs(abs(val22)-1.d0).le.tol .and. abs(val23).le.tol) then
      stat2 = .true.
    end if
    !
    test_stat = .false.
    if (stat1 .and. stat2) then
      test_stat = .true.
    end if
  end subroutine test_inner_product_dis
!
!> @brief subroutine for testing inner_product_dis routine in qd_helper.f90
  subroutine test_expectation_value_dis(test_stat)
    use global_m
    use qd_helper_m
    use linear_algebra_helper_m
    use matrix_generator_m
    implicit none
    ! io variables
    logical, intent(out)              :: test_stat
    ! internal variables
    complex(8), dimension(2,2)        :: mat1
    complex(8), dimension(2,2)        :: mat2
    complex(8), dimension(2,2)        :: eig1
    complex(8), dimension(2,2)        :: eig2
    double precision, dimension(2)    :: eigval1
    double precision, dimension(2)    :: eigval2
    complex(8)                        :: val11
    complex(8)                        :: val12
    complex(8)                        :: val21
    complex(8)                        :: val22
    logical                           :: stat1
    logical                           :: stat2
    !
    mat1 = pauli_matrices(1)
    mat2 = pauli_matrices(2)
    !
    call diagonalize_matrix(2, mat1, eig1, eigval1)
    call diagonalize_matrix(2, mat2, eig2, eigval2)
    !
    val11 = expectation_value_dis(eig1(:,1), mat1)
    val12 = expectation_value_dis(eig1(:,2), mat1)
    val21 = expectation_value_dis(eig2(:,1), mat2)
    val22 = expectation_value_dis(eig2(:,2), mat2)
    !
    stat1 = .false.
    if (abs(real(val11)-eigval1(1)).le.tol .and. abs(real(val12)-eigval1(2)).le.tol) then
      stat1 = .true.
    end if
    stat2 = .false.
    if (abs(real(val21)-eigval2(1)).le.tol .and. abs(real(val22)-eigval2(2)).le.tol) then
      stat2 = .true.
    end if
    !
    test_stat = .false.
    if (stat1 .and. stat2) then
      test_stat = .true.
    end if
  end subroutine test_expectation_value_dis
!> @brief subroutine for testing factorial function in math_helper.f90
  subroutine test_factorial(test_stat)
    use global_m
    use math_helper_m
    implicit none
    ! io variables
    logical, intent(out)          :: test_stat
    ! internal variables
    integer, parameter            :: n = 5
    !
    test_stat = .false.
    if (abs(factorial(n)-120.d0) .le. tol) then
      test_stat = .true.
    end if
  end subroutine test_factorial
!> @brief subroutine for testing spin-spin interaction function in ising_model.f90
  subroutine test_spin_spin_interaction(test_stat)
    use global_m
    use ising_model_m
    use matrix_generator_m
    use linear_algebra_helper_m
    implicit none
    ! io variables
    logical, intent(out)                    :: test_stat
    ! internal variables
    complex(8), dimension(2,2)              :: sz
    complex(8), dimension(2,2)              :: sy
    complex(8), dimension(2,2)              :: sx
    complex(8), dimension(2**2,2**2)        :: spinz_int_2spin
    complex(8), dimension(2**2,2**2)        :: spinx_int_2spin
    complex(8), dimension(2**2,2**2)        :: spiny_int_2spin
    complex(8), dimension(2**2,2**2)        :: spinz_int_2spin_exp
    complex(8), dimension(2**2,2**2)        :: spinx_int_2spin_exp
    complex(8), dimension(2**2,2**2)        :: spiny_int_2spin_exp
    complex(8), dimension(2**3,2**3)        :: spinz_int_3spinobc
    complex(8), dimension(2**3,2**3)        :: spinz_int_3spinpbc
    complex(8), dimension(2**3,2**3)        :: spinz_int_3spinobc_exp
    complex(8), dimension(2**3,2**3)        :: spinz_int_3spinpbc_exp
    complex(8), dimension(2**2,2**2)        :: diff1
    complex(8), dimension(2**2,2**2)        :: diff4
    complex(8), dimension(2**2,2**2)        :: diff5
    complex(8), dimension(2**3,2**3)        :: diff2
    complex(8), dimension(2**3,2**3)        :: diff3
    logical                                 :: stat1
    logical                                 :: stat2
    logical                                 :: stat3
    logical                                 :: stat4
    logical                                 :: stat5
    !
    sz = pauli_matrices(3)
    !
    spinz_int_2spin = spin_spin_interaction(sz, 2, "o", 1.d0)
    spinz_int_2spin_exp = reshape((/cmplx(1.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), &
                                  cmplx(0.d0,0.d0), cmplx(-1.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), &
                                  cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(-1.d0,0.d0), cmplx(0.d0,0.d0), &
                                  cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(1.d0,0.d0)/), &
                        shape(spinz_int_2spin_exp), order=(/2,1/))
    !
    diff1 = spinz_int_2spin_exp-spinz_int_2spin
    call if_null_c(diff1, stat1)
    !
    spinz_int_3spinobc = spin_spin_interaction(sz, 3, "o", 1.d0)
    spinz_int_3spinobc_exp = reshape((/ cmplx(2.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), &
                                       cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), &
                                       !
                                       cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), &
                                       cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), &
                                       !
                                       cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(-2.d0,0.d0), cmplx(0.d0,0.d0), &
                                       cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), &
                                       !
                                       cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), &
                                       cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), &
                                       !
                                       cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), &
                                       cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), &
                                       !
                                       cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), &
                                       cmplx(0.d0,0.d0), cmplx(-2.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), &
                                       !
                                       cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), &
                                       cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), &
                                       !
                                       cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), &
                                       cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(2.d0,0.d0) &
                                     /), shape(spinz_int_3spinobc_exp), order=(/2,1/))
       !
       spinz_int_3spinpbc = spin_spin_interaction(sz, 3, "p", 1.d0)
       spinz_int_3spinpbc_exp = reshape((/ cmplx(3.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), &
                                       cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), &
                                       !
                                       cmplx(0.d0,0.d0), cmplx(-1.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), &
                                       cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), &
                                       !
                                       cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(-1.d0,0.d0), cmplx(0.d0,0.d0), &
                                       cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), &
                                       !
                                       cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(-1.d0,0.d0), &
                                       cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), &
                                       !
                                       cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), &
                                       cmplx(-1.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), &
                                       !
                                       cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), &
                                       cmplx(0.d0,0.d0), cmplx(-1.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), &
                                       !
                                       cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), &
                                       cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(-1.d0,0.d0), cmplx(0.d0,0.d0), &
                                       !
                                       cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), &
                                       cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(3.d0,0.d0) &
                                     /), shape(spinz_int_3spinpbc_exp), order=(/2,1/))
    !
    diff2 = spinz_int_3spinobc_exp - spinz_int_3spinobc
    diff3 = spinz_int_3spinpbc_exp - spinz_int_3spinpbc
    call if_null_c(diff2, stat2)
    call if_null_c(diff3, stat3)
    !
    sx = pauli_matrices(1)
    sy = pauli_matrices(2)
    spinx_int_2spin = spin_spin_interaction(sx, 2, "o", 1.d0)
    spiny_int_2spin = spin_spin_interaction(sy, 2, "o", 1.d0)
    !
    spinx_int_2spin_exp = reshape((/ cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(1.d0,0.d0), &
                                     cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(1.d0,0.d0), cmplx(0.d0,0.d0), &
                                     cmplx(0.d0,0.d0), cmplx(1.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), &
                                     cmplx(1.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0) &
                                  /), shape(spinx_int_2spin_exp), order=(/2,1/))
    spiny_int_2spin_exp = reshape((/ cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(-1.d0,0.d0), &
                                     cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(1.d0,0.d0), cmplx(0.d0,0.d0), &
                                     cmplx(0.d0,0.d0), cmplx(1.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), &
                                     cmplx(-1.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0) &
                                  /), shape(spiny_int_2spin_exp), order=(/2,1/))

    !
    diff4 = spinx_int_2spin_exp-spinx_int_2spin
    call if_null_c(diff4, stat4)
    diff5 = spiny_int_2spin_exp-spiny_int_2spin
    call if_null_c(diff5, stat5)
    !
    test_stat = .false.
    if (stat1 .and. stat2 .and. stat3 .and. stat4 .and. stat5) then
      test_stat = .true.
    end if
  end subroutine test_spin_spin_interaction
!> @brief subroutine for testing zeeman interaction function in ising_model.f90
  subroutine test_zeeman_term(test_stat)
    use global_m
    use ising_model_m
    use linear_algebra_helper_m
    use matrix_generator_m
    implicit none
    ! io variables
    logical, intent(out)                    :: test_stat
    ! internal variables
    complex(8), dimension(2,2)              :: sx
    complex(8), dimension(2,2)              :: sy
    complex(8), dimension(2,2)              :: sz
    complex(8), dimension(2**2,2**2)        :: spinz_zeeman_2spin
    complex(8), dimension(2**2,2**2)        :: spinx_zeeman_2spin
    complex(8), dimension(2**2,2**2)        :: spiny_zeeman_2spin
    complex(8), dimension(2**2,2**2)        :: spinz_zeeman_2spin_exp
    complex(8), dimension(2**2,2**2)        :: spiny_zeeman_2spin_exp
    complex(8), dimension(2**2,2**2)        :: spinx_zeeman_2spin_exp
    complex(8), dimension(2**2,2**2)        :: diff1
    complex(8), dimension(2**2,2**2)        :: diff2
    complex(8), dimension(2**2,2**2)        :: diff3
    logical                                 :: stat1
    logical                                 :: stat2
    logical                                 :: stat3
    complex(8), parameter                   :: h = cmplx(1.d0, 0.d0)
    !
    sx = pauli_matrices(1)
    sy = pauli_matrices(2)
    sz = pauli_matrices(3)
    !
    spinz_zeeman_2spin = zeeman_term(sz, 2, h)
    spiny_zeeman_2spin = zeeman_term(sy, 2, h)
    spinx_zeeman_2spin = zeeman_term(sx, 2, h)
    !
    spinz_zeeman_2spin_exp = reshape((/ cmplx(2.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), &
                                        cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), &
                                        cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0),&
                                        cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(-2.d0,0.d0) &
                                      /), shape(spinz_zeeman_2spin_exp), order=(/2,1/))
    spinx_zeeman_2spin_exp = reshape((/ cmplx(0.d0,0.d0), cmplx(1.d0,0.d0), cmplx(1.d0,0.d0), cmplx(0.d0,0.d0), &
                                        cmplx(1.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(1.d0,0.d0), &
                                        cmplx(1.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(1.d0,0.d0),&
                                        cmplx(0.d0,0.d0), cmplx(1.d0,0.d0), cmplx(1.d0,0.d0), cmplx(0.d0,0.d0) &
                                      /), shape(spinz_zeeman_2spin_exp), order=(/2,1/))
    spiny_zeeman_2spin_exp = reshape((/ cmplx(0.d0,0.d0), cmplx(0.d0,-1.d0), cmplx(0.d0,-1.d0), cmplx(0.d0,0.d0), &
                                        cmplx(0.d0,1.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,-1.d0), &
                                        cmplx(0.d0,1.d0), cmplx(0.d0,0.d0), cmplx(0.d0,0.d0), cmplx(0.d0,-1.d0),&
                                        cmplx(0.d0,0.d0), cmplx(0.d0,1.d0), cmplx(0.d0,1.d0), cmplx(0.d0,0.d0) &
                                      /), shape(spinz_zeeman_2spin_exp), order=(/2,1/))
    !
    diff1 = spinz_zeeman_2spin_exp - spinz_zeeman_2spin
    diff2 = spinx_zeeman_2spin_exp - spinx_zeeman_2spin
    diff3 = spiny_zeeman_2spin_exp - spiny_zeeman_2spin
    !
    call if_null_c(diff1, stat1)
    call if_null_c(diff2, stat2)
    call if_null_c(diff3, stat3)
    !
    test_stat = .false.
    if (stat1 .and. stat2 .and. stat3) then
      test_stat = .true.
    end if
  end subroutine test_zeeman_term
!> @brief subroutine for tesing magnetization operator function in ising_model.f90
  subroutine test_magnetization_operator(test_stat)
    use global_m
    use linear_algebra_helper_m
    use matrix_generator_m
    use ising_model_m
    use qd_helper_m
    implicit none
    ! io variables
    logical, intent(out)                     :: test_stat
    ! internal variables
    integer, parameter                       :: n=2
    complex(8), dimension(2,2)               :: sz
    complex(8), allocatable, dimension(:,:)  :: mag_op
    complex(8), dimension(2**n)              :: state
    complex(8)                               :: mag_val
    !
    sz = pauli_matrices(3)
    mag_op = magnetization_operator(sz, n)
    state = cmplx(0.d0, 0.d0)
    state(1) = cmplx(1.d0, 0.d0)
    mag_val = expectation_value_dis(state, mag_op)
    !
    test_stat = .false.
    if (abs(abs(mag_val)-n*1.d0) .le. tol) then
      test_stat = .true.
    end if
    !write(*,*) mag_val
    !write(*,*)  mag_op(1,:)
    !write(*,*)  mag_op(2,:)
    !write(*,*)  mag_op(3,:)
    !write(*,*)  mag_op(4,:)
  end subroutine test_magnetization_operator
!> @brief subroutine for testing transverse_ising1d_hamiltonian function in ising_model.f90
  subroutine test_transverse_ising1d_hamiltonian(test_stat)
    use global_m
    use linear_algebra_helper_m
    use matrix_generator_m
    use ising_model_m
    implicit none
    ! io variables
    logical, intent(out)                  :: test_stat
    ! internal variables
    complex(8), dimension(2,2)            :: sz
    complex(8), dimension(2,2)            :: sx
    integer, parameter                    :: n=2
    double precision                      :: j
    complex(8)                            :: h
    complex(8), dimension(2**n,2**n)      :: hamiltonian
    complex(8), dimension(2**n,2**n)      :: spin_intobc
    complex(8), dimension(2**n,2**n)      :: spin_intpbc
    complex(8), dimension(2**n,2**n)      :: field_int
    complex(8), dimension(2**n,2**n)      :: dif1
    complex(8), dimension(2**n,2**n)      :: dif2
    complex(8), dimension(2**n,2**n)      :: dif3
    logical                               :: stat1
    logical                               :: stat2
    logical                               :: stat3
    !
    sz = (hbar/2.d0)*pauli_matrices(3)
    sx = (hbar/2.d0)*pauli_matrices(1)
    !
    j = exp(1.d0)
    h = cmplx(0.d0,0.d0)
    hamiltonian = transverse_ising1d_hamiltonian(n, sz, sx, j, h, "o")
    spin_intobc = spin_spin_interaction(sz, n, "o", j)
    dif1 = hamiltonian - spin_intobc
    call if_null_c(dif1, stat1)
    !
    j = pi
    h = cmplx(0.d0,0.d0)
    hamiltonian = transverse_ising1d_hamiltonian(n, sz, sx, j, h, "p")
    spin_intpbc = spin_spin_interaction(sz, n, "p", j)
    dif2 = hamiltonian - spin_intpbc
    call if_null_c(dif2, stat2)
    !
    j = 0.d0
    h = cmplx(pi,exp(1.d0))
    hamiltonian = transverse_ising1d_hamiltonian(n, sz, sx, j, h, "p")
    field_int = zeeman_term(sx, n, h)
    dif3 = hamiltonian - field_int
    call if_null_c(dif3, stat3)
    !
    test_stat = .false.
    if (stat1 .and. stat2 .and. stat3) then
      test_stat = .true.
    end if
  end subroutine test_transverse_ising1d_hamiltonian
!> @brief subroutine for testing bch_c subroutine in math_helper.f90
  subroutine test_bch_c(test_stat)
    use global_m
    use math_helper_m
    use matrix_generator_m
    use linear_algebra_helper_m
    implicit none
    ! io variables
    logical, intent(out)              :: test_stat
    ! internal variables
    complex(8), dimension(2,2)                     :: sx
    complex(8), dimension(2,2)                     :: sy
    complex(8), dimension(2,2)                     :: sz
    complex(8), allocatable, dimension(:,:)        :: mat0
    complex(8), allocatable, dimension(:,:)        :: mat1
    complex(8), allocatable, dimension(:,:)        :: mat2
    complex(8), allocatable, dimension(:,:)        :: mat3
    complex(8), dimension(2,2)                     :: exp0
    complex(8), dimension(2,2)                     :: exp1
    complex(8), dimension(2,2)                     :: exp2
    complex(8), dimension(2,2)                     :: exp3
    complex(8), dimension(2,2)                     :: dif0
    complex(8), dimension(2,2)                     :: dif1
    complex(8), dimension(2,2)                     :: dif2
    complex(8), dimension(2,2)                     :: dif3
    logical                                        :: stat0
    logical                                        :: stat1
    logical                                        :: stat2
    logical                                        :: stat3
    !
    sx = pauli_matrices(1)
    sy = pauli_matrices(2)
    sz = pauli_matrices(3)
    !
    call bch_c(sx, sz, 0, mat0)
    exp0 = sz
    dif0 = exp0 - mat0
    call if_null_c(dif0, stat0)
    !
    call bch_c(sx, sz, 1, mat1)
    exp1 = exp0 + (-2.d0*iota*hbar)*sy
    dif1 = exp1 - mat1
    call if_null_c(dif1, stat1)
    !
    call bch_c(sx, sz, 2, mat2)
    exp2 = exp1 - 0.5d0*((2.d0*iota*hbar)**2.d0)*sz
    dif2 = exp2 - mat2
    call if_null_c(dif2, stat2)
    !
    test_stat = .false.
    if (stat0 .and. stat1 .and. stat2) then
      test_stat = .true.
    end if
    !
  end subroutine test_bch_c
!> @brief subroutine for testing crank-nicolson evolution in qd_helper.f90
  subroutine test_crank_nicolson_evolution(test_stat)
    use global_m
    use qd_helper_m
    use ising_model_m
    use matrix_generator_m
    use linear_algebra_helper_m
    implicit none
    ! io variables
    logical, intent(out)                            :: test_stat
    ! internal variables
    double precision                                :: n_real
    integer                                         :: n
    double precision, parameter                     :: j=1.d0
    complex(8), parameter                           :: h=cmplx(0.d0,0.d0)
    complex(8), dimension(2,2)                      :: sz
    complex(8), dimension(2,2)                      :: sx
    complex(8), allocatable, dimension(:,:)         :: ham
    complex(8), allocatable, dimension(:,:)         :: eig_vect
    double precision, allocatable, dimension(:)     :: eig_val
    complex(8), allocatable, dimension(:)           :: psi0
    complex(8), allocatable, dimension(:)           :: psit
    complex(8)                                      :: norm0
    complex(8)                                      :: normt
    double precision                                :: dt
    integer                                         :: dimm
    logical                                         :: stat1
    logical                                         :: stat2
    complex(8), allocatable, dimension(:)           :: dif
    !
    call random_number(n_real)
    n = floor(n_real*10.d0)+1
    sz = pauli_matrices(3)
    sx = pauli_matrices(1)
    !
    dimm = 2**n
    ham = transverse_ising1d_hamiltonian(n, sz, sx, j, h, "o")
    !
    allocate(eig_vect(dimm,dimm))
    allocate(eig_val(dimm))
    allocate(psi0(dimm), psit(dimm))
    allocate(dif(dimm))
    !
    call diagonalize_matrix(2**n, ham, eig_vect, eig_val)
    psi0 = eig_vect(:,1)
    !
    call random_number(dt)
    psit = crank_nicolson_evolution(psi0, ham, dt)
    dif = psit-psi0
    if (sum(abs(dif)).le.tol) then
      stat1 = .false.
    end if
    !
    stat2 = .false.
    norm0=inner_product_dis(psi0,psi0)
    normt=inner_product_dis(psit,psit)
    if (abs(real(norm0)-real(normt)).le.tol) then
      stat2 = .true.
    end if
    !
    test_stat = .false.
    if (.not.stat1 .and. stat2) then
      test_stat = .true.
    end if
    !
    write(*,*) norm0
    write(*,*) normt
    !
    deallocate(ham, eig_vect, eig_val)
    deallocate(psi0, psit)
  end subroutine test_crank_nicolson_evolution
!> @brief subroutine for testing vi_truncation_index in simulate_helper.f90
  subroutine test_vi_truncation_index(test_stat)
    use global_m
    use simulate_helper_m
    use matrix_generator_m
    implicit none
    ! io variables
    logical, intent(out)                  :: test_stat
    ! internal variables
    integer, parameter                    :: n=5
    complex(8), dimension(2,2)            :: sz
    complex(8), dimension(2,2)            :: sx
    double precision, parameter           :: j=0.3d0
    complex(8), parameter                 :: h=cmplx(1.0d0,0.0d0)
    integer                               :: val
    !
    sx = (hbar/2.d0)*pauli_matrices(1)
    sz = (hbar/2.d0)*pauli_matrices(3)
    val = vi_truncation_index(n, j, h, sz, sx)
    !
    test_stat = .false.
    if (val .gt. 0) then
      test_stat = .true.
    end if
    !
    write(*,*) val
  end subroutine test_vi_truncation_index
end module test_driver_m
