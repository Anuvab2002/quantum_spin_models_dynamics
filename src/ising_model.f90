!> @file ising_model.f90
!> @brief provides different routines for quantum isning model calculations
!> @author ap
module ising_model_m
contains
!> @brief function for calculating the spin-spin interaction terms of the Hamiltonian
!> @param[in]      s            spin matrix of given direction
!> @param[in]      n            size of the chain
!> @param[in]      bc           boundary condition
!> @param[in]      j            coupling strength
!> @return         spin_int     spin-spin interaction term
  function spin_spin_interaction(s, n, boundary, j)result(spin_int)
    use linear_algebra_helper_m
    use matrix_generator_m
    implicit none
    ! io variables
    complex(8), dimension(:,:), intent(in)    :: s
    integer, intent(in)                       :: n
    character(len=1), intent(in)              :: boundary
    double precision, intent(in)              :: j
    complex(8), allocatable, dimension(:,:)   :: spin_int
    ! internal variable
    integer                                   :: dim
    integer                                   :: dimm
    complex(8), allocatable, dimension(:,:)   :: op_s
    integer                                   :: idx
    complex(8), allocatable, dimension(:,:)   :: i1
    complex(8), allocatable, dimension(:,:)   :: i2
    integer                                   :: p
    integer                                   :: q
    complex(8), allocatable, dimension(:,:)   :: temp1
    complex(8), allocatable, dimension(:,:)   :: temp2
    complex(8), allocatable, dimension(:,:)   :: spin_int_obc
    complex(8), allocatable, dimension(:,:)   :: spin_int_pbc
    !
    dim = size(s,1)
    if (size(s,2).ne.dim) then
      spin_int = cmplx(0.d0,0.d0)
      stop "Execution error! Spin-spin interaction calculation failed! Provid a square spin matrix."
    end if
    dim = size(s,1)
    dimm = dim**n
    !
    allocate(spin_int(dimm, dimm))
    spin_int = cmplx(0.d0, 0.d0)
    allocate(spin_int_obc(dim**n,dim**n))
    allocate(spin_int_pbc(dim**n,dim**n))
    spin_int_obc = cmplx(0.d0, 0.d0)
    spin_int_pbc = cmplx(0.d0, 0.d0)
    !
    allocate(op_s(dim**2,dim**2))
    call kron_product(s, s, op_s)
    !
    allocate(temp2(dimm,dimm))
    do idx = 1,n-1
      temp2 = cmplx(0.d0, 0.d0)
      !
      p = idx-1
      q = n-idx-1
      i1 = identity_matrix_complex(dim**p)
      i2 = identity_matrix_complex(dim**q)
      allocate(temp1(dim**(idx+1),dim**(idx+1)))
      call kron_product(i1, op_s, temp1)
      call kron_product(temp1, i2, temp2)
      deallocate(i1,i2)
      deallocate(temp1)
      spin_int_obc = spin_int_obc + j*temp2
    end do
    !
  select case(boundary)
  case("o")
    spin_int = spin_int_obc
    deallocate(temp2)
  case("p")
    temp2 = cmplx(0.0d0, 0.0d0)
    allocate(temp1(dim**(n-1),dim**(n-1)))
    i1 = identity_matrix_complex(dim**(n-2))
    call kron_product(s, i1, temp1)
    call kron_product(temp1, s, temp2)
    spin_int_pbc = spin_int_obc + j*temp2
    deallocate(i1)
    deallocate(temp1, temp2)
    spin_int = spin_int_pbc
  end select
  deallocate(spin_int_obc, spin_int_pbc)
  end function spin_spin_interaction
!> @brief function for calculating the zeeman term of the Hamiltonian
!> @param[in]       s           the spin matrix
!> @param[in]       n           number of spin sites in the chain
!> @param[in]       h           strenght of the spin-field coupling
!> return           field_int   the zeeman term
  function zeeman_term(s, n, h)result(field_int)
    use linear_algebra_helper_m
    use matrix_generator_m
    implicit none
    ! io variables
    complex(8), dimension(:,:), intent(in)      :: s
    integer, intent(in)                         :: n
    complex(8), intent(in)                      :: h
    complex, allocatable, dimension(:,:)        :: field_int
    ! internal variables
    integer                                     :: dim
    integer                                     :: idx
    integer                                     :: p
    integer                                     :: q
    complex(8), allocatable, dimension(:,:)     :: i1
    complex(8), allocatable, dimension(:,:)     :: i2
    complex(8), allocatable, dimension(:,:)     :: temp1
    complex(8), allocatable, dimension(:,:)     :: temp2
    !
    dim = size(s,1)
    if (size(s,2).ne.dim) then
      field_int = cmplx(0.d0,0.d0)
      stop "Execution error! Zeeman term calculation failed! Provid a square spin matrix."
    end if
    allocate(field_int(dim**n,dim**n))
    field_int = cmplx(0.d0, 0.d0)
    !
    allocate(temp2(dim**n,dim**n))
    do idx = 1,n
      temp2 = cmplx(0.d0, 0.d0)
      !
      p = idx-1
      q = n-idx
      !
      i1 = identity_matrix_complex(dim**p)
      i2 = identity_matrix_complex(dim**q)
      !
      allocate(temp1(dim**(idx),dim**(idx)))
      call kron_product(i1, s, temp1)
      call kron_product(temp1, i2, temp2)
      !
      deallocate(i1, i2, temp1)
      field_int = field_int + h*temp2
    end do
    deallocate(temp2)
  end function zeeman_term
!> @brief function for calculating magnetization operator
!> @param[in]       s         spin matrix
!> @param[in]       n         number of spin chains
!> @return           mag       magnetization operator
  function magnetization_operator(s, n)result(mag)
    use linear_algebra_helper_m
    use matrix_generator_m
    implicit none
    ! io variables
    complex(8), dimension(:,:), intent(in)      :: s
    integer, intent(in)                         :: n
    complex(8), allocatable, dimension(:,:)     :: mag
    ! interal variables
    integer                                     :: dim
    integer                                     :: idx
    integer                                     :: p
    integer                                     :: q
    complex(8), allocatable, dimension(:,:)     :: i1
    complex(8), allocatable, dimension(:,:)     :: i2
    complex(8), allocatable, dimension(:,:)     :: temp1
    complex(8), allocatable, dimension(:,:)     :: temp2
    !
    dim = size(s,1)
    if (size(s,2).ne.dim) then
      mag = cmplx(0.d0,0.d0)
      stop "Execution error! Magnetization operator calculation failed. Provid a square spin matrix."
    end if
    !
    allocate(mag(dim**n,dim**n))
    mag = cmplx(1.d0, 0.d0)
    !
    allocate(temp2(dim**n,dim**n))
    do idx = 1,n-1
      temp2 = cmplx(0.d0, 0.d0)
      !
      p = idx-1
      q = n-idx
      !
      i1 = identity_matrix_complex(dim**p)
      i2 = identity_matrix_complex(dim**q)
      !
      allocate(temp1(dim**(idx),dim**(idx)))
      call kron_product(i1, s, temp1)
      call kron_product(temp1, i2, temp2)
      !
      deallocate(i1, i2, temp1)
      mag = mag + temp2
    end do
    deallocate(temp2)
  end function magnetization_operator
end module ising_model_m
