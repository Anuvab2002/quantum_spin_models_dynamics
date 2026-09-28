!> @file global.f90
!> @brief provides global parameters
!> @author ap
module global_m
  double precision, parameter           :: zero=1.0d-8! zero
  double precision, parameter           :: tol=1.0d-5! tolerance
  complex(8), parameter                 :: iota=(0.d0,1.d0)! iota
  double precision, parameter           :: pi=4.d0*atan(1.d0)! pi
  double precision, parameter           :: hbar=1.d0! hbar
!--------------------------- Spin chain parameters ----------------------------!
  integer, parameter                    :: site=3! number of spin sites in the chain
  character(len=1), parameter           :: bc="o"! boundary condition: o(open)/p(perodic)
  double precision, parameter           :: j_z=1.d0! z-directional spin-spin coupling strength
  double precision, parameter           :: j_x=1.d0! x-directional spin-spin coupling strength
  double precision, parameter           :: j_y=1.d0! y-directional spin-spin coupling strength
!---------------------------- Dynamical parameters ----------------------------!
  integer, parameter                    :: n_bch=5! order of BCH calculation

end module global_m
