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
  character(len=1), parameter           :: bc="p"! boundary condition: o(open)/p(perodic)
  double precision, parameter           :: jz=1.d0! z-directional spin-spin coupling strength
  double precision, parameter           :: jx=1.d0! x-directional spin-spin coupling strength
  double precision, parameter           :: jy=1.d0! y-directional spin-spin coupling strength
  complex(8), parameter                 :: h0=cmplx(0.5d0,0.d0)! constant parameter for spin-field coupling strength
!---------------------------- Dynamical parameters ----------------------------!
  integer, parameter                    :: n_bch=10! order of BCH calculation
  double precision, parameter           :: tfinal=1.d0! total time of dynamics
  double precision, parameter           :: step=0.001d0! size of the time-step
  integer, parameter                    :: ramp_cut=5! cut-off number for ramp function
  integer, parameter                    :: nosc=50! total number of oscillation
  double precision, parameter           :: freq=2.d0*pi*nosc/tfinal! frequency of oscillation of spin-field coupling constant
end module global_m
