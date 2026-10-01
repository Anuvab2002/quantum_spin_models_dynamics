!> @file run_driver.f90
!> @brief provides routines for final time evolution calculations
!> @author ap
module run_driver_m
contains
!> @brief subroutine for time evolution of 1-d transverse field isng model
  subroutine time_evol_transverse_ising_1d()
    use global_m
    use qd_helper_m
    use dqc_probe_m
    use ising_model_m
    use simulate_helper_m
    use matrix_generator_m
    use linear_algebra_helper_m
    implicit none
    ! internal variables
    double precision                    :: t
    complex(8)                          :: ht
    !
    integer, parameter                  :: dimm=2**site
    !
    complex(8), dimension(2,2)          :: sx
    complex(8), dimension(2,2)          :: sy
    complex(8), dimension(2,2)          :: sz
    !
    complex(8), dimension(dimm,dimm)    :: ham0
    complex(8), dimension(dimm,dimm)    :: v_s
    complex(8), dimension(dimm,dimm)    :: v_i
    complex(8), dimension(dimm,dimm)    :: eig_vect
    double precision, dimension(dimm)   :: eig_val
    complex(8), dimension(dimm,dimm)    :: mag_s
    complex(8), dimension(dimm,dimm)    :: mag_i
    complex(8), dimension(dimm,dimm)    :: dmagdt_i
    !
    complex(8), dimension(dimm)         :: psi_old
    complex(8), dimension(dimm)         :: psi0
    complex(8), dimension(dimm)         :: psi_new
    !
    integer                             :: idx
    integer                             :: ntime
    !
    double precision                    :: norm
    double precision                    :: norm0
    double precision                    :: le1
    complex(8)                          :: mag_val
    complex(8)                          :: dmagdt_val
    !
    sx = (hbar*0.5d0)*pauli_matrices(1)
    sy = (hbar*0.5d0)*pauli_matrices(2)
    sz = (hbar*0.5d0)*pauli_matrices(3)
    !
    t = 0.d0
    ham0 = spin_spin_interaction(sz, site, bc, jz)
    call diagonalize_matrix(dimm, ham0, eig_vect, eig_val)
    psi_old = eig_vect(:,dimm)
    psi0 = psi_old
    norm0 = abs(inner_product_dis(psi0, psi0))
    !
    mag_s = magnetization_operator(sz, site)
    !
    open(unit=100, file="ht.dat", status="replace", action="write")
    open(unit=101, file="norm_3spinp_1d.dat", status="replace", action="write")
    open(unit=102, file="le1_3spinp_1d.dat", status="replace", action="write")
    open(unit=104, file="mag_3spinp_1d.dat", status="replace", action="write")
    open(unit=105, file="dmagdt_3spinp_1d.dat", status="replace", action="write")
    open(unit=106, file="mag_dmagdt_3spinp_1d.dat", status="replace", action="write")
    ntime = floor(tfinal/step)
    do idx = 1,ntime
      v_s = cmplx(0.d0,0.d0)
      v_i = cmplx(0.d0,0.d0)
      mag_i = cmplx(0.d0,0.d0)
      psi_new = cmplx(0.d0,0.d0)
      !
      t = t+step
      !
      ht = h_t(h0, freq, t, ramp_cut)
      v_s = zeeman_term(sx, site, ht)
      v_i = operator_interaction(v_s, ham0, t, n_bch)
      psi_new = interaction_wavefunction(psi_old, ham0, v_s, t, step, n_bch)
      !
      norm = abs(inner_product_dis(psi_new, psi_new))
      !
      le1 = le(psi_new, psi0)
      !
      mag_i = operator_interaction(mag_s, ham0, t, n_bch)
      mag_val = expectation_value_dis(psi_new, mag_i)
      !
      dmagdt_i = heisenberg_eom(mag_i, ham0)
      dmagdt_val = expectation_value_dis(psi_new, dmagdt_i)
      !
      psi_old = psi_new
      write(100,*)  t, real(ht)
      write(101,*)  t, (norm-norm0)*(10.d0**12.d0)
      write(102,*)  t, le1
      write(104,*)  t, real(mag_val)
      write(105,*)  t, real(dmagdt_val)
      write(106,*)  real(mag_val), real(dmagdt_val)
    end do
    !
    close(100)
    close(101)
    close(102)
    close(104)
    close(105)
  end subroutine time_evol_transverse_ising_1d
end module run_driver_m
