program H2_disp_v2
   use iso_fortran_env
   implicit none

   integer,  parameter :: dp    = real64
   real(dp), parameter :: pi    = 3.14159265358979323846264338327950288419_dp
   
   ! Physical constants
   real(dp), parameter :: kB_SI       = 1.380649e-23_dp      ! J/K
   real(dp), parameter :: hbar_SI     = 1.054571817e-34_dp   ! J s
   real(dp), parameter :: eV_to_J     = 1.602176634e-19_dp   ! J/eV
   real(dp), parameter :: amu_to_kg   = 1.66054e-27_dp 
   real(dp), parameter :: mass_H2_amu = 2.016_dp
   real(dp), parameter :: mass_H2_kg  = mass_H2_amu * amu_to_kg  
   
   ! Supercell parameters for phonon calculation
   integer,  parameter :: nc      = 7
   integer,  parameter :: n_basis = 2                   ! 2 atoms per unit cell (HCP)
   integer,  parameter :: n_unit  = nc**3               ! Number of unit cells = 343
   integer,  parameter :: n       = n_unit * n_basis    ! Total atoms = 686 
   
   ! Lattice parameters (Angstrom)
   real(dp), parameter :: lattice_param = 3.79_dp      ! Lattice constant a
   real(dp), parameter :: c_lat         = 6.19_dp      ! c lattice parameter
   real(dp), parameter :: disp_step     = 0.0001_dp    ! Finite displacement epsilon
   
   real(dp), parameter :: sg_alpha  = 1.713_dp    
   real(dp), parameter :: sg_beta   = 2.9615_dp      ! In 1/Angstrom
   real(dp), parameter :: sg_gamma  = 0.03546_dp   
   real(dp), parameter :: sg_rc     = 4.403_dp       ! Cutoff for damping function in Angstrom
   real(dp), parameter :: sg_C6     = 84178.9_dp     ! In K*Angstrom^6
   real(dp), parameter :: sg_C8     = 417858.4_dp  
   real(dp), parameter :: sg_C9     = 147037.3_dp   
   real(dp), parameter :: sg_C10    = 2617496.9_dp 
   real(dp), parameter :: cutoff    = 12.0_dp        ! Cutoff for intermolecular interaction

   integer,  parameter :: n1_seg   = 64                 ! Number of q points along path in reciprocal space
   integer,  parameter :: n2_seg   = 32
   integer,  parameter :: n3_seg   = 56
   integer,  parameter :: n4_seg   = 28 
   integer,  parameter :: n_points = n1_seg + n2_seg + n3_seg + n4_seg

   real(dp)    :: r0(3,n)                            ! Equilibrium positions
   real(dp)    :: r(3,n)                             ! Current positions (for displacements)
   real(dp)    :: f(3,n), f_plus(3,n), f_minus(3,n) ! Forces
   real(dp)    :: phi(3, 3, n_basis, n_basis, n_unit) ! Force constant matrix
   real(dp)    :: b_vec(3,3)                         ! Reciprocal lattice vectors
   integer     :: i, j, k

   call setup_structure()
   call calculate_phi_matrix()
   call calculate_dispersion()
   print *, 'Results saved to dispersion2.dat'


contains

   subroutine setup_structure()
      ! Setup equilibrium HCP structure and reciprocal lattice vectors
      real(dp) :: avec(3,3), basis(3,2)
      integer :: ix, iy, iz, idx, ib
      
      ! Primitive lattice vectors (Equations 2-4 from my note)
      avec(:,1) = [lattice_param, 0.0_dp, 0.0_dp]
      avec(:,2) = [0.5_dp*lattice_param, sqrt(3.0_dp)/2.0_dp*lattice_param, 0.0_dp]
      avec(:,3) = [0.0_dp, 0.0_dp, c_lat]
      
      ! Reciprocal lattice vectors (Equations 6-8 from my note)
      ! b1 = (2pi / sqrt(3)a) * (sqrt(3)x - y)
      b_vec(:,1) = (2.0_dp*pi / (sqrt(3.0_dp)*lattice_param)) * &
                   [sqrt(3.0_dp), -1.0_dp, 0.0_dp]
      ! b2 = (4pi / sqrt(3)a) * y
      b_vec(:,2) = (4.0_dp*pi / (sqrt(3.0_dp)*lattice_param)) * &
                   [0.0_dp, 1.0_dp, 0.0_dp]
      ! b3 = 2pi/c * z
      b_vec(:,3) = [0.0_dp, 0.0_dp, 2.0_dp*pi/c_lat]
      
      ! Basis atoms (positions within unit cell)
      ! Atom 1 at (0,0,0)
      basis(:,1) = [0.0_dp, 0.0_dp, 0.0_dp]
      ! Atom 2 at (1/3, 1/3, 1/2) in fractional coordinates
      basis(:,2) = (avec(:,1) + avec(:,2))/3.0_dp + avec(:,3)/2.0_dp
      
      ! Build supercell positions
      idx = 0
      do iz = 0, nc-1
         do iy = 0, nc-1
            do ix = 0, nc-1
               do ib = 1, n_basis
                  idx = idx + 1
                  r0(:, idx) = real(ix,dp)*avec(:,1) + &
                               real(iy,dp)*avec(:,2) + &
                               real(iz,dp)*avec(:,3) + &
                               basis(:,ib)
               end do
            end do
         end do
      end do
      
   end subroutine setup_structure

   subroutine construct_dynamical_matrix(q_cart, Dmat)
      ! Construct dynamical matrix D(q) using Equation (12) from my note
      ! D_ka,k'a'(q) = (1/sqrt(m*m)) * Sum_l Phi_ka,k'a'(0,l) * exp(i*q*R_l)
      
      real(dp), intent(in) :: q_cart(3)           ! q-vector in Cartesian coordinates
      complex(dp), intent(out) :: Dmat(6, 6)      ! 6x6 dynamical matrix
      
      integer :: l, k1, k2, a1, a2, row, col
      integer :: ix_c, iy_c, iz_c, central_cell
      integer :: i_source, i_target
      real(dp) :: r_diff(3), phase
      complex(dp) :: c_phase, im_unit
      
      im_unit = cmplx(0.0_dp, 1.0_dp, kind=dp)
      Dmat = cmplx(0.0_dp, 0.0_dp, kind=dp)
      
      ! Find central unit cell (same as in calculate_phi_matrix)
      ix_c = nc / 2
      iy_c = nc / 2
      iz_c = nc / 2
      central_cell = iz_c * nc * nc + iy_c * nc + ix_c + 1
      
      ! Loop over all unit cells in the supercell
      do l = 1, n_unit
         do k1 = 1, n_basis  ! Source atom (kappa) - in central cell
            do k2 = 1, n_basis  ! Target atom (kappa') - in cell l
               
               ! Calculate position difference R_{l,k2} - R_{central,k1}
               i_source = (central_cell-1)*n_basis + k1  ! Source in central cell
               i_target = (l-1)*n_basis + k2             ! Target in cell l
               r_diff = r0(:, i_target) - r0(:, i_source)
               
               ! Calculate phase factor exp(i*q*R)
               phase = dot_product(q_cart, r_diff)
               c_phase = exp(im_unit * phase)
               
               ! Build dynamical matrix elements
               do a1 = 1, 3  ! alpha (row)
                  do a2 = 1, 3  ! alpha' (col)
                     row = (k1-1)*3 + a1
                     col = (k2-1)*3 + a2
                     
                     ! D = (1/m) * Sum( Phi * exp(iqR) )
                     Dmat(row, col) = Dmat(row, col) + &
                        (phi(a1, a2, k2, k1, l) / mass_H2_kg) * c_phase
                  end do
               end do
            end do
         end do
      end do
      
   end subroutine construct_dynamical_matrix

   subroutine calculate_dispersion()
      ! Calculate phonon dispersion along Γ-K-M-Γ-A path
      
      complex(dp), allocatable :: Dmat(:,:)
      real(dp), allocatable :: w_eig(:)
      real(dp), allocatable :: q_path(:,:), dist(:)
      integer, parameter :: n_segments = 4
      real(dp) :: q_frac(3), q_cart(3)
      integer :: ip, j, info
      complex(dp), allocatable :: work(:)
      real(dp), allocatable :: rwork(:)
      integer :: lwork
      
      allocate(Dmat(6,6), w_eig(6))
      allocate(q_path(3, n_points), dist(n_points))

      call generate_path(q_path, dist)

      open(unit=20, file='dispersion2.dat', status='replace')
      write(20, '(A)') '# Distance  Energy_meV(1..6)'

      do ip = 1, n_points
         if (mod(ip, 50) == 0 .or. ip == 1) then
            print *, 'q-point ', ip, '/', n_points
         end if
         
         q_frac = q_path(:, ip)
         
         ! Convert fractional q to Cartesian
         q_cart = q_frac(1)*b_vec(:,1) + q_frac(2)*b_vec(:,2) + q_frac(3)*b_vec(:,3)

         call construct_dynamical_matrix(q_cart, Dmat)

         lwork = 20
         allocate(work(lwork), rwork(3*6-2))
         call zheev('N', 'U', 6, Dmat, 6, w_eig, work, lwork, rwork, info)
         deallocate(work, rwork)
         
         if (info /= 0) then
            print *, 'ERROR: Diagonalisation failed at q-point', ip, 'info=', info
            stop
         end if
         
         ! Convert eigenvalues to phonon energies in meV
         ! w_eig has units [K/Å²]/[kg] from D = Phi/mass 
         ! ω^2 [rad^2/s^2] = w_eig x k_B[J/K] / (10^-20 m^2/Å^2)
         ! ω [rad/s] = sqrt(w_eig x k_B x 10^20)
         ! E [meV] = hbar x ω x 1000 / eV_to_J
         write(20, '(F10.4, 6ES16.6)') dist(ip), &
            (hbar_SI * sqrt(max(w_eig(j), 0.0_dp) * kB_SI * 1.0e20_dp) * 1000.0_dp / eV_to_J, j=1,6)
      end do
      
      close(20)
      deallocate(Dmat, w_eig, q_path, dist)
      
   end subroutine calculate_dispersion

! ------------------------------------------------------------------
    ! Arguments:
    !   m   (in)    : Dimension of the matrix
    !   A   (inout) : The matrix to diagonalise. 
    !   w   (out)   : Vector containing the m eigenvalues (sorted lowest to highest)
    ! ------------------------------------------------------------------
   subroutine diag_hermitian(m, A, w)
      integer, intent(in)        :: m
      complex(dp), intent(inout) :: A(m, m)
      real(dp), intent(out)      :: w(m)

      integer :: info, lwork
      complex(dp), allocatable :: work(:)
      real(dp), allocatable    :: rwork(:)
      complex(dp) :: work_query(1)

      allocate(rwork(3*m - 2))
        
      call zheev('N', 'U', m, A, m, w, work_query, -1, rwork, info)

      if (info /= 0) then
         print *, "Error: Workspace query failed in ZHEEV."
         stop
      end if

      lwork = int(real(work_query(1)))
      allocate(work(lwork))

      call zheev('N', 'U', m, A, m, w, work, lwork, rwork, info)

      if (info /= 0) then
         print *, "Error: Diagonalisation failed! Info =", info
         stop
      end if
        
   end subroutine diag_hermitian

   subroutine generate_path(path, d_out)
      ! Generate k-point path for HCP Brillouin zone
      ! Path: Γ(0,0,0) -> K(2/3,1/3,0) -> M(1/2,0,0) -> Γ -> A(0,0,1/2)
      ! Following Equations (15-19) from note
      
      real(dp), intent(out) :: path(:,:)    ! (3, n_total) array
      real(dp), intent(out) :: d_out(:)     ! Distance along path
      
      real(dp) :: k_G(3), k_K(3), k_M(3), k_A(3)
      integer :: i, n
      
      ! High symmetry points in fractional coordinates
      k_G = [0.0_dp, 0.0_dp, 0.0_dp]        ! Γ
      k_K = [2.0_dp/3.0_dp, 1.0_dp/3.0_dp, 0.0_dp]  ! K
      k_M = [0.5_dp, 0.0_dp, 0.0_dp]        ! M
      k_A = [0.0_dp, 0.0_dp, 0.5_dp]        ! A
      
      n = 0
      d_out(1) = 0.0_dp
      
      ! Segment 1: Gamma -> K (Λ line)
      do i = 1, n1_seg
         n = n + 1
         path(:,n) = k_G + (real(i-1,dp)/real(n1_seg-1,dp)) * (k_K - k_G)
         if (n > 1) d_out(n) = d_out(n-1) + norm2(path(:,n)-path(:,n-1))
      end do
      
      ! Segment 2: K -> M (T line)
      do i = 1, n2_seg
         n = n + 1
         path(:,n) = k_K + (real(i,dp)/real(n2_seg,dp)) * (k_M - k_K)
         d_out(n) = d_out(n-1) + norm2(path(:,n)-path(:,n-1))
      end do
      
      ! Segment 3: M -> Γ (Σ line)
      do i = 1, n3_seg
         n = n + 1
         path(:,n) = k_M + (real(i,dp)/real(n3_seg,dp)) * (k_G - k_M)
         d_out(n) = d_out(n-1) + norm2(path(:,n)-path(:,n-1))
      end do
      
      ! Segment 4: Γ-> A (Δ line)
      do i = 1, n4_seg
         n = n + 1
         path(:,n) = k_G + (real(i,dp)/real(n4_seg,dp)) * (k_A - k_G)
         d_out(n) = d_out(n-1) + norm2(path(:,n)-path(:,n-1))
      end do
      
   end subroutine generate_path

   subroutine calculate_phi_matrix()
      ! Calculate force constant matrix Phi using finite displacement method
      ! Phi(beta, alpha, k_target, k_source, l_unit) = -(F_plus - F_minus) / (2*epsilon)
      ! This implements Equation (14) from my note
      
      integer :: k_source, alpha, i_target, k_target, l_idx
      integer :: ix_c, iy_c, iz_c, central_cell, i_source
      
      phi = 0.0_dp
      
      ! Find central unit cell
      ix_c = nc / 2
      iy_c = nc / 2
      iz_c = nc / 2
      central_cell = iz_c * nc * nc + iy_c * nc + ix_c + 1
      
      print *, 'Displacing atoms in central unit cell:', central_cell
      print *, 'Central cell position: (', ix_c, ',', iy_c, ',', iz_c, ')'
      
      ! Loop over basis atoms in the CENTRAL unit cell
      do k_source = 1, n_basis
         ! Global index of source atom in central cell
         i_source = (central_cell - 1) * n_basis + k_source
         
         do alpha = 1, 3
            r = r0
            r(alpha, i_source) = r(alpha, i_source) + disp_step
            call force_vpimd(r, f_plus)
            r = r0
            r(alpha, i_source) = r(alpha, i_source) - disp_step
            call force_vpimd(r, f_minus)
            
            ! Calculate force constants for all target atoms
            do i_target = 1, n
               ! Map linear index to (unit cell, basis atom)
               k_target = mod(i_target - 1, n_basis) + 1
               l_idx = (i_target - 1) / n_basis + 1
               
               ! Phi_beta,alpha (k_target, k_source, l_idx) = -(F+ - F-) / (2*eps)
               phi(:, alpha, k_target, k_source, l_idx) = &
                  -(f_plus(:, i_target) - f_minus(:, i_target)) / (2.0_dp * disp_step)
            end do
         end do
      end do
      
   end subroutine calculate_phi_matrix

   subroutine force_vpimd(r, f)
      real(dp), intent(in)  :: r(:,:)
      real(dp), intent(out) :: f(:,:)
      real(dp), parameter   :: T_ref = 5.4_dp 
      
      ! Persistent Data (Saved between calls)
      real(dp), allocatable, save :: coeffs(:,:) 
      real(dp), allocatable, save :: R_val(:)
      real(dp), save :: r_min, dr, dr_inv, r_cut_sq
      integer, save :: n_bins
      logical, save :: initialised = .false.

      integer :: n_atoms, i, j, idx, unit_num, ios, n_lines
      real(dp) :: rij(3), r2, r_mag, delta, f_mag
      real(dp) :: c1, c2, c3
      character(len=128) :: line
      real(dp), allocatable :: V_dat(:), V_d2(:), u(:)
      real(dp) :: val_r, val_v
      real(dp) :: yp0, ypn, sig, p, h
      real(dp) :: ya, yb, y2a, y2b

      n_atoms = size(r, 2)
      f = 0.0_dp

      ! ----------------------------------------------------------------
      ! PART 1: INITIALISATION (Executed only on the first call)
      ! ----------------------------------------------------------------
      if (.not. initialised) then
         print *, '--- Initialising Potential from Vpimd.dat ---'
          
         unit_num = 10
         open(unit=unit_num, file='Vpimd.dat', status='old', iostat=ios)
         if (ios /= 0) stop 'Error: Could not open Vpimd.dat'

         ! Count valid data lines
         n_lines = 0
         do
            read(unit_num, '(A)', iostat=ios) line
            if (ios /= 0) exit
            if (len_trim(line) == 0 .or. line(1:1)=='#' .or. line(1:1)=='@') cycle
            n_lines = n_lines + 1
         end do
         rewind(unit_num)
          
         n_bins = n_lines - 1
         allocate(R_val(0:n_bins), V_dat(0:n_bins))
         allocate(V_d2(0:n_bins), u(0:n_bins))
         allocate(coeffs(0:n_bins, 4))
          
         ! Read data and convert units
         i = 0
         do
            read(unit_num, '(A)', iostat=ios) line
            if (ios /= 0) exit
            if (len_trim(line) == 0 .or. line(1:1)=='#' .or. line(1:1)=='@') cycle ! Check for empty and non data lines
             
            read(line, *) val_r, val_v
             
            R_val(i) = val_r
            V_dat(i) = val_v * T_ref ! Convert (kB*5.4K) -> Kelvin
             
            i = i + 1
         end do
         close(unit_num)

         r_min = R_val(0)
         dr = R_val(1) - R_val(0)
         dr_inv = 1.0_dp / dr
         r_cut_sq = R_val(n_bins)**2

         print *, '  R_min (A):', r_min
         print *, '  dR (A):', dr
         print *, '  R_cut (A):', sqrt(r_cut_sq)

         ! Taken from numerical recipe page 109 for cubic spline
         V_d2(0) = 0.0_dp; u(0) = 0.0_dp
          
         do i = 1, n_bins-1
            sig = 0.5_dp
            p = sig * V_d2(i-1) + 2.0_dp
            V_d2(i) = (sig - 1.0_dp) / p
            u(i) = (V_dat(i+1) - V_dat(i)) / dr - &
                   (V_dat(i) - V_dat(i-1)) / dr
            u(i) = (6.0_dp * u(i) / (2.0_dp*dr) - sig*u(i-1)) / p
         end do
          
         V_d2(n_bins) = 0.0_dp
         do i = n_bins-1, 0, -1
            V_d2(i) = V_d2(i) * V_d2(i+1) + u(i)
         end do
          
         ! Now we splint
         h = dr
         do i = 0, n_bins-1
            ya = V_dat(i); yb = V_dat(i+1)
            y2a = V_d2(i); y2b = V_d2(i+1)
             
            coeffs(i, 1) = ya
            coeffs(i, 2) = (yb-ya)/h - (h/6.0_dp)*(2.0_dp*y2a+y2b) 
            coeffs(i, 3) = y2a / 2.0_dp
            coeffs(i, 4) = (y2b-y2a) / (6.0_dp*h) ! See my second note for details
         end do
         coeffs(n_bins, :) = 0.0_dp 
          
         deallocate(V_dat, V_d2, u)
         initialised = .true.
      end if

      ! ----------------------------------------------------------------
      ! PART 2: FORCE CALCULATION
      ! ----------------------------------------------------------------
      f = 0.0_dp

      do i = 1, n_atoms
         do j = i+1, n_atoms
            rij = r(:,j) - r(:,i) ! vector direction i -> j
            r2 = dot_product(rij, rij)
            
            if (r2 < r_cut_sq .and. r2 > 1.0e-8_dp) then
               r_mag = sqrt(r2)
               
               ! Find bin index
               idx = int((r_mag - r_min) * dr_inv)
               
               if (idx >= 0 .and. idx < n_bins) then
                  delta = r_mag - (r_min + real(idx, dp) * dr)
                   
                  c1 = coeffs(idx, 2) ! Linear coefficient
                  c2 = coeffs(idx, 3) ! Quadratic coefficient
                  c3 = coeffs(idx, 4) ! Cubic coefficient
                   
                  ! Force Magnitude F = dV/dr (applied along rij vector)
                  ! Force on i = -grad_i V = -V' * grad_i r = -V' * (-rij/r) = V' * rij/r
                  f_mag = (c1 + delta * (2.0_dp*c2 + delta * 3.0_dp*c3))
                   
                  ! F_vec = (F_mag / r) * r_vec
                  f_mag = f_mag / r_mag
                  rij = rij * f_mag
                   
                  f(:,i) = f(:,i) + rij
                  f(:,j) = f(:,j) - rij
               end if
            end if
         end do
      end do
   end subroutine force_vpimd

end program H2_disp_v2








