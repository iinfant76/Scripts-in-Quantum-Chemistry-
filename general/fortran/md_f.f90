      SUBROUTINE make_bond_matrix_f(coords, n_atoms, r)
      integer*4 :: n_atoms, i_atm, j_atm, i_ax
!f2py intent(in) n_atoms
      real*8, dimension(n_atoms, 3)  :: coords
      real*8, dimension(n_atoms,n_atoms) :: r, d
!f2py intent(in) coords
!f2py intent(out) r 

      ! Initialize some variables 
      d = 0.0
      r = 0.0

      do i_atm = 1, n_atoms
         do j_atm = 1, n_atoms
            do i_ax = 1, 3  
               d(i_atm, j_atm) = d(i_atm, j_atm)+(coords(i_atm, i_ax) & 
                   - coords(j_atm, i_ax))**2
            enddo 
            r(i_atm, j_atm) = sqrt(d(i_atm, j_atm)) 
         enddo
      enddo

      END

      SUBROUTINE make_angle_matrix_f(r_mtx, n_atoms, teta)
      integer*4 :: n_atoms, i, j, k 
!f2py intent(in) n_atoms
      real*8, dimension(n_atoms, n_atoms) :: r_mtx 
      real*8, dimension(n_atoms,n_atoms,n_atoms) :: teta
!f2py intent(in) r_mtx
!f2py intent(out) teta  
      ! Initialize some variables 
      teta = 0.0

      do k = 1, n_atoms
        do j = 1, n_atoms
           do i = 1, n_atoms 
             if (i.ne.j .and. i.ne.k) then  
                teta(i,j,k) = acosd( ( r_mtx(i,j)**2 + r_mtx(j,k)**2   &
                    - r_mtx(i,k)**2 ) / ( 2 * r_mtx(i,j) * r_mtx(i,k)))  
             endif 
           enddo
        enddo
      enddo 
     
      end 


!def rdf(fn, atoms_i, atoms_j, dr, rmax, start, stride): 
    ! Retrieve number of lines in trajectory file
!    cmd = "wc -l {}".format(fn)
!    l = subprocess.check_output(cmd.split()).decode()
!    n_lines = int(l.split()[0]) ! Number of lines in traj file

    ! Read number of atoms in the molecule. It is usually the first line in a xyz file  
!    with open(fn) as f:
!       l = f.readline()
!       n_atoms = int(l.split()[0])

    ! Calculate total number of frames in the trajectory file  
!    n_frames = int(int(n_lines)/(n_atoms+2))
    
    ! Create grid of r values to compute the pair correlation function g_ij 
!    r_grid = np.arange(0, rmax, dr)
    
    ! Compute the volume in a concentric sphere of size dr at a distance r from the origin  
!    volume = 4 * np.pi * np.power(r_grid, 2) * dr ! elemental volume dV = 4*pi*r^2*dr
!    vol_tot = np.sum(volume) 
    
    ! Read the order of atoms from the first frame of the trajectory. This order is unchanged. 
!    atoms = np.genfromtxt(fn, skip_header = 2, skip_footer=(int(n_lines) - (n_atoms + 2)), usecols=0, dtype=str)

    ! Find the indexes in the bond_matrix of the atomic types involved in the calculation of g_ij
!    index_i = np.where( atoms == atoms_i )
!    index_j = np.where( atoms == atoms_j )

!    for iframe in range(start, n_frames, stride): 
        ! Read coordinates from iframe
!        coords = np.genfromtxt(fn, skip_header = (2 + (n_atoms + 2) * (iframe - 1)), 
!                    skip_footer=(int(n_lines) - ((n_atoms + 2) * iframe)), usecols=(1,2,3))
        ! Compute bond distance matrix for iframe
!        bond_mtx = make_bond_matrix(n_atoms, coords)
        ! Slice the bond_matrix with only the atom types
!        sliced_mtx = bond_mtx[np.ix_(index_i[0], index_j[0])]
        ! Count the number of atoms within r and r+dr  
!        rho_ab = np.stack( np.where( (sliced_mtx > r_grid[idr]) & (sliced_mtx < r_grid[idr+1]) )[0].size for idr in range(r_grid.size-1))
        ! Compute the total density of pair of atoms
!        rho_tot = np.sum(rho_ab) 
!        rho_bulk = rho_tot / vol_tot  
!        density = counts / volume[1:] ! skip the first element dV which is 0 
!        g_ab = rho_ab / rho_bulk 
        ! Store density in g_ij
!        if (iframe == start):
!            g_ij = g_ab
!        else:
!            g_ij = np.column_stack((g_ij,g_ab))

    ! Average over all frames
!    g_ij_av = np.average(g_ij, axis=1)

!    return r_grid[1:], g_ij_av 



