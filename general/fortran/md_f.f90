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

