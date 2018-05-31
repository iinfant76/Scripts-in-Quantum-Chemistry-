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

      SUBROUTINE make_angle_matrix_f(r_mtx,coords,n_atoms,bond_max1, &
                                     bond_max2, teta)
      integer*4 :: n_atoms, i, j, k, axes 
!f2py intent(in) n_atoms
      real*8 :: bond_max1, bond_max2
      real*8, dimension(n_atoms, 3)  :: coords 
      real*8, dimension(n_atoms, n_atoms) :: r_mtx 
      real*8, dimension(3) :: uij, vjk
      real*8, dimension(n_atoms,n_atoms,n_atoms) :: teta
!f2py intent(in) coords 
!f2py intent(in) r_mtx
!f2py intent(in) bond_max1
!f2py intent(in) bond_max2 
!f2py intent(out) teta  
      ! Initialize some variables 
      teta = 0.0
      d = 0.0 
      bond_min = 0.5
      do k = 1, n_atoms
        do j = 1, n_atoms
           if (r_mtx(j, k).gt.bond_min .and.   & 
                     r_mtx(j, k).lt.bond_max2 ) then 
              do i = 1, n_atoms 
                 if ( r_mtx(i, j).gt.bond_min .and. & 
                                  r_mtx(i, j).lt.bond_max1) then 
                    if (i.ne.j .and. i.lt.k .and. j.ne.k) then  
                       do axes = 1, 3
                          uij(axes) = coords(i, axes) - coords(j, axes)
                          vjk(axes) = coords(j, axes) - coords(k, axes) 
                       enddo
                       teta(i,j,k) = acos( dot_product(uij, vjk) /  & 
                              ( r_mtx(i,j) * r_mtx(j,k)))  
                    endif 
                 endif 
              enddo
           endif 
        enddo
      enddo 
     
      end 

