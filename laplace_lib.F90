module laplace_lib
implicit none

private
public laplace, jacobian

contains

subroutine laplace(inl, jnl, x, side_x, r, side_r, phi, nx, nr, nx_loc, f_temp)
    integer, intent(in):: inl, jnl, side_x, nx, side_r, nr, nx_loc
    real(8), intent(in):: x(nx), r(nr), phi(nr,nx_loc)
    real(8), intent(out):: f_temp
    real(8):: dx_p = 0.d0, dx_m = 0.d0
    real(8):: dr_p = 0.d0, dr_m = 0.d0
    real(8):: dphi_dx = 0.d0, dphi_dr = 0.d0
    integer:: i

!-----------------------------------------------------------------------
!*******************POISSON equation***********************
!-----------------------------------------------------------------------
    if (nx > 1) then
        ! X-dir Center
        if (side_x .eq. 0) then
            dx_p  = x(jnl+1) - x(jnl)
            dx_m = x(jnl)   - x(jnl-1)

            ! term: d^2(phi)/dx^2 (x-term)
            dphi_dx = 2.d0*phi(inl,jnl+1)/(dx_p &
            * (dx_p+dx_m)) &
            - 2.d0*phi(inl,jnl)/(dx_p*dx_m) &
            + 2.d0*phi(inl,jnl-1)/(dx_m &
            * (dx_p+dx_m))

        ! X-dir left (vacuum)
        else if (side_x < 0) then
            dx_p = x(jnl+1) - x(jnl)

            ! BC is E_x = 0
            dphi_dx = 2.d0 * (phi(inl,jnl+1) - phi(inl,jnl)) / dx_p**2.d0

        ! X-dir right (vacuum) is fixed phi = 0
        else if (side_x > 0) then
            dx_m = x(jnl) - x(jnl-1)

            dphi_dx = 2.d0 * (phi(inl,jnl-1) - phi(inl,jnl)) / dx_m**2.d0
        end if
    end if

    if (nr > 1) then
        ! r-dir Center
        if (side_r .eq. 0) then
            dr_p  = r(inl+1) - r(inl)
            dr_m = r(inl)   - r(inl-1)

            ! term: d^2(phi)/dr^2 (r-term)
            dphi_dr = 2.d0*phi(inl+1,jnl)/(dr_p * (dr_p+dr_m)) &
            - 2.d0*phi(inl,jnl)/(dr_p*dr_m) &
            + 2.d0*phi(inl-1,jnl)/(dr_m * (dr_p+dr_m)) &
            + (phi(inl+1,jnl) - phi(inl-1,jnl)) / ( r(inl) * (dr_p + dr_m) )

        ! r-dir left (vacuum)
        else if (side_r < 0) then
            dr_p = r(inl+1) - r(inl)

            ! BC is E_r = 0
            dphi_dr = 2.d0 * (phi(inl+1,jnl) - phi(inl,jnl)) / dr_p**2.d0

        ! r-dir right (vacuum) is fixed phi = 0
        else if (side_r > 0) then
            dr_m = r(inl) - r(inl-1)

            dphi_dr = 2.d0 * (phi(inl-1,jnl) - phi(inl,jnl)) / dr_m**2.d0
        end if
    end if
    
    !if (rank == 1) then
    !write(*,'(4I5)') inl,jnl,side_x,side_r
    !write(*,'(2F5.1)') dr_p, dr_m
    !write(*,'(5F5.1)') phi(inl,jnl+1), phi(inl,jnl), phi(inl,jnl+1)
    !write(*,'(3F5.1)') dphi_dx, dphi_dr, max(dx_p,dx_m,dr_p,dr_m)
    !read(*,*) i
    !end if
    f_temp = (dphi_dx + dphi_dr) * max(dx_p,dx_m,dr_p,dr_m)

end subroutine laplace

subroutine jacobian(i_local, j_local, i_global, j_global, z, side_z, r, side_r, phi, &
                    nz, nr, nz_loc, neqn, node_global, b_temp, cols, A_temp)
    integer, intent(in):: i_local, j_local, i_global, j_global, side_z, side_r, &
                          nz, nr, nz_loc, neqn, node_global(nr,nz)
    integer, intent(inout):: cols(5)
    real(8), intent(in):: z(nz), r(nr), b_temp
    real(8), intent(inout):: phi(nr,nz_loc), A_temp(1,5)
    logical:: zero_perturb
    real(8):: perturb, b_pert
    real(8):: temporary_Storage
    integer:: I, J, K, width, k_start, k_stop, cols_idx
    integer, dimension(5,2):: stencil

    ! initialize
    cols_idx = 0
    temporary_Storage = 0
    Perturb = 0.0001
    stencil = 0
    width = 0
    
    k_start = 1
    k_stop  = 5
    if (nz .eq. 1) k_start = 3
    if (nr .eq. 1) k_stop  = 3
    
    do K = k_start, k_stop
        if ((K .eq. 1) .and. (i_global .ne. 1)) then
            if (node_global(i_global-1,j_global) > 0) then
                width = width + 1
                stencil(width,1) = -1
                stencil(width,2) =  0
            end if
        else if ((K .eq. 2) .and. (i_global .ne. nr)) then
            if (node_global(i_global+1,j_global) > 0) then
                width = width + 1
                stencil(width,1) = 1
                stencil(width,2) = 0
            end if
        else if (K .eq. 3) then
            width = width + 1
            stencil(width,1) = 0
            stencil(width,2) = 0
        else if ((K .eq. 4) .and. (j_global .ne. nz)) then
            if (node_global(i_global,j_global+1) > 0) then
                width = width + 1
                stencil(width,1) = 0
                stencil(width,2) = 1
            end if
        else if ((K .eq. 5) .and. (j_global .ne. 1)) then
            if (node_global(i_global,j_global-1) > 0) then
                width = width + 1
                stencil(width,1) =  0
                stencil(width,2) = -1
            end if
        end if
    end do
    
    DO k=1, width
        I = i_local + stencil(k,1)
        J = j_local + stencil(k,2)
        
        zero_perturb = .false.
!-----------------------------------------------------------------------
!*******************JACOBIAN SEtUP*************************
!-----------------------------------------------------------------------
        cols_idx = cols_idx + 1
        temporary_Storage = phi(I,J)
        if (Abs(phi(I,J)) > 0) then
            phi(I,J) = phi(I,J) + &
                phi(I,J)*perturb
        else
            zero_perturb = .true.
            phi(I,J) = perturb
        end if
        call laplace(i_local, j_local, z, side_z, r, side_r, phi, &
                     nz, nr, nz_loc, b_pert)
        if (.not. zero_perturb) then
            phi(I,J) = temporary_Storage
        else
            phi(I,J) = 1
        end if
        
        cols(cols_idx) = node_global(i_global + stencil(k,1), &
                                     j_global + stencil(k,2)) - 1
        A_temp(1,cols_idx) = (b_pert - b_temp)/(phi(I,J)*perturb)
        
        if (Zero_Perturb) then
            phi(I,J) = temporary_Storage
        end if
    end DO
end subroutine jacobian
end module laplace_lib


