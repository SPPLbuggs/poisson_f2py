    module laplace_lib
    use properties
    implicit none

    contains

    subroutine laplace(inl, jnl, side_z, side_r, b_temp)
    integer, intent(in):: inl, jnl, side_z, side_r
    real(dp), intent(out):: b_temp
    real(dp):: dz_p = 0, dz_m = 0
    real(dp):: dr_p = 0, dr_m = 0
    real(dp):: dphi_dz = 0, dphi_dr = 0

!-----------------------------------------------------------------------
!*******************POISSON equation***********************
!-----------------------------------------------------------------------
    if (nz > 1) then
        ! X-dir Center
        if (side_z .eq. 0) then
            dz_p  = z(jnl+1) - z(jnl)
            dz_m = z(jnl)   - z(jnl-1)

            ! term: d^2(phi)/dz^2 (x-term)
            dphi_dz = 2.0_dp*phi(inl,jnl+1)/(dz_p &
            * (dz_p+dz_m)) &
            - 2.0_dp*phi(inl,jnl)/(dz_p*dz_m) &
            + 2.0_dp*phi(inl,jnl-1)/(dz_m &
            * (dz_p+dz_m))

        ! X-dir left (vacuum)
        else if (side_z < 0) then
            dz_p = z(jnl+1) - z(jnl)

            ! BC is E_x = 0
            dphi_dz = 2.0_dp * (phi(inl,jnl+1) - phi(inl,jnl)) / dz_p**2.0_dp

        ! X-dir right (vacuum) is fixed phi = 0
        else if (side_z > 0) then
            dz_m = z(jnl) - z(jnl-1)

            dphi_dz = 2.0_dp * (phi(inl,jnl-1) - phi(inl,jnl)) / dz_m**2.0_dp
        end if
    end if

    if (nr > 1) then
        ! r-dir Center
        if (side_r .eq. 0) then
            dr_p  = r(inl+1) - r(inl)
            dr_m = r(inl)   - r(inl-1)

            ! term: d^2(phi)/dr^2 (r-term)
            dphi_dr = 2.0_dp*phi(inl+1,jnl)/(dr_p * (dr_p+dr_m)) &
            - 2.0_dp*phi(inl,jnl)/(dr_p*dr_m) &
            + 2.0_dp*phi(inl-1,jnl)/(dr_m * (dr_p+dr_m)) &
            + (phi(inl+1,jnl) - phi(inl-1,jnl)) / ( r(inl) * (dr_p + dr_m) )

        ! r-dir left (vacuum)
        else if (side_r < 0) then
            dr_p = r(inl+1) - r(inl)

            ! BC is E_r = 0
            dphi_dr = 2.0_dp * (phi(inl+1,jnl) - phi(inl,jnl)) / dr_p**2.0_dp

        ! r-dir right (vacuum) is fixed phi = 0
        else if (side_r > 0) then
            dr_m = r(inl) - r(inl-1)

            dphi_dr = 2.0_dp * (phi(inl-1,jnl) - phi(inl,jnl)) / dr_m**2.0_dp
        end if
    end if
    
    !if (rank == 1) then
    !write(*,'(4I5)') inl,jnl,side_z,side_r
    !write(*,'(2F5.1)') dr_p, dr_m
    !write(*,'(5F5.1)') phi(inl,jnl+1), phi(inl,jnl), phi(inl,jnl+1)
    !write(*,'(3F5.1)') dphi_dz, dphi_dr, max(dz_p,dz_m,dr_p,dr_m)
    !read(*,*) i
    !end if
    b_temp = (dphi_dz + dphi_dr) * max(dz_p,dz_m,dr_p,dr_m)

    end subroutine laplace

    subroutine jacobian(i_local, j_local, i_global, j_global, side_z, side_r, &
                    b_temp, cols, A_temp)
    integer, intent(in):: i_local, j_local, i_global, j_global, side_z, side_r
    integer, intent(inout):: cols(5)
    real(dp), intent(in):: b_temp
    real(dp), intent(inout):: A_temp(1,5)
    logical:: zero_perturb
    real(dp):: perturb, b_pert
    real(dp):: temporary_Storage
    integer:: I, J, K, width, k_start, k_stop, cols_idz
    integer, dimension(5,2):: stencil

    ! initialize
    cols_idz = 0
    temporary_Storage = 0
    Perturb = 1e-4_dp
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
        cols_idz = cols_idz + 1
        temporary_Storage = phi(I,J)
        if (Abs(phi(I,J)) > 0) then
            phi(I,J) = phi(I,J) + &
                phi(I,J)*perturb
        else
            zero_perturb = .true.
            phi(I,J) = perturb
        end if
        call laplace(i_local, j_local, side_z, side_r, b_pert)
        if (.not. zero_perturb) then
            phi(I,J) = temporary_Storage
        else
            phi(I,J) = 1
        end if
        
        cols(cols_idz) = node_global(i_global + stencil(k,1), &
                                     j_global + stencil(k,2)) - 1
        A_temp(1,cols_idz) = (b_pert - b_temp)/(phi(I,J)*perturb)
        
        if (Zero_Perturb) then
            phi(I,J) = temporary_Storage
        end if
    end DO
    end subroutine jacobian
    end module laplace_lib

