    module mod
    use properties
    use laplace_lib
    use petsc_lib
    implicit none
    
    !f2py integer:: neqn, nz, nr, nz_loc
    !f2py integer, allocatable :: type_z(:,:), type_r(:,:), global_idx(:,:)
    !f2py integer, allocatable :: local_idx(:,:), node_global(:,:)
    !f2py real(8), allocatable :: z(:), r(:), phi(:,:)
    
    contains
    subroutine run
    integer:: i_global, j_global, i_local, j_local, node, cols(5)
    real(dp):: A_temp(1,5), b_temp, soln
    
    call petsc_initialize
    
    ! assemble A and b
    do node = Istart, Iend-1
        i_local  = local_idx(node+1, 1)
        j_local  = local_idx(node+1, 2)
        i_global = global_idx(node+1,1)
        j_global = global_idx(node+1,2)
        
        b_temp = 0.0_dp
        A_temp = 0.0_dp
        cols = -5
        
        call laplace(i_local, j_local, type_z(i_global,j_global), &
                     type_r(i_global,j_global), b_temp)
        call jacobian(i_local, j_local, i_global, j_global, &
                      type_z(i_global,j_global), type_r(i_global,j_global), &
                      b_temp, cols, A_temp)
        
        call MatSetValues(A, 1, node, 5, cols, A_temp, insert_values, ierr)
        call VecSetValues(b, 1, node, -b_temp, insert_values, ierr)
    end do
    
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
    
    call vecassemblybegin(b, ierr)
    call vecassemblyend(b, ierr)
    
    !call VecView(b,PETSC_VIEWER_STDOUT_WORLD,ierr)
    !call MatView(A,PETSC_VIEWER_STDOUT_WORLD,ierr)
    !read(*,*) j_global

    ! create linear solver
    call KSPCreate(comm,ksp,ierr)
    call KSPSetOperators(ksp,A,A,ierr)
    call KSPSetFromOptions(ksp,ierr)
    !call KSPGetPC(ksp,pc,ierr)
    !call KSPsetType(ksp,KSPCG,ierr)
    !call PCSetType(pc,PCMG,ierr)
    call KSPSetTolerances(ksp,1.d-7,1.d-7,PETSC_DEFAULT_REAL,50000,ierr)
    
    ! solve system
    call KSPSolve(ksp,b,x,ierr)
    
    do node = Istart, Iend-1
        call VecGetValues(x, 1, node, soln, ierr)
        i_local = local_idx(node + 1, 1)
        j_local = local_idx(node + 1, 2)
        phi(i_local, j_local) = phi(i_local, j_local) + soln
    end do
  
    call MatMult(A,x,u,ierr)
    call VecAYPX(u,none,b,ierr)
    call VecNorm(u,NORM_2,norm,ierr)
    call KSPGetIterationNumber(ksp,its,ierr)
    if (rank .eq. 0) then
        if (norm .gt. 1.e-12) then
            write(6,100) norm,its
        else
            write(6,110) its
        endif
    endif
100 format('Norm of error ',e11.4,' iterations ',i5)
110 format('Norm of error < 1.e-12 iterations ',i5)
    
    call petsc_finalize
    end subroutine
    
    end module
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
