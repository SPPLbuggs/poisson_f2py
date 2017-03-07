    module poisson_solver
    use laplace_lib
    implicit none

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscksp.h>
#include <petsc/finclude/petscpc.h>

    public main

    MPI_Comm comm
    PetscMPIInt rank,tasks
    Mat A
    Vec b,u,x
    KSP ksp
    PC pc
    PetscInt Istart, Iend, its
    PetscReal  norm
    PetscErrorCode ierr
    PetscScalar none
    
    contains

    subroutine main(neqn, z, type_z, r, type_r, phi, local_idx, global_idx, &
    node_global, nz, nr,nz_loc)

    integer, intent(in):: neqn, nz, nr, nz_loc, type_z(nr,nz), type_r(nr,nz), &
    global_idx(nz*nr,2),local_idx(nr*nz,2), node_global(nr,nz)
    real(8), intent(in):: z(nz_loc), r(nr)
    real(8), intent(inout):: phi(nr,nz_loc)
    integer:: i_global, j_global, i_local, j_local, node, cols(5), &
              loc_cols, loc_rows
    real(8):: A_temp(1,5), b_temp, soln
    
    call PetscInitialize(petsc_null_character, ierr)
    comm = PETSC_COMM_WORLD
    
    none = -1.0
    
    call MPI_Comm_rank(comm,rank,ierr)
    call MPI_Comm_size(comm,tasks,ierr)
    
    ! create parallel matrix
    ! parallel partitioning determined at runtime
    call MatCreate(comm,A,ierr)
    
    loc_cols = neqn
    loc_rows = 0
    do j_global = max(1,nz/tasks*rank+1), min(nz,nz/tasks*(rank+1))
        do i_global = 1,nr
            if (node_global(i_global, j_global) > 0) loc_rows = loc_rows+1
        end do
    end do
            
    call MatSetSizes(A,loc_rows,loc_rows,neqn,neqn,ierr)
    call MatSetUp(A,ierr)
    
    call MatSetFromOptions(A,ierr)
    call MatSeqAIJSetPreallocation(A, 5, petsc_null_integer, ierr)
    call MatSetOption(A, mat_ignore_zero_entries, petsc_true, ierr)
    
    ! find parallel partitioning range
    call MatGetOwnershipRange(A,Istart,Iend,ierr)
    
    ! create parallel vectors
    ! parallel partitioning determined at runtime
    call VecCreateMPI(comm,loc_rows,neqn,u,ierr)
    call VecCreateMPI(comm,loc_rows,neqn,b,ierr)
    call VecCreateMPI(comm,loc_rows,neqn,x,ierr)
    call VecSetFromOptions(u,ierr)
    call VecSetFromOptions(b,ierr)
    call VecSetFromOptions(x,ierr)
    call VecSetOption(b, vec_ignore_negative_indices, petsc_true, ierr)   
    
    ! assemble A and b
    do node = Istart, Iend-1
        i_local  = local_idx(node+1, 1)
        j_local  = local_idx(node+1, 2)
        i_global = global_idx(node+1,1)
        j_global = global_idx(node+1,2)
        
        b_temp = 0.d0
        A_temp = 0.d0
        cols = -5

        call laplace(i_local, j_local, z, type_z(i_global,j_global), r, &
                     type_r(i_global,j_global), phi, nz, nr, nz_loc, b_temp)
        call jacobian(i_local, j_local, i_global, j_global, z, &
                      type_z(i_global,j_global), r, type_r(i_global,j_global), &
                      phi, nz, nr, nz_loc, neqn, node_global, b_temp, &
                      cols, A_temp)
        
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
    
    call KSPDestroy(ksp,ierr)
    call VecDestroy(u,ierr)
    call VecDestroy(x,ierr)
    call VecDestroy(b,ierr)
    call MatDestroy(A,ierr)
    
    call PetscFinalize(ierr)
    end subroutine
    
    end module
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
