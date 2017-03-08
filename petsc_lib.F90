    module petsc_lib
    use properties
    implicit none

    contains 
    subroutine petsc_initialize
    
    integer:: i_global, j_global, loc_cols, loc_rows

    call PetscInitialize(petsc_null_character, ierr)
    comm = PETSC_COMM_WORLD
        
    none = -1.0
        
    call MPI_Comm_rank(comm,rank,ierr)
    call MPI_Comm_size(comm,tasks,ierr)
        
    ! create parallel matrix
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
    call VecCreateMPI(comm,loc_rows,neqn,u,ierr)
    call VecCreateMPI(comm,loc_rows,neqn,b,ierr)
    call VecCreateMPI(comm,loc_rows,neqn,x,ierr)
    call VecSetFromOptions(u,ierr)
    call VecSetFromOptions(b,ierr)
    call VecSetFromOptions(x,ierr)
    call VecSetOption(b, vec_ignore_negative_indices, petsc_true, ierr)
    
    end subroutine

    subroutine petsc_finalize
    
    call KSPDestroy(ksp,ierr)
    call VecDestroy(u,ierr)
    call VecDestroy(x,ierr)
    call VecDestroy(b,ierr)
    call MatDestroy(A,ierr)
    call PetscFinalize(ierr)

    end subroutine
    end module
