    module properties
    implicit none
    !integer, parameter:: dp=selected_real_kind(15)
    
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscksp.h>
#include <petsc/finclude/petscpc.h>

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
    
    integer:: neqn, nz, nr, nz_loc
    integer, allocatable :: type_z(:,:), type_r(:,:), global_idx(:,:), &
                            local_idx(:,:), node_global(:,:)
    real(8), allocatable :: z(:), r(:), phi(:,:)

    end module properties
