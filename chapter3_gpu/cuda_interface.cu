#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#include <petscmat.h>
#include <petscksp.h>
#include <petsctime.h>

extern "C" void petsc_solve_cuda_(Mat *A, Vec *b, Vec *y) {
    // Correct C syntax trace print with argument pointers casted to void*
    printf("\n[GPU Bridge Trace] Entered petsc_solve_cuda_! A=%p, b=%p, y=%p\n", (void*)A, (void*)b, (void*)y);
    fflush(stdout);

    PetscLogDouble t_start, t_gpu_copy_start, t_gpu_copy_end;
    PetscLogDouble t_setup_start, t_setup_end;
    PetscLogDouble t_solve_start, t_solve_end;
    PetscLogDouble t_gpu_back_start, t_gpu_back_end, t_end;

    Mat A_gpu = NULL;
    Vec b_gpu = NULL, y_gpu = NULL;
    KSP ksp = NULL;
    PC pc = NULL;
    MPI_Comm comm;

    // Get the MPI communicator from the matrix
    PetscCallAbort(PETSC_COMM_WORLD, PetscTime(&t_start));
    PetscCallAbort(PETSC_COMM_WORLD, PetscObjectGetComm((PetscObject)*A, &comm));

    PetscCallAbort(comm, PetscTime(&t_gpu_copy_start));

    // 1. Convert CPU Matrix to GPU Matrix (handles both Seq and MPI automatically)
    PetscCallAbort(comm, MatConvert(*A, MATAIJCUSPARSE, MAT_INITIAL_MATRIX, &A_gpu));

    // 2. Duplicate vectors layout and set types to VECCUDA
    PetscCallAbort(comm, VecDuplicate(*b, &b_gpu));
    PetscCallAbort(comm, VecSetType(b_gpu, VECCUDA));
    PetscCallAbort(comm, VecDuplicate(*y, &y_gpu));
    PetscCallAbort(comm, VecSetType(y_gpu, VECCUDA));

    // 3. Copy host data to GPU
    PetscCallAbort(comm, VecCopy(*b, b_gpu));
    PetscCallAbort(comm, VecCopy(*y, y_gpu)); // Includes warm-start initial guess

    cudaDeviceSynchronize();
    PetscCallAbort(comm, PetscTime(&t_gpu_copy_end));

    // 4. Configure KSP Solver on GPU
    PetscCallAbort(comm, PetscTime(&t_setup_start));
    PetscCallAbort(comm, KSPCreate(comm, &ksp));
    PetscCallAbort(comm, KSPSetOperators(ksp, A_gpu, A_gpu));

    PetscCallAbort(comm, KSPSetType(ksp, KSPGMRES));
    PetscCallAbort(comm, KSPGMRESSetRestart(ksp, 30));

    PetscCallAbort(comm, KSPGetPC(ksp, &pc));
    PetscCallAbort(comm, PCSetType(pc, PCJACOBI));

    PetscCallAbort(comm, KSPSetTolerances(ksp, 1e-7, PETSC_DEFAULT, PETSC_DEFAULT, 1000));
    PetscCallAbort(comm, KSPSetFromOptions(ksp));
    PetscCallAbort(comm, KSPSetUp(ksp));

    cudaDeviceSynchronize();
    PetscCallAbort(comm, PetscTime(&t_setup_end));

    // 5. Solve on GPU
    PetscCallAbort(comm, PetscTime(&t_solve_start));
    PetscCallAbort(comm, KSPSolve(ksp, b_gpu, y_gpu));
    cudaDeviceSynchronize();
    PetscCallAbort(comm, PetscTime(&t_solve_end));

    // 6. Copy GPU result vector back to CPU vector *y
    PetscCallAbort(comm, PetscTime(&t_gpu_back_start));
    PetscCallAbort(comm, VecCopy(y_gpu, *y));
    cudaDeviceSynchronize();
    PetscCallAbort(comm, PetscTime(&t_gpu_back_end));

    // 7. Performance metrics
    PetscInt its;
    PetscReal rnorm;
    PetscCallAbort(comm, KSPGetIterationNumber(ksp, &its));
    PetscCallAbort(comm, KSPGetResidualNorm(ksp, &rnorm));

    // 8. Clean up GPU structures (original Mat *A and Vec *y are untouched)
    PetscCallAbort(comm, KSPDestroy(&ksp));
    PetscCallAbort(comm, VecDestroy(&b_gpu));
    PetscCallAbort(comm, VecDestroy(&y_gpu));
    PetscCallAbort(comm, MatDestroy(&A_gpu));

    PetscCallAbort(comm, PetscTime(&t_end));

    PetscInt rank;
    MPI_Comm_rank(comm, &rank);
    if (rank == 0) {
        printf("\n==================== [GPU SOLVER REPORT] ====================\n");
        printf(" MPI Comm Size         : (Dynamic Parallel Active)\n");
        printf(" Iterations Taken      : %d\n", (int)its);
        printf(" Final Residual Norm   : %e\n", (double)rnorm);
        printf("-------------------------------------------------------------\n");
        printf(" CPU->GPU Copy Time    : %f sec\n", t_gpu_copy_end - t_gpu_copy_start);
        printf(" Solver Setup Time     : %f sec\n", t_setup_end - t_setup_start);
        printf(" KSP GPU Solve Time    : %f sec\n", t_solve_end - t_solve_start);
        printf(" GPU->CPU Copy Time    : %f sec\n", t_gpu_back_end - t_gpu_back_start);
        printf("-------------------------------------------------------------\n");
        printf(" Total Execution Time  : %f sec\n", t_end - t_start);
        printf("=============================================================\n\n");
        fflush(stdout);
    }
}
