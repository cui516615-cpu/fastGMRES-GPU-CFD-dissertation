#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <petscmat.h>
#include <petscksp.h>
#include <petsctime.h>

int main(int argc, char **args) {
    PetscCall(PetscInitialize(&argc, &args, NULL, NULL));
    printf("[GPU Sandbox] Robust Benchmarking Started...\n");

    int num_rows_plus_one, nnz, x_length;
    FILE *f;

    f = fopen("fluidity_matrix_row.dat", "rb");
    fread(&num_rows_plus_one, sizeof(int), 1, f);
    int *row_ptr = (int*)malloc(num_rows_plus_one * sizeof(int));
    fread(row_ptr, sizeof(int), num_rows_plus_one, f);
    fclose(f);

    f = fopen("fluidity_matrix_col.dat", "rb");
    fread(&nnz, sizeof(int), 1, f);
    int *col_idx = (int*)malloc(nnz * sizeof(int));
    fread(col_idx, sizeof(int), nnz, f);
    fclose(f);

    f = fopen("fluidity_matrix_val.dat", "rb");
    fread(&nnz, sizeof(int), 1, f);
    double *values = (double*)malloc(nnz * sizeof(double));
    fread(values, sizeof(double), nnz, f);
    fclose(f);

    f = fopen("fluidity_rhs.dat", "rb");
    fread(&x_length, sizeof(int), 1, f);
    double *rhs_val = (double*)malloc(x_length * sizeof(double));
    fread(rhs_val, sizeof(double), x_length, f);
    fclose(f);

    PetscInt nrows = num_rows_plus_one - 1;
    if (row_ptr[0] == 1) {
        for(int i = 0; i <= nrows; i++) row_ptr[i] -= 1;
        for(int i = 0; i < nnz; i++) col_idx[i] -= 1;
    }

    Mat A; Vec b, x; KSP ksp; PC pc;

    MatCreateSeqAIJWithArrays(PETSC_COMM_SELF, nrows, nrows, row_ptr, col_idx, values, &A);
    MatSetType(A, MATSEQAIJCUSPARSE); 
    MatSetFromOptions(A);

    VecCreateSeqWithArray(PETSC_COMM_SELF, 1, nrows, rhs_val, &b);
    VecSetType(b, VECSEQCUDA);

    VecCreate(PETSC_COMM_SELF, &x);
    VecSetSizes(x, nrows, nrows);
    VecSetType(x, VECSEQCUDA);

    KSPCreate(PETSC_COMM_SELF, &ksp);
    KSPSetOperators(ksp, A, A);
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCJACOBI);
    KSPSetType(ksp, KSPCG);
    KSPSetTolerances(ksp, 1.0e-7, PETSC_DEFAULT, PETSC_DEFAULT, 10000);
    KSPSetFromOptions(ksp); 


    int n_warmup = 5;
    int n_runs = 20;
    double times[20];
    double sum = 0.0, mean = 0.0, variance = 0.0, std_dev = 0.0;
    double min_time = 1e9, max_time = 0.0;

    printf("[GPU Sandbox] Performing %d Warm-up runs to initialize CUDA context...\n", n_warmup);
    for(int i = 0; i < n_warmup; i++) {
        VecSet(x, 0.0); // Reset initial guess!
        KSPSolve(ksp, b, x);
    }

    printf("[GPU Sandbox] Performing %d Benchmark runs...\n", n_runs);
    for(int i = 0; i < n_runs; i++) {
        VecSet(x, 0.0); // Reset initial guess to ensure identical computational work
        
        PetscLogDouble t1, t2;
        PetscTime(&t1);
        KSPSolve(ksp, b, x);  // STRICTLY timing KSPSolve only
        PetscTime(&t2);
        
        times[i] = t2 - t1;
        sum += times[i];
        if (times[i] < min_time) min_time = times[i];
        if (times[i] > max_time) max_time = times[i];
    }

    /* Calculate Statistics */
    mean = sum / n_runs;
    for(int i = 0; i < n_runs; i++) {
        variance += (times[i] - mean) * (times[i] - mean);
    }
    std_dev = sqrt(variance / n_runs);

    PetscInt its;
    KSPGetIterationNumber(ksp, &its);

    printf("\n======================================================\n");
    printf("📊 [Robust GPU Performance Benchmark]\n");
    printf("   - Iterations        : %d\n", its);
    printf("   - Mean Solve Time   : %.4f seconds\n", mean);
    printf("   - Min (Best) Time   : %.4f seconds\n", min_time);
    printf("   - Max Time          : %.4f seconds\n", max_time);
    printf("   - Std Deviation     : %.4f seconds\n", std_dev);
    printf("======================================================\n\n");

    KSPDestroy(&ksp); MatDestroy(&A); VecDestroy(&b); VecDestroy(&x);
    free(row_ptr); free(col_idx); free(values); free(rhs_val);

    PetscCall(PetscFinalize());
    return 0;
}
