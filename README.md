# fastGMRES-GPU-CFD-dissertation

Supporting code for the MSc dissertation:

**A Study of Sparse Linear Solvers for Unstructured-Grid CFD: GPU Acceleration and Randomized Flexible GMRES**

**Bo Cui**  
Department of Mathematics  
The University of Manchester  
2026

---

## Overview

This repository contains supporting code for my MSc dissertation on the numerical
solution of large sparse linear systems arising from unstructured-grid computational
fluid dynamics (CFD).

The dissertation investigates two complementary aspects of sparse iterative solvers:

1. **GPU acceleration of sparse linear algebra**, using PETSc and CUDA on linear
   systems extracted from the open-source CFD framework *Fluidity*.

2. **Numerical evaluation of fastGMRES**, a flexible GMRES framework using an
   inner randomized sketched GMRES iteration.

The repository also contains scripts used to profile the sequential fastGMRES
implementation and to construct a model-based estimate of its potential performance
on GPU hardware.

---

## Dissertation Structure

The code in this repository is organised according to the main computational
experiments in the dissertation.

### Chapter 3 — Fluidity Matrix Extraction and GPU Benchmarking

Chapter 3 develops a non-intrusive workflow for extracting sparse linear systems
from *Fluidity* and reconstructing them in a standalone PETSc/CUDA environment.

The extracted Pressure Poisson Equation (PPE) system was solved using equivalent
CG--Jacobi configurations on a single CPU process and an NVIDIA GPU.

Over twenty timed solves, the mean PETSc `KSPSolve` time decreased from

- **CPU:** 5.0936 s
- **GPU:** 0.1475 s

corresponding to a measured solve-phase speedup of approximately **34.53×**.

The corresponding code will be placed in:

```text
chapter3_gpu/
```

This directory contains scripts and source files related to:

- Fluidity CSR matrix extraction;
- binary `.dat` matrix loading;
- Fortran-to-C indexing conversion;
- PETSc sparse matrix construction;
- PETSc CUDA matrix and vector types;
- CG--Jacobi CPU/GPU benchmarking;
- generation of the GPU performance figure used in the dissertation.

---

### Chapter 5 — Numerical Evaluation of fastGMRES

Chapter 5 evaluates fastGMRES on non-symmetric sparse linear systems.

The first test problem is a non-symmetric momentum matrix extracted from the
`FluAtm` branch of *Fluidity*.

The study is then extended to ten public non-symmetric sparse matrices from the
SuiteSparse Matrix Collection.

Three Krylov methods are compared:

- restarted GMRES(50);
- BiCGSTAB;
- fastGMRES.

Each public matrix is tested under two configurations:

1. without an external preconditioner;
2. with ILU(0) right preconditioning.

This produces twenty public benchmark cases in total.

Across these twenty cases, the convergence counts reported in the dissertation are:

| Solver | No preconditioner | ILU(0) | Overall |
|---|---:|---:|---:|
| GMRES(50) | 6/10 | 8/10 | 14/20 |
| BiCGSTAB | 5/10 | 7/10 | 12/20 |
| fastGMRES | 8/10 | 9/10 | 17/20 |

The corresponding MATLAB code will be placed in:

```text
chapter5_fastgmres/
```

This directory contains code for:

- the fastGMRES implementation;
- the detailed FluAtm momentum-matrix experiment;
- GMRES(50), BiCGSTAB and fastGMRES comparisons;
- the ten-matrix SuiteSparse benchmark;
- ILU(0) right-preconditioned experiments;
- runtime performance profiles;
- residual and convergence analysis.

---

### Chapter 6 — GPU Performance Estimate for fastGMRES

Chapter 6 investigates the possible performance of a future GPU implementation of
fastGMRES.

The sequential MATLAB implementation was profiled on the FluAtm momentum system.

The measured baseline was approximately:

```text
CPU fastGMRES runtime : 54.02 s
Total counted SpMVs   : 9203
True relative residual: approximately 9.84e-8
```

A separate CPU microbenchmark estimated an average sparse matrix-vector product
cost of approximately:

```text
3.32e-3 s
```

The estimated cumulative SpMV time was therefore approximately 30.6 s, corresponding
to about 56.6% of the measured sequential fastGMRES runtime.

The measured GPU result from Chapter 3 was then used as an empirical acceleration
reference in an Amdahl-law model:

```text
S = 1 / ((1 - p) + p / 34.53)
```

where `p` is the assumed fraction of the sequential fastGMRES workload that can be
effectively accelerated on a GPU.

The dissertation considers three illustrative scenarios:

| Scenario | Assumed effectively GPU-accelerated fraction | Estimated speedup |
|---|---:|---:|
| Conservative | 56.6% | 2.22× |
| Moderate | 70% | 3.12× |
| Optimistic | 80% | 4.48× |

These values are **model-based projections rather than measured GPU fastGMRES
results**.

The corresponding profiling code will be placed in:

```text
chapter6_profiling/
```

---

## Repository Structure

The repository is organised as follows:

```text
fastGMRES-GPU-CFD-dissertation/
│
├── README.md
├── .gitignore
│
├── chapter3_gpu/
│   ├── README.md
│   ├── matrix_loading/
│   ├── petsc_cuda/
│   └── plotting/
│
├── chapter5_fastgmres/
│   ├── README.md
│   ├── fastgmres.m
│   ├── fluidity_case/
│   ├── public_matrix_benchmark/
│   └── plotting/
│
├── chapter6_profiling/
│   ├── README.md
│   └── profiling/
│
└── results/
    └── selected_summary_files/
```

The exact directory structure may be refined as the supporting code is cleaned and
uploaded.

---

## Software Environment

The experiments in the dissertation used two different computational environments.

### PETSc/CUDA experiments

The GPU experiments used:

```text
Operating system : Ubuntu 22.04.5 LTS
CPU              : AMD EPYC 7713
GPU              : NVIDIA A100-PCIE-40GB
PETSc            : 3.23.7
CUDA Toolkit     : 13.1
Compiler         : GCC / gfortran
```

PETSc CUDA matrix and vector types used in the experiments include:

```text
MATSEQAIJCUSPARSE
VECSEQCUDA
```

The sparse GPU backend uses NVIDIA cuSPARSE.

### MATLAB experiments

The sequential fastGMRES experiments used:

```text
Operating system : Windows 11
CPU              : Intel Core i7-9750H
Memory           : 16 GB
MATLAB           : R2026a
```

No explicit GPU acceleration or MATLAB Parallel Computing Toolbox was used for the
Chapter 5 fastGMRES benchmark.

---

## Public Benchmark Matrices

The ten public matrices used in Chapter 5 are:

```text
NACA12
ML_Laplace
vas_stokes_1M
Transport
cavity10
ns3Da
atmosmodl
atmosmodd
GT01R
poisson3Da
```

These matrices were obtained from the **SuiteSparse Matrix Collection**.

The original SuiteSparse matrix files are not redistributed in this repository.
Users wishing to reproduce the benchmark should obtain the corresponding matrices
from the SuiteSparse Matrix Collection.

For experimental consistency, the public benchmark does not use dataset-specific
right-hand sides. Instead, each test constructs a reproducible consistent system:

```matlab
x_true = randn(n,1);
b = A*x_true;
```

using a fixed random seed for each matrix.

---

## Convergence Criterion

For the MATLAB benchmarks, solver convergence is determined using the independently
evaluated true relative residual:

```text
||b - A*x||_2 / ||b||_2 <= 1e-7
```

A solve is classified as unsuccessful if this criterion is not satisfied before the
specified iteration or runtime limit is reached.

---

## fastGMRES Configuration

The detailed FluAtm case and the public benchmark use different solver limits.

For the public benchmark, the principal fastGMRES settings are:

```text
Tolerance                  : 1e-7
Maximum outer iterations   : 150
Maximum inner dimension m  : 80
Sketch dimension s         : 162
Arnoldi truncation t       : 0
Conditioning threshold     : 1e15
Maximum runtime            : 600 s
Random seed                : rng('default')
```

The same parameter configuration is used across the ten public matrices.

---

## Preconditioning

The extended public benchmark additionally uses zero-fill incomplete LU
factorisation:

```text
A ≈ L*U
```

as a common ILU(0) right preconditioner.

The same ILU(0) preconditioned operator is used for GMRES(50), BiCGSTAB and
fastGMRES.

ILU setup time is recorded separately and is not included in the iterative solver
runtime used in the runtime performance profiles.

---

## Data Availability

Large matrix and simulation data files are not necessarily included directly in this
repository.

In particular:

- SuiteSparse matrices should be downloaded from the original SuiteSparse Matrix
  Collection;
- some *Fluidity*-generated binary matrix files may be omitted because of file size,
  redistribution, or project-data considerations;
- scripts describe the expected matrix filenames and formats where required.

The extracted *Fluidity* matrices use binary files containing:

```text
CSR row pointers
CSR column indices
non-zero matrix values
right-hand side vector
```

---

## Reproducing the Experiments

Detailed instructions will be provided in the README file associated with each
chapter directory.

The intended workflow is approximately:

### Chapter 5 public benchmark

```text
1. Download the required SuiteSparse matrices.
2. Place the matrix files in the benchmark directory.
3. Ensure fastgmres.m is available on the MATLAB path.
4. Run the public benchmark script.
5. Inspect the generated convergence and runtime summaries.
6. Generate the performance profiles.
```

### Chapter 6 profiling experiment

```text
1. Place the FluAtm matrix data in the profiling directory.
2. Ensure fastgmres.m is available.
3. Run the fastGMRES profiling script.
4. Record the baseline runtime and SpMV count.
5. Run the kernel microbenchmarks.
6. Use the reported values in the Amdahl-law performance model.
```

---

## Important Notes

The GPU and fastGMRES experiments were performed on different matrices and in
different software environments.

Therefore, the GPU speedup estimated for fastGMRES in Chapter 6 should not be
interpreted as a directly measured GPU result.

The estimated range of approximately **2.2--4.5×** is a model-based projection based
on:

- measured CPU fastGMRES profiling;
- microbenchmark estimates of SpMV cost;
- the measured PETSc/CUDA acceleration obtained in Chapter 3;
- Amdahl's law.

A native C++/CUDA implementation of fastGMRES would be required to verify this
estimate experimentally.

---

## Related Software and References

The dissertation makes use of or discusses the following software and numerical
methods:

- *Fluidity*
- PETSc
- NVIDIA CUDA
- cuSPARSE
- GMRES
- Flexible GMRES
- BiCGSTAB
- randomized/sketched GMRES
- fastGMRES
- ILU(0)
- SuiteSparse Matrix Collection

The mathematical formulation, references, experimental methodology and complete
discussion are provided in the dissertation.

---

## Code Availability and Attribution

This repository contains supporting research code associated with my MSc
dissertation.

Some algorithms implemented or adapted in this repository are based on methods
described in the academic literature. Their inclusion here does not imply authorship
of the underlying algorithms.

In particular, fastGMRES refers to the Flexible GMRES / randomized sGMRES framework
studied in:

> Stefan Güttel and John W. Pearson,  
> *Stabilizing randomized GMRES through flexible GMRES*.

Relevant references and citations are provided in the dissertation and should be
consulted when using or discussing these algorithms.

---

## Dissertation

**Bo Cui**

*A Study of Sparse Linear Solvers for Unstructured-Grid CFD: GPU Acceleration and
Randomized Flexible GMRES*

MSc Dissertation  
Department of Mathematics  
The University of Manchester  
2026

---

## Author

**Bo Cui**

Department of Mathematics  
The University of Manchester

GitHub repository:

```text
https://github.com/cui516615-cpu/fastGMRES-GPU-CFD-dissertation
```

---

## License

No open-source licence has currently been assigned to this repository.

Please contact the author before redistributing or reusing substantial portions of
the repository.

Third-party software and algorithms remain subject to their respective licences,
copyrights and citation requirements.
