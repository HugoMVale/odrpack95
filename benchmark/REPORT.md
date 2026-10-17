# BLAS run-time comparison: reference BLAS vs OpenBLAS

## Setup

- Compiler: gfortran (MinGW, Windows), release build, `fortran_std=none`.
- Builds: `builddir-refblas` (`-Dblas_option=refblas`) and `builddir-openblas` (`-Dblas_option=openblas`), both with `-Dbuild_benchmark=true`.
- Benchmark: `benchmark/bench_odr.f90`. Explicit linear and exponential models, n = 1000, 10000, 100000, with free and bounded parameters. Seeded synthetic data, so every build fits identical problems. Only the `odr` call is timed.
- Runner: `benchmark/run_benchmark.ps1 -Reps 10`. One warm-up run per build, then 10 recorded repetitions per case. Run order alternates between repetitions.
- Configurations: reference BLAS; OpenBLAS with default threads (`openblas_mt`); OpenBLAS pinned to one thread (`openblas_1t`).
- Raw data: `benchmark/results/raw.csv`. Medians: `benchmark/results/summary.csv`.

## Results (median wall time, seconds)

| Model     |      n | Bound   |      RefBLAS |  OpenBLAS MT |  OpenBLAS 1T | Saving MT (%) | Saving 1T (%) |
| :-------- | -----: | :------ | -----------: | -----------: | -----------: | ------------: | ------------: |
| linear    |   1000 | free    |     0.001140 |     0.001113 |     0.001167 |           2.4 |          -2.3 |
| linear    |   1000 | bounded |     0.001056 |     0.000999 |     0.001008 |           5.4 |           4.5 |
| linear    |  10000 | free    |     0.012450 |     0.011870 |     0.011715 |           4.7 |           5.9 |
| linear    |  10000 | bounded |     0.011898 |     0.011050 |     0.012091 |           7.1 |          -1.6 |
| linear    | 100000 | free    |     0.139472 |     0.139803 |     0.134874 |          -0.2 |           3.3 |
| linear    | 100000 | bounded |     0.137150 |     0.138056 |     0.132430 |          -0.7 |           3.4 |
| expo      |   1000 | free    |     0.006810 |     0.006491 |     0.006357 |           4.7 |           6.6 |
| expo      |   1000 | bounded |     0.006744 |     0.006420 |     0.006263 |           4.8 |           7.1 |
| expo      |  10000 | free    |     0.111309 |     0.106928 |     0.108243 |           3.9 |           2.8 |
| expo      |  10000 | bounded |     0.110339 |     0.107597 |     0.107984 |           2.5 |           2.1 |
| expo      | 100000 | free    |     2.166439 |     2.202090 |     2.120389 |          -1.6 |           2.1 |
| expo      | 100000 | bounded |     2.158921 |     2.214665 |     2.117968 |          -2.6 |           1.9 |
| **Total** |        |         | **4.863728** | **4.947082** | **4.760489** |      **-1.7** |       **2.1** |

Every case converged with `info = 1` in all builds. The largest relative difference in fitted parameters between reference BLAS and OpenBLAS is 3.1e-6, for the exponential model at n = 100000.

## Decision

The total saving from OpenBLAS is -1.7% (default threads) and 2.1% (one thread). Both are below the 10% threshold, so the self-contained reference BLAS is the recommended default.

## Caveats

- Per-case differences are of the same order as run-to-run noise. The summary reports medians only, and the spread is not reported.
- The bundled BLAS interfaces cover only level-1 routines. OpenBLAS gains would likely be larger for level-2 or level-3 routines, which this library does not use.
- Results come from one Windows machine with synthetic data and two model shapes. Other platforms or problem shapes could differ.
- The 1e-6 parameter difference in the exponential model comes from changed floating-point operation order in BLAS, not a correctness failure. Convergence status matches in every case.
