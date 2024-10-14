# magma 2.7.0

#### Incompatible Dependancies

 - magma requires CUDA less than vesion 12.0 to support texture mapping

 - CUDA 11.7.0 can't be used with GCC > 11.0

 The nvcc flag '-allow-unsupported-compiler' can be used to override this version check

#### magma 2.7.1

Magma 2.7.1 add support for CUDA 12!  2023-02-22

#### CMake Output

```
-- Found CUDA
--     CUDA_CUDART_LIBRARY: CUDA::cudart
--     compile for CUDA arch 6.1 (Pascal)
--     compile for CUDA arch 7.5 (Turing)
```

```
-- Found BLAS: /app/software/FlexiBLAS/3.2.1-GCC-12.2.0/lib/libflexiblas.so
-- Looking for Fortran cheev
-- Looking for Fortran cheev - found
-- Found LAPACK: /app/software/FlexiBLAS/3.2.1-GCC-12.2.0/lib/libflexiblas.so;-lm;-ldl
```


