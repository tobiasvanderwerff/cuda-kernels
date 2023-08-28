Contents: 

```
---
 |-- `histogram`: Parallel histogram implementation.
 |
 |-- `infoli`: A parallelized version of C code that implements a model of a biological
 |   neural network (the inferior-olivary nucleus). The C code was written by someone
 |   else, and the challenge was to make it run faster by parallelizing it using CUDA. 
 |
 |-- `parallel-prefix-sum`: Optimized version of the parallel prefix sum algorithm,
 |    which uses balanced trees using a hierarchical approach, and has O(N) runtime. 
 |    While some other implementations such as Kogge-Stone run in O(N*logN), this version
 |    is work-efficient (which, should be noted, doesn't mean it's necessarily fastest). 
 |
 |-- `sgemm-tiled`: Optimized SGEMM kernel, using tiling. Makes use of common techniques
 |    such as coalesced access and using shared memory to minimize global memory reads.
 |
 |-- `sgemm`: Single-Precision General Matrix Multiplication (SGEMM) kernel.
 |
 |-- `sum-reduction`: Parallel sum reduction, i.e. computes the sum of an array of values.
```
