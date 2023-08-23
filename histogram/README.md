# Optimizations

Note that optimizations are cumulative, i.e. optimization 2 also includes optimization 1.

## Naive implementation
[Commit link](https://github.com/tobiasvanderwerff/cuda-kernels/commit/bb35ec6f14a345a729e291bee342d29160795c38)

**Runtime over 5 runs**: 14.504ms


## Optimization 1: coalesced global memory access
[Commit link](https://github.com/tobiasvanderwerff/cuda-kernels/commit/f50db76b0433cad50bb8d39becad545500bacc40)
https://github.com/tobiasvanderwerff/cuda-kernels/commit/f50db76b0433cad50bb8d39becad545500bacc40

**Description**: Ensure that threads within a warp are accessing global memory
contiguously among threads, by using a strided access pattern. E.g. thread i accesses
idx i, and thread i+1 accesses idx i+1 of global memory. This is faster because of
so-called DRAM bursing, which allows obtaining several contiguous global memory elements
and sharing them among threads, using only a single read.

**Runtime over 5 runs**: 4.25ms
**Speed-up vs naive**:   3.4x  


## Optimization 2: use more threads and less work per thread

**Description**: The naive implementation started with each thread computing a fixed
number of output elements: ELEMENTS_PER_THREAD. However, this does not make full
utilization of the available GPU processors, since the later threads end up doing
nothing, since each previous thread is already taking care of ELEMENTS_PER_THREAD
elements. Therefore, I modified the stride in such a way that as many threads as
possible are performing useful work.

**Runtime over 5 runs**: 0.283ms
**Speed-up vs naive**: 51.3x


