# Optimizations

## Naive implementation
**Runtime over 5 runs**: 14.504ms


## Optimization 1: coalesced global memory access
**Description**: Ensure that threads within a warp are accessing global memory
contiguously among threads, by using a strided access pattern. E.g. thread i accesses
idx i, and thread i+1 accesses idx i+1 of global memory. This is faster because of
so-called DRAM bursing, which allows obtaining several contiguous global memory elements
and sharing them among threads, using only a single read.

**Runtime over 5 runs**: 4.25ms
**Speed-up vs naive**:   3.4x  
