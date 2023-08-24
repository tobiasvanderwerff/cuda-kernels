# Optimizations

Note that optimizations are cumulative, e.g. optimization 2 also includes optimization 1.
All performance measurements are done using an input array with 1,000,000 values (unless
mentioned otherwise), randomly generated.

Hardware used: NVIDIA Geforce GTX 770

## Naive implementation
[Commit link](https://github.com/tobiasvanderwerff/cuda-kernels/commit/bb35ec6f14a345a729e291bee342d29160795c38)

### Performance
**Runtime over 5 runs**: 14.504ms


## Optimization 1: coalesced global memory access
[Commit link](https://github.com/tobiasvanderwerff/cuda-kernels/commit/f50db76b0433cad50bb8d39becad545500bacc40)

### Description
Ensure that threads within a warp are accessing global memory
contiguously among threads, by using a strided access pattern. E.g. thread i accesses
idx i, and thread i+1 accesses idx i+1 of global memory. This is faster because of
so-called DRAM bursing, which allows obtaining several contiguous global memory elements
and sharing them among threads, using only a single read.

### Performance
**Runtime over 5 runs**: 4.25ms  
**Speed-up vs naive**:   3.4x  


## Optimization 2: use more threads and less work per thread
[Commit link](https://github.com/tobiasvanderwerff/cuda-kernels/commit/e92e441b615b08bd3134c62269547c51107835f8)

### Description
The naive implementation started with each thread computing a fixed
number of output elements: ELEMENTS_PER_THREAD. However, this does not make full
utilization of the available GPU processors, since the later threads end up doing
nothing, since each previous thread is already taking care of ELEMENTS_PER_THREAD
elements. Therefore, I modified the stride in such a way that as many threads as
possible are performing useful work.

### Performance
**Runtime over 5 runs**: 0.283ms  
**Speed-up vs naive**: 51.3x


## Optimization 3: shared memory (privatization)
[Commit link](https://github.com/tobiasvanderwerff/cuda-kernels/commit/09ea898ad203ede867e10dc53da9120f515dd058)

### Description
Instead of all threads trying to write directly to global memory (which requires atomic
adds, which means many threads will be stalled), each thread block maintains a local
version of the histogram (an example of the
[privatization](https://en.wikipedia.org/wiki/Privatization_(computer_programming))
technique). Each thread in a block first writes to the block-specific histogram, leading
to less traffic and higher throughput on those local histograms. Then, all blocks merge
their partial histograms into the global memory histogram. This leads to fewer total
writes to global memory and less contention among threads. The reduction of writes to
global memory is by a factor of THREADS_PER_BLOCK, e.g. 256. Note that the actual
increase in throughput due to fewer threads waiting for atomic operations to finish
depends on the number of bins in the histogram.

Downside is that it is necessary to set the number of bins as a constant since it has to
be known at compile time.

It's interesting to note from the performance that it doesn't lead to an improvement
compared to optimization 2 (17x vs 51x). It seems that the additional work required for
initializing shared memory, using two barrier synchronizations, and doing the final
shared->global transfer, lead to some noticeable overhead. I imagine that there would be
a bigger advantage to the shared memory approach where there is more traffic on
individual bins, i.e. when the total number of bins is reduced. This can be seen from
the additional runs shown below, where optimization 3 _does_ outperform opt2 when the
bins are reduced to 32.

```
-- num_elements=100,000,000, num_bins=4096  
opt2: 13.17ms  
opt3: 52.75ms  

-- num_elements=100,000,000, num_bins=32  
opt2: 30.94ms
opt3: 15.49ms
```

The conclusion seems to be that whether to prefer optimization 2 or 3 is
situation-dependent, depending in particular on the number of bins of the output
histogram.

### Performance
**Runtime over 5 runs**: 0.819ms  
**Speed-up vs naive**: 17.7x
