/******************************************************************************
 *cr
 *cr            (C) Copyright 2010 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ******************************************************************************/

#define BLOCK_SIZE 512

__global__ void reduction(float *out, float *in, unsigned size)
{
    /********************************************************************
    Load a segment of the input vector into shared memory
    Traverse the reduction tree
    Write the computed sum to the output vector at the correct index
    ********************************************************************/

    __shared__ float partial_sum[2*BLOCK_SIZE];

    const unsigned int tx = threadIdx.x, bx = blockIdx.x, bd = blockDim.x;
    const unsigned int start = 2 * bx * bd;

    // Let each thread load 2 elements into shared memory (coalesced),
    // making sure to fill up the entire shared array.
    partial_sum[tx] = start + tx < size ? in[start + tx] : 0;
    partial_sum[tx + bd] = start + tx + bd < size ? in[start + tx + bd] : 0;
    __syncthreads();

    // Traverse the reduction tree until all elements are added together
    // This implementation is slightly inefficient because after only a few iterations
    // of this loop, only a small percentage of threads in each warp is actually doing
    // work.
    for (unsigned int stride = 1; stride <= BLOCK_SIZE; stride *= 2) {
        __syncthreads(); // wait for all threads in the block before continuing
        if (tx % stride == 0)  {
            partial_sum[2*tx] += partial_sum[2*tx + stride];
        }
    }

    // The final result is the first element of the shared memory array
    if (tx == 0)  // Avoid race conditions
        out[bx] = partial_sum[0];
}
