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
    // Element 1
    if (start + tx < size) {
        partial_sum[tx] = in[start + tx];
    } else {
        partial_sum[tx] = 0;
    }
    // Element 2
    if (start + tx + bd < size) {
        partial_sum[tx + bd] = in[start + tx + bd];
    } else {
        partial_sum[tx + bd] = 0;
    }
    __syncthreads();

    // Traverse the reduction tree until all elements are added together
    // Ensure minimal amount of thread divergence by keeping adjacent threads active.
    for (unsigned int stride = BLOCK_SIZE; stride >= 1; stride /= 2)  {
        __syncthreads();
        if (tx < stride) {
            partial_sum[tx] += partial_sum[tx + stride];
        }
    }

    // The final result is the first element of the shared memory array
    if (tx == 0)  // Avoid race conditions
        out[bx] = partial_sum[0];
}
