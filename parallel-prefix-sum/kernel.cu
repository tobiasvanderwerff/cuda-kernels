/******************************************************************************
 *cr
 *cr            (C) Copyright 2010 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ******************************************************************************/

/* 
 * O(N) implementation of parallel prefix scan, using a balanced tree.
 * See https://developer.nvidia.com/gpugems/gpugems3/part-vi-gpu-computing/chapter-39-parallel-prefix-sum-scan-cuda
 * for implementation details.
 */

#define BLOCK_SIZE 512

#include <math.h>
#include "support.h"

__global__ void partial_scan(float *out, float *in, float *sums, unsigned int in_size) {
    /* Perform parallel prefix scan on a block-by-block basis. */

    __shared__ float temp[2*BLOCK_SIZE];

    const unsigned int tx = threadIdx.x, bx = blockIdx.x, bd = blockDim.x;
    const unsigned int start = 2 * bx * bd;
    const unsigned int n = 2 * BLOCK_SIZE;

    // Load the input array from global into shared memory
    temp[tx]    = (start+tx < in_size) ? in[start+tx] : 0;
    temp[tx+bd] = (start+tx+bd < in_size) ? in[start+tx+bd] : 0;

    // TODO (optimization): no. of threads within a warp are reduced by a factor of 2
    // each iteration, which is inefficient. Fix this.
    // TODO (optimization): address bank conflicts by using padding in 'temp'. See the
    // Nvidia doc mentioned at the top of this file for an example using a
    // CONFLICT_FREE_OFFSET macro. Something to keep in mind: It is mentioned in this video
    // (https://youtu.be/CZgM3DEBplE?t=830) that bank conflicts are not always a big
    // deal, since shared memory reads are so fast, and the scheduler can also optimize
    // it away by switching to another warp while the banks figure out the requested
    // data

    // Up-sweep (reduce) phase
    for (int i = 0; i < log2((double)n); i++) {
        __syncthreads();
        int pw = powf(2, i);
        if (tx % pw == 0) {
            int ix1 = 2*tx + pw - 1;
            int ix2 = 2*tx + 2*pw - 1;
            temp[ix2] += temp[ix1];
        }
    }

    __syncthreads();
    if (tx == 0) {
        // Store the total sum
        int ix = (start+n < in_size) ? n-1 : in_size % n - 1;
        sums[bx] = temp[ix];
        // Set last element to 0 for down-sweep
        temp[n-1] = 0;
    }

    // Down-sweep phase
    for (int i = log2((double)n)-1; i >= 0; i--) {
        __syncthreads();
        int pw = powf(2, i);
        if (tx % pw == 0) {
            int ix1 = 2*tx + pw - 1;
            int ix2 = 2*tx + 2*pw - 1;
            float t = temp[ix1];
            temp[ix1] = temp[ix2];
            temp[ix2] += t;
        }
    }

    // Write the result
    __syncthreads();
    if (start+tx < in_size)
        out[start+tx] = temp[tx];
    if (start+tx+bd < in_size)
        out[start+tx+bd] = temp[tx+bd];
}


__global__ void uniform_add(float* out, float* in, float* incr, int n) {
    /* Add values in 'incr' to 'in' on a block-by-block basis. */
    const int bx = blockIdx.x, bd = blockDim.x, tx = threadIdx.x;
    const int i = bd * bx + tx;
    if (i < n) {
        out[i] = in[i] + incr[bx];
    }
}


/******************************************************************************
Setup and invoke your kernel(s) in this function. You may also allocate more
GPU memory if you need to
*******************************************************************************/
void preScan(float *out, float *in, unsigned int in_size)
{
    int nblocks;
    float *sums, *incr;

    nblocks = (in_size-1) / (BLOCK_SIZE*2) + 1;

    cudaError_t cuda_ret;
    cuda_ret = cudaMalloc((void**)&sums, nblocks*sizeof(float));
    if (cuda_ret != cudaSuccess) FATAL("Failed to allocate memory for 'sums'");
    cuda_ret = cudaMalloc((void**)&incr, nblocks*sizeof(float));
    if (cuda_ret != cudaSuccess) FATAL("Failed to allocate memory for 'incr'");

    unsigned int gridDim, blockDim;

    // Calculate prefix sum on a block-by-block basis
    gridDim = nblocks;
    blockDim = BLOCK_SIZE;
    partial_scan<<< gridDim, blockDim >>>(out, in, sums, in_size);

    cudaDeviceSynchronize();

    // Calculate increments
    gridDim = (nblocks-1)/BLOCK_SIZE + 1;
    blockDim = BLOCK_SIZE;
    partial_scan<<< gridDim, blockDim  >>>(incr, sums, sums, nblocks);

    cudaDeviceSynchronize();

    // Add increments to block results
    gridDim = nblocks;
    blockDim = 2*BLOCK_SIZE;
    uniform_add<<< gridDim, blockDim >>>(out, out, incr, in_size);

    cudaDeviceSynchronize();

    cudaFree(sums);
    cudaFree(incr);
}

