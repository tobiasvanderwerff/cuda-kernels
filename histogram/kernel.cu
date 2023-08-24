/******************************************************************************
 *cr
 *cr            (C) Copyright 2010 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ******************************************************************************/
#define NTHREADS_HOST 256
#define NUM_BINS 4096

// Define your kernels in this file you may use more than one kernel if you
// need to

__global__ void histogram_kernel_interleaved_shared(unsigned int* input, unsigned int* bins,
             unsigned int num_elements, unsigned int num_bins) {
    __shared__ unsigned int Sbins[NUM_BINS];  // private copy of the histogram

    const unsigned int tx = threadIdx.x, bx = blockIdx.x, 
                       bdim = blockDim.x, gdim = gridDim.x;
    const unsigned int start = bx*bdim + tx;
    const unsigned int stride = bdim*gdim;

    // Set shared memory values
    for (int i = tx; i < NUM_BINS; i += bdim)
        Sbins[i] = 0;
    __syncthreads();

    // Update histogram values
    int i = start;
    while (i < num_elements) {
        unsigned int val = input[i];
        if (val < NUM_BINS) {
            atomicAdd(&(Sbins[val]), 1);
        }
        i += stride;
    }

    // Merge partial histogram into global memory
    __syncthreads();
    for (int i = tx; i < NUM_BINS; i += bdim) {
        atomicAdd(&(bins[i]), Sbins[i]);
    }
}


__global__ void convert_kernel(unsigned int *bins32, uint8_t *bins8,
            unsigned int num_bins) {
	
    const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i < num_bins) {
        bins8[i] = (bins32[i] <= 255) ? (uint8_t)bins32[i] : 255u;
    }
}

/******************************************************************************
Setup and invoke your kernel(s) in this function. You may also allocate more
GPU memory if you need to
*******************************************************************************/
void histogram(unsigned int* input, uint8_t* bins, unsigned int num_elements,
        unsigned int num_bins) {
    
    // Create 32 bit bins
    cudaError_t cuda_ret;
    unsigned int* bins32;
    // Allocate memory
    cuda_ret = cudaMalloc((void**) &bins32, num_bins*sizeof(unsigned int));
    if (cuda_ret != cudaSuccess) FATAL("Failed to set bins32 elements to 0");
    cudaDeviceSynchronize();
    // Fill allocated memory with 0 values
    cuda_ret = cudaMemset(bins32, 0, num_bins*sizeof(unsigned int));
    if (cuda_ret != cudaSuccess) FATAL("Failed to set bins32 elements to 0");
    cudaDeviceSynchronize();

    Timer timer;
    startTime(&timer);
	
    // Launch histogram kernel using 32-bit bins
    int nthreads = NTHREADS_HOST;
    int nblocks = (num_elements-1)/nthreads + 1;
    histogram_kernel_interleaved_shared<<<nblocks, nthreads>>>(input, bins32, num_elements, num_bins); 
	
    cudaDeviceSynchronize();

    // Convert 32-bit bins into 8-bit bins
    nthreads = NTHREADS_HOST;
    nblocks = (num_bins-1)/nthreads + 1;
    convert_kernel<<<nthreads, nblocks>>>(bins32, bins, num_bins);
	
    stopTime(&timer); printf("\nActual time interleaved shared= %f s\n", elapsedTime(timer));

    // Free allocated device memory
    cudaFree(bins32);
}


