/******************************************************************************
 *cr
 *cr            (C) Copyright 2010 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ******************************************************************************/

#include <stdio.h>

#define TILE_SIZE 16

__global__ void mysgemm(int m, int n, int k, const float *A, const float *B, float* C) {

    /********************************************************************
     *
     * Compute C = A x B
     *   where A is a (m x k) matrix
     *   where B is a (k x n) matrix
     *   where C is a (m x n) matrix
     *
     * Use shared memory for tiling
     *
     ********************************************************************/

    __shared__ float tileA[TILE_SIZE][TILE_SIZE];
    __shared__ float tileB[TILE_SIZE][TILE_SIZE];

    // Load reused variables into registers for fast access
    const unsigned int tx = threadIdx.x;
    const unsigned int ty = threadIdx.y;

    // Figure out the position to fill, i.e. C[row][col]
    const unsigned int row = blockIdx.y * blockDim.y + ty;
    const unsigned int col = blockIdx.x * blockDim.x + tx;

    float sum = 0;
    for (int ti = 0; ti < (k-1)/TILE_SIZE + 1; ++ti) {  // ti is the phase number
      // Collaboratively allocate memory
      // Assign 0s to obviate a control branch when calculating the partial result,
      // avoiding a branch divergence.
      if ((TILE_SIZE*ti + tx < k) && (row < m)) {
        tileA[ty][tx] = A[row*k + TILE_SIZE*ti + tx];   // fixed row
      } else {
        tileA[ty][tx] = 0;
      }
      if ((TILE_SIZE*ti + ty < k) && (col < n)) {
        tileB[ty][tx] = B[(TILE_SIZE*ti + ty)*k + col]; // fixed column
      } else {
        tileB[ty][tx] = 0;
      }
      __syncthreads();
      // Calculate partial result
      if (row < m && col < n) {
        for (int ki = 0; ki < TILE_SIZE; ++ki) {
          /* if (TILE_SIZE*ti + ki < k) */  // avoid this branch divergence
          sum += tileA[ty][ki] * tileB[ki][tx];
        }
      }
      __syncthreads();
    }
    // Write final result
    if (row < m && col < n)
      C[row*k + col] = sum;
}

void basicSgemm(char transa, char transb, int m, int n, int k, float alpha, const float *A, int lda, const float *B, int ldb, float beta, float *C, int ldc)
{
    if ((transa != 'N') && (transa != 'n')) {
	printf("unsupported value of 'transa'\n");
    	return;
    }

    if ((transb != 'N') && (transb != 'n')) {
	printf("unsupported value of 'transb'\n");
	return;
    }

    if ((alpha - 1.0f > 1e-10) || (alpha - 1.0f < -1e-10)) {
	printf("unsupported value of alpha\n");
	return;
    }

    if ((beta - 0.0f > 1e-10) || (beta - 0.0f < -1e-10)) {
	printf("unsupported value of beta\n");
	return;
    }

    // Initialize thread block and kernel grid dimensions ---------------------

    const unsigned int BLOCK_SIZE = TILE_SIZE;

    dim3 gridDim((n-1)/BLOCK_SIZE + 1, (m-1)/BLOCK_SIZE + 1, 1);
    dim3 blockDim(BLOCK_SIZE, BLOCK_SIZE, 1);


    // Invoke CUDA kernel -----------------------------------------------------
    
    mysgemm<<<gridDim, blockDim>>>(m, n, k, A, B, C);


}


