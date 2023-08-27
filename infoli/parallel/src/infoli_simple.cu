
#include "computations/axon.h"
#include "computations/dendrite.h"
#include "computations/soma.h"
#include "infoli_simple.h"
#include "util/io.h"

/**
 * Allocates memory for the cellPtr
 *
 * CellPtr is a 2-D array of size 2*CellCount
 * and containd cell states used in every simulation
 * step accordingly
 */
CellState **allocCellPtr(int cellCount) {
  CellState **cellPtr;
  printf("Allocating Neuron State Memory\n");
  cellPtr = (CellState **)malloc(2 * sizeof(CellState *));
  if (cellPtr == NULL) {
    printf("Error: Couldn't malloc for cellPtr\n");
    exit(EXIT_FAILURE);
  }

  cellPtr[0] = (CellState *)malloc(cellCount * sizeof(CellState));
  cellPtr[1] = (CellState *)malloc(cellCount * sizeof(CellState));
  if ((!cellPtr[0]) || (!cellPtr[1])) {
    printf("Error: Couldn't malloc the array for CellStates\n");
    exit(EXIT_FAILURE);
  }
  return cellPtr;
}

CellState *allocAndCopyCellPtrCUDA(int cellCount, CellState** cellPtr) {
  // Note that instead of returning a 2D array (double pointer) with dimensions
  // arr[2][cellCount], we return a 1D array (single pointer) of size
  // 2*cellCount. This is because it's annoying to allocate a double pointer on
  // the CUDA device, and I couldn't figure out how to do it properly.

  CellState* cellPtr_d;
  cudaError_t err;

  const size_t rowsize = size_t(cellCount) * sizeof(CellState);
  err = cudaMalloc((void**) &cellPtr_d, 2*rowsize);
  cudaSuccessOrExit(err);
  // Copy first row
  err = cudaMemcpy(cellPtr_d, cellPtr[0], rowsize, cudaMemcpyHostToDevice);
  cudaSuccessOrExit(err);
  // Copy second row
  err = cudaMemcpy(cellPtr_d+cellCount, cellPtr[1], rowsize, cudaMemcpyHostToDevice);
  cudaSuccessOrExit(err);

  return cellPtr_d;
}

/**
 * Allocates memory for the cellParamsPtr
 *
 * CellCompParams struct is used to update a cell state.
 * It contains the current flowing to this specific cell
 * voltages from communicating cells , and pointers to the previous
 * cell state and to the next cell state
 * We allocate cellCount CellCompParams for all the core's cells
 * Pointer to previous and next states are pointers to CellPtr[0] and
 * CellPtr[1] elements.
 */
CellCompParams *allocCellParams(int cellCount) {
  CellCompParams *cellParamsPtr;
  cellParamsPtr = (CellCompParams *)malloc(cellCount * sizeof(CellCompParams));
  if (cellParamsPtr == NULL) {
    printf("Error: Couldn't malloc for cellParamsPtr\n");
    exit(EXIT_FAILURE);
  }
  return cellParamsPtr;
}

CellCompParams *allocAndCopyCellParamsCUDA(int cellCount, CellCompParams* cellParamsPtr) { 
  cudaError_t err;
  CellCompParams* cellParamsPtr_d;

  const unsigned int memsize = cellCount * sizeof(CellCompParams);
  err = cudaMalloc((void**) &cellParamsPtr_d, memsize);
  cudaSuccessOrExit(err);
  err = cudaMemcpy(cellParamsPtr_d, cellParamsPtr, memsize, cudaMemcpyHostToDevice);
  cudaSuccessOrExit(err);

  return cellParamsPtr_d;
}

/**
 * Sets the initial state of the simulation
 * Does this by reading from the file specified in conFile
 */
void init(const char *conFile, CellCompParams *cellParamsPtr, CellState **cellPtr, int cellCount) {
  // initial amount of neighbours for each cell is 0 (bug fix in case the cell
  // stays isolated)
  for (int i = 0; i < cellCount; i++) {
    cellParamsPtr[i].total_amount_of_neighbours = 0;
  }

  mod_prec cond_value = CONDUCTANCE;
  if (ALLTOALL == 0) {
    FILE *pConFile = fopen(conFile, "r");
    if (pConFile == NULL) {
      printf("Error: Couldn't open file %s\n", conFile);
      exit(EXIT_FAILURE);
    }
    // reading into a temporary float, because fscanf apparently doesn't
    // read into mod_prec
    float temp_cond_value;
    fscanf(pConFile, "%f ", &temp_cond_value);
    fclose(pConFile);
    cond_value = (mod_prec)temp_cond_value;
  }

  // handle connectivity file parsing so that each cell knows what it needs
  printf("Reading Network Connectivity Data\n");
  for (int line_counter = 0; line_counter < cellCount; line_counter++) {
    for (int i = 0; i < cellCount; i++) {
      // this connection is considered not existing if conductance = 0
      if (cond_value != 0) {
        // part of the code handling RECEIVING and noting which of my cells
        // needs input from which other cells, from ANY core

        CellCompParams *currParams = &cellParamsPtr[i];
        int n = currParams->total_amount_of_neighbours;

        // if this is the first neighbour, initialize buffers
        if (currParams->total_amount_of_neighbours == 0) {
          currParams->neighVdend = NULL;
          currParams->neighConductances = NULL;
          currParams->neighId = NULL;
        }

        currParams->neighId = allocIntArr(currParams->neighId, n);
        // which cell sends this voltage to us (GLOBAL ID)
        currParams->neighId[n] = line_counter;
        currParams->neighConductances = allocModArr(currParams->neighConductances, n);
        // what conductance we use to calculate its impact
        currParams->neighConductances[n] = cond_value;
        // allocate space for storing this voltage
        currParams->neighVdend = allocModArr(currParams->neighVdend, n);
        currParams->total_amount_of_neighbours++;
      }
    }
  }
  printf("int: %d, dend: %d, soma: %d, axon: %d, cellstate: %d\n", sizeof(int), sizeof(Dend), sizeof(Soma),
         sizeof(Axon), sizeof(CellState));
  // Initialise cellPtr[0] with appropriate values
  initState(cellPtr[0], cellCount);
  if (G_CAL_FROM_FILE) {
    readGCalFromFile(cellPtr[0], cellCount);
  }

  // Initialize g_CaL
  int seedvar = time(NULL);
  srand(seedvar);

  for (int i = 0; i < cellCount; i++) {
    cellPtr[1][i].soma.g_CaL = cellPtr[0][i].soma.g_CaL;
    if (RAND_INIT) {
      cellPtr[0][i].soma.g_CaL = 0.6 + (0.2 * (rand() % 100) / 100);
      cellPtr[1][i].soma.g_CaL = cellPtr[0][i].soma.g_CaL;
    }
  }

  // random initialization process
  if (RAND_INIT) {
    int initSteps;
    for (int i = 0; i < cellCount; i++) {
      initSteps = rand() % (int)ceil(100 / DELTA);
      // make it odd, so that the final state is in prevCellState
      initSteps = initSteps | 0x00000001;

      CellCompParams *currParams = &cellParamsPtr[i];
      for (int j = 0; j < initSteps; j++) {
        // Arrange inputs
        currParams->iAppIn = 0;  // No stimulus
        currParams->prevCellState = &cellPtr[j % 2][i];
        currParams->newCellState = &cellPtr[(j % 2) ^ 1][i];

        compDendrite(currParams, 1);
        compSoma(currParams);
        compAxon(currParams);
      }
    }
  }
}

void copyDeviceToHost(
    CellCompParams *cellParamsPtr_d, CellState *cellPtr_d, 
    CellCompParams *cellParamsPtr_h, CellState **cellPtr_h,
    int cellCount) {
  // Copy device->host
  cudaMemcpy(cellParamsPtr_h, cellParamsPtr_d, cellCount * sizeof(CellCompParams), cudaMemcpyDeviceToHost);
  cudaMemcpy(cellPtr_h[0], cellPtr_d, cellCount * sizeof(CellState), cudaMemcpyDeviceToHost);
  cudaMemcpy(cellPtr_h[1], cellPtr_d + cellCount, cellCount * sizeof(CellState), cudaMemcpyDeviceToHost);
  cudaDeviceSynchronize();
}


void copyHostToDevice(
    CellCompParams *cellParamsPtr_d, CellState *cellPtr_d, 
    CellCompParams *cellParamsPtr_h, CellState **cellPtr_h,
    int cellCount) {
  cudaMemcpy(cellParamsPtr_d, cellParamsPtr_h, cellCount * sizeof(CellCompParams), cudaMemcpyHostToDevice);
  cudaMemcpy(cellPtr_d, cellPtr_h[0], cellCount * sizeof(CellState), cudaMemcpyHostToDevice);
  cudaMemcpy(cellPtr_d + cellCount, cellPtr_h[1], cellCount * sizeof(CellState), cudaMemcpyHostToDevice);
  cudaDeviceSynchronize();
}

__global__ void update_cells(CellCompParams *cellParamsPtr, CellState *cellPtr, int cellCount, mod_prec iApp, int simArrayId) {
  /* Update all cells in parallel. */

  const int targetCell = blockIdx.x*blockDim.x + threadIdx.x;

  // The reason we divide by 4 instead of 3 is that this way, all sequences of 4
  // elements will be in the same block, and can synchronize. This is because 4
  // will always be a multiple of the block size (assuming the block size is a
  // power of 2!), whereas 3 is not. If you divide by 3, you could have a
  // situation where e.g. the first 2 elements of a sequence are in block i, but
  // the last element is in block i+1 (which means you can't synchronize them).
  if (targetCell / 4 < cellCount) {
    CellCompParams *currParams = &cellParamsPtr[targetCell / 4];
    // TODO: move to shared memory?
    // TODO: try coalesced access (ie strided)

    /* we simulate a hardcoded input pulse here
      * that differs from step to step
      */
    if (targetCell % 4 == 0) {
      currParams->iAppIn = iApp;
      currParams->prevCellState = &cellPtr[simArrayId*cellCount + targetCell/4];
      currParams->newCellState = &cellPtr[(simArrayId ^ 1)*cellCount + targetCell/4];
    }

    __syncthreads();
    switch (targetCell % 4) {
      case 0:
        compDendrite(currParams, 0);
        break;
      case 1:
        compSoma(currParams);
        break;
      case 2:
        compAxon(currParams);
        break;
      case 3:
        break;
    }
  }
}

/**
 * MAIN Loop
 *
 * Iterates over all time steps in the simulation and performs the various
 * calculations This is the actual iteration
 */
// void performSimulation(CellCompParams *cellParamsPtr, CellState *cellPtr, int cellCount, int totalSimSteps) {
void performSimulation(CellCompParams *cellParamsPtr_d, CellState *cellPtr_d, 
                       int cellCount, int totalSimSteps) {

  const unsigned int gridDimComm = CEILDIV(cellCount*cellCount, CUDA_BLOCK_SIZE);
  const unsigned int gridDimCell = CEILDIV(4*cellCount, CUDA_BLOCK_SIZE);

  for (int simStep = 0; simStep < totalSimSteps; simStep++) {
    int simArrayId = simStep % 2;
    mod_prec iApp = ((simStep >= 20000) && (simStep < 20500 - 1)) ? 6 : 0;

    /* Perform_Communication() performs the inter core
     * core dendrite communication with connected cells
     * See definition for more details
     */
    communicationStep<<<gridDimComm, CUDA_BLOCK_SIZE>>>(cellParamsPtr_d, cellPtr_d + simArrayId*cellCount, cellCount);

    cudaDeviceSynchronize();

    update_cells<<<gridDimCell, CUDA_BLOCK_SIZE>>>(cellParamsPtr_d, cellPtr_d, cellCount, iApp, simArrayId);

    cudaDeviceSynchronize();
  }
}

/**
 * The simulation
 *
 * Starts and times the simulation
 */
// void simulate(CellCompParams *cellParamsPtr, CellState *cellPtr, int cellCount) {
void simulate(CellCompParams *cellParamsPtr_d, CellState *cellPtr_d, 
              CellCompParams *cellParamsPtr_h, CellState **cellPtr_h,
              int cellCount) {

  printf("Beginning Execution\n");
  if (ALLTOALL == 1) {
    printf("All-to-All Connectivity Benchmark Active\n");
  }
  timestamp_t t0 = getTimeStamp();

  int totalSimSteps = (int)(SIMTIME / DELTA);

  // Perform simulation
  copyHostToDevice(cellParamsPtr_d, cellPtr_d, cellParamsPtr_h, cellPtr_h, cellCount);
  performSimulation(cellParamsPtr_d, cellPtr_d, cellCount, totalSimSteps);
  copyDeviceToHost(cellParamsPtr_d, cellPtr_d, cellParamsPtr_h, cellPtr_h, cellCount);

  timestamp_t t1 = getTimeStamp();
  printf(
      "Execution completed: %ld ms of brain time of %ld neurons in %ld "
      "simulation steps executed\n",
      SIMTIME, cellCount, totalSimSteps);
  printf("Application Execution Time : %lld usecs\n", t1 - t0);

  if (PRINT_STATE) {
    char resultFileName[BUFF_SIZE];
    sprintf(resultFileName, PRINT_STATE_LOCATION);
    // simStep%2 here should refer to the cellPtr which has the last state of
    // the network that we calculated
    // printState(cellPtr_h + (totalSimSteps % 2)*cellCount, resultFileName, cellCount);
    printState(cellPtr_h[totalSimSteps%2], resultFileName, cellCount);
  }
}