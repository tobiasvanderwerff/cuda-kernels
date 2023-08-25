
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

  FILE *pConFile = fopen(conFile, "r");
  if (pConFile == NULL) {
    printf("Error: Couldn't open file %s\n", conFile);
    exit(EXIT_FAILURE);
  }
  // handle connectivity file parsing so that each cell knows what it needs
  printf("Reading Network Connectivity Data\n");
  for (int line_counter = 0; line_counter < cellCount; line_counter++) {
    for (int i = 0; i < cellCount; i++) {
      mod_prec cond_value;
      if (ALLTOALL == 0) {
        // reading into a temporary float, because fscanf apparently doesn't
        // read into mod_prec
        float temp_cond_value;
        fscanf(pConFile, "%f ", &temp_cond_value);
        cond_value = (mod_prec)temp_cond_value;

      } else {
        cond_value = CONDUCTANCE;
      }
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

  fclose(pConFile);
}

/**
 * MAIN Loop
 *
 * Iterates over all time steps in the simulation and performs the various
 * calculations This is the actual iteration
 */
void performSimulation(CellCompParams *cellParamsPtr, CellState **cellPtr, int cellCount, int totalSimSteps) {
  for (int simStep = 0; simStep < totalSimSteps; simStep++) {
    mod_prec iApp;
    int simArrayId = simStep % 2;
    if ((simStep >= 20000) && (simStep < 20500 - 1)) {
      iApp = 6;
    } else {
      iApp = 0;
    }

    /* Perform_Communication() performs the inter core
     * core dendrite communication with connected cells
     * See definition for more details
     */
    communicationStep(cellParamsPtr, cellPtr[simArrayId], cellCount);

    for (int targetCell = 0; targetCell < cellCount; targetCell++) {
      CellCompParams *currParams = &cellParamsPtr[targetCell];

      /* we simulate a hardcoded input pulse here
       * that differs from step to step
       */
      currParams->iAppIn = iApp;
      currParams->prevCellState = &cellPtr[simArrayId][targetCell];
      currParams->newCellState = &cellPtr[simArrayId ^ 1][targetCell];

      compDendrite(currParams, 0);
      compSoma(currParams);
      compAxon(currParams);
    }
  }
}

/**
 * The simulation
 *
 * Starts and times the simulation
 */
void simulate(CellCompParams *cellParamsPtr, CellState **cellPtr, int cellCount) {
  printf("Beginning Execution\n");
  if (ALLTOALL == 1) {
    printf("All-to-All Connectivity Benchmark Active\n");
  }
  timestamp_t t0 = getTimeStamp();

  int totalSimSteps = (int)(SIMTIME / DELTA);

  performSimulation(cellParamsPtr, cellPtr, cellCount, totalSimSteps);

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
    printState(cellPtr[totalSimSteps % 2], resultFileName, cellCount);
  }
}
