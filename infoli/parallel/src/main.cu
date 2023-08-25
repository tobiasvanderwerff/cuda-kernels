/*
 *
 * Copyright (c) 2012, Neurasmus B.V., The Netherlands,
 * web: www.neurasmus.com email: voltage@neurasmus.com
 *
 * Any use or reproduction in whole or in parts is prohibited
 * without the written consent of the copyright owner.
 *
 * All Rights Reserved.
 *
 *
 * Author: Sebastian Isaza
 * Created: 19-01-2012
 * Modified: 07-08-2012
 *
 * Description: Top source file of the Inferior Olive model, originally written
 * in Matlab by Jornt De Gruijl. It contains the implementation of all
 * functions. The main function allocates the necessary memory, initializes the
 * system state and runs the model calculations.
 *
 */

/*
 * we assume that dim1 refers to the number of ROWS
 * and dim2 refers to COLUMNS (cell network size).
 */
#include <stdio.h>
#include <stdlib.h>

#include "infoli_simple.h"

#define CONFIG_FILE_NAME "config/cellConnections.txt"

int main(int argc, char *argv[]) {
  /* Process command line arguments
   * Case argc = 2 then a one-pulse input is stimulated.
   */
  if (argc == 1) {
    printf("Error: Too few arguments.\nUsage: ./infoli_simple <Net_Size>\n");
    exit(EXIT_FAILURE);
  } else if (argc > 2) {
    printf("Error: Too many arguments.\nUsage: /infoli_simple <Net_Size>\n");
    exit(EXIT_FAILURE);
  }

  int cellCount = atoi(argv[1]);

  // Allocate space for network
  CellState **cellPtr = allocCellPtr(cellCount);
  CellCompParams *cellParamsPtr = allocCellParams(cellCount);

  // CUDA Setup 
  CellCompParams* cellParamsPtr_d = allocCellParamsCUDA(cellCount);
  CellState* cellPtr_d = allocCellPtrCUDA(cellCount);

  cudaDeviceSynchronize();


  // Init network
  char conFile[BUFF_SIZE];
  sprintf(conFile, CONFIG_FILE_NAME);
  // TODO <TVDW>: Go over 'init' again after doing 'simulate'
  init(conFile, cellParamsPtr, cellPtr, cellCount);

  // Perform simulation
  simulate(cellParamsPtr, cellPtr, cellCount);

  // Clean up
  free(cellPtr[0]);
  free(cellPtr[1]);
  free(cellPtr);
  free(cellParamsPtr);

  cudaFree(cellParamsPtr_d);
  cudaFree(cellPtr_d);

  return 0;
}
