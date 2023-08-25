/*
 *
 * Copyright (c) 2012, Neurasmus B.V., The Netherlands,
 * web: www.neurasmus.com email: info@neurasmus.com
 *
 * Any use reproduction in whole or in parts is prohibited
 * without the written consent of the copyright owner.
 *
 * All Rights Reserved.
 *
 *
 * Author: Sebastian Isaza
 * Created: 10-04-2012
 * Modified: 06-06-2012
 *
 * Description : Top header file of the Inferior Olive model. It contains the
 * constant model conductances, the data structures that hold the cell state and
 * the function prototypes.
 *
 */

#pragma once

#include <stdio.h>
#include <stdlib.h>

#include "cell.h"

/*** MACROS ***/
//	flag enabling dumping of the final states of the neurons (search for it
// in results folder)
#define PRINT_STATE 1
#define PRINT_STATE_LOCATION "results/lastStateDump.txt"

// make it zero to facilitate debugging , 0 for debugging / 1 for random states
#define RAND_INIT 0

// in ms, for when no input file is provided , time in msec when you do not have
// an input file
#define SIMTIME 1000
//	2 integer, the dot, 2 decimals and the delimiter
#define IAPP_MAX_CHARS 6
//	flag enabling or disabling user defined initial settings of gcal
#define G_CAL_FROM_FILE 0
#define ALLTOALL 0
// Cell properties , biological properties for the cells ( irrevelant for now )

// Conductance for neighbors' coupling
#define CONDUCTANCE 0.004

/*** TYPEDEFS AND STRUCTS***/
CellState **allocCellPtr(int cellCount);
CellState *allocAndCopyCellPtrCUDA(int cellCount, CellState** cellPtr);
CellCompParams *allocCellParams(int cellCount);
CellCompParams *allocAndCopyCellParamsCUDA(int cellCount, CellCompParams* cellParamsPtr); 

void init(const char *conFile, CellCompParams *cellParamsPtr, CellState **cellPtr, int cellCount);
// void simulate(CellCompParams *cellParamsPtr, CellState *cellPtr, int cellCount);
void simulate(CellCompParams *cellParamsPtr, CellState *cellPtr, 
              CellCompParams *cellParamsPtr_h, CellState **cellPtr_h,
              int cellCount);