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
#pragma once

#include "../cell.h"
#include "util.h"

int readFileLine(FILE *pInFile, mod_prec *iAppArray, int cellCount);

void readGCalFromFile(CellState *cellPtr, int cellCount);

void printState(const CellState *cellPtr, const char *filename, int cellCount);
