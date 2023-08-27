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

#include <stdio.h>
#include <stdlib.h>

#define MIN(A, B) A > B ? B : A
#define CEILDIV(A, B) (A-1)/B + 1
#define BUFF_SIZE 100

// CUDA block size
#define CUDA_BLOCK_SIZE 512

// typedef double mod_prec;
// BE VERY CAREFUL TO CHECK ALL SCANFS TO BE SURE YOU SCAN FOR DOUBLE-POINT
// ACCURACY, KNOWN ISSUE WITH COND VALUES) AND MPI_TYPES
typedef float mod_prec;
typedef unsigned long long timestamp_t;

timestamp_t getTimeStamp();

int *allocIntArr(int *pointer, int existing_slots);
void cudaSuccessOrExit(cudaError_t err, const char* file, int line);

mod_prec *allocModArr(mod_prec *pointer, int existing_slots);
