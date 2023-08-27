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
#include "../util/util.h"

// Dendritic conductances (mS/cm2)

// Potassium gate conductance (35)
#define G_K_CA 35
// High-threshold Ca gate conductance (4.5)
#define G_CAH 4.5
// Dendrite leak conductance (0.015)
#define G_LD 0.016
// H current gate conductance (1.5) (0.15 in SCHWEIGHOFER 2004)
#define G_H 0.125

typedef struct DendCurrVoltParams {
  mod_prec iApp;
  mod_prec iC;
  mod_prec vDend;
  mod_prec vSoma;
  mod_prec q, r, s;
  mod_prec newVDend;
  mod_prec I_CaH;
} DendCurrVoltParams;

void compDendrite(CellCompParams *cellParamsPtr, int randomness);
__global__ void compDendriteCUDA(CellCompParams *cellParamsPtr, int randomness);
