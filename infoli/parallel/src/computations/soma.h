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

// Somatic conductances (mS/cm2)

// Na gate conductance (=90 in Schweighofer code, 70 in paper) 120 too little
#define G_NA_S 150
// K delayed rectifier gate conductance (alternative value: 18)
#define G_KDR_S 9.0
// Voltage-dependent (fast) potassium
#define G_K_S 5
// Leak conductance (0.015)
#define G_LS 0.016

typedef struct SomaCurrVoltParams {
  mod_prec g_CaL;
  mod_prec vSoma;
  mod_prec vDend;
  mod_prec vAxon;
  mod_prec k, l, m, h, n, x_s;
  mod_prec newVSoma;
} SomaCurrVoltParams;

void compSoma(CellCompParams *cellParamsPtr);
__global__ void compSomaCUDA(CellCompParams *cellParamsPtr);
