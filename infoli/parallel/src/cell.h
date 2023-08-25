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

#include "util/util.h"

// Cell morphology
#define P1 0.25     // Cell surface ratio soma/dendrite (0.2)
#define P2 0.15     // Cell surface ratio axon(hillock)/soma (0.1)
#define G_INT 0.13  // Cell internal conductance (0.13)
// Reversal potentials
#define V_NA 55   // Na reversal potential (55)
#define V_K -75   // K reversal potential
#define V_CA 120  // Ca reversal potential (120)
#define V_H -43   // H current reversal potential
#define V_L 10    // leak current

// Capacitance
#define C_M 1
// 0.05 milli sec = 50 micro sec
#define DELTA 0.05

typedef struct Dend {
  mod_prec V_dend;
  mod_prec Hcurrent_q;
  mod_prec Calcium_r;
  mod_prec Potassium_s;
  mod_prec I_CaH;
  mod_prec Ca2Plus;
} Dend;

typedef struct Soma {
  mod_prec g_CaL;
  mod_prec V_soma;
  mod_prec Sodium_m;
  mod_prec Sodium_h;
  mod_prec Calcium_k;
  mod_prec Calcium_l;
  mod_prec Potassium_n;
  mod_prec Potassium_p;
  mod_prec Potassium_x_s;
} Soma;

typedef struct Axon {
  mod_prec V_axon;
  mod_prec Sodium_m_a;
  mod_prec Sodium_h_a;
  mod_prec Potassium_x_a;
} Axon;

/* struct that represents a cell
 *
 */
typedef struct CellState {
  int cellID;
  struct Dend dend;
  struct Soma soma;
  struct Axon axon;
} CellState;

typedef struct CellCompParams {
  mod_prec iAppIn;
  mod_prec *neighVdend;
  mod_prec *neighConductances;
  int *neighId;
  CellState *prevCellState;
  CellState *newCellState;
  int index_of_neighVdend;
  int total_amount_of_neighbours;
} CellCompParams;

/*** FUNCTION PROTOTYPES ***/

void communicationStep(CellCompParams *params, const CellState *cells, int cellCount);

void initState(CellState *cellPtr, int cellCount);
