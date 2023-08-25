#pragma once

#include "../cell.h"
#include "../util/util.h"

// Axon hillock conductances (mS/cm2)

// Na gate conductance (according to literature: 100 to 200 times as big as
// somatic conductance)
#define G_NA_A 240
// Na (resurgent) gate conductance
#define G_NA_R 0
// K voltage-dependent
#define G_K_A 20
// Leak conductance
#define G_LA 0.016

typedef struct AxonCurrVoltParams {
  mod_prec vSoma;
  mod_prec vAxon;
  mod_prec m_a, h_a, x_a;
  mod_prec newVAxon;
} AxonCurrVoltParams;

void compAxon(CellCompParams *cellParamsPtr);
