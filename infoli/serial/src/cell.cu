//
// Created by niels on 5/31/21.
//

#include "cell.h"

void initState(CellState *cellPtr, int cellCount) {
  CellState initState;
  // Initial dendritic parameters
  initState.dend.V_dend = -60;
  initState.dend.Calcium_r = 0.0112788;    // High-threshold calcium
  initState.dend.Potassium_s = 0.0049291;  // Calcium-dependent potassium
  initState.dend.Hcurrent_q = 0.0337836;   // H current
  initState.dend.Ca2Plus = 3.7152;         // Calcium concentration
  initState.dend.I_CaH = 0.5;              // High-threshold calcium current
  // Initial somatic parameters
  initState.soma.g_CaL =
      0.68;  // default arbitrary value but it should be randomized per cell
  initState.soma.V_soma = -60;
  initState.soma.Sodium_m = 1.0127807;  // Sodium (artificial)
  initState.soma.Sodium_h = 0.3596066;
  initState.soma.Potassium_n = 0.2369847;  // Potassium (delayed rectifier)
  initState.soma.Potassium_p = 0.2369847;
  initState.soma.Potassium_x_s = 0.1;    // Potassium (voltage-dependent)
  initState.soma.Calcium_k = 0.7423159;  // Low-threshold calcium
  initState.soma.Calcium_l = 0.0321349;
  // Initial axonal parameters
  initState.axon.V_axon = -60;
  // sisaza: Sodium_m_a doesn't have a state, therefore this assignment
  // doesn'thave any effect
  initState.axon.Sodium_m_a = 0.003596066;  // Sodium (thalamocortical)
  initState.axon.Sodium_h_a = 0.9;
  initState.axon.Potassium_x_a = 0.2369847;  // Potassium (transient)

  initState.cellID = 0;  // Compute the cellID of the first cell in this core

  // Copy init state to all cell states
  for (int j = 0; j < cellCount; j++) {
    memcpy(&cellPtr[j], &initState, sizeof(CellState));
    initState.cellID++;  // next cell's ID is increased by 1
  }
}

/* The function where the magic happens: using the structures created during
 * initialization, exchange necessary dendritic voltage
 */
void communicationStep(CellCompParams *params, const CellState *cells,
                       int cellCount) {
  for (int j = 0; j < cellCount; j++) {
    CellCompParams param = params[j];
    for (int i = 0; i < param.total_amount_of_neighbours; ++i) {
      param.neighVdend[i] = cells[param.neighId[i]].dend.V_dend;
    }
  }
}