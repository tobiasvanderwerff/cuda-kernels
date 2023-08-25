#include "dendrite.h"

mod_prec dendHCurr(mod_prec prevV_dend, mod_prec prevHcurrent_q) {
  // Update dendritic H current component
  mod_prec q_inf = 1 / (1 + exp((prevV_dend + 80) / 4));
  mod_prec tau_q = 1 / (exp(-0.086 * prevV_dend - 14.6) + exp(0.070 * prevV_dend - 1.87));
  mod_prec dq_dt = (q_inf - prevHcurrent_q) / tau_q;
  mod_prec q_local = DELTA * dq_dt + prevHcurrent_q;
  return q_local;
}

mod_prec dendCaCurr(mod_prec prevV_dend, mod_prec prevCalcium_r) {
  // Update dendritic high-threshold Ca current component
  mod_prec alpha_r = 1.7 / (1 + exp(-(prevV_dend - 5) / 13.9));
  mod_prec beta_r = 0.02 * (prevV_dend + 8.5) / (exp((prevV_dend + 8.5) / 5) - 1);
  mod_prec r_inf = alpha_r / (alpha_r + beta_r);
  mod_prec tau_r = 5 / (alpha_r + beta_r);
  mod_prec dr_dt = (r_inf - prevCalcium_r) / tau_r;
  mod_prec r_local = DELTA * dr_dt + prevCalcium_r;
  return r_local;
}

mod_prec dendKCurr(mod_prec prevPotassium_s, mod_prec prevCa2Plus) {
  // Update dendritic Ca-dependent K current component
  mod_prec alpha_s = min((0.00002 * prevCa2Plus), 0.01);
  mod_prec beta_s = 0.015;
  mod_prec s_inf = alpha_s / (alpha_s + beta_s);
  mod_prec tau_s = 1 / (alpha_s + beta_s);
  mod_prec ds_dt = (s_inf - prevPotassium_s) / tau_s;
  mod_prec s_local = DELTA * ds_dt + prevPotassium_s;
  return s_local;
}

// Consider merging cal into kCurr since cal's output doesn't go to currVolt but
// to kCurr
mod_prec dendCal(mod_prec prevCa2Plus, mod_prec prevI_CaH) {
  // update Calcium concentration
  mod_prec dCa_dt = -3 * prevI_CaH - 0.075 * prevCa2Plus;
  mod_prec Ca2Plus_local = DELTA * dCa_dt + prevCa2Plus;
  return Ca2Plus_local;
}

mod_prec icNeighbours(mod_prec *neighVdend, mod_prec *neighConductances, mod_prec prevV_dend, int neighbors) {
  mod_prec I_c = 0;
  for (int i = 0; i < neighbors; i++) {
    mod_prec V = prevV_dend - neighVdend[i];
    mod_prec f = 0.8 * exp(-1 * pow(V, 2) / 100) + 0.2;  // SCHWEIGHOFER 2004 VERSION
    mod_prec cond = neighConductances[i];
    I_c += (cond * f * V);
  }

  return I_c;
}

mod_prec dendCICa(mod_prec prevV_dend, mod_prec r) {
  // DENDRITIC CURRENTS
  // Inward high-threshold Ca current I_CaH
  mod_prec I_CaH = G_CAH * r * r * (prevV_dend - V_CA);
  return I_CaH;
}

mod_prec dendCurrVolt(const DendCurrVoltParams params) {
  // Get inputs
  mod_prec I_c = params.iC;
  mod_prec I_app = params.iApp;
  mod_prec prevV_dend = params.vDend;
  mod_prec prevV_soma = params.vSoma;
  mod_prec q = params.q;
  mod_prec s = params.s;
  mod_prec I_CaH = params.I_CaH;

  // DENDRITIC CURRENTS

  // Soma-dendrite interaction current I_sd
  mod_prec I_sd = (G_INT / (1 - P1)) * (prevV_dend - prevV_soma);
  // Outward Ca-dependent K current I_K_Ca
  mod_prec I_K_Ca = G_K_CA * s * (prevV_dend - V_K);
  // Leakage current I_ld
  mod_prec I_ld = G_LD * (prevV_dend - V_L);
  // Inward anomalous rectifier I_h
  mod_prec I_h = G_H * q * (prevV_dend - V_H);

  mod_prec dVd_dt = (-(I_CaH + I_sd + I_ld + I_K_Ca + I_c + I_h) + I_app) / C_M;

  return DELTA * dVd_dt + prevV_dend;
}

void compDendrite(CellCompParams *cellParamsPtr, const int randomness) {
  // only dendNew is modified, so that one is a pointer
  Dend *dendNew = &cellParamsPtr->newCellState->dend;
  Dend dendPrev = cellParamsPtr->prevCellState->dend;

  // Set up parameters
  mod_prec prevV_dend = dendPrev.V_dend;
  mod_prec prevHcurrent_q = dendPrev.Hcurrent_q;
  mod_prec prevCalcium_r = dendPrev.Calcium_r;
  mod_prec prevPotassium_s = dendPrev.Potassium_s;
  mod_prec prevCa2Plus = dendPrev.Ca2Plus;
  mod_prec prevI_CaH = dendPrev.I_CaH;

  dendNew->Hcurrent_q = dendHCurr(prevV_dend, prevHcurrent_q);
  dendNew->Calcium_r = dendCaCurr(prevV_dend, prevCalcium_r);
  dendNew->Potassium_s = dendKCurr(prevPotassium_s, prevCa2Plus);
  dendNew->Ca2Plus = dendCal(prevCa2Plus, prevI_CaH);

  /* ANDREAS change here the last parameter of icNeighbours, we no longer use
   * number_of_neighbours instead we have an index in struct cellParameters
   * which tells us where the cell's neighbouring voltages end in the array
   */

  // Set up params for final computation
  DendCurrVoltParams params;
  // in random initialization mode, cells run in closed circuit mode so
  // neighboring current equals zero
  if (randomness == 1) {
    params.iC = 0;
  } else {
    params.iC = icNeighbours(cellParamsPtr->neighVdend, cellParamsPtr->neighConductances, dendPrev.V_dend,
                             cellParamsPtr->total_amount_of_neighbours);
  }

  dendNew->I_CaH = dendCICa(prevV_dend, dendNew->Calcium_r);

  params.iApp = cellParamsPtr->iAppIn;
  params.vDend = dendPrev.V_dend;
  params.newVDend = dendNew->V_dend;
  params.vSoma = cellParamsPtr->prevCellState->soma.V_soma;
  params.q = dendNew->Hcurrent_q;
  params.r = dendNew->Calcium_r;
  params.s = dendNew->Potassium_s;
  params.I_CaH = dendNew->I_CaH;

  // Final computation
  dendNew->V_dend = dendCurrVolt(params);
}
