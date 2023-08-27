#include "axon.h"

CUDA_CALLABLE_FUNCTION mod_prec axonSodiumH(mod_prec prevV_axon, mod_prec prevSodium_h_a) {
  // Update axonal Na components
  // NOTE: current has shortened inactivation to account for high
  // firing frequencies in axon hillock
  mod_prec h_inf_a = 1 / (1 + (exp((-60 - prevV_axon) / -5.8)));
  mod_prec tau_h_a = 1.5 * exp((-40 - prevV_axon) / 33);
  mod_prec dh_dt_a = (h_inf_a - prevSodium_h_a) / tau_h_a;
  mod_prec h_a_local = prevSodium_h_a + DELTA * dh_dt_a;
  // Put result
  return h_a_local;
}

CUDA_CALLABLE_FUNCTION mod_prec axonSodiumM(mod_prec prevV_axon) {
  // keep separate variable for readability sake
  mod_prec m_inf_a = 1 / (1 + (exp((-30 - prevV_axon) / 5.5)));
  return m_inf_a;
}

CUDA_CALLABLE_FUNCTION mod_prec axonPotassium(mod_prec prevV_axon, mod_prec prevPotassium_x_a) {
  // D'ANGELO 2001 -- Voltage-dependent potassium
  mod_prec alpha_x_a = 0.13 * (prevV_axon + 25) / (1 - exp(-(prevV_axon + 25) / 10));
  mod_prec beta_x_a = 1.69 * exp(-0.0125 * (prevV_axon + 35));
  mod_prec alpha_beta_x_a = alpha_x_a + beta_x_a;
  mod_prec x_inf_a = alpha_x_a / alpha_beta_x_a;
  mod_prec tau_x_a = 1 / alpha_beta_x_a;
  mod_prec dx_dt_a = (x_inf_a - prevPotassium_x_a) / tau_x_a;
  mod_prec x_a_local = 0.05 * dx_dt_a + prevPotassium_x_a;
  // Put result
  return x_a_local;
}

CUDA_CALLABLE_FUNCTION mod_prec axonCurrVolt(const AxonCurrVoltParams params) {
  // Get inputs
  mod_prec prevV_soma = params.vSoma;
  mod_prec prevV_axon = params.vAxon;
  mod_prec m_a = params.m_a;
  mod_prec h_a = params.h_a;
  mod_prec x_a = params.x_a;

  // AXONAL CURRENTS
  // Sodium
  mod_prec I_Na_a = G_NA_A * m_a * m_a * m_a * h_a * (prevV_axon - V_NA);
  // Leak
  mod_prec I_la = G_LA * (prevV_axon - V_L);
  // Soma-axon interaction current I_sa
  mod_prec I_sa = (G_INT / P2) * (prevV_axon - prevV_soma);
  // Potassium (transient)
  mod_prec I_K_a = G_K_A * pow(x_a, 4) * (prevV_axon - V_K);
  mod_prec dVa_dt = (-(I_K_a + I_sa + I_la + I_Na_a)) / C_M;

  return DELTA * dVa_dt + prevV_axon;
}

// update somatic components
// SCHWEIGHOFER:
CUDA_CALLABLE_FUNCTION void compAxon(CellCompParams *cellParamsPtr) {
  // only axonNew is modified, so that one is a pointer
  Axon *axonNew = &cellParamsPtr->newCellState->axon;
  Axon axonPrev = cellParamsPtr->prevCellState->axon;

  // Set up parameters
  mod_prec prevV_axon = axonPrev.V_axon;
  mod_prec prevSodium_h_a = axonPrev.Sodium_h_a;
  mod_prec prevPotassium_x_a = axonPrev.Potassium_x_a;

  // Compute
  axonNew->Sodium_h_a = axonSodiumH(prevV_axon, prevSodium_h_a);
  axonNew->Sodium_m_a = axonSodiumM(prevV_axon);
  axonNew->Potassium_x_a = axonPotassium(prevV_axon, prevPotassium_x_a);

  // Set up params for final computation
  AxonCurrVoltParams params;
  params.vSoma = cellParamsPtr->prevCellState->soma.V_soma;
  params.vAxon = axonPrev.V_axon;
  params.m_a = axonNew->Sodium_m_a;
  params.h_a = axonNew->Sodium_h_a;
  params.x_a = axonNew->Potassium_x_a;

  // Final computation
  axonNew->V_axon = axonCurrVolt(params);
}