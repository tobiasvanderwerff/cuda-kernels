#include "soma.h"

#define CUDA_CALLABLE_FUNCTION __host__ __device__

CUDA_CALLABLE_FUNCTION mod_prec somaCalciumK(mod_prec prevV_soma, mod_prec prevCalcium_k) {
  mod_prec k_inf = (1 / (1 + exp(-1 * (prevV_soma + 61) / 4.2)));
  mod_prec tau_k = 1;
  mod_prec dk_dt = (k_inf - prevCalcium_k) / tau_k;
  mod_prec k_local = DELTA * dk_dt + prevCalcium_k;
  return k_local;
}

CUDA_CALLABLE_FUNCTION mod_prec somaCalciumL(mod_prec prevV_soma, mod_prec prevCalcium_l) {
  mod_prec l_inf = (1 / (1 + exp((prevV_soma + 85.5) / 8.5)));
  mod_prec tau_l = ((20 * exp((prevV_soma + 160) / 30) / (1 + exp((prevV_soma + 84) / 7.3))) + 35);
  mod_prec dl_dt = (l_inf - prevCalcium_l) / tau_l;
  mod_prec l_local = DELTA * dl_dt + prevCalcium_l;
  return l_local;
}

CUDA_CALLABLE_FUNCTION mod_prec somaSodiumM(mod_prec prevV_soma) {
  // RAT THALAMOCORTICAL SODIUM:
  mod_prec m_inf = 1 / (1 + (exp((-30 - prevV_soma) / 5.5)));
  mod_prec m_local = m_inf;
  return m_local;
}

CUDA_CALLABLE_FUNCTION mod_prec somaSodiumH(mod_prec prevV_soma, mod_prec prevSodium_h) {
  // RAT THALAMOCORTICAL SODIUM:
  mod_prec h_inf = 1 / (1 + (exp((-70 - prevV_soma) / -5.8)));
  mod_prec tau_h = 3 * exp((-40 - prevV_soma) / 33);
  mod_prec dh_dt = (h_inf - prevSodium_h) / tau_h;
  mod_prec h_local = prevSodium_h + DELTA * dh_dt;
  return h_local;
}

CUDA_CALLABLE_FUNCTION mod_prec somaPotassiumN(mod_prec prevV_soma, mod_prec prevPotassium_n) {
  // NEOCORTICAL
  mod_prec n_inf = 1 / (1 + exp((-3 - prevV_soma) / 10));
  mod_prec tau_n = 5 + (47 * exp(-(-50 - prevV_soma) / 900));
  mod_prec dn_dt = (n_inf - prevPotassium_n) / tau_n;
  mod_prec n_local = DELTA * dn_dt + prevPotassium_n;
  return n_local;
}

CUDA_CALLABLE_FUNCTION mod_prec somaPotassiumP(mod_prec prevV_soma, mod_prec prevPotassium_p) {
  // NEOCORTICAL
  mod_prec p_inf = 1 / (1 + exp((-51 - prevV_soma) / -12));
  mod_prec tau_p = 5 + (47 * exp(-(-50 - prevV_soma) / 900));
  mod_prec dp_dt = (p_inf - prevPotassium_p) / tau_p;
  mod_prec p_local = DELTA * dp_dt + prevPotassium_p;
  return p_local;
}

CUDA_CALLABLE_FUNCTION mod_prec somaPotassiumX(mod_prec prevV_soma, mod_prec prevPotassium_x_s) {
  // Voltage-dependent (fast) potassium
  mod_prec alpha_x_s = 0.13 * (prevV_soma + 25) / (1 - exp(-(prevV_soma + 25) / 10));
  mod_prec beta_x_s = 1.69 * exp(-0.0125 * (prevV_soma + 35));
  mod_prec x_inf_s = alpha_x_s / (alpha_x_s + beta_x_s);
  mod_prec tau_x_s = 1 / (alpha_x_s + beta_x_s);
  mod_prec dx_dt_s = (x_inf_s - prevPotassium_x_s) / tau_x_s;
  mod_prec x_s_local = 0.05 * dx_dt_s + prevPotassium_x_s;
  return x_s_local;
}

CUDA_CALLABLE_FUNCTION mod_prec somaCurrVolt(const SomaCurrVoltParams params) {
  // Get inputs
  mod_prec g_CaL = params.g_CaL;
  mod_prec prevV_dend = params.vDend;
  mod_prec prevV_soma = params.vSoma;
  mod_prec prevV_axon = params.vAxon;
  mod_prec k = params.k;
  mod_prec l = params.l;
  mod_prec m = params.m;
  mod_prec h = params.h;
  mod_prec n = params.n;
  mod_prec x_s = params.x_s;

  // SOMATIC CURRENTS

  // Dendrite-soma interaction current I_ds
  mod_prec I_ds = (G_INT / P1) * (prevV_soma - prevV_dend);
  // Inward low-threshold Ca current I_CaL
  mod_prec I_CaL = g_CaL * k * k * k * l * (prevV_soma - V_CA);  // k^3
  // Inward Na current I_Na_s
  mod_prec I_Na_s = G_NA_S * m * m * m * h * (prevV_soma - V_NA);
  // Leakage current I_ls
  mod_prec I_ls = G_LS * (prevV_soma - V_L);
  // Outward delayed potassium current I_Kdr
  mod_prec I_Kdr_s = G_KDR_S * n * n * n * n * (prevV_soma - V_K);  // SCHWEIGHOFER
  // I_K_s
  mod_prec I_K_s = G_K_S * pow(x_s, 4) * (prevV_soma - V_K);
  // Axon-soma interaction current I_as
  mod_prec I_as = (G_INT / (1 - P2)) * (prevV_soma - prevV_axon);

  mod_prec dVs_dt = (-(I_CaL + I_ds + I_as + I_Na_s + I_ls + I_Kdr_s + I_K_s)) / C_M;

  return DELTA * dVs_dt + prevV_soma;
}

// update somatic components
// SCHWEIGHOFER:
CUDA_CALLABLE_FUNCTION void _compSoma(CellCompParams *cellParamsPtr) {
  // only somaNew is modified, so that one is a pointer
  Soma *somaNew = &cellParamsPtr->newCellState->soma;
  Soma somaPrev = cellParamsPtr->prevCellState->soma;

  // Set up parameters
  mod_prec prevV_soma = somaPrev.V_soma;
  mod_prec prevCalcium_k = somaPrev.Calcium_k;
  mod_prec prevCalcium_l = somaPrev.Calcium_l;
  mod_prec prevSodium_h = somaPrev.Sodium_h;
  mod_prec prevPotassium_n = somaPrev.Potassium_n;
  mod_prec prevPotassium_p = somaPrev.Potassium_p;
  mod_prec prevPotassium_x_s = somaPrev.Potassium_x_s;

  // Compute
  somaNew->Calcium_k = somaCalciumK(prevV_soma, prevCalcium_k);
  somaNew->Calcium_l = somaCalciumL(prevV_soma, prevCalcium_l);

  somaNew->Sodium_m = somaSodiumM(prevV_soma);
  somaNew->Sodium_h = somaSodiumH(prevV_soma, prevSodium_h);

  somaNew->Potassium_n = somaPotassiumN(prevV_soma, prevPotassium_n);
  somaNew->Potassium_p = somaPotassiumP(prevV_soma, prevPotassium_p);

  somaNew->Potassium_x_s = somaPotassiumX(prevV_soma, prevPotassium_x_s);

  // Set up params for final computation
  SomaCurrVoltParams params;
  params.g_CaL = somaPrev.g_CaL;
  params.vSoma = somaPrev.V_soma;
  params.vDend = cellParamsPtr->prevCellState->dend.V_dend;
  params.vAxon = cellParamsPtr->prevCellState->axon.V_axon;
  params.k = somaNew->Calcium_k;
  params.l = somaNew->Calcium_l;
  params.m = somaNew->Sodium_m;
  params.h = somaNew->Sodium_h;
  params.n = somaNew->Potassium_n;
  params.x_s = somaNew->Potassium_x_s;

  // Final computation
  somaNew->V_soma = somaCurrVolt(params);
}

void compSoma(CellCompParams *cellParamsPtr) {
  _compSoma(cellParamsPtr);
}

__global__ void compSomaCUDA(CellCompParams *cellParamsPtr) {
  _compSoma(cellParamsPtr);
}
