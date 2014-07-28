// This file was automatically generated by Kmisc::registerFunctions()

#include <R.h>
#include <Rinternals.h>

#include <R_ext/Rdynload.h>

SEXP Algorithm5(SEXP restrictTable_,SEXP numDrivers_,SEXP numGenes_,SEXP typeCBN_,SEXP birthRate_, SEXP s_, SEXP death_,SEXP mu_,SEXP initSize_,SEXP sampleEvery_,SEXP detectionSize_,SEXP finalTime_,SEXP initSize_species_,SEXP initSize_iter_,SEXP seed_gsl_,SEXP verbose_,SEXP speciesFS_,SEXP ratioForce_,SEXP typeFitness_,SEXP maxram_,SEXP mutatorGenotype_,SEXP initMutant_,SEXP maxWallTime_,SEXP keepEvery_,SEXP alpha_,SEXP sh_,SEXP K_,SEXP endTimeEvery_,SEXP finalDrivers_);

R_CallMethodDef callMethods[]  = {
  {"C_Algorithm5", (DL_FUNC) &Algorithm5, 29},
  {NULL, NULL, 0}
};

void R_init_OncoSimulR(DllInfo *info) {
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}

