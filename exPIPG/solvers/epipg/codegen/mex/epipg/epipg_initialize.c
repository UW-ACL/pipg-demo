/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * epipg_initialize.c
 *
 * Code generation for function 'epipg_initialize'
 *
 */

/* Include files */
#include "epipg_initialize.h"
#include "_coder_epipg_mex.h"
#include "epipg_data.h"
#include "rt_nonfinite.h"
#include "timeKeeper.h"

/* Function Declarations */
static void epipg_once(void);

/* Function Definitions */
static void epipg_once(void)
{
  mex_InitInfAndNan();
  savedTime_not_empty_init();
}

void epipg_initialize(void)
{
  mexFunctionCreateRootTLS();
  emlrtClearAllocCountR2012b(emlrtRootTLSGlobal, false, 0U, NULL);
  emlrtEnterRtStackR2012b(emlrtRootTLSGlobal);
  if (emlrtFirstTimeR2012b(emlrtRootTLSGlobal)) {
    epipg_once();
  }
}

/* End of code generation (epipg_initialize.c) */
