/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * epipg_terminate.c
 *
 * Code generation for function 'epipg_terminate'
 *
 */

/* Include files */
#include "epipg_terminate.h"
#include "_coder_epipg_mex.h"
#include "epipg_data.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void epipg_atexit(void)
{
  mexFunctionCreateRootTLS();
  emlrtEnterRtStackR2012b(emlrtRootTLSGlobal);
  emlrtLeaveRtStackR2012b(emlrtRootTLSGlobal);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
  emlrtExitTimeCleanup(&emlrtContextGlobal);
}

void epipg_terminate(void)
{
  emlrtLeaveRtStackR2012b(emlrtRootTLSGlobal);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/* End of code generation (epipg_terminate.c) */
