/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_epipg_mex.c
 *
 * Code generation for function '_coder_epipg_mex'
 *
 */

/* Include files */
#include "_coder_epipg_mex.h"
#include "_coder_epipg_api.h"
#include "epipg_data.h"
#include "epipg_initialize.h"
#include "epipg_terminate.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs,
                 const mxArray *prhs[])
{
  mexAtExit(&epipg_atexit);
  /* Module initialization. */
  epipg_initialize();
  /* Dispatch the entry-point. */
  unsafe_epipg_mexFunction(nlhs, plhs, nrhs, prhs);
  /* Module termination. */
  epipg_terminate();
}

emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLSR2021a(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1,
                           NULL);
  return emlrtRootTLSGlobal;
}

void unsafe_epipg_mexFunction(int32_T nlhs, mxArray *plhs[5], int32_T nrhs,
                              const mxArray *prhs[26])
{
  const mxArray *outputs[5];
  int32_T b_nlhs;
  /* Check for proper number of arguments. */
  if (nrhs < 26) {
    emlrtErrMsgIdAndTxt(emlrtRootTLSGlobal,
                        "EMLRT:runTime:TooFewInputsConstants", 9, 4, 5, "epipg",
                        4, 5, "epipg", 4, 5, "epipg");
  }
  if (nrhs != 26) {
    emlrtErrMsgIdAndTxt(emlrtRootTLSGlobal, "EMLRT:runTime:WrongNumberOfInputs",
                        5, 12, 26, 4, 5, "epipg");
  }
  if (nlhs > 5) {
    emlrtErrMsgIdAndTxt(emlrtRootTLSGlobal,
                        "EMLRT:runTime:TooManyOutputArguments", 3, 4, 5,
                        "epipg");
  }
  /* Call the function. */
  epipg_api(prhs, nlhs, outputs);
  /* Copy over outputs to the caller. */
  if (nlhs < 1) {
    b_nlhs = 1;
  } else {
    b_nlhs = nlhs;
  }
  emlrtReturnArrays(b_nlhs, &plhs[0], &outputs[0]);
}

/* End of code generation (_coder_epipg_mex.c) */
