/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * tic.c
 *
 * Code generation for function 'tic'
 *
 */

/* Include files */
#include "tic.h"
#include "rt_nonfinite.h"
#include "timeKeeper.h"
#include "emlrt.h"

/* Function Definitions */
void tic(void)
{
  emlrtTimespec t;
  emlrtClockGettimeMonotonic(&t);
  timeKeeper(t);
}

/* End of code generation (tic.c) */
