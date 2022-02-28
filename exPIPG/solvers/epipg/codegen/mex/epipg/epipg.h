/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * epipg.h
 *
 * Code generation for function 'epipg'
 *
 */

#pragma once

/* Include files */
#include "rtwtypes.h"
#include "emlrt.h"
#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Function Declarations */
void epipg(real_T q[800], real_T u[398], real_T v[796], real_T omg, real_T N,
           real_T M, real_T *exitflag, real_T *solve_time);

/* End of code generation (epipg.h) */
