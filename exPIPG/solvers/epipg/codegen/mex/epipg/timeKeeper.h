/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * timeKeeper.h
 *
 * Code generation for function 'timeKeeper'
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
void b_timeKeeper(real_T *outTime_tv_sec, real_T *outTime_tv_nsec);

void savedTime_not_empty_init(void);

void timeKeeper(const emlrtTimespec newTime);

/* End of code generation (timeKeeper.h) */
