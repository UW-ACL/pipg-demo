/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * epipg.c
 *
 * Code generation for function 'epipg'
 *
 */

/* Include files */
#include "epipg.h"
#include "rt_nonfinite.h"
#include "tic.h"
#include "toc.h"
#include "blas.h"
#include "mwmathutil.h"
#include <stddef.h>
#include <stdio.h>
#include <string.h>

/* Function Definitions */
void epipg(real_T q[800], real_T u[398], real_T v[796], real_T omg, real_T N,
           real_T M, real_T *exitflag, real_T *solve_time)
{
  static const real_T a[16] = {-1.0, -0.0, -0.0, -0.0, -0.2, -1.0, -0.0, -0.0,
                               -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.2, -1.0};
  static const real_T b_a[16] = {1.0, 0.2, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                                 0.0, 0.0, 1.0, 0.2, 0.0, 0.0, 0.0, 1.0};
  static const real_T e_a[16] = {1.0, 0.0, 0.0, 0.0, 0.2, 1.0, 0.0, 0.0,
                                 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.2, 1.0};
  static const real_T c_a[8] = {0.020000000000000004, 0.0, 0.2, 0.0, 0.0,
                                0.020000000000000004, 0.0, 0.2};
  static const real_T d_a[8] = {
      -0.020000000000000004, -0.2, -0.0, -0.0, -0.0, -0.0,
      -0.020000000000000004, -0.2};
  static const real_T g_a[8] = {0.020000000000000004, 0.2, 0.0, 0.0, 0.0, 0.0,
                                0.020000000000000004, 0.2};
  static const int8_T h_a[16] = {1, 0, 0, 0, 0, 1, 0, 0,
                                 0, 0, 1, 0, 0, 0, 0, 1};
  static const int8_T b_varargin_1[4] = {6, 1, 6, 1};
  static const int8_T f_a[4] = {1, 0, 0, 1};
  static const int8_T varargin_1[4] = {-6, -1, -6, -1};
  ptrdiff_t incx_t;
  ptrdiff_t n_t;
  real_T x[1198];
  real_T q2[800];
  real_T qp[800];
  real_T vp[796];
  real_T u2[398];
  real_T up[398];
  real_T absx;
  real_T alf;
  real_T b_varargin_2;
  real_T bet;
  real_T conv_test_dual;
  real_T sig;
  real_T sig1;
  real_T varargin_2;
  int32_T i;
  int32_T i1;
  int32_T i2;
  int32_T i3;
  int32_T k;
  int32_T qp_tmp;
  int32_T t;
  boolean_T conv_flag;
  boolean_T exitg1;
  /*                                                                                                                 %
   * main variables (modified) */
  /*               % constants */
  /*  barebones epipg which is not verbose and doesn't compute cost at each
   * iteration */
  /*  test for dual infeasibility is turned off */
  /*  additional copies of primal and dual variables */
  /*  start epipg timing */
  tic();
  sig = 1.0;
  /*  power iteration */
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 15000)) {
    sig1 = sig;
    i = (int32_T)(N - 1.0);
    for (t = 0; t < i; t++) {
      for (i1 = 0; i1 < 4; i1++) {
        i2 = i1 << 2;
        i3 = i1 << 1;
        v[t + 199 * i1] = (q[(t + 200 * i1) + 1] -
                           (((q[t] * b_a[i2] + q[t + 200] * b_a[i2 + 1]) +
                             q[t + 400] * b_a[i2 + 2]) +
                            q[t + 600] * b_a[i2 + 3])) -
                          (u[t] * c_a[i3] + u[t + 199] * c_a[i3 + 1]);
      }
    }
    for (i = 0; i < 4; i++) {
      q[200 * i] = 0.0;
      i1 = i << 2;
      q[200 * i] += v[0] * a[i1];
      q[200 * i] += v[199] * a[i1 + 1];
      q[200 * i] += v[398] * a[i1 + 2];
      q[200 * i] += v[597] * a[i1 + 3];
    }
    for (i = 0; i < 2; i++) {
      u[199 * i] = 0.0;
      i1 = i << 2;
      u[199 * i] += v[0] * d_a[i1];
      u[199 * i] += v[199] * d_a[i1 + 1];
      u[199 * i] += v[398] * d_a[i1 + 2];
      u[199 * i] += v[597] * d_a[i1 + 3];
    }
    i = (int32_T)((N - 1.0) + -1.0);
    for (t = 0; t < i; t++) {
      for (i1 = 0; i1 < 4; i1++) {
        i2 = i1 << 2;
        q[(t + 200 * i1) + 1] =
            v[t + 199 * i1] -
            (((v[t + 1] * e_a[i2] + v[t + 200] * e_a[i2 + 1]) +
              v[t + 399] * e_a[i2 + 2]) +
             v[t + 598] * e_a[i2 + 3]);
      }
      for (i1 = 0; i1 < 2; i1++) {
        i2 = (t + 199 * i1) + 1;
        u[i2] = 0.0;
        i3 = i1 << 2;
        u[i2] += v[t + 1] * d_a[i3];
        u[i2] += v[t + 200] * d_a[i3 + 1];
        u[i2] += v[t + 399] * d_a[i3 + 2];
        u[i2] += v[t + 598] * d_a[i3 + 3];
      }
    }
    q[(int32_T)N - 1] = v[(int32_T)N - 2];
    q[(int32_T)N + 199] = v[(int32_T)N + 197];
    q[(int32_T)N + 399] = v[(int32_T)N + 396];
    q[(int32_T)N + 599] = v[(int32_T)N + 595];
    i = 0;
    i1 = 0;
    for (i2 = 0; i2 < 800; i2++) {
      q2[i2] = q[i1 + 200 * i];
      i++;
      if (i > 3) {
        i = 0;
        i1++;
      }
    }
    i = 0;
    i1 = 0;
    for (i2 = 0; i2 < 398; i2++) {
      u2[i2] = u[i1 + 199 * i];
      i++;
      if (i > 1) {
        i = 0;
        i1++;
      }
    }
    memcpy(&x[0], &q2[0], 800U * sizeof(real_T));
    memcpy(&x[800], &u2[0], 398U * sizeof(real_T));
    n_t = (ptrdiff_t)1198;
    incx_t = (ptrdiff_t)1;
    sig = dnrm2(&n_t, &x[0], &incx_t);
    if (muDoubleScalarAbs(sig1 - sig) < 0.001) {
      exitg1 = true;
    } else {
      for (i = 0; i < 800; i++) {
        q[i] /= sig;
      }
      for (i = 0; i < 398; i++) {
        u[i] /= sig;
      }
      k++;
    }
  }
  sig *= 1.1;
  /*  pipg step size */
  alf = 2.0 / (muDoubleScalarSqrt(4.0 * omg * sig + 1.0) + 1.0);
  bet = omg * alf;
  /*  pipg diagnostic variables */
  /*  conv_test_primal_prev = coder.nullcopy(0.0);     */
  conv_flag = false;
  *exitflag = 0.0;
  /*   0 -> Maximum iterations reached; desired solve or infeasibility detection
   * accuracy not met */
  /*   1 -> Primal and Dual feasible (solved to desired accuracy) */
  /*  -1 -> Primal infeasible */
  /*  -2 -> Dual infeasible (not available) */
  /*  pipg iterations */
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 50000)) {
    /*  copy of the previous iterate, used for monitoring convergence */
    /*  proj grad: t = 1 */
    qp[0] = 0.0;
    qp[200] = -1.0;
    qp[400] = 0.0;
    qp[600] = 0.0;
    /*  projection onto initial value set */
    for (i = 0; i < 2; i++) {
      i1 = i << 1;
      i2 = i << 2;
      up[199 * i] = muDoubleScalarMax(
          -2.0, muDoubleScalarMin(
                    2.0, u[199 * i] -
                             alf * ((u[0] * (real_T)f_a[i1] +
                                     u[199] * (real_T)f_a[i1 + 1]) -
                                    (((v[0] * g_a[i2] + v[199] * g_a[i2 + 1]) +
                                      v[398] * g_a[i2 + 2]) +
                                     v[597] * g_a[i2 + 3]))));
    }
    /*  proj grad on u */
    /*  proj grad: t in [2, M-1] */
    i = (int32_T)((M - 1.0) + -1.0);
    for (t = 0; t < i; t++) {
      for (i1 = 0; i1 < 4; i1++) {
        i2 = i1 << 2;
        qp_tmp = (t + 200 * i1) + 1;
        qp[qp_tmp] = muDoubleScalarMax(
            varargin_1[i1],
            muDoubleScalarMin(
                b_varargin_1[i1],
                q[qp_tmp] -
                    alf * (((((q[t + 1] * (real_T)h_a[i2] +
                               q[t + 201] * (real_T)h_a[i2 + 1]) +
                              q[t + 401] * (real_T)h_a[i2 + 2]) +
                             q[t + 601] * (real_T)h_a[i2 + 3]) +
                            v[t + 199 * i1]) -
                           (((v[t + 1] * e_a[i2] + v[t + 200] * e_a[i2 + 1]) +
                             v[t + 399] * e_a[i2 + 2]) +
                            v[t + 598] * e_a[i2 + 3]))));
      }
      for (i1 = 0; i1 < 2; i1++) {
        i2 = i1 << 1;
        i3 = i1 << 2;
        qp_tmp = (t + 199 * i1) + 1;
        up[qp_tmp] = muDoubleScalarMax(
            -2.0,
            muDoubleScalarMin(
                2.0,
                u[qp_tmp] -
                    alf * ((u[t + 1] * (real_T)f_a[i2] +
                            u[t + 200] * (real_T)f_a[i2 + 1]) -
                           (((v[t + 1] * g_a[i3] + v[t + 200] * g_a[i3 + 1]) +
                             v[t + 399] * g_a[i3 + 2]) +
                            v[t + 598] * g_a[i3 + 3]))));
      }
    }
    /*  proj grad t in [M, N-1] */
    i = (int32_T)((N - 1.0) + (1.0 - M));
    for (t = 0; t < i; t++) {
      sig = M + (real_T)t;
      qp[(int32_T)sig - 1] = 4.0;
      qp[(int32_T)sig + 199] = 0.0;
      qp[(int32_T)sig + 399] = 4.0;
      qp[(int32_T)sig + 599] = 0.0;
      for (i1 = 0; i1 < 2; i1++) {
        i2 = i1 << 1;
        i3 = i1 << 2;
        qp_tmp = ((int32_T)sig + 199 * i1) - 1;
        up[qp_tmp] = muDoubleScalarMax(
            -2.0,
            muDoubleScalarMin(
                2.0, u[qp_tmp] -
                         alf * ((u[(int32_T)sig - 1] * (real_T)f_a[i2] +
                                 u[(int32_T)sig + 198] * (real_T)f_a[i2 + 1]) -
                                (((v[(int32_T)sig - 1] * g_a[i3] +
                                   v[(int32_T)sig + 198] * g_a[i3 + 1]) +
                                  v[(int32_T)sig + 397] * g_a[i3 + 2]) +
                                 v[(int32_T)sig + 596] * g_a[i3 + 3]))));
      }
    }
    qp[(int32_T)N - 1] = 4.0;
    qp[(int32_T)N + 199] = 0.0;
    qp[(int32_T)N + 399] = 4.0;
    qp[(int32_T)N + 599] = 0.0;
    for (i = 0; i < 800; i++) {
      q2[i] = 2.0 * qp[i] - q[i];
    }
    for (i = 0; i < 398; i++) {
      u2[i] = 2.0 * up[i] - u[i];
    }
    i = (int32_T)N;
    for (t = 0; t <= i - 2; t++) {
      sig = u2[t];
      sig1 = u2[t + 199];
      for (i1 = 0; i1 < 4; i1++) {
        i2 = i1 << 2;
        i3 = i1 << 1;
        qp_tmp = t + 199 * i1;
        vp[qp_tmp] = v[qp_tmp] +
                     bet * ((q2[(t + 200 * i1) + 1] -
                             (((q2[t] * b_a[i2] + q2[t + 200] * b_a[i2 + 1]) +
                               q2[t + 400] * b_a[i2 + 2]) +
                              q2[t + 600] * b_a[i2 + 3])) -
                            (sig * c_a[i3] + sig1 * c_a[i3 + 1]));
      }
    }
    for (i1 = 0; i1 < 800; i1++) {
      q[i1] = -0.89999999999999991 * q[i1] + 1.9 * qp[i1];
    }
    for (i1 = 0; i1 < 398; i1++) {
      u[i1] = -0.89999999999999991 * u[i1] + 1.9 * up[i1];
    }
    for (i1 = 0; i1 < 796; i1++) {
      v[i1] = -0.89999999999999991 * v[i1] + 1.9 * vp[i1];
    }
    if (muDoubleScalarRem((real_T)k + 1.0, 15.0) == 0.0) {
      sig1 = conv_test_dual;
      /*  conv_test_primal_prev = conv_test_primal; */
      sig = 0.0;
      conv_test_dual = 0.0;
      for (qp_tmp = 0; qp_tmp <= i - 2; qp_tmp++) {
        varargin_2 = 0.0;
        absx = muDoubleScalarAbs(q[qp_tmp] - qp[qp_tmp]);
        if (muDoubleScalarIsNaN(absx) || (absx > 0.0)) {
          varargin_2 = absx;
        }
        absx = muDoubleScalarAbs(q[qp_tmp + 200] - qp[qp_tmp + 200]);
        if (muDoubleScalarIsNaN(absx) || (absx > varargin_2)) {
          varargin_2 = absx;
        }
        absx = muDoubleScalarAbs(q[qp_tmp + 400] - qp[qp_tmp + 400]);
        if (muDoubleScalarIsNaN(absx) || (absx > varargin_2)) {
          varargin_2 = absx;
        }
        absx = muDoubleScalarAbs(q[qp_tmp + 600] - qp[qp_tmp + 600]);
        if (muDoubleScalarIsNaN(absx) || (absx > varargin_2)) {
          varargin_2 = absx;
        }
        b_varargin_2 = 0.0;
        absx = muDoubleScalarAbs(u[qp_tmp] - up[qp_tmp]);
        if (muDoubleScalarIsNaN(absx) || (absx > 0.0)) {
          b_varargin_2 = absx;
        }
        absx = muDoubleScalarAbs(u[qp_tmp + 199] - up[qp_tmp + 199]);
        if (muDoubleScalarIsNaN(absx) || (absx > b_varargin_2)) {
          b_varargin_2 = absx;
        }
        sig =
            muDoubleScalarMax(muDoubleScalarMax(sig, varargin_2), b_varargin_2);
        varargin_2 = 0.0;
        absx = muDoubleScalarAbs(v[qp_tmp] - vp[qp_tmp]);
        if (muDoubleScalarIsNaN(absx) || (absx > 0.0)) {
          varargin_2 = absx;
        }
        absx = muDoubleScalarAbs(v[qp_tmp + 199] - vp[qp_tmp + 199]);
        if (muDoubleScalarIsNaN(absx) || (absx > varargin_2)) {
          varargin_2 = absx;
        }
        absx = muDoubleScalarAbs(v[qp_tmp + 398] - vp[qp_tmp + 398]);
        if (muDoubleScalarIsNaN(absx) || (absx > varargin_2)) {
          varargin_2 = absx;
        }
        absx = muDoubleScalarAbs(v[qp_tmp + 597] - vp[qp_tmp + 597]);
        if (muDoubleScalarIsNaN(absx) || (absx > varargin_2)) {
          varargin_2 = absx;
        }
        conv_test_dual = muDoubleScalarMax(conv_test_dual, varargin_2);
      }
      varargin_2 = 0.0;
      absx = muDoubleScalarAbs(q[(int32_T)N - 1] - qp[(int32_T)N - 1]);
      if (muDoubleScalarIsNaN(absx) || (absx > 0.0)) {
        varargin_2 = absx;
      }
      absx = muDoubleScalarAbs(q[(int32_T)N + 199] - qp[(int32_T)N + 199]);
      if (muDoubleScalarIsNaN(absx) || (absx > varargin_2)) {
        varargin_2 = absx;
      }
      absx = muDoubleScalarAbs(q[(int32_T)N + 399] - qp[(int32_T)N + 399]);
      if (muDoubleScalarIsNaN(absx) || (absx > varargin_2)) {
        varargin_2 = absx;
      }
      absx = muDoubleScalarAbs(q[(int32_T)N + 599] - qp[(int32_T)N + 599]);
      if (muDoubleScalarIsNaN(absx) || (absx > varargin_2)) {
        varargin_2 = absx;
      }
      sig = muDoubleScalarMax(sig, varargin_2);
      if ((sig <= 0.005) && (conv_test_dual <= 0.005)) {
        conv_flag = true;
        *exitflag = 1.0;
        mexPrintf("PRIMAL DUAL FEASIBLE\nPIPG converged in %.0f iterations.\n",
                  (real_T)k + 1.0);
        exitg1 = true;
      } else if (muDoubleScalarAbs(conv_test_dual / sig1 - 1.0) <= 0.05) {
        conv_flag = true;
        *exitflag = -1.0;
        mexPrintf("PRIMAL INFEASIBLE\nPIPG converged in %.0f iterations.\n",
                  (real_T)k + 1.0);
        exitg1 = true;
      } else {
        k++;
      }
    } else {
      k++;
    }
  }
  *solve_time = toc();
  /*  end epipg timing */
  if (!conv_flag) {
    mexPrintf("INCONCLUSIVE\nMaximum iterations reached.\n");
  }
}

/* End of code generation (epipg.c) */
