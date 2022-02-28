/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_epipg_api.c
 *
 * Code generation for function '_coder_epipg_api'
 *
 */

/* Include files */
#include "_coder_epipg_api.h"
#include "epipg.h"
#include "epipg_data.h"
#include "rt_nonfinite.h"

/* Function Declarations */
static void b_emlrt_marshallIn(const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               real_T y[800]);

static const mxArray *b_emlrt_marshallOut(const real_T u[398]);

static void c_emlrt_marshallIn(const mxArray *u, const char_T *identifier,
                               real_T y[398]);

static const mxArray *c_emlrt_marshallOut(const real_T u);

static void d_emlrt_marshallIn(const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               real_T y[398]);

static const mxArray *d_emlrt_marshallOut(const real_T u[796]);

static void e_emlrt_marshallIn(const mxArray *v, const char_T *identifier,
                               real_T y[796]);

static void emlrt_marshallIn(const mxArray *q, const char_T *identifier,
                             real_T y[800]);

static const mxArray *emlrt_marshallOut(const real_T u[800]);

static void f_emlrt_marshallIn(const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               real_T y[796]);

static real_T g_emlrt_marshallIn(const mxArray *omg, const char_T *identifier);

static real_T h_emlrt_marshallIn(const mxArray *u,
                                 const emlrtMsgIdentifier *parentId);

static void i_emlrt_marshallIn(const mxArray *src,
                               const emlrtMsgIdentifier *msgId,
                               real_T ret[800]);

static void j_emlrt_marshallIn(const mxArray *src,
                               const emlrtMsgIdentifier *msgId,
                               real_T ret[398]);

static void k_emlrt_marshallIn(const mxArray *src,
                               const emlrtMsgIdentifier *msgId,
                               real_T ret[796]);

static real_T l_emlrt_marshallIn(const mxArray *src,
                                 const emlrtMsgIdentifier *msgId);

/* Function Definitions */
static void b_emlrt_marshallIn(const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               real_T y[800])
{
  real_T dv[800];
  int32_T i;
  int32_T i1;
  i_emlrt_marshallIn(emlrtAlias(u), parentId, dv);
  for (i = 0; i < 4; i++) {
    for (i1 = 0; i1 < 200; i1++) {
      y[i1 + 200 * i] = dv[i + (i1 << 2)];
    }
  }
  emlrtDestroyArray(&u);
}

static const mxArray *b_emlrt_marshallOut(const real_T u[398])
{
  static const int32_T iv[2] = {2, 199};
  const mxArray *m;
  const mxArray *y;
  real_T dv[398];
  real_T *pData;
  int32_T b_i;
  int32_T i;
  int32_T i1;
  y = NULL;
  for (i = 0; i < 199; i++) {
    i1 = i << 1;
    dv[i1] = u[i];
    dv[i1 + 1] = u[i + 199];
  }
  m = emlrtCreateNumericArray(2, (const void *)&iv[0], mxDOUBLE_CLASS, mxREAL);
  pData = emlrtMxGetPr(m);
  i = 0;
  for (b_i = 0; b_i < 199; b_i++) {
    i1 = b_i << 1;
    pData[i] = dv[i1];
    pData[i + 1] = dv[i1 + 1];
    i += 2;
  }
  emlrtAssign(&y, m);
  return y;
}

static void c_emlrt_marshallIn(const mxArray *u, const char_T *identifier,
                               real_T y[398])
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  d_emlrt_marshallIn(emlrtAlias(u), &thisId, y);
  emlrtDestroyArray(&u);
}

static const mxArray *c_emlrt_marshallOut(const real_T u)
{
  const mxArray *m;
  const mxArray *y;
  y = NULL;
  m = emlrtCreateDoubleScalar(u);
  emlrtAssign(&y, m);
  return y;
}

static void d_emlrt_marshallIn(const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               real_T y[398])
{
  real_T dv[398];
  int32_T i;
  int32_T i1;
  j_emlrt_marshallIn(emlrtAlias(u), parentId, dv);
  for (i = 0; i < 2; i++) {
    for (i1 = 0; i1 < 199; i1++) {
      y[i1 + 199 * i] = dv[i + (i1 << 1)];
    }
  }
  emlrtDestroyArray(&u);
}

static const mxArray *d_emlrt_marshallOut(const real_T u[796])
{
  static const int32_T iv[2] = {4, 199};
  const mxArray *m;
  const mxArray *y;
  real_T dv[796];
  real_T *pData;
  int32_T b_i;
  int32_T i;
  int32_T i1;
  y = NULL;
  for (i = 0; i < 199; i++) {
    i1 = i << 2;
    dv[i1] = u[i];
    dv[i1 + 1] = u[i + 199];
    dv[i1 + 2] = u[i + 398];
    dv[i1 + 3] = u[i + 597];
  }
  m = emlrtCreateNumericArray(2, (const void *)&iv[0], mxDOUBLE_CLASS, mxREAL);
  pData = emlrtMxGetPr(m);
  i = 0;
  for (b_i = 0; b_i < 199; b_i++) {
    i1 = b_i << 2;
    pData[i] = dv[i1];
    pData[i + 1] = dv[i1 + 1];
    pData[i + 2] = dv[i1 + 2];
    pData[i + 3] = dv[i1 + 3];
    i += 4;
  }
  emlrtAssign(&y, m);
  return y;
}

static void e_emlrt_marshallIn(const mxArray *v, const char_T *identifier,
                               real_T y[796])
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  f_emlrt_marshallIn(emlrtAlias(v), &thisId, y);
  emlrtDestroyArray(&v);
}

static void emlrt_marshallIn(const mxArray *q, const char_T *identifier,
                             real_T y[800])
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  b_emlrt_marshallIn(emlrtAlias(q), &thisId, y);
  emlrtDestroyArray(&q);
}

static const mxArray *emlrt_marshallOut(const real_T u[800])
{
  static const int32_T iv[2] = {4, 200};
  const mxArray *m;
  const mxArray *y;
  real_T dv[800];
  real_T *pData;
  int32_T b_i;
  int32_T i;
  int32_T i1;
  y = NULL;
  for (i = 0; i < 200; i++) {
    i1 = i << 2;
    dv[i1] = u[i];
    dv[i1 + 1] = u[i + 200];
    dv[i1 + 2] = u[i + 400];
    dv[i1 + 3] = u[i + 600];
  }
  m = emlrtCreateNumericArray(2, (const void *)&iv[0], mxDOUBLE_CLASS, mxREAL);
  pData = emlrtMxGetPr(m);
  i = 0;
  for (b_i = 0; b_i < 200; b_i++) {
    i1 = b_i << 2;
    pData[i] = dv[i1];
    pData[i + 1] = dv[i1 + 1];
    pData[i + 2] = dv[i1 + 2];
    pData[i + 3] = dv[i1 + 3];
    i += 4;
  }
  emlrtAssign(&y, m);
  return y;
}

static void f_emlrt_marshallIn(const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               real_T y[796])
{
  real_T dv[796];
  int32_T i;
  int32_T i1;
  k_emlrt_marshallIn(emlrtAlias(u), parentId, dv);
  for (i = 0; i < 4; i++) {
    for (i1 = 0; i1 < 199; i1++) {
      y[i1 + 199 * i] = dv[i + (i1 << 2)];
    }
  }
  emlrtDestroyArray(&u);
}

static real_T g_emlrt_marshallIn(const mxArray *omg, const char_T *identifier)
{
  emlrtMsgIdentifier thisId;
  real_T y;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = h_emlrt_marshallIn(emlrtAlias(omg), &thisId);
  emlrtDestroyArray(&omg);
  return y;
}

static real_T h_emlrt_marshallIn(const mxArray *u,
                                 const emlrtMsgIdentifier *parentId)
{
  real_T y;
  y = l_emlrt_marshallIn(emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static void i_emlrt_marshallIn(const mxArray *src,
                               const emlrtMsgIdentifier *msgId, real_T ret[800])
{
  static const int32_T dims[2] = {4, 200};
  real_T(*r)[800];
  int32_T i;
  emlrtCheckBuiltInR2012b(emlrtRootTLSGlobal, msgId, src,
                          (const char_T *)"double", false, 2U,
                          (void *)&dims[0]);
  r = (real_T(*)[800])emlrtMxGetData(src);
  for (i = 0; i < 800; i++) {
    ret[i] = (*r)[i];
  }
  emlrtDestroyArray(&src);
}

static void j_emlrt_marshallIn(const mxArray *src,
                               const emlrtMsgIdentifier *msgId, real_T ret[398])
{
  static const int32_T dims[2] = {2, 199};
  real_T(*r)[398];
  int32_T i;
  emlrtCheckBuiltInR2012b(emlrtRootTLSGlobal, msgId, src,
                          (const char_T *)"double", false, 2U,
                          (void *)&dims[0]);
  r = (real_T(*)[398])emlrtMxGetData(src);
  for (i = 0; i < 398; i++) {
    ret[i] = (*r)[i];
  }
  emlrtDestroyArray(&src);
}

static void k_emlrt_marshallIn(const mxArray *src,
                               const emlrtMsgIdentifier *msgId, real_T ret[796])
{
  static const int32_T dims[2] = {4, 199};
  real_T(*r)[796];
  int32_T i;
  emlrtCheckBuiltInR2012b(emlrtRootTLSGlobal, msgId, src,
                          (const char_T *)"double", false, 2U,
                          (void *)&dims[0]);
  r = (real_T(*)[796])emlrtMxGetData(src);
  for (i = 0; i < 796; i++) {
    ret[i] = (*r)[i];
  }
  emlrtDestroyArray(&src);
}

static real_T l_emlrt_marshallIn(const mxArray *src,
                                 const emlrtMsgIdentifier *msgId)
{
  static const int32_T dims = 0;
  real_T ret;
  emlrtCheckBuiltInR2012b(emlrtRootTLSGlobal, msgId, src,
                          (const char_T *)"double", false, 0U, (void *)&dims);
  ret = *(real_T *)emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

void epipg_api(const mxArray *const prhs[26], int32_T nlhs,
               const mxArray *plhs[5])
{
  real_T q[800];
  real_T v[796];
  real_T u[398];
  real_T M;
  real_T N;
  real_T exitflag;
  real_T omg;
  real_T solve_time;
  /* Marshall function inputs */
  emlrt_marshallIn(emlrtAliasP(prhs[0]), "q", q);
  c_emlrt_marshallIn(emlrtAliasP(prhs[1]), "u", u);
  e_emlrt_marshallIn(emlrtAliasP(prhs[2]), "v", v);
  omg = g_emlrt_marshallIn(emlrtAliasP(prhs[10]), "omg");
  N = g_emlrt_marshallIn(emlrtAliasP(prhs[12]), "N");
  M = g_emlrt_marshallIn(emlrtAliasP(prhs[13]), "M");
  /* Invoke the target function */
  epipg(q, u, v, omg, N, M, &exitflag, &solve_time);
  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(q);
  if (nlhs > 1) {
    plhs[1] = b_emlrt_marshallOut(u);
  }
  if (nlhs > 2) {
    plhs[2] = c_emlrt_marshallOut(exitflag);
  }
  if (nlhs > 3) {
    plhs[3] = d_emlrt_marshallOut(v);
  }
  if (nlhs > 4) {
    plhs[4] = c_emlrt_marshallOut(solve_time);
  }
}

/* End of code generation (_coder_epipg_api.c) */
