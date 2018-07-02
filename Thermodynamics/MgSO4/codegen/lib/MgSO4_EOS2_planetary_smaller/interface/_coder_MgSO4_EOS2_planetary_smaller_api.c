/*
 * _coder_MgSO4_EOS2_planetary_smaller_api.c
 *
 * Code generation for function '_coder_MgSO4_EOS2_planetary_smaller_api'
 *
 */

/* Include files */
#include "tmwtypes.h"
#include "_coder_MgSO4_EOS2_planetary_smaller_api.h"
#include "_coder_MgSO4_EOS2_planetary_smaller_mex.h"

/* Variable Definitions */
emlrtCTX emlrtRootTLSGlobal = NULL;
emlrtContext emlrtContextGlobal = { true,/* bFirstTime */
  false,                               /* bInitialized */
  131435U,                             /* fVersionInfo */
  NULL,                                /* fErrorFunction */
  "MgSO4_EOS2_planetary_smaller",      /* fFunctionName */
  NULL,                                /* fRTCallStack */
  false,                               /* bDebugMode */
  { 2045744189U, 2170104910U, 2743257031U, 4284093946U },/* fSigWrd */
  NULL                                 /* fSigMem */
};

/* Function Declarations */
static real_T b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId);
static real_T c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId);
static real_T emlrt_marshallIn(const emlrtStack *sp, const mxArray *m, const
  char_T *identifier);
static const mxArray *emlrt_marshallOut(const real_T u);

/* Function Definitions */
static real_T b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId)
{
  real_T y;
  y = c_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static real_T c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId)
{
  real_T ret;
  static const int32_T dims = 0;
  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 0U, &dims);
  ret = *(real_T *)mxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static real_T emlrt_marshallIn(const emlrtStack *sp, const mxArray *m, const
  char_T *identifier)
{
  real_T y;
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = b_emlrt_marshallIn(sp, emlrtAlias(m), &thisId);
  emlrtDestroyArray(&m);
  return y;
}

static const mxArray *emlrt_marshallOut(const real_T u)
{
  const mxArray *y;
  const mxArray *m0;
  y = NULL;
  m0 = emlrtCreateDoubleScalar(u);
  emlrtAssign(&y, m0);
  return y;
}

void MgSO4_EOS2_planetary_smaller_api(const mxArray * const prhs[3], const
  mxArray *plhs[3])
{
  real_T m;
  real_T P;
  real_T T;
  real_T rho;
  real_T Cp;
  real_T alpha;
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;

  /* Marshall function inputs */
  m = emlrt_marshallIn(&st, emlrtAliasP((const mxArray *)prhs[0]), "m");
  P = emlrt_marshallIn(&st, emlrtAliasP((const mxArray *)prhs[1]), "P");
  T = emlrt_marshallIn(&st, emlrtAliasP((const mxArray *)prhs[2]), "T");

  /* Invoke the target function */
  MgSO4_EOS2_planetary_smaller(m, P, T, &rho, &Cp, &alpha);

  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(rho);
  plhs[1] = emlrt_marshallOut(Cp);
  plhs[2] = emlrt_marshallOut(alpha);
}

void MgSO4_EOS2_planetary_smaller_atexit(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtEnterRtStackR2012b(&st);
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
  MgSO4_EOS2_planetary_smaller_xil_terminate();
}

void MgSO4_EOS2_planetary_smaller_initialize(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtClearAllocCountR2012b(&st, false, 0U, 0);
  emlrtEnterRtStackR2012b(&st);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

void MgSO4_EOS2_planetary_smaller_terminate(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/* End of code generation (_coder_MgSO4_EOS2_planetary_smaller_api.c) */
