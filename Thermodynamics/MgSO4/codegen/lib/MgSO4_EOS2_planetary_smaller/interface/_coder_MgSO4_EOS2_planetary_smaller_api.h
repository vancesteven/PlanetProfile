/*
 * _coder_MgSO4_EOS2_planetary_smaller_api.h
 *
 * Code generation for function '_coder_MgSO4_EOS2_planetary_smaller_api'
 *
 */

#ifndef _CODER_MGSO4_EOS2_PLANETARY_SMALLER_API_H
#define _CODER_MGSO4_EOS2_PLANETARY_SMALLER_API_H

/* Include files */
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include <stddef.h>
#include <stdlib.h>
#include "_coder_MgSO4_EOS2_planetary_smaller_api.h"

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

/* Function Declarations */
extern void MgSO4_EOS2_planetary_smaller(real_T m, real_T P, real_T T, real_T
  *rho, real_T *Cp, real_T *alpha);
extern void MgSO4_EOS2_planetary_smaller_api(const mxArray * const prhs[3],
  const mxArray *plhs[3]);
extern void MgSO4_EOS2_planetary_smaller_atexit(void);
extern void MgSO4_EOS2_planetary_smaller_initialize(void);
extern void MgSO4_EOS2_planetary_smaller_terminate(void);
extern void MgSO4_EOS2_planetary_smaller_xil_terminate(void);

#endif

/* End of code generation (_coder_MgSO4_EOS2_planetary_smaller_api.h) */
