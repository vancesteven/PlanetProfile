/*
 * _coder_MgSO4_EOS2_planetary_smaller_mex.cpp
 *
 * Code generation for function '_coder_MgSO4_EOS2_planetary_smaller_mex'
 *
 */

/* Include files */
#include "_coder_MgSO4_EOS2_planetary_smaller_api.h"
#include "_coder_MgSO4_EOS2_planetary_smaller_mex.h"

/* Function Declarations */
static void c_MgSO4_EOS2_planetary_smaller_(int32_T nlhs, mxArray *plhs[3],
  int32_T nrhs, const mxArray *prhs[3]);

/* Function Definitions */
static void c_MgSO4_EOS2_planetary_smaller_(int32_T nlhs, mxArray *plhs[3],
  int32_T nrhs, const mxArray *prhs[3])
{
  int32_T n;
  const mxArray *inputs[3];
  const mxArray *outputs[3];
  int32_T b_nlhs;
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;

  /* Check for proper number of arguments. */
  if (nrhs != 3) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 3, 4,
                        28, "MgSO4_EOS2_planetary_smaller");
  }

  if (nlhs > 3) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 28,
                        "MgSO4_EOS2_planetary_smaller");
  }

  /* Temporary copy for mex inputs. */
  for (n = 0; n < nrhs; n++) {
    inputs[n] = prhs[n];
  }

  /* Call the function. */
  MgSO4_EOS2_planetary_smaller_api(inputs, outputs);

  /* Copy over outputs to the caller. */
  if (nlhs < 1) {
    b_nlhs = 1;
  } else {
    b_nlhs = nlhs;
  }

  emlrtReturnArrays(b_nlhs, plhs, outputs);

  /* Module termination. */
  MgSO4_EOS2_planetary_smaller_terminate();
}

void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs, const mxArray
                 *prhs[])
{
  mexAtExit(MgSO4_EOS2_planetary_smaller_atexit);

  /* Initialize the memory manager. */
  /* Module initialization. */
  MgSO4_EOS2_planetary_smaller_initialize();

  /* Dispatch the entry-point. */
  c_MgSO4_EOS2_planetary_smaller_(nlhs, plhs, nrhs, prhs);
}

emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1);
  return emlrtRootTLSGlobal;
}

/* End of code generation (_coder_MgSO4_EOS2_planetary_smaller_mex.cpp) */
