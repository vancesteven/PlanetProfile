/*
 * TensorGriddedInterp.cpp
 *
 * Code generation for function 'TensorGriddedInterp'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "MgSO4_EOS2_planetary_smaller.h"
#include "TensorGriddedInterp.h"
#include "ppval.h"
#include "spline.h"

/* Function Definitions */
double TensorGriddedInterp(const double varargin_1[14], const double varargin_2
  [20], const double varargin_3[15], const double varargin_4[4200], double
  varargin_5, double varargin_6, double varargin_7)
{
  double ppk_breaks[15];
  static double ppk_coefs[15680];
  double vkj[280];
  double varargout_2[280];
  double b_ppk_breaks[20];
  double b_ppk_coefs[1064];
  double b_vkj[14];
  double varargout_1[14];
  double c_ppk_breaks[14];
  double c_ppk_coefs[52];
  spline(varargin_3, varargin_4, ppk_breaks, ppk_coefs);
  ppval(ppk_breaks, ppk_coefs, varargin_7, vkj);
  memcpy(&varargout_2[0], &vkj[0], 280U * sizeof(double));
  b_spline(varargin_2, varargout_2, b_ppk_breaks, b_ppk_coefs);
  b_ppval(b_ppk_breaks, b_ppk_coefs, varargin_6, b_vkj);
  memcpy(&varargout_1[0], &b_vkj[0], 14U * sizeof(double));
  c_spline(varargin_1, varargout_1, c_ppk_breaks, c_ppk_coefs);
  return c_ppval(c_ppk_breaks, c_ppk_coefs, varargin_5);
}

/* End of code generation (TensorGriddedInterp.cpp) */
